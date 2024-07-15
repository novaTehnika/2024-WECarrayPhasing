% study_hydElecPTO_checkValves.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/8/2024
%
% PURPOSE/DESCRIPTION:
% This script performs parameter variation study
% using the model contained in sys_hydElecPTO.m and solved by
% sim_hydElecPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_hydElecPTO.m script is called.
%
% This specific script studies the size of check valves at the WEC-driven
% pump.
%
% This script is set up to be run as part of a SLURM job array. The
% following lines are required before this script is called:
%   iVar = ${SLURM_ARRAY_TASK_ID};
%   SS=1;
%
% FILE DEPENDENCY:
% ./Hydraulic-Electric PTO/
%   initialConditionDefault_hydElecPTO.m
%   parameters_hydElecPTO.m
%   sim_hydElecPTO.m
%   stateIndex_hydElecPTO.m
%   sys_hydElecPTO.m
% ./WEC model/
%   flapModel.m
%   hydroStaticTorque.m
%   parameters_WECmodel.m
% ./WEC model/WECdata
%   nemohResults_vantHoff2009_20180802.mat
%   vantHoffTFCoeff.mat
% ./Solvers/
%   deltaE_NI.m
%   deltaV_NI.m
%   ode1.m
% ./Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
% ./Components/Pipeline
%   flowR.m
%   lineCap.m
%   pipelineNPi.m
% ./Utilities/
%   startParPool.m
%   statsTimeVar_cdf.m
%   get_current_git_hash.m
%   leadingZeros.m
%
% UPDATES:
% 7/8/2024 - Created from study_hydElecPTO_arrayHPaccum.m and 
% study_parPTO_checkValves.m in the repository 
% https://github.com/novaTehnika/2023-DynPTOModelDesignStudies
%
% Copyright (C) 2024  Jeremy W. Simmons II
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program. If not, see <https://www.gnu.org/licenses/>.
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear
% clc
addpath('WEC model')
addpath(['WEC model' filesep 'WECdata'])
addpath('Hydraulic-Electric PTO')
addpath('Components')
addpath(['Components' filesep 'Pipeline'])
addpath('Sea States')
addpath('Solvers')
addpath('Utilities')
[git_hash_string, git_status_string] = get_current_git_hash();

%% %%%%%%%%%%%%   SIMULATION PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation timeframe
par.tstart = 0; %[s] start time of simulation
par.tend = 500; %[s] end time of simulation

par.Tramp = 250; % [s] excitation force ramp period
par.TrampWEC = min(25,par.Tramp); % [s] excitation force ramp period

% Solver parameters
par.solver = 'fixed time'; % 'variable time' OR 'fixed time'
switch par.solver
    case 'fixed time'
        par.MaxStep = 1e-5;             % [s] time step size
        par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
        if mod(par.downSampledStepSize,par.MaxStep)
            warning('down-sampled time step is not an integer multiple of the maximum step size')
        end
    case 'variable time'
        par.odeSolverRelTol = 1e-4; % Rel. error tolerance parameter for ODE solver
        par.odeSolverAbsTol = 1e-4; % Abs. error tolerance parameter for ODE solver
        par.MaxStep = 5e-2;             % [s] max step size for variable timestep solver
        par.stepSizeWECramp = 1e-4;     % [s] step size for solver during WEC ramp
        par.downSampledStepSize = 1e-2; % [s] specifies time step for data output
        if mod(par.downSampledStepSize,par.MaxStep)
            warning('down-sampled time step is not an integer multiple of the maximum step size')
        end
end

% Sea State and Wave construction parameters
switch 1
    case 1
        load('Sea States/SSdata_Bull2017WEPrize.mat')
    case 2
        load('Sea States/SSdata_HumboltBay_1D.mat')
end
SS = 2;
par.wave.Hs = Hs(SS);
par.wave.Tp = Tp(SS);
par.wave.waveDirection = 0; % [rad]
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% number of WECs and their positions
par.NumWECs = 1;

 % location perpendicular to 0 degree wave direction
WECspacing = 30; % [m]
par.WEC.y = (0:par.NumWECs-1)*WECspacing;

 % location parallel to 0 degree wave direction
par.WEC.x = [0, 10];

% load parameters (must come after NumWECs and WEC pos. are specified)
par = parameters_hydElecPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

%% Special modifications to base parameters
par.control.p_nom = 8e6; % [Pa]
% par.w_c = (2500)*2*pi/60; % [(rpm) -> rad/s] Charge pump speed
par.control.p_l_nom = 0.5e6; % [Pa]

par.motor.D = (1000)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad] Motor displacement per WEC

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nVar = 50;
kv = logspace(log10(0.85e-4),log10(1e-3),nVar);% [m^3/s/Pa] valve coefficient for high-pressure outlet check valve
X = 1.5; % proportion between low and high-pressure check valves

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% Set study variables
% change design parameter
par.kvWECout = kv(iVar);
par.kvWECin = X*kv(iVar);

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define state indices
par.iy = stateIndex_hydElecPTO(par);

% Define initial conditions
y0 = initialConditionDefault_hydElecPTO(par); % default ICs, provides 'y0'

% run simulation
ticSIM = tic;
out = sim_hydElecPTO(y0,par);
toc(ticSIM)

% Calculate metrics
eff_wecPump = mean(out.power.P_wp) ...
                    /mean(out.power.P_WEC);
p_min_wp = min([out.p_a;out.p_b]);

if ~saveSimData
    clear out
end

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolobj = gcp('nocreate'); delete(poolobj);

%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datetime("now",'format','yyyy-MM-dd''T''HH:mm'); % time in ISO8601

% Save data
filename = ['data_hydElecPTO_checkValves', ...
            '_',char(datetime("now",'Format','yyyyMMdd')), ...
            '_',num2str(SS,leadingZeros(999)), ...
            '_',num2str(iVar,leadingZeros(nVar))];
save(filename,'-v7.3')

return
