% study_hydElecPTO_arrayHPaccum.m script m-file
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
% This specific script studies the total high-pressure accumulator volume
% and number of WECs within an array.
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
% 7/8/2024 - Created from study_parPTO_accum_wActiveRV.m in the
% repository https://github.com/novaTehnika/2023-DynPTOModelDesignStudies
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
par.tend = 2000; %[s] end time of simulation

par.Tramp = 250; % [s] excitation force ramp period
par.TrampWEC = min(25,par.Tramp); % [s] excitation force ramp period

% Solver parameters
par.solver = 'fixed time'; % 'variable time' OR 'fixed time'
switch par.solver
    case 'fixed time'
        par.MaxStep = 5e-5;             % [s] time step size
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

% load parameters
par = parameters_hydElecPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of WECs
NumWECs = 1:10;
nVar1 = numel(NumWECs);

% total accumulator volume
nVar2 = 10;
VperWEC = 1e-3*logspace(log10(10),log10(2e3),nVar2);% [L->m^3] total accumulator volume

% nominal pressure
nVar3 = 20;
p_nom = linspace(10e6,30e6,nVar3); % [Pa]

[meshVar.NumWECs, meshVar.VperWEC, meshVar.p_nom] = ...
                                    meshgrid(NumWECs,VperWEC,p_nom);
NumWECs_mesh = meshVar.NumWECs(:);
VperWEC_mesh = meshVar.VperWEC(:);
p_nom_mesh = meshVar.p_nom(:);

nVar = length(NumWECs_mesh);

saveSimData = 1; % save simulation data (1) or just output variables (0)

%% Special modifications to base parameters
par.control.p_l_nom = 0.5e6; % [Pa]

D_m_base = (1000)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad] Motor displacement per WEC

%% Set study variables
% nominal pressure
par.control.p_h_nom = p_nom_mesh(iVar); % [Pa]

% number of WECs and their positions
par.NumWECs = NumWECs_mesh(iVar);

 % location perpendicular to 0 degree wave direction
WECspacing = 30; % [m]
par.WEC.y = (0:par.NumWECs-1)*WECspacing; 

 % location parallel to 0 degree wave direction
  % wavelength for specified energy period
    Te_design = 9.86; % [s]
    theta_design = 0;
    K_design = fzero(@(K) K*par.WEC.g*tanh(K*par.WEC.H) - (2*pi/Te_design)^2, ...
                                        (2*pi/Te_design)^2/par.WEC.g);
    k_design = K_design*cos(theta_design); % x component of wave number vector
    waveLen_design = 2*pi/k_design;

  % phase periodicity
  phase_periodicity = pi;

  %
  rng(3); % seed random number generator
  phase_spacing = phase_periodicity/par.NumWECs;
  lin_spacing = phase_spacing/2/pi*waveLen_design;

  x = ((1:par.NumWECs)-1)*lin_spacing - lin_spacing/2;
  par.WEC.x = x(randperm(par.NumWECs));

% update WEC model parameters based on modified num. WECs and pos.
    par = parameters_WECmodel(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

% accumulator volume
par.Vc_h = NumWECs_mesh(iVar)*VperWEC_mesh(iVar);

% hydraulic motor displacement
par.motor.D = D_m_base*par.NumWECs;

% Pressure relief valve
 % high-pressure outlet of WEC-driven pump
    maxPressure = 30e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = (100)*1e-3*par.NumWECs; % [(L/s) -> m^3/s]
    par.hPRV.p_crack = maxPressure - margin;
    par.hPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

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
% max rate of change in pressure
% 97th percentile ratof change
% power loss from valve
% power loss through PRVs
% permeate production
% power loss from pump/motor and power generated for normalization
% 
PP_WEC = mean(out.power.P_WEC);
PP_wp = mean(out.power.P_wp);
PP_mLoss = mean(out.power.P_mLoss);
PP_gen = mean(out.power.P_gen);
PP_hPRV = mean(out.power.P_hPRV);
dpdt_max = max(abs(out.dydt(:,par.iy.p_h)));

if ~saveSimData
    clear out
end

%% %%%%%%%%%%%%   End Computations  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

poolobj = gcp('nocreate'); delete(poolobj);

%% %%%%%%%%%%%%   Save Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
timeStamp = datetime("now",'format','yyyy-MM-dd''T''HH:mm'); % time in ISO8601

% Save data
filename = ['data_hydElecPTO_arrayHPaccum', ...
            '_',char(datetime("now",'Format','yyyyMMdd')), ...
            '_',num2str(SS,leadingZeros(999)), ...
            '_',num2str(iVar,leadingZeros(nVar))];
save(filename,'-v7.3')

return
