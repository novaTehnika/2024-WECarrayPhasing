% run_hydElecPTO.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This script serves as a shell for running a single simulation
% using the model contained in sys_parPTO.m and solved by 
% sim_hydElecPTO.m.
% The parameter initiallization functions are called within this
% script before the sim_parPTO.m script is called.
%
% FILE DEPENDENCY:
% ./Hydraulic-Electric PTO/
%   initialConditionDefault_hydElecPTO
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
%
% UPDATES:
% 6/7/2024 - Created from run_parPTO.m.
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
clear
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
par.WEC.nw = 1000; % num. of frequency components for harmonic superposition 
par.wave.rngSeedPhase = 3; % seed for the random number generator

% load parameters
par = parameters_parPTO(par,...
    'nemohResults_vantHoff2009_20180802.mat','vantHoffTFCoeff.mat');

%% Special modifications to base parameters
% par.Sro = 3700; % [m^3]
% par.D_WEC = 0.23;         % [m^3/rad] flap pump displacement

% Operating parameters
% p_ro_nom = 1e6*[4.0000 4.9435 8.0000 5.2661 8.0000 7.1052]; % [Pa]
% w_c = [2500 2500 3000 3000 3000 3000]*2*pi/60; % [(rpm) -> rad/s]

par.control.p_ro_nom = 6.11e6; % [Pa]
par.w_c = (2500)*2*pi/60; % [(rpm) -> rad/s] Charge pump speed
par.duty_sv = 0.0;

% Configuration
par.ERUconfig.present = 1;

par.rvConfig.included = 0; % RO inlet valve is 1 - present, 0 - absent
par.rvConfig.active = (0)*par.rvConfig.included; % RO inlet valve is 1 - active, 0 - passive

par.plConfig.included = 0;


%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define state indices
par.iy = stateIndex_hydElecPTO(par);

% Define initial conditions
y0 = initialConditionDefault_hydElecPTO(par); % default ICs, provides 'y0'

% Simulation
ticSIM = tic;
out = sim_hydElecPTO(y0,par);
toc(ticSIM)

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WEC pos., vel., and Torque
bottomEdge = 1;
leftEdge = 2;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];


ax(1) = subplot(3,1,1);
plot(out.t,out.waveElev)
xlabel('time (s)')
ylabel('elevation (m)')
title('Wave Elevation')

ax(2) = subplot(3,1,2);
xlabel('time (s)')
hold on

yyaxis left
plot(out.t,out.theta)
xlabel('time (s)')
ylabel('position (rad)')
% ylim([-pi/2 pi/2])

yyaxis right
plot(out.t,out.theta_dot)
ylabel('angular velocity (rad/s)')
% ylim(10*[-pi/2 pi/2])

ax(3) = subplot(3,1,3);
hold on
plot(out.t,1e-6*out.T_wave)
plot(out.t,1e-6*out.T_pto)
plot(out.t,1e-6*out.T_rad)
plot(out.t,1e-6*out.T_hydroStatic)
ylabel('Torque (MNm)')

legend('T_{w}','T_{PTO}','T_{r}','T_{h}')
      
linkaxes(ax,'x');

sgtitle('WEC Behaviour')

%% Pressure
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 9;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(3,1,1);
plot(out.t,1e-6*out.p_hin)
hold on
plot(out.t,1e-6*out.p_hout)
plot(out.t,1e-6*out.p_ro)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{hin}','p_{hout}','p_{ro}')

ax(2) = subplot(3,1,2);
plot(out.t,1e-6*out.p_lin)
hold on
plot(out.t,1e-6*out.p_lout)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{lin}','p_{lout}')

ax(3) = subplot(3,1,3);
plot(out.t,1e-3*out.dydt(:,par.iy.p_ro))
hold on

addpath('Utilities')
dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,par.iy.p_ro)));
dpdt_97 = dist_dpdt.xi(find(dist_dpdt.f > 0.97,1,'first'));
plot(out.t([1 end]),1e-3*dpdt_97*[1 1],'-.k')
p(1) = plot(out.t([1 end]),1e-3*dpdt_97*[-1 -1],'-.k');
p(1).HandleVisibility='off';

dist_dpdt = statsTimeVar_cdf(out.t,abs(out.dydt(:,par.iy.p_ro)));
dpdt_99 = dist_dpdt.xi(find(dist_dpdt.f > 0.99,1,'first'));
plot(out.t([1 end]),1e-3*dpdt_99*[1 1],'-.k')
p(2) = plot(out.t([1 end]),1e-3*dpdt_99*[-1 -1],'--k');
p(2).HandleVisibility='off';

plot(out.t([1 end]),1e-3*out.par.control.dpdt_ROmax*[1 1],'-r')
p(3) = plot(out.t([1 end]),1e-3*out.par.control.dpdt_ROmax*[-1 -1],'-r');
p(3).HandleVisibility='off';

xlabel('Time (s)')
ylabel('Rate of change in pressure (kPa/s)')
legend('dpdt_{ro}','+-97th p-tile |dpdt_{ro}|','+-99th p-tile |dpdt_{ro}|','target limit')

linkaxes(ax,'x');

sgtitle('Pressures')

%% Flow rates
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(3,1,1);
plot(out.t,1e3*60*out.q_hwp(:))
hold on
plot(out.t,1e3*60*out.q_hin(:))
plot(out.t,1e3*60*out.q_hout(:))
plot(out.t,1e3*60*out.q_pm(:))
plot(out.t,1e3*60*out.q_hinPRV(:))
plot(out.t,1e3*60*out.q_rv(:))
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{hwp}','q_{hin}','q_{hout}','q_{pm}','q_{hinPRV}','q_{rv}')

ax(2) = subplot(3,1,2);
plot(out.t,1e3*60*out.q_rv(:))
hold on
plot(out.t,1e3*60*out.q_ERUfeed(:))
plot(out.t,1e3*60*out.q_feed(:))
plot(out.t,1e3*60*out.q_perm(:))
plot(out.t,1e3*60*out.q_roPRV(:))
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{rv}','q_{ERUfeed}','q_{feed}','q_{perm}','q_{roPRV}')

ax(3) = subplot(3,1,3);
plot(out.t,1e3*60*out.q_lwp(:))
hold on
plot(out.t,1e3*60*out.q_pm(:))
plot(out.t,1e3*60*out.q_c(:))
plot(out.t,1e3*60*out.q_ERUfeed(:))
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{lwp}','q_{pm}','q_{c}','q_{ERUfeed}')

linkaxes(ax,'x')

sgtitle('Flow rates')


%% Controller behavior
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 9;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(5,1,1);
plot(out.t,1e-6*out.p_hout)
hold on
plot(out.t,1e-6*out.p_ro)
plot(out.t,1e-6*out.control.p_filt)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{hout}','p_{ro}','p_{filt}')

ax(2) = subplot(5,1,2);
hold on
yyaxis left
plot(out.t,60/(2*pi)*out.control.w_pm)
plot(out.t,60/(2*pi)*out.w_pm)
ylabel('Shaft speed (rpm)')
yyaxis right
plot(out.t,1e-3*out.Tgen)
ylabel('Torque (kNm)')
legend('nominal','actual','Generator')

ax(3) = subplot(5,1,3);
hold on
plot(out.t,1e3*60*out.q_pm)
plot(out.t,1e3*60*out.q_hout)
plot(out.t,1e3*60*out.q_rv)
ylabel('Flow rate (Lpm)')
legend('q_{pm}','q_{hout}','q_{rv}')

ax(4) = subplot(5,1,4);
hold on
plot(out.t,out.control.errInt_p_filt)
legend('pressure control')
ylabel('Error integral')
xlabel('Time (s)')

ax(5) = subplot(5,1,5);
hold on
plot(out.t,sqrt(1e3)*1e3*out.control.kv_ideal)
plot(out.t,sqrt(1e3)*1e3*out.control.kv_rv)
legend('kv_ideal','kv_rv')
ylabel('valve coefficient (L/s/kPa^{1/2})')
xlabel('Time (s)')

linkaxes(ax,'x')

sgtitle('Controller behaviour')

%% WEC-driven pump
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 6;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(2,1,1);
plot(out.t,1e-6*out.p_a)
hold on
plot(out.t,1e-6*out.p_b)
plot(out.t,1e-6*out.p_hin,'k-')
plot(out.t,1e-6*out.p_lout,'k--')
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{a}','p_{b}','p_{hin}','p_{lout}')

ax(2) = subplot(2,1,2);
plot(out.t,1e3*60*out.q_hwp,'k-')
hold on
plot(out.t,1e3*60*out.q_lwp,'k--')
plot(out.t,1e3*60*out.q_sv,'b')
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{hin}','q_{lout}','q_{sv}')

linkaxes(ax,'x');

sgtitle('WEC-driven Pump Dynamics')

%% Pipelines
bottomEdge = 1;
leftEdge = 3;
width = 7.5; % one column: 3+9/16, two column: 7.5
height = 8;
fontSize = 8;
lineWidth = 1;

fig = figure;
fig.Units = 'inches';
fig.Position = [leftEdge bottomEdge width height ];

ax(1) = subplot(4,1,1);
plot(out.t,1e-6*out.p_hin)
hold on
plot(out.t,1e-6*out.p_hout)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{hin}','p_{hout}')

ax(2) = subplot(4,1,2);
plot(out.t,1e3*60*out.q_hin,'k-')
hold on
plot(out.t,1e3*60*out.q_hout,'k--')
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{hin}','q_{hout}')

ax(3) = subplot(4,1,3);
plot(out.t,1e-6*out.p_lin)
hold on
plot(out.t,1e-6*out.p_lout)
xlabel('Time (s)')
ylabel('Pressure (MPa)')
legend('p_{lin}','p_{lout}')

ax(4) = subplot(4,1,4);
plot(out.t,1e3*60*out.q_lin,'k-')
hold on
plot(out.t,1e3*60*out.q_lout,'k--')
xlabel('Time (s)')
ylabel('Flow rate (Lpm)')
legend('q_{lin}','q_{lout}')

linkaxes(ax,'x');

sgtitle('Pipeline Dynamics')
