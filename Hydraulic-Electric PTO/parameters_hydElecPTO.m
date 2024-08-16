function par = parameters_hydElecPTO(par,filenameCoeff,filenameRadSS)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_hydElecPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This function loads default parameters into the "par" structure. This
% includes calling a simular function to load parameters for the WEC model.
% File names for the hydrodynamic data for the WEC are passed through that
% function. Not all parameters can be modified after this function runs
% without issue. Changing parameters in any design study script should be
% check for affecting the other parameters.
%
% This function is for a use with sys_hydElecPTO.m.
%
% FILE DEPENDENCY:
% ../WEC model/
%   parameters_WECmodel.m
%
% UPDATES:
% 6/7/2024 - created from parameters_parPTO.m.
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

    % fluid and entrianed gas properties
    par.rho = 0.860*1000; % [kg/m3] density of working fluid
    par.mu = (32)*1e-6*par.rho; % [(cSt) -> Pa-s]  Dynamic (absolute) viscosity
    par.beta = 1.8e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    % WEC parameters
    par = parameters_WECmodel(par,filenameCoeff,filenameRadSS);

    % WEC-pump
    par.theta_max = pi/2; % [rad] maximum stroke (plus and minus) from upright poistion
     % pumping chamber
    par.D_WEC = 0.15;         % [m^3/rad] flap pump displacement
    V_wecPumpTotal = 2*par.theta_max*par.D_WEC;
    par.V_wecDead = V_wecPumpTotal/2*(0.2); % [m^3] dead volume attached to each port (a and b) of the WEC-driven pump
    par.eta_v_WEC = 1;
    par.eta_m_WEC = 0.9;                % Flap pump mechanical efficiency

     % check valves
      % inlet check valves (low-pressure)
    par.kvWECin = 1.5^2*2.1082e-04;
    par.pc_WECin = 0.5e5; % [Pa] Cracking pressure
    par.dp_WECin = 0.5e5; % [Pa] margin between cracking pressure and fully open condition
    
     % outlet check valves (high-pressure)
    par.kvWECout = 2.1082e-04;
    par.pc_WECout = 0.5e5; % [Pa] Cracking pressure
    par.dp_WECout = 0.5e5; % [Pa] margin between cracking pressure and fully open condition

    % power control unit
      % pump/motor
    par.motor.D = (1000)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement
    par.motor.w_max = (3600)/60*2*pi; % [(rpm) -> rad/s] maximum speed of motor
    par.motor.w_min = (1)/60*2*pi; % [(rpm) -> rad/s] minimum speed of motor

      % efficiency model coeffs. (McCandlish and Dory model)
       % Axial piston pump (data from Danfoss APP 43/1700)
    par.motor.C_s = 1.529827313040421e-09;
    par.motor.V_r = 1.103;
    par.motor.C_v = 9.272715872908728e+03;
    par.motor.C_f = 0.047594936708861;
    
      % generator
    par.eta_g = 0.9;  % efficiency of electric generator
       % rotor inertia estimated by NEMA MG-1 (14.46)
    % n_poles = 3; % [qty.] number of poles of induction motor
    % genPowerRating = 100; % [hp] power rating of induction motor
    % par.Jpm = ( 0.02*2^n_poles*genPowerRating^(1.35-.05*n_poles/2) ) ...
    %             * 0.0421401101; % [lb ft^2 --> kg m^2] rotor inertia
    
    % Accumulators
    par.f = 1e-2; % fraction of dead volume of working fluid compared to 
                  % charge volume
     % LPA
    par.Vc_l = (2000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_l = 0.15e6; % [Pa] charge pressure
     % HPA
    par.Vc_h = (5000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_h = 4e6; % [Pa] charge pressure

    % Contoller Parameters
    par.control.p_l_nom = 0.5e6; % [Pa]
    par.control.p_h_nom = 20e6; % [Pa]

     % Signal filtering
    par.control.tau_pfilt = 0.01; % [s] time constant for LPF for pressure signal
    
     % PI control of pressure using w_m
    par.control.w_m_ctrl.max = par.motor.w_max;
    par.control.w_m_ctrl.min = par.motor.w_min;
    par.control.w_m_ctrl.kp = 5e-4;
    par.control.w_m_ctrl.ki = 0*5e-6;
    
    % Charging system (Intake & Boost pump)
    par.cn = 5.5;
    par.cq = -5e6;
    par.w_c = (2500)*2*pi/60; % [(rpm) -> rad/s]
    par.eta_c = 0.7;  % pumping efficiency of pressure boost pump
    par.eta_m = 0.9;  % efficiency of charge pump motor

    % Pressure relief valves
     % inlet to low-pressure pipeline/outlet of charge pump
    maxPressure = 10e5; % [Pa]
    margin = 1e4; % [Pa]
    maxFlow = (100)*1e-3; % [(L/s) -> m^3/s]
    par.lPRV.p_crack = maxPressure - margin;
    par.lPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

     % high-pressure outlet of WEC-driven pump
    maxPressure = 35e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = (100)*1e-3; % [(L/s) -> m^3/s]
    par.hPRV.p_crack = maxPressure - margin;
    par.hPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;


end