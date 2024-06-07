function par = parameters_parPTO(par,filenameCoeff,filenameRadSS)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_parPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% This function loads default parameters into the "par" structure. This
% includes calling a simular function to load parameters for the WEC model.
% File names for the hydrodynamic data for the WEC are passed through that
% function. Not all parameters can be modified after this function runs
% without issue. Changing parameters in any design study script should be
% check for affecting the other parameters.
%
% This function is for a use with sys_parPTO.m.
%
% FILE DEPENDENCY:
% ../WEC model/
%   parameters_WECmodel.m
%
% UPDATES:
% 11/2/2023 - created from parameters_refPTO.m.
%
% Copyright (C) 2023  Jeremy W. Simmons II
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
    par.rho = 1023; % [kg/m3] density of air
    par.mu = 9.4e-4; % [Pa-s]  Dynamic (absolute) viscosity
    par.beta = 2.2e9; % [Pa]  Bulk Modulus of air free fluid
    par.p_vap = 0.037e5; % [Pa] vapour pressure of seawater
    par.R = 0.0001; % [-] fraction  Baseline fraction of air by volume entrained in the hydraulic fluid at atm
    par.p_o = 101.3e3; % [Pa]  Atmospheric pressure (reference)
    par.gamma = 1.4; % [-]ratio of specific heats for air
    
    % WEC parameters
    par = parameters_WECmodel(par,filenameCoeff,filenameRadSS);

    % WEC-pump
    par.theta_max = pi/2; % [rad] maximum stroke (plus and minus) from upright poistion
     % pumping chamber
    par.D_WEC = 0.23;         % [m^3/rad] flap pump displacement
    V_wecPumpTotal = 2*par.theta_max*par.D_WEC;
    par.V_wecDead = V_wecPumpTotal/2*(0.2); % [m^3] dead volume attached to each port (a and b) of the WEC-driven pump
    par.eta_v_WEC = 1;
    par.eta_m_WEC = 0.9;                % Flap pump mechanical efficiency

     % switching valve
    par.duty_sv = 0;
    par.T_sv = 0.5; % [s] switching period (1/switching freq)
    par.tr_sv = 0.05; % transition ratio (frac. of cyc. transitioning 1way)

    dp_rated = 1e6; % [Pa] 
    q_rated = 300e-3; % [m^3/s]
    par.kv_sv = q_rated/sqrt(dp_rated);
    par.dp_svlin = 1e4; % [Pa] linearization margin
    par.kv_svlin = par.kv_sv*sqrt(par.dp_svlin)/par.dp_svlin;

     % check valves
      % inlet check valves (low-pressure)
    par.kvWECin = 1.5*0.32e-3;
    par.pc_WECin = 1e5; % [Pa] Cracking pressure
    par.dp_WECin = 1e5; % [Pa] margin between cracking pressure and fully open condition
    
     % outlet check valves (high-pressure)
    par.kvWECout = 0.32e-3;
    par.pc_WECout = 1e5; % [Pa] Cracking pressure
    par.dp_WECout = 1e5; % [Pa] margin between cracking pressure and fully open condition

    % RO module parameters
    par.p_perm = par.p_o;
    par.p_osm = 2.7e6;
    par.Sro = 3700; % [m^3]
    par.Aperm = 2.57e-12; % [m^3/(N-s)] permeabiity coefficient (Yu and Jenne,2018)
    par.Y = 0.25;

    % ERU
    par.ERUconfig.present = 1;
    par.eta_ERUv = 0.95;
    par.eta_ERUm = 0.95;
    
    % power control unit
      % pump/motor
    par.D_pm = (1000)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement
    par.w_pm_max = (1750)/60*2*pi; % [(rpm) -> rad/s] maximum speed of motor
    par.w_pm_min = (1)/60*2*pi; % [(rpm) -> rad/s] minimum speed of motor

      % efficiency model coeffs. (McCandlish and Dory model)
       % Axial piston pump (data from Danfoss APP 43/1700)
    par.APP.C_s = 3.0554e-10;
    par.APP.V_r = 1.103;
    par.APP.C_v = 7.1755e5;
    par.APP.C_f = 0.0259;
    
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
     % LPA, onshore
    par.Vc_lin = (2000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_lin = 0.15e6; % [Pa] charge pressure
     % LPA, offshore/at inelt of WEC-driven pump
    par.Vc_lout = (8000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_lout = 0.15e6; % [Pa] charge pressure
     % HPA at outlet of WEC-driven pump and inlet of high-pressure pipeline
     % (offshore)
    par.Vc_hin = (5000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_hin = 4e6; % [Pa] charge pressure
     % HPA at outlet of high-pressure pipeline (onshore)
    par.Vc_hout = (5000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_hout = 4e6; % [Pa] charge pressure
     % HPA at inlet to RO module
    par.Vc_ro = (5000)*1e-3; % [(L) -> m^3] gas volume at charge pressure
    par.pc_ro = 4e6; % [Pa] charge pressure

    % Pipelines
    par.plConfig.included = 1; % pipelines are 1 - present, 0 - absent
     % LP pipeline
    LineID = 1;
    par.L_line(LineID) = 500; % [m] length of LP pipeline
    par.d_line(LineID) = 0.2; % [m] diameter of LP pipeline
    par.A_line(LineID) = pi/4*par.d_line(LineID)^2; % crosssectional flow area
    par.n_seg(LineID) = 3; % minimum of 2
    par.I(LineID) = par.rho*(par.L_line(LineID)/par.n_seg(LineID))/par.A_line(LineID);
     % HP pipeline
    LineID = 2;
    par.L_line(LineID) = 500; % [m] length of HP pipeline
    par.d_line(LineID) = 0.1; % [m] diameter of HP pipeline
    par.A_line(LineID) = pi/4*par.d_line(2)^2; % crosssectional flow area
    par.n_seg(LineID) = 3; % minimum of 2
    par.I(LineID) = par.rho*(par.L_line(LineID)/par.n_seg(LineID))/par.A_line(LineID);

    % Contoller Parameters
    par.control.p_h_nom = 6e6; % [Pa]
    par.control.p_ro_nom = par.control.p_h_nom; % [Pa] (not actually used as control ref.)
    par.control.p_l_nom = 0.5e6; % [Pa] (not actually used as control ref.)
    par.control.p_ro_max = 8.3e6; % [Pa]
    par.control.p_ro_min = max(3e6,par.pc_ro); % [Pa]

     % Signal filtering
    par.control.tau_pfilt = 0.01; % [s] time constant for LPF for pressure signal
    
     % PI control of pressure (p_h or p_ro) using w_pm
    par.control.w_pm_ctrl.max = par.w_pm_max;
    par.control.w_pm_ctrl.min = par.w_pm_min;
    par.control.w_pm_ctrl.kp = 5e-4;
    par.control.w_pm_ctrl.ki = 0*5e-6;

    % RO inlet valve for pressure ripple reduction
    par.rvConfig.included = 1; % RO inlet valve is 1 - present, 0 - absent
    par.rvConfig.active = (0)*par.rvConfig.included; % RO inlet valve is 1 - active, 0 - passive
    % dp_rated = 1e5; % [Pa] 
    % q_rated = 50e-3; % [m^3/s]
    % par.kv_rv = q_rated/sqrt(dp_rated);
    par.kv_rv = (4.7)/1e3/sqrt(1e3); % [(L/s/kPa^0.5) -> m^3/s/Pa^0.5]
    par.control.dpdt_ROmax = (10)*6894.76;
    
    % Charging system (Intake & Boost pump)
    par.cn = 5.5;
    par.cq = -5e6;
    par.w_c = (2500)*2*pi/60; % [(rpm) -> rad/s]
    par.eta_c = 0.7;  % pumping efficiency of pressure boost pump
    par.eta_m = 0.9;  % efficiency of charge pump motor
    % par.p_c = .65e6;

    % Pressure relief valves
     % inlet to low-pressure pipeline/outlet of charge pump
    maxPressure = 10e5; % [Pa]
    margin = 1e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.lPRV.p_crack = maxPressure - margin;
    par.lPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;

     % high-pressure outlet of WEC-driven pump
    maxPressure = 20e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = (100)*1e-3; % [(L/s) -> m^3/s]
    par.hPRV.p_crack = maxPressure - margin;
    par.hPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;
    
     % high-pressure inlet to RO module
    maxPressure = 8.3e6; % [Pa]
    margin = 5e4; % [Pa]
    maxFlow = 100e-3; % [m^3/s]
    par.roPRV.p_crack = maxPressure - margin;
    par.roPRV.C = (maxPressure^(3/2) ...
                 - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;


end