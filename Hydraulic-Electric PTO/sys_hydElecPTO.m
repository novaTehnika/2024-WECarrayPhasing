function [dydt, nonState, control] = sys_hydElecPTO(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_hydElecPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for an array of wave energy converters
% that share a common hydraulic-electric power take-off.
%
% FILE DEPENDENCY:
% ../Hydraulic-Electric PTO/
%   stateIndex_hydElecPTO.m
% ../WEC model/
%   flapModel.m
% ../Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
% ../Components/Pipeline
%   flowR.m
%   lineCap.m
%   pipelineNPi.m
%
% UPDATES:
% 6/7/2024 - created from sys_parPTO.m.
% 7/3/2024 - set low-pressure branch pressure to fixed
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

%% Calculate control input and the change in controller states (if any)

%% PI control of p_load using w_m
 % Error
err_p = y(par.iy.p_filt) - par.control.p_load_nom;
 % Feedforward
w_ff = par.control.w_m_ctrl.min;
 % Control signal
w_m_nom = w_ff ...
    + (par.control.w_m_ctrl.kp*err_p ...
    + par.control.w_m_ctrl.ki*y(par.iy.errInt_p_filt));
nomAboveMax = w_m_nom > par.control.w_m_ctrl.max;
nomBelowMin = w_m_nom < par.control.w_m_ctrl.min;

control.w_m = w_m_nom ...
          + nomAboveMax*(par.control.w_m_ctrl.max - w_m_nom) ...
          + nomBelowMin*(par.control.w_m_ctrl.min - w_m_nom);
              
%% deriviatives for filtered signal and error integrals (w/
 % anti-wind-up)
dydt_control = [    % filtered signal
                (y(par.iy.p_h) - y(par.iy.p_filt))...
                /par.control.tau_pfilt;

                    % error integral for pressure control
                (y(par.iy.errInt_p_filt) < 1e12 ...
                & y(par.iy.errInt_p_filt) > -1e12) ...
                *err_p];

%% Calculate non-state variables like coeffs., forces, and flow rates
% Accumulator capacitance
nonState.C_l = capAccum(y(par.iy.p_l),par.pc_l,par.Vc_l,par.f,par);
nonState.C_h = capAccum(y(par.iy.p_h),par.pc_h,par.Vc_h,par.f,par);

% WEC-driven pump
nonState.q_totalIn = 0;
nonState.q_totalOut = 0;
for iWEC = 1:par.NumWECs
    % pumping chamber capacitance
    V_a = par.V_wecDead + par.D_WEC*(par.theta_max - y(par.iy.theta(iWEC)));
    nonState.C_a(iWEC) = deadVCap(y(par.iy.p_a(iWEC)),V_a,par);
    V_b = par.V_wecDead + par.D_WEC*(par.theta_max + y(par.iy.theta(iWEC)));
    nonState.C_b(iWEC) = deadVCap(y(par.iy.p_b(iWEC)),V_b,par);

    % check valves
    nonState.q_aIn(iWEC) = flowCV(y(par.iy.p_l) - y(par.iy.p_a(iWEC)), ...
                     par.kvWECin,par.pc_WECin,par.dp_WECin);
    nonState.q_aOut(iWEC) = flowCV(y(par.iy.p_a(iWEC)) - y(par.iy.p_h), ...
                     par.kvWECout,par.pc_WECout,par.dp_WECout);
    nonState.q_bIn(iWEC) = flowCV(y(par.iy.p_l) - y(par.iy.p_b(iWEC)), ...
                     par.kvWECin,par.pc_WECin,par.dp_WECin);
    nonState.q_bOut(iWEC) = flowCV(y(par.iy.p_b(iWEC)) - y(par.iy.p_h), ...
                     par.kvWECout,par.pc_WECout,par.dp_WECout);

    nonState.q_totalIn = nonState.q_totalIn ...
                        + nonState.q_aIn(iWEC) + nonState.q_bIn(iWEC);
    nonState.q_totalOut = nonState.q_totalOut ...
                        + nonState.q_aOut(iWEC) + nonState.q_bOut(iWEC);
    
     % Reaction torque on WEC
    delta_p_wp = y(par.iy.p_b(iWEC))-y(par.iy.p_a(iWEC));
    WECpumpPumping = y(par.iy.theta_dot(iWEC))*delta_p_wp < 0;
    WECpumpMotoring = ~WECpumpPumping;
    
    nonState.T_pto(iWEC) = par.D_WEC*delta_p_wp...
                * (WECpumpPumping/par.eta_m_WEC ...
                + WECpumpMotoring*par.eta_m_WEC);
end

% House power pump/motor
delta_p_m = y(par.iy.p_l) - y(par.iy.p_h);
nonState.pmPumping = control.w_m*(delta_p_m) >= 0;
nonState.pmMotoring = control.w_m*(delta_p_m) < 0;

volLoss_m = par.APP.C_s*(delta_p_m/(par.mu*abs(control.w_m))) ...
                   + (delta_p_m/par.beta)*(par.APP.V_r + 1);
nonState.eta_v_m = nonState.pmPumping*(1 - volLoss_m) ...
                  + nonState.pmMotoring/(1 + volLoss_m);
nonState.q_m = (nonState.pmPumping*nonState.eta_v_m ...
                + nonState.pmMotoring/nonState.eta_v_m) ...
                *par.D_m*control.w_m;

mechLoss_m = par.APP.C_v*par.mu*abs(control.w_m)/delta_p_m + par.APP.C_f;
nonState.eta_m_m = nonState.pmPumping/(1 + mechLoss_m) ...
                  + nonState.pmMotoring*(1 - mechLoss_m);
nonState.T_m = (nonState.pmPumping/nonState.eta_m_m ...
                + nonState.pmMotoring*nonState.eta_m_m)...
                *par.D_m*delta_p_m;

% Charge Pump
dP_SO = (y(par.iy.p_l) - par.p_o) - par.cn*par.w_c^2; % difference between shut-off pressure and current pressure differential
nonState.q_c = (dP_SO < 0)* sqrt(dP_SO/par.cq);

% Pressure relief valves
nonState.q_lPRV = flowPRV(y(par.iy.p_l),par.lPRV.p_crack,par.lPRV.C);
nonState.q_hPRV = flowPRV(y(par.iy.p_h),par.hPRV.p_crack,par.hPRV.C);

%% Calculate the hydrodynamics WEC state derivatives and output nonstate
% variables like forces and wave elevation
for iWEC = 1:par.NumWECs
    [dydt_WEC(iWEC,:), nonState.torqueWEC(iWEC), nonState.waveElev(iWEC)] = ...
                    flapModel(t,y(par.iy.WEC(iWEC,:)),nonState.T_pto(iWEC),par,iWEC);
end

%% State derivatives
dydt = zeros(par.iy.ny,1);

% WEC and WEC-driven pump(s)
for iWEC = 1:par.NumWECs
    dydt(par.iy.p_a(iWEC)) = 1/nonState.C_a(iWEC) ...
                    *(par.D_WEC*y(par.iy.theta_dot(iWEC)) ...
                    + nonState.q_aIn(iWEC) - nonState.q_aOut(iWEC));
    dydt(par.iy.p_b(iWEC)) = 1/nonState.C_b(iWEC) ...
                    *(-par.D_WEC*y(par.iy.theta_dot(iWEC)) ...
                    + nonState.q_bIn(iWEC) - nonState.q_bOut(iWEC));

    dydt(par.iy.theta(iWEC)) = dydt_WEC(iWEC,1); % angular velocity
    dydt(par.iy.theta_dot(iWEC)) = dydt_WEC(iWEC,2); % angular acceleration
    dydt(par.iy.rad(iWEC,:)) = dydt_WEC(iWEC,3:end); % radiation damping states for WEC model
end

% low-pressure accumulator
% dydt(par.iy.p_l) = 1/nonState.C_l*(nonState.q_c + nonState.q_m ...
%                 - nonState.q_totalIn + nonState.q_hPRV);
dydt(par.iy.p_l) = 0;

% high-pressure accumulators
dydt(par.iy.p_h) = 1/nonState.C_h*(nonState.q_totalOut - nonState.q_m ...
                - nonState.q_hPRV);

% Control system
dydt(par.iy.control) = dydt_control;

end
