function [dydt, nonState, control] = sys_parPTO(t,y,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys_parPTO.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a parallel-type wave energy PTO 
% with long high- and low-pressure pipelines and an option for a throttling
% valve for enhanced pressure ripple filtering.
%
% FILE DEPENDENCY:
% ../Parallel-type PTO/
%   stateIndex_refPTO.m
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
% 11/2/2023 - created from sys_refPTO.m.
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

%% Calculate control input and the change in controller states (if any)

%% PI control of p_hout using w_pm
 % Error
err_p = y(par.iy.p_filt) - par.control.p_ro_nom;
 % Feedforward
w_ff = 0*par.control.w_pm_ctrl.min;
 % Control signal
w_pm_nom = w_ff ...
    + (par.control.w_pm_ctrl.kp*err_p ...
    + par.control.w_pm_ctrl.ki*y(par.iy.errInt_p_filt));
nomAboveMax = w_pm_nom > par.control.w_pm_ctrl.max;
nomBelowMin = w_pm_nom < par.control.w_pm_ctrl.min;
p_hAbovep_ro = par.rvConfig.included*(y(par.iy.p_hout) > y(par.iy.p_ro)) ...
            + ~par.rvConfig.included;
control.w_pm = p_hAbovep_ro*(w_pm_nom ...
          + nomAboveMax*(par.control.w_pm_ctrl.max - w_pm_nom) ...
          + nomBelowMin*(par.control.w_pm_ctrl.min - w_pm_nom)) ...
          + ~p_hAbovep_ro*par.control.w_pm_ctrl.min;

%% Feedforward control of active RO inlet valve
% determine ideal valve coefficient to satisfy dpdt limit
C_ro = capAccum(y(par.iy.p_ro),par.pc_ro,par.Vc_ro,par.f,par);
q_perm = par.Sro*par.Aperm*(y(par.iy.p_ro) - par.p_perm - par.p_osm);
q_ro = q_perm/par.Y ...
            *(1-par.ERUconfig.present*par.eta_ERUv^2*(1-par.Y));
dp = y(par.iy.p_hout) - y(par.iy.p_ro);
control.kv_ideal = (sign(dp)*C_ro*par.control.dpdt_ROmax + q_ro) ...
        /sqrt(abs(dp));
% satisfy conditions related to limitations, control objectives, 
% and configuration
p_roAbovepc_ro = y(par.iy.p_ro) > par.pc_ro;
rvOpen = ~p_roAbovepc_ro | ~par.rvConfig.active; % cond.s for valve full open
control.kv_rv = ~rvOpen*max(0,min(par.kv_rv,control.kv_ideal)) ...
               + rvOpen*par.kv_rv;
              
%% deriviatives for filtered signal and error integrals (w/
 % anti-wind-up)
dydt_control = [    % filtered signal
                (y(par.iy.p_ro) - y(par.iy.p_filt))...
                /par.control.tau_pfilt;

                    % error integral for pressure control
                (y(par.iy.errInt_p_filt) < 1e12 ...
                & y(par.iy.errInt_p_filt) > -1e12) ...
                *err_p];

%% Calculate non-state variables like coeffs., forces, and flow rates
% Accumulator capacitance
nonState.C_lin = capAccum(y(par.iy.p_lin),par.pc_lin,par.Vc_lin,par.f,par);
nonState.C_lout = capAccum(y(par.iy.p_lout),par.pc_lout,par.Vc_lout,par.f,par);
nonState.C_hin = capAccum(y(par.iy.p_hin),par.pc_hin,par.Vc_hin,par.f,par);
nonState.C_hout = capAccum(y(par.iy.p_hout),par.pc_hout,par.Vc_hout,par.f,par);
nonState.C_ro = capAccum(y(par.iy.p_ro),par.pc_ro,par.Vc_ro,par.f,par);

% WEC-driven pump
 % pumping chamber capacitance
V_a = par.V_wecDead + par.D_WEC*(par.theta_max - y(par.iy.theta));
nonState.C_a = deadVCap(y(par.iy.p_a),V_a,par);
V_b = par.V_wecDead + par.D_WEC*(par.theta_max + y(par.iy.theta));
nonState.C_b = deadVCap(y(par.iy.p_b),V_b,par);

 % Switching valve flow
dp = y(par.iy.p_a) - y(par.iy.p_b);
lin_lowdp = (abs(dp) < par.dp_svlin);
nonState.q_sv = areaFracPWM(t,par.duty_sv,par.T_sv,par.tr_sv) ...
                *(lin_lowdp*par.kv_svlin*dp ...
                + ~lin_lowdp*sign(dp)*par.kv_sv*sqrt(abs(dp)));
 % check valves
nonState.q_ain = flowCV(y(par.iy.p_lout) - y(par.iy.p_a), ...
                 par.kvWECin,par.pc_WECin,par.dp_WECin);
nonState.q_aout = flowCV(y(par.iy.p_a) - y(par.iy.p_hin), ...
                 par.kvWECout,par.pc_WECout,par.dp_WECout);
nonState.q_bin = flowCV(y(par.iy.p_lout) - y(par.iy.p_b), ...
                 par.kvWECin,par.pc_WECin,par.dp_WECin);
nonState.q_bout = flowCV(y(par.iy.p_b) - y(par.iy.p_hin), ...
                 par.kvWECout,par.pc_WECout,par.dp_WECout);

 % Reaction torque on WEC
delta_p_wp = y(par.iy.p_b)-y(par.iy.p_a);
WECpumpPumping = y(par.iy.theta_dot)*delta_p_wp < 0;
WECpumpMotoring = ~WECpumpPumping;

nonState.T_pto = par.D_WEC*delta_p_wp...
            * (WECpumpPumping/par.eta_m_WEC ...
            + WECpumpMotoring*par.eta_m_WEC);

% House power pump/motor
delta_p_pm = y(par.iy.p_lin) - y(par.iy.p_hout);
nonState.pmPumping = control.w_pm*(delta_p_pm) >= 0;
nonState.pmMotoring = control.w_pm*(delta_p_pm) < 0;

volLoss_pm = par.APP.C_s*(delta_p_pm/(par.mu*abs(control.w_pm))) ...
                   + (delta_p_pm/par.beta)*(par.APP.V_r + 1);
nonState.eta_v_pm = nonState.pmPumping*(1 - volLoss_pm) ...
                  + nonState.pmMotoring/(1 + volLoss_pm);
nonState.q_pm = (nonState.pmPumping*nonState.eta_v_pm ...
                + nonState.pmMotoring/nonState.eta_v_pm) ...
                *par.D_pm*control.w_pm;

mechLoss_pm = par.APP.C_v*par.mu*abs(control.w_pm)/delta_p_pm + par.APP.C_f;
nonState.eta_m_pm = nonState.pmPumping/(1 + mechLoss_pm) ...
                  + nonState.pmMotoring*(1 - mechLoss_pm);
nonState.Tpm = (nonState.pmPumping/nonState.eta_m_pm ...
                + nonState.pmMotoring*nonState.eta_m_pm)...
                *par.D_pm*delta_p_pm;

% Reverse osmosis module
nonState.q_perm = par.Sro*par.Aperm*(y(par.iy.p_ro) - par.p_perm - par.p_osm);
nonState.q_feed = nonState.q_perm/par.Y;

 % RO inlet valve/"ripple filter"
dp = y(par.iy.p_hout) - y(par.iy.p_ro);
nonState.q_rv = par.rvConfig.included ...
               * control.kv_rv*sqrt(abs(dp))*sign(dp);

 % ERU
nonState.q_brine = nonState.q_feed - nonState.q_perm;
nonState.q_ERUfeed = par.ERUconfig.present ...
                     *(par.eta_ERUv)^2*nonState.q_brine;

% Charge Pump
dP_SO = (y(par.iy.p_lin) - par.p_o) - par.cn*par.w_c^2; % difference between shut-off pressure and current pressure differential
nonState.q_c = (dP_SO < 0)* sqrt(dP_SO/par.cq);

% Pressure relief valves
nonState.q_linPRV = flowPRV(y(par.iy.p_lin),par.lPRV.p_crack,par.lPRV.C);
nonState.q_hinPRV = flowPRV(y(par.iy.p_hin),par.hPRV.p_crack,par.hPRV.C);
nonState.q_roPRV = flowPRV(y(par.iy.p_ro),par.roPRV.p_crack,par.roPRV.C);

%% Calculate the hydrodynamics WEC state derivatives and output nonstate
% variables like forces and wave elevation
[dydt_WEC, nonState.torqueWEC, nonState.waveElev] = ...
                    flapModel(t,y(par.iy.WEC),nonState.T_pto,par);

%% Calculate states of the pipelines
if par.plConfig.included
    [q_lin, q_lout, dydt_LPPL, ~] = pipelineNPi(t,y(par.iy.LPPL), ...
                            y(par.iy.p_lin),y(par.iy.p_lout),par,1,0);
    [q_hin, q_hout, dydt_HPPL, ~] = pipelineNPi(t,y(par.iy.HPPL), ...
                            y(par.iy.p_hin),y(par.iy.p_hout),par,2,0);
else
    q_lin = 0;
    q_lout = 0;
    q_hin = 0;
    q_hout = 0;
end

%% State derivatives
dydt = zeros(par.iy.ny,1);

% WEC-driven pump
dydt(par.iy.p_a) = 1/nonState.C_a*(par.D_WEC*y(par.iy.theta_dot) - nonState.q_sv ...
                + nonState.q_ain - nonState.q_aout);
dydt(par.iy.p_b) = 1/nonState.C_b*(-par.D_WEC*y(par.iy.theta_dot) + nonState.q_sv ...
                + nonState.q_bin - nonState.q_bout);

% low-pressure accumulators
dydt_p_lin = 1/(nonState.C_lin ...
                + ~par.plConfig.included*nonState.C_lout) ...
                *(nonState.q_c - q_lin ...
                + nonState.q_pm - nonState.q_ERUfeed ...
                - nonState.q_linPRV + nonState.q_roPRV);
dydt_p_lout = 1/(nonState.C_lout ...
                + ~par.plConfig.included*nonState.C_lin) ...
                *(q_lout ...
                - nonState.q_ain - nonState.q_bin ...
                + nonState.q_hinPRV);

dydt(par.iy.p_lin) = dydt_p_lin + ~par.plConfig.included*dydt_p_lout;
dydt(par.iy.p_lout) = dydt_p_lout + ~par.plConfig.included*dydt_p_lin;

% high-pressure accumulators
dydt_p_hin = 1/(nonState.C_hin ...
                + ~par.plConfig.included*(nonState.C_hout ...
                + ~par.rvConfig.included*nonState.C_ro)) ...
                *(nonState.q_aout + nonState.q_bout ...
                - q_hin ...
                - nonState.q_hinPRV);
dydt_p_hout = 1/(nonState.C_hout ...
                + ~par.plConfig.included*nonState.C_hin ...
                + ~par.rvConfig.included*nonState.C_ro) ...
                *(q_hout - nonState.q_pm - nonState.q_rv);
dydt_p_ro = 1/(nonState.C_ro ...
                + ~par.rvConfig.included*(nonState.C_hout ...
                + ~par.plConfig.included*nonState.C_hin)) ...
                *(nonState.q_rv - nonState.q_feed + nonState.q_ERUfeed ...
                - nonState.q_roPRV);

dydt(par.iy.p_hin) = dydt_p_hin + ~par.plConfig.included*(dydt_p_hout ...
                                + ~par.rvConfig.included*dydt_p_ro);
dydt(par.iy.p_hout) = dydt_p_hout + ~par.rvConfig.included*dydt_p_ro ...
                                + ~par.plConfig.included*dydt_p_hin;
dydt(par.iy.p_ro) = dydt_p_ro + ~par.rvConfig.included*(dydt_p_hout ...
                                + ~par.plConfig.included*dydt_p_hin);

% Control system
dydt(par.iy.control) = dydt_control;

% WEC
dydt(par.iy.theta) = dydt_WEC(1); % angular velocity
dydt(par.iy.theta_dot) = dydt_WEC(2); % angular acceleration
dydt(par.iy.rad) = dydt_WEC(3:end); % radiation damping states for WEC model

% Pipeline
if par.plConfig.included
    dydt(par.iy.LPPL) = dydt_LPPL;
    dydt(par.iy.HPPL) = dydt_HPPL;
end

end
