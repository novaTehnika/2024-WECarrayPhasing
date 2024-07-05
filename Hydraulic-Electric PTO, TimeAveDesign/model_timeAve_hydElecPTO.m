function  varargout = model_timeAve_hydElecPTO(x,par,iPTO,outputConfig)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_timeAve_hydElecPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/3/2024
%
% PURPOSE:
% This function implements a model for a wave energy-to-electric power
% system.
% 
% The model is a simple, static model with that includes two-way coupling
% with the time-averaged simulation results of a WEC; the coupling is set
% up such that the reaction force from the PTO is a function of the average
% WEC speed (or power absortion) and the average WEC speed (or power
% absorption) is function of the reaction torque from the PTO.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 7/3/2024 - adapted from model_timeAvePTO.m in 
% 2021-TimeAvePTOarchetectureStudy repository.
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

    power = @(T_c) interp1(par.T_c_data(par.SS,:),...
                par.PP_w_data(par.SS,:),T_c,'spline');

    switch iPTO
        case 2
            duty = x(2);
        case 1
            duty = 1;
    end

    p_h = x(1);
    dp_w = p_h - par.p_l;
    T_c = (duty*par.D_w)*dp_w/par.eta_w;
    PP_w = power(T_c);
    q_w = PP_w*par.eta_w/dp_w;

    % House power pump/motor
     % determine speed
    w_m = (q_w/par.motor.D - par.motor.C_s*(dp_w/par.mu))/ ...
          (dp_w/par.beta*(par.motor.V_r + 1));
    w_m = min(w_m,par.motor.w_max);

     % flow
    q_m = par.motor.D*w_m*(1 + par.motor.C_s*dp_w/(par.mu*w_m) ...
                       + dp_w/par.beta*(par.motor.V_r + 1));

     % torque
    T_m = par.D_m*dp_w*(1 - par.motor.C_v*par.mu*w_m/dp_w + par.motor.C_f);

    % Elec. power
    PP_gen = par.eta_gen*T_m*w_m;

    % Power lost to PRV
    PP_PRV = dp_w*(q_w-q_m);

    switch outputConfig
        case 1
            varargout = {PP_gen};
        case 2
            % constraints c <= 0, ceq = 0
            c(1) = PP_c - PP_gen;
            c(2) = p_h - par.p_h_bnds(2);
            c(3) = par.p_h_bnds(1) - p_h;
            c(4) = T_c - par.T_c_data(par.SS,end);
            c(5) = par.motor.w_min - w_m;
            ceq = [];
            varargout = {c,ceq};
        case 3
            varargout = {T_c, PP_w, PP_gen, PP_PRV};
    end

end