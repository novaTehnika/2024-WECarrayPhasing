function par = parameters_timeAve_hydElecPTO()
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_timeAve_hydElecPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/3/2024
%
% PURPOSE:
% This function outputs parameters for a model of a wave energy-to-electric
% power system.
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
% 7/3/2024 - adapted from parameters_timeAvePTO.m in 
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

    % fluid properties
    par.rho = 0.860*1000; % [kg/m3] density of working fluid
    par.mu = (32)*1e-6*par.rho; % [(cSt) -> Pa-s]  Dynamic (absolute) viscosity
    par.beta = 1.8e9; % [Pa]  Bulk Modulus of air free fluid

     % Charge pump
    par.p_l = 0.3e6; % [Pa] charge pressure

     % hydraulic motor
    par.motor.D = (500)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]  Motor displacement
    par.motor.w_max = (3600)/60*2*pi; % [(rpm) -> rad/s] maximum speed of motor
    par.motor.w_min = (0)/60*2*pi; % [(rpm) -> rad/s] minimum speed of motor

      % efficiency model coeffs. (McCandlish and Dory model)
       % Axial piston pump (data from Danfoss APP 43/1700)
    par.motor.C_s = 1.529827313040421e-09;
    par.motor.V_r = 1.103;
    par.motor.C_v = 9.272715872908728e+03;
    par.motor.C_f = 0.047594936708861;

     % Generator
    par.eta_gen = 0.9; % [-] elec. generator efficiency

     % WEC-driven pump
    par.eta_w = 0.9; % [-] WEC-driven pump efficiency

end