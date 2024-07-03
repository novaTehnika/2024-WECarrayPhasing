function par = parameters_timeAvePTO()
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_timeAvePTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/31/2021
%
% PURPOSE:
% This function outputs parameters for models of a wave-powered reverse 
% osmosis system.
%
% The models are simple, static models with that include two-way coupling 
% with the time-averaged simulation results of a WEC; the coupling is set 
% up such that the reaction force from the PTO is a function of the average
% WEC speed (or power absortion) and the average WEC speed (or power 
% absorption) is function of the reaction torque from the PTO.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 12/31/2021 - created.
%
% Copyright (C) 2022  Jeremy W. Simmons II
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

     % RO module
    par.p_osm =  2.275e6; % [Pa] osmotic pressure of feed water
    par.A_w = 2.57e-12; % [m^3/(N-s)] permeabiity coefficient (Yu and Jenne,2018)
    par.Y = 0.25;

     % Charge pump
    par.p_c = 0.3e6; % [Pa] charge pressure
    par.eta_c = 0.7; % [-] charge pump efficiency
    par.eta_m = 0.9; % [-] elec. motor efficiency

     % Pump/motor
    par.eta_pm = 0.9; % [-] pump/motor efficiency
    par.eta_gen = 0.9; % [-] elec. generator efficiency

     % WEC-driven pump
    par.eta_w = 0.9; % [-] WEC-driven pump efficiency

end