% run_PTOsizing.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% July 24, 2023
%
% PURPOSE:
% The purpose of this script is to perform a design for a wave-powered 
% reverse osmosis system operating in a known distibution of sea conditions
% given a PTO architecture, WEC-driven pump displacement, and install RO 
% membrane area. The design is specified in the script under the section
% heading "Variables".
% 
% The models used in the design are simple, static models with that 
% include two-way coupling with the time-averaged simulation results of a
% WEC; the coupling is set up such that the reaction force from the PTO is   
% a function of the average WEC speed (or power absortion) and the average 
% WEC speed (or power absorption) is function of the reaction torque from 
% the PTO.
%
% In the design, the WEC-driven pump displacement and RO module
% membrane area are varied across a grid of values. For each set
% of these design variables, an optimization is performed to select the
% nominal operating pressure of the system (either the RO feed pressure or 
% the pressure at the outlet of the WEC_driven pump depending on the PTO 
% architecture) and the switching duty (if applicable). This routine is
% performed for each sea state.
% 
% The optimization of the operating pressure and switching duty is a
% nonlinear, constrained optimization which seeks to maximize the permeate
% production subject to the following constraints:
% 1) the electrical power production meets or exceeds the electrical demand
% 2) the operating pressure at the RO module is within a prescribe range
% If the design does not meet these constraints, the value of zero is
% recorded as the permeate production rate for that design.
% 
% Once the results are obtained for the specified grid of design
% parameters, design parameters are selected for each sea state to give 
% optimal combinations based on four possible configurations:
% A) the WEC-driven pump and RO module size are fixed across all sea 
% conditions
% B) the WEC-driven pump is variable across all sea states while the active 
% RO module size is fixed
% C) the active RO module size is variable across all sea states while the
% WEC-driven pump displacement is fixed.
% D) Both the WEC-driven pump and active RO module size are variable across
% all sea states
% The design combinations are organized based on the largest allowable
% value for the two design parameters; the design parameters are selected 
% for each sea state from the set of designs that have design parameters 
% values are less than or equal to the maximum allowable value, when 
% variable by sea state. Otherwise, when the design parameter is fixed
% across the sea states, the set includes only designs that have the same 
% value of the fixed design parameter. The design with greatest permeate
% production is selected from the set of available designs in each sea
% state.
%
%
% FILE DEPENDENCIES:
% PTOsizing_multiSS.m
% zero2nan.m
% parameters_timeAvePTO.m
% loadColors.m
% data_coulombPTO_dampingStudy_24-Aug-2022_1_slim.mat
% data_coulombPTO_dampingStudy_6SS_20230724_slim.mat
%
% UPDATES
% 07/24/2023 - created from study_PTOsizing.m.
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
clear
clc

%% %%%%%%%%%%%%   SIMULATION/DESIGN PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize parameter structure and get misc. base parameters
par = parameters_timeAvePTO();

% bounds on pressures in system
bounds.p_f_bnds = [4e6 8e6]; % [Pa/Pa] Bounds for feed pressure
bounds.p_w_bnds = [4e6 30e6]; % [Pa/Pa] Bounds for pressure at WEC driven pump
bounds.D_bnds = [0.1 1]; % [-] bounds for valve switching duty

% WEC: load time averaged results for WEC performance
switch 1
    case 1
        filename_WECpowerCurve = 'data_coulombPTO_dampingStudy_20220927_slim.mat';
        SSset = 1:114;
    case 2
        filename_WECpowerCurve = 'data_coulombPTO_dampingStudy_6SS_20230724_slim.mat';
        SSset = 1:6;
end
par.SSset = SSset;
load(filename_WECpowerCurve)
par.T_c_data = T_c_data; % [Nm] Torque applied to WEC by PTO
par.PP_w_data = PP_w_data; % [W] Power transmitted by WEC to PTO
par.weight = weight;
par.Hs = Hs;
par.Tp = Tp;
clearvars T_c_data PP_w_data weight Hs Tp

%% %%%%%%%%%%%%   Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Available PTO configurations
PTOarray = [1 1 3 3 4 1 1 3 3 4];
design_case = [1 2 1 2 1 3 4 3 4 3];

% Specified design
iPTO = 1;
 % WEC-driven pump displacement
D_w = 0.23;
f_D_w = 0.01;
inc_D_w = (0.01)*D_w;
 % Total Ro membrane area;
S_ro = 3700;
f_S_ro = 0.01;
inc_S_ro = 37;
 % ERU configuration
ERUconfig = 1; % 0-w/o ERU; 1-w/ ERU

% Create array of displacement and memebrane area if these are variable.
 % WEC-driven pump displacment
if design_case(iPTO) == 2 || design_case(iPTO) == 4
    min_D_w = inc_D_w*(floor(f_D_w*D_w/inc_D_w)+(mod(D_w*f_D_w,inc_D_w)>0));
    D_wArray  = min_D_w:inc_D_w:D_w; % [m^3/s] displacement
else
    D_wArray = D_w;
end

 % membrane area in RO module
if design_case(iPTO) == 3 || design_case(iPTO) == 4
    min_S_ro = inc_S_ro*(floor(f_S_ro*S_ro/inc_S_ro)+(mod(S_ro*f_S_ro,inc_S_ro)>0));
    S_roArray  = min_S_ro:inc_S_ro:S_ro; % [m^3/s] displacement
else
    S_roArray = S_ro;
end

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = PTOsizing_multiSS(D_wArray,S_roArray, ...
                bounds,PTOarray(iPTO),design_case(iPTO), ...
                ERUconfig,par);

%% %%%%%%%%%%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['data_PTOsizing_',char(datetime("now",'Format','yyyyMMddHHmmss'))];
save(filename,'-v7.3')
