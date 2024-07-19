% run_timeAvePTOsizing.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5 July 2024
%
% PURPOSE:
% The purpose of this script is to determine optimal operating parameters
% for a hydraulic-electric PTO architectures for wave-powered electric 
% power production operating in a known distibution of sea conditions,
% WEC-driven pump displacement, and hydraulic motor displacement. The 
% design is specified in the script under the section heading "Variables".
% 
% The model is a simple, static model with that includes two-way coupling
% with the time-averaged simulation results of a WEC; the coupling is set
% up such that the reaction force from the PTO is a function of the average
% WEC speed (or power absortion) and the average WEC speed (or power
% absorption) is function of the reaction torque from the PTO.
% 
% The optimization of the operating pressure and switching duty is a
% nonlinear, constrained optimization which seeks to maximize the electric
% power production subject to the following constraints:
% 1) the operating pressure at the RO module is within a prescribe range
% 2) the minimum speed of the hydrualic motor is maintained
% If the design does not meet these constraints, the value of zero is
% recorded for electric power production rate for that design and sea
% condition.
% 
% Once the results are obtained for the specified grid of design
% parameters, design parameters are selected for each sea state to give 
% optimal combinations based on four possible configurations:
% A) the WEC-driven pump and hydraulic motor displacement are fixed across
% all sea conditions.
% B) the WEC-driven pump displacement is variable across all sea states 
% while the hydraulic motor displacement is fixed.
% C) the hydraulic motor displacement is variable across all sea states
% while the WEC-driven pump displacement is fixed.
% D) Both the WEC-driven pump and hydraulic motor displacement are variable
% across all sea states
% The design combinations are organized based on the largest allowable
% value for the two design parameters; the design parameters are selected 
% for each sea state from the set of designs that have design parameters 
% values less than or equal to the maximum allowable value, when variable
% by sea state. Otherwise, when the design parameter is fixed across the
% sea states, the set includes only designs that have the same value of the
% fixed design parameter. The design with greatest electric power
% production is selected from the set of available designs in each sea
% state.
%
%
% FILE DEPENDENCIES:
% PTOsizing_multiSS.m
% model_timeAve_hydElecPTO
% parameters_timeAve_hydElecPTO.m
% studyData_coulombPTO_dampingStudy_24-Aug-2022_1_slim.mat
% studyData_coulombPTO_dampingStudy_6SS_20230724_slim.mat
% loadColors.m
% maxRC.m
% zero2nan.m
%
% UPDATES
% 7/5/2024 - adapted from 2021-TimeAvePTOarchetectureStudy repository.
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
addpath('Hydraulic-Electric PTO, TimeAveDesign')
addpath('Utilities')
[git_hash_string, git_status_string] = get_current_git_hash();

%% %%%%%%%%%%%%   SIMULATION/DESIGN PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize parameter structure and get misc. base parameters
par = parameters_timeAve_hydElecPTO();

% bounds on pressures in system
bounds.p_h_bnds = [4e6 30e6]; % [Pa/Pa] Bounds for system pressure
bounds.D_bnds = [0.1 1]; % [-] bounds for valve switching duty

% WEC: load time averaged results for WEC performance
switch 2
    case 1
        filename_WECpowerCurve = ...
            'studyData_coulombPTO_dampingStudy_20220927_slim.mat';
    case 2
        filename_WECpowerCurve = ...
            'studyData_coulombPTO_dampingStudy_6SS_20230724_slim.mat';
end
load(filename_WECpowerCurve)

par.T_c_data = T_c_data; % [Nm] Torque applied to WEC by PTO
par.PP_w_data = PP_w_data; % [W] Power transmitted by WEC to PTO
par.weight = weight;
par.Hs = Hs;
par.Tp = Tp;

SSset = 1:numel(Hs);
par.SSset = SSset;

clearvars T_c_data PP_w_data weight Hs Tp

%% %%%%%%%%%%%%   Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Available PTO configurations
PTOarray = [1 1 1 1 2 2 2 2];
design_case = [1 2 3 4 1 2 3 4];

% Specified design
iPTO = 1;
 % WEC-driven pump displacement
D_w = 0.1;
f_D_w = 0.01;
inc_D_w = (0.01)*D_w;
 % motor displacement
D_m = (500)*1e-6/(2*pi); % [(cc/rev) -> m^3/rad]
f_D_m = 0.01;
inc_D_m = (0.01)*D_m;

% Create array of displacement and motor displacement if these are variable.
 % WEC-driven pump displacment
if design_case(iPTO) == 2 || design_case(iPTO) == 4
    min_D_w = inc_D_w*(floor(f_D_w*D_w/inc_D_w)+(mod(D_w*f_D_w,inc_D_w)>0));
    D_wArray  = min_D_w:inc_D_w:D_w;
else
    D_wArray = D_w;
end

 % motor displacement
if design_case(iPTO) == 3 || design_case(iPTO) == 4
    min_D_m = inc_D_m*(floor(f_D_m*D_m/inc_D_m)+(mod(D_m*f_D_m,inc_D_m)>0));
    D_mArray  = min_D_m:inc_D_m:D_m;
else
    D_mArray = D_m;
end

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = PTOsizing_multiSS(D_wArray,D_mArray, ...
                bounds,PTOarray(iPTO),design_case(iPTO),par);

%% %%%%%%%%%%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = ['data_PTOsizing_',char(datetime("now",'Format','yyyyMMddHHmmss'))];
save(filename,'-v7.3')
