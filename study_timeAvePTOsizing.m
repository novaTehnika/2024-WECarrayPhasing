% study_timeAvePTOsizing.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5 July 2024
%
% PURPOSE:
% The purpose of this script is to perform a design study sizing a 
% hydraulic-electric PTO architectures for wave-powered electric power
% production operating in a known distibution of sea conditions.
% 
% The model is a simple, static model with that includes two-way coupling
% with the time-averaged simulation results of a WEC; the coupling is set
% up such that the reaction force from the PTO is a function of the average
% WEC speed (or power absortion) and the average WEC speed (or power
% absorption) is function of the reaction torque from the PTO.
%
% In the design study, the WEC-driven pump displacement and hydraulic motor
% displacement are varied using a grid of values. For each set of these
% design variables, an optimization is performed to select the nominal 
% operating pressure of the system and the switching duty (if applicable).
% This routine is performed for each sea state.
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
switch 1
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

%% %%%%%%%%%%%%   Study Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% WEC-driven pump displacment
nD_w = 201; % Size of array for pump displacement
D_wArray = logspace(log10(0.01),log10(0.5),nD_w); % [m^3/rad] displacement

% motor displacement
nD_m = 201; % Size of array for motor displacement
D_mArray = logspace(log10(1e2),log10(5e3),nD_m)* ...
                                    1e-6/(2*pi); % [(cc/rev) -> m^3/rad]

% Specify PTO configurations
 % Available PTO configurations
PTOarray = [1 1 1 1 2 2 2 2];
design_case = [1 2 3 4 1 2 3 4];

iPTO = 1;

%% %%%%%%%%%%%%   COLLECT DATA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data(iPTO) = PTOsizing_multiSS(D_wArray,D_mArray, ...
             bounds,PTOarray(iPTO),design_case(iPTO),par);

%% %%%%%%%%%%%%   SAVE DATA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = ['data_PTOsizing_contour_',char(datetime("now",'Format','yyyyMMdd'))];
save(filename,'-v7.3')

return

%% %%%%%%%%%%%%   PLOTTING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lineWidth = 1.5;

loadColors
levels = [50 100 125 150 165];
color1 = maroon;
color2 = gold; 
nlines = length(levels);
colorWeight = linspace(0,1,nlines);
co = zeros(nlines,3);
for i=1:nlines
    co(i,:) = colorWeight(i)*color2 + colorWeight(nlines-(i-1))*color1; 
end

f = figure;
colormap(f,co)

[M,c1] = contour(data.D_m*2*pi/1e-6,data.D_w,1e-3*data.PP_genTotal, ...
                 levels,'-','ShowText','on')
c1.LineWidth = lineWidth;
xlabel('Motor displacement (cc/rev)')
ylabel('Pump displacement (m^3/rad)')
title(['Yearly Average Elec. Power Production (kW): PTO ',num2str(iPTO)])