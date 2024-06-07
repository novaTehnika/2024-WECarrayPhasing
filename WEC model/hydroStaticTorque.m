function T_hydroStatic = hydroStaticTorque(theta, h_wave, par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hydroStaticTorque.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/29/2021
%
% PURPOSE/DESCRIPTION:
% This function calculates the hydrostatic torque on an oscilating
% surge wave converter. The submerged volume of the flap is approximated
% assuming the depth of submergance under the wave does not change 
% significantly due to the angle of the flap and that the flap is a thin 
% plate. The result of this function is the sum of torque due to the weight
% of the flap and the torque due to the weight of the displaced water.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 6/29/2021 - created from buoyTorque.m developed in 2018 
% by the author. The file is rewritten to take a 
% structure of parameters as an input instead of using hard coded 
% variables.
% 8/23/2022 - added position measurement for center of mass
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


%% Determine buoyant torque from submerged volume

% Height of top edge of flap from hinge
h_top = par.WEC.L_flap*cos(theta);

% Length of submerged part of plate (by centerline) centroids 
% from thin plate assumption
L_sub = (h_top > par.WEC.h+h_wave)*(par.WEC.h+h_wave)/cos(theta) + ...
    (h_top <= par.WEC.h+h_wave)*par.WEC.L_flap; 

% Calculate the centroid of the buoyant force, position perpendicular to
% gravity
centroid_x = L_sub/2*sin(theta);

% calculate bouyant torque
T_buoy = par.WEC.T*L_sub*par.WEC.W*par.WEC.rho*par.WEC.g*centroid_x;

%% Determine torque due to off center weight
% Calculate position of the center of mass of the flap, perpendicular to
% gravity
cm_x = par.WEC.L_cm*sin(theta);

% calculate torque due to weight
T_weight = par.WEC.m*par.WEC.g*cm_x;

%% Determine total effective torque
T_hydroStatic = T_buoy - T_weight;

end