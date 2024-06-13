function y0 = initialConditionDefault_hydElecPTO(par)
% initialConditionDefault_hydElecPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This script loads default intial conditions for the hydElecPTO model.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 6/7/2024 - Created from initialConditionDefault_parPTO.m.
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

for iWEC = 1:par.NumWECs
    y0(par.iy.theta(iWEC),1) = 0.1*pi/2; % [rad]
    y0(par.iy.theta_dot(iWEC),1) = 0; % [rad/s]
    y0(par.iy.rad(iWEC,:),1) = zeros(numel(par.iy.rad(iWEC,:)),1);
    
    y0(par.iy.p_a(iWEC),1) = 1e6; % [Pa]
    y0(par.iy.p_b(iWEC),1) = 1e6; % [Pa]
end

y0(par.iy.p_l,1) = par.control.p_l_nom; % [Pa]
y0(par.iy.p_h,1) = par.control.p_h_nom; % [Pa]

y0(par.iy.p_filt,1) = par.control.p_h_nom;
y0(par.iy.errInt_p_filt,1) = 0;

end