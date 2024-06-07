function iy = stateIndex_hydElecPTO(par)
% stateIndex_hydElecPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This script loads the state indices for the hydElecPTO model. The system
% includes shared low and high-pressure accumulators, a control system, and
% multiple WECs with WEC drive pumps. The WEC driven pumps have
% compressible fluid in the pumping chambers.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 6/7/2024 - Created from stateIndex_parPTO.m.
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


iy.p_l = 1;
iy.p_h = 2;

iy.p_filt = 3;
iy.errInt_p_filt = 4;
iy.control = [iy.p_filt; iy.errInt_p_filt];

% first WEC
iy.p_a = 5;
iy.p_b = 6;
iy.theta = 7;
iy.theta_dot = 8;
iy.rad = (1:par.WEC.ny_rad) + iy.theta_dot; % state vector indices for radiation damping states for WEC model
iy.WEC = [iy.theta, iy.theta_dot, iy.rad];

% remaining WECs
ny_WEC = 2+2+par.WEC.ny_rad;
for iWEC = 2:par.NumWECs
iy.p_a = [iy.p_a, iy.p_a(iWEC-1)+ny_WEC*(iWEC-1)];
iy.p_b = [iy.p_b, iy.p_b(iWEC-1)+ny_WEC*(iWEC-1)];
iy.theta = [iy.theta, iy.theta(iWEC-1)+ny_WEC*(iWEC-1)];
iy.theta_dot = [iy.theta_dot, iy.theta_dot(iWEC-1)+ny_WEC*(iWEC-1)];
iy.rad = [iy.rad; iy.rad(iWEC-1,:)+ny_WEC*(iWEC-1)];
iy.WEC = [iy.WEC; iy.WEC(iWEC-1,:)+(2+par.ny_rad)*(iWEC-1)];
end

iy.ny = iy.rad(end);

end