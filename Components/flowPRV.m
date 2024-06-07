function q = flowPRV(p,p_crack,C)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flowPRV.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/3/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the flow through a pressure relief valve. A poppet-style valve
% with linear return spring is assumed. This yeilds an analytical solution 
% for flow as a function of pressure. The flow constant is calculated from
% a specified flow rate at a given pressure, along with a margin that
% specifies the difference between the given pressure and the cracking 
% pressure. The following code gives an example for calculating the
% parameters of the model.
%
%     maxPressure = 8.3e6; % [Pa]
%     margin = 5e4; % [Pa]
%     maxFlow = 100e-3; % [m^3/s]
%     p_crack = maxPressure - margin;
%     C = (maxPressure^(3/2) ...
%                  - (maxPressure-margin)*maxPressure^(1/2))/maxFlow;
%
% FILE DEPENDENCY: 
% NA
%
% UPDATES:
% 5/3/2023 - created from implimentation in sys_seriesPTO.m
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

q = 1/C*max(0,p.^(3/2)-p_crack*p.^(1/2));

end