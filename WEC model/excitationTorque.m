function T_e = excitationTorque(t,par,iWEC)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% excitationTorque.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/6/2023
%
% PURPOSE/DESCRIPTION:
% This function calculates the excitation torque on an oscilating
% surge wave converter.
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 12/6/2023 - created from flapModel.m
% 6/7/2024 - modified to account for multiple WECs
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

T_e = sum( par.WEC.F_amp.*sqrt(2*par.wave.S_w.*par.WEC.dw) ...
               .*sin(par.WEC.w*t + par.WEC.phi(:,iWEC) + par.WEC.phi_e) );

end