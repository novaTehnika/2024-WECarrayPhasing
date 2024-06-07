function y0 = initialConditionDefault_parPTO(par)
% initialConditionDefault_parPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% This script loads default intial conditions for the parPTO model.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 11/2/2023 - Created from initialConditionDefault_refPTO.m.
% 12/6/2023 - Changed from script to function
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

y0(par.iy.theta,1) = 0.1*pi/2; % [rad]
y0(par.iy.theta_dot,1) = 0; % [rad/s]
y0(par.iy.rad,1) = zeros(numel(par.iy.rad),1);

y0(par.iy.p_a,1) = 1e6; % [Pa]
y0(par.iy.p_b,1) = 1e6; % [Pa]

y0(par.iy.p_lin,1) = par.control.p_l_nom; % [Pa]
y0(par.iy.p_lout,1) = par.control.p_l_nom; % [Pa]
y0(par.iy.p_hin,1) = par.control.p_h_nom; % [Pa]
y0(par.iy.p_hout,1) = par.control.p_h_nom; % [Pa]
y0(par.iy.p_ro,1) = par.control.p_h_nom; % [Pa]

y0(par.iy.p_filt,1) = par.control.p_h_nom;
y0(par.iy.errInt_p_filt,1) = 0;

y0(par.iy.LPPL,1) = y0(par.iy.p_lout)*mod(2:2*par.n_seg(1),2)';
y0(par.iy.HPPL,1) = y0(par.iy.p_hin)*mod(2:2*par.n_seg(2),2)';

end