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
% This script loads the state indices for the hydElecPTO model
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

iy.p_a = 1;
iy.p_b = 2;

iy.p_lin = 3;
iy.p_lout = 4;
iy.p_hin = 5;
iy.p_hout = 6;
iy.p_ro = 7;

iy.p_filt = 8;
iy.errInt_p_filt = 9;
iy.control = [iy.p_filt; iy.errInt_p_filt];

iy.theta = 10;
iy.theta_dot = 11;
iy.rad = (1:par.WEC.ny_rad) + iy.theta_dot; % state vector indices for radiation damping states for WEC model
iy.WEC = [iy.theta iy.theta_dot iy.rad];
iy.LPPL = (1:2*par.n_seg(1)-1) + iy.rad(end);
iy.qLP = (1:2:2*par.n_seg(1)-1) + (iy.LPPL(1)-1);
iy.pLP = (2:2:2*par.n_seg(1)-1) + (iy.LPPL(1)-1);

iy.HPPL = (1:2*par.n_seg(2)-1) + iy.LPPL(end);
iy.qHP = (1:2:2*par.n_seg(2)-1) + (iy.HPPL(1)-1);
iy.pHP = (2:2:2*par.n_seg(2)-1) + (iy.HPPL(1)-1);

iy.ny = iy.HPPL(end);

end