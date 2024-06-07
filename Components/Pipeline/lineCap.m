function C = lineCap(p,lineID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lineCap.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/3/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the capacitance of a single pressure node in a segmented
% pipline. Isothermal compression of the entrained air is assumed.
%
% FILE DEPENDENCY: 
% NA
%
% UPDATES:
% 5/3/2023 - created from implimentation in pipelineNPi.m
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

% calculate capacitance in pipeline lump

V_line = par.A_line(lineID)*par.L_line(lineID)/par.n_seg(lineID); % volume of each lump

% calculate effective bulk modulus 

 % via Cho method (as arranged in Yudell, 2017)
%         beta_eff = par.beta* ...
%         (((p/par.p_o)^(1/par.gamma)*exp((par.p_o-p)/par.beta)+par.R) / ...
%         (par.R/par.gamma*par.beta/p+(p/par.p_o)^(1/par.gamma)*...
%         exp((par.p_o-p)/par.beta))); 
 
 % isothermal bulk modulus
beta_eff = par.beta/(1 + par.beta*(par.R*par.p_o/p^2));

% calculate capacitance
C = V_line/beta_eff;

end