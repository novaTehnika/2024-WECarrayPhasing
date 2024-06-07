function C = capAccum(p,pc,Vc,f,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% capAccum.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/3/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the combined capacitance of an accumulator and a dead volume.
% The dead volume allows this function to work for pressures (p) below the
% charge pressure (pc) of the accumulator. The capacitance of the working 
% fluid is assumed to be insignificant at pressures above the charge 
% pressure. As such, accumulated working fluid volume in the accumulator is
% neglected.
%
% The dead volume is defined as a fraction (f) of the total of the 
% accumulator (Vc).
%
% Isothermal compression is assumed for the charged gas in the accumulator
% and the entrained air in the working fluid.
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

% calculate effective bulk modulus of dead volume
beta_eff = par.beta/(1 + par.beta*(par.R*par.p_o/p^2));
% calculate combined capacitance of accumulator and dead volume 
C = (p>pc) * Vc*pc/p^2 ...
    + f*Vc/beta_eff;

end