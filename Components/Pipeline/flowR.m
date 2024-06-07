function R = flowR(q,lineID,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flowR.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/3/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the flow resistance of a single resistive element in a 
% segmented pipline. The Blassius corolation is assumed for turbulent flow.
% A linear interpolation between laminar and turbulent regimes is used 
% between specified Reynolds numbers.
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

% parameters of the pipe friction model with linear transisiton in 
% transitional region of laminar and turbulent regimes
Re1 = 2300; % transition form laminar to transitional flow
f1 = 64/Re1; % friction factor at transisiton
Re2 = 4500; % transition form transitional to laminar flow
f2  = 0.316*Re2^-0.25; % friction factor at transisiton

% calculate the Re of flow (set floor to avoid zero result)
Re = max(4*par.rho*abs(q)/(par.mu*pi*par.d_line(lineID)),1e-3);

% Calculate darcy friction factor based on Reynolds number
if Re < Re1  % Laminar flow regime
    f = 64/max(Re,0.01);
elseif Re < Re2 % transitional flow 
    f = f1 + (f2-f1)/(Re2 - Re1)*(Re - Re1);
else % turbulent flow regime
    f = 0.316*Re^-.25; % Blasius correlation
end
% calculate flow resistance coefficient
R = f*Re*par.mu*(par.L_line(lineID)/par.n_seg(lineID)) ...
              /(2*par.d_line(lineID)^2*par.A_line(lineID));

end