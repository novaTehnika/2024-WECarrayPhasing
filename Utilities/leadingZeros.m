% leadingZeros.m script function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 08/01/2023
%
% PURPOSE/DESCRIPTION:
% This script creates a format string for the function num2str() for
% integers on the order of the arguement.
%
% FILE DEPENDENCY: N/A
%
% UPDATES:
% 08/01/2023 - Created from implimentation in 
% study_refPTO_chargePumpAccum_wPassiveRV.m
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
function str = leadingZeros(x)
    str = ['%0',num2str(floor(log10(x))),'.f'];
end