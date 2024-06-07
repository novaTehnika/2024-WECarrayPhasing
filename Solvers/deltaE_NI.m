function deltaE = deltaE_NI(p1,p2,cap,dp_abs)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% deltaE_NI.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 5/3/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the change in stored energy for a given intial and final 
% pressure using numerical intergration.  This function does this for a 
% component with a capacitance given by cap().
%
% FILE DEPENDENCY: 
% NA
%
% UPDATES:
% 5/3/2023 - created
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

% integration step size and direction
deltap = p2-p1;
dp = sign(deltap)*dp_abs;

% initialize
p = p1;
deltaE = 0;

% perform numerical integration
while ((deltap>0) && (p<p2)) || ((deltap<0) && (p>p2))
    C = cap(p);
    dE = p*C*dp;
    p = p + dp;
    deltaE = deltaE + dE;
end

end