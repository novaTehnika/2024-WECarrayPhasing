function q = flowCV(deltap,kv,pcrack,pmargin)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flowCV.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/22/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the flow through a check valve (CV).
%
% The flow is assumed to be linear with pressure drop over some margin
% above the cracking pressure and with the square root of above that.
%
% FILE DEPENDENCY:
% NA
%
% UPDATES:
% 6/22/2023 - created.
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

q = max(0, ...
    min(kv*sqrt(pmargin+pcrack)/pmargin*(deltap-pcrack), ...
        sign(deltap)*kv*sqrt(abs(deltap))));

end