function X = areaFracPWM(t,duty,T_sw,tr)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% areaFracPWM.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/21/2023
%
% PURPOSE/DESCRIPTION:
% Calculate the flow coefficient, as a fraction of the maximum, of an 
% on/off valve on a fixed swithcing period with a trapozoidal profile in 
% time.
%
% The duty is the portion of time the valve is open.
%
% The valve is closed for the first portion of the cycle. At 1-duty the
% valve begins to transition to the open position. The valve has
% transistioned to the closed position at the end of the cycle.
% 
% The calculation follows from the function of a line with a slope and
% intercept that switches halfway through the open phase. The lines
% impliment the transition of the valve. The result is then bounded between
% zero and one.
%
% FILE DEPENDENCY:
% NA
%
% UPDATES:
% 6/21/2023 - created.
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

% calculate the current time as fraction within switching period
t_frac = mod(t,T_sw)/T_sw;

% determine phase of cycle
flag = (t_frac > 1-duty/2);

% calculate fraction of maximum flow coefficient
X = max(0,min(1,(1-2*flag)*1/tr.*(t_frac - (1 - duty*(~flag)))));

% if duty is one force valve always open
flag = (duty==1);
X = flag + X*(~flag);

end