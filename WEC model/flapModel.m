function [dydt, torqueFlap, waveElev] = flapModel(t,y, ptoTorque, par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flapModel.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/29/2021
%
% PURPOSE/DESCRIPTION:
% Calculate the state derivatives for a wave energy converter
%
% FILE DEPENDENCY: 
%
% UPDATES:
% 6/29/2021 - created from flapModel.m developed in 2018 by the 
% author. The file is rewritten to exclude the use of
% persistance variables and initialization in favor of calculating
% parameters outside the ODE solver. The file is also written for use in
% a variable timestep ODE solver.
% 7/9/2021 - added an exponential ramp in the excitation force so that the
% simulation does not experience an impulse at t = 0. 
% 08/02/2022 - added ramp period to par struct and included tstart so that
% ramp ends at t=0;
%
% Copyright (C) 2022  Jeremy W. Simmons II
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

iyrad = 3:2+par.WEC.ny_rad; % state vector indices for radiation damping states for WEC model

% Wave excitation force
waveElev = waveElevation(t,par);%par.WEC.waveElevation(t);
    
% Excitation Torque
t0 = par.tstart-par.Tramp;
torqueFlap.wave = ramp(t-t0,par.TrampWEC)*excitationTorque(t,par);%par.WEC.excitationTorque(t);
    
% Hydrostatic torque
torqueFlap.hydroStatic = -hydroStaticTorque(y(1), waveElev, par);

% Radiation torque
torqueFlap.radiation = -y(iyrad(1));
dydt(iyrad) =  par.WEC.A_rad*y(iyrad) + par.WEC.B_rad*y(2);

% flap motion
thetaLim = 8/9*pi/2;
endDamping = 1e9;
dydt(1) = y(2);
dydt(2) = 1/(par.WEC.I + par.WEC.I_inf) ...
                *(torqueFlap.wave + torqueFlap.radiation ...
                + torqueFlap.hydroStatic + ptoTorque - endDamping*y(2)* ...
                ( (y(1)>thetaLim)*(y(2)>0) ... 
                + (y(1)<-thetaLim)*(y(2)<0)));

end
