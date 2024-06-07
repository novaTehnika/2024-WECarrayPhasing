function [tout, yout, exitCode] = ode1(F,t0,dt,tfinal,y0,downSampleRate)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ode1.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 
%
% PURPOSE/DESCRIPTION:
% Created from ode1() written by MathWorks
%
% FILE DEPENDENCY: 
% parameters_WECmodel_V02x00.m
%
% UPDATES:
% 05/12/2023 - add down sampling with an error check for the downsamping
% rate as being an integer value.
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
    exitCode = 1;
    % parse optional options
    if exist('downSampleRate','var')
        if mod(downSampleRate,1)
            warning('The down-sampling rate must be an integer value.')
        end
    else
        downSampleRate = 1;
    end
    
    % Generate time array and determine the number of solver time steps
    tout = t0:dt*downSampleRate:tfinal;
    nt = ceil((tout(end)-t0)/dt);

    % Initialize state vector array for output
    yout = zeros(length(tout),length(y0));

    % Specify intial conditions
    t = t0;
    y = y0;
    yout(1,:) = y0;
    
    % Loop through time, starting from second time step
    itds = 2;
    for it = 2:nt
        dydt = F(t,y);      % Calculate the state derivitive vector at 
                            % previous time

        % handle errors
        if any(isnan(dydt),'all')
            exitCode = 2;
            warning('error: (exit code 2) states resulted in NaN for dydt')
            return
        elseif any(imag(dydt),'all')
            exitCode = 3;
            warning('error: (exit code 3) states resulted in imaginary value for dydt')
            return
        end

        t = t + dt;         % Increment time from previous time
        y = y + dydt.*dt;   % Perform first order numerical integration for
                            % state at next time

        yout(itds,:) = y;   % Store the system state at the next time step 
                            % for output


        itds = itds + ~(mod(it-1,downSampleRate));  % increment the current
                                                    % output time step if
                                                    % the curent solver
                                                    % time step matches an
                                                    % output time step
    end
      
end    