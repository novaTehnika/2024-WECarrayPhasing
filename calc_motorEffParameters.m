% calc_motorEffParameters.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/17/2024
%
% PURPOSE/DESCRIPTION:
% This script calculates the coefficients for the Wilson pump and motor
% efficiency model, given one data point for volumetric efficiency and one
% data point for mechanical efficiency. An assumption is made for the
% volume fraction coefficient (V_r) and for mechanical efficiency at a
% lower pressure (500kPa).
%
% FILE DEPENDENCY:
%
% UPDATES:
% 7/17/2024 - Adapted from calc_pumpEffParameters.m in the repository at
% https://github.com/novaTehnika/2021Q2-PTOmodeling.
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
clear

rho = 0.860*1000; % [kg/m3] density of working fluid
mu = (32)*1e-6*rho; % [(cSt) -> Pa-s]  Dynamic (absolute) viscosity
beta = 1.8e9; % [Pa]  Bulk Modulus of air free fluid

% Volumetric efficiency loss coefficients
P = 30e6; % [Pa] Pressure for flow rate data point
wq = 2*pi*3600/60; % [rpm -> rad/s] Speed at flow rate data point
eta_v1 = 0.95;

V_r = 1.103; % assumed volume fraction
C_s = mu*wq/P*(1/eta_v1 - 1 - P/beta*(V_r+1));

% Mechanical efficiency loss coefficients
P1 = 30e6; % [Pa] Pressure for torque data point
P2 = 500e3; % [Pa] Pressure for assumed mechanical efficiency data point
wm = 2*pi*3600/60; % [rpm -> rad/s] Speed for torque data point

eta_m1 = 0.95; % P1*D/T1; % efficiency at torque data point
eta_m2 = 0.8*eta_m1; % assumed efficiency at lower pressure

A = [mu*wm/P1, 1; mu*wm/P2, 1];
b = [1 - eta_m1; 1 - eta_m2];
x = A\b;
C_v = x(1);
C_f = x(2);


%% Efficiency as a motor
p = linspace(1e6,P,20);
w = linspace((200)*2*pi/60,wm,20);
for j = 1:length(p)
    for k = 1:length(w)
        eta_v(j,k) = 1/(1 + C_s*(p(j)/(mu*abs(w(k))))...
            + (p(j)/beta)*(V_r + 1));         % Pump volumetric efficiency
        eta_m(j,k) = (1 - C_v*mu*abs(w(k))/...
            p(j) - C_f);                      % Pump mechanical efficiency
    end
end


figure
surf(w,p*1e-6,eta_v)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficieny')
title('Volumetric Efficiency in Motoring Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficiency')
title('Mechanical Efficiency in Motoring Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_v.*eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficiency')
title('Total Efficiency in Motoring Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

%% Efficiency as a pump
p = linspace(1e6,P,20);
w = linspace((200)*2*pi/60,wm,20);
for j = 1:length(p)
    for k = 1:length(w)
        eta_v(j,k) = 1 - C_s*(p(j)/(mu*abs(w(k))))...
            - (p(j)/beta)*(V_r + 1);         % Pump volumetric efficiency
        eta_m(j,k) = 1/(1 + C_v*mu*abs(w(k))/...
            p(j) + C_f);                      % Pump mechanical efficiency
    end
end


figure
surf(w,p*1e-6,eta_v)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficiency')
title('Volumetric Efficiency in Pumping Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficiency')
title('Mechanical Efficiency in Pumping Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])

figure
surf(w,p*1e-6,eta_v.*eta_m)
xlabel('Angular Velocity (rad/s)')
ylabel('Pressure (MPa)')
zlabel('Efficiency')
title('Total Efficiency in Pumping Mode')
xlim([0 w(end)])
ylim(1e-6*[0 p(end)])