% waveElev.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This script tests a model of a wave field that provides a amplitude and
% phase parameters as a function of frequency for a trigometric. The model
% implements the Pierson-Moskowitz wave spectrum and is parameterized by
% the significant wave height, peak period, and direction of travel for the
% waves.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 6/8/2023 - Created from initialConditionDefault_seriesPTO.m.
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
par.wave.Hs = 2.5;
par.wave.Tp = 8.7;

par.wave.rngSeedPhase = 3; % set the seed for the random number generator
par.WEC.nw = 100;


H = 10; % [m] mean water depth
g = 9.81; % [m/s^2] gravity const.
waveDirection = (0)*pi/180; % [(deg)->rad] direction of wave travel wrt x-axis

psi = 0; % [rad] phase of frequency component at [x,y] = [0,0]


%%
% determine new freqencies based on equal energy method
 % determine the area under the curve of the wave spectrum for
 % discretization
  % high density grid of frequencies evenly spaced
w = linspace(0.1,10,1e4); 
dw = w(2) - w(1); % grid spacing of frequencies
  % wave spectrum calculated for high-density even grid
S_w = PiersonSpec(w,par);
  % area for each bin
a = trapz(w,S_w)/(par.WEC.nw+1);
A = cumtrapz(w,S_w);
        
 % determine the frequencies of the equal area grid using obj. function
par.WEC.w = zeros(par.WEC.nw,1);
par.WEC.dw = zeros(par.WEC.nw,1);
wmax = w(end);
Atarget = 0;
options = optimset('Display','off');
for iw = 1:par.WEC.nw
    if iw == 1 % set wmin to min frequency in range considered
        wmin = w(1);
    else % set wmin to right bound of previous bin
        wmin = par.WEC.w(iw-1) + par.WEC.dw(iw-1)/2;
    end
    Atarget = a*iw;

    par.WEC.w(iw) = ...
        fzero(@(w_ea) Atarget - interp1(w,A,w_ea,'spline'),...
        (wmax - wmin)/10,options);

	par.WEC.dw(iw) = a/PiersonSpec(par.WEC.w(iw),par);

end

% calculate the wave spectrum for the equal area frequencies
par.wave.S_w = PiersonSpec(par.WEC.w,par);

% generate random phases for each frequency component for the wave elevation
% seed the random number generator
rng(par.wave.rngSeedPhase); 
% generate uniformly distributed random numbers [-pi,pi]
par.wave.phi = 2*pi*(rand(par.WEC.nw,1)-0.5);

a = sqrt(2*par.wave.S_w(:).*par.WEC.dw); % [m] wave amplitude
w = par.WEC.w;
nw = numel(w);
%%

% Wave number
for iw = 1:nw
    K(iw,1) = fzero(@(K) K*g*tanh(K*H) - w(iw)^2, w(iw)^2/g);
    k(iw,1) = K(iw)*cos(waveDirection); % x component of wave number vector
    l(iw,1) = K(iw)*sin(waveDirection); % y component of wave number vector
end


t = 0:0.1:100;
nt = numel(t);
x = 0:0.5:20;
nx = numel(x);
y = 0;
ny = numel(y);

z_surf = zeros(nx,ny,nt);
for ix = 1:nx
    for iy = 1:ny
        for it = 1:nt
            z_surf(ix,iy,it) = sum(a.*cos(k*x(ix) + l*y(iy) - w*t(it) + par.wave.phi));
        end
    end
end

%%
fig = figure('Units','inches',"Position",[1,1,3,2]);
clearvars F
for it = 1:nt
    plot(x,z_surf(:,:,it),'k-',[x(1),x(end)],[-H,-H],'LineWidth',2)
%     surf(y,x,z_surf(:,:,it),'EdgeColor','none')
%     view(2)
ylim([-11,3])
    F(it) = getframe(fig);
end

%%
fig = figure('Units','inches',"Position",[1,1,3,2]);
movie(fig,F,1,1/(t(2)-t(1)))

%% test wave elevation function as implimented in WEC model
t = 0:0.1:1200;

x = [0 30];
y = [0 15];
par.WEC.phi(:,1) = par.wave.phi - (k*x(1) + l*y(1));
par.WEC.phi(:,2) = par.wave.phi - (k*x(2) + l*y(2));

for it = 1:numel(t)
    waveElev_tseries(1,it) = waveElevation(t(it),par,1);
    waveElev_tseries(2,it) = waveElevation(t(it),par,2);
end


figure
plot(t,waveElev_tseries(1,:),'LineWidth',1.5)
hold on

plot(t,waveElev_tseries(2,:),'LineWidth',1.5)

% plot(t,(waveElev_tseries(1,:) + waveElev_tseries(2,:))./2,'k-','LineWidth',1.5)

xlim([0 100])

xlabel('time (s)')
ylabel('elevation (m)')

legend('at WEC A','at WEC B','average')

title('Wave Elevation at WECs Placed Out of Phase')
%%
function S_w = PiersonSpec(w,par)
    % Based on Falnes (2002) "Ocean Waves and Oscillating Systems:..."
    S_w = 10*pi^5*par.wave.Hs^2/par.wave.Tp^4./w.^5 ... 
        .*exp(-20*pi^4/par.wave.Tp^4./w.^4)/(2*pi);
end