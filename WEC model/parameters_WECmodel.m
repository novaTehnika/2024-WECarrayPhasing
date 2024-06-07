function par = parameters_WECmodel(par,filenameCoeff,filenameRadSS)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters_WECmodel.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/29/2021
%
% PURPOSE/DESCRIPTION:
% This function loads prepares parameters used in the hydrodynamcic WEC
% model. Parameters are loaded from two files. One contains hydrodynamic
% data that are outputs of boundary element method code and the other
% contains general parameters and coefficients of a radiation damping model
% calculated seperately from the BEM code results. This funciton calculates
% the wave specrum based on the significant wave height and peak period,
% and prepares the matrices for the radiation damping, state-space model.
%
% FILE DEPENDENCY: 
% ./WECdata/
%   nemohResults_vantHoff2009_20180802.mat
%   vantHoffTFCoeff.mat
%
% UPDATES:
% 6/29/2021 - created from buoyTorque.m developed in 2018 
% by the author. The file is rewritten to take a 
% structure of parameters as an input instead of using hard coded 
% variables.
% 7/9/2021 - implemented the equal-area (equal-energy) method of 
% generating wave elevation timeseries and excitation force. This method 
% involves interpolation of the excitation force. Additions include the use
% of nested functions for the wave spectrum and obj. func.
% 8/23/2022 - added position measurement for center of mass and corrected
% moment of inertia calculation to include effect of parallel axis (value
% for MOI specified is about center of mass).
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

% Load saved WEC parameters
% load('nemohResults_vantHoff2009_20180802.mat','OSWC_vantHoff2009')
load(filenameCoeff,'OSWC_vantHoff2009')
par.imprt.WEC.w = OSWC_vantHoff2009.w;          % [rad/s] frequency
par.imprt.WEC.F_amp = OSWC_vantHoff2009.F_amp;  % [Nm/m] frequency dependent excitation torque amplitude coefficient
par.imprt.WEC.phi_e = OSWC_vantHoff2009.phi_e;  % [rad] frequency dependent phase of excitation torque
par.imprt.WEC.I_inf = OSWC_vantHoff2009.A_inf(1);  % [kgm^2] added mass at infinite frequency
par.WEC.I_inf = par.imprt.WEC.I_inf;
% load('vantHoffTFCoeff_v2.mat','KradDen','KradNum')
load(filenameRadSS,'KradDen','KradNum')
par.WEC.a = KradDen';     % transfer function denominator coefficients for raditation damping derived to approximate the hydrodynamic response
par.WEC.b = KradNum';     % transfer function numerator coefficients for raditation damping derived to approximate the hydrodynamic response
par.WEC.L_flap = 11;      % [m] length of flap (hinge to top edge)
par.WEC.T = 2;            % [m] flap thickness
par.WEC.h = 8.9;          % [m] mean water depth to hinge
par.WEC.H = 11;           % [m] mean water depth to sea bed
par.WEC.W = 18;           % [m] flap width
par.WEC.L_cm = 5;         % [m] distance of center of mass from hinge
par.WEC.m = 127000;       % [kg] mass of flap
I = 1.85e6;     % [kgm^2] mom. of inertia of flap about center of mass
par.WEC.I = I + par.WEC.m*par.WEC.L_cm^2; % [kgm^2] mom. of inertia of flap about hinge
par.WEC.rho = 1025;       % [kg/m^3] density of sea water
par.WEC.g = 9.81;         % [m/s^2] gravitational acceleration at sea level

% construct matrices for radiation damping state-space model
 % number of states describing radiation damping memory
par.WEC.ny_rad = length(par.WEC.a)-1; 

% Cast linear radiation model to Observable Cononical Form
A1 = [ eye(par.WEC.ny_rad-1); ...
       zeros(1,par.WEC.ny_rad-1)] ; 
par.WEC.A_rad  = [-par.WEC.a(2:end) A1];
par.WEC.B_rad = [zeros(par.WEC.ny_rad-length(par.WEC.b),1); ...
                 par.WEC.b];
             
             
% determine new freqencies based on equal energy method
 % determine the area under the curve of the wave spectrum for
 % discretization
  % high density grid of frequencies evenly spaced
w = linspace(par.imprt.WEC.w(1),par.imprt.WEC.w(end),1e6); 
  % wave spectrum calculated for high-density even grid
S_w = PiersonSpec(w,par);
  % area for each bin
a = trapz(w,S_w)/(par.WEC.nw+1);

  % calculate the cumulative area under the curve
A = cumtrapz(w,S_w);

 % determine the frequencies for an equal area grid by finding upper and
 % lower bounds for the bins
par.WEC.w = zeros(par.WEC.nw,1);
par.WEC.dw = zeros(par.WEC.nw,1);
w_lb = zeros(par.WEC.nw,1); % lower bound for each bin
w_ub = zeros(par.WEC.nw,1); % upper bound for each bin
  % min freq. set to a significant value
wmin = w(find(A >= 0.01*a,1,'first')); 
  % max freq. set to a significant value
wmax = w(find(A <= (A(end) - 0.01*a),1,'last')); % max freq. set to a significant value
for iw = 1:par.WEC.nw
    % set lower bound
    if iw == 1
        w_lb(iw) = wmin;
    else
        w_lb(iw) = w_ub(iw-1);
    end
    % find upper bound to satisfy area increment
    if iw == par.WEC.nw
        w_ub(iw) = wmax;
    else
       w_ub(iw) = w(find(A >= a*iw,1,'first'));
    end
	par.WEC.dw(iw) = w_ub(iw) - w_lb(iw);
    par.WEC.w(iw) =  (w_ub(iw) + w_lb(iw))/2;
end

% calculate the wave spectrum for the equal area frequencies
par.wave.S_w = PiersonSpec(par.WEC.w,par);

% interpolate frequency dependent coefficients
 % frequency dependent excitation torque amplitude coefficient
par.WEC.F_amp = ...
    interp1(par.imprt.WEC.w,par.imprt.WEC.F_amp,par.WEC.w,'linear');

 % frequency dependent phase of excitation torque
par.WEC.phi_e = ...
    interp1(par.imprt.WEC.w,par.imprt.WEC.phi_e,par.WEC.w,'linear');


% generate random phases for each frequency component for the wave elevation
 % seed the random number generator
rng(par.wave.rngSeedPhase); 

 % generate uniformly distributed random numbers [-pi,pi]
par.wave.phi = 2*pi*(rand(par.WEC.nw,1)-0.5);


% determine phase of wave elevation freq. components that reflect the
% position of each WEC
 % generate wave numbers as a function of frequency
for iw = 1:nw
    par.K(iw,1) = fzero(@(K) ...
                    K*par.WEC.g*tanh(K*par.WEC.H) - par.WEC.w(iw)^2, ...
                                        par.WEC.w(iw)^2/par.WEC.g);
    par.k(iw,1) = K(iw)*cos(par.wave.waveDirection); % x component of wave number vector
    par.l(iw,1) = K(iw)*sin(par.wave.waveDirection); % y component of wave number vector
end

% add phase due to position
for iWEC = 1:par.NumWECs
    par.WEC.phi(iWEC,:) = par.wave.phi - (k*par.WEC.x(iWEC) + l*par.WEC.y(iWEC));
end

% regular waves (supersedes previous calculations)
if par.WEC.nw == 1
    par.wave.S_w = 0.5*par.wave.Hs^2;
    par.WEC.dw = 1;
    par.WEC.w = 2*pi/par.wave.Tp;
    par.WEC.F_amp = ...
        interp1(par.imprt.WEC.w,par.imprt.WEC.F_amp,par.WEC.w,'linear');
    par.WEC.phi_e = ...
        interp1(par.imprt.WEC.w,par.imprt.WEC.phi_e,par.WEC.w,'linear');

    par.K = fzero(@(K) K*par.WEC.g*tanh(K*par.WEC.H) - par.WEC.w^2, ...
                                        par.WEC.w^2/par.WEC.g);
    par.k = K*cos(par.wave.waveDirection); % x component of wave number vector
    par.l = K*sin(par.wave.waveDirection); % y component of wave number vector
    par.WEC.phi = - (k*par.WEC.x(iWEC) + l*par.WEC.y(iWEC));
end

%% %%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%
    function S_w = PiersonSpec(w,par)
        % Based on Falnes (2002) "Ocean Waves and Oscillating Systems:..."
        S_w = 10*pi^5*par.wave.Hs^2/par.wave.Tp^4./w.^5 ... 
            .*exp(-20*pi^4/par.wave.Tp^4./w.^4)/(2*pi);
    end

end