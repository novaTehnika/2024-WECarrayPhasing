function [out, exitCode] = sim_hydElecPTO(y0,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_hydElecPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 6/7/2024
%
% PURPOSE/DESCRIPTION:
% This script executes the set-up, solution, and basic post-processing for
% the hydElecPTO model.
%
% FILE DEPENDENCY:
% ../Hydraulic-Electric PTO/
%   sys_hydElecPTO.m
%   stateIndex_hydElecPTO.m
% ../WEC model/
%   flapModel.m
%   hydroStaticTorque.m
% ../Solvers/
%   deltaE_NI.m
%   deltaV_NI.m
%   ode1.m
% ../Components/
%   areaFracPWM.m
%   capAccum.m
%   deadVCap.m
%   flowCV.m
%   flowPRV.m
% ../Components/Pipeline
%   flowR.m
%   lineCap.m
%   pipelineNPi.m
% ./Utilities/
%   startParPool.m
%
% UPDATES:
% 6/7/2024 - Created from sim_parPTO.m.
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
exitCode = 1;

%% %%%%%%%%%%%%   SOLUTION   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solve for states of the dynamic system

switch par.solver
    case 'fixed time'
        dt = par.MaxStep;
        downSampleRate = floor(par.downSampledStepSize/dt);

     % Run solver
        tspan = [par.tstart-par.Tramp par.tend];   % time interval
        ticODE = tic;
        [t, y, exitCode] = ode1(@(t,y) sys(t,y,par), ...
                                  tspan(1),dt,tspan(2),y0,downSampleRate);
        toc(ticODE)

    case 'variable time'
 % Solver options
    options = odeset('RelTol',1e-4,...
                     'AbsTol',1e-4,...
                     'MaxOrder',3);
    options.MaxStep = par.MaxStep;
    options.InitialStep = par.stepSizeWECramp;
    dt = par.stepSizeWECramp;
    downSampleRate = floor(par.downSampledStepSize/dt);

 % Run 1st order solver through ramp period for WEC
    tspan = [par.tstart-par.Tramp ...
             par.tstart-(par.Tramp-par.TrampWEC)];   % time interval
    ticODE = tic;
    [t, y] = ode1(@(t,y) sys(t,y,par), ...
                                tspan(1),dt,tspan(2),y0,downSampleRate);
    toc(ticODE)

 % Run variable order solver for time after ramp period for WEC
    tspan = [tspan(2) par.tend];   % time interval
    ticODE = tic;
    [t, y] = ode15s(@(t,y) sys(t,y,par), ...
                                tspan,y(end-1,:)',options);
    toc(ticODE)
end

handle errors
if exitCode == 1 % normal operation, no error in solver
    % check for negative pressure values
    if any(y(:,[par.iy.p_a, par.iy.p_b, par.iy.p_l]) < 0,'all')
        warning('Negative pressures detected.')
        exitCode = 4;
        out = [];
        return
    end
else % error, return empty
    out = [];
    return
end

%% %%%%%%%%%%%%   POST-PROCESS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Parameters
    out.par = par;
    
    % Select desired time indices
    itVec = find(t >= par.tstart);
    
    % Extract system states from simulation results post ramp
    out.t = t(itVec);
    out.y = y(itVec,:);
    
    out.p_l = y(itVec,par.iy.p_l);
    out.p_h = y(itVec,par.iy.p_h);
    
    out.control.p_filt = y(itVec,par.iy.p_filt);
    out.control.errInt_p_filt = y(itVec,par.iy.errInt_p_filt);
    
    for iWEC = 1:par.NumWECs
        out.p_a(:,iWEC) = y(itVec,par.iy.p_a(iWEC));
        out.p_b(:,iWEC) = y(itVec,par.iy.p_b(iWEC));
        out.theta(:,iWEC) = y(itVec,par.iy.theta(iWEC)); % [rad] position of the WEC
        out.theta_dot(:,iWEC) = y(itVec,par.iy.theta_dot(iWEC)); % [rad/s] angular velocity of the WEC
    end

    % Post-process non-state variables and state derivatives
    syspost = @(t,y,par) sysPost(t,y,par);

    nt_ramp = itVec(1)-1;
    % startParPool
    ntVec = numel(itVec);
    dydt = zeros(ntVec,par.iy.ny);
    for it = ntVec:-1:1
        it_star = it+nt_ramp;
        t_star = t(it_star);
        y_star = y(it_star,:)';

        [dydt(it,:), nonState(it), control(it)] = ...
                            syspost(t_star,y_star,par);
        
        % Move WEC/WEC-driven pump results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        for iWEC = 1:par.NumWECs
            % System input: wave elevation
            temp(it,iWEC).waveElev = nonState(it).waveElev(iWEC);

            % forces on WEC
            temp(it,iWEC).T_pto = nonState(it).T_pto(iWEC);
            temp(it,iWEC).T_hydroStatic = nonState(it).torqueWEC(iWEC).hydroStatic;
            temp(it,iWEC).T_wave = nonState(it).torqueWEC(iWEC).wave;
            temp(it,iWEC).T_rad = nonState(it).torqueWEC(iWEC).radiation;

            % WEC-driven pump flow
            temp(it,iWEC).q_aIn = nonState(it).q_aIn(iWEC);
            temp(it,iWEC).q_aOut = nonState(it).q_aOut(iWEC);
            temp(it,iWEC).q_bIn = nonState(it).q_bIn(iWEC);
            temp(it,iWEC).q_bOut = nonState(it).q_bOut(iWEC);
        end

    end
    
    % State derivatives
    out.dydt = dydt;
    
    % WEC/WEC-driven pump variables
    for iWEC = 1:par.NumWECs
        % System input: wave elevation
        out.waveElev(:,iWEC) = [temp(:,iWEC).waveElev]';

        % forces on WEC
        out.T_pto(:,iWEC) = [temp(:,iWEC).T_pto]';
        out.T_hydroStatic(:,iWEC) = [temp(:,iWEC).T_hydroStatic]';
        out.T_wave(:,iWEC) = [temp(:,iWEC).T_wave]';
        out.T_rad(:,iWEC) = [temp(:,iWEC).T_rad]';

        % WEC-driven pump flow
        out.q_aIn(:,iWEC) = [temp(:,iWEC).q_aIn]';
        out.q_aOut(:,iWEC) = [temp(:,iWEC).q_aOut]';
        out.q_bIn(:,iWEC) = [temp(:,iWEC).q_bIn]';
        out.q_bOut(:,iWEC) = [temp(:,iWEC).q_bOut]';

        out.q_hwp(:,iWEC) = out.q_aOut(:,iWEC) + out.q_bOut(:,iWEC);
        out.q_lwp(:,iWEC) = out.q_aIn(:,iWEC) + out.q_bIn(:,iWEC);
    end

    % Control signals
    out.control.w_m = [control(:).w_m]';

    % torque on pump/motor and generator shafts
    out.T_m = [nonState(:).T_m]';
    out.T_gen = out.T_m;

    % pump/motor flow
    out.w_m = out.control.w_m;
    out.q_m = [nonState(:).q_m]';

    % Charge Pump
    out.q_c = [nonState(:).q_c]';

     % Pressure relief valves at system pressure nodes
    out.q_lPRV = [nonState(:).q_lPRV]';
    out.q_hPRV = [nonState(:).q_hPRV]';

%% Post-process analysis

    %% Energy analysis

    % WEC-driven pump(s)
    for iWEC = 1:par.NumWECs

        out.power.P_WEC(:,iWEC) = -out.T_pto(:,iWEC).*out.theta_dot(:,iWEC);
        out.power.P_wp(:,iWEC) = out.p_h.*out.q_hwp(:,iWEC) - out.p_l.*out.q_lwp(:,iWEC);
        out.power.P_wpLoss(:,iWEC) = out.power.P_WEC(:,iWEC) - out.power.P_wp(:,iWEC);

     % check valve rectifier
        out.power.P_aIn(:,iWEC) = out.q_aIn(:,iWEC).*(out.p_l-out.p_a(:,iWEC));
        out.power.P_aOut(:,iWEC) = out.q_aOut(:,iWEC).*(out.p_a(:,iWEC)-out.p_h);
        out.power.P_bIn(:,iWEC) = out.q_bIn(:,iWEC).*(out.p_l-out.p_b(:,iWEC));
        out.power.P_bOut(:,iWEC) = out.q_bOut(:,iWEC).*(out.p_b(:,iWEC)-out.p_h);
    end

    % Pump/motor and generator
    out.power.P_mLoss = (out.p_l - out.p_h).*out.q_m ...
                       - out.T_m.*out.w_m;
    out.power.P_gen = -((-out.T_gen.*out.w_m > 0).*par.eta_g ...
                    + (-out.T_gen.*out.w_m < 0)./par.eta_g) ...
                    .*out.T_gen.*out.w_m;
    out.power.P_genLoss = -out.T_gen.*out.w_m - out.power.P_gen;

    % Charge pump
    out.power.P_cElec = 1/(par.eta_c*par.eta_m)*out.q_c.*(out.p_l-out.par.p_o);
    out.power.P_cLoss = out.power.P_cElec - out.q_c.*(out.p_l-out.par.p_o);

    % Pressure relief valves
    out.power.P_lPRV = out.q_lPRV.*(out.p_l-out.par.p_o);
    out.power.P_hPRV = out.q_hPRV.*(out.p_h-out.p_l);

    %% Energy Balance

    % Change in available potential energy in WEC-driven pump chambers
    dp = 1e2;
    cap = @(p,V) deadVCap(p,V,par);
    Va = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    Vb = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max + theta);
    deltaE_a = zeros(par.NumWECs,1);
    deltaE_b = zeros(par.NumWECs,1);
    for iWEC = 1:par.NumWECs
         % Chamber 'a'
        cap1 = @(p) cap(p,Va(out.theta(1,iWEC)));
        capend = @(p) cap(p,Va(out.theta(end,iWEC)));
        deltaE_a(iWEC) = deltaE_NI(out.par.p_o,out.p_a(end,iWEC),capend,dp) ...
                    - deltaE_NI(out.par.p_o,out.p_a(1,iWEC),cap1,dp);
    
         % Chamber 'b'
        cap1 = @(p) cap(p,Vb(out.theta(1,iWEC)));
        capend = @(p) cap(p,Vb(out.theta(end,iWEC)));
        deltaE_b(iWEC) = deltaE_NI(out.par.p_o,out.p_b(end,iWEC),capend,dp) ...
                    - deltaE_NI(out.par.p_o,out.p_b(1,iWEC),cap1,dp);
    end

    % Change in available potential energy in accumulators
    dp = 1e2;
     % Low-pressure (charge pump outlet and inlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_l,par.Vc_l,par.f,par);
    deltaE_l = deltaE_NI(out.p_l(1),out.p_l(end),cap,dp);

     % High-pressure(outlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_h,par.Vc_h,par.f,par);
    deltaE_h = deltaE_NI(out.p_h(1),out.p_h(end),cap,dp);

    % Total change in stored energy in the system
    deltaE_sys = sum(deltaE_a + deltaE_b) ...
        + deltaE_l + deltaE_h;

    % Power flow at boundaries
    P_in = sum(out.power.P_WEC,2);
    P_out = out.power.P_gen - out.power.P_cElec;
    P_loss = sum(out.power.P_wpLoss,2) ...
            + sum(out.power.P_aIn + out.power.P_aOut,2) ...
            + sum(out.power.P_bIn + out.power.P_bOut,2) ...
            + out.power.P_mLoss + out.power.P_genLoss ...
            + out.power.P_cLoss ...
            + out.power.P_lPRV + out.power.P_hPRV;
    P_bnds = P_in - P_out - P_loss;

    % Total balance of energy: energy added minus change in energy stored
    out.Ebal = trapz(out.t,P_bnds) - deltaE_sys;
    out.Ebal_error = out.Ebal/trapz(out.t,P_in);

    %% Mass Balance
    % Equivalent change in fluid volume in WEC-driven pump chambers
    dp = 1e2;
    cap = @(p,V) deadVCap(p,V,par);
    Va = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    Vb = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    deltaV_a = zeros(1,par.NumWECs);
    deltaV_b = zeros(1,par.NumWECs);
    for iWEC = 1:par.NumWECs
         % Chamber 'a'
        cap1 = @(p) cap(p,Va(out.theta(1,iWEC)));
        capend = @(p) cap(p,Va(out.theta(end,iWEC)));
        deltaV_a(iWEC) = deltaV_NI(out.par.p_o,out.p_a(end,iWEC),capend,dp) ...
                    - deltaV_NI(out.par.p_o,out.p_a(1,iWEC),cap1,dp);

         % Chamber 'b'
        cap1 = @(p) cap(p,Vb(out.theta(1,iWEC)));
        capend = @(p) cap(p,Vb(out.theta(end,iWEC)));
        deltaV_b(iWEC) = deltaV_NI(out.par.p_o,out.p_a(end,iWEC),capend,dp) ...
                    - deltaV_NI(out.par.p_o,out.p_a(1,iWEC),cap1,dp);
    end

    % Change in fluid volume in accumulators
    dp = 1e2;
     % Low-pressure (charge pump outlet and inlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_l,par.Vc_l,par.f,par);
    deltaV_l = deltaV_NI(out.p_l(1),out.p_l(end),cap,dp);

     % High-pressure(outlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_h,par.Vc_h,par.f,par);
    deltaV_h = deltaV_NI(out.p_h(1),out.p_h(end),cap,dp);

    % Total change in stored volume in the system
    out.deltaV_total = sum(deltaV_a + deltaV_b,1) ...
                     + deltaV_l + deltaV_h;

    % Flow at boundaries
    qbnds = out.q_c - out.q_lPRV;

    % Total balance of flow: flow in minus change in volume stored
    out.Vbal = trapz(out.t,qbnds) - out.deltaV_total;
    out.Vbal_error = out.Vbal/trapz(out.t,out.q_m);

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ , ~ ] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sysPost(t,y,par)
    	[dydt, nonState, control] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sys_prototype(t,y,par)
        [dydt, nonState, control] = sys_hydElecPTO(t,y,par);
    end
end
