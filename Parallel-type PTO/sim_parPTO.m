function [out, exitCode] = sim_parPTO(y0,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_parfPTO.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/2/2023
%
% PURPOSE/DESCRIPTION:
% This script executes the set-up, solution, and basic post-processing for
% the parPTO model.
%
% FILE DEPENDENCY:
% ../Parallel-type PTO/
%   sys_parPTO.m
%   stateIndex_parPTO.m
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
% 11/2/2023 - Created from sim_refPTO.m.
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

% handle errors
if exitCode == 1 % normal operation, no error in solver
    % check for negative pressure values
    if any(y(:,[par.iy.p_a, par.iy.p_b, ...
            par.iy.p_lin, par.iy.p_lout]) < 0,'all')
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
    
    out.p_a = y(itVec,par.iy.p_a);
    out.p_b = y(itVec,par.iy.p_b);
    
    out.p_lin = y(itVec,par.iy.p_lin);
    out.p_lout = y(itVec,par.iy.p_lout);
    out.p_hin = y(itVec,par.iy.p_hin);
    out.p_hout = y(itVec,par.iy.p_hout);
    out.p_ro = y(itVec,par.iy.p_ro);
    
    out.control.p_filt = y(itVec,par.iy.p_filt);
    out.control.errInt_p_filt = y(itVec,par.iy.errInt_p_filt);
    
    out.theta = y(itVec,par.iy.theta); % [rad] position of the WEC
    out.theta_dot = y(itVec,par.iy.theta_dot); % [rad/s] angular velocity of the WEC

    out.yLP = y(itVec,par.iy.LPPL);
    out.qLP = y(itVec,par.iy.qLP);
    out.pLP = y(itVec,par.iy.pLP);
    out.q_lin = out.qLP(:,1);
    out.q_lout = out.qLP(:,end);
    
    out.yHP = y(itVec,par.iy.HPPL);
    out.qHP = y(itVec,par.iy.qHP);
    out.pHP = y(itVec,par.iy.pHP);
    out.q_hin = out.qHP(:,1);
    out.q_hout = out.qHP(:,end);
            
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
        
        % Move WEC torque results up a level in nonState stucture so that
        % they can be used like arrays in assiging to the output structure
        temp(it).T_hydroStatic = nonState(it).torqueWEC.hydroStatic;
        temp(it).T_wave = nonState(it).torqueWEC.wave;
        temp(it).T_rad = nonState(it).torqueWEC.radiation;

        [~ ,~, ~, pLsoln_LP(it)] = ...
                        pipelineNPi([],y_star(par.iy.LPPL), ...
                        y_star(par.iy.p_lin),y_star(par.iy.p_lout),par,1,1);
        [~ ,~, ~, pLsoln_HP(it)] = ...
                        pipelineNPi([],y_star(par.iy.HPPL), ...
                        y_star(par.iy.p_hin),y_star(par.iy.p_hout),par,2,1);
    end
    
    % State derivatives
    out.dydt = dydt;
    
     % System input: wave elevation
    out.waveElev = [nonState(:).waveElev]';
    
     % forces on WEC
    out.T_pto = [nonState(:).T_pto]';
    out.T_hydroStatic = [temp(:).T_hydroStatic]';
    out.T_wave = [temp(:).T_wave]';
    out.T_rad = [temp(:).T_rad]';
          
     % Control signals
    out.control.w_pm = [control(:).w_pm]';
    out.control.kv_rv = [control(:).kv_rv]';
    out.control.kv_ideal = [control(:).kv_ideal]';

     % torque on pump/motor and generator shafts
    out.Tpm = [nonState(:).Tpm]';
    out.Tgen = out.Tpm;
    
     % WEC-driven pump flow
    out.q_ain = [nonState(:).q_ain]';
    out.q_aout = [nonState(:).q_aout]';
    out.q_bin = [nonState(:).q_bin]';
    out.q_bout = [nonState(:).q_bout]';
    out.q_sv = [nonState(:).q_sv]';

    out.q_hwp = [nonState(:).q_aout]' + [nonState(:).q_bout]';
    out.q_lwp = [nonState(:).q_ain]' + [nonState(:).q_bin]';

     % pump/motor flow
    out.w_pm = out.control.w_pm;
    out.q_pm = [nonState(:).q_pm]';
    
     % RO performance
    out.q_perm = [nonState(:).q_perm]';
    out.q_brine = [nonState(:).q_brine]';
    out.q_feed = out.q_perm + out.q_brine;

     % RO inlet valve
    out.q_rv = [nonState(:).q_rv]';

     % ERU flow rate
    out.q_ERUfeed = [nonState(:).q_ERUfeed]';

     % Charge Pump
    out.q_c = [nonState(:).q_c]';

     % Pressure relief valves at system pressure nodes
    out.q_linPRV = [nonState(:).q_linPRV]';
    out.q_hinPRV = [nonState(:).q_hinPRV]';
    out.q_roPRV = [nonState(:).q_roPRV]';

%% Post-process analysis

    %% Energy analysis
    % WEC-driven pump
    out.power.P_WEC = -out.T_pto.*out.theta_dot;
    out.power.P_wp = out.p_hin.*out.q_hwp - out.p_lout.*out.q_lwp;
    out.power.P_wpLoss = out.power.P_WEC - out.power.P_wp;
    
     % switching valve
    out.power.P_sv = out.q_sv.*(out.p_a-out.p_b);
    
     % check valve rectifier
    out.power.P_ain = out.q_ain.*(out.p_lout-out.p_a);
    out.power.P_aout = out.q_aout.*(out.p_a-out.p_hin);
    out.power.P_bin = out.q_bin.*(out.p_lout-out.p_b);
    out.power.P_bout = out.q_bout.*(out.p_b-out.p_hin);

    % Pump/motor and generator
    out.power.P_pmLoss = (out.p_lin - out.p_hout).*out.q_pm ...
                       - out.Tpm.*out.w_pm;
    out.power.P_gen = -((-out.Tgen.*out.w_pm > 0).*par.eta_g ...
                    + (-out.Tgen.*out.w_pm < 0)./par.eta_g) ...
                    .*out.Tgen.*out.w_pm;
    out.power.P_genLoss = -out.Tgen.*out.w_pm - out.power.P_gen;

    % RO ripple control valve
    out.power.P_rv = out.q_rv.*(out.p_hout-out.p_ro);

    % Charge pump
    out.power.P_cElec = 1/(par.eta_c*par.eta_m)*out.q_c.*(out.p_lin-out.par.p_o);
    out.power.P_cLoss = out.power.P_cElec - out.q_c.*(out.p_lin-out.par.p_o);

    % ERU
    dp_ERUfeed = out.p_ro - out.p_lin;
    P_ERUfeed = out.q_ERUfeed.*dp_ERUfeed;
    P_ERUbrine = out.q_brine.*(out.p_ro - out.par.p_o);
    PbalERU = par.ERUconfig.present ...
                *(1/(par.eta_ERUv*par.eta_ERUm)*P_ERUfeed ...
                - par.eta_ERUv*par.eta_ERUm*P_ERUbrine);
    out.power.P_ERULoss = PbalERU + P_ERUbrine - P_ERUfeed;
    out.power.P_ERUelec = ((PbalERU > 0)./par.eta_m ...
                        + (PbalERU < 0).*par.eta_m) ...
                        .*PbalERU;
    out.power.P_ERUelecLoss = out.power.P_ERUelec - PbalERU;

    % Pressure relief valves
    out.power.P_linPRV = out.q_linPRV.*(out.p_lin-out.par.p_o);
    out.power.P_hinPRV = out.q_hinPRV.*(out.p_hin-out.p_lout);
    out.power.P_roPRV = out.q_roPRV.*(out.p_ro-out.p_lin);

    % Pipeline losses
    out.power.P_LPPL = [pLsoln_LP(:).PPfric]';
    out.power.P_HPPL = [pLsoln_HP(:).PPfric]';

    % Electrical Energy Storage
    out.power.P_battery = out.power.P_gen ...
                        - out.power.P_cElec - out.power.P_ERUelec;
    out.power.deltaE_battery = trapz(out.t,out.power.P_battery);

    %% Energy Balance
    % Change in available potential energy in WEC-driven pump chambers
    dp = 1e2;
    cap = @(p,V) deadVCap(p,V,par);
    Va = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max - theta);
    Vb = @(theta) par.V_wecDead + par.D_WEC*(par.theta_max + theta);

     % Chamber 'a'
    cap1 = @(p) cap(p,Va(out.theta(1)));
    capend = @(p) cap(p,Va(out.theta(end)));
    deltaE_a = deltaE_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaE_NI(out.par.p_o,out.p_a(1),cap1,dp);

     % Chamber 'b'
    cap1 = @(p) cap(p,Vb(out.theta(1)));
    capend = @(p) cap(p,Vb(out.theta(end)));
    deltaE_b = deltaE_NI(out.par.p_o,out.p_b(end),capend,dp) ...
                - deltaE_NI(out.par.p_o,out.p_b(1),cap1,dp);

    % Change in available potential energy in accumulators
    dp = 1e2;
     % Low-pressure pipeline inlet (charge pump outlet)
    cap = @(p) capAccum(p,par.pc_lin,par.Vc_lin,par.f,par);
    deltaE_lin = deltaE_NI(out.p_lin(1),out.p_lin(end),cap,dp);

     % Low-pressure pipeline outlet (inlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_lout,par.Vc_lout,par.f,par);
    deltaE_lout = deltaE_NI(out.p_lout(1),out.p_lout(end),cap,dp);

     % High-pressure pipeline inlet (outlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_hin,par.Vc_hin,par.f,par);
    deltaE_hin = deltaE_NI(out.p_hin(1),out.p_hin(end),cap,dp);

     % High-pressure pipeline outlet
    cap = @(p) capAccum(p,par.pc_hout,par.Vc_hout,par.f,par);
    deltaE_hout = deltaE_NI(out.p_hout(1),out.p_hout(end),cap,dp);

     % High-pressure inlet of RO module
    cap = @(p) capAccum(p,par.pc_ro,par.Vc_ro,par.f,par);
    deltaE_ro = deltaE_NI(out.p_ro(1),out.p_ro(end),cap,dp);

    % Change in available potential energy in pipelines
    dp = 1e2;
     % Low-pressure pipeline
    cap = @(p) lineCap(p,1,par);
    pPL = out.pLP;
    deltaE_PL = zeros(size(pPL(1,:)));
    for iyp = 1:numel(deltaE_PL)
        deltaE_PL(iyp) = deltaE_NI(pPL(1,iyp),pPL(end,iyp),cap,dp);
    end
    deltaPE_LPPL = sum(deltaE_PL);
     % High-pressure pipeline
    cap = @(p) lineCap(p,2,par);
    pPL = out.pHP;
    deltaE_PL = zeros(size(pPL(1,:)));
    for iyp = 1:numel(deltaE_PL)
        deltaE_PL(iyp) = deltaE_NI(pPL(1,iyp),pPL(end,iyp),cap,dp);
    end
    deltaPE_HPPL = sum(deltaE_PL);

    % Change in kinetic energy in pipelines
     % Low-pressure pipeline
    LineID = 1;
    deltaKE_LPPL = 1/2*par.I(LineID)/par.A_line(LineID)^2* ...
                    sum((signed_sq(out.qLP(end,:)) ...
                       - signed_sq(out.qLP(1,:))));
     % High-pressure pipeline
    LineID = 2;
    deltaKE_HPPL = 1/2*par.I(LineID)/par.A_line(LineID)^2* ...
                    sum((signed_sq(out.qHP(end,:)) ...
                       - signed_sq(out.qHP(1,:))));

    % Total change in stored energy in the system
    deltaE_sys = deltaE_a + deltaE_b ...
        + deltaE_lin + deltaE_lout ...
        + deltaE_hin + deltaE_hout ...
        + deltaE_ro ...
        + deltaPE_LPPL + deltaKE_LPPL ...
        + deltaPE_HPPL + deltaKE_HPPL ...
        + out.power.deltaE_battery;

    % Power flow at boundaries
    P_in = out.power.P_WEC;
    P_out = out.q_perm.*(out.p_ro - out.par.p_o);
    P_loss = out.power.P_wpLoss + out.power.P_sv ...
            + out.power.P_ain + out.power.P_aout ...
            + out.power.P_bin + out.power.P_bout ...
            + out.power.P_pmLoss + out.power.P_genLoss ...
            + out.power.P_rv ...
            + out.power.P_cLoss ...
            + out.power.P_ERULoss + out.power.P_ERUelecLoss ...
            + out.power.P_linPRV + out.power.P_hinPRV + out.power.P_roPRV ...
            + out.power.P_LPPL + out.power.P_HPPL;
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

     % Chamber 'a'
    cap1 = @(p) cap(p,Va(out.theta(1)));
    capend = @(p) cap(p,Va(out.theta(end)));
    deltaV_a = deltaV_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaV_NI(out.par.p_o,out.p_a(1),cap1,dp);

     % Chamber 'b'
    cap1 = @(p) cap(p,Vb(out.theta(1)));
    capend = @(p) cap(p,Vb(out.theta(end)));
    deltaV_b = deltaV_NI(out.par.p_o,out.p_a(end),capend,dp) ...
                - deltaV_NI(out.par.p_o,out.p_a(1),cap1,dp);

    % Change in fluid volume in accumulators
    dp = 1e2;
     % Low-pressure pipeline inlet (charge pump outlet)
    cap = @(p) capAccum(p,par.pc_lin,par.Vc_lin,par.f,par);
    deltaV_lin = deltaV_NI(out.p_lin(1),out.p_lin(end),cap,dp);

     % Low-pressure pipeline outlet (inlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_lout,par.Vc_lout,par.f,par);
    deltaV_lout = deltaV_NI(out.p_lout(1),out.p_lout(end),cap,dp);

     % High-pressure pipeline inlet (outlet of WEC-driven pump)
    cap = @(p) capAccum(p,par.pc_hin,par.Vc_hin,par.f,par);
    deltaV_hin = deltaV_NI(out.p_hin(1),out.p_hin(end),cap,dp);

     % High-pressure pipeline outlet
    cap = @(p) capAccum(p,par.pc_hout,par.Vc_hout,par.f,par);
    deltaV_hout = deltaV_NI(out.p_hout(1),out.p_hout(end),cap,dp);

     % High-pressure inlet of RO module
    cap = @(p) capAccum(p,par.pc_ro,par.Vc_ro,par.f,par);
    deltaV_ro = deltaV_NI(out.p_ro(1),out.p_ro(end),cap,dp);

    % Equivalent change in fluid volume in pipelines
    dp = 1e2;
     % Low-pressure pipeline
    cap = @(p) lineCap(p,1,par);
    pPL = out.pLP;
    deltaV_PL = zeros(size(pPL(1,:)));
    for iyp = 1:size(pPL,2)
        deltaV_PL(iyp) = deltaV_NI(pPL(1,iyp),pPL(end,iyp),cap,dp);
    end
    deltaV_LPPL = sum(deltaV_PL);
     % High-pressure pipeline
    cap = @(p) lineCap(p,2,par);
    pPL = out.pHP;
    deltaV_PL = zeros(size(pPL(1,:)));
    for iyp = 1:size(pPL,2)
        deltaV_PL(iyp) = deltaV_NI(pPL(1,iyp),pPL(end,iyp),cap,dp);
    end
    deltaV_HPPL = sum(deltaV_PL);

    % Total change in stored volume in the system
    out.deltaV_total = deltaV_a + deltaV_b ...
                     + deltaV_lin + deltaV_lout ...
                     + deltaV_hin + deltaV_hout ...
                     + deltaV_ro ...
                     + deltaV_LPPL + deltaV_HPPL;

    % Flow at boundaries
    qbnds = out.q_c ...
            - (out.q_perm + out.q_brine + out.q_linPRV);

    % Total balance of flow: flow in minus change in volume stored
    out.Vbal = trapz(out.t,qbnds) - out.deltaV_total;
    out.Vbal_error = out.Vbal/trapz(out.t,out.q_perm);

%% %%%%%%%%%%%%   FUNCTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = sys(t,y,par)
    	[dydt, ~ , ~ ] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sysPost(t,y,par)
    	[dydt, nonState, control] = sys_prototype(t,y,par);
    end

    function [dydt, nonState, control] = sys_prototype(t,y,par)
        [dydt, nonState, control] = sys_parPTO(t,y,par);
    end
end
