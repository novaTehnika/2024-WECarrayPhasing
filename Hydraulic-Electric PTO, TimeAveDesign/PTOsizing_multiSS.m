function data = PTOsizing_multiSS(D_wArray,S_roArray,bounds,iPTO, ...
                                  design_case,ERUconfig,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTOsizing_multiSS.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 12/31/2021
%
% PURPOSE:
% The purpose of this script is to perform a design study of PTO
% architectures for wave-powered reverse osmosis operating in a known 
% distibution of sea conditions. The script is set up to
% perform the same study for four different PTO architectures. 
% 
% The models used in the design study are simple, static models with that 
% include two-way coupling with the time-averaged simulation results of a
% WEC; the coupling is set up such that the reaction force from the PTO is   
% a function of the average WEC speed (or power absortion) and the average 
% WEC speed (or power absorption) is function of the reaction torque from 
% the PTO.
%
% In the design study, the WEC-driven pump displacement and RO module
% permeabilty coefficient are varied using a grid of values. For each set
% of these design variables, an optimization is performed to select the
% nominal operating pressure of the system (either the RO feed pressure or 
% the pressure at the outlet of the WEC_driven pump depending on the PTO 
% architecture) and the switching duty (if applicable). This routine is
% performed for each sea state.
% 
% The optimization of the operating pressure and switching duty is a
% nonlinear, constrained optimization which seeks to maximize the permeate
% production subject to the following constraints:
% 1) the electrical power production meets or exceeds the electrical demand
% 2) the operating pressure at the RO module is within a prescribe range
% If the design does not meet these constraints, the value of zero is
% recorded as the permeate production rate for that design.
% 
% Once the results are obtained for the specified grid of design
% parameters, design parameters are selected for each sea state to give 
% optimal combinations based on four possible configurations:
% A) the WEC-driven pump and RO module size are fixed across all sea 
% conditions
% B) the WEC-driven pump is variable across all sea states while the active 
% RO module size is fixed
% C) the active RO module size is variable across all sea states while the
% WEC-driven pump displacement is fixed.
% D) Both the WEC-driven pump and active RO module size are variable across
% all sea states
% The design combinations are organized based on the largest allowable
% value for the two design parameters; the design parameters are selected 
% for each sea state from the set of designs that have design parameters 
% values are less than or equal to the maximum allowable value, when 
% variable by sea state. Otherwise, when the design parameter is fixed
% across the sea states, the set includes only designs that have the same 
% value of the fixed design parameter. The design with greatest permeate
% production is selected from the set of available designs in each sea
% state.
%
%
% FILE DEPENDENCY: 
% maxRC.m
% zero2nan.m
% model_timeAvePTO.m
% 
% UPDATES:
% 12/31/2021 - created.
% 08/22/2022 - Corrected constraint on WEC-driven torque to be included in 
% constraint function in the optimization.
% 06/15/2023 - add optional ERU (ERUconfig=0 -> w/o ERU; ERUconfig=1 ->
% w/ ERU). Values between 0 and 1 effectively set an efficiency of the ERU.
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

    %% get size of design space
    nD_w = length(D_wArray);
    nS_ro = length(S_roArray);

    %% Options for optimization routine fmincon()
    % options = optimset('Display','iter','PlotFcns',@optimplotfval);
    options = optimoptions('fmincon');
    options.Display = 'off';
    options.MaxFunEvals = 1e4;
    options.MaxIter = 1e2;
%     options.TolFun = 1e-17;
%     options.TolCon = 1e-17;
    options.OptimalityTolerance = 1e-20;
    options.StepTolerance = 1e-10;

    nSS = length(par.SSset);
    
    %% Set bounds and intial guesses for operational variable to be optimized
    
    % Scaling factors
    pScale = 1e7; % [Pa] Scaling factor for pressure variables
    dutyScale = 1; % [-] Scaling factor for switching duty
    
    % Bounds
    p_f_bnds = bounds.p_f_bnds/pScale; % [Pa/Pa] Bounds for feed pressure
    p_w_bnds = bounds.p_w_bnds/pScale; % [Pa/Pa] Bounds for pressure at WEC driven pump
    D_bnds = bounds.D_bnds/dutyScale; % [-] bounds for valve switching duty
    
    par.p_f_bnds = p_f_bnds*pScale; % [Pa] Bounds on feed pressure to be 
                                     % passed to contraint function nonlcon()
    
    % Inital values
    p_f_o = (6e6)/pScale;
    p_w_o = (8e6)/pScale;
    duty_o = (0.9)/dutyScale;
    
    %% Build arguments for fmincon()
    switch iPTO
        case 1
            lb = p_f_bnds(1);
            ub = p_f_bnds(2);
            x0 = p_f_o;
            x_scale = pScale;
        case 2
            lb = [p_f_bnds(1); D_bnds(1)];
            ub = [p_f_bnds(2); D_bnds(2)];
            x0 = [p_f_o; duty_o];
            x_scale = [pScale; dutyScale];
        case 3
            lb = p_w_bnds(1);
            ub = p_w_bnds(2);
            x0 = p_w_o;
            x_scale = pScale;
        case 4
            lb = [p_w_bnds(1); D_bnds(1)];
            ub = [p_w_bnds(2); D_bnds(2)];
            x0 = [p_w_o; duty_o];
            x_scale = [pScale; dutyScale];
    end
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    %% Initialize arrays for record of design study results
    q_perm = zeros(nD_w, nS_ro, nSS); % [m^3/s] mean permeate flow rate
    D_w = zeros(nD_w, nS_ro); % [m^3/rad] WEC-driven pump displacement
    S_ro = zeros(nD_w, nS_ro); % [m^3/Pa.s] Permeate coefficient
    T_c = zeros(nD_w, nS_ro, nSS); % [Nm] Torque applied to WEC by PTO
    PP_w = zeros(nD_w, nS_ro, nSS); % [W] Power transmitted to PTO by WEC
    PP_gen = zeros(nD_w, nS_ro, nSS); % [W] Elec. power generated
    PP_c = zeros(nD_w, nS_ro, nSS); % [W] Elec. power consumed by charge pump
    c = zeros(nD_w, nS_ro, nSS,4); % Constraint variables (c < 0)
    feasible = zeros(nD_w,nS_ro,nSS);
    p_i = zeros(nD_w, nS_ro, nSS); % [Pa] Operating pressure
    duty = zeros(nD_w, nS_ro, nSS); % [-] switching valve duty
    
    %% Execute grid study with embedded optimization routine for operational 
    % parameters to obtain the optimal performance of each set of design
    % parameters.
    
    % loop through WEC-pump displacement
    parfor iD_w = 1:nD_w
        param = par; % pass parameter structure to par-for loop
        
        % Set WEC-driven pump displacement
        param.D_w = D_wArray(iD_w);
        S_roArray_par = S_roArray;
        
        % loop through RO module size
        for iS_ro = 1:nS_ro
            
            % Set RO membrane area
            param.S_ro = S_roArray_par(iS_ro);
            
            % loop through sea states
            for iSS = 1:nSS
                
                % Set sea state
                param.SS = par.SSset(iSS);
                
                % Record design variables for the current case
                D_w(iD_w,iS_ro) = D_wArray(iD_w);
                S_ro(iD_w,iS_ro) = S_roArray(iS_ro);
                
                % Build the objective function for the current case
                obj = @(x) 1/model_timeAvePTO(x.*x_scale,param,iPTO,ERUconfig,1);
                
                % Build the contraint function for the current case
                nonlcon = @(x) model_timeAvePTO(x.*x_scale,param,iPTO,ERUconfig,2);
                
                % Execute optimization
                x = fmincon(obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
                
                % Record optimization result for operational parameters
                p_i(iD_w,iS_ro,iSS) = x(1)*pScale;
                if iPTO == 2 || iPTO == 4
                    duty(iD_w,iS_ro,iSS) = x(2)*dutyScale;
                end
                
                % Calculate and record contraint variables and record
                % feasiblity
                [temp,~] = nonlcon(x);
                c(iD_w,iS_ro,iSS,:) = temp;
                feasible(iD_w,iS_ro,iSS) = temp(1) <= 0 ...
                                  &  temp(2) <= 0 ...
                                  &  temp(3) <= 0 ...
                                  &  temp(4) <= 0;
    
                % Calculate and record performance variables
                [q_perm(iD_w,iS_ro,iSS),...
                 T_c(iD_w,iS_ro,iSS),...
                 PP_w(iD_w,iS_ro,iSS),...
                 PP_gen(iD_w,iS_ro,iSS),...
                 PP_c(iD_w,iS_ro,iSS)] = model_timeAvePTO(x.*x_scale,param,iPTO,ERUconfig,3);
    
                % modify permeate production rate with weight for sea state and
                % feasibility of result. Weight is given as percentage.
                q_perm(iD_w,iS_ro,iSS) = feasible(iD_w,iS_ro,iSS)...
                                            *par.weight(iSS)/100 ...
                                            *q_perm(iD_w,iS_ro,iSS);
    
            end
        end
    end
    
    clearvars temp
    
    %% Post-Process
    % combine designs across sea states based on the four possible ways in
    % which the system may have adjustability between sea states (whether a 
    % design parameter is fixed or not) with the desing combination organized 
    % based on the maximum allowable value for each design variable

    % Initialize structure array recording the set of design combinations
    clearvars design
    design.q_perm = [];
    design.iD_w = [];
    design.iS_ro = [];
    
    % Initialize arrays recording contribution to total production from 
    % selection for each sea state and the combined total production
    q_permTotal = zeros(nD_w,nS_ro);
    
    % Perform selection for each sea state and design case
    
    % loop though max displacement
    for iD_w = nD_w:-1:1
        % loop though max RO size
        for iS_ro = nS_ro:-1:1
            % loop through sea states
            for iSS = nSS:-1:1
                % Mask production results for the current sea state and design 
                % parameter combinations based on whether the the design 
                % combination fits the constraints of the design case for the
                % current combination of maximum values for the design
                % parameters. Then, find the maximum production value and the 
                % coresponding (index to the) design parameters        

                switch design_case
                    case 1 
                        % Option A: fixed displacment pump, fixed RO size
                        temp.q_perm = q_perm(iD_w,iS_ro,iSS);
                        temp.iD_w = iD_w;
                        temp.iS_ro = iS_ro;
    
                    case 2 
                        % Option B: variable displacement pump, fixed RO
                        [temp.q_perm, temp.iD_w, ~] = ...
                                        maxRC(q_perm(:,:,iSS) ...
                                        .* (D_w(:,:) <= D_wArray(iD_w) ...
                                        & S_ro(:,:) == S_roArray(iS_ro)));
                        temp.iS_ro = iS_ro;

                    case 3 
                        % Option C: fixed displacement pump and variable RO
                        [temp.q_perm, ~, temp.iS_ro] = ...
                                        maxRC(q_perm(:,:,iSS) ...
                                        .* (D_w(:,:) == D_wArray(iD_w) ...
                                         & S_ro(:,:) <= S_roArray(iS_ro)));
                        temp.iD_w = iD_w;
    
                    case 4 
                        % Option D: variable displacement pump and RO size
                        [temp.q_perm, temp.iD_w, temp.iS_ro] = ...
                                        maxRC(q_perm(:,:,iSS) ...
                                        .* (D_w(:,:) <= D_wArray(iD_w) ...
                                        & S_ro(:,:) <= S_roArray(iS_ro)));
                end
                design(iD_w,iS_ro,iSS).q_perm = temp.q_perm;
                design(iD_w,iS_ro,iSS).D_w = ...
                                        D_wArray(temp.iD_w);
                design(iD_w,iS_ro,iSS).S_ro = ...
                                        S_roArray(temp.iS_ro);
                design(iD_w,iS_ro,iSS).p_i = ...
                                        p_i(temp.iD_w,temp.iS_ro,iSS);
                design(iD_w,iS_ro,iSS).duty = ...
                                        duty(temp.iD_w,temp.iS_ro,iSS);
                design(iD_w,iS_ro,iSS).PP_w = ...
                                        PP_w(temp.iD_w,temp.iS_ro,iSS);
                design(iD_w,iS_ro,iSS).PP_gen = ...
                                        PP_gen(temp.iD_w,temp.iS_ro,iSS);
                design(iD_w,iS_ro,iSS).PP_c = ...
                                        PP_c(temp.iD_w,temp.iS_ro,iSS);
                design(iD_w,iS_ro,iSS).feasible = ...
                                        feasible(temp.iD_w,temp.iS_ro,iSS);

            end
            % Sum contributions to total production results
            q_permTotal(iD_w,iS_ro) = sum([design(iD_w,iS_ro,:).q_perm]);
    
        end
    end
    
    %% Package results
     % Results given displacement and membrane area
    data.par = par;
    data.D_w = D_w; % given displacement
    data.S_ro = S_ro; % given membrane area
    data.p_i = p_i;
    data.duty = duty;
    data.feasible = feasible;
    data.PP_w = PP_w;
    data.PP_gen = PP_gen;
    data.PP_c = PP_c;

     % Results of optimal performance for the given PTO architecture
    data.q_permTotal = q_permTotal;
    data.design = design;

end