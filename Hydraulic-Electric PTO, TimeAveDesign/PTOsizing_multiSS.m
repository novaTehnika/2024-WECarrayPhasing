function data = PTOsizing_multiSS(D_wArray,D_mArray,bounds,iPTO, ...
                                  design_case,par)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PTOsizing_multiSS.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/3/2024
%
% PURPOSE:
% The purpose of this script is to perform a design study sizing a 
% hydraulic-electric PTO architectures for wave-powered electric power
% production operating in a known distibution of sea conditions.
% 
% The model is a simple, static model with that includes two-way coupling
% with the time-averaged simulation results of a WEC; the coupling is set
% up such that the reaction force from the PTO is a function of the average
% WEC speed (or power absortion) and the average WEC speed (or power
% absorption) is function of the reaction torque from the PTO.
%
% In the design study, the WEC-driven pump displacement and hydraulic motor
% displacement are varied using a grid of values. For each set of these
% design variables, an optimization is performed to select the nominal 
% operating pressure of the system and the switching duty (if applicable).
% This routine is performed for each sea state.
% 
% The optimization of the operating pressure and switching duty is a
% nonlinear, constrained optimization which seeks to maximize the electric
% power production subject to the following constraints:
% 1) the operating pressure at the RO module is within a prescribe range
% 2) the minimum speed of the hydrualic motor is maintained
% If the design does not meet these constraints, the value of zero is
% recorded for electric power production rate for that design and sea
% condition.
% 
% Once the results are obtained for the specified grid of design
% parameters, design parameters are selected for each sea state to give 
% optimal combinations based on four possible configurations:
% A) the WEC-driven pump and hydraulic motor displacement are fixed across
% all sea conditions.
% B) the WEC-driven pump displacement is variable across all sea states 
% while the hydraulic motor displacement is fixed.
% C) the hydraulic motor displacement is variable across all sea states
% while the WEC-driven pump displacement is fixed.
% D) Both the WEC-driven pump and hydraulic motor displacement are variable
% across all sea states
% The design combinations are organized based on the largest allowable
% value for the two design parameters; the design parameters are selected 
% for each sea state from the set of designs that have design parameters 
% values less than or equal to the maximum allowable value, when variable
% by sea state. Otherwise, when the design parameter is fixed across the
% sea states, the set includes only designs that have the same value of the
% fixed design parameter. The design with greatest electric power
% production is selected from the set of available designs in each sea
% state.
%
%
% FILE DEPENDENCY: 
% model_timeAve_hydElecPTO.m
% maxRC.m
% zero2nan.m
% 
% UPDATES:
% 7/3/2024 - adapted from PTOsizing_multiSS.m in 
% 2021-TimeAvePTOarchetectureStudy repository.
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

    %% get size of design space
    nD_w = length(D_wArray);
    nD_m = length(D_mArray);

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
    
    %% Set bounds and intial guesses for operational variable to be
    % optimized
    
    % Scaling factors
    pScale = 1e7; % [Pa] Scaling factor for pressure variables
    dutyScale = 1; % [-] Scaling factor for switching duty
    
    % Bounds
    p_h_bnds = bounds.p_h_bnds/pScale; % [Pa/Pa] Bounds for system pressure
    D_bnds = bounds.D_bnds/dutyScale; % [-] bounds for valve switching duty
    
    par.p_h_bnds = p_h_bnds*pScale; % [Pa] Bounds on feed pressure to be 
                                    % passed to contraint fnc nonlcon()
    
    % Inital values
    p_h_o = (10e6)/pScale;
    duty_o = (0.9)/dutyScale;
    
    %% Build arguments for fmincon()
    switch iPTO
        case 1
            lb = p_h_bnds(1);
            ub = p_h_bnds(2);
            x0 = p_h_o;
            x_scale = pScale;
        case 2
            lb = [p_h_bnds(1); D_bnds(1)];
            ub = [p_h_bnds(2); D_bnds(2)];
            x0 = [p_h_o; duty_o];
            x_scale = [pScale; dutyScale];
    end
    
    % Linear constraints
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    
    %% Initialize arrays for record of design study results
    D_w = zeros(nD_w, nD_m); % [m^3/rad] WEC-driven pump displacement
    D_m = zeros(nD_w, nD_m); % [m^3/Pa.s] Motor displacement
    T_c = zeros(nD_w, nD_m, nSS); % [Nm] Torque applied to WEC by PTO
    PP_w = zeros(nD_w, nD_m, nSS); % [W] Power transmitted to PTO by WEC
    PP_gen = zeros(nD_w, nD_m, nSS); % [W] Elec. power generated
    PP_PRV = zeros(nD_w, nD_m, nSS); % [W] Power lost to PRV
    c = zeros(nD_w, nD_m, nSS,4); % Constraint variables (c < 0)
    feasible = zeros(nD_w,nD_m,nSS);
    p_h = zeros(nD_w, nD_m, nSS); % [Pa] Operating pressure
    duty = zeros(nD_w, nD_m, nSS); % [-] switching valve duty
    
    %% Execute grid study with embedded optimization routine for
    % operational  parameters to obtain the optimal performance of each set
    % of design parameters.
    
    % loop through WEC-pump displacement
    parfor iD_w = 1:nD_w
        % pass information to par-for loop
        param = par;
        D_mArray_par = D_mArray;

        % Set WEC-driven pump displacement
        param.D_w = D_wArray(iD_w);
        
        % loop through motor displacement
        for iD_m = 1:nD_m
            
            % Set motor displacement
            param.motor.D = D_mArray_par(iD_m);
            
            % loop through sea states
            for iSS = 1:nSS
                
                % Set sea state
                param.SS = par.SSset(iSS);
                
                % Record design variables for the current case
                D_w(iD_w,iD_m) = D_wArray(iD_w);
                D_m(iD_w,iD_m) = D_mArray(iD_m);
                
                % Build the objective function for the current case
                obj = ...
                  @(x) 1/model_timeAve_hydElecPTO(x.*x_scale,param,iPTO,1);
                
                % Build the contraint function for the current case
                nonlcon = ...
                  @(x) model_timeAve_hydElecPTO(x.*x_scale,param,iPTO,2);
                
                % Execute optimization
                x = fmincon(obj,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
                
                % Record optimization result for operational parameters
                p_h(iD_w,iD_m,iSS) = x(1)*pScale;
                if iPTO == 2
                    duty(iD_w,iD_m,iSS) = x(2)*dutyScale;
                end
                
                % Calculate and record contraint variables and record
                % feasiblity
                [temp,~] = nonlcon(x);
                c(iD_w,iD_m,iSS,:) = temp;
                feasible(iD_w,iD_m,iSS) = temp(1) <= 0 ...
                                  &  temp(2) <= 0 ...
                                  &  temp(3) <= 0 ...
                                  &  temp(4) <= 0;
    
                % Calculate and record performance variables
                [T_c(iD_w,iD_m,iSS),...
                 PP_w(iD_w,iD_m,iSS),...
                 PP_gen(iD_w,iD_m,iSS),...
                 PP_PRV(iD_w,iD_m,iSS)] = ...
                 model_timeAve_hydElecPTO(x.*x_scale,param,iPTO,3);
    
                % modify electric power production with weight for sea
                % state and feasibility of result. Weight is given as
                % percentage.
                PP_gen(iD_w,iD_m,iSS) = feasible(iD_w,iD_m,iSS)...
                                            *par.weight(iSS)/100 ...
                                            *PP_gen(iD_w,iD_m,iSS);
    
            end
        end
    end
    
    clearvars temp
    
    %% Post-Process
    % combine designs across sea states based on the four possible ways in
    % which the system may have adjustability between sea states (whether a
    % design parameter is fixed or not) with the design combination
    % organized based on the maximum allowable value for each design
    % variable

    % Initialize structure array recording the set of design combinations
    clearvars design
    design.PP_gen = [];
    design.iD_w = [];
    design.iD_m = [];
    
    % Initialize arrays recording contribution to total production from
    % selection for each sea state and the combined total production
    PP_genTotal = zeros(nD_w,nD_m);
    
    % Perform selection for each sea state and design case
    
    % loop though max displacement
    for iD_w = nD_w:-1:1
        % loop though max motor displacement
        for iD_m = nD_m:-1:1
            % loop through sea states
            for iSS = nSS:-1:1
                % Mask production results for the current sea state and
                % design parameter combinations based on whether the the
                % design combination fits the constraints of the design
                % case for the current combination of maximum values for
                % the design parameters. Then, find the maximum production
                % value and the coresponding (index to the) design
                % parameters        

                switch design_case
                    case 1 
                        % Option A: fixed displacment pump
                        %           fixed displacment motor
                        temp.PP_gen = PP_gen(iD_w,iD_m,iSS);
                        temp.iD_w = iD_w;
                        temp.iD_m = iD_m;
    
                    case 2 
                        % Option B: variable displacement pump
                        %           fixed displacment motor
                        [temp.PP_gen, temp.iD_w, ~] = ...
                                        maxRC(PP_gen(:,:,iSS) ...
                                        .* (D_w(:,:) <= D_wArray(iD_w) ...
                                        & D_m(:,:) == D_mArray(iD_m)));
                        temp.iD_m = iD_m;

                   case 3 
                       % Option C: fixed displacement pump
                       %           variable displacment motor
                        [temp.PP_gen, ~, temp.iD_m] = ...
                                        maxRC(PP_gen(:,:,iSS) ...
                                        .* (D_w(:,:) == D_wArray(iD_w) ...
                                         & D_m(:,:) <= D_mArray(iD_m)));
                        temp.iD_w = iD_w;
    
                    case 4 
                        % Option D: variable displacement pump and motor
                        [temp.PP_gen, temp.iD_w, temp.iD_m] = ...
                                        maxRC(PP_gen(:,:,iSS) ...
                                        .* (D_w(:,:) <= D_wArray(iD_w) ...
                                        & D_m(:,:) <= D_mArray(iD_m)));

                end
                design(iD_w,iD_m,iSS).D_w = ...
                                        D_wArray(temp.iD_w);
                design(iD_w,iD_m,iSS).D_m = ...
                                        D_mArray(temp.iD_m);
                design(iD_w,iD_m,iSS).p_i = ...
                                        p_h(temp.iD_w,temp.iD_m,iSS);
                design(iD_w,iD_m,iSS).duty = ...
                                        duty(temp.iD_w,temp.iD_m,iSS);
                design(iD_w,iD_m,iSS).PP_w = ...
                                        PP_w(temp.iD_w,temp.iD_m,iSS);
                design(iD_w,iD_m,iSS).PP_gen = ...
                                        PP_gen(temp.iD_w,temp.iD_m,iSS);
                design(iD_w,iD_m,iSS).PP_PRV = ...
                                        PP_PRV(temp.iD_w,temp.iD_m,iSS);
                design(iD_w,iD_m,iSS).feasible = ...
                                        feasible(temp.iD_w,temp.iD_m,iSS);

            end
            % Sum contributions to total production results
            PP_genTotal(iD_w,iD_m) = sum([design(iD_w,iD_m,:).PP_gen]);
    
        end
    end
    
    %% Package results
     % Results given pump and motor displacement
    data.par = par;
    data.D_w = D_w; % given displacement
    data.D_m = D_m; % given membrane area
    data.p_h = p_h;
    data.duty = duty;
    data.feasible = feasible;
    data.PP_w = PP_w;
    data.PP_gen = PP_gen;
    data.PP_PRV = PP_PRV;

     % Results of optimal performance for the given PTO architecture
    data.PP_genTotal = PP_genTotal;
    data.design = design;

end