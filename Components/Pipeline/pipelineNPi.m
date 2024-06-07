function [q_in, q_out, dydt_pLine, pLsoln] = ...
         pipelineNPi(t,y,p_in,p_out,par,lineID,postProcess)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pipelineNPi.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/25/2021
%
% PURPOSE/DESCRIPTION:
% 
%
% FILE DEPENDENCY: NA
%
% UPDATES:
% 3/25/2021 - created from implimentation in simPTO_V01x05.m
% 11/2/2021 - added lineID specification for par.n_seg for "load arrays for
% easier indexing" (fixing a mistake). Added lineID to parameters
% 5/3/2023 - lineCap() and flowR() moved to separate function m-files
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

    % load arrays for easier indexing
    q = y(1:2:2*par.n_seg(lineID)-1);
    p = [y(2:2:2*par.n_seg(lineID)-1); p_out];

    % output boundary conditions to system
    q_in = q(1);
    q_out = q(end);

    % intialize dy/dt vector
    dydt_pLine = zeros(2*par.n_seg(lineID)-1,1); 

    % dq/dt at first pi-lump
    dydt_pLine(1) = (p_in - p(1) - flowR(q(1),lineID,par)*q(1))...
                /par.I(lineID);

    % dp/dt and dq/dt at remaining pi-lumps
    for j = 1:par.n_seg(lineID)-1 
        % dp/dt
        dydt_pLine(2*j) = (q(j) - q(j+1)) / lineCap(p(j),lineID,par); 
        % dq/dt
        dydt_pLine(2*j+1) = (p(j) - p(j+1) ...
                        - flowR(q(j+1),lineID,par)*q(j+1))...
                        /par.I(lineID); 
    end

    % Calculate friction loss
    if postProcess
        
        % intialize sketch variable
        PPfric = zeros(par.n_seg(lineID),1);
        
        % calc. friction loss in each segment
        for i = 1:par.n_seg(lineID)
            PPfric(i) = flowR(q(i),lineID,par)*q(i)^2;
        end

        % calc. total friction loss
        pLsoln.PPfric = sum(PPfric);
        
    else
        pLsoln = [];
        
    end

end