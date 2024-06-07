% paretoFront2D.m function m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 7/19/2023
%
% PURPOSE/DESCRIPTION:
% This script serves finds indicies for the individuals on a 2D pareto 
% front. The objective data must be in a singel dimensional array. This
% function assumes the objective is to minimize the objective values.
%
% FILE DEPENDENCY:
%
% UPDATES:
% 7/19/2023 - Created from analysis in postprocess_datafiles_accum_wRV.m
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
function non_dominated = paretoFront2D(obj1,obj1_minORmax,...
                                       obj2,obj2_minORmax)
    % Check that dimensions of the objective data match
    sz = size(obj1);
    if ~(sz(1) == 1 || sz(2) == 1)
        error('The objectives are not 1D.')
    end

    % Check that the objective data are one-dimensional
    if ~(sum(size(obj1) == size(obj2)) == 2  ...
            || sum(size(obj1') == size(obj2)) == 2)
        error('The dimensions of objectives do not match.')
    end
    
    % Check that whether each objective should be minimized or maximized
    errorMinMax = @() error("Specify 'min' or 'max' for each objective.");
    if contains(obj1_minORmax,'min','IgnoreCase',true)
        obj1 = -obj1;
    elseif contains(obj1_minORmax,'max','IgnoreCase',true)
        % do nothing
    else
        errorMinMax()
    end
    if contains(obj2_minORmax,'min','IgnoreCase',true)
        obj2 = -obj2;
    elseif contains(obj2_minORmax,'max','IgnoreCase',true)
        % do nothing
    else
        errorMinMax()
    end

    % Number of individuals
    n = numel(obj1);

    % Initialize dominance matrix
    dominance = false(n);
    
    % Find and mark dominating individuals
    for i = 1:n
        for j = 1:n
            if obj1(i) >= obj1(j) && obj2(i) >= obj2(j) ...
               && (obj1(i) > obj1(j) || obj2(i) > obj2(j)) ...
               || isnan(obj2(j))
                dominance(i,j) = true;
            end
        end
    end
    
    % Identify non-dominated individuals
    non_dominated = find(sum(dominance, 1) == 0);
end