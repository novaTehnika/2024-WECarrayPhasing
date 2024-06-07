function dist = statsTimeVar_cdf(t,x)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statsTimeVar_cdf.m function m-file
% AUTHORS: 
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 3/29/2021
%
% PURPOSE/DESCRIPTION:
% Estimate the cumulative density function using kernel density estimation 
% for a variable coming from a variable time step simulation. Then
% calculate the  maximum value (from the original data). 
%
% FILE DEPENDENCY:
% ksdensity.m - Built-in MATLAB function
%
% UPDATES:
% 3/29/2021 - created 
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

    dt = [0; diff(t(:))];
    weight_varTime = dt/t(end);
    
    [dist.f, dist.xi] = ksdensity(x,...
                    'Weights',weight_varTime,...
                    'NumPoints',100,...
                    'Function','cdf',... %'Support','positive',...
                    'BoundaryCorrection','reflection');  
                
%     dxi = dist.xi(2) - dist.xi(1);
%     dist.mu = sum(dist.xi.*dist.f)*dxi;
%     dist.var = sum(dist.xi.^2.*dist.f)*dxi - dist.mu^2;
%     dist.sigma = sqrt(dist.var);
    dist.max = max(x);

end
