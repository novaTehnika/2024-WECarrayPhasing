% parSafeStartSlurm.m script m-file
% AUTHORS:
% Jeremy Simmons (email: simmo536@umn.edu)
% University of Minnesota
% Department of Mechanical Engineering
%
% CREATION DATE:
% 11/30/2023
%
% PURPOSE/DESCRIPTION:
% This script executes creates and assigns a folder specific to the SLURM
% job being executed. THis avoids issues of multiple instances  of MATLAB
% creating identically named jobs.
%
% This code was adapted from that given at:
% https://www2.physics.ox.ac.uk/it-services/ ...
% running-multiple-matlab-jobs-at-the-same-time (30 Nov 2023)
%
% FILE DEPENDENCY:
%
% UPDATES:
% 11/30/2023 - Created
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

cl = parcluster(); %create a cluster object

slurmid = getenv('SLURM_JOB_ID'); %get the unique job ID as seen by Hydra

storage_folder = strcat('.',filesep,'matlabtemp',filesep,slurmid);

mkdir(storage_folder); %make a temporary folder for this session

cl.JobStorageLocation = storage_folder; %and tell the cluster to use it

parpool(cl,nWorkers); %create the parallel pool for this cluster