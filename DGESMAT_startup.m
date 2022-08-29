function DGESMAT_startup()
% DGESMAT_STARTUP  Startup file for DGESMAT
%   MAKE adds paths of the DGESMAT to Matlab.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


file_path = mfilename('fullpath');
tmp = strfind(file_path, 'DGESMAT_startup');
file_path = file_path(1 : (tmp(end)-1));

% Folder for all source files recursively
addpath(genpath([file_path 'src']));

% Folder for all external files recursively
% addpath(genpath([file_path 'external']));

% Folder for all pseudopotential files
addpath([file_path 'ppdata/default']);


end