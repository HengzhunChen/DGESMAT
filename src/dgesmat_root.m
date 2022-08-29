function path = dgesmat_root()
% DGESMAT_ROOT Root directory of DGESMAT.
%
%   S = dgesmat_root returns a string that is the absolute path of the
%   directory where the DGESMAT software is installed.
%  
%   See also matlabroot.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

path = mfilename('fullpath');
path = path(1 : end-16);

end