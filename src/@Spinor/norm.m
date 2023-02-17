function beta = norm(X, varargin)
% SPINOR/NORM Norm function for wave function
%
%    beta = norm(X, options) returns the norm of the wave functions.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

beta = norm(X.wavefun, varargin{:});

end