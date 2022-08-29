function beta = sum(X, varargin)
% SPINOR/SUM Sum function for wave function
%
%    beta = sum(X) returns the sum of the wave functions.
%
%    beta = sum(X, dim) returns the sum of the wave functions along the 
%    dim dimension.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

beta = sum(X.wavefun, varargin{:});

end