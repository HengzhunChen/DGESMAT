function beta = vecnorm(X, varargin)
% SPINOR/VECNORM Vecnorm function for wave function
%
%    beta = vecnorm(X, options) returns the vector norm of wave functions
%    according to options.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

beta = vecnorm(X.wavefun, varargin{:});

end