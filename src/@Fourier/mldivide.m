function x = mldivide(F, y)
% FOURIER/MLDIVIDE Operator overload performing inverse Fourier transform
%
%    Usage: x = F \ y;
%    y can be numerical data or Spinor object.
%
%    See also Fourier, Fourier/mtimes, Fourier/ctranspose.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

x = F' * y;

end