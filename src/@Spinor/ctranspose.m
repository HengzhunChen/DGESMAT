function Y = ctranspose(X)
% SPINOR/CTRANPOSE Operator overload. Ctranspose function for wavefunction.
%
%    Y = ctranspose(X) returns a Spinor object with the conjugate 
%    transpose value of the wave function of X.
%
%    See aslo Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Y = X;
Y.trans = ( X.trans == 0 );

end