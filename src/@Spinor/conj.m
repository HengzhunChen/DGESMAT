function Xconj = conj(X)
% SPINOR/CONJ Conjugate function for wave function.
%
%    Xconj = conj(X) returns a Spinor object with the conjugate value of 
%    the wave function of X.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Xconj = X;
Xconj.wavefun = conj(X.wavefun);

end