function precResidual = KerkerPrecond(scf, residual)
% SCF/KERKERPRECOND computes the KerKer precondition used in mixing. 
%
%    NOTE:
%    From the point of view of the elliptic preconditioner
%
%    (-\Delta + 4 * pi * b) r_p = - \Delta r
%
%    The Kerker preconditioner in Fourier space is
%
%    k^2 / (k^2 + 4 * pi * b)
%
%    See also SCF, SCF/AndersonMix.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% Here we choose KerkerB to be a fixed number
% FIXME hard coded
KerkerB = 0.08;

F = scf.eigSol.fft;
gkkFine = F.gkkFine;

y = F * residual;
idxnz = gkkFine ~= 0;
y(idxnz) = y(idxnz) .* ...
    gkkFine(idxnz) ./ (gkkFine(idxnz) + 4 * pi * KerkerB);
x = F' * y;

precResidual = real(x);

end