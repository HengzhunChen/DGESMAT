function precResidual = KerkerPrecond(scfDG, residual)
% SCFDG/KERKERPRECOND computes the KerKer precondition used in mixing. 
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
%    See also SCFDG, SCFDG/AndersonMix.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


F = scfDG.fft;

numElem = scfDG.numElem;
numGridFine = scfDG.domain.numGridFine;

% NOTE Fixed KerkerB parameter
% FIXME hard code
KerkerB = 0.08;

% convert residual from 3-dim array into 1-dim vector over global domain
tempVec = ElemVecToGlobal(residual, numGridFine, numElem);

Y = F * tempVec;
% Do not touch the zero frequency 
% Procedure taken from VASP
gkkFine = F.gkkFine;
idxnz = gkkFine ~= 0;
Y(idxnz) = Y(idxnz) .* gkkFine(idxnz) ./ (gkkFine(idxnz) + 4 * pi * KerkerB);
tempVec = real(F' * Y);

% convert tempVec from 1-dim vector into 3-dim array over global domain
precResidual = GlobalVecToElem(tempVec, numGridFine, numElem);


end