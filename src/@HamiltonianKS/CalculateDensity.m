function [totalCharge, density] = CalculateDensity(HamKS, psi, occrate)
% HAMILTONIANKS/CALCULATEDENSITY calculate density of HamKS.
% 
%    [totalCharge, density] = CalculateDensity(HamKS, psi, occrate)
%    computes density over grid and its sum totalCharge from wavefun psi
%    and occupation rate occrate.
%
%    See also HamiltonianKS, Spinor, SCF/CalculateOccupationRate.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% only restricted spin case
numSpin = HamKS.numSpin;
vol = HamKS.domain.Volume();
ntotFine = HamKS.domain.NumGridTotalFine();
numOcc = HamKS.numOccupiedState;

density = zeros(ntotFine, 1);
psiFine = CoarseToFine(psi, HamKS.fft); 

% FIXME factor to be simplified
occrate = reshape(occrate, 1, []);
fac = numSpin * occrate;  % fac should be row vector
density = density + sum( abs(psiFine).^2 .* fac, 2 );
       
totalCharge = sum(density);
% scale the density
density = density .* (numSpin * numOcc * ntotFine) ./ (vol * totalCharge);

% Double check (can be neglected)
totalCharge = sum( density .* vol ./ ntotFine );


end