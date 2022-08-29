function outputOcc = PopulateMPoccup(Tsigma, MPsmearingOrder, eigvals, fermiMu)
% POPULATEMPOCCUP subroutine for MP (and GB) type smearing 
%
%    outputOcc = PopulateMPoccup(Tsigma, MPsmearingOrder, eigvals, fermiMu) 
%    fills up the the outputOcc occupations using the input eigvals, 
%    according to the Methfessel-Paxton recipe. 
%
%    See also SCFDG/CalculateOccupationRate.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

x = (eigvals - fermiMu) ./ Tsigma;

t = MPoccupations(x, MPsmearingOrder);

idx1 = t < 0.0;
t(idx1) = 0.0;

idx2 = t > 1.0;
t(idx2) = 1.0;

outputOcc = t;

end