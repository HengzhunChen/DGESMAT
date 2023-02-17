function y = MPentropy(x, order)
% MPENTROPY subroutine for MP (and GB) type smearing 
%
%    See also SCFDG/CalculateHarrisEnergy.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

Avec = [  1.0 / sqrt(pi), ...
         -1.0 / (4.0 * sqrt(pi)), ...
          1.0 / (32.0 * sqrt(pi)), ...
         -1.0 / (384 * sqrt(pi)) ];

y = 0.5 * Avec(order) .* Hermite_poly(x, 2*order) .* exp(- x.^2);

end