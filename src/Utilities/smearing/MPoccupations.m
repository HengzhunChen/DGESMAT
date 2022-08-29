function y = MPoccupations(x, order)
% MPOCCUPATIONS subroutine for MP (and GB) type smearing
%
%    See also MPoccupResidual.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

Avec = [  1.0 / sqrt(pi), ...
         -1.0 / (4.0 * sqrt(pi)), ...
          1.0 / (32.0 * sqrt(pi)), ...
         -1.0 / (384 * sqrt(pi)) ];
     
y = 0.5 * (1.0 - erf(x));

for m = 1 : order
    y = y + Avec(m) * Hermite_poly(x, 2*order - 1) .* exp(- x.^2);
end

end