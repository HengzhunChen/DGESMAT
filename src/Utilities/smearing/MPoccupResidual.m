function y = MPoccupResidual(Tsigma, MPSmearingOrder, eigvals, fermiMu, numSolve)
% MPOCCUPRESIDUAL subroutine for MP (and GB) type smearing
%
%    This computes the residual of (\sum_i f_i ) - n_e used for computing 
%    the Fermi level.
%
%    See also SCFDG/CalculateOccupationRate.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

x = (eigvals - fermiMu) ./ Tsigma;
t = MPoccupations(x, MPSmearingOrder);

idx1 = t < 0.0;
t(idx1) = 0.0;
idx2 = t > 1.0;
t(idx2) = 1.0;

y = sum(t) - numSolve;

end