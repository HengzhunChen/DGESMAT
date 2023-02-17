function [vMix, dfMat, dvMat] = AndersonMix(scf, iter, vOld, vNew, dfMat, dvMat)
% SCF/ANDERSONMIX Anderson mixing to iter-th SCF iteration.
%
%    NOTE: This function is written for potential mixing. Density mixing 
%    is also available but need some post processing, which should be 
%    considered outside here.
%
%    Anderson mixing is used for accelerating the SCF iteration, which is 
%    a special case for Broyden mixing.
%    The mixing can be expressed by
%
%    vnew = vin - beta* Jinv * r
%
%    where r = vout - vin, and Jinv is an approximation to the inverse
%    of the Jacobian.
% 
%    For Anderson mixing, Jinv is defined to be
%
%    Jinv = I + (dv - df)*pinv(df)
%
%    This expression can be derived from  min_Jinv || Jinv - I ||_F
%                                         s.t. dv =  Jinv*df  
%
%    See also SCF.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% mixMaxDim is the number of columns in dfMat and dvMat. This is the the 
% mixing memory. Only information from the previous mixdim iterations are used.
mixMaxDim = scf.mixing.mixMaxDim;
mixStepLength = scf.mixing.mixStepLength;
mixType = scf.mixing.mixType;

% number of iterations used, iter should start from 1
iterused = min(iter-1, mixMaxDim);
% the current position of dfMat, dvMat
ipos = iter - 1 - floor((iter-2) / mixMaxDim) * mixMaxDim;
% the next position of dfMat, dvMat
inext = iter - floor((iter-1) / mixMaxDim) * mixMaxDim;

res = vOld - vNew;  % residual
vOpt = vOld;  % Optimal input potential in Anderson mixing
resOpt = res;  % Optimal residual in Anderson mixing


if iter > 1
    dfMat(:, ipos) = res - dfMat(:, ipos);
    dvMat(:, ipos) = vOld -  dvMat(:, ipos);
    
    % calculate pseudo-inverse
    rcond = 1e-12;  % FIXME Magic number
    
    gammas = pinv(dfMat(:, 1:iterused), rcond) * res; 
    
    % update vOpt, resOpt
    vOpt = vOpt - dvMat(:, 1:iterused) * gammas;
    resOpt = resOpt - dfMat(:, 1:iterused) * gammas;
end

if mixType == "kerker+anderson"
    precResOpt = scf.KerkerPrecond(resOpt);
elseif mixType == "anderson"
    precResOpt = resOpt;
else
    error('Invalid mixing type.');
end

% Update dfMat, dvMat, vMix
dfMat(:, inext) = res;
dvMat(:, inext) = vOld;

vMix = vOpt - mixStepLength * precResOpt;

end