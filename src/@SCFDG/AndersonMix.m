function [vMix, dfMat, dvMat] = AndersonMix(scfDG, iter, vOld, vNew, dfMat, dvMat)
% SCFDG/ANDERSONMIX Anderson mixing to iter-th SCF iteration.
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
%    See also SCFDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


% *********************************************************************
% Initialize
% *********************************************************************

numElem = scfDG.numElem;
numElemTotal = prod(numElem);

% mixMaxDim is the number of columns in dfMat and dvMat. This is the the 
% mixing memory. Only information from the previous mixdim iterations are used.
mixMaxDim = scfDG.mix.mixMaxDim;
mixStepLength = scfDG.mix.mixStepLength;
mixType = scfDG.mix.mixType;

% number of iterations used, iter should start from 1
iterused = min(iter-1, mixMaxDim);
% the current position of dfMat, dvMat
ipos = iter - 1 - floor((iter-2) / mixMaxDim) * mixMaxDim;
% the next position of dfMat, dvMat
inext = iter - floor((iter-1) / mixMaxDim) * mixMaxDim;

vMix   = cell(numElemTotal, 1);
res    = cell(numElemTotal, 1);  % residual
vOpt   = cell(numElemTotal, 1);  % optimal input potential in Anderson mixing
resOpt = cell(numElemTotal, 1);  % optimal residual in Anderson mixing 


% *********************************************************************
% Anderson mixing
% *********************************************************************

for elemIdx = 1 : numElemTotal
    res{elemIdx} = vOld{elemIdx} - vNew{elemIdx};
    vOpt{elemIdx} = vOld{elemIdx};
    resOpt{elemIdx} = res{elemIdx};
    if iter > 1
        dfMat{elemIdx}(:, ipos) = res{elemIdx} - dfMat{elemIdx}(:, ipos);
        dvMat{elemIdx}(:, ipos) = vOld{elemIdx} - dvMat{elemIdx}(:, ipos);
    end
end

% For iter == 1, Anderson mixing is the same as simple mixing.
if iter > 1
    nrow = iterused;
    
    % normal matrix FTF = F^T * F
    FTF = zeros(nrow, nrow);
    % right hand side FTv = F^T * vout
    FTv = zeros(nrow, 1);
    
    for elemIdx = 1 : numElemTotal
        df = dfMat{elemIdx};
        Res = res{elemIdx};
        
        FTv = FTv + df(:, 1 : nrow)' * Res;
        FTF = FTF + df(:, 1 : nrow)' * df(:, 1 : nrow);
    end
    
    % FIXME Magic number for pseudo-inverse
    rcond = 1e-12;
    
    % FTv = pinv(FTF) * res;
    [FTv, ~] = lsqr(FTF, FTv, rcond);
    
    % update vOpt, resOpt
    % FTv = Y^{\dagger} r as in the usual notation
    for elemIdx = 1 : numElemTotal
        vOpt{elemIdx} = vOpt{elemIdx} - dvMat{elemIdx}(:, 1: nrow) * FTv;        
        resOpt{elemIdx} = resOpt{elemIdx} - dfMat{elemIdx}(:, 1: nrow) * FTv;
    end
end

% preconditioned optimal residual in Anderson mixing
if mixType == "kerker+anderson"
    precResOpt = scfDG.KerkerPrecond(resOpt);
elseif mixType == "anderson"
    precResOpt = resOpt;
else
    error('Invalid mixing type.');
end

% update dfMat, dvMat, vMix
for elemIdx = 1 : numElemTotal
    dfMat{elemIdx}(:, inext) = res{elemIdx};
    dvMat{elemIdx}(:, inext) = vOld{elemIdx};
    
    vMix{elemIdx} = vOpt{elemIdx} - mixStepLength * precResOpt{elemIdx};
end


end