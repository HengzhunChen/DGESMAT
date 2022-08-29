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

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% *********************************************************************
% Initialize
% *********************************************************************

numElem = scfDG.numElem;

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

vMix   = cell(numElem);
res    = cell(numElem);  % residual
vOpt   = cell(numElem);  % optimal input potential in Anderson mixing
resOpt = cell(numElem);  % optimal  residual in Anderson mixing 


% *********************************************************************
% Anderson mixing
% *********************************************************************

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            res{i, j, k} = vOld{i, j, k} - vNew{i, j, k};
            
            vOpt{i, j, k} = vOld{i, j, k};
            
            resOpt{i, j, k} = res{i, j, k};
            
            if iter > 1
                dfMat{i, j, k}(:, ipos) = res{i, j, k} - dfMat{i, j, k}(:, ipos);
                dvMat{i, j, k}(:, ipos) = vOld{i, j, k} - dvMat{i, j, k}(:, ipos);
            end
        end
    end
end

% For iter == 1, Anderson mixing is the same as simple mixing.
if iter > 1
    nrow = iterused;
    
    % normal matrix FTF = F^T * F
    FTF = zeros(nrow, nrow);
    % right hand side FTv = F^T * vout
    FTv = zeros(nrow, 1);
    
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                df = dfMat{i, j, k};
                Res = res{i, j, k};
                
                FTv = FTv + df(:, 1 : nrow)' * Res;
                FTF = FTF + df(:, 1 : nrow)' * df(:, 1 : nrow);
            end
        end
    end
    
    % FIXME Magic number for pseudo-inverse
    rcond = 1e-12;
    
    % FTv = pinv(FTF) * res;
    [FTv, ~] = lsqr(FTF, FTv, rcond);
    
    % update vOpt, resOpt
    % FTv = Y^{\dagger} r as in the usual notation
    
    for k = 1 : numElem(3)
        for j = 1 : numElem(2)
            for i = 1 : numElem(1)
                vOpt{i, j, k} = vOpt{i, j, k} - dvMat{i, j, k}(:, 1: nrow) * FTv;
                
                resOpt{i, j, k} = resOpt{i, j, k} - dfMat{i, j, k}(:, 1: nrow) * FTv;
            end
        end
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
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            dfMat{i, j, k}(:, inext) = res{i, j, k};
            dvMat{i, j, k}(:, inext) = vOld{i, j, k};
            
            vMix{i, j, k} = vOpt{i, j, k} - mixStepLength * precResOpt{i, j, k};
        end
    end
end


end