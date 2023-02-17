function HamDG = CalculateDGMatrix(HamDG)
% HAMILTONIANDG/CALCULATEDGMATRIX calculate the DG Hamiltonian matrix.
%
%    HamDG = CalculateDGMatrix(HamDG) calculates
%    (1) HamDG.elemBasisIdx, a cell whose size is numElemTotal, containing 
%        index of each basis function in the total basis functions array.
%    (2) HamDG.elemBasisInvIdx, a cell in the length of the total basis 
%        functions array, containing the element index each basis function 
%        belonging to.
%    (3) HamDG.sizeHMat, number of rows/columns of DG Hamiltonian matrix.
%    (4) HamDG.vnlCoef, HamDG.vnlDrvCoef, i.e., coefficients for the 
%        nonlocal pseudopotential. Save the coef and its derivatives in 
%        vnlCoef and vnlDrvCoef cell structures over all elements, for the 
%        later use of force computation and a posteriori error estimation.
%    (5) HamDG.HMat, DG Hamiltonian matrix, a cell with size
%        numElemTotal-by-numElemTotal, only those non-zero blocks will be
%        defined and used. 
%
%     See also HamiltonianDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


DIM = dimDef();
numElem = HamDG.numElem;
numElemTotal = prod(numElem);
numAtom = length(HamDG.atomList);

% here use the LGL grid
numLGLGrid = HamDG.grid.numLGLGridElem;

% Jump of the value of the basis, and average of the derivative of the
% basis function, each of 6 types describing the different faces along the
% X/Y/Z directions. L/R: left/right.
basisJumpXL = cell(numElemTotal, 1);
basisJumpXR = cell(numElemTotal, 1);
basisJumpYL = cell(numElemTotal, 1);
basisJumpYR = cell(numElemTotal, 1);
basisJumpZL = cell(numElemTotal, 1);
basisJumpZR = cell(numElemTotal, 1);

DbasisAverageXL = cell(numElemTotal, 1);
DbasisAverageXR = cell(numElemTotal, 1);
DbasisAverageYL = cell(numElemTotal, 1);
DbasisAverageYR = cell(numElemTotal, 1);
DbasisAverageZL = cell(numElemTotal, 1);
DbasisAverageZR = cell(numElemTotal, 1);


% The derivative of basisLGL along x,y,z directions
DbasisDimX = cell(numElemTotal, 1);
DbasisDimY = cell(numElemTotal, 1);
DbasisDimZ = cell(numElemTotal, 1);

DMat = HamDG.grid.DMat;

% Integration weights
LGLWeight2D = HamDG.grid.LGLWeight2D;
LGLWeight3D = HamDG.grid.LGLWeight3D;
LGLWeight2Dx = reshape(LGLWeight2D{1}, [], 1);
LGLWeight2Dy = reshape(LGLWeight2D{2}, [], 1);
LGLWeight2Dz = reshape(LGLWeight2D{3}, [], 1);
LGLWeight3D = reshape(LGLWeight3D, [], 1);

% pseudopotential
HamDG.vnlCoef = cell(numElemTotal, 1);
for d = 1 : DIM
    HamDG.vnlDrvCoef{d} = cell(numElemTotal, 1);
end

pseudoListElem = HamDG.pseudoListElem;
vnlCoef = cell(numElemTotal, 1);
vnlDrvCoefX = cell(numElemTotal, 1);
vnlDrvCoefY = cell(numElemTotal, 1);
vnlDrvCoefZ = cell(numElemTotal, 1); 


% Clear the DG matrix
HMat = cell(numElemTotal, numElemTotal);


% *********************************************************************
% Initial setup
% *********************************************************************

basisLGL = HamDG.basisLGL;

% compute the indices for all basis functions
elemBasisIdx = cell(numElemTotal, 1);
elemBasisInvIdx = cell(numElemTotal, 1);
count = 0;

for elemIdx = 1 : numElemTotal
    [~, n] = size(basisLGL{elemIdx});
    idxStart = count + 1;
    idxEnd = count + n;
    elemBasisIdx{elemIdx} = idxStart : idxEnd;
    for m = idxStart : idxEnd
        elemBasisInvIdx{m} = elemIdx;
    end
    count = count + n;
end

HamDG.elemBasisIdx = elemBasisIdx;
HamDG.elemBasisInvIdx = elemBasisInvIdx;
HamDG.sizeHMat = count;


% **********************************************************************
% Element-wise preparation
% **********************************************************************

% NOTE: parfor is used here
parfor elemIdx = 1 : numElemTotal

% *********************************************************************
% Compute the local derivatives
% *********************************************************************

%
% compute derivatives on each local element
%
    basis = basisLGL{elemIdx};
    [m, n] = size(basis);
    numBasis = n;
    DbasisX = zeros(m, n);
    DbasisY = zeros(m, n);
    DbasisZ = zeros(m, n);
   
    for g = 1 : numBasis
        DbasisX(:, g) = DiffPsi(DMat, numLGLGrid, basis(:, g), 1);
        DbasisY(:, g) = DiffPsi(DMat, numLGLGrid, basis(:, g), 2);
        DbasisZ(:, g) = DiffPsi(DMat, numLGLGrid, basis(:, g), 3);
    end
    
    DbasisDimX{elemIdx} = DbasisX;
    DbasisDimY{elemIdx} = DbasisY;
    DbasisDimZ{elemIdx} = DbasisZ;

%    
% compute the average of derivatives and jumps of values
%    
    % x-direction
    numGridFace = numLGLGrid(2) * numLGLGrid(3);
    valL = zeros(numGridFace, numBasis);
    valR = zeros(numGridFace, numBasis);
    drvL = zeros(numGridFace, numBasis);
    drvR = zeros(numGridFace, numBasis);
        
    % form jumps and averages from volume to face.
    for g = 1 : numBasis
        % 0.5 comes from average
        %
        % {{a}} = 1/2 (a_L + a_R)
        %
        DbasisXg = DbasisX(:, g);
        DbasisXg = reshape(DbasisXg, numLGLGrid);
        drvLg = +0.5 * DbasisXg(1, :, :);
        drvL(:, g) = reshape(drvLg, [], 1);
        drvRg = +0.5 * DbasisXg(end, :, :);
        drvR(:, g) = reshape(drvRg, [], 1);
        
        % 1.0, -1.0 comes from jump with different normal vectors
        %
        % [[a]] = -(1.0) a_L + (1.0) a_R
        %
        basisg = basis(:, g);
        basisg = reshape(basisg, numLGLGrid);
        valLg = -1.0 * basisg(1, :, :);
        valL(:, g) = reshape(valLg, [], 1);
        valRg = +1.0 * basisg(end, :, :);
        valR(:, g) = reshape(valRg, [], 1);
    end
    
    basisJumpXL{elemIdx} = valL;
    basisJumpXR{elemIdx} = valR;
    DbasisAverageXL{elemIdx} = drvL;
    DbasisAverageXR{elemIdx} = drvR;
    
    
    % y-direction
    numGridFace = numLGLGrid(1) * numLGLGrid(3);
    valL = zeros(numGridFace, numBasis);
    valR = zeros(numGridFace, numBasis);
    drvL = zeros(numGridFace, numBasis);
    drvR = zeros(numGridFace, numBasis);
        
    % form jumps and averages from volume to face.
    for g = 1 : numBasis
        % 0.5 comes from average
        %
        % {{a}} = 1/2 (a_L + a_R)
        %
        DbasisYg = DbasisY(:, g);
        DbasisYg = reshape(DbasisYg, numLGLGrid);
        drvLg = +0.5 * DbasisYg(:, 1, :);
        drvL(:, g) = reshape(drvLg, [], 1);
        drvRg = +0.5 * DbasisYg(:, end, :);
        drvR(:, g) = reshape(drvRg, [], 1);
        
        % 1.0, -1.0 comes from jump with different normal vectors
        %
        % [[a]] = -(1.0) a_L + (1.0) a_R
        %
        basisg = basis(:, g);
        basisg = reshape(basisg, numLGLGrid);
        valLg = -1.0 * basisg(:, 1, :);
        valL(:, g) = reshape(valLg, [], 1);
        valRg = +1.0 * basisg(:, end, :);
        valR(:, g) = reshape(valRg, [], 1);
    end
    
    basisJumpYL{elemIdx} = valL;
    basisJumpYR{elemIdx} = valR;
    DbasisAverageYL{elemIdx} = drvL;
    DbasisAverageYR{elemIdx} = drvR;
    
    
    % z-direction
    numGridFace = numLGLGrid(1) * numLGLGrid(2);
    valL = zeros(numGridFace, numBasis);
    valR = zeros(numGridFace, numBasis);
    drvL = zeros(numGridFace, numBasis);
    drvR = zeros(numGridFace, numBasis);
        
    % form jumps and averages from volume to face.
    for g = 1 : numBasis
        % 0.5 comes from average
        %
        % {{a}} = 1/2 (a_L + a_R)
        %
        DbasisZg = DbasisZ(:, g);
        DbasisZg = reshape(DbasisZg, numLGLGrid);
        drvLg = +0.5 * DbasisZg(:, :, 1);
        drvL(:, g) = reshape(drvLg, [], 1);
        drvRg = +0.5 * DbasisZg(:, :, end);
        drvR(:, g) = reshape(drvRg, [], 1);
        
        % 1.0, -1.0 comes from jump with different normal vectors
        %
        % [[a]] = -(1.0) a_L + (1.0) a_R
        %
        basisg = basis(:, g);
        basisg = reshape(basisg, numLGLGrid);
        valLg = -1.0 * basisg(:, :, 1);
        valL(:, g) = reshape(valLg, [], 1);
        valRg = +1.0 * basisg(:, :, end);
        valR(:, g) = reshape(valRg, [], 1);
    end
    
    basisJumpZL{elemIdx} = valL;
    basisJumpZR{elemIdx} = valR;
    DbasisAverageZL{elemIdx} = drvL;
    DbasisAverageZR{elemIdx} = drvR;


% *********************************************************************
% Nonlocal pseudopotential term, Part I
% Compute the coefficients for the nonlocal pseudopotential 
% *********************************************************************

% Compute the coefficient (i.e., the inner product of the nonlocal 
% pseduopotential and basis functions in the form of <phi|l> ) for 
% nonlocal pesudopotential projectors locally.
%
% Also get the inner product of the form <D_{x,y,z} phi | l> for
% nonlocal pseudopotential projectors locally

    coefMap     = containers.Map('KeyType', 'double', 'ValueType', 'any');
    coefDrvXMap = containers.Map('Keytype', 'double', 'ValueType', 'any');
    coefDrvYMap = containers.Map('Keytype', 'double', 'ValueType', 'any');
    coefDrvZMap = containers.Map('Keytype', 'double', 'ValueType', 'any');

    pseudoMap = pseudoListElem{elemIdx};

    % Loop over atoms, regardless of whether this atom belongs
    % to this element or not.
    atomIdxCell = keys(pseudoMap);
    atomIdxList = cell2mat(atomIdxCell);  % from cell to number
    for atomIdx = atomIdxList
        pseudo = pseudoMap(atomIdx);
        vnlList = pseudo.vnlList;
        
        if isempty(vnlList)
            coefMap(atomIdx) = [];
            coefDrvXMap(atomIdx) = [];
            coefDrvYMap(atomIdx) = [];
            coefDrvZMap(atomIdx) = [];
        else
            % using local inner product computation
            idx = vnlList.idx;
            val = vnlList.val;
            LGLval = LGLWeight3D(idx) .* val;
            coef = basis(idx, :)' * LGLval;
            coefDrvX = DbasisX(idx, :)' * LGLval;
            coefDrvY = DbasisY(idx, :)' * LGLval;
            coefDrvZ = DbasisZ(idx, :)' * LGLval;

            coefMap(atomIdx) = coef;
            coefDrvXMap(atomIdx) = coefDrvX;
            coefDrvYMap(atomIdx) = coefDrvY;
            coefDrvZMap(atomIdx) = coefDrvZ;
        end  % non-empty

    end % end of loop for atoms

    % Save coef and its derivatives in vnlCoef and vnlDrvCoef
    % structures, for the use of constructing DG matrix, force 
    % computation and a posteriori error estimation.
    vnlCoef{elemIdx} = coefMap;
    vnlDrvCoefX{elemIdx} = coefDrvXMap;
    vnlDrvCoefY{elemIdx} = coefDrvYMap;
    vnlDrvCoefZ{elemIdx} = coefDrvZMap;    
end    

HamDG.vnlCoef = vnlCoef;
HamDG.vnlDrvCoef{1} = vnlDrvCoefX;
HamDG.vnlDrvCoef{2} = vnlDrvCoefY;
HamDG.vnlDrvCoef{3} = vnlDrvCoefZ;


% *********************************************************************
% Diagonal part of the DG Hamiltonian matrix:  
% 1) Laplacian 
% 2) Local potential
% 3) Intra-element part of boundary terms
% *********************************************************************

penaltyAlpha = HamDG.penaltyAlpha;
vtotLGL = HamDG.vtotLGL;
HMatDiag = cell(numElemTotal, 1);

% NOTE: parfor is used here
parfor elemIdx = 1 : numElemTotal
    basis = basisLGL{elemIdx};    
    
    % 
    % Laplacian part
    %
    DbasisX = DbasisDimX{elemIdx};
    DbasisY = DbasisDimY{elemIdx};
    DbasisZ = DbasisDimZ{elemIdx};
    
    % Dphi_i * w * Dphi_j
    localMat = 0.5 * (...
        DbasisX' * (LGLWeight3D .* DbasisX) + ...
        DbasisY' * (LGLWeight3D .* DbasisY) + ...
        DbasisZ' * (LGLWeight3D .* DbasisZ) );     

                
    %
    % Local potential part
    %
    vtot = vtotLGL{elemIdx};
    localMat = localMat + basis' * ((LGLWeight3D .* vtot) .* basis);


    %
    % x-direction: intra-element part of the boundary term
    %            
    valL = basisJumpXL{elemIdx};
    valR = basisJumpXR{elemIdx};
    drvL = DbasisAverageXL{elemIdx};
    drvR = DbasisAverageXR{elemIdx};
    
    intByPartTerm = -0.5 * (...
        drvL' * (LGLWeight2Dx .* valL) + valL' * (LGLWeight2Dx .* drvL) + ...
        drvR' * (LGLWeight2Dx .* valR) + valR' * (LGLWeight2Dx .* drvR) );
    
    penaltyTerm = penaltyAlpha * (...
        valL' * (LGLWeight2Dx .* valL) + valR' * (LGLWeight2Dx .* valR) );
    
    localMat = localMat + intByPartTerm + penaltyTerm;

    
    %
    % y-direction: intra-element part of the boundary term
    %            
    valL = basisJumpYL{elemIdx};
    valR = basisJumpYR{elemIdx};
    drvL = DbasisAverageYL{elemIdx};
    drvR = DbasisAverageYR{elemIdx};

    intByPartTerm = -0.5 * (...
        drvL' * (LGLWeight2Dy .* valL) + valL' * (LGLWeight2Dy .* drvL) + ...
        drvR' * (LGLWeight2Dy .* valR) + valR' * (LGLWeight2Dy .* drvR) );

    penaltyTerm = penaltyAlpha * (...
        valL' * (LGLWeight2Dy .* valL) + valR' * (LGLWeight2Dy .* valR) );

    localMat = localMat + intByPartTerm + penaltyTerm;
    

    %
    % z-direction: intra-element part of the boundary term
    %            
    valL = basisJumpZL{elemIdx};
    valR = basisJumpZR{elemIdx};
    drvL = DbasisAverageZL{elemIdx};
    drvR = DbasisAverageZR{elemIdx};
    
    intByPartTerm = -0.5 * (...
        drvL' * (LGLWeight2Dz .* valL) + valL' * (LGLWeight2Dz .* drvL) + ...
        drvR' * (LGLWeight2Dz .* valR) + valR' * (LGLWeight2Dz .* drvR) );

    penaltyTerm = penaltyAlpha * (...
        valL' * (LGLWeight2Dz .* valL) + valR' * (LGLWeight2Dz .* valR) );

    localMat = localMat + intByPartTerm + penaltyTerm;


    HMatDiag{elemIdx} = localMat;
end

% add localMat to HamDG.HMat
for elemIdx = 1 : numElemTotal
    HMat{elemIdx, elemIdx} = HMatDiag{elemIdx};
end


% *********************************************************************
% Nonlocal pseudopotential term, Part II
% Update the nonlocal potential part of the matrix
% *********************************************************************

% loop over atoms
for atomIdx = 1 : numAtom
    vnlWeight = HamDG.vnlWeightMap(atomIdx);
    vnlCoef = HamDG.vnlCoef;

    % skip those atoms that have not nonlocal pseudopotential
    if isempty(vnlWeight)
        continue;
    end

    % Loop over element 1
    for elemIdx1 = 1 : numElemTotal
        coefMap1 = vnlCoef{elemIdx1};

        if isKey(coefMap1, atomIdx)
            coef1 = coefMap1(atomIdx);

            % Skip the calculation if there is no adaptive local basis 
            % function.
            if isempty(coef1)
                continue;
            end

            % Loop over element 2
            for elemIdx2 = 1 : numElemTotal
                coefMap2 = vnlCoef{elemIdx2};

                % compute the contribution to HMat{elemIdx1, elemIdx2}
                if isKey(coefMap2, atomIdx)
                    coef2 = coefMap2(atomIdx);

                    % Skip the calculation if there is no adaptive local 
                    % basis function.
                    if isempty(coef2)
                        continue;
                    end

                    % check size consistency
                    [~, numProjector1] = size(coef1);
                    [~, numProjector2] = size(coef2);
                    if numProjector1 ~= numProjector2 || ...
                       numProjector1 ~= length(vnlWeight)
                        [i1, j1, k1] = ElemIdxToKey(elemIdx1, numElem);
                        [i2, j2, k2] = ElemIdxToKey(elemIdx2, numElem);
                        msg = "Error in assembling the nonlocal pseudopotential part" + ...
                            " of the DG matrix. Atom number " + num2str(atomIdx) + ...
                            " Element 1: " + num2str([i1, j1, k1]) + ...
                            ", Element 2: " + num2str([i2, j2, k2]);
                        error(msg);
                    end

                    % outer product with the weight of the nonlocal 
                    % pseudopotential to form local matrix with
                    % size (numBasis1, numBasis2)
                    localMat = coef1 * (vnlWeight .* coef2');

                    % add to HamDG.HMat
                    if isempty(HMat{elemIdx1, elemIdx2})
                        HMat{elemIdx1, elemIdx2} = zeros(size(localMat));
                    end
                    HMat{elemIdx1, elemIdx2} = HMat{elemIdx1, elemIdx2} + localMat;

                end  % found atomIdx in element 2
            end  % end of elemIdx2
        end  % found atomIdx in element 1
    end  % end of elemIdx1
end  % end of atom list
       

% *********************************************************************
% Boundary terms, Part II
% Update the inter-element boundary part of the matrix
% *********************************************************************

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            
            %
            % x-direction
            %
            
            % elemIdxL is the previous element, elemIdxR is the current element
            if i == 1
                pre1 = numElem(1);
            else
                pre1 = i - 1;
            end
            elemIdxL = ElemKeyToIdx(pre1, j, k, numElem);
            elemIdxR = ElemKeyToIdx(i, j, k, numElem);
                        
            % NOTE that notation can be very confusing here:
            % The left element (elemIdxL) contributes to the right face (XR),
            % and the right element (elemIdxR) contributes to the left face
            % (XL).
            valL = basisJumpXR{elemIdxL};
            valR = basisJumpXL{elemIdxR};
            drvL = DbasisAverageXR{elemIdxL};
            drvR = DbasisAverageXL{elemIdxR};
            
            intByPartTerm = -0.5 * (...
                drvL' * (LGLWeight2Dx .* valR) + valL' * (LGLWeight2Dx .* drvR) );
            
            penaltyTerm = penaltyAlpha * (valL' * (LGLWeight2Dx .* valR));
            
            localMat = intByPartTerm + penaltyTerm;

            localMatTran = localMat';
            
            if isempty(HMat{elemIdxL, elemIdxR})
                HMat{elemIdxL, elemIdxR} = zeros(size(localMat));
                HMat{elemIdxR, elemIdxL} = zeros(size(localMatTran));
            end
                                        
            % add (elemIdxL, elemIdxR) to HamDG.HMat
            HMat{elemIdxL, elemIdxR} = HMat{elemIdxL, elemIdxR} + localMat;
            % add (elemIdxR, elemIdxL) to HamDG.HMat
            HMat{elemIdxR, elemIdxL} = HMat{elemIdxR, elemIdxL} + localMatTran;                 
            
            
            %
            % y-direction
            %
            
            % elemIdxL is the previous element, elemIdxR is the current element
            if j == 1
                pre2 = numElem(2);
            else
                pre2 = j - 1;
            end
            elemIdxL = ElemKeyToIdx(i, pre2, k, numElem);
            elemIdxR = ElemKeyToIdx(i, j, k, numElem);
                        
            % NOTE that notation can be very confusing here:
            % The left element (elemIdxL) contributes to the right face (YR),
            % and the right element (elemIdxR) contributes to the left face
            % (YL).
            valL = basisJumpYR{elemIdxL};
            valR = basisJumpYL{elemIdxR};
            drvL = DbasisAverageYR{elemIdxL};
            drvR = DbasisAverageYL{elemIdxR};
            
            intByPartTerm = -0.5 * (...
                drvL' * (LGLWeight2Dy .* valR) + valL' * (LGLWeight2Dy .* drvR) );

            penaltyTerm = penaltyAlpha * (valL' * (LGLWeight2Dy .* valR));

            localMat = intByPartTerm + penaltyTerm;

            localMatTran = localMat';
            
            if isempty(HMat{elemIdxL, elemIdxR})
                HMat{elemIdxL, elemIdxR} = zeros(size(localMat));
                HMat{elemIdxR, elemIdxL} = zeros(size(localMatTran));
            end
            
            % add (elemIdxL, elemIdxR) to HamDG.HMat
            HMat{elemIdxL, elemIdxR} = HMat{elemIdxL, elemIdxR} + localMat;
            % add (elemIdxR, elemIdxL) to HamDG.HMat
            HMat{elemIdxR, elemIdxL} = HMat{elemIdxR, elemIdxL} + localMatTran;                 

            
            %
            % z-direction
            %
            
            % elemIdxL is the previous element, elemIdxR is the current element
            if k == 1
                pre3 = numElem(3);
            else
                pre3 = k - 1;
            end
            elemIdxL = ElemKeyToIdx(i, j, pre3, numElem);
            elemIdxR = ElemKeyToIdx(i, j, k, numElem);
                        
            % NOTE that notation can be very confusing here:
            % The left element (elemIdxL) contributes to the right face (ZR),
            % and the right element (elemIdxR) contributes to the left face
            % (ZL).
            valL = basisJumpZR{elemIdxL};
            valR = basisJumpZL{elemIdxR};
            drvL = DbasisAverageZR{elemIdxL};
            drvR = DbasisAverageZL{elemIdxR};
            
            intByPartTerm = -0.5 * (...
                drvL' * (LGLWeight2Dz .* valR) + valL' * (LGLWeight2Dz .* drvR) );

            penaltyTerm = penaltyAlpha * (valL' * (LGLWeight2Dz .* valR));

            localMat = intByPartTerm + penaltyTerm;

            localMatTran = localMat';
            
            if isempty(HMat{elemIdxL, elemIdxR})
                HMat{elemIdxL, elemIdxR} = zeros(size(localMat));
                HMat{elemIdxR, elemIdxL} = zeros(size(localMatTran));
            end
            
            % add (elemIdxL, elemIdxR) to HamDG.HMat
            HMat{elemIdxL, elemIdxR} = HMat{elemIdxL, elemIdxR} + localMat;
            % add (elemIdxR, elemIdxL) to HamDG.HMat
            HMat{elemIdxR, elemIdxL} = HMat{elemIdxR, elemIdxL} + localMatTran;               
            
        end
    end
end

HamDG.HMat = HMat;


end