function HamDG = CalculateDGMatrix(HamDG)
% HAMILTONIANDG/CALCULATEDGMATRIX calculate the DG Hamiltonian matrix.
%
%    HamDG = CalculateDGMatrix(HamDG) calculates
%    (1) HamDG.elemBasisIdx, a cell whose size is numElem, containing 
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

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% FIXME: check all numBasis = 0 and see whether it should be numBasisTotal
% = 0

DIM = dimDef();
numElem = HamDG.numElem;
numAtom = length(HamDG.atomList);

% here use the LGL grid
numLGLGrid = HamDG.grid.numLGLGridElem;

% Jump of the value of the basis, and average of the derivative of the
% basis function, each of size 6 describing the different faces along the
% X/Y/Z directions. L/R: left/right.
XL = 1;
XR = 2;
YL = 3;
YR = 4;
ZL = 5;
ZR = 6;
NUM_FACE = 6;

basisJump = cell(NUM_FACE, 1);
DbasisAverage = cell(NUM_FACE, 1);

% The derivative of basisLGL along x,y,z directions
Dbasis = cell(DIM, 1);
for d = 1 : DIM
    Dbasis{d} = cell(numElem);
end

% Integration weights
LGLWeight2D = HamDG.grid.LGLWeight2D;
LGLWeight3D = HamDG.grid.LGLWeight3D;


% Clear the DG matrix
HMat = cell(prod(numElem), prod(numElem));


% *********************************************************************
% Initial setup
% *********************************************************************

% compute the indices for all basis functions
elemBasisIdx = cell(numElem);
elemBasisInvIdx = cell(prod(numElem), 1);
count = 0;

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            [~, n] = size(HamDG.basisLGL{i, j, k});
            idxStart = count + 1;
            idxEnd = count + n;
            elemBasisIdx{i, j, k} = idxStart : idxEnd;
            for m = idxStart : idxEnd
                elemBasisInvIdx{m} = [i, j, k];
            end
            count = count + n;
        end
    end
end

HamDG.elemBasisIdx = elemBasisIdx;
HamDG.elemBasisInvIdx = elemBasisInvIdx;
HamDG.sizeHMat = count;



% *********************************************************************
% Compute the local derivatives
% *********************************************************************

% compute derivatives on each local element
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            basis = HamDG.basisLGL{i, j, k};
            [m, n] = size(basis);
            numBasis = n;
            DbasisX = zeros(m, n);
            DbasisY = zeros(m, n);
            DbasisZ = zeros(m, n);
            
            for g = 1 : numBasis
                DbasisX(:, g) = HamDG.DiffPsi(numLGLGrid, basis(:, g), 1);
                DbasisY(:, g) = HamDG.DiffPsi(numLGLGrid, basis(:, g), 2);
                DbasisZ(:, g) = HamDG.DiffPsi(numLGLGrid, basis(:, g), 3);
            end
            
            Dbasis{1}{i, j, k} = DbasisX;
            Dbasis{2}{i, j, k} = DbasisY;
            Dbasis{3}{i, j, k} = DbasisZ;
        end
    end
end
    
% compute the avergae of derivatives and jumps of values
for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            basis = HamDG.basisLGL{i, j, k};
            [~, numBasis] = size(basis);
            
            % x-direction
            numGridFace = numLGLGrid(2) * numLGLGrid(3);
            valL = zeros(numGridFace, numBasis);
            valR = zeros(numGridFace, numBasis);
            drvL = zeros(numGridFace, numBasis);
            drvR = zeros(numGridFace, numBasis);
            
            DbasisX = Dbasis{1}{i, j, k};
            
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
            
            basisJump{XL}{i, j, k} = valL;
            basisJump{XR}{i, j, k} = valR;
            DbasisAverage{XL}{i, j, k} = drvL;
            DbasisAverage{XR}{i, j, k} = drvR;
            
            
            % y-direction
            numGridFace = numLGLGrid(1) * numLGLGrid(3);
            valL = zeros(numGridFace, numBasis);
            valR = zeros(numGridFace, numBasis);
            drvL = zeros(numGridFace, numBasis);
            drvR = zeros(numGridFace, numBasis);
            
            DbasisY = Dbasis{2}{i, j, k};
            
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
            
            basisJump{YL}{i, j, k} = valL;
            basisJump{YR}{i, j, k} = valR;
            DbasisAverage{YL}{i, j, k} = drvL;
            DbasisAverage{YR}{i, j, k} = drvR;
            
            
            % z-direction
            numGridFace = numLGLGrid(1) * numLGLGrid(2);
            valL = zeros(numGridFace, numBasis);
            valR = zeros(numGridFace, numBasis);
            drvL = zeros(numGridFace, numBasis);
            drvR = zeros(numGridFace, numBasis);
            
            DbasisZ = Dbasis{3}{i, j, k};
            
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
            
            basisJump{ZL}{i, j, k} = valL;
            basisJump{ZR}{i, j, k} = valR;
            DbasisAverage{ZL}{i, j, k} = drvL;
            DbasisAverage{ZR}{i, j, k} = drvR;
            
        end
    end
end


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

HamDG.vnlCoef = cell(numElem);
for d = 1 : DIM
    HamDG.vnlDrvCoef{d} = cell(numElem);
end

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            coefMap     = containers.Map('KeyType', 'double', 'ValueType', 'any');
            coefDrvXMap = containers.Map('Keytype', 'double', 'ValueType', 'any');
            coefDrvYMap = containers.Map('Keytype', 'double', 'ValueType', 'any');
            coefDrvZMap = containers.Map('Keytype', 'double', 'ValueType', 'any');

            pseudoMap = HamDG.pseudoListElem{i, j, k};
            basis = HamDG.basisLGL{i, j, k};
            DbasisX = Dbasis{1}{i, j, k};
            DbasisY = Dbasis{2}{i, j, k};
            DbasisZ = Dbasis{3}{i, j, k};

            [~, numBasis] = size(basis);

            % Loop over atoms, regardless of whether this atom belongs
            % to this element or not.
            for atomIdxcell = keys(pseudoMap)
                atomIdx = atomIdxcell{1};  % from cell to number
                pseudo = pseudoMap(atomIdx);
                vnlList = pseudo.vnlList;
                
                if isempty(vnlList(1).idx)
                    coefMap(atomIdx) = [];
                    coefDrvXMap(atomIdx) = [];
                    coefDrvYMap(atomIdx) = [];
                    coefDrvZMap(atomIdx) = [];
                else
                    coef = zeros(numBasis, length(vnlList));
                    coefDrvX = zeros(numBasis, length(vnlList));
                    coefDrvY = zeros(numBasis, length(vnlList));
                    coefDrvZ = zeros(numBasis, length(vnlList));

                    % Loop over projectors
                    % implementation using local inner product computation

                    for g = 1 : length(vnlList)
                        vnl = vnlList(g);
                        idx = vnl.idx;
                        val = vnl.val(:, 1);

                        weight = reshape(LGLWeight3D, [], 1);
                        % Loop over basis function
                        for b = 1 : numBasis
                            coef(b, g) = sum(weight(idx) .* basis(idx, b) .* val);
                            coefDrvX(b, g) = sum(weight(idx) .* DbasisX(idx, b) .* val);
                            coefDrvY(b, g) = sum(weight(idx) .* DbasisY(idx, b) .* val);
                            coefDrvZ(b, g) = sum(weight(idx) .* DbasisZ(idx, b) .* val);
                        end    
                    end % end for (g)

                    coefMap(atomIdx) = coef;
                    coefDrvXMap(atomIdx) = coefDrvX;
                    coefDrvYMap(atomIdx) = coefDrvY;
                    coefDrvZMap(atomIdx) = coefDrvZ;
                end  % non-empty

            end % end of loop for atoms

            % Save coef and its derivatives in vnlCoef and vnlDrvCoef
            % structures, for the use of constructing DG matrix, force 
            % computation and a posteriori error estimation.
            HamDG.vnlCoef{i, j, k} = coefMap;
            HamDG.vnlDrvCoef{1}{i, j, k} = coefDrvXMap;
            HamDG.vnlDrvCoef{2}{i, j, k} = coefDrvYMap;
            HamDG.vnlDrvCoef{3}{i, j, k} = coefDrvZMap;    
        end
    end
end    


% *********************************************************************
% Diagonal part of the DG Hamiltonian matrix:  
% 1) Laplacian 
% 2) Local potential
% 3) Intra-element part of boundary terms
% *********************************************************************

for k = 1 : numElem(3)
    for j = 1 : numElem(2)
        for i = 1 : numElem(1)
            basis = HamDG.basisLGL{i, j, k};
            [~, numBasis] = size(basis);
            localMat = zeros(numBasis, numBasis);
            
            % In all matrix assembly process, only compute the lower
            % triangular matrix use symmetry later.
            
            % 
            % Laplacian part
            %
            DbasisX = Dbasis{1}{i, j, k};
            DbasisY = Dbasis{2}{i, j, k};
            DbasisZ = Dbasis{3}{i, j, k};
            
            % Dphi_i * w * Dphi_j
            for a = 1 : numBasis
                % for b = 1 : numBasis
                for b = 1 : a
                    localMat(a, b) = ...
                        0.5 * sum(DbasisX(:, a) .* LGLWeight3D(:) .* DbasisX(:, b)) + ...
                        0.5 * sum(DbasisY(:, a) .* LGLWeight3D(:) .* DbasisY(:, b)) + ...
                        0.5 * sum(DbasisZ(:, a) .* LGLWeight3D(:) .* DbasisZ(:, b));
                end
            end
                        
            %
            % Local potential part
            %
            vtot = HamDG.vtotLGL{i, j, k};
            for a = 1 : numBasis
                for b = 1 : a
                    localMat(a, b) = localMat(a, b) + ...
                        sum(basis(:, a) .* LGLWeight3D(:) .* vtot .* basis(:, b));
                end
            end
                        
            %
            % x-direction: intra-element part of the boundary term
            %            
            valL = basisJump{XL}{i, j, k};
            valR = basisJump{XR}{i, j, k};
            drvL = DbasisAverage{XL}{i, j, k};
            drvR = DbasisAverage{XR}{i, j, k};
            
            for a = 1 : numBasis
                for b = 1 : a
                    % integration by part term
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{1}(:) .* valL(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{1}(:) .* drvL(:, b)) ...
                        -0.5 * sum(drvR(:, a) .* LGLWeight2D{1}(:) .* valR(:, b)) ...
                        -0.5 * sum(valR(:, a) .* LGLWeight2D{1}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{1}(:) .* valL(:, b)) + ...
                        HamDG.penaltyAlpha * sum(valR(:, a) .* LGLWeight2D{1}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end

            %
            % y-direction: intra-element part of the boundary term
            %            
            valL = basisJump{YL}{i, j, k};
            valR = basisJump{YR}{i, j, k};
            drvL = DbasisAverage{YL}{i, j, k};
            drvR = DbasisAverage{YR}{i, j, k};
            
            for a = 1 : numBasis
                for b = 1 : a
                    % integration by part term
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{2}(:) .* valL(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{2}(:) .* drvL(:, b)) ...
                        -0.5 * sum(drvR(:, a) .* LGLWeight2D{2}(:) .* valR(:, b)) ...
                        -0.5 * sum(valR(:, a) .* LGLWeight2D{2}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{2}(:) .* valL(:, b)) + ...
                        HamDG.penaltyAlpha * sum(valR(:, a) .* LGLWeight2D{2}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end

            %
            % z-direction: intra-element part of the boundary term
            %            
            valL = basisJump{ZL}{i, j, k};
            valR = basisJump{ZR}{i, j, k};
            drvL = DbasisAverage{ZL}{i, j, k};
            drvR = DbasisAverage{ZR}{i, j, k};
            
            for a = 1 : numBasis
                for b = 1 : a
                    % integration by part term
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{3}(:) .* valL(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{3}(:) .* drvL(:, b)) ...
                        -0.5 * sum(drvR(:, a) .* LGLWeight2D{3}(:) .* valR(:, b)) ...
                        -0.5 * sum(valR(:, a) .* LGLWeight2D{3}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{3}(:) .* valL(:, b)) + ...
                        HamDG.penaltyAlpha * sum(valR(:, a) .* LGLWeight2D{3}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end

            %
            % Symmetrize the diagonal part of the DG matrix
            %
            diagVec = diag(localMat);
            localMat = localMat + localMat';
            localMat = localMat - diag(diagVec);
            
            % add localMat to HamDG.HMat
            key1 = i + (j-1) * numElem(1) + (k-1) * numElem(1) * numElem(2);
            key2 = key1;
            if isempty(HMat{key1, key2})
                HMat{key1, key2} = zeros(size(localMat));
            end
            HMat{key1, key2} = HMat{key1, key2} + localMat;
            
        end
    end
end


% *********************************************************************
% Nonlocal pseudopotential term, Part II
% Update the nonlocal potential part of the matrix
% *********************************************************************

% loop over atoms
for atomIdx = 1 : numAtom
    vnlWeight = HamDG.vnlWeightMap(atomIdx);

    % skip those atoms that have not nonlocal pseudopotential
    if isempty(vnlWeight)
        continue;
    end

    % Loop over element 1
    for k1 = 1 : numElem(3)
      for j1 = 1 : numElem(2)
        for i1 = 1 : numElem(1)
            coefMap1 = HamDG.vnlCoef{i1, j1, k1};

            if isKey(coefMap1, atomIdx)
                coef1 = coefMap1(atomIdx);

                % Skip the calculation if there is no adaptive
                % local basis function.
                if isempty(coef1)
                    continue;
                end

                % Loop over element 2
                for k2 = 1 : numElem(3)
                  for j2 = 1 : numElem(2)
                    for i2 = 1 : numElem(1)
                        coefMap2 = HamDG.vnlCoef{i2, j2, k2};

                        % compute the contribution to
                        % HamDG.HMat{key1, key2}

                        if isKey(coefMap2, atomIdx)
                            coef2 = coefMap2(atomIdx);

                            % Skip the calculation if there is 
                            % no adaptive local basis function.
                            if isempty(coef2)
                                continue;
                            end

                            [~, numProjector1] = size(coef1);
                            [~, numProjector2] = size(coef2);

                            % check size consistency
                            if numProjector1 ~= numProjector2 || ...
                                numProjector1 ~= length(vnlWeight)
                                msg = "Error in assembling the nonlocal pseudopotential part of the DG matrix. " + ...
                                    "Atom number " + num2str(atomIdx) + ...
                                    " Element 1: " + num2str([i1, j1, k1]) + ", Element 2: " + num2str([i2, j2, k2]);
                                error(msg);
                            end

                            % outer product with the weight of
                            % the nonlocal pseudopotential to
                            % form local matrix
                            %
                            localMat = coef1 * diag(vnlWeight) * coef2';
                            % size (numBasis1, numBasis2)

                            % add to HamDG.HMat
                            key1 = i1 + (j1 - 1) * numElem(1) + (k1 - 1) * numElem(1) * numElem(2);
                            key2 = i2 + (j2 - 1) * numElem(1) + (k2 - 1) * numElem(1) * numElem(2);

                            if isempty(HMat{key1, key2})
                                HMat{key1, key2} = zeros(size(localMat));
                            end

                            HMat{key1, key2} = HMat{key1, key2} + localMat;

                        end  % found atomIdx in element 2
                    end  % end of i2
                  end  % end of j2
                end  % end of k2
            end  % found atomIdx in element 1
        end  % end of i1
      end  % end of j1
    end  % end of k1
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
            
            % keyL is the previous element, keyR is the current element
            if i == 1
                pre1 = numElem(1);
            else
                pre1 = i - 1;
            end
                        
            % NOTE that notation can be very confusing here:
            % The left element (keyL) contributes to the right face (XR),
            % and the right element (keyR) contributes to the left face
            % (XL).
            valL = basisJump{XR}{pre1, j, k};
            valR = basisJump{XL}{i, j, k};
            drvL = DbasisAverage{XR}{pre1, j, k};
            drvR = DbasisAverage{XL}{i, j, k};
            
            [~, numBasisL] = size(valL);
            [~, numBasisR] = size(valR);
            
            localMat = zeros(numBasisL, numBasisR);
            
            % inter-element part of the boudary term
            for a = 1 : numBasisL
                for b = 1 : numBasisR
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{1}(:) .* valR(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{1}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{1}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end
            
            keyL = pre1 + (j-1) * numElem(1) + (k-1) * numElem(1) * numElem(2);
            keyR = i + (j-1) * numElem(1) + (k-1) * numElem(1) * numElem(2);
            
            localMatTran = localMat';
            
            if isempty(HMat{keyL, keyR})
                HMat{keyL, keyR} = zeros(size(localMat));
                HMat{keyR, keyL} = zeros(size(localMatTran));
            end
                                        
            % add (keyL, keyR) to HamDG.HMat
            HMat{keyL, keyR} = HMat{keyL, keyR} + localMat;
            % add (keyR, keyL) to HamDG.HMat
            HMat{keyR, keyL} = HMat{keyR, keyL} + localMatTran;                 
            
            
            %
            % y-direction
            %
            
            % keyL is the previous element, keyR is the current element
            if j == 1
                pre2 = numElem(2);
            else
                pre2 = j - 1;
            end
                        
            % NOTE that notation can be very confusing here:
            % The left element (keyL) contributes to the right face (YR),
            % and the right element (keyR) contributes to the left face
            % (YL).
            valL = basisJump{YR}{i, pre2, k};
            valR = basisJump{YL}{i, j, k};
            drvL = DbasisAverage{YR}{i, pre2, k};
            drvR = DbasisAverage{YL}{i, j, k};
            
            [~, numBasisL] = size(valL);
            [~, numBasisR] = size(valR);
            
            localMat = zeros(numBasisL, numBasisR);
            
            % inter-element part of the boudary term
            for a = 1 : numBasisL
                for b = 1 : numBasisR
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{2}(:) .* valR(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{2}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{2}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end
            
            keyL = i + (pre2 - 1) * numElem(1) + (k - 1) * numElem(1) * numElem(2);
            keyR = i + (j - 1) * numElem(1) + (k - 1) * numElem(1) * numElem(2);
            
            localMatTran = localMat';
            
            if isempty(HMat{keyL, keyR})
                HMat{keyL, keyR} = zeros(size(localMat));
                HMat{keyR, keyL} = zeros(size(localMatTran));
            end
            
            % add (keyL, keyR) to HamDG.HMat
            HMat{keyL, keyR} = HMat{keyL, keyR} + localMat;
            % add (keyR, keyL) to HamDG.HMat
            HMat{keyR, keyL} = HMat{keyR, keyL} + localMatTran;                 

            
            %
            % z-direction
            %
            
            % keyL is the previous element, keyR is the current element
            if k == 1
                pre3 = numElem(3);
            else
                pre3 = k - 1;
            end
                        
            % NOTE that notation can be very confusing here:
            % The left element (keyL) contributes to the right face (ZR),
            % and the right element (keyR) contributes to the left face
            % (ZL).
            valL = basisJump{ZR}{i, j, pre3};
            valR = basisJump{ZL}{i, j, k};
            drvL = DbasisAverage{ZR}{i, j, pre3};
            drvR = DbasisAverage{ZL}{i, j, k};
            
            [~, numBasisL] = size(valL);
            [~, numBasisR] = size(valR);
            
            localMat = zeros(numBasisL, numBasisR);
            
            % inter-element part of the boudary term
            for a = 1 : numBasisL
                for b = 1 : numBasisR
                    intByPartTerm = ...
                        -0.5 * sum(drvL(:, a) .* LGLWeight2D{3}(:) .* valR(:, b)) ...
                        -0.5 * sum(valL(:, a) .* LGLWeight2D{3}(:) .* drvR(:, b));
                    
                    penaltyTerm = ...
                        HamDG.penaltyAlpha * sum(valL(:, a) .* LGLWeight2D{3}(:) .* valR(:, b));
                    
                    localMat(a, b) = localMat(a, b) + intByPartTerm + penaltyTerm;
                end
            end
            
            keyL = i + (j - 1) * numElem(1) + (pre3 - 1) * numElem(1) * numElem(2);
            keyR = i + (j - 1) * numElem(1) + (k - 1) * numElem(1) * numElem(2);
            
            localMatTran = localMat';
            
            if isempty(HMat{keyL, keyR})
                HMat{keyL, keyR} = zeros(size(localMat));
                HMat{keyR, keyL} = zeros(size(localMatTran));
            end
            
            % add (keyL, keyR) to HamDG.HMat
            HMat{keyL, keyR} = HMat{keyL, keyR} + localMat;
            % add (keyR, keyL) to HamDG.HMat
            HMat{keyR, keyL} = HMat{keyR, keyL} + localMatTran;               
            
        end
    end
end

HamDG.HMat = HMat;


end