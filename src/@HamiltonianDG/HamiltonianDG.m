classdef HamiltonianDG
    % HAMILTONIANDG Main class of DGDFT for storing and assembling the 
    %    DG matrix.
    %
    %    H = HamiltonianDG() returns a HamiltonianDG object which was 
    %    initialized by the global variable esdfParam. 
    %
    %    See also Domain, Atom, Fourier, PeriodTable.

    %  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
    %                     Fudan University
    %  This file is distributed under the terms of the MIT License.    
    
    properties (SetAccess = public)
        domain
        domainElem
        grid
        atomList
        fft
        
        numElem
        numSpin
        numExtraState
        numOccupiedState
        
        pseudoType
        pseudoListElem  % pseudo info of each atom over an extended element
        vnlCoef
        vnlDrvCoef
        vnlWeightMap
        
        pseudoCharge
        density
        gradDensity
        atomDensity
        densityLGL
        
        vext
        vhart
        vxc
        epsxc
        vtot
        vtotLGL
        
        XCType
        
        basisUniformFine
        basisLGL
        elemBasisIdx  % index of a basis in total basis function array
        elemBasisInvIdx  % element index a basis belonging to
        HMat  % DG Hamiltonian matrix, stores in cell
        sizeHMat
        sparseHMat  % HMat stores in sparse matrix
        hasConvertSparse  % whether HMat save in sparseHMat 
        penaltyAlpha
        
        eigVal
        occupationRate
        eigvecCoef
        
    end
    
    methods
        function H = HamiltonianDG(varargin)
            H.grid = struct(...
                'uniformGrid', [], ...
                'uniformGridFine', [], ...
                ...
                'uniformGridElem', [], ...
                'uniformGridElemFine', [], ...
                'numUniformGridElem', [], ...
                'numUniforomGridElemFine', [], ...
                ...
                'numLGLGridElem', [], ...
                'LGLGridElem', [], ...
                'LGLWeight1D', [], ...
                'LGLWeight2D', [], ...
                'LGLWeight3D', [], ...
                'DMat', [], ... % differentiation matrix on the LGL grid
                ...
                'LGLToUniformMat', [], ...
                'LGLToUniformMatFine', [], ...
                'LGLToUniformGaussMatFine', [] ...
                );
                            
            switch (nargin)
                case 0
                    H.domain = Domain();
                    H.domainElem = Domain();
                    H = Setup(H);
                    return;
                otherwise
                    error('Wrong number of arguments');
            end            
        end
        
        function val = NumBasisTotal(HamDG)
            val = HamDG.sizeHMat;
        end
        function val = NumStateTotal(HamDG)
            val = HamDG.numExtraState + HamDG.numOccupiedState;
        end
        
        function rhoUniform = InterpLGLToUniform(HamDG, rhoLGL)
            % transform data over LGL grid to uniform fine grid
            rhoUniform = InterpLGLToUniform(...
                HamDG.grid.LGLToUniformMatFine, ...
                HamDG.grid.numLGLGridElem, ...
                HamDG.grid.numUniformGridElemFine, ...
                rhoLGL);
        end
        
        function sparseH = ElemMatToSparse(HamDG)
            % transform an element-wise matrix into a sparse matrix
            % note: here does not change hasConvertSparse
            sparseH = ElemMatToSparse(...
                HamDG.HMat,  ...
                HamDG.numElem, ...
                HamDG.sizeHMat);
            sparseH = (sparseH + sparseH') / 2;
        end

    end
end