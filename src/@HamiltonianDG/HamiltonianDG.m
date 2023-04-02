classdef HamiltonianDG
    % HAMILTONIANDG Main class of DGDFT for storing and assembling the 
    %    DG matrix.
    %
    %    H = HamiltonianDG(esdfParam) returns a HamiltonianDG object which 
    %    is initialized by the ESDFInputParam object esdfParam. 
    %
    %    See also Domain, Atom, Fourier, PeriodTable, ESDFInputParam.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.    
    
    properties (SetAccess = public)
        domain
        domainElem
        grid
        atomList
        fft
        ecutWavefun
        
        numElem
        numSpin
        numExtraState
        numOccupiedState

        XCType
        VDWType

        pseudoType
        pseudoListElem  % pseudo info of each atom over an extended element
        vnlCoef
        vnlDrvCoef
        vnlWeightMap
        
        pseudoCharge  % Gaussian compensation charge
        density
        gradDensity
        atomDensity
        densityLGL
        
        vext
        vLocalSR
        vhart
        vxc
        epsxc
        vtot
        vtotLGL
        
        EVdw
        Eself
        EIonSR
        Eext
        
        forceVdw
        forceIonSR
        forceExt
        
        basisLGL
        elemBasisIdx  % index of a basis in total basis function array
        elemBasisInvIdx  % element index a basis belonging to
        HMat  % DG Hamiltonian matrix, stores in cell
        sizeHMat
        penaltyAlpha
        
        eigVal
        occupationRate
        eigvecCoef

        userOption
        
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
                'numUniformGridElemFine', [], ...
                ...
                'LGLGridElem', [], ...
                'numLGLGridElem', [], ...
                'LGLWeight1D', [], ...
                'LGLWeight2D', [], ...
                'LGLWeight3D', [], ...
                'DMat', [], ... % differentiation matrix on the LGL grid
                ...
                'LGLToUniformMat', [], ...
                'LGLToUniformMatFine', [] ...
                );

            H.userOption = struct(...
                ...
                );
                            
            switch (nargin)
                case 1
                    esdfParam = varargin{1};
                    H.domain = Domain();
                    H.domainElem = Domain();
                    H = Setup(H, esdfParam);
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
                
        function sparseH = ElemMatToSparse(HamDG)
            % transform an element-wise matrix into a sparse matrix
            sparseH = ElemMatToSparse(...
                        HamDG.HMat, HamDG.numElem, HamDG.sizeHMat);
            sparseH = (sparseH + sparseH') / 2;
        end

    end
end