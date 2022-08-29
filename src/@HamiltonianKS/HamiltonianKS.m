classdef HamiltonianKS
    % HAMILTONIANKS class of Hamiltonian for planewave basis 
    %    diagonalization method.
    %
    %    H = HamiltonianKS() returns an empty HamiltonianKS object.
    %
    %    H = HamiltonianKS(esdfParam) returns a HamiltonianKS object
    %    according to parameters in esdfParam.
    %
    %    H = HamiltonianKS(domain, atomList) returns a HamiltonianKS 
    %    object with respect to domain and atomList.
    %
    %    H = HamiltonianKS(domain, atomList, ptable) returns a
    %    HamiltonianKS object with respect to domain, atomList and 
    %    PeriodTale ptable.
    %
    %    H = HamiltonianKS(domain, atomList, ptable, fft) returns a
    %    HamiltonianKS object with respect to domain, atomList, 
    %    PeriodTable ptable and Fourier fft.
    %  
    %    See also Domain, Atom, Fourier, PeriodTable, xcRef.

    %  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
    %                     Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = public)
        domain
        atomList
        ptable
        fft
        
        numSpin
        numExtraState
        numOccupiedState
        
        eigVal
        occupationRate
                
        vLocalSR
        vext
        vhart
        vxc
        epsxc
        vtot

        pseudoCharge
        gaussianCharge  % used if isUseVLocal
        density
        gradDensity
        atomDensity
        
        EVdw
        Eself
        EIonSR
        Eext
        
        forceVdw
        forceIonSR
        forceExt
        
        pseudoType
        XCType
        isHybrid
        isEXXActive        
        
        pseudoList  % pseudo info in each atom
        ecutFilter  % used in CheFSI
        hybrid
        exx
    end        
    
    
    methods
        function H = HamiltonianKS(varargin)
            
            H.pseudoList = struct(...
                'pseudoCharge', [], ...
                'vLocalSR', [], ...  % used if isUseVLocal
                'vnlList', [] ...
                );

            H.hybrid = struct(...
                'ScreenMu', 0.106, ...  
                'ExxFraction', 0.25, ...  
                'DFType', "", ...
                'DFKmeansWFType', "", ...
                'DFKmeansWFAlpha', [], ...
                'DFKmeansTolerance', [], ...
                'DFKmeansMaxIter', [], ...
                'DFNumMu', [], ...
                'DFNumGaussianRandom', [], ...
                'DFTolerance', [], ...
                'exxDivergenceType', [], ...
                'exxDiv', [] ...
                );
            
            H.ecutFilter = struct(...
                'isApplyFilter', [], ...
                'wfnCutoff', [] ...
                );
            
            H.exx = struct(...
                'phiEXX', [], ...
                'vexxProj', [], ...
                'exxgkk', [] ...
                );
            
            switch (nargin)
                case 0
                    H.domain = Domain();
                    return;
                case 1
                    esdfParam = varargin{1};
                    
                    H.domain = esdfParam.basic.domain;
                    H.atomList = esdfParam.basic.atomList;
                    ptable = PeriodTable();
                    H.ptable = ptable;
                    H.fft = Fourier(esdfParam.basic.domain);
                    H = Setup(H);
                    return;                    
                case 2
                    domain = varargin{1};
                    atomList = varargin{2};
                    
                    H.domain = domain;
                    H.atomList = atomList;                    
                    ptable = PeriodTable();
                    H.ptable = ptable;
                    H.fft = Fourier(domain);
                    H = Setup(H);
                    return;
                case 3
                    domain = varargin{1};
                    atomList = varargin{2};
                    ptable = varargin{3};
                    
                    H.domain = domain;
                    H.atomList = atomList;
                    H.ptable = ptable;
                    H.fft = Fourier(domain);
                    H = Setup(H);
                    return;
                case 4
                    domain = varargin{1};
                    atomList = varargin{2};
                    ptable = varargin{3};
                    fft = varargin{4};
                    
                    H.domain = domain;
                    H.atomList = atomList;
                    H.ptable = ptable;
                    H.fft = fft;
                    H = Setup(H);
                    return;
                otherwise
                    error('Wrong number of arguments');
            end
        end
        
        function val = NumStateTotal(H)
            val = H.numExtraState + H.numOccupiedState;
        end
        
    end
end