classdef HamiltonianKS
    % HAMILTONIANKS class of Hamiltonian for planewave basis 
    %    diagonalization method.
    %
    %    H = HamiltonianKS() returns an empty HamiltonianKS object.
    %
    %    H = HamiltonianKS(esdfParam) returns a HamiltonianKS object
    %    according to parameters in ESDFInputParam object esdfParam.
    %
    %    H = HamiltonianKS(esdfParam, domain, atomList) returns a 
    %    HamiltonianKS object with respect to ESDFInputParam object 
    %    esdfParam, domain and atomList.
    %
    %    H = HamiltonianKS(esdfParam, domain, atomList, ptable) returns a
    %    HamiltonianKS object with respect to ESDFInputParam object 
    %    esdfParam, domain, atomList and PeriodTale object ptable.
    %
    %    H = HamiltonianKS(esdfParam, domain, atomList, ptable, fft) 
    %    returns a HamiltonianKS object with respect to ESDFInputParam 
    %    object esdfParam domain, atomList, PeriodTable ptable and 
    %    Fourier fft.
    %  
    %    See also Domain, Atom, Fourier, PeriodTable, xcRef.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = public)
        domain
        atomList
        ptable
        fft
        ecutWavefunction
        
        numSpin
        numExtraState
        numExtraElectron
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
        VDWType
        
        pseudoList  % pseudo info in each atom
        ecutFilter  % used in CheFSI

        userOption
        isDGDFT  % whether this is a subproblem of DGDFT
    end        
    
    
    methods
        function H = HamiltonianKS(varargin)
            
            H.pseudoList = struct(...
                'pseudoCharge', [], ...
                'vLocalSR', [], ...  % used if isUseVLocal
                'vnlList', [] ...
                );
            
            H.ecutFilter = struct(...
                'isApplyFilter', [], ...
                'wfnCutoff', [] ...
                );
            
            H.userOption = struct(...
                'isUseVLocal', [] ...
                );
            
            switch (nargin)
                case 0
                    H.domain = Domain();
                    return;
                case 1
                    esdfParam = varargin{1};                    
                    H.domain = esdfParam.basic.domain;
                    H.atomList = esdfParam.basic.atomList;
                    ptable = PeriodTable(esdfParam);
                    H.ptable = ptable;
                    H.fft = Fourier(esdfParam.basic.domain);
                    H = Setup(H, esdfParam);
                    return;                    
                case 3
                    esdfParam = varargin{1};
                    domain = varargin{2};
                    atomList = varargin{3};
                    
                    H.domain = domain;
                    H.atomList = atomList;                    
                    ptable = PeriodTable(esdfParam);
                    H.ptable = ptable;
                    H.fft = Fourier(domain);
                    H = Setup(H, esdfParam);
                    return;
                case 4
                    esdfParam = varargin{1};
                    domain = varargin{2};
                    atomList = varargin{3};
                    ptable = varargin{4};
                    
                    H.domain = domain;
                    H.atomList = atomList;
                    H.ptable = ptable;
                    H.fft = Fourier(domain);
                    H = Setup(H, esdfParam);
                    return;
                case 5
                    esdfParam = varargin{1};
                    domain = varargin{2};
                    atomList = varargin{3};
                    ptable = varargin{4};
                    fft = varargin{5};
                    
                    H.domain = domain;
                    H.atomList = atomList;
                    H.ptable = ptable;
                    H.fft = fft;
                    H = Setup(H, esdfParam);
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