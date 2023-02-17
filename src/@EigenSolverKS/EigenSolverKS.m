classdef EigenSolverKS
    % EIGENSOLVERKS class of eigensolvers for planewave-basis Kohn Sham DFT 
    %    in the global domain or extended element.
    %
    %    eigSol = EigenSolverKS() returns an EigenSolverKS object. 
    %
    %    eigSol = EigenSolverKS(hamKS, psi) returns an EigenSolverKS object
    %    with respect to HamiltonianKS object hamKS and Spinor object psi.
    %
    %    Options of eigensolver:
    %        - LOBPCG (default)
    %        - eigs
    %        - PPCG
    %        - CheFSI
    %    
    %    See also HamiltonianKS, Spinor, Fourier, SCF, SCFDG.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.
    
    properties (SetAccess = public)
       eigVal
       resVal
       hamKS
       psi
       fft
    end
        
    methods 
        function eigSol = EigenSolverKS(varargin)
            switch (nargin)
                case 0
                    eigSol.hamKS = HamiltonianKS();
                    eigSol.fft = Fourier();
                    eigSol.psi = Spinor();
                    return;
                case 2
                    hamKS = varargin{1};
                    psi = varargin{2};
                    
                    eigSol.hamKS = hamKS;
                    eigSol.psi = psi;
                    eigSol.fft = hamKS.fft;
                    eigSol.eigVal = zeros(psi.NumStateTotal(), 1);
                    eigSol.resVal = zeros(psi.NumStateTotal(), 1);
                otherwise
                    error('Wrong number of arguments');
            end
        end
    end
       
end