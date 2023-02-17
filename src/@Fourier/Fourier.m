classdef Fourier
    % FOURIER class for fft    
    %    fft = Fourier() returns a Fourier class object.
    %
    %    Initialize() and InitializeFine() are used to initialize the 
    %    fourier structure over different grids of a domain.
    %
    %    See also Domain.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = public)
        domain
        isInitialized
        isInitializedFine

        % bool variables indicate for direction of fft
        forward
        
        % Laplacian operator related
        gkk
        ik
        TeterPrecond
        gkkFine
        ikFine
        TeterPrecondFine
        % Teter Preconditioner used for convergence rate of eigensolvers
        % (Teter, Payne, and Allan 1989)
                        
        % index array for mapping a coarse grid to a fine grid
        idxFineGrid
    end
    
    methods
        function fft = Fourier(varargin)
            switch (nargin)
                case 0            
                    fft.isInitialized = false;
                    fft.isInitializedFine = false;

                    fft.forward = 1;

                    fft.domain = Domain();
                case 1
                    domain = varargin{1};
                    fft = Initialize(fft, domain);
                    fft = InitializeFine(fft, domain);
                    
                    fft.forward = 1;
                otherwise
                    error('Wrong number of arguments');
            end
        end
    end
    
end