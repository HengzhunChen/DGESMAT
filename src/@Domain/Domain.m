classdef Domain
    % DOMAIN class for computational domain
    %    dm = Domain() returns a Domain class object
    %
    %    NOTE:
    %    numGrid used for wavefun and numGridFine used for density
    %    starting position mainly used for subdomain, 
    %    i.e., element in DGDFT
    %

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.
    
    properties (SetAccess = public)
        length          % length in each direction
        posStart        % starting position
        numGrid         % number of coarse grids points in each direction
        numGridFine     % number of fine grids points in each direction
    end
    
    methods
        function domain = Domain()
            domain.length = [0, 0, 0];
            domain.posStart = [0, 0, 0];
            domain.numGrid = [0, 0, 0];
            domain.numGridFine = [0, 0, 0];         
        end
    end
end