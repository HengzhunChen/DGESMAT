classdef PTEntry
    % PTENTRY structure for each entry of PeriodTable.
    %
    %    See also PeriodTable, ReadUPF.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = public)
        params     % basic parameters of the chemical element
        samples    % radial grid, vLocal and nonlocal projectors and their derivatives
        cutoffs    % cutoff value for different samples
        
        nlweights  % weight of the nonlocal pseudopotential, i.e., eigenvalues of DIJ matrix
        nltypes    % type of nonlocal projectors
    end

    methods
        function PTentry = PTEntry()
            PTentry.params = struct(...
                'Znuc', [], ...
                'Mass', [], ...
                'Zval', [], ...  % number of valence charge
                'RGaussian', [] ...  % effective radius of the Gaussain charge
                );
            PTentry.samples = struct(...
                'rad', [], ...
                'vLocalSR', [], ...
                'drv_vLocalSR', [], ...
                'rhoAtom', [], ...
                'drv_rhoAtom', [], ...
                'nonlocal', [] ...
                );
            PTentry.cutoffs = struct(...
                'rad', [], ...
                'vLocalSR', [], ...  % cutoff for short range part of vLocal
                'drv_vLocalSR', [], ...
                'rhoAtom', [], ...
                'drv_rhoAtom', [], ...
                'nonlocal', [] ...
                );
        end
    end 
         
end