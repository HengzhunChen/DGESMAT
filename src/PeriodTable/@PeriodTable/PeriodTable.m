classdef PeriodTable
    % PERIODTABLE class for potential information, especially for pseudo
    %    potential for each type of atom being used.
    %
    %    ptable = PeriodTable(esdfParam) returns a PeriodTable object 
    %    according to info from esdfParam.
    %
    %    See also ESDFInputParam, ReadUPF, PTEntry, Atom.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.
    
    
    properties (SetAccess = private)        
        ptNLtype = struct(...
            'L0', 0, ...
            'L1', 1, ...
            'L2', 2, ...
            'L3', 3, ...
            'SpinorOrbit_L1', -1, ...
            'SpinorOrbit_L2', -2, ...
            'SpinorOrbit_L3', -3 ...
            );  % type of nonlocal projectors
            
        pteMap     % map, from atomic number to PTEntry
        splineMap  % map, from atomic number to spline coefficients for pseudopotentials
                   % struct of spline see also Setup() function of PeiodTable

        userOption = struct(...
            ...
            );
    end

    methods
        function PT = PeriodTable(esdfParam)
            PT.pteMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
            PT.splineMap = containers.Map('KeyType', 'double', 'ValueType', 'any');

            PT = Setup(PT, esdfParam);
        end
        
        function flag = IsNonlocal(PT, atomNum)
            % Whether the atom type has nonlocal pseudopotential
            if ~isempty(PT.pteMap(atomNum).samples.nonlocal)  
                flag = true;
            else
                flag = false;
            end
        end
        
        function val = RcutVLocalSR(PT, atomNum)
            % cutoff radius for VLocal short range part in the real space
            val = PT.pteMap(atomNum).cutoffs.vLocalSR;
        end
        
        function val = RcutRhoAtom(PT, atomNum)
            % cutoff radius for model atomic density in the real space. 
            % This is only used for constructing initial charge density,  
            % and does not need to be very accurate 
            val = PT.pteMap(atomNum).cutoffs.rhoAtom;
        end
        
        function val = RcutNonlocal(PT, atomNum)
            % cutoff radius for the nonlocal pseudopotential in the real 
            % space.
            % If there are multiple pseudopotentials, Rcut should be the
            % maximum radius so that ALL nonlocal pseudopotentials are
            % accurate.
            if PT.IsNonlocal(atomNum)
                val = max(PT.pteMap(atomNum).cutoffs.nonlocal);
            else
                val = 0.0;
            end
        end
        
        function val = Mass(PT, atomNum)
            val = PT.pteMap(atomNum).params.Mass;
        end
        
        function val = Zval(PT, atomNum)
            val = PT.pteMap(atomNum).params.Zval;
        end
        
        function val = RGaussian(PT, atomNum)
            val = PT.pteMap(atomNum).params.RGaussian;
        end
        
    end
end