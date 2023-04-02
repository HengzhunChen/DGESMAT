classdef Element
    % ELEMENT structure records basic information for chemical element in 
    % Periodic Table.
    %
    %    See also ElementTable.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.
    
    properties (SetAccess = public)
        z
        symbol
        config
        mass
    end
    methods
        function elem = Element(zz, s, c, m)
            elem.z = zz;
            elem.symbol = s;
            elem.config = c;
            elem.mass = m;
        end
    end
end