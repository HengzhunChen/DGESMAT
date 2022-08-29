classdef Spinor
    % SPINOR class of wave function for the global domain or extended 
    %    element.
    %
    %    psi = Spinor() returns an empty Spinor class object.
    %
    %    psi = Spinor(domain, numStateList, data) returns a Spinor object 
    %    with corresponding input, data can be a number or matrix,
    %    numStateList is a list for number of states of each spin component.
    %
    %    NOTE: component means Two-component spin, spin-up and spin-down 
    %    or Four-component spin, large/small spin-up/down.
    %
    %    See also Domain.

    %  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
    %                     Fudan University
    %  This file is distributed under the terms of the MIT License.
    
    
    properties (SetAccess = private)
        domain     % mesh should be used here for general cases 
        
        numComponent
        numStateList  % list for number of states of each spin component
        
        trans  % indicator for transpose        
    end
    
    properties (SetAccess = public)
        % local data of the wavefunction
        wavefun
    end
    
    methods
        function X = Spinor(varargin)
            switch (nargin)
                case 0
                    X.domain = Domain();
                    return;
                case 3
                    domain = varargin{1};
                    numStateList = varargin{2};
                    data = varargin{3};

                    X.domain = domain;
                    X.numStateList = numStateList;
                    X.numComponent = length(numStateList);

                    X.trans = 0;

                    if length(data) > 1
                        X.wavefun = data;
                    else
                        % data is 1-by-1 double
                        X.wavefun = data * ...
                            ones(domain.NumGridTotal(), sum(numStateList));
                    end 

                otherwise
                    error('Wrong number of arguments');
            end
        end
        
        % -------------------- Inquiries ---------------------------
        
        function val = NumGridTotal(Psi)
            % NOTE: this may not equal to Psi.domain.NumGridTotal()
            [val, ~] = size(Psi.wavefun);
        end
        
        function val = NumStateTotal(Psi)
            val = sum(Psi.numStateList);
        end
        
    end
        
end