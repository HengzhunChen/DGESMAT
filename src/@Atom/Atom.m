classdef Atom
    % Atom class for atom information
    %
    %    atom = Atom() returns an empty Atom class object.
    %
    %    atom = Atom(type, pos, vel, force) returns an Atom object with
    %    atomic number type, position pos, velocity vel and force force.
    %
    %    See also ESDFInputParam, ESDFReadInput, PeriodTable. 

    %  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
    %                     Fudan University
    %  This file is distributed under the terms of the MIT License.

    properties (SetAccess = public)
        type       % Atomic number
        pos        % Position
        vel        % Velocity
        force      % Force
    end
    
    methods
        function a = Atom(varargin)
            switch (nargin)
                case 0
                    return;
                case 4
                    t = varargin{1};
                    p = varargin{2};
                    v = varargin{3};
                    f = varargin{4};            
                    a.type = t;
                    a.pos = p;
                    a.vel = v;
                    a.force = f;
               otherwise
                    error('Wrong number of arguments');
            end
        end
    end
end