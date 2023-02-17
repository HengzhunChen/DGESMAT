function X = subsasgn(X, S, Y)
% SPINOR/SUBSASGN Operator overload. Subsasgn function for wave function
%
%    X(idx) = Y assign the sub wave function.
%
%    See also Spinor.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

switch S(1).type
    case '()'
        idx1 = S(1).subs{1};
        idx2 = S(1).subs{2};
        if isa(Y, 'Spinor')
            X.wavefun(idx1, idx2) = Y.wavefun;
        else
            X.wavefun(idx1, idx2) = Y;
        end
    otherwise
        X = builtin('subsasgn', X, S, Y);
end

end