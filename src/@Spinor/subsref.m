function Y = subsref(X, S)
% SPINOR/SUBSREF Subref function for wave function
%
%    Y = X(idx) returns the sub wave function.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

if numel(S) > 1
    Y = builtin('subsref', X, S);
    return;
end

switch S(1).type
    case '()'
        rows = S(1).subs{1};
        cols = S(1).subs{2};
        Y = X;
        Y.wavefun = X.wavefun(rows, cols);
    otherwise
        Y = builtin('subsref', X, S);
end

end