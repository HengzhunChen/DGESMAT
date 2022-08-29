function Y = horzcat(varargin)
% SPINOR/HORZCAT Operator overload. Horzcat function for wave function
%
%    Y = horzcat(X1, X2, ...) returns the horizontal cat of wave function
%    from X1, X2, ...
%
%    Y = [X1, X2, ...] returns the horizontal cat of wave functions.
%
%    NOTE: If there are multi-components in Xs to be horizontal 
%    concatenated, output Y will have only one component.
%
%    See also Spinor.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

X = varargin{1};

nXcols = cellfun(@ncols, varargin);
Ypsi = zeros(nrows(X), sum(nXcols));

idx = 0;
for j = 1 : length(varargin)
    idx = idx(end) + (1 : nXcols(j));
    Ypsi(:, idx) = varargin{j}.wavefun;
end

[~, numState] = size(Ypsi);
Y = Spinor(X.domain, numState, Ypsi);
Y.trans = X.trans;

end