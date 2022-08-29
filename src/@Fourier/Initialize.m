function F = Initialize(F, domain)
% FOURIER/INITIALIZE Initialize the Fourier variables of coarse grid over 
%   domain.
%
%    See also Fourier.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


dim = dimDef();

if F.isInitialized
    error('Fourier has been initialized.');
end

F.domain = domain;


% Fourier grid over coarse grid
numGrid = domain.numGrid;
length = domain.length;

KGrid = cell(dim, 1);
for d = 1 : dim
    ng = numGrid(d);
    grid = (0:ng-1) - ng *( (0:ng-1) > ng/2 );
    KGrid{d} = grid * 2 * pi ./ length(d);  % row vector
end
   

[I, J, K] = ndgrid(KGrid{1}, KGrid{2}, KGrid{3});
kkxyz = [I(:), J(:), K(:)];
F.gkk = sum(kkxyz .* kkxyz, 2);

F.ik = cell(dim, 1);
F.ik{1} = 1i .* kkxyz(:, 1);
F.ik{2} = 1i .* kkxyz(:, 2);
F.ik{3} = 1i .* kkxyz(:, 3);

% Teter Preconditioner
% used for convergence rate of eigensolvers
a = F.gkk;
b = 27 + a .*( 18 + a.*(12 + a.*8) );
F.TeterPrecond = b ./ ( b + 16 .* a.^4 );


% mark Fourier to be initialized over coarse grid
F.isInitialized = true;

end