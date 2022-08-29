function F = InitializeFine(F, domain)
% FOURIER/INITIALIZEFINE Initialize the Fourier variables of fine grid 
%    over domain.
%
%    See also Fourier.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


dim = dimDef();

if F.isInitializedFine
    error('Fourier has been initialized.');
end

F.domain = domain;


% construct Fourier grid over fine grid
numGridFine = domain.numGridFine;
length = domain.length;

KGridFine = cell(dim, 1);
for d = 1 : dim
    ng = numGridFine(d);
    grid = (0:ng-1) - ng *( (0:ng-1) > ng/2 );
    KGridFine{d} = grid * 2 * pi ./ length(d);  % row vector
end
   

[I, J, K] = ndgrid(KGridFine{1}, KGridFine{2}, KGridFine{3});
kkxyz = [I(:), J(:), K(:)];
F.gkkFine = sum(kkxyz .* kkxyz, 2);

F.ikFine = cell(dim, 1);
F.ikFine{1} = 1i .* kkxyz(:, 1);
F.ikFine{2} = 1i .* kkxyz(:, 2);
F.ikFine{3} = 1i .* kkxyz(:, 3);


% Teter Preconditioner
a = F.gkkFine;
b = 27 + a .*( 18 + a.*(12 + a.*8) );
F.TeterPrecondFine = b ./ ( b + 16 .* a.^4 );


% compute the index for mapping coarse to fine grid
F.idxFineGrid = zeros(domain.NumGridTotal(), 1);

ng = domain.numGrid;
ngF = domain.numGridFine;

iF = [1: floor(ng(1)/2), floor(ngF(1)/2)+1, (floor(ng(1)/2)+2 :ng(1)) + ngF(1)-ng(1)];
jF = [1: floor(ng(2)/2), floor(ngF(2)/2)+1, (floor(ng(2)/2)+2 :ng(2)) + ngF(2)-ng(2)];
kF = [1: floor(ng(3)/2), floor(ngF(3)/2)+1, (floor(ng(3)/2)+2 :ng(3)) + ngF(3)-ng(3)];

[IF, JF, KF] = ndgrid(iF, jF, kF);
idxF = IF(:) + (JF(:)-1) .* ngF(1) + (KF(:)-1) .* ngF(1) .* ngF(2);
F.idxFineGrid = idxF(:);


% mark Fourier to be initialized over fine grid
F.isInitializedFine = true;

end