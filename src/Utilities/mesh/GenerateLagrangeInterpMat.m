function Mat = GenerateLagrangeInterpMat(newGrids, oldGrids, domainLength)
% GENERATELAGRANGEINTERPMAT generates the transfer matrix of Barycentric
%    Lagrange interpolation from old grid to new grid in three dimension.
%
%    Mat = GenerateLagrangeInterpMat(newGrids, oldGrids, domainLength)
%    generates the transfer matrix from old grid to new grid by Lagrange
%    interpolation. newGrids, oldGrids and Mat are all cells in length 
%    three containing information in each dimension, domainLength is also 
%    a vector in length three.
%
%    NOTE: The Lagrange polynomials involved in the transfer matrix is
%    computed using the Barycentric method. For more information see
%
%    [J.P. Berrut and L.N. Trefethen, Barycentric Lagrange Interpolation,
%    SIAM Rev. 2004]
%
%    See also GenerateFourierInterpMat, InterpByTransferMat.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


DIM = dimDef();
Mat = cell(DIM, 1);

numGridOld = zeros(DIM, 1);
numGridNew = zeros(DIM, 1);

for d = 1 : DIM
    numGridOld(d) = length(oldGrids{d});
    numGridNew(d) = length(newGrids{d});
end

EPS = 1e-13;  % small stabilization parameter

for d = 1 : DIM
    oldGrid = oldGrids{d};
    newGrid = newGrids{d};
    
    % stabilization constant factor, according to Berrut and Trefethen
    stableFac = 0.25 * domainLength(d);

    % ! Here row vector or column vector is important for computation
        
    lambda = ( oldGrid - oldGrid' ) ./ stableFac;
    lambda(1 : numGridOld(d)+1 : end) = 1;  % fix diagonal value
    lambda = prod(lambda);
    lambda = 1 ./ lambda;  % size is [1, numGridOld(d)]
    
    % denominator
    denom = repmat(lambda, numGridNew(d), 1) ./ (newGrid' - oldGrid + EPS);
    denom = sum(denom, 2);  % size is [numGridNew(d), 1]
    
    localMat = (lambda ./ denom) .* (1 ./ (newGrid' - oldGrid + EPS));
    % size is [numGridNew(d), numGridOld(d)]
    
    Mat{d} = localMat;
end


% pseudo code used for reference
% for i = 1 : numGridOld(d)
%     lambda(i) = 1.0;
%     for j = 1 : numGridOld(d)
%         if j ~= i 
%             lambda(i) = lambda(i) * (oldGrids{d}(i) - oldGrids{d}(j)) / stableFac; 
%         end
%     end
%     lambda(i) = 1.0 / lambda(i);
%     for j = 1 : numGridNew(d)
%         denom(j) = denom(j) + lambda(i) / ( newGrids{d}(j) - oldGrids{d}(i) + EPS );
%     end
% end
% 
% for i = 1 : numGridOld(d)
%     for j = 1 : numGridNew(d)
%         localMat(j, i) = (lambda(i) / ( newGrids{d}(j) - oldGrids{d}(i) + EPS )) / denom(j); 
%     end
% end


end