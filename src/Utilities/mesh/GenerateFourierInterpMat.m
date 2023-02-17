function Mat = GenerateFourierInterpMat(newGrid, oldGrid, domainLength)
% GENERATEFOURIERINTERPMAT generates the transfer matrix of Fourier
%    interpolation from old grid to new grid in three dimension.
%
%    Mat = GenerateFourierInterpMat(newGrid, oldGrid, domainLength)
%    generates the transfer matrix from old grid to new grid by Fourier
%    interpolation. newGrid, oldGrid and Mat are all cells in length three
%    containing information in each dimension, domainLength is also a
%    vector in length three.
%
%    See also GenerateLagrangeInterpMat, InterpByTransferMat.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


DIM = dimDef();
Mat = cell(DIM, 1);

numGridOld = zeros(DIM, 1);
numGridNew = zeros(DIM, 1);

for d = 1 : DIM
    numGridOld(d) = length(oldGrid{d});
    numGridNew(d) = length(newGrid{d});
end

for d = 1 : DIM
    ng = numGridOld(d);
    grid = (0:ng-1) - ng *( (0:ng-1) > ng/2 );
    KGrid = grid * 2 * pi ./ domainLength(d);  % row vector
    
    tempMat = newGrid{d}' - oldGrid{d};
    tempTns = repmat(tempMat, [1, 1, numGridOld(d)]); 
    % size [numGridNew(d), numGridOld(d), numGridOld(d)]
    tempKGrid = repmat(KGrid, [numGridNew(d), 1, numGridOld(d)]);
    tempKGrid = permute(tempKGrid, [1, 3, 2]);  
    % size [numGridNew(d), numGridOld(d), numGridOld(d)], i-th page is all KGrid(i)
    Mat{d} = sum( cos(tempTns .* tempKGrid), 3 ) ./ numGridOld(d);  
    % size [numGridNew(d), numGridOld(d)]
end


% pseudo code used for reference
% for i = 1 : newGrid(d)
%     for j = 1 : oldGrid(d)
%         Mat(i, j) = 0.0;
%         for k = 1 : oldGrid(d)
%             Mat(i, j) = Mat(i, j) + ...
%                 cos( KGrid(k) * ( newGrid{d}(i) - oldGrid{d}(j) ) ) / numGridOld(d);
%         end
%     end
% end


end