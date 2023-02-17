function valueNew = InterpByTransferMat(transferMat, numGridOld, numGridNew, valueOld)
% INTERPBYTRANSFERMAT performs the 3-dim interpolation from old grid to 
%    new grid by corresponding transform matrix.
%
%    valueNew = InterpByTransferMat(transferMat, numGridOld, numGridNew, valueOld)
%    computes the new value over the new grid by applying the transform 
%    matrix to value over old grid. Note that valueNew and valueOld are 
%    tensor in three dimension and numGridOld, numGridNew are both vector 
%    with length three, transferMat is a cell in length three with
%    interpolation matrix in three dimensions.
%
%    See also GenerateFourierInterpMat.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

% x-direction
vecOldX = reshape(valueOld, numGridOld(1), numGridOld(2)*numGridOld(3));
vecNewX = transferMat{1} * vecOldX;
temp1 = reshape(vecNewX, numGridNew(1), numGridOld(2), numGridOld(3));

% y-direction
temp2 = permute(temp1, [2, 1, 3]);
vecOldY = reshape(temp2, numGridOld(2), numGridNew(1)*numGridOld(3));
vecNewY = transferMat{2} * vecOldY;
temp2 = reshape(vecNewY, numGridNew(2), numGridNew(1), numGridOld(3));
temp2 = permute(temp2, [2, 1, 3]);

% z-direction
temp3 = permute(temp2, [3, 2, 1]);
vecOldZ = reshape(temp3, numGridOld(3), numGridNew(1)*numGridNew(2));
vecNewZ = transferMat{3} * vecOldZ;
temp3 = reshape(vecNewZ, numGridNew(3), numGridNew(2), numGridNew(1));
temp3 = permute(temp3, [3, 2, 1]);

valueNew = reshape(temp3, [], 1);

end