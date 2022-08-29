function rhoLGL = InterpPeriodicUniformFineToLGL(PeriodicUniformFineToLGLMat, numUniformGridFine, numLGLGrid, rhoUniform)
% INTERPPERIODICUNIFORMFINETOLGL transform data over uniform fine grid to
%    LGL grid.
%
%    rhoLGL = InterpPeriodicUniformFineToLGL(PeriodicUniformFineToLGLMat, 
%    numUniformGridFine, numLGLGrid, rhoUniform) transforms the data over
%    over uniform fine grid rhoUniform to data over LGL grid rhoLGL, by 
%    using transform matrix PeriodicUniformFineToLGLMat, number of uniform 
%    fine grid in 3-dim numUniformGridFine(vector), number of LGL grid in 
%    3-dim numLGLGrid(vector).
%
%    See also SCFDG/InterpPeriodicUniformFineToLGL, SCFDG/Setup.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% x-direction
vecUniformX = reshape(rhoUniform, ...
        numUniformGridFine(1), numUniformGridFine(2)*numUniformGridFine(3));
vecLGLX = PeriodicUniformFineToLGLMat{1} * vecUniformX;
temp1 = reshape(vecLGLX, ...
        numLGLGrid(1), numUniformGridFine(2), numUniformGridFine(3));

% y-direction
temp2 = permute(temp1, [2, 1, 3]);
vecUniformY = reshape(temp2, ...
        numUniformGridFine(2), numLGLGrid(1)*numUniformGridFine(3));
vecLGLY = PeriodicUniformFineToLGLMat{2} * vecUniformY;
temp2 = reshape(vecLGLY, ...
        numLGLGrid(2), numLGLGrid(1), numUniformGridFine(3));
temp2 = permute(temp2, [2, 1, 3]);

% z-direction
temp3 = permute(temp2, [3, 2, 1]);
vecUniformZ = reshape(temp3, ...
        numUniformGridFine(3), numLGLGrid(1)*numLGLGrid(2));
vecLGLZ = PeriodicUniformFineToLGLMat{3} * vecUniformZ;
temp3 = reshape(vecLGLZ, ...
        numLGLGrid(3), numLGLGrid(2), numLGLGrid(1));
temp3 = permute(temp3, [3, 2, 1]);

rhoLGL = reshape(temp3, [], 1);

end