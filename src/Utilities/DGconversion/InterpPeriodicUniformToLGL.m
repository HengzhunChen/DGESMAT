function psiLGL = InterpPeriodicUniformToLGL(PeriodicUniformToLGLMat, numUniformGrid, numLGLGrid, psiUniform)
% INTERPPERIODICUNIFORMTOLGL transform data over uniform grid to LGL grid.
%
%    psiLGL = InterpPeriodicUniformToLGL(PeriodicUniformToLGLMat, 
%    numUniformGrid, numLGLGrid, psiUniform) transforms the data over
%    uniform grid psiUniform to data over LGL grid psiLGL, by using 
%    transform matrix PeriodicUniformToLGLMat, number of uniform grid in 
%    3-dim numUniformGrid(vector), number of LGL grid in 3-dim 
%    numLGLGrid(vector).
%
%    See also SCFDG/InterpPeriodicUniformToLGL, SCFDG/Setup.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% x-direction
vecUniformX = reshape(psiUniform, ...
        numUniformGrid(1), numUniformGrid(2)*numUniformGrid(3));
vecLGLX = PeriodicUniformToLGLMat{1} * vecUniformX;
temp1 = reshape(vecLGLX, ...
        numLGLGrid(1), numUniformGrid(2), numUniformGrid(3));

% y-direction
temp2 = permute(temp1, [2, 1, 3]);
vecUniformY = reshape(temp2, ...
        numUniformGrid(2), numLGLGrid(1)*numUniformGrid(3));
vecLGLY = PeriodicUniformToLGLMat{2} * vecUniformY;
temp2 = reshape(vecLGLY, ...
        numLGLGrid(2), numLGLGrid(1), numUniformGrid(3));
temp2 = permute(temp2, [2, 1, 3]);

% z-direction
temp3 = permute(temp2, [3, 2, 1]);
vecUniformZ = reshape(temp3, ...
        numUniformGrid(3), numLGLGrid(1)*numLGLGrid(2));
vecLGLZ = PeriodicUniformToLGLMat{3} * vecUniformZ;
temp3 = reshape(vecLGLZ, ...
        numLGLGrid(3), numLGLGrid(2), numLGLGrid(1));
temp3 = permute(temp3, [3, 2, 1]);

psiLGL = reshape(temp3, [], 1);
 

end