function rhoUniform = InterpLGLToUniform(LGLToUniformMatFine, numLGLGrid, numUniformGridFine, rhoLGL)
% INTERPLGLTOUNIFORM transform data over LGL grid to uniform fine grid.
%
%    rhoUniform = InterpLGLToUniform(LGLToUniformMatFine, numLGLGrid, 
%    numUniformGridFine, rhoLGL) transfroms the data over LGL grid rhoLGL 
%    to data over uniform fine grid rhoUniform, by using transform matrix
%    LGLToUniformMatFine, number of LGL grid in 3-dim numLGLGrid(vector)
%    and number of uniform fine grid in 3-dim numUniformGridFine(vector).
%
%    See also HamiltonianDG/InterpLGLToUniform, HamiltonianDG/Setup.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% x-direction
vecLGLX = reshape(rhoLGL, ...
        numLGLGrid(1), numLGLGrid(2)*numLGLGrid(3));
vecUniformX = LGLToUniformMatFine{1} * vecLGLX;
temp1 = reshape(vecUniformX, ...
        numUniformGridFine(1), numLGLGrid(2), numLGLGrid(3));

% y-direction
temp2 = permute(temp1, [2, 1, 3]);
vecLGLY = reshape(temp2, ...
        numLGLGrid(2), numUniformGridFine(1)*numLGLGrid(3));
vecUniformY = LGLToUniformMatFine{2} * vecLGLY;
temp2 = reshape(vecUniformY, ...
        numUniformGridFine(2), numUniformGridFine(1), numLGLGrid(3));
temp2 = permute(temp2, [2, 1, 3]);

% z-direction
temp3 = permute(temp2, [3, 2, 1]);
vecLGLZ = reshape(temp3, ...
        numLGLGrid(3), numUniformGridFine(1)*numUniformGridFine(2));
vecUniformZ = LGLToUniformMatFine{3} * vecLGLZ;
temp3 = reshape(vecUniformZ, ...
        numUniformGridFine(3), numUniformGridFine(2), numUniformGridFine(1));
temp3 = permute(temp3, [3, 2, 1]);

rhoUniform = reshape(temp3, [], 1);


end