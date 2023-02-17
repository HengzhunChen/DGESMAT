function FT = ctranspose(F)
% FOURIER/CTRANSPOSE Ctranspse function for Fourier class
%    Use to change the direction of the unitary Fourier transform, forward 
%    means from real space to reciprocal space.
%
%    Usage: FT = F';
%
%    See also Fourier.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

if nargin == 1
    if F.forward
        F.forward = 0;
    elseif ~F.forward
        F.forward = 1;
    else
        error('Not sure which direction the original transform goes');
    end
    FT = F;    
else
    error('Fourier: Invalid syntax')
end

end