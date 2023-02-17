function y = mtimes(F, x)
% FOURIER/MTIMES Perform the Fourier transform between real space and
%    reciprocal space.
%
%    Usage: y = F * x;
%    x can be numerical data or Spinor Object.
%
%    See also Fourier, Fourier/ctranpose

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


if nargin == 2
    numGrid = F.domain.numGrid;
    numGridFine = F.domain.numGridFine;
    ntot = F.domain.NumGridTotal();
    ntotFine = F.domain.NumGridTotalFine();
    vol = F.domain.Volume();

    if isa(x, 'numeric')
        [nrows, ~] = size(x);
        if nrows == ntot
            if ~F.forward
                y = ifft3Execute(x, numGrid, vol);
            elseif F.forward
                y = fft3Execute(x, numGrid, vol);
            else
                error('Fourier: something wrong with the FFT configuration');
            end
        elseif nrows == ntotFine
            if ~F.forward
                y = ifft3Execute(x, numGridFine, vol);
            elseif F.forward
                y = fft3Execute(x, numGridFine, vol);
            else
                error('Fourier: something wrong with the FFT configuration');
            end
        else
            error('the number of rows in x does not match with the Fourier object');
        end
    elseif isa(x, 'Spinor')
        psi = x;
        x = psi.wavefun;
        [nrows, ~] = size(x);
        if nrows == ntot
            if ~F.forward
                y = ifft3Execute(x, numGrid, vol);
            elseif F.forward
                y = fft3Execute(x, numGrid, vol);
            else
                error('Fourier: something wrong with the FFT configuration');
            end
        elseif nrows == ntotFine
            if ~F.forward
                y = ifft3Execute(x, numGridFine, vol);
            elseif F.forward
                y = fft3Execute(x, numGridFine, vol);
            else
                error('Fourier: something wrong with the FFT configuration');
            end
        else
            error('the number of rows in x does not match with the Fourier object');
        end
        psi.wavefun = y;
        y = psi;
    end
else
    error(' Fourier Syntax: y=F*x ');
end
        
end