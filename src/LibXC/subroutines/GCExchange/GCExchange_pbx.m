function [epsx, v1x, v2x] = GCExchange_pbx(rho, grho2)
% GCEXCHANGE_PBX subroutine to calculate PBE type gradient correction 
%    for exchange part of exchange correlation.
%
%    See also xcRef.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.

c1 = 0.75/pi;
c2 = 3.093667726280136;  % c2 = (3*pi^2)^(1/3)
c5 = 8/3;

k = 0.804;
% mu = 0.21951; 
mu = 0.219514972764517;  % mu = beta * (pi^2 / 3)

agrho = sqrt(grho2);
kf = c2 * rho.^(1/3);
dsg = 0.5 ./ kf;
s1 = agrho .* dsg ./ rho;

% Energy
f2 = 1 + s1.*s1*mu/k;
fx = k - k./f2;

exunif = -c1 * kf;
epsx = exunif .* fx;

% Potential
dxunif = exunif / 3;

dfx = mu * s1 ./ f2.^2;
v1x = epsx + dxunif.*fx - exunif.*dfx*c5.*s1;
v2x = exunif .* dfx .* dsg ./ agrho;

end