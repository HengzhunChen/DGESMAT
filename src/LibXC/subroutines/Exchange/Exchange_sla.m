function [epsx, vx] = Exchange_sla(rs)
% EXCHANGE_SLA subroutine to calculate SLA type exchange functional
%
%    See also xcRef.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


alpha = 2/3;
f = -9/8 * (3/2/pi)^(2/3);
falpha = f * alpha;
% falpha = -0.458165293283143;

vx = 4/3 * falpha ./ rs;
epsx = falpha ./ rs;

end