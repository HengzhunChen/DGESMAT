function [out_r, out_v] = spline_rad(r, v, even)
% SPLINE_RAD Interpolating odd/even functions along the radius direction
%
%    [out_r, out_v] = spline_rad(r, v, even) interpolates radial grid r
%    with value over grid v according to flag even. r Must be non-negative 
%    and can contain zero. even == 1 for even function and even == 0 for 
%    odd function. After interpolation, new grid out_r contains the same 
%    number of positive and negative points and does not contain zero, and 
%    try to avoid the singular behavior near r = 0. out_v is the value over 
%    new grid out_r.
%
%    See also spline.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% avoid r=0 and extrapolation beyond the given grid
r_positive = r(r > 0); 
rmin = min(r_positive);  
rmax = max(r_positive);
dstep = 0.001;

% extend the positive grid to negative part 
% [-reverse(rtemp), 0, rtemp]
rtemp = rmin : dstep : rmax;  % row vector
rtemp = rtemp(1 : end-1);
out_r = [-fliplr(rtemp), rtemp];

% interpolation by cubic spline
[slpb, slpc, slpd] = spline(r, v);
vtemp = seval(rtemp, r, v, slpb, slpc, slpd);

vtemp = reshape(vtemp, 1, []);
if even
    out_v = fliplr(vtemp);
else
    out_v = - fliplr(vtemp);
end

out_v = [out_v, vtemp];


end