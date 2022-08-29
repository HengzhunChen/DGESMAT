function factor = esdf_convfac(from, to)
% ESDF_CONVFAC finds the conversion factor between physical units.
% 
%    See also esdf_unit, esdf_get.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global phy_var phy_name phy_unit

% find index of the from and to units
phy_name = esdf_reduce(phy_name);
idx_from = find(phy_name == esdf_reduce(from));
idx_to   = find(phy_name == esdf_reduce(to));

% check the units are recognized
if isempty(idx_from)
   error('Units not recognized in input file : %s', esdf_reduce(from));
end
if isempty(idx_to)
    error('Units not recognized in Program : %s', esdf_reduce(to));
end

% check from and to units are of the same physical variable
if phy_var(idx_from) ~= phy_var(idx_to)
    error('Dimensions do not match : %s vs %s', esdf_reduce(from), esdf_reduce(to));
end

% set the converiosn factor
factor = phy_unit(idx_from) / phy_unit(idx_to);
    
end