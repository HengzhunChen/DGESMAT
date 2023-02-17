function num_grid = NumGridTotal(domain)
% DOMAIN/NUMGRIDTOTAL NumGridTotal function for Domain class
%    
%    See also Domain.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


num_grid = domain.numGrid(1) * domain.numGrid(2) * domain.numGrid(3);

end