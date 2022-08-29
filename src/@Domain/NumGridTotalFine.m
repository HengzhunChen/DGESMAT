function num_grid = NumGridTotalFine(domain)
% DOMAIN/NUMGRIDTOTALFINE NumGridTotalFine function for Domain class
%    
%    See also Domain.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


num_grid = domain.numGridFine(1) * domain.numGridFine(2) * domain.numGridFine(3);

end