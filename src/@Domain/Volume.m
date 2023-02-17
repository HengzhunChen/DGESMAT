function volume = Volume(domain)
% DOMAIN/VOLUME Volume function for Domain class
%    
%    See also Domain.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


volume = domain.length(1) * domain.length(2) * domain.length(3);
    
end
