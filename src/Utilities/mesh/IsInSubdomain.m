function flag = IsInSubdomain(r, sub_domain, length_global)
% ISINSUBDOMAIN Check whether a point is in a domain 
% 
%    flag = IsInSubdomain(r, domain, length_global) returns the flag that
%    indicates whether a point r is in the domain sub_domain, which is a 
%    part of a global domain with periodic boundary conditions and length
%    length_global.
%
%    See also Domain.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

flag = true;
posStart = sub_domain.posStart;
length_elem = sub_domain.length;

% shift for periodic boundary conditions
shift_start = mod( posStart, length_global );
shift_r = mod( r, length_global );

for i = 1 : dimDef()    
    % Case 1 of the buffer interval %
    if shift_start(i) + length_elem(i) > length_global(i)
        if shift_r(i) > shift_start(i) + length_elem(i) - length_global(i) && ...
           shift_r(i) < shift_start(i)
            flag = false;
        end
    else    
    % Case 2 of the buffer interval %
        if shift_r(i) < shift_start(i) || ...
           shift_r(i) > shift_start(i) + length_elem(i)
            flag = false;
        end
    end
end

end