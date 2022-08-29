function t = Smoother(x)
% SMOOTHER Smoother function used for periodizing the potential in the 
%    extended element.
%
%    See also SCFDG/Setup.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


idx1 = x <= 0;
idx2 = (0 < x) & (x < 1);
idx3 = x >= 1;

t = x;
t(idx1) = 1;
t(idx3) = 0;

z = x;
z(idx2) = -1 ./ x(idx2) + 1 ./ (1 - x(idx2));
idxz1 = (z < 0) & idx2;
idxz2 = (z >= 0) & idx2;

t(idxz1) = 1 ./ ( exp(z(idxz1)) + 1 );
t(idxz2) = exp(-z(idxz2)) ./ ( exp(-z(idxz2)) + 1 );

% used for reference
% if x <= 0
%     t = 1.0;
% elseif x >= 1
%     t = 0.0;
% else
%     z = -1 / x + 1 / (1 - x);
%     if z < 0
%         t = 1 / (exp(z) + 1);
%     else
%         t = exp(-z) / (exp(-z) + 1);
%     end
% end


end