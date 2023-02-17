function reduced_label = esdf_reduce(label)
% ESDF_REDUCE returns the reduced form of the input label.
%    Reduce the string to lower case and remove punctuation.
%
%    Note: Only global variable tokenList will be reduced, lineList will
%    not be reduced since it contains data separated by punctuation.
%    
%    See also esdf_init, esdf_get.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

labelstr = lower(label);
labelstr = erase(labelstr, '-');
labelstr = erase(labelstr, '.');
reduced_label = labelstr;

end