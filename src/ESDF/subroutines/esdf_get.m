function out = esdf_get(label, default, dunit)
% ESDF_GET is used to setup the value of esdf parameters. It returns the 
%          value attached to the "label", separated into two types, 
%          numerical type in double, string type in string.
%
%    out = esdf_get(label, default) returns the value attached to "label"
%    from the input, if no value is assigned from input, returns the
%    default value.
%
%    out = esdf_get(label, default, dunit) returns the value attached to
%    "label" (physical quantity) with default unit "dunit".
%
%    See also ESDFReadInput, ESDFInputParam.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


global tokenList kw_label;

if nargin == 2  % without unit of physical quantities
    % check whether the label is defined
    % reduce the string to lower case and remove punctuation
    labelstr = esdf_reduce(label);
    index = find(kw_label == labelstr, 1);
    if isempty(index)
        error('Label "%s" not recognized in keyword list', label);
    end
    
    % set to default
    out = default;        
    
    % read corresponding value of label in inputFile if there exists
    token_idx = find(esdf_reduce(tokenList(:, 1)) == esdf_reduce(label));
    if ~isempty(token_idx)
        kw_value = tokenList(token_idx, 2);
        % if keyword value is numerical
        if isnumeric(default)
            kw_value = str2double(kw_value);
            if isnan(kw_value)
                error('Unable to parse "%s" in esdf_get', label);
            end
        end
        out = kw_value;
    end
    
elseif nargin == 3  % physical quantity with unit
    % check whether the label is defined
    % reduce the string to lower case and remove punctuation
    labelstr = esdf_reduce(label);
    
    index = find(kw_label == labelstr, 1);
    if isempty(index)
        error('Label "%s" not recognized in keyword list', label);
    end
    
    % set to default
    out = default;
    
    % read corresponding value of label in inputFile if there exists
    token_idx = find(esdf_reduce(tokenList(:, 1)) == esdf_reduce(label));
    if ~isempty(token_idx)
        kw_value = tokenList(token_idx(1), 2);
        % if keyword value is numerical
        if isnumeric(default)
            kw_value = str2double(kw_value);
            if isnan(kw_value)
                error('Unable to parse "%s" in esdf_get', label);
            end
        end
        out = kw_value; 
        
        % deal with physical variable unit
        if length(tokenList(token_idx(1), :)) == 3  % line with unit
            iunit = esdf_reduce(tokenList(token_idx(1), 3));
            out = esdf_convfac(iunit, dunit) * out;
        end
    end
end

end