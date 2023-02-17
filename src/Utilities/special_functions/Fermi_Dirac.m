function f = Fermi_Dirac(eigVal, Efermi, Tbeta)
% FERMI_DIRAC calculates Fermi-Dirac distribution.
%    
%    f = FERMIDIRAC(eigVal, Efermi, Tbeta) returns the Fermi-Dirac 
%    distribution with Fermi energy Efermi, inverse of temperature in 
%    atomic unit Tbeta at energy eigVal.
%
%    The Fermi-Dirac distribution is defined as follows,
%       f(e) = \frac{1}{ exp{ (e-e_{Fermi}) * T_{beta} } + 1 }.
%

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

f = 1 ./ ( 1 + exp( (eigVal - Efermi) * Tbeta ) );

end