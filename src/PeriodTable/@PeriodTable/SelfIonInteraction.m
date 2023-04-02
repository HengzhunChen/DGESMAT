function Eself = SelfIonInteraction(PT, type)
% PERIODTABLE/SELFIONINERACTION computes the self ionic interaction energy
%    according to atom type (atomic number).
%  
%    See also PeriodTable, CalculateIonSelfEnergyAndForce.
%
%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.

RGaussian = PT.RGaussian(type);
Zval = PT.Zval(type);
Eself = Zval * Zval / ( sqrt(2*pi) * RGaussian );

end