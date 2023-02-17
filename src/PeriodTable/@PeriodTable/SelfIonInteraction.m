function Eself = SelfIonInteraction(PT, type)
% PERIODTABLE/SELFIONINERACTION computes the self ionic interaction energy
%    according to atom type (atomic number).
%  
%    See also PeriodTable, HamiltonianKS/CalculateIonSelfEnergyAndForce,
%    Atom.
%
%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


if ~PT.userOption.isUseVLocal
    Eself = PT.pteMap(type).params.Eself;
else
    RGaussian = PT.RGaussian(type);
    Zion = PT.Zion(type);
    Eself = Zion * Zion / ( sqrt(2*pi) * RGaussian );
end

end