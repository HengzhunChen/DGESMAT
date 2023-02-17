function PrintState(scfDG)
% SCFDG/PRINTSTATE print state information for SCF iteration of DG method. 
%
%    See also SCFDG.

%  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
%                          Fudan University
%  This file is distributed under the terms of the MIT License.


hamDG = scfDG.hamDG;

InfoPrint(0, '\n Eigenvalues in the global domain.\n');
for i = 1 : length(hamDG.eigVal)
    InfoPrint(0, ... 
        "band#    = ", i, ...
        "eigval   = ", hamDG.eigVal(i), ...
        "occrate  = ", hamDG.occupationRate(i));
end
InfoPrint(0, '\n');

InfoPrint(0, '\n');
InfoPrint(0, 'NOTE:  Ecor  = Exc - EVxc - Ehart -Eself \n');
InfoPrint(0, '       Etot  = Ekin + Ecor \n');
InfoPrint(0, '       Efree = Etot + Entropy \n \n');

InfoPrint(0, "EfreeHarris       = ",  scfDG.EfreeHarris, "[au]");
InfoPrint(0, "EfreeSecondOrder  = ",  scfDG.EfreeSecondOrder, "[au]");
InfoPrint(0, "Etot              = ",  scfDG.Etot, "[au]");
InfoPrint(0, "Efree             = ",  scfDG.Efree, "[au]");
InfoPrint(0, "Ekin              = ",  scfDG.Ekin, "[au]");
InfoPrint(0, "Ehart             = ",  scfDG.Ehart, "[au]");
InfoPrint(0, "EVxc              = ",  scfDG.EVxc, "[au]");
InfoPrint(0, "Exc               = ",  scfDG.Exc, "[au]"); 
InfoPrint(0, "Evdw              = ",  scfDG.Evdw, "[au]"); 
InfoPrint(0, "Eself             = ",  scfDG.Eself, "[au]");
InfoPrint(0, "Ecor              = ",  scfDG.Ecor, "[au]");
InfoPrint(0, "Fermi             = ",  scfDG.fermi, "[au]");
    

end