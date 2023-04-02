classdef ElementTable
    % ElementTable records basic information of most chemical elements in 
    % Periodic table.
    %
    %    See also Element, ReadUPF, PeriodTable.

    %  Copyright (c) 2022-2023 Hengzhun Chen and Yingzhou Li, 
    %                          Fudan University
    %  This file is distributed under the terms of the MIT License.
        
    properties (SetAccess = private)
        etable
        zmap
    end
    
    methods
        function elemTable = ElementTable()
            etable = Element.empty();
            etable(end+1) = Element(1,"H","1s1",1.00794);
            etable(end+1) = Element(2,"He","1s2",4.00260);
            etable(end+1) = Element(3, "Li","1s2 2s1",     6.941);
            etable(end+1) = Element(4, "Be","1s2 2s2",     9.01218);
            etable(end+1) = Element(5, "B", "1s2 2s2 2p1",10.811);
            etable(end+1) = Element(6, "C", "1s2 2s2 2p2",12.0107);
            etable(end+1) = Element(7, "N", "1s2 2s2 2p3",14.00674);
            etable(end+1) = Element(8, "O", "1s2 2s2 2p4",15.9994);
            etable(end+1) = Element(9, "F", "1s2 2s2 2p5",18.9884);
            etable(end+1) = Element(10,"Ne","1s2 2s2 2p6",20.1797);

            etable(end+1) = Element(11,"Na","[Ne] 3s1",    22.98977);
            etable(end+1) = Element(12,"Mg","[Ne] 3s2",    24.3050);
            etable(end+1) = Element(13,"Al","[Ne] 3s2 3p1",26.98154);
            etable(end+1) = Element(14,"Si","[Ne] 3s2 3p2",28.0855);
            etable(end+1) = Element(15,"P", "[Ne] 3s2 3p3",30.97376);
            etable(end+1) = Element(16,"S", "[Ne] 3s2 3p4",32.066);
            etable(end+1) = Element(17,"Cl","[Ne] 3s2 3p5",35.4527);
            etable(end+1) = Element(18,"Ar","[Ne] 3s2 3p6",39.948);

            etable(end+1) = Element(19,"K", "[Ar] 4s1",39.0983);
            etable(end+1) = Element(20,"Ca","[Ar] 4s2",40.078);
            etable(end+1) = Element(21,"Sc","[Ar] 3d1 4s2",44.95591);
            etable(end+1) = Element(22,"Ti","[Ar] 3d2 4s2",47.867);
            etable(end+1) = Element(23,"V", "[Ar] 3d3 4s2",50.9415);
            etable(end+1) = Element(24,"Cr","[Ar] 3d5 4s1",51.9961);
            etable(end+1) = Element(25,"Mn","[Ar] 3d5 4s2",54.93805);
            etable(end+1) = Element(26,"Fe","[Ar] 3d6 4s2",55.845);
            etable(end+1) = Element(27,"Co","[Ar] 3d7 4s2",58.9332);
            etable(end+1) = Element(28,"Ni","[Ar] 3d8 4s2",58.6934);
            etable(end+1) = Element(29,"Cu","[Ar] 3d10 4s1",63.546);
            etable(end+1) = Element(30,"Zn","[Ar] 3d10 4s2",65.39);
            etable(end+1) = Element(31,"Ga","[Ar] 3d10 4s2 4p1",69.723);
            etable(end+1) = Element(32,"Ge","[Ar] 3d10 4s2 4p2",72.61);
            etable(end+1) = Element(33,"As","[Ar] 3d10 4s2 4p3",74.9216);
            etable(end+1) = Element(34,"Se","[Ar] 3d10 4s2 4p4",78.96);
            etable(end+1) = Element(35,"Br","[Ar] 3d10 4s2 4p5",79.904);
            etable(end+1) = Element(36,"Kr","[Ar] 3d10 4s2 4p6",83.80);

            etable(end+1) = Element(37,"Rb","[Kr] 5s1",85.4678);
            etable(end+1) = Element(38,"Sr","[Kr] 5s2",87.62);
            etable(end+1) = Element(39,"Y" ,"[Kr] 4d1 5s2",88.90585);
            etable(end+1) = Element(40,"Zr","[Kr] 4d2 5s2",91.224);
            etable(end+1) = Element(41,"Nb","[Kr] 4d4 5s1",92.90638);
            etable(end+1) = Element(42,"Mo","[Kr] 4d5 5s1",95.94);
            etable(end+1) = Element(43,"Tc","[Kr] 4d5 5s2",98.0);
            etable(end+1) = Element(44,"Ru","[Kr] 4d7 5s1",101.07);
            etable(end+1) = Element(45,"Rh","[Kr] 4d8 5s1",102.9055);
            etable(end+1) = Element(46,"Pd","[Kr] 4d10",106.42);
            etable(end+1) = Element(47,"Ag","[Kr] 4d10 5s1",107.8682);
            etable(end+1) = Element(48,"Cd","[Kr] 4d10 5s2",112.411);
            etable(end+1) = Element(49,"In","[Kr] 4d10 5s2 5p1",114.818);
            etable(end+1) = Element(50,"Sn","[Kr] 4d10 5s2 5p2",118.710);
            etable(end+1) = Element(51,"Sb","[Kr] 4d10 5s2 5p3",121.760);
            etable(end+1) = Element(52,"Te","[Kr] 4d10 5s2 5p4",127.60);
            etable(end+1) = Element(53,"I" ,"[Kr] 4d10 5s2 5p5",126.90447);
            etable(end+1) = Element(54,"Xe","[Kr] 4d10 5s2 5p6",131.29);

            etable(end+1) = Element(55,"Cs","[Xe] 6s1",132.90545);
            etable(end+1) = Element(56,"Ba","[Xe] 6s2",137.327);
            etable(end+1) = Element(57,"La","[Xe] 5d1 6s2",138.9055);
            etable(end+1) = Element(58,"Ce","[Xe] 4f1 5d1 6s2",140.116);
            etable(end+1) = Element(59,"Pr","[Xe] 4f3 6s2",140.90765);
            etable(end+1) = Element(60,"Nd","[Xe] 4f4 6s2",144.24);
            etable(end+1) = Element(61,"Pm","[Xe] 4f5 6s2",145.0);
            etable(end+1) = Element(62,"Sm","[Xe] 4f6 6s2",150.36);
            etable(end+1) = Element(63,"Eu","[Xe] 4f7 6s2",151.964);
            etable(end+1) = Element(64,"Gd","[Xe] 4f7 5d1 6s2",157.25);
            etable(end+1) = Element(65,"Tb","[Xe] 4f9 6s2",158.92534);
            etable(end+1) = Element(66,"Dy","[Xe] 4f10 6s2",162.50);
            etable(end+1) = Element(67,"Ho","[Xe] 4f11 6s2",164.93032);
            etable(end+1) = Element(68,"Er","[Xe] 4f12 6s2",167.26);
            etable(end+1) = Element(69,"Tm","[Xe] 4f13 6s2",168.93421);
            etable(end+1) = Element(70,"Yb","[Xe] 4f14 6s2",173.04);
            etable(end+1) = Element(71,"Lu","[Xe] 4f14 5d1 6s2",174.967);
            etable(end+1) = Element(72,"Hf","[Xe] 4f14 5d2 6s2",178.49);
            etable(end+1) = Element(73,"Ta","[Xe] 4f14 5d3 6s2",180.9479);
            etable(end+1) = Element(74,"W" ,"[Xe] 4f14 5d4 6s2",183.84);
            etable(end+1) = Element(75,"Re","[Xe] 4f14 5d5 6s2",186.207);
            etable(end+1) = Element(76,"Os","[Xe] 4f14 5d6 6s2",190.23);
            etable(end+1) = Element(77,"Ir","[Xe] 4f14 5d7 6s2",192.217);
            etable(end+1) = Element(78,"Pt","[Xe] 4f14 5d9 6s1",195.078);
            etable(end+1) = Element(79,"Au","[Xe] 4f14 5d10 6s1",196.96655);
            etable(end+1) = Element(80,"Hg","[Xe] 4f14 5d10 6s2",200.59);
            etable(end+1) = Element(81,"Tl","[Xe] 4f14 5d10 6s2 6p1",204.3833);
            etable(end+1) = Element(82,"Pb","[Xe] 4f14 5d10 6s2 6p2",207.2);
            etable(end+1) = Element(83,"Bi","[Xe] 4f14 5d10 6s2 6p3",208.98038);
            etable(end+1) = Element(84,"Po","[Xe] 4f14 5d10 6s2 6p4",209.0);
            etable(end+1) = Element(85,"At","[Xe] 4f14 5d10 6s2 6p5",210.0);
            etable(end+1) = Element(86,"Rn","[Xe] 4f14 5d10 6s2 6p6",222.0);

            etable(end+1) = Element(87,"Fr","[Rn] 7s1",223.0);
            etable(end+1) = Element(88,"Ra","[Rn] 7s2",226.0);
            etable(end+1) = Element(89,"Ac","[Rn] 6d1 7s2",227.0);
            etable(end+1) = Element(90,"Th","[Rn] 6d2 7s2",232.0381);
            etable(end+1) = Element(91,"Pa","[Rn] 5f2 6d1 7s2",231.03588);
            etable(end+1) = Element(92,"U" ,"[Rn] 5f3 6d1 7s2",238.0289);
            etable(end+1) = Element(93,"Np","[Rn] 5f4 6d1 7s2",237.0);
            etable(end+1) = Element(94,"Pu","[Rn] 5f5 6d1 7s2",244.0);

            zmap = containers.Map([etable.symbol], [etable.z]);
            elemTable.etable = etable;
            elemTable.zmap = zmap;
        end
        
        function atom_number = z(elemTable, symbol)
            atom_number = elemTable.zmap(symbol);
        end
        
        function atom_symbol = symbol(elemTable, z)
            atom_symbol = elemTable.etable(z).symbol;
        end
        
        function config = configuration(elemTable, varargin)
            if isnumeric(varargin{1})
                z = varargin{1};
                config = elemTable.etable(z).config;
            else
                symbol = varargin{1};
                idx = elemTable.zmap(symbol);
                config = elemTable.etable(idx).config;
            end
        end
        
        function m = mass(elemTable, varargin)
            if isnumeric(varargin{1})
                z = varargin{1};
                m = elemTable.etable(z).mass;
            else
                symbol = varargin{1};
                idx = elemTable.zmap(symbol);
                m = elemTable.etable(idx).mass;
            end
        end
                
        function s = size(elemTable)
            s = length(elemTable.etable);
        end

    end
    
end