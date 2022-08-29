function esdf_unit()
% ESDF_UNIT constructs factors for variations of the units.
%
%    Shorthands of physical units
%        m - mass 
%        l - length 
%        t - time 
%        e - energy 
%        f - force 
%        p - pressure
%        c - charge 
%        d - dipole 
%        mom - mom inert 
%        ef - efield
%
%    See also esdf_get.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


global phy_var phy_name phy_unit;

phy = {
    "m",   "kg",           1.0;
    "m",   "g",            1.0e-3;
    "m",   "amu",          1.66054e-27;
    "l",   "m",            1.0;
    "l",   "nm",           1.0e-9;
    "l",   "ang",          1.0e-10;
    "l",   "bohr",         0.529177e-10;
    "t",   "s",            1.0;
    "t",   "ns",           1.0e-9;
    "t",   "ps",           1.0e-12;
    "t",   "fs",           1.0e-15;
    "e",   "j",            1.0;
    "e",   "erg",          1.0e-7;
    "e",   "ev",           1.60219e-19;
    "e",   "mev",          1.60219e-22;
    "e",   "ry",           2.17991e-18;
    "e",   "mry",          2.17991e-21;
    "e",   "hartree",      4.35982e-18;
    "e",   "kcal/mol",     6.94780e-21;
    "e",   "mhartree",     4.35982e-21;
    "e",   "kj/mol",       1.6606e-21;
    "e",   "hz",           6.6262e-34;
    "e",   "thz",          6.6262e-22;
    "e",   "cm-1",         1.986e-23;
    "e",   "cm^-1",        1.986e-23;
    "e",   "cm**-1",       1.986e-23;
    "f",   "N",            1.0;
    "f",   "ev/ang",       1.60219e-9;
    "f",   "ry/bohr",      4.11943e-8;
    "l",   "cm",           1.0e-2;
    "p",   "pa",           1.0;
    "p",   "mpa",          1.0e6;
    "p",   "gpa",          1.0e9;
    "p",   "atm",          1.01325e5;
    "p",   "bar",          1.0e5;
    "p",   "mbar",         1.0e11;
    "p",   "ry/bohr**3",   1.47108e13;
    "p",   "ev/ang**3",    1.60219e11;
    "c",   "c",            1.0;
    "c",   "e",            1.602177e-19;
    "d",   "C*m",          1.0;
    "d",   "D",            3.33564e-30;
    "d",   "debye",        3.33564e-30;
    "d",   "e*bohr",       8.47835e-30;
    "d",   "e*ang",        1.602177e-29;
    "mom", "kg*m**2",      1.0;
    "mom", "ry*fs**2",     2.1799e-48;
    "ef",  "v/m",          1.0;
    "ef",  "v/nm",         1.0e9;
    "ef",  "v/ang",        1.0e10;
    "ef",  "v/bohr",       1.8897268e10;
    "ef",  "ry/bohr/e",    2.5711273e11;
    "ef",  "har/bohr/e",   5.1422546e11;
    "e",   "k",            1.38066e-23;
    "b",   "t",            1.0;
    "b",   "ry/mu_bohr",   2.350499e5;
    "b",   "g",            1.0e4;
};

phy_var  = string(phy(:, 1));
phy_name = string(phy(:, 2));
phy_unit = cell2mat(phy(:, 3));

end