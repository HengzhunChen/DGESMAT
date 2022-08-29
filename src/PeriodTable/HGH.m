function tempEntry = HGH(Znuc)
% HGH generates HGH type pseudopotential.
%
%    tempEntry = HGH(Zunc) produces HGH pseudopotential according to atomic
%    number Znuc and stores in the structure of PTEntry.
%
%    NOTE: 
%    To generate more types of HGH with different exchange-correlation
%    functionals, use information from
%    https://www.cp2k.org/static/potentials
%
%    See also PeriodTable, PTEntry.

%  Copyright (c) 2022 Hengzhun Chen and Yingzhou Li, 
%                     Fudan University
%  This file is distributed under the terms of the MIT License.


% TODO: some elements in periodic table not being supported here

InfoPrint([0, 1], 'Use HGH type pseudopotential');
InfoPrint([0, 1], 'Znuc %d\n', Znuc);


% ----------------------------------------------------------------------
% Initialize  
% ----------------------------------------------------------------------

% IMPORTANT: Everything is zero but r variables initialized to 1
% to avoid overflow /underflow 
Zion = 0;
mass = 0;

% Local pseudopotential
rloc = 1;
C1   = 0;
C2   = 0;
C3   = 0;
C4   = 0;

% Regular nonlocal pseudopotential
r0   = 1;
h011 = 0;
h022 = 0;
h033 = 0;
r1   = 1;
h111 = 0;
h122 = 0;
h133 = 0;
r2   = 1;
h211 = 0;
h222 = 0;
h233 = 0;
r3   = 1;
h311 = 0;

% Spin-orbit coupling
k111 = 0;
k122 = 0;
k133 = 0;
k211 = 0;
k222 = 0;
k233 = 0; 
k311 = 0;

% New element template
% if(Znuc==)
% Zion = ;
% mass = ;
% rloc    = ;
% C1      = ;
% r0      = ;
% h011    = ;
% h022    = ;
% r1      = ;
% h111    = ;
% h122    = ;
% k111    = ;
% k122    = ;
% k133    = ;
% r2      = ;
% h211    = ;
% h222    = ;
% k211    = ;
% k222    = ;
%
% rhocut  = ;
% wavcut  = ;
% end 


% ----------------------------------------------------------------------
% Detail info for each element in PeriodicTable
% ----------------------------------------------------------------------

% H
if(Znuc==1)
	Zion = 1;
	mass = 1.00794;
	rloc    = 0.200000;
	C1      = -4.180237;
	C2      = 0.725075;
	%
	rhocut  = 2.0;
	wavcut  = 2.0;
end 

% Li (semicore)
if(Znuc==3)
	Zion = 3;
	mass = 6.941;
	rloc    = 0.400000;
	C1      = -14.034868;
	C2      = 9.553476;
	C3      = -1.766488;
	C4      = 0.084370;
	%
	rhocut  = 3.5;
	wavcut  = 3.5;
end 

% Li (no semicore)
% if(Znuc==3)
	% Zion = 1;
	% mass = 6.941;
	% rloc    = 0.787553;
	% C1      = -1.892612;
	% C2      = 0.286060;
	% r0      = 0.666375;
	% h011    = 1.858811;
	% r1      = 1.079306;
	% h111    = -0.005895;
	% k111    = 0.000019;  
	% %
	% rhocut  = 7.5;
	% wavcut  = 7.5;
% end 

% B
if(Znuc==5)
	Zion = 3;
	mass = 10.81; 
	rloc = 0.433930;
	C1 = -5.578642;
	C2 = 0.804251;
	r0 = 0.373843;
	h011 = 6.233928;
	r1 = 0.360393;
	k111 = 0.000878;  
	rhocut = 3.5;
	wavcut = 3.0;
end

% C
if(Znuc==6)
	Zion = 4;
	mass = 12.011;
	rloc = 0.348830;
	C1 = -8.513771;
	C2 = 1.228432;
	r0 = 0.304553;
	h011 = 9.522842;
	r1 = 0.232677;
	k111 = 0.004104;  
	rhocut = 3.0;
	wavcut = 3.0;
end

% N
if(Znuc==7)
	Zion = 5;
	mass = 14.007;
	rloc = 0.289179;
	C1 = -12.234820;
	C2 = 1.766407;
	r0 = 0.256605;
	h011 = 13.552243;
	r1 = 0.270134;
	k111 = 0.003131;  
	rhocut = 2.5;
	wavcut = 2.5;
end

% O. New
if(Znuc==8)
	Zion = 6;
	mass = 15.9994;
	rloc    = 0.247621;
	C1      = -16.580318;
	C2      = 2.395701;
	r0      = 0.221786;
	h011    = 18.266917;
	r1      = 0.256829;
	k111    = 0.004476;
	%
	rhocut  = 2.0;
	wavcut  = 2.0;
end 

% F
if(Znuc==9)
	Zion = 7;
	mass = 18.9984032;
	rloc    = 0.218525;
	C1      = -21.307361;
	C2      = 3.072869;
	r0      = 0.195567;
	h011    = 23.584942;
	r1      = 0.174268;
	k111    = 0.015106;

	rhocut = 2.0;
	wavcut = 2.0;
end 

% Na (No semicore)
if(Znuc==11)
	Zion = 1;
	mass = 22.989768;
	rloc = 0.885509;
	C1 = -1.238867;
	r0 = 0.661104;
	h011 = 1.847271;
	h022 = 0.582004;
	r1 = 0.857119;
	h111 = 0.471133;
	k111 = 0.002623;
	rhocut = 6.5;
	wavcut = 6.5;
end

% Al 
if(Znuc==13)
	Zion = 3;
	mass = 26.9815386;
	rloc = 0.450000;
	C1   = -8.491351;
	r0   = 0.460104;
	h011 = 5.088340;
	h022 = 2.679700;
	r1   = 0.536744;
	h111 = 2.193438;
	k111 = 0.006154;
	k122 = 0.003947;
	rhocut = 3.5;
	wavcut = 3.5;
end

% Al: local potential only
% if(Znuc==13)
	% Zion = 3;
	% mass = 26.9815386;
	% rloc = 0.450000;
	% C1   = -8.491351;
	% r0   = 0.460104;
	% rhocut = 3.5;
	% wavcut = 3.5;
% end

%Si
if(Znuc==14)
	Zion = 4;
	mass = 28.0855;
	rloc = 0.440000;
	C1 = -7.336103;
	r0 = 0.422738;
	h011 = 5.906928;
	h022 = 3.258196;
	r1 = 0.484278;
	h111 = 2.727013;
	k111 = 0.000373;
	k122 = 0.014437;
	rhocut = 3.5;
	wavcut = 3.5;
end

% P
if(Znuc==15) 
	Zion = 5;
	mass = 30.973762;
	rloc = 0.430000; 
	C1 = -6.654220;
	r0 = 0.389803;
	h011 = 6.842136;
	h022 = 3.856693;
	r1 = 0.440796;
	h111 = 3.282606;
	k111 = 0.002544;
	k122 = 0.017895;
	rhocut = 3.5;
	wavcut = 3.5;
end

% Ti (WITH semicore)
if(Znuc==22) 
	Zion = 12;
	mass = 47.867;
	rloc = 0.380000; 
	C1 = 7.548789;
C2 = -0.588377;
	r0 = 0.334235;
	h011 = 6.925740;
	h022 = -3.142005;
	r1 = 0.242416;
	h111 = 5.079086;
h122 = -6.284281;
	k111 = 0.122395;
	k122 = 0.057447;
	r2      = 0.242947;
	h211    = -9.125896;
	k211    = 0.005822;
	rhocut = 3.5;
	wavcut = 3.5;
end

% Copper (With semicore) 
if(Znuc==29)
    Zion = 11;
    mass = 63.546;
    rloc = 0.530000;

    r0   = 0.423734;
    h011 = 3.888050;
    h022 = 3.276584;
    h033 = 2.290091;

    r1   =  0.572177;
    h111 =  1.751272;
    h122 =  0.374943;
    k111 = -0.024067;
    k122 =  0.076481;

    r2   =  0.266143;
    h211 = -12.676957;
    k211 =  0.010489;

    rhocut = 3.5;
    wavcut = 3.5;
end

% Ga
if(Znuc==31)
	Zion    = 3;
	mass    = 69.723;
	rloc    = 0.560000;
	r0      = 0.610791;
	h011    = 2.369325;
	h022    = -0.249015;
	h033    = -0.551796;
	r1      = 0.704596;
	h111    = 0.746305;
	h122    = -0.513132;
	k111    = 0.029607;
	k122    = -0.000873;
	r2      = 0.982580;
	h211    = 0.075437;
	k211    = 0.001486;
	rhocut = 4.0;
	wavcut = 5.5;
end

% As 
if(Znuc==33)
	Zion    = 5;
	mass    = 74.92160;
	rloc    = 0.520000;
	r0      = 0.456400;
	h011    = 4.560761;
	h022    = 1.692389;
	h033    = -1.373804;
	r1      = 0.550562;
	h111    = 1.812247;
	h122    = -0.646727;
	k111    = 0.052466;
	k122    = 0.020562;
	r2      = 0.685283;
	h211    = 0.312373;
	k211    = 0.004273;
	rhocut = 4.0;
	wavcut = 4.0;
end

% Se
if(Znuc==34)
	Zion    = 6;
	mass    = 78.96;
	rloc    = 0.510000;
	r0      = 0.432531;
	h011    = 5.145131;
	h022    = 2.052009;
	h033    = -1.369203;
	r1      = 0.472473;
	h111    = 2.858806;
	h122    = -0.590671;
	k111    = 0.062196;
	k122    = 0.064907;
	r2      = 0.613420;
	h211    = 0.434829;
	k211    = 0.005784;
	rhocut = 3.0;
	wavcut = 3.0;
end

% Sr (WITH semicore)
if(Znuc==38)
	Zion    = 10;
	mass    = 87.62;
	rloc    = 0.480000;
	C1      = 5.571455;
	C2      = -1.079963;
	r0      = 0.275441;
	h011    = 9.995135;
	h022    = 9.336679;
	r1      = 0.302243;
	h111    = 3.169126;
	h122    = 4.049231;
	k111    = -0.576265;
	k122    = 0.990062;
	r2      = 0.502045;
	h211    = 0.438728;
	k211    = 0.008991;
	rhocut = 4.0;
	wavcut = 3.0;
end

% La (WITH semi-core)
if(Znuc==57)
	Zion = 11;
	mass = 138.90;
	rloc    = 0.535000;
	C1      = 19.909308;
	C2      = -1.474830;
	r0      = 0.551775;
	h011    = 1.293272;
	h022    = -1.121819;
	r1      = 0.476308;
	h111    = 1.172527;
	h122    = -0.828810;
    h133    = 0.029857;
	k111    = 0.524623;
	k122    = -0.030901;
	k133    = 0.142077;
	r2      = 0.626672;
	h211    = 0.328377;
	k211    = 0.020900;
	r3      = 0.299310;
	h311    = -18.269439;
	k311    = 0.007193;
	%
	rhocut  = 4.5;
	wavcut  = 3.5;
end 

% Ce
if(Znuc==58)
	Zion    = 12;
	mass    = 140.116;
	rloc    = 0.535000;
	C1      = 18.847470;
	C2      = -0.765636;
	r0      = 0.521790;
	h011    = 1.321616;
	h022    = -1.700444;
	r1      = 0.470324;
	h111    = 0.972641;
	h122    = -1.451337;
	h133    = 0.000000;
	k111    = 0.463710;
	k122    = 0.090257;
	k133    = 0.012566;
	r2      = 0.703593;
	h211    = 0.074241;
	k211    = 0.013265;
	r3      = 0.306717;
	h311    = -17.214790;
	k311    = 0.007568;
	%
	rhocut  = 4.0;
	wavcut  = 3.5;
end 

% Yb (WITH semi-core)
if(Znuc==70)
	Zion = 24;
	mass = 173.054;
	rloc    = 0.500000;
	C1      = 17.357144;
	C2      = -1.773916;
	r0      = 0.402309;
	h011    = 2.120771;
	h022    = -4.802990;
	r1      = 0.414358;
	h111    = -0.923212;
	h122    = -1.678540;
	k111    = -0.349470;
	k122    = 2.074295;
	k133    = -1.814297;
	r2      = 0.444025;
	h211    = -0.889967;
	k211    = 0.070076;
	r3      = 0.238298;
	h311    = -29.932854;
	k311    = 0.025718;
	%
	rhocut  = 4.0;
	wavcut  = 3.0;
end 

% Lu (WITH semi-core)

if(Znuc==71)
	Zion = 25;
	mass = 174.9668;
	rloc    = 0.497000;
	C1      = 17.037053;
	C2      = -1.661610;
	r0      = 0.391206;
	h011    = 2.184678;
	h022    = -5.432346;
	r1      = 0.393896;
	h111    = -0.719819;
	h122    = -2.723799;
	k111    = 0.152450;
	k122    = 1.395416;
	k133    = -1.238744;
	r2      = 0.436518;
	h211    = -1.173245;
	k211    = 0.072130;
	r3      = 0.232629;
	h311    = -31.852262;
	k311    = 0.028006;
	%
	rhocut  = 3.5;
	wavcut  = 3.0;
end 

% Pt (WITHOUT Semi-core)
% if(Znuc==78)
% 	Zion    = 10;
% 	mass    = 195.084;
% 	rloc    = 0.616000;
% 	C1      = 11.027417;
% 	r0      = 0.520132;
% 	h011    = 2.447430;
% 	h022    = 2.640360;
% 	r1      = 0.658976;
% 	h111    = 0.408453;
% 	h122    = 1.647716;
% 	k111    = -0.763296;
% 	k122    = 1.065883;
% 	r2      = 0.451243;
% 	h211    = -4.552295;
% 	h222    = -2.102396;
% 	k211    = 0.146912;
% 	k222    = -0.169306;
% 
% 	rhocut = 4.0;
% 	wavcut = 3.5;
% end

% Pt (WITH Semi-core)
if(Znuc==78)
	Zion    = 18;
	mass    = 195.084;
	rloc    = 0.500000;
	C1      = 5.445832;
	C2      = 1.156382;
	r0      = 0.409942;
	h011    = 2.994366;
	h022    = -7.448772;
	h033    = 4.243095;
	r1      = 0.398652;
	h111    = -0.225181;
	h122    = -3.776974;
	k111    = 1.017060;
	k122    = -0.348213;
	k133    = -0.331919;
	r2      = 0.367964;
	h211    = 0.632067;
	h222    = -5.755431;
	k211    = 0.226472;
	k222    = -0.114346;

	rhocut = 3.5;
	wavcut = 3.5;
end

% Au (WITH Semi-core)
if(Znuc==79)
	Zion    = 11;
	mass    = 196.966569;
	rloc    = 0.590000;
	C1      = 11.604428;
	r0      = 0.521180;
	h011    = 2.538614;
	h022    = 2.701113;
	r1      = 0.630613;
	h111    = 0.394853;
	h122    = 2.057831;
	k111    = -0.960055;
	k122    = 1.296571;
	r2      = 0.440706; 
	h211    = -4.719070;
	h222    = -1.650429;
	k211    = 0.148484;
	k222    = -0.169493;

	rhocut = 3.5;
	wavcut = 3.5;
end

% Tl (WITH semi-core)
if(Znuc==81)
	Zion    = 13;
	mass    = 204.3833;
	rloc    = 0.550000;
	C1      = 7.301886;
	r0      = 0.502423;
	h011    = 3.326560;
	h022    = 4.341390;
	r1      = 0.572016;
	h111    = 1.272807;
	h122    = 2.992206;
	k111    = 0.012233;
	k122    = 0.031664;
	k133    = 1.019164;
	r2      = 0.393185;
	h211    = -3.200652;
	h222    = -3.008296;
	k211    = 0.186849;
	k222    = -0.170651;

	rhocut = 3.5;
	wavcut = 3.5;
end

% Hg (WITH Semi-core)
if(Znuc==80)
	Zion    = 12;
	mass    = 200.59;
	rloc    = 0.570000;
	C1      = 2.134572;
	r0      = 0.521802;
	h011    = 3.293920;
	h022    = 4.661001;
	r1      = 0.621648;
	h111    = 2.100960;
	h122    = 1.689988;
	k111    = 0.084989;
	k122    = 0.072771;
	k133    = 0.653348;
	r2      = 0.401894;
	h211    = -1.669886;
	h222    = -2.473265;
	k211    = 0.155759;
	k222    = -0.122282;

	rhocut = 3.5;
	wavcut = 3.5;
end

% Bi
if(Znuc==83)
	Zion    = 5;
	mass    = 208.98040;
	rloc    = 0.605000;
	C1      = 6.679437;
	r0      = 0.678858;
	h011    = 1.377634;
	h022    = -0.513697;
	h033    = -0.471028;
	r1      = 0.798673;
	h111    = 0.655578;
	h122    = -0.402932;
	k111    = 0.305314;
	k122    = -0.023134;
	r2      = 0.934683;
	h211    = 0.378476;
	k211    = 0.029217;

	rhocut = 4.5;
	wavcut = 4.5;
end


% Derived quantities
h012 = -1/2*sqrt(3/5)*h022;
h013 = 1/2 * sqrt(5/21) * h033;
h023 = -1/2 * sqrt(100/63) * h033;
h112 = -1/2 * sqrt(5/7) * h122;
h113 = 1/6 * sqrt(35/11) * h133;
h123 = -1/6 * 14 / sqrt(11) * h133;
h212 = -1/2 * sqrt(7/9) * h222;
h213 = 1/2  * sqrt(63/143) * h233;
h223 = -1/2 * 18 / sqrt(143) * h233;

k112 = -1/2 * sqrt(5/7) * k122;
k113 = 1/6 * sqrt(35/11) * k133;
k123 = -1/6 * 14 / sqrt(11) * k133;
k212 = -1/2 * sqrt(7/9) * k222;
k213 = 1/2  * sqrt(63/143) * k233;
k223 = -1/2 * 18 / sqrt(143) * k233;


% ----------------------------------------------------------------------
% Functional form for pseudopotentials
% ----------------------------------------------------------------------

Vlocfn = @(r) ...
	(-Zion./r.*erf(r./(sqrt(2)*rloc)) + ...
	exp(-1/2*(r./rloc).^2).*(C1+C2*(r/rloc).^2+C3*(r/ ...
	rloc).^4+C4*(r./rloc).^6));

rho0fn = @(r) ...
	(1/(4*pi).*(-sqrt(2/pi)*Zion./rloc.^3.*exp(-1/2*(r./rloc).^2) ...
	+ exp(-1/2*(r./rloc).^2).* ...
	((-6*C2+3*C1)/rloc.^2 + (-20*C3+7*C2-C1).*r.^2/rloc.^4 ...
	+ (-42*C4+11*C3-C2).*r.^4/rloc.^6 + ...
	(15*C4-C3).*r.^6/rloc.^8 - C4.*r.^8./rloc.^10)));

drho0fn = @(r) ...
	(1/(4*pi).*(sqrt(2/pi)*Zion*r./rloc.^5.*exp(-1/2*(r./rloc).^2) ...
	+ exp(-1/2*(r./rloc).^2).* ...
	((-40*C3+20*C2-5*C1)*r/rloc.^4 + ...
	(-168*C4+64*C3-11*C2+C1).*r.^3/rloc.^6 + ...
	(132*C4-17*C3+C2).*r.^5/rloc.^8 + ...
	(-23*C4+C3).*r.^7/rloc.^10 + C4.*r.^9./ ...
	rloc.^12)));


l=0; i=1; rl=r0;
p01fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=0; i=2; rl=r0;
p02fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=0; i=3; rl=r0;
p03fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=1; i=1; rl=r1;
p11fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=1; i=2; rl=r1;
p12fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=1; i=3; rl=r1;
p13fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=2; i=1; rl=r2;
p21fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=2; i=2; rl=r2;
p22fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=2; i=3; rl=r2;
p23fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));

l=3; i=1; rl=r3;
p31fn = @(r) ...
	(sqrt(2) * r.^(l+2*(i-1)).*exp(-r.^2/(2*rl^2))) ./...
	(rl.^(l+(4*i-1)/2).*sqrt(gamma(l+(4*i-1)/2)));


% Diagonalize the h coefficient matrices to obtain regular pseudopotential
% l=0 channel
tmp = [h011 h012 h013;
       h012 h022 h023;
       h013 h023 h033];
[V0,D0] = eig(tmp);
D0 = diag(D0);  [D0, idx] = sort(D0,'descend'); V0 = V0(:,idx);

p0afn = @(r) (p01fn(r)*V0(1,1) + p02fn(r)*V0(2,1) + p03fn(r)*V0(3,1));
h0a = D0(1);

p0bfn = @(r) (p01fn(r)*V0(1,2) + p02fn(r)*V0(2,2) + p03fn(r)*V0(3,2));
h0b = D0(2);

p0cfn = @(r) (p01fn(r)*V0(1,3) + p02fn(r)*V0(2,3) + p03fn(r)*V0(3,3));
h0c = D0(3);

% l=1 channel
tmp = [h111 h112 h113;
       h112 h122 h123;
       h113 h123 h133];
[V1,D1] = eig(tmp);
D1 = diag(D1);  [D1, idx] = sort(D1,'descend'); V1 = V1(:,idx);

p1afn = @(r) (p11fn(r)*V1(1,1) + p12fn(r)*V1(2,1) + p13fn(r)*V1(3,1));
h1a   = D1(1);

p1bfn = @(r) (p11fn(r)*V1(1,2) + p12fn(r)*V1(2,2) + p13fn(r)*V1(3,2));
h1b   = D1(2); 

p1cfn = @(r) (p11fn(r)*V1(1,3) + p12fn(r)*V1(2,3) + p13fn(r)*V1(3,3));
h1c   = D1(3);

% l=2 channel. 
tmp = [h211 h212 h213;
       h212 h222 h223;
       h213 h223 h233];

[V2,D2] = eig(tmp);
D2 = diag(D2);  [D2, idx] = sort(D2,'descend'); V2 = V2(:,idx);

p2afn = @(r) (p21fn(r)*V2(1,1) + p22fn(r)*V2(2,1) + p23fn(r)*V2(3,1));
h2a   = D2(1);

p2bfn = @(r) (p21fn(r)*V2(1,2) + p22fn(r)*V2(2,2) + p23fn(r)*V2(3,2));
h2b   = D2(2);

p2cfn = @(r) (p21fn(r)*V2(1,3) + p22fn(r)*V2(2,3) + p23fn(r)*V2(3,3));
h2c   = D2(3); 

% l=3 is special
p3afn = @(r) p31fn(r);
h3a   = h311; 


% Diagonalize the k coefficient matrices to obtain pseudopotential for
% spin-orbit coupling

% l=1 channel
tmp = [k111 k112 k113;
       k112 k122 k123;
       k113 k123 k133];
[Vk1,Dk1] = eig(tmp);
Dk1 = diag(Dk1);  [Dk1, idx] = sort(Dk1,'descend'); Vk1 = Vk1(:,idx);

p1kafn = @(r) (p11fn(r)*Vk1(1,1) + p12fn(r)*Vk1(2,1) + p13fn(r)*Vk1(3,1));
k1a   = Dk1(1);

p1kbfn = @(r) (p11fn(r)*Vk1(1,2) + p12fn(r)*Vk1(2,2) + p13fn(r)*Vk1(3,2));
k1b   = Dk1(2); 

p1kcfn = @(r) (p11fn(r)*Vk1(1,3) + p12fn(r)*Vk1(2,3) + p13fn(r)*Vk1(3,3));
k1c   = Dk1(3);

% l=2 channel. 
tmp = [k211 k212 k213;
       k212 k222 k223;
       k213 k223 k233];

[Vk2,Dk2] = eig(tmp);
Dk2 = diag(Dk2);  [Dk2, idx] = sort(Dk2,'descend'); Vk2 = Vk2(:,idx);

p2kafn = @(r) (p21fn(r)*Vk2(1,1) + p22fn(r)*Vk2(2,1) + p23fn(r)*Vk2(3,1));
k2a   = Dk2(1);

p2kbfn = @(r) (p21fn(r)*Vk2(1,2) + p22fn(r)*Vk2(2,2) + p23fn(r)*Vk2(3,2));
k2b   = Dk2(2);

p2kcfn = @(r) (p21fn(r)*Vk2(1,3) + p22fn(r)*Vk2(2,3) + p23fn(r)*Vk2(3,3));
k2c    = Dk2(3); 

% l=3 is special
p3kafn = @(r) p31fn(r);
k3a    = k311; 


% ---------------------------------------------------------------------
% sample
% ---------------------------------------------------------------------

stp = 0.0005;
r = (-10 + stp/2) : stp : 10;

Vloc = Vlocfn(r);
p0a = p0afn(r);
p0b = p0bfn(r);
p0c = p0cfn(r);
p1a = p1afn(r);
p1b = p1bfn(r);
p1c = p1cfn(r);
p2a = p2afn(r);
p2b = p2bfn(r);
p2c = p2cfn(r);
p3a = p3afn(r);

p1ka = p1kafn(r);
p1kb = p1kbfn(r);
p1kc = p1kcfn(r);
p2ka = p2kafn(r);
p2kb = p2kbfn(r);
p2kc = p2kcfn(r);
p3ka = p3kafn(r);

p0a_pp = csape(r, p0a);
p0b_pp = csape(r, p0b);
p0c_pp = csape(r, p0c);
p1a_pp = csape(r, p1a);
p1b_pp = csape(r, p1b);
p1c_pp = csape(r, p1c);
p2a_pp = csape(r, p2a);
p2b_pp = csape(r, p2b);
p2c_pp = csape(r, p2c);
p3a_pp = csape(r, p3a);

p1ka_pp = csape(r, p1ka);
p1kb_pp = csape(r, p1kb);
p1kc_pp = csape(r, p1kc);
p2ka_pp = csape(r, p2ka);
p2kb_pp = csape(r, p2kb);
p2kc_pp = csape(r, p2kc);
p3ka_pp = csape(r, p3ka);


% Pseudo-charge and its derivatives
rho = rho0fn(r); 
drho = drho0fn(r); 

% Regular nonlocal pseudopotential and their derivatives
p0a_dpp = fnder(p0a_pp, 1);
dp0a = fnval(p0a_dpp, r);

p0b_dpp = fnder(p0b_pp, 1);
dp0b = fnval(p0b_dpp, r);

p0c_dpp = fnder(p0c_pp, 1);
dp0c = fnval(p0c_dpp, r);

p1a_dpp = fnder(p1a_pp, 1);
dp1a = fnval(p1a_dpp, r);

p1b_dpp = fnder(p1b_pp, 1);
dp1b = fnval(p1b_dpp, r);

p1c_dpp = fnder(p1c_pp, 1);
dp1c = fnval(p1c_dpp, r);

p2a_dpp = fnder(p2a_pp, 1);
dp2a = fnval(p2a_dpp, r);

p2b_dpp = fnder(p2b_pp, 1);
dp2b = fnval(p2b_dpp, r);

p2c_dpp = fnder(p2c_pp, 1);
dp2c = fnval(p2c_dpp, r);

p3a_dpp = fnder(p3a_pp, 1);
dp3a = fnval(p3a_dpp, r);

% Nonlocal pseudopotential for spin-orbit coupling and their
% derivatives
p1ka_dpp = fnder(p1ka_pp, 1);
dp1ka = fnval(p1ka_dpp, r);

p1kb_dpp = fnder(p1kb_pp, 1);
dp1kb = fnval(p1kb_dpp, r);

p1kc_dpp = fnder(p1kc_pp, 1);
dp1kc = fnval(p1kc_dpp, r);

p2ka_dpp = fnder(p2ka_pp, 1);
dp2ka = fnval(p2ka_dpp, r);

p2kb_dpp = fnder(p2kb_pp, 1);
dp2kb = fnval(p2kb_dpp, r);

p2kc_dpp = fnder(p2kc_pp, 1);
dp2kc = fnval(p2kc_dpp, r);

p3ka_dpp = fnder(p3ka_pp, 1);
dp3ka = fnval(p3ka_dpp, r);

gd = find(r>0);
Es = 1/2* sum(4*pi*r(gd).^2 .* Vloc(gd) .* rho(gd) * stp);


% ---------------------------------------------------------------------
% save data into PTEntry() structure
% ---------------------------------------------------------------------

tempEntry = PTEntry();

tempEntry.params.Znuc = Znuc;
tempEntry.params.Mass = mass;
tempEntry.params.Zion = Zion;
tempEntry.params.Eself = Es;

tempEntry.samples.rad = r(:);
tempEntry.samples.pseudoCharge = -rho(:); 
tempEntry.samples.drv_pseudoCharge = -drho(:);

tempEntry.cutoffs.rad = rhocut;
tempEntry.cutoffs.pseudoCharge = rhocut;
tempEntry.cutoffs.drv_pseudoCharge = rhocut;

% nonlocal part
spl = zeros(numel(r), 0);
wgt = zeros(1,0);
typ = zeros(1,0);
cut = zeros(1,0);

cnt = 1;

if(abs(h0a) > 1e-10)
	spl(:, 2*cnt-1) = p0a(:);
	spl(:, 2*cnt) = dp0a(:);
    wgt(cnt) = h0a;  typ(cnt) = 0;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h0b) > 1e-10)
	spl(:, 2*cnt-1) = p0b(:);   
	spl(:, 2*cnt) = dp0b(:);  
 	wgt(cnt) = h0b;  typ(cnt) = 0;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h0c) > 1e-10)
	spl(:, 2*cnt-1) = p0c(:);   
	spl(:, 2*cnt) = dp0c(:);  
	wgt(cnt) = h0c;  typ(cnt) = 0;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h1a) > 1e-10)
	spl(:, 2*cnt-1) = p1a(:);   
	spl(:, 2*cnt) = dp1a(:);  
	wgt(cnt) = h1a;  typ(cnt) = 1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h1b) > 1e-10)
	spl(:, 2*cnt-1) = p1b(:);   
	spl(:, 2*cnt) = dp1b(:);  
	wgt(cnt) = h1b;  typ(cnt) = 1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h1c) > 1e-10)
	spl(:, 2*cnt-1) = p1c(:);   
	spl(:, 2*cnt) = dp1c(:);  
	wgt(cnt) = h1c;  typ(cnt) = 1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h2a) > 1e-10)
	spl(:, 2*cnt-1) = p2a(:);   
	spl(:, 2*cnt) = dp2a(:);  
	wgt(cnt) = h2a;  typ(cnt) = 2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h2b) > 1e-10)
	spl(:, 2*cnt-1) = p2b(:);   
	spl(:, 2*cnt) = dp2b(:);  
	wgt(cnt) = h2b;  typ(cnt) = 2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h2c) > 1e-10)
	spl(:, 2*cnt-1) = p2c(:);   
	spl(:, 2*cnt) = dp2c(:);  
	wgt(cnt) = h2c;  typ(cnt) = 2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(h3a) > 1e-10)
	spl(:, 2*cnt-1) = p3a(:);   
	spl(:, 2*cnt) = dp3a(:);  
	wgt(cnt) = h3a;  typ(cnt) = 3;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k1a) > 1e-10)
	spl(:, 2*cnt-1) = p1ka(:);   
	spl(:, 2*cnt) = dp1ka(:);  
	wgt(cnt) = k1a;  typ(cnt) = -1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k1b) > 1e-10)
	spl(:, 2*cnt-1) = p1kb(:);   
	spl(:, 2*cnt-1) = dp1kb(:);  
	wgt(cnt) = k1b;  typ(cnt) = -1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k1c) > 1e-10)
	spl(:, 2*cnt-1) = p1kc(:);   
	spl(:, 2*cnt) = dp1kc(:);  
	wgt(cnt) = k1c;  typ(cnt) = -1;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k2a) > 1e-10)
	spl(:, 2*cnt-1) = p2ka(:);   
	spl(:, 2*cnt) = dp2ka(:);  
	wgt(cnt) = k2a;  typ(cnt) = -2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k2b) > 1e-10)
	spl(:, 2*cnt-1) = p2kb(:);   
	spl(:, 2*cnt) = dp2kb(:);  
	wgt(cnt) = k2b;  typ(cnt) = -2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k2c) > 1e-10)
	spl(:, 2*cnt-1) = p2kc(:);   
	spl(:, 2*cnt) = dp2kc(:);  
	wgt(cnt) = k2c;  typ(cnt) = -2;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

if(abs(k3a) > 1e-10)
	spl(:, 2*cnt-1) = p3ka(:);   
	spl(:, 2*cnt) = dp3ka(:);  
	wgt(cnt) = k3a;  typ(cnt) = -3;  cut(cnt) = wavcut;  cnt = cnt + 1;
end

tempEntry.samples.nonlocal = spl;
tempEntry.cutoffs.nonlocal = cut;
tempEntry.NLweights = wgt;
tempEntry.NLtypes = typ;


% ----------------------------------------------------------------------
% check normalization condition
% ----------------------------------------------------------------------

InfoPrint(0, 'Check for normalization condition for element %3d\n', Znuc);
InfoPrint(0, 'Total number of pseudopotentials terms %3d\n', cnt-1);

% Local pseudopotential
pk = find(r>0 & r<rhocut); 
InfoPrint(0, 'int rho         = %15.10e\n', sum(4*pi*r(pk).^2.*(-rho(pk))*stp) );

% Nonlocal pseudopotential
pk = find(r>0 & r<wavcut);
if( abs(h0a) > 1e-10 )
	InfoPrint(0, 'Abs: int p0a^2  = %15.10e\n', sum(r(pk).^2.*p0a(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p0a^2  = %15.10e\n', ...
		sum(r(pk).^2.*p0a(pk).^2*stp) / (sum(r.^2.*p0a.^2*stp)/2) );
end
if( abs(h0b) > 1e-10 )
	InfoPrint(0, 'Abs: int p0b^2  = %15.10e\n', sum(r(pk).^2.*p0b(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p0b^2  = %15.10e\n', ...
		sum(r(pk).^2.*p0b(pk).^2*stp) / (sum(r.^2.*p0b.^2*stp)/2) );
end
if( abs(h0c) > 1e-10 )
	InfoPrint(0, 'Abs: int p0c^2  = %15.10e\n', sum(r(pk).^2.*p0c(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p0c^2  = %15.10e\n', ...
		sum(r(pk).^2.*p0c(pk).^2*stp) / (sum(r.^2.*p0c.^2*stp)/2) );
end
if( abs(h1a) > 1e-10 )
	InfoPrint(0, 'Abs: int p1a^2  = %15.10e\n', sum(r(pk).^2.*p1a(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1a^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1a(pk).^2*stp) / (sum(r.^2.*p1a.^2*stp)/2) );
end
if( abs(h1b) > 1e-10 )
	InfoPrint(0, 'Abs: int p1b^2  = %15.10e\n', sum(r(pk).^2.*p1b(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1b^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1b(pk).^2*stp) / (sum(r.^2.*p1b.^2*stp)/2) );
end
if( abs(h1c) > 1e-10 )
	InfoPrint(0, 'Abs: int p1c^2  = %15.10e\n', sum(r(pk).^2.*p1c(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1c^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1c(pk).^2*stp) / (sum(r.^2.*p1c.^2*stp)/2) );
end
if( abs(h2a) > 1e-10 )
	InfoPrint(0, 'Abs: int p2a^2  = %15.10e\n', sum(r(pk).^2.*p2a(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2a^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2a(pk).^2*stp) / (sum(r.^2.*p2a.^2*stp)/2) );
end
if( abs(h2b) > 1e-10 )
	InfoPrint(0, 'Abs: int p2b^2  = %15.10e\n', sum(r(pk).^2.*p2b(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2b^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2b(pk).^2*stp) / (sum(r.^2.*p2b.^2*stp)/2) );
end
if( abs(h2c) > 1e-10 )
	InfoPrint(0, 'Abs: int p2c^2  = %15.10e\n', sum(r(pk).^2.*p2c(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2c^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2c(pk).^2*stp) / (sum(r.^2.*p2c.^2*stp)/2) );
end
if( abs(h3a) > 1e-10 )
	InfoPrint(0, 'Abs: int p3a^2  = %15.10e\n', sum(r(pk).^2.*p3a(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p3a^2  = %15.10e\n', ...
		sum(r(pk).^2.*p3a(pk).^2*stp) / (sum(r.^2.*p3a.^2*stp)/2) );
end

if( abs(k1a) > 1e-10 )
	InfoPrint(0, 'Abs: int p1ka^2  = %15.10e\n', sum(r(pk).^2.*p1ka(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1ka^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1ka(pk).^2*stp) / (sum(r.^2.*p1ka.^2*stp)/2) );
end
if( abs(k1b) > 1e-10 )
	InfoPrint(0, 'Abs: int p1kb^2  = %15.10e\n', sum(r(pk).^2.*p1kb(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1kb^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1kb(pk).^2*stp) / (sum(r.^2.*p1kb.^2*stp)/2) );
end
if( abs(k1c) > 1e-10 )
	InfoPrint(0, 'Abs: int p1kc^2  = %15.10e\n', sum(r(pk).^2.*p1kc(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p1kc^2  = %15.10e\n', ...
		sum(r(pk).^2.*p1kc(pk).^2*stp) / (sum(r.^2.*p1kc.^2*stp)/2) );
end
if( abs(k2a) > 1e-10 )
	InfoPrint(0, 'Abs: int p2ka^2  = %15.10e\n', sum(r(pk).^2.*p2ka(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2ka^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2ka(pk).^2*stp) / (sum(r.^2.*p2ka.^2*stp)/2) );
end
if( abs(k2b) > 1e-10 )
	InfoPrint(0, 'Abs: int p2kb^2  = %15.10e\n', sum(r(pk).^2.*p2kb(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2kb^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2kb(pk).^2*stp) / (sum(r.^2.*p2kb.^2*stp)/2) );
end
if( abs(k2c) > 1e-10 )
	InfoPrint(0, 'Abs: int p2kc^2  = %15.10e\n', sum(r(pk).^2.*p2kc(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p2kc^2  = %15.10e\n', ...
		sum(r(pk).^2.*p2kc(pk).^2*stp) / (sum(r.^2.*p2kc.^2*stp)/2) );
end
if( abs(k3a) > 1e-10 )
	InfoPrint(0, 'Abs: int p3ka^2  = %15.10e\n', sum(r(pk).^2.*p3ka(pk).^2*stp) );
	InfoPrint(0, 'Rel: int p3ka^2  = %15.10e\n', ...
		sum(r(pk).^2.*p3ka(pk).^2*stp) / (sum(r.^2.*p3ka.^2*stp)/2) );
end

InfoPrint(0, '\n');


end