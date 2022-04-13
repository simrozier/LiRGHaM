(* ::Package:: *)

(* ::Text:: *)
(*Producing plots for the linear response of the MW halo to the LMC.*)


(* ::Section::Closed:: *)
(*Useful tools and general parameters*)


(* ::Input::Initialization:: *)
$SetHistoryLength=1;
(*Generic options for the plots.*)
fontname = "Book Antiqua";
fontsize = 13;
fontsizesmall = 8;
basestyleplot = {FontFamily -> fontname,FontSize->fontsize};
labelstyleplot = {FontFamily -> fontname, FontSize -> fontsize};
labelstyleplot2 = {FontFamily -> fontname, FontSize -> 1.5fontsize};
imagesize = 350;
(*Function defining the precision of some numbers*)
choprad[r_?NumericQ]:= r 10^10 // Round // #/10^10 & //N ; 
(*Function removing the point from integers. See https://mathematica.stackexchange.com/questions/7463/is-this-the-most-efficient-way-to-round-approximate-integers-to-integers-while-l*)
IntegerChop=With[{r=Round[#]},r+Chop[#-r]]&;
(**)
(*Function replacing a list of resonance numbers with a String where "{" \[Rule] "l", "," \[Rule] "u", "}" \[Rule] "r" *)
toLURString[resListf_]:=StringReplace[Module[{i=0},StringReplace[ToString[resListf],(("\""/;(++i;False))|", "/;(EvenQ[i])):>","]],{","->"u","{"->"l","}"->"r"}];
(**)
(*Building a function which takes a complex output from Julia and transforms it into a MMA usable table.*)
importDataFromJulia[string_]:=Block[{data, expoTabData},
data = Import[string, "Table", "TextDelimiters"->" "];
data = Map[StringReplace[ToString[#],{"im" -> " *I"}]&,data,{2}];
data = Map[StringSplit[#," "]&,data,{2}];
data = Map[Flatten,data,{1}];
expoTabData = Map[StringCases[#,{y___~~"e"~~x___ ->x}]&,data,{2} ];
expoTabData = Map[ToExpression,expoTabData,{3}];
expoTabData = Map[Cases[#,x___->If[x<-50,-50,x]]&,expoTabData,{2}];
Quiet@(data = MapIndexed[StringReplace[#, {y___~~"e"~~x___ ->y~~"e"~~ToString[expoTabData[[#2[[1]],#2[[2]]]][[1]]]}]&, data, {2}]);
data = Map[StringReplace[ToString[#],{"e" -> "*10^"}]&,data,{2}];
data = Map[Partition[#,4]&,data];
data = Map[#[[1]]<>#[[2]]<>#[[3]]<>#[[4]]&,data,{2}];
data = Map[ToExpression[#]&,data,{2}];
Transpose[data]];
(**)
(*Turn the green from "RedGreenSplit" into blue.*)
greenToBlue[color_]:=Block[{redBase,greenBase,fR,fG},
redBase = 0.3;
greenBase = 0.45;
fR[x_]:=(0.97-redBase)/(0.95-0.5)(x-0.5)+redBase;
fG[x_]:=(0.97-greenBase)/(0.95-0.5)(x-0.5)+greenBase;
If[color[[1]]<1.,RGBColor[{fR[#[[2]]],fG[#[[2]]],1.}]&@color,color]];
(*Define redBlueSplit*)
redBlueSplit = Function[x,If[x<0.5,Blend[{RGBColor[0.05, 0.05, 1],RGBColor[{1,1,1}]},2x],Blend[{RGBColor[1., 0.0054530410000000005`, 0.0054530410000000005`],RGBColor[{1,1,1}]},2(1-x)]]];
pinkCyanSplit = Function[x,If[x<0.5,Blend[{RGBColor[{Rational[123, 256], Rational[203, 256], Rational[255, 256]}],RGBColor[{1,1,1}]},2x],Blend[{RGBColor[{Rational[255, 256], Rational[7, 16], Rational[153, 256]}],RGBColor[{1,1,1}]},2(1-x)]]];


(* ::Section::Closed:: *)
(*Prerequisites*)


(* ::Subsection::Closed:: *)
(*Units*)


(* ::Text:: *)
(*To go from real to theoretical units, divide by the following numbers. To recover the true values from the theoretical ones, multiply by them.*)


(* ::Input::Initialization:: *)
(*Scale radius of the Hernquist sphere for the Milky Way halo.*)
aMW = 40.8; (*kpc*)
(*Milky Way halo mass*)
MMW = 1.57*10^12; (*Msun*)
(*Newton's constant.*)
G = 4.30*10^-6;(*(km/s)^2 kpc/Msun*)


(* ::Input::Initialization:: *)
(*Milky Way halo virial radius*)
RvirMW = 280; (*kpc*)


(* ::Input::Initialization:: *)
(*Velocity unit*)
vUnit = Sqrt[G MMW / aMW]; (*Convertion from km/s to G = MMW = aMW = 1.*)
(*Time unit*)
tUnit = aMW / vUnit; (*Convertion from s kpc / km to G = MMW = aMW = 1.*)
(*Density unit*)
\[Rho]Unit = MMW / aMW^3;


(* ::Input::Initialization:: *)
(*List of times*)
nTStepsIntegration = 621;
nTStepsOrbit = 20;
tmin = 0.;
(*Time step of the integration.*)
\[CapitalDelta]tIntegration = 0.01 / tUnit;
\[CapitalDelta]tOrbit = 10*\[CapitalDelta]tIntegration;
(*Max time of the final orbit.*)
tmaxIntegration = nTStepsIntegration * \[CapitalDelta]tIntegration;
tmaxOrbit = nTStepsOrbit *\[CapitalDelta]tOrbit;


(* ::Subsection::Closed:: *)
(*LMC's density*)


(* ::Input::Initialization:: *)
(*Normalised Hernquist density*)
\[Rho]HernquistNormalised[r_]=1/(2Pi r (1 + r)^3);


(* ::Input::Initialization:: *)
(*Off-centered LMC density.*)
MLMC = 1.8*10^11 ;
aLMC = 20;
\[Rho]LMC[r_,rcen_] =( MLMC / aLMC^3 / \[Rho]Unit)  \[Rho]HernquistNormalised[Norm[r - rcen] / (aLMC / aMW)];


(* ::Subsection::Closed:: *)
(*LMC' s orbit*)


(* ::Input:: *)
(*(*We use rp and vp from Salem+2015 and simply integrate back the orbit in the MW potential.*)*)


(* ::Input::Initialization:: *)
\[Psi]HERN[r_]=- 1 /(r + 1);
d\[Psi]HERNdr[r_] = D[\[Psi]HERN[r],r];
(*Angular frequency of circular orbits*)
  \[CapitalOmega]HERN[rad_]=Sqrt[1/rad D[\[Psi]HERN[rad],rad]];
(*Epicycle frequency*)
\[Kappa]HERN[rad_]=Sqrt[rad D[\[CapitalOmega]HERN[rad]^2,rad] + 4 \[CapitalOmega]HERN[rad]^2]//Simplify;


(* ::Input:: *)
(*Clear[lpfgintg]*)


(* ::Input::Initialization:: *)
(*Orbit integrator. Leap Frog, Herquist potential. Initial state given by rp and vp.*)
lpfgintg[tmin_,tmax_,\[CapitalDelta]t_,rp_,vp_]:=lpfgintg[tmin,tmax,\[CapitalDelta]t,rp,vp]=Block[{Eloc,Lloc,x0,y0,vx0,vy0,x,y,r,vx,vy,v,vxh,vyh,vh,ax,ay,i,xvtab},
Eloc=vp^2/2 + \[Psi]HERN[rp];
Lloc=rp vp;
x0=rp;
y0=0;
vx0=0;
vy0=Sqrt[2(Eloc-\[Psi]HERN[rp])];
(*Initialization*)
       x=x0;
y=y0;
r = Sqrt[x^2 + y^2];
vx=vx0;
vy=vy0;
v = Sqrt[vx^2 + vy^2];
ax=-(x/r) d\[Psi]HERNdr[r];
ay=-(y/r) d\[Psi]HERNdr[r];
xvtab=Table[0,{Floor[(tmax-tmin)/\[CapitalDelta]t]},{4}];
xvtab[[1]]={x0,y0,vx0,vy0};
For[i=1,i<Floor[(tmax-tmin)/\[CapitalDelta]t],i++,
(*Leap-frog integrator. ATTENTION: error in the B&T*)
vxh=vx+ax \[CapitalDelta]t/2;
vyh=vy+ay \[CapitalDelta]t/2;
x=x+vxh \[CapitalDelta]t;
y=y+vyh \[CapitalDelta]t;
r = Sqrt[x^2 + y^2];
ax=-(x/r) d\[Psi]HERNdr[r];
ay=-(y/r) d\[Psi]HERNdr[r];
vx=vxh+ax \[CapitalDelta]t/2;
vy = vyh + ay \[CapitalDelta]t/2;
xvtab[[i+1]]={x,y,vx,vy};];
xvtab];


(* ::Input::Initialization:: *)
(*Pericenter values of the LMC.*)
rpLMC = 48; (*kpc*)
vpLMC = 340; (*km/s*)


(* ::Subsection::Closed:: *)
(*Integrator for the Kepler problem*)


(* ::Input::Initialization:: *)
(*Orbit integrator. Leap Frog, Herquist potential. Initial state given by rp and vp.*)
lpfgintgAccKep[tmin_,tmax_,\[CapitalDelta]t_,rp_,vp_]:=(*lpfgintgAcc[tmin,tmax,\[CapitalDelta]t,rp,vp]=*)Block[{Eloc,Lloc,x10,x20,y10,y20,vx10,vx20,vy10,vy20,x1,x2,y1,y2, x, y,r,vx1, vx2, vy1,vy2,v,vx1h,vy1h, vx2h, vy2h,vh,ax1, ax2,ay1, ay2,i,xvtab},
Eloc=vp^2/2 -G MMW/rp;
Lloc=rp vp;
x10=0;
y10=0;
vx10=0;
vy10=0;
x20=rp;
y20=0;
vx20=0;
vy20=vp;
(*Initialization*)
x1=x10;
y1=y10;
       x2=x20;
y2=y20;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
vx1=vx10;
vx2=vx20;
vy1=vy10;
vy2=vy20;
v = Sqrt[(vx2-vx1)^2 + (vy2-vy1)^2];
ax1=x/r MLMC/MMW/r^2;
ay1=y/r MLMC/MMW/r^2;
ax2=-(x/r) /r^2;
ay2=-(y/r)  /r^2;
xvtab=Table[0,{Floor[(tmax-tmin)/\[CapitalDelta]t]},{8}];(*x1, y1, vx1, vy1, x2, y2, vx2, vy2*)
xvtab[[1]]={x10,y10,vx10,vy10, x20,y20,vx20,vy20};
For[i=1,i<Floor[(tmax-tmin)/\[CapitalDelta]t],i++,
(*Leap-frog integrator. ATTENTION: error in the B&T*)
vx1h=vx1+ax1 \[CapitalDelta]t/2;
vy1h=vy1+ay1 \[CapitalDelta]t/2;
vx2h=vx2+ax2 \[CapitalDelta]t/2;
vy2h=vy2+ay2 \[CapitalDelta]t/2;
x1=x1+vx1h \[CapitalDelta]t;
y1=y1+vy1h \[CapitalDelta]t;
x2=x2+vx2h \[CapitalDelta]t;
y2=y2+vy2h \[CapitalDelta]t;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
ax1=x/r MLMC/MMW/r^2;
ay1=y/r MLMC/MMW/r^2;
ax2=-(x/r) /r^2;
ay2=-(y/r)  /r^2;
vx1=vx1h+ax1 \[CapitalDelta]t/2;
vy1 = vy1h + ay1 \[CapitalDelta]t/2;
vx2=vx2h+ax2 \[CapitalDelta]t/2;
vy2= vy2h + ay2 \[CapitalDelta]t/2;
xvtab[[i+1]]={x1,y1,vx1,vy1, x2, y2, vx2, vy2};];
xvtab];


(* ::Subsection::Closed:: *)
(*LMC' s orbit with total acceleration*)


(* ::Input:: *)
(*(*We use rp and vp from Salem+2015 and simply integrate back the orbit in the MW potential.*)*)


(* ::Input::Initialization:: *)
\[Psi]HERN[r_]=- 1 /(r + 1);
\[Psi]HERNSat[r_]=\[Psi]HERN[r * aMW / aLMC] * MLMC / MMW * aMW / aLMC;
\[Rho]HERN[r_]= 1 /(2Pi)/r/(r + 1)^3;
\[Rho]HERNSat[r_]=\[Rho]HERN[r * aMW / aLMC] * MLMC / MMW * (aMW / aLMC)^3;
d\[Psi]HERNdr[r_] = D[\[Psi]HERN[r],r];
d\[Psi]HERNSatdr[r_]=d\[Psi]HERNdr[r * aMW / aLMC] * MLMC / MMW * (aMW / aLMC)^2;
(*Angular frequency of circular orbits*)
  \[CapitalOmega]HERN[rad_]=Sqrt[1/rad D[\[Psi]HERN[rad],rad]];
\[CapitalOmega]HERNSat[rad_]=\[CapitalOmega]HERN[rad * aMW / aLMC] * Sqrt[MLMC / MMW] * (aMW / aLMC)^(3/2);
(*Epicycle frequency*)
\[Kappa]HERN[rad_]=Sqrt[rad D[\[CapitalOmega]HERN[rad]^2,rad] + 4 \[CapitalOmega]HERN[rad]^2]//Simplify;
\[Kappa]HERNSat[rad_]=\[Kappa]HERN[rad * aMW / aLMC] * Sqrt[MLMC / MMW] * (aMW / aLMC)^(3/2);
(*Pericenter values of the LMC.*)
rpLMC = 48; (*kpc*)
vpLMC = 340; (*km/s*)


(* ::Input::Initialization:: *)
(*Force of two interpenetrating spheres.*)
testForce[rLMC_]:=Quiet[-2 Pi NIntegrate[d\[Psi]HERNdr[Sqrt[R^2 + z^2]] * z / Sqrt[R^2 + z^2] * R *\[Rho]HERNSat[Sqrt[R^2 + (z - rLMC)^2]],{R,0,Infinity},{z,-Infinity,Infinity}]];
forceIntp = Interpolation[Table[{x,testForce[x]},{x,0.1,20.0,0.1}]];


(* ::Input::Initialization:: *)
(*Orbit integrator. Leap Frog, Herquist potential. Initial state given by rp and vp.*)
lpfgintgAccBack[tmin_,tmax_,\[CapitalDelta]t_,rp_,vp_]:=(*lpfgintgAcc[tmin,tmax,\[CapitalDelta]t,rp,vp]=*)Block[{Eloc,Lloc,x10,x20,y10,y20,vx10,vx20,vy10,vy20,x1,x2,y1,y2, x, y,r,vx1, vx2, vy1,vy2,v,vx1h,vy1h, vx2h, vy2h,vh,ax1, ax2,ay1, ay2,i,xvtab},
x10=0;
y10=0;
vx10=0;
vy10=0;
x20=rp;
y20=0;
vx20=0;
vy20=vp;
(*Initialization*)
x1=x10;
y1=y10;
       x2=x20;
y2=y20;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
vx1=vx10;
vx2=vx20;
vy1=vy10;
vy2=vy20;
v = Sqrt[(vx2-vx1)^2 + (vy2-vy1)^2];
ax1=- x / r  forceIntp[r];
ay1=-(y/r)   forceIntp[r];
ax2=x/r   forceIntp[r]MMW/MLMC;
ay2=y/r   forceIntp[r]MMW/MLMC;
xvtab=Table[0,{Floor[(tmax-tmin)/\[CapitalDelta]t]},{8}];(*x1, y1, vx1, vy1, x2, y2, vx2, vy2*)
xvtab[[1]]={x10,y10,vx10,vy10, x20,y20,vx20,vy20};
For[i=1,i<Floor[(tmax-tmin)/\[CapitalDelta]t],i++,
(*Leap-frog integrator. ATTENTION: error in the B&T*)
vx1h=vx1-ax1 \[CapitalDelta]t/2;
vy1h=vy1-ay1 \[CapitalDelta]t/2;
vx2h=vx2-ax2 \[CapitalDelta]t/2;
vy2h=vy2-ay2 \[CapitalDelta]t/2;
x1=x1-vx1h \[CapitalDelta]t;
y1=y1-vy1h \[CapitalDelta]t;
x2=x2-vx2h \[CapitalDelta]t;
y2=y2-vy2h \[CapitalDelta]t;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
ax1=- x / r  forceIntp[r];
ay1=-(y/r)   forceIntp[r];
ax2=x/r   forceIntp[r]MMW/MLMC;
ay2=y/r   forceIntp[r]MMW/MLMC;
vx1=vx1h-ax1 \[CapitalDelta]t/2;
vy1 = vy1h - ay1 \[CapitalDelta]t/2;
vx2=vx2h-ax2 \[CapitalDelta]t/2;
vy2= vy2h - ay2 \[CapitalDelta]t/2;
xvtab[[i+1]]={x1,y1,vx1,vy1, x2, y2, vx2, vy2};];
xvtab];


(* ::Subsection::Closed:: *)
(*LMC' s orbit with cusp acceleration*)


(* ::Input:: *)
(*(*We use rp and vp from Salem+2015.*)*)


(* ::Input::Initialization:: *)
\[Psi]HERN[r_]=- 1 /(r + 1);
\[Psi]HERNSat[r_]=\[Psi]HERN[r * aMW / aLMC] * MLMC / MMW * aMW / aLMC;
\[Rho]HERN[r_]= 1 /(2Pi)/r/(r + 1)^3;
\[Rho]HERNSat[r_]=\[Rho]HERN[r * aMW / aLMC] * MLMC / MMW * (aMW / aLMC)^3;
d\[Psi]HERNdr[r_] = D[\[Psi]HERN[r],r];
d\[Psi]HERNSatdr[r_]=d\[Psi]HERNdr[r * aMW / aLMC] * MLMC / MMW * (aMW / aLMC)^2;
(*Angular frequency of circular orbits*)
  \[CapitalOmega]HERN[rad_]=Sqrt[1/rad D[\[Psi]HERN[rad],rad]];
\[CapitalOmega]HERNSat[rad_]=\[CapitalOmega]HERN[rad * aMW / aLMC] * Sqrt[MLMC / MMW] * (aMW / aLMC)^(3/2);
(*Epicycle frequency*)
\[Kappa]HERN[rad_]=Sqrt[rad D[\[CapitalOmega]HERN[rad]^2,rad] + 4 \[CapitalOmega]HERN[rad]^2]//Simplify;
\[Kappa]HERNSat[rad_]=\[Kappa]HERN[rad * aMW / aLMC] * Sqrt[MLMC / MMW] * (aMW / aLMC)^(3/2);
(*Pericenter values of the LMC.*)
rpLMC = 48; (*kpc*)
vpLMC = 340; (*km/s*)


(* ::Input::Initialization:: *)
(*Force applied by the LMC on the cusp only.*)
forceCusp[rLMC_]:=Quiet[-2 Pi NIntegrate[d\[Psi]HERNSatdr[Sqrt[R^2 +( z-rLMC)^2]] * (z-rLMC) / Sqrt[R^2 + (z-rLMC)^2] * R *\[Rho]HERN[Sqrt[R^2 + z^2]],{R,0,rLMC},{z,-Sqrt[rLMC^2-R^2],Sqrt[rLMC^2-R^2]}]];
forceCuspIntp = Interpolation[Table[{x,forceCusp[x]},{x,0.1,20.0,0.1}]];


(* ::Input::Initialization:: *)
(*Orbit integrator. Leap Frog, Herquist potential. Initial state given by rp and vp.*)
lpfgintgCuspAccBack[tmin_,tmax_,\[CapitalDelta]t_,rp_,vp_]:=(*lpfgintgAcc[tmin,tmax,\[CapitalDelta]t,rp,vp]=*)Block[{Eloc,Lloc,x10,x20,y10,y20,vx10,vx20,vy10,vy20,x1,x2,y1,y2, x, y,r,vx1, vx2, vy1,vy2,v,vx1h,vy1h, vx2h, vy2h,vh,ax1, ax2,ay1, ay2,i,xvtab},
x10=0;
y10=0;
vx10=0;
vy10=0;
x20=rp;
y20=0;
vx20=0;
vy20=vp;
(*Initialization*)
x1=x10;
y1=y10;
       x2=x20;
y2=y20;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
vx1=vx10;
vx2=vx20;
vy1=vy10;
vy2=vy20;
v = Sqrt[(vx2-vx1)^2 + (vy2-vy1)^2];
ax1= x / r  (*forceCuspIntp[r]*)d\[Psi]HERNSatdr[r];
ay1=y/r   (*forceCuspIntp[r]*)d\[Psi]HERNSatdr[r];
ax2=x/r   forceIntp[r]MMW/MLMC;
ay2=y/r   forceIntp[r]MMW/MLMC;
xvtab=Table[0,{Floor[(tmax-tmin)/\[CapitalDelta]t]},{8}];(*x1, y1, vx1, vy1, x2, y2, vx2, vy2*)
xvtab[[1]]={x10,y10,vx10,vy10, x20,y20,vx20,vy20};
For[i=1,i<Floor[(tmax-tmin)/\[CapitalDelta]t],i++,
(*Leap-frog integrator. ATTENTION: error in the B&T*)
vx1h=vx1-ax1 \[CapitalDelta]t/2;
vy1h=vy1-ay1 \[CapitalDelta]t/2;
vx2h=vx2-ax2 \[CapitalDelta]t/2;
vy2h=vy2-ay2 \[CapitalDelta]t/2;
x1=x1-vx1h \[CapitalDelta]t;
y1=y1-vy1h \[CapitalDelta]t;
x2=x2-vx2h \[CapitalDelta]t;
y2=y2-vy2h \[CapitalDelta]t;
x=x2-x1;
y=y2-y1;
r = Sqrt[x^2 + y^2];
ax1= x / r  (*forceCuspIntp[r]*)d\[Psi]HERNSatdr[r];
ay1=y/r   (*forceCuspIntp[r]*)d\[Psi]HERNSatdr[r];
ax2=x/r   forceIntp[r]MMW/MLMC;
ay2=y/r   forceIntp[r]MMW/MLMC;
vx1=vx1h-ax1 \[CapitalDelta]t/2;
vy1 = vy1h - ay1 \[CapitalDelta]t/2;
vx2=vx2h-ax2 \[CapitalDelta]t/2;
vy2= vy2h - ay2 \[CapitalDelta]t/2;
xvtab[[i+1]]={x1,y1,vx1,vy1, x2, y2, vx2, vy2};];
xvtab];


(* ::Subsection::Closed:: *)
(*Clutton - Brock basis*)


(* ::Input::Initialization:: *)
(*Definition of the bi-orthogonal basis from Clutton-Brock 1973. 
This basis is constructed from ultraspherical polynomials, but unlike that of Hernquist & Ostriker 1992, the first term corresponds to a cored density.
Some modifications are made to account for:
- a change in the definition of Gegenbauer polynomials between 1973 and today.
- a prefactor which normalises the basis.*)
(**)
lmax=6;
nmode=80;
(*Proportionality factor between potential and density*)
KCB[n_,l_]:=KCB[n,l]=N@(4n(n+2l+2)+(2l+1)(2l+3));
(*Normalisation term.*)
ICB[n_,l_]:=ICB[n,l]=Pi 2^(1-2l)(N@Gamma[n+2l]/N[n!])/N[(n+l)]/N@Gamma[l]^2;
(*Intermediate variable*)
\[Chi][x_]= ((x)^2 - 1)/((x)^2 + 1);
(*DEFINITION OF THE POTENTIAL*)
(*Analytic version of the radial functions. Different from Clutton-Brock 1973: the Gegenbauer polynomials have shifted indexes l\[Rule]l+1.
The prefactor (as compared to CB73) is required to guarantee the normalisation of the basis.*)
UbasisCBa[l_,n_,r_,Rbasisf_]:=-Sqrt[4Pi 2^(2 l +3)/KCB[n,l]/ICB[n,l+1]/Rbasisf] (r/Rbasisf)^l/((1 + (r/Rbasisf)^2)^(l + 1/2)) GegenbauerC[n,l+1,Evaluate@\[Chi][r,Rbasisf]];
(*(*Compiled version of the basis. It requires to define the variable usableUpoly. This formulation ensures that only the function we want is compiled. See https://mathematica.stackexchange.com/questions/180377/efficient-compilation-of-a-list-of-functions.*)
usableUCB =(Table[{lcomp==lt&&ncomp==nt,UbasisCBa[lt,nt,rcomp,Rbcomp]},{lt,0,lmax},{nt,0,nmode}]//Flatten[#,1]&//Flatten)//Join[#,{True,0.}]&//Apply[Which,#]&;
UbasisCB = Compile[{{lcomp,_Integer},{ncomp,_Integer},{rcomp,_Real},{Rbcomp,_Real}},
(*ATTENTION TO THE INDEX OFFSETS THAT COME FROM THE FACT THAT l STARTS AT 0.*)
usableUCB
,
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];*)
(*dPsi/dr*)
dUbasisCBDRa[l_,n_,r_,Rbasisf_]=D[UbasisCBa[l,n,r,Rbasisf],r];
(*usableDUCBDR =(Table[{lcomp==lt&&ncomp==nt,dUbasisCBDRa[lt,nt,rcomp,Rbcomp]},{lt,0,lmax},{nt,0,nmode}]//Flatten[#,1]&//Flatten)//Join[#,{True,0.}]&//Apply[Which,#]&;
dUbasisCBDR = Compile[{{lcomp,_Integer},{ncomp,_Integer},{rcomp,_Real},{Rbcomp,_Real}},
usableDUCBDR,
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];*)
(**)
(*DEFINITION OF THE DENSITY*)
(*Analytic version of the basis. Different from Clutton-Brock 1973: the Gegenbauer polynomials have shifted indexes l\[Rule]l+1.
The prefactor (as compared to CB73) is required to guarantee the normalisation of the basis.*)
DbasisCoeff[n_,l_]:=DbasisCoeff[n,l]=Sqrt[(2^(2 l +3)KCB[n,l])/(4 Pi ICB[n,l+1])] ;
DbasisCBa[l_,n_,r_,Rbasisf_]:= With[{x=N[r/Rbasisf]},1/Rbasisf^2.5  N@DbasisCoeff[n,l](x)^(l )/((1 + (x)^2)^(l + 5/2)) GegenbauerC[n,l+1,N@\[Chi][x]]];
(*(*Compiled version of the basis. It requires to define the variable usableUpoly. This formulation ensures that only the function we want is compiled. See https://mathematica.stackexchange.com/questions/180377/efficient-compilation-of-a-list-of-functions.*)
usableDCB =(Table[{lcomp==lt&&ncomp==nt,DbasisCBa[lt,nt,rcomp,Rbcomp]},{lt,0,lmax},{nt,0,nmode}]//Flatten[#,1]&//Flatten)//Join[#,{True,0.}]&//Apply[Which,#]&;
DbasisCB = Compile[{{lcomp,_Integer},{ncomp,_Integer},{rcomp,_Real},{Rbcomp,_Real}},
(*ATTENTION TO THE INDEX OFFSETS THAT COME FROM THE FACT THAT l STARTS AT 0.*)
usableDCB
,
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];*)


(* ::Input::Initialization:: *)
(*.*)
mmax=6;
(*Compiled version of the spherical harmonics.*)
usableY=(Table[{lcomp==lt&&mcomp==mt,SphericalHarmonicY[lt,mt,\[Theta]comp,\[Phi]comp]},{lt,0,lmax},{mt,-lt,lt,1}]//Flatten[#,1]&//Flatten)//Join[#,{True,0.}]&//Apply[Which,#]&;
SphericalHarmonicYC=Compile[{{lcomp,_Integer},{mcomp,_Integer},{\[Theta]comp,_Real},{\[Phi]comp,_Real}},
usableY
,
CompilationTarget->"C",CompilationOptions->{"ExpressionOptimization" -> True,"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True},RuntimeOptions->{"CatchMachineOverflow"->False , "CatchMachineUnderflow"->False,"CatchMachineIntegerOverflow"->False,"CompareWithTolerance"->False,"EvaluateSymbolically"->False}];


(* ::Input::Initialization:: *)
Ub[Ubasisf_] := 
  Ub[Ubasisf] = 
   Which[
   ToString[Ubasisf] == ToString[Hernquist], UbasisPOLY,
 ToString[Ubasisf] == ToString[Clutton], UbasisCB,
ToString[Ubasisf] == ToString[Bessel], UbasisISOa
    ];
Db[Ubasisf_] := 
  Db[Ubasisf] = 
   Which[
   ToString[Ubasisf] == ToString[Hernquist], DbasisPOLY,
 ToString[Ubasisf] == ToString[Clutton], DbasisCB,
 ToString[Ubasisf] == ToString[CluttonA], DbasisCBa,
ToString[Ubasisf] == ToString[Bessel], DbasisISOa
    ];
nStart[Ubasisf_] := 
  nStart[Ubasisf] = 
   Which[
   ToString[Ubasisf] == ToString[Hernquist], 0,
 ToString[Ubasisf] == ToString[Clutton], 0,
 ToString[Ubasisf] == ToString[CluttonA], 0,
ToString[Ubasisf] == ToString[Bessel], 1
    ];


(* ::Subsection::Closed:: *)
(*Bessel basis*)


(* ::Input::Initialization:: *)
(* Memorization of BesselJZero. *)
UbasisZero[l_,n_] = BesselJZero[l-0.5,n];
(*Amplitude of the basis elements.*)
AmplitudeUbasis[l_,n_,Rbasisf_]= -(Sqrt[8\[Pi] /Rbasisf])/(UbasisZero[l,n]*Abs[SphericalBesselJ[l,UbasisZero[l,n]]]);
(*Radial component of the potential elements*)
UbasisISOa[l_,n_,r_,Rbasisf_] = AmplitudeUbasis[l,n,Rbasisf]*SphericalBesselJ[l,UbasisZero[l,n]*r/Rbasisf];
(*Amplitude of the density basis element.*)
AmplitudeDbasis[l_,n_,Rbasisf_]=-AmplitudeUbasis[l,n,Rbasisf]/(4Pi )(UbasisZero[l,n]/Rbasisf)^2;
(*Radial component of the density elements.*)
DbasisISOa[l_,n_,r_,Rbasisf_]=AmplitudeDbasis[l,n,Rbasisf]SphericalBesselJ[l,UbasisZero[l,n]*r/Rbasisf];
dUbasisISOa[l_, n_, r_, Rbasisf_] = D[UbasisISOa[l, n, r, Rbasisf], r];


(* ::Subsection::Closed:: *)
(*Deprojection*)


(* ::Input::Initialization:: *)
(*Deprojection of the vector*)
deprojectionLMCResponse[Ubasisf_,Rbasisf_,vector_,r_, \[Theta]_, \[Phi]_,m_,l_,n_]:= Block[{vectloc,Dbloc,sumM0, sumOtherM},
Dbloc=Db[Ubasisf];
vectloc=vector;
sumM0 = Sum[Sum[vectloc[[l1 + 1,1,n1+1 - nStart[Ubasisf]]]N[Dbloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,0,\[Theta],\[Phi]],{l1,0,l}];
sumOtherM = Sum[Sum[vectloc[[l1 + 1, m1 + 1, n1+1 - nStart[Ubasisf]]]N[Dbloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,m1,\[Theta],\[Phi]],{l1,0,l},{m1,1, Min[l1,m]}];
sumM0 + 2 Re[sumOtherM]];
(*Deprojection of the vector with a single l*)
deprojectionLMCResponseSingleL[Ubasisf_,Rbasisf_,vector_,r_, \[Theta]_, \[Phi]_,m_,l_,n_]:= Block[{vectloc,Dbloc,sumM0, sumOtherM},
Dbloc=Db[Ubasisf];
vectloc=vector;
sumM0 = Sum[Sum[vectloc[[l1 + 1,1,n1+1 - nStart[Ubasisf]]]N[Dbloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,0,\[Theta],\[Phi]],{l1,{l}}];
sumOtherM = Sum[Sum[vectloc[[l1 + 1, m1 + 1, n1+1 - nStart[Ubasisf]]]N[Dbloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,m1,\[Theta],\[Phi]],{l1,{l}},{m1,1, Min[l1,m]}];
sumM0 + 2 Re[sumOtherM]];


(* ::Input::Initialization:: *)
(*Deprojection of the vector potential with a single l*)
deprojectionPotentialSingleL[Ubasisf_,Rbasisf_,vector_,r_, \[Theta]_, \[Phi]_,m_,l_,n_]:= Block[{vectloc,Ubloc,sumM0, sumOtherM},
Ubloc=Ub[Ubasisf];
vectloc=vector;
sumM0 = Sum[Sum[vectloc[[l1 + 1,1,n1+1 - nStart[Ubasisf]]]N[Ubloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,0,\[Theta],\[Phi]],{l1,{l}}];
sumOtherM = Sum[Sum[vectloc[[l1 + 1, m1 + 1, n1+1 - nStart[Ubasisf]]]N[Ubloc[l1,n1,r,Rbasisf]],{n1,nStart[Ubasisf],n}]SphericalHarmonicYC[l1,m1,\[Theta],\[Phi]],{l1,{l}},{m1,1, Min[l1,m]}];
sumM0 + 2 Re[sumOtherM]];


(* ::Input::Initialization:: *)
(*Function giving the full resonance list for a given value of l and n1bound.*)
fullResonances[l_,n1bound_]:= (Join@@Table[({#,lloc}&/@Range[-n1bound,n1bound]),{lloc,l,1,-2}])~Join~If[EvenQ[l],({#,0}&/@Range[1,n1bound]),{}];


(* ::Input::Initialization:: *)
harmonicFactors[l_, m_, phiMin_, phiMax_, nPoints_] :=harmonicFactors[l, m, phiMin, phiMax, nPoints] =Table[ SphericalHarmonicYC[l,m,Pi/2,\[Phi]],{\[Phi],phiMin,phiMax,(phiMax - phiMin)/nPoints}];


(* ::Section:: *)
(*figLMCFullResponse*)


(*Data directory*)
nbdir = "/Users/rozier/Documents/PostDoc/sphericalresponse/Time_Response/Julia_code_published/LiRGHaM/plots/";
newDataDirectory = nbdir<>"../source/MW_LMC_Response/Data/";


(* ::Subsection:: *)
(*Parameters*)


(* ::Input::Initialization:: *)
Rbasis = 11; (*Scale radius of the basis elements*) 
RbasisS = "11.0"; (*String version of the scale radius (keeping the ".0")*)
RbasisString[Rbasis]=RbasisS;  (*Wrapping the correspondence in a function*)
beta = 0.01;  (*Anisotropy parameter of the MW model.*)
betaS = "0.01"; (*String version of the anisotropy parameter*)
BetaString[beta]=betaS; (*Wrapping the correspondence in a function*)
nmax = 3; (*Maximum order of the radial basis functions.*)
uBasis = "Clutton"; (*Type of basis functions.*)
ellMax = 1; (*Maximum harmonic ell.*)


(* ::Subsection:: *)
(*Individual panels and bar*)


(*Files with values of the density.*)
responseDensityRadialFileName[l_,m_, resonanceList_,beta_,nmax_, Ubasisf_, Rbasisf_] := newDataDirectory<>"ResponseDensityRadial_To_LMC_Ell"<>ToString[l]<>"_m"<>ToString[m]<>"_resonances"<>toLURString[resonanceList]<>"_beta"<>BetaString[beta]<>"_nmax"<>ToString[nmax]<>"_Ubasis"<>Ubasisf<>"_Rbasis"<>Rbasisf<>".txt";
responseDensityRadial[l_,m_, resonanceList_, beta_, nmax_, Ubasisf_, Rbasisf_]:=responseDensityRadial[l, m, resonanceList, beta, nmax, Ubasisf, Rbasisf] =  Transpose@importDataFromJulia[responseDensityRadialFileName[l, m, resonanceList,beta, nmax,Ubasisf, RbasisString[Rbasisf]]];
responseDensityRadial[l_,m_, {},beta_, nmax_, Ubasisf_, Rbasisf_] = Table[0.,{nmax + 1},{21}];
(*Files with values of the density.*)
bareResponseDensityRadialFileName[l_,m_, resonanceList_,beta_,nmax_, Ubasisf_, Rbasisf_] := newDataDirectory<>"BareResponseDensityRadial_To_LMC_Ell"<>ToString[l]<>"_m"<>ToString[m]<>"_resonances"<>toLURString[resonanceList]<>"_beta"<>BetaString[beta]<>"_nmax"<>ToString[nmax]<>"_Ubasis"<>Ubasisf<>"_Rbasis"<>Rbasisf<>".txt";
bareResponseDensityRadial[l_,m_, resonanceList_, beta_,nmax_, Ubasisf_, Rbasisf_]:=bareResponseDensityRadial[l, m, resonanceList,beta,nmax, Ubasisf, Rbasisf] = If[OddQ[l-m],Table[0.,{nmax + 1},{21}], Transpose@importDataFromJulia[bareResponseDensityRadialFileName[l, m, resonanceList,beta,nmax,Ubasisf, RbasisString[Rbasisf]]]];
bareResponseDensityRadial[l_,m_, {},beta_, nmax_, Ubasisf_, Rbasisf_] = Table[0.,{nmax + 1},{21}];


(* ::Input::Initialization:: *)
(*Loading the response.*)
Do[responseDensityRadial[l,m,fullResonances[l, 2], beta, nmax, uBasis, Rbasis],{l,0,ellMax},{m,l,0,-2}];


(* ::Input::Initialization:: *)
(*Full trajectory of the LMC in (x,y) coordinates.*)
LMCFullTrajectory = {-#[[1]],#[[2]]}&/@lpfgintg[0,tmaxIntegration,\[CapitalDelta]tIntegration,rpLMC / aMW,vpLMC / vUnit][[408;;608;;1,{1,2}]];
LMCTrajectory = LMCFullTrajectory[[1;;201;;10]];


(* ::Input::Initialization:: *)
(*Plotting the response*)
(*Scales*)
rmax=300/aMW;
xmin = -rmax;
xmax = rmax;
ymin = -rmax;
ymax = rmax;
zmaxP = 0.6 ;
(*Ticks*)
xticks=Join[Table[{u,""},{u,-300/aMW,xmax,100/aMW}],Table[{u,Style[ToString[IntegerChop[u*aMW]],fontname,fontsize]},{u,-300/aMW,xmax,100/aMW}]];
yticks=Join[Table[{u,""},{u,-300/aMW,ymax,100/aMW}],Table[{u,Style[ToString[IntegerChop[u*aMW]],fontname,fontsize]},{u,-300/aMW,ymax,100/aMW}]];
zticks=Join[Table[{u,""},{u,-zmaxP,zmaxP,zmaxP/6}],Table[{u,Style[ToString[IntegerChop[u]],fontname,0.8fontsize]},{u,-zmaxP,zmaxP,zmaxP/3}]];
(*Axes labels*)
ax1 = Style["x (kpc)",fontname,fontsize,Black];
ax2 = Style["y (kpc)",fontname,fontsize,Black];
ax3 =Row[{Superscript[Style["\[Rho]",fontname,fontsize,Black],Style["s",fontname,fontsizesmall,Black]],Subscript[ Style["/\[Rho]",fontname,fontsize,Black],Style["MW", fontname, fontsizesmall,Black]]}];
(*Axes styles*)
axstyle = {fontname,fontsize,Black};
(*plotting the perturber and the response*)
Block[{r,\[Theta],\[Phi], z0, t,lmax,l,n,npoints,rMinGrid,rMaxGrid,phiMinGrid,phiMaxGrid,phiRot, tTable,nmaxloc,npointsR,npoints\[Phi],ltab, resonancesList},
nmaxloc = nmax;
(*Grid of points in the orbital plane of the LMC.*)
npointsR = 200;
npoints\[Phi]= 200;
rMinGrid =0.001;
rMaxGrid = 10.401;
phiMinGrid = -0.999 Pi;
phiMaxGrid = 0.999 Pi;
(*Rotation applied to the whole plot, to have the Milky Way disk roughly in the horizontal plane.*)
phiRot = Pi/5;
(*All values of ell.*)
ltab = Range[0,ellMax];
(*All resonances for each ell.*)
resonancesList ={fullResonances[#,2]&/@ltab};
pLMCResponse= Table[
Block[{density,tabdensity,MWcirc,LMCcirc,time,lvalue, ell,x,y,totalTable,harmonicFact,respdens,tensProd},
(*We sum all values of the response coming from different ell,m harmonics. Here, we account for the \[PlusMinus]m symmetry. The values of the "radial density" coming from the data files are multiplied by the "harmonic factors" to recover the full density of each harmonic.*)
totalTable = Table[Sum[
(harmonicFact = harmonicFactors[ltab[[i]],m, phiMinGrid,phiMaxGrid,npoints\[Phi]];
respdens =responseDensityRadial[ltab[[i]],m,resonancesList[[p,i]] ,beta,nmaxloc, uBasis, Rbasis][[All,t]];
tensProd = TensorProduct[respdens, harmonicFact];
If[m==0, tensProd, 2 Re[tensProd]]),{i,1,ellMax + 1},{m,ltab[[i]],0,-2}],{p,1,1}];
(*Dividing the overdensity by the background to get the relative overdensity. Applying a slight smoothing to polish the plot.*)
tabdensity =Map[MovingAverage[#,10]&, Table[r=rMinGrid + (ir - 1)(rMaxGrid - rMinGrid)/npointsR;
\[Phi]=Mod[phiMinGrid + (iphi - 1)(phiMaxGrid - phiMinGrid)/npoints\[Phi]+Pi-phiRot,2Pi]-Pi;
{x,y}=FromPolarCoordinates[{r,\[Phi]}];{x,y,Min[(Re@ totalTable[[1,ir,iphi]])/(\[Rho]HernquistNormalised[r] ),0.999zmaxP]},{iphi,1,npoints\[Phi] + 1},{ir,1,npointsR + 1}]]//Flatten[#,1]&;
density = ListDensityPlot[tabdensity,PlotRange->{{xmin,xmax},{ymin,ymax},{-zmaxP,zmaxP}},ColorFunction->(redBlueSplit[Rescale[#1,{-zmaxP,zmaxP}]]&),ColorFunctionScaling->False,(*MaxRecursion\[Rule]2,PlotPoints\[Rule]50,*)PerformanceGoal->"Quality",FrameLabel->{ax1,ax2},LabelStyle->labelstyleplot,Frame->True,FrameTicks->{xticks,yticks},PlotLegends->BarLegend[{(redBlueSplit[Rescale[#1,{-zmaxP,zmaxP}]]&),{-zmaxP,zmaxP}},LegendLabel->ax3,Ticks->zticks],FrameTicksStyle->Directive[fontname,fontsize,Black],ImageSize->Medium];
(*Adding a circle for the MW.*)
MWcirc = Graphics[{Thick,Circle[{0,0},1.0]}];
(*Adding a circle for the LMC.*)
LMCcirc = Graphics[{Dashed,Thick,Circle[{{Cos[phiRot],Sin[phiRot]},{-Sin[phiRot],Cos[phiRot]}} . LMCTrajectory[[t]],aLMC/aMW]}];
(*Displaying the time.*)
time = Graphics[Text[Row[{Subscript[Style["t",FontFamily->fontname, 1.2fontsize],Style[ToString[t-1],FontFamily->fontname, 1.2fontsizesmall]],Style[" = ",FontFamily->fontname, 1.2fontsize],Style[ToString[IntegerChop@N[(t-1) /10]],FontFamily->fontname,1.2 fontsize],Style[" Gyr",FontFamily->fontname,1.2fontsize]}],{- 0.5 rmax,0.8 rmax},BaseStyle->Directive[FontFamily->fontname]],Background->Transparent];
Show[density,MWcirc,LMCcirc(*,time*)]],
{t,1,21,1}]];


(* ::Input::Initialization:: *)
(*Exporting the plots.*)
Table[Export[nbdir<>"figResponseMW"<>ToString[ (i-1)]<>".pdf",Rasterize[pLMCResponse[[i,1]],RasterSize->2000],"PDF"],{i,1,21}];
Export[nbdir<>"barResponseMW.pdf",pLMCResponse[[1,2,1]],"PDF"]


(* ::Input::Initialization:: *)
(*Exporting a movie.*)
Quiet@Export["FiducialResponse.gif",pLMCResponse[[All,1]],"GIF","DisplayDurations"->0.5] ;
