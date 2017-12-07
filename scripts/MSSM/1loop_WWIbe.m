(* ::Package:: *)

Quit[]


path = "/Users/jamesmckay/Documents/Programs/Mass_builder/";
kappa=1/(16\[Pi]^2);
sw=Sin[\[Theta]];
cw=Cos[\[Theta]];

(* Check that counter-term couplings agree with Ibe et al. (2013) *)
g=e/sw;
gammaZ = -e^2(-5/3);
mB0[m1_,m2_]= \[Delta];

mB1[m1_,m2_]=-\[Delta]/2;

mB21[m1_,m2_] =\[Delta]/3;

mB22[m1_,m2_]=-p^2 (mB1[m1,m2]+mB21[m1,m2])-(p^2/4)*mB0[m1,m2]-(1/4)*(m1^2-m2^2)*(mB0[m1,m2]+2*mB1[m1,m2]);

PI[m1_,m2_]=-p^2(mB1[m1,m2]+mB21[m1,m2]);

(* below we have in order from left to right, Leptons (e,mu,tau) , up type quarks and down type quarks (neutrinos have 0 charge so are not here) *)
(*SEIbeFermionLepton =  e^2*(1/(2*Pi^2))*( (-1)^2*(3*PI[mf,mf]) + (2/3)^2*(2*PI[mf,mf]+PI[mt,mt]) + (1/3)^2*(3*PI[mf,mf]));

SEIbeFermionLepton =  e^2*(1/(2*Pi^2))*( (-1)^2*(3*PI[mf,mf]) + 3(2/3)^2*(3*PI[mf,mf]) + 3(1/3)^2*(3*PI[mf,mf]));*)


SEIbeFermionLepton =  e^2*(1/(2*Pi^2))*( (-1)^2*(3*PI[mf,mf]) + 3(2/3)^2*(2*PI[mf,mf]+PI[mt,mt]) + 3(1/3)^2*(3*PI[mf,mf]));


SEIbeBoson = -  e^2*(3/(4*Pi^2))*( mB22[mw,mw]+ p^2/18  ) -e^2*p^2*(1/(4*Pi^2))*( mB0[mw,mw]);
SEIbeWino = e^2*(1/(2*Pi^2))*( PI[MChi,MChi] );


SEIbeTotal = FullSimplify[SEIbeWino+SEIbeBoson+SEIbeFermionLepton];
deltaIbe=FullSimplify[Coefficient[SEIbeTotal,\[Delta],1]];
deltaZgamma = - (e^2/(16*Pi^2))*(32/9 *Ng - 5/3)p^2;


eq1=FullSimplify[deltaZgamma+deltaIbe]
Solve[{eq1==0},{Ng}]





deltaZgamma/.Ng->3
