#include <lal/SphericalHarmonics.h>
#include <math.h>

static COMPLEX16 Complex(REAL8 real, REAL8 imag){
	COMPLEX16 z = real + imag * I;
	return z;
 }
/* All the hGO_l_m fucntions contain the 3.5PN quasi-circular spinning corrections as given in Henry et al, arXiv:2209.00374v2 */
//H22

static COMPLEX16 hGO_2_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT, UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 kappa1 = 1.0; /*for black holes kappa and lambda is 1*/
    REAL8 kappa2 = 1.0;
    REAL8 lambda1 = 1.0;
    REAL8 lambda2 = 1.0;
    REAL8 delta = sqrt(1-4*Nu);
    
    if(vpnorder == 0){
        return(mass/r + pow(PhiDOT,2)*pow(r,2) + Complex(0,2)*PhiDOT*r*rDOT - 
        pow(rDOT,2));
    }

    else if(vpnorder == 2){
        return ((21*pow(mass,2)*(-10 + Nu) - 27*(-1 
        + 3*Nu)*pow(r,2)*(PhiDOT*r - Complex(0,1)*rDOT)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,3) + mass*r*((11 + 156*Nu)*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,10)*(5 
        + 27*Nu)*PhiDOT*r*rDOT - 3*(15 + 32*Nu)*pow(rDOT,2)))/(42.*pow(r,2)) );
    }

    else if(vpnorder == 3){
        return ((pow(mass,2)*(Complex(0,-1)*rDOT*
        ((3 + 3*delta - 8*Nu)*S1z + (3 - 3*delta - 8*Nu)*S2z) + 
       PhiDOT*r*((-3 - 3*delta + 5*Nu)*S1z + (-3 + 3*delta + 5*Nu)*S2z)))/(3.*pow(r,2)) /* (<--This is the general orbit term) */
      /* (This is the quasi-circular limit of the general orbit term-->) */ - ((-4*((1 + delta - Nu)*S1z + S2z - (delta + Nu)*S2z)*pow(x,2.5))/3.) 
       + ((-4*(S1z + delta*S1z + S2z - delta*S2z - Nu*(S1z + S2z))*pow(x,2.5))/3.)); /* (<--This is Quentins quasi-circular term) */
       
    }

    else if(vpnorder == 4){
        return ((6*pow(mass,3)*(3028 + 1267*Nu + 158*pow(Nu,2)) 
        + 9*(83 - 589*Nu + 1111*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + pow(mass,2)*r*((-11891 - 36575*Nu + 13133*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,8)*(-773 - 3767*Nu + 2852*pow(Nu,2))*PhiDOT*r*rDOT - 6*(-619 
        + 2789*Nu + 934*pow(Nu,2))*pow(rDOT,2)) - 3*mass*pow(r,2)*(2*(-835 - 19*Nu
        + 2995*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) + Complex(0,6)*(-433 - 721*Nu 
        + 1703*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT + 6*(-33 + 1014*Nu 
        + 232*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) + Complex(0,4)*(-863 
        + 1462*Nu + 2954*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        - 3*(-557 + 664*Nu + 1712*pow(Nu,2))*pow(rDOT,4)))/(1512.*pow(r,3))
         + (3*pow(mass,3)*(S1z*(4*Nu*S2z + (1 + delta - 2*Nu)*S1z*kappa1) - 
       (-1 + delta + 2*Nu)*pow(S2z,2)*kappa2))/(4.*pow(r,3)) 
       - ((kappa1*(1 + delta - 2*Nu)*pow(S1z,2) + S2z*(4*Nu*S1z 
       - kappa2*(-1 + delta + 2*Nu)*S2z))*pow(x,3)) + ((kappa1*(1 + delta - 2*Nu)*pow(S1z,2) 
       + S2z*(4*Nu*S1z - kappa2*(-1 + delta + 2*Nu)*S2z))*pow(x,3)));
    }

    else if(vpnorder == 5){
        return ((pow(mass,2)*Nu*(2*mass*(Complex(0,-702)*PhiDOT*r + rDOT) 
        + 3*r*(Complex(0,-316)*pow(PhiDOT,3)*pow(r,3) 
        - 847*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,184)*PhiDOT*r*pow(rDOT,2) 
        - 122*pow(rDOT,3))))/(105.*pow(r,3)) + ((2*(56*delta*Nu*(-S1z + S2z) 
        + 101*Nu*(S1z + S2z) + 132*pow(Nu,2)*(S1z + S2z) - 80*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/63.) );
    }

    else if(vpnorder == 6){
        return ((4*pow(mass,4)*(-8203424 + 2180250*pow(Nu,2) 
        + 592600*pow(Nu,3) + 15*Nu*(-5503804 + 142065*pow(M_PI,2))) 
        - 2700*(-507 + 6101*Nu - 25050*pow(Nu,2) 
        + 34525*pow(Nu,3))*pow(r,4)*pow(PhiDOT*r - Complex(0,1)*rDOT,3)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,5) + pow(mass,3)*r*(pow(PhiDOT,2)*(337510808 
        - 198882000*pow(Nu,2) + 56294600*pow(Nu,3) + Nu*(183074880 
        - 6392925*pow(M_PI,2)))*pow(r,2) + Complex(0,110)*PhiDOT*(-5498800 
        - 785120*pow(Nu,2) + 909200*pow(Nu,3) + 3*Nu*(-1849216 
        + 38745*pow(M_PI,2)))*r*rDOT + 2*(51172744 - 94929000*pow(Nu,2) 
        - 5092400*pow(Nu,3) + 45*Nu*(2794864 + 142065*pow(M_PI,2)))*pow(rDOT,2)) 
        - 20*pow(mass,2)*pow(r,2)*((-986439 + 1873255*Nu - 9961400*pow(Nu,2) 
        + 6704345*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) + Complex(0,4)*(-273687 
        - 978610*Nu - 4599055*pow(Nu,2) + 2783005*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + (-181719 + 19395325*Nu + 8237980*pow(Nu,2) 
        + 2612735*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,8)*(-234312 + 1541140*Nu + 1230325*pow(Nu,2) 
        + 1828625*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) - 3*(-370268 + 1085140*Nu 
        + 2004715*pow(Nu,2) + 1810425*pow(Nu,3))*pow(rDOT,4)) 
        + 300*mass*pow(r,3)*(4*(12203 - 36427*Nu - 27334*pow(Nu,2) 
        + 149187*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) + Complex(0,2)*(44093 - 68279*Nu 
        - 295346*pow(Nu,2) + 541693*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT 
        + 2*(27432 - 202474*Nu + 247505*pow(Nu,2) 
        + 394771*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,2)*(97069 - 383990*Nu - 8741*pow(Nu,2) 
        + 1264800*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        + (-42811 + 53992*Nu + 309136*pow(Nu,2) - 470840*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,2)*(51699 - 252256*Nu + 131150*pow(Nu,2) 
        + 681160*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) - 3*(16743 - 75104*Nu + 26920*pow(Nu,2) 
        + 207200*pow(Nu,3))*pow(rDOT,6)))/(3.3264e6*pow(r,4)) + (((4*(1 + delta)*(-7 + 9*kappa1) 
        - 7*(9 + 17*delta)*Nu - 9*(15 + 7*delta)*kappa1*Nu + 12*(7 - 17*kappa1)*pow(Nu,2))*pow(S1z,2) + 
       2*S1z*(Complex(0,-42)*(1 + delta - 2*Nu) - 84*(1 + delta - Nu)*M_PI + Nu*(-271 + 288*Nu)*S2z) + 
       S2z*(12*(7 - 17*kappa2)*pow(Nu,2)*S2z + 4*(-1 + delta)*(Complex(0,21) + 42*M_PI + 7*S2z - 9*kappa2*S2z) + 
          Nu*(168*(Complex(0,1) + M_PI) + 7*delta*(17 + 9*kappa2)*S2z - 9*(7 + 15*kappa2)*S2z)))*pow(x,4))/63. );
    }
    else if(vpnorder==7){

        return(((3318*pow(Nu,3)*(S1z + S2z) + Nu*(-504*((7 + delta)*kappa1 - 3*(3 + delta)*lambda1)*pow(S1z,3) - 
          1008*pow(S1z,2)*(3*kappa1*M_PI - 3*(1 + delta)*S2z + 2*(1 + delta)*kappa1*S2z) + 
          S1z*(17387 + 20761*delta + 1008*S2z*(6*M_PI + (-1 + delta)*(-3 + 2*kappa2)*S2z)) + 
          S2z*(17387 - 20761*delta + 504*S2z*(-6*kappa2*M_PI + (-7 + delta)*kappa2*S2z - 3*(-3 + delta)*lambda2*S2z))) + 
          2*(2809*(1 + delta)*S1z + 756*(1 + delta)*kappa1*M_PI*pow(S1z,2) + 756*(1 + delta)*(kappa1 - lambda1)*pow(S1z,3) - 
          (-1 + delta)*S2z*(2809 + 756*S2z*(-(lambda2*S2z) + kappa2*(M_PI + S2z)))) - 
          2*pow(Nu,2)*(708*delta*(-S1z + S2z) + (S1z + S2z)*(4427 + 1008*(kappa1*pow(S1z,2) + S2z*(-2*S1z + kappa2*S2z)))))*pow(x,4.5))/756.);
    }

    else{
        return 0;
    } 
    
}
// hQC_l_m() functions contains only the hereditary terms at particular PN order.

static COMPLEX16 hQC_2_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    double EulerGamma = 0.5772156649015329;
   // keeping only the hereditary terms
    // if(vpnorder == 0){
    //   return(2*x);
    // }

    /* else if(vpnorder == 2){
      return((-5.095238095238095 + (55*Nu)/21.)*pow(x,2));
    } */

    /* else */ if(vpnorder == 3){
        return(4*M_PI*pow(x,2.5));
    }

    /* else if(vpnorder == 4){
       return((-2.874338624338624 - (1069*Nu)/108. + (2047*pow(Nu,2))/756.)*
       pow(x,3));
    } */

    else if(vpnorder == 5){
        //return((Complex(0,-48)*Nu - (214*M_PI)/21. + (68*Nu*M_PI)/21.)*pow(x,3.5));
        return((2*(-107 + 34*Nu)*M_PI*pow(x,3.5))/21.);
    }

    else if(vpnorder == 6){
        return((pow(x,4)*(-27392*EulerGamma + 
         M_PI*(Complex(0,13696) + 35*(64 + 41*Nu)*M_PI) - 13696*log(16*x)))/
         1680.);
    }

    else if(vpnorder == 7){
        return((Complex(0,176.9506172839506)*Nu - 
         Complex(0,8.605291005291006)*pow(Nu,2) - (2173*M_PI)/378. - 
         (2495*Nu*M_PI)/189. + (80*pow(Nu,2)*M_PI)/27.)*pow(x,4.5));
    }

    else{
        return 0;
    }

}

static COMPLEX16 hl_2_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_2_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_2_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_2_m_2(Nu,vpnorder,x))* cpolar(1,-2*Phi);
    }
}

static COMPLEX16 hl_2_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_2_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,2) * conj((hGO_2_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) + hQC_2_m_2(Nu,vpnorder,x))) * cpolar(1,2*Phi); 
    }
}

//H21

 static COMPLEX16 hGO_2_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;
    

    if(vpnorder == 1){
        return (Complex(0,0.6666666666666666)*delta*mass*PhiDOT );
    }

    else if(vpnorder == 2){
        return ((Complex(0,-0.5)*pow(mass,2)*((1 + delta)*S1z + (-1 + delta)*S2z))/pow(r,2)
         - (Complex(0,-0.5)*(S1z - S2z + delta*(S1z + S2z))*pow(x,2))
         + (Complex(0,-0.5)*(S1z - S2z + delta*(S1z + S2z))*pow(x,2)));
    }

    else if(vpnorder == 3){
        return ((Complex(0,0.023809523809523808)*delta*mass*PhiDOT*(4*mass*(-9 + 11*Nu) 
        + r*((19 - 24*Nu)*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,2)*(83 + 2*Nu)*PhiDOT*r*rDOT + 2*(-33 + 10*Nu)*pow(rDOT,2))))/r );
    }

    else if(vpnorder == 4){
        return ((Complex(0,0.011904761904761904)*pow(mass,2)*
     (2*mass*((77 + 17*Nu + 11*delta*(7 + Nu))*S1z + 
          (-77 - 17*Nu + 11*delta*(7 + Nu))*S2z) + 
       r*(Complex(0,-2)*PhiDOT*r*rDOT*
           (147*(1 + delta)*S1z + (-83 + 13*delta)*Nu*S1z + 
             147*(-1 + delta)*S2z + (83 + 13*delta)*Nu*S2z) + 
          pow(rDOT,2)*((105*(1 + delta) - 4*(13 + 15*delta)*Nu)*S1z + 
             (-105 + 15*delta*(7 - 4*Nu) + 52*Nu)*S2z) + 
          4*pow(PhiDOT,2)*pow(r,2)*
           ((-21 - 21*delta + 66*Nu + 4*delta*Nu)*S1z + 
             (21 - 21*delta - 66*Nu + 4*delta*Nu)*S2z))))/pow(r,3)
              - (Complex(0,0.023809523809523808)*((-7 + 205*Nu + delta*(-7 + 33*Nu))*S1z + 
             (7 - 205*Nu + delta*(-7 + 33*Nu))*S2z)*pow(x,3)) + (Complex(0,0.023809523809523808)*pow(x,3)*( 
             (-7 + 205*Nu + delta*(-7 + 33*Nu))*S1z + (7 - 205*Nu + delta*(-7 + 33*Nu))*S2z )));
    }

    else if(vpnorder == 5){
        return ((Complex(0,0.0013227513227513227)*delta*mass*PhiDOT*(10*pow(mass,2)*(31 
        - 205*Nu + 111*pow(Nu,2)) - 2*mass*r*((-197 + 5*Nu 
        + 660*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) + Complex(0,1)*(-3167 - 5278*Nu 
        + 201*pow(Nu,2))*PhiDOT*r*rDOT + 8*(202 + 587*Nu - 177*pow(Nu,2))*pow(rDOT,2)) 
        + 3*pow(r,2)*((152 - 692*Nu + 333*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,2)*(308 - 1607*Nu + 111*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        - 3*(75 - 560*Nu + 68*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,2)*(-265 + 526*Nu + 18*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        + (-241 + 550*Nu - 264*pow(Nu,2))*pow(rDOT,4))))/pow(r,2)
        + (Complex(0,-0.08333333333333333)*pow(x,3.5)*(2*(-((1 + delta)*(5 + kappa1)) 
        + 2*(6 + delta + (4 + 3*delta)*kappa1)*Nu)*pow(S1z,2) + S1z*(Complex(0,-3) - Complex(0,3)*delta + 6*(1 + delta)*M_PI 
        - 16*delta*Nu*S2z - Complex(0,3)*(1 + delta)*log(16)) + S2z*(6*(-1 + delta)*M_PI - 2*(-1 + delta)*(5 + kappa2)*S2z 
        + 4*(-6 + delta - 4*kappa2 + 3*delta*kappa2)*Nu*S2z - Complex(0,3)*(-1 + delta)*(1 + log(16))))));
    }

    else if(vpnorder == 6){
        return ((delta*pow(mass,2)*Nu*PhiDOT*(mass*(195*PhiDOT*r 
        - Complex(0,946)*rDOT) + 9*r*(270*pow(PhiDOT,3)*pow(r,3) 
        - Complex(0,483)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - 580*PhiDOT*r*pow(rDOT,2) + Complex(0,42)*pow(rDOT,3))))/(315.*pow(r,2)) 
        +(Complex(0,-0.0006613756613756613)*(-4*pow(Nu,2)*(-1365*(S1z - S2z)
         + 179*delta*(S1z + S2z)) + 6*(208*(1 + delta)*S1z + 63*(1 + delta)*kappa1*pow(S1z,3) + (-1 + delta)*S2z*(208 + 63*kappa2*pow(S2z,2))) - 
      Nu*(378*(3 + delta)*kappa1*pow(S1z,3) + 378*(1 + delta)*(-2 + kappa1)*pow(S1z,2)*S2z + S1z*(8351 + 7027*delta + 378*(-1 + delta)*(-2 + kappa2)*pow(S2z,2)) + 
        S2z*(-8351 + 7027*delta + 378*(-3 + delta)*kappa2*pow(S2z,2))))*pow(x,4)));
    }

    else if(vpnorder == 7){

        return(Complex(0,0.0003968253968253968)*pow(x,4.5)*(10*(6*(1 + delta)*(63 + 26*kappa1) - (84*(26 + 17*delta) + (1069 + 757*delta)*kappa1)*Nu + 
        (2844 + 312*delta + (1409 + 36*delta)*kappa1)*pow(Nu,2))*pow(S1z,2) + 
        S2z*(30*(14 - 431*Nu + delta*(-14 + 87*Nu))*M_PI + 10*(6*(-1 + delta)*(63 + 26*kappa2) + (2184 - 1428*delta + 1069*kappa2 - 757*delta*kappa2)*Nu + 
           (-2844 + 312*delta - 1409*kappa2 + 36*delta*kappa2)*pow(Nu,2))*S2z + 
        Complex(0,3)*(66 - 66*delta - 2957*Nu + 4405*delta*Nu - 20*(14 - 14*delta - 431*Nu + 87*delta*Nu)*log(2))) + 
        3*S1z*(Complex(0,-1)*(66 + 66*delta - 2957*Nu - 280*log(2)) + 5*
         (-28*(1 + delta)*M_PI + 2*(431 + 87*delta)*Nu*M_PI + delta*Nu*(Complex(0,881) + 16*(7 + 23*Nu)*S2z) - Complex(0,4)*(431*Nu + delta*(-14 + 87*Nu))*log(2)))));
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hQC_2_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 1){
        return(Complex(0,0.6666666666666666)*delta*pow(x,1.5));
    } */

    /* else if(vpnorder == 3){
        return((Complex(0,-0.40476190476190477)*delta + 
     Complex(0,0.47619047619047616)*delta*Nu)*pow(x,2.5));
    } */

    /* else */ if(vpnorder == 4){
        return((2*delta*pow(x,3)*(Complex(0,1)*M_PI + log(4)))/3.);
    }

    /* else if(vpnorder == 5){
        return((Complex(0,-0.2275132275132275)*delta - 
     Complex(0,2.693121693121693)*delta*Nu + 
     Complex(0,0.3134920634920635)*delta*pow(Nu,2))*pow(x,3.5));
    } */

    else if(vpnorder == 6){
        return(M_PI*(Complex(0,-0.40476190476190477)*delta*pow(x,4) + 
        Complex(0,0.14285714285714285)*delta*Nu*pow(x,4)) + 
        ((-17*delta*pow(x,4))/21. + (2*delta*Nu*pow(x,4))/7.)*log(2));
    }

    else{

        return 0;
    }

}


static COMPLEX16 hl_2_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_2_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R) * (hGO_2_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_2_m_1(Nu,vpnorder,x)) * cpolar(1,-1*Phi);
    }

}

static COMPLEX16 hl_2_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_2_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R) * pow(-1,2) * conj((hGO_2_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_2_m_1(Nu,vpnorder,x))) * cpolar(1,1*Phi);
    }
    
}

//H33
static COMPLEX16 hGO_3_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder, REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;


    if(vpnorder == 1){
        return ((sqrt(0.11904761904761904)*delta*(2*r*pow(Complex(0,1)*PhiDOT*r - rDOT,3) 
        + mass*(Complex(0,-7)*PhiDOT*r + 4*rDOT)))/(2.*r) 
        );
    } 

    else if(vpnorder == 3){
        return ((sqrt(0.11904761904761904)*delta*(6*(-5 
        + 19*Nu)*pow(r,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,4)*(Complex(0,1)*PhiDOT*r 
        + rDOT) + 2*pow(mass,2)*(Complex(0,-3)*(-101 + 43*Nu)*PhiDOT*r 
        + (-109 + 86*Nu)*rDOT) + 3*mass*r*(Complex(0,-12)*(1 + 4*Nu)*pow(PhiDOT,3)*pow(r,3) 
        + 6*(14 + 31*Nu)*pow(PhiDOT,2)*pow(r,2)*rDOT + Complex(0,3)*(33 
        + 62*Nu)*PhiDOT*r*pow(rDOT,2) 
        - 4*(8 + 17*Nu)*pow(rDOT,3))))/(36.*pow(r,2)));
    }    

    else if(vpnorder == 4){
        return ((Complex(0,-0.125)*sqrt(0.11904761904761904)*pow(mass,2)*
     (4*mass*(-1 + 5*Nu)*((1 + delta)*S1z + (-1 + delta)*S2z) + 
       r*(2*pow(rDOT,2)*(6*(1 + delta)*S1z - 5*(5 + 3*delta)*Nu*S1z + 
             (-6 + delta*(6 - 15*Nu) + 25*Nu)*S2z) + 
          pow(PhiDOT,2)*pow(r,2)*
           (-24*(1 + delta)*S1z + (119 + 33*delta)*Nu*S1z + 
             (24 - 119*Nu + 3*delta*(-8 + 11*Nu))*S2z) + 
          Complex(0,2)*PhiDOT*r*rDOT*
           (-18*(1 + delta)*S1z + (77 + 39*delta)*Nu*S1z + 
             (18 - 77*Nu + 3*delta*(-6 + 13*Nu))*S2z))))/pow(r,3) 
             - (Complex(0,-0.375)*sqrt(1.0714285714285714)*
            (-4*S1z + 19*Nu*(S1z - S2z) + 4*S2z - 4*delta*(S1z + S2z) + 5*delta*Nu*(S1z + S2z))*pow(x,3))
             + (Complex(0,-0.375)*sqrt(0.04285714285714286)*pow(x,3)*( 5*(-4 + 19*Nu + delta*(-4 + 5*Nu))*S1z + 5*(4 - 19*Nu + delta*(-4 + 5*Nu))*S2z )));
    } 

    else if(vpnorder == 5){
        return ((delta*(30*(183 - 1579*Nu + 3387*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(Complex(0,-1)*PhiDOT*r + rDOT,5) 
        + 10*pow(mass,3)*(Complex(0,-1)*(26473 - 27451*Nu + 9921*pow(Nu,2))*PhiDOT*r 
        + 4*(623 - 732*Nu + 1913*pow(Nu,2))*rDOT) 
        + 2*pow(mass,2)*r*(Complex(0,-11)*(-5353 - 13493*Nu 
        + 4671*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) + (-75243 
        - 142713*Nu + 192821*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,220)*(-256 + 781*Nu + 840*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) 
        - 10*(-756 + 8238*Nu + 7357*pow(Nu,2))*pow(rDOT,3)) 
        + 3*mass*pow(r,2)*(Complex(0,2)*(-7633 + 9137*Nu 
        + 28911*pow(Nu,2))*pow(PhiDOT,5)*pow(r,5) - 4*(-8149 
        + 1576*Nu + 43533*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        - Complex(0,2)*(-9297 - 19517*Nu 
        + 64839*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        - 32*(-1288 + 3667*Nu + 4056*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - Complex(0,5)*(-9851 + 17954*Nu + 40968*pow(Nu,2))*PhiDOT*r*pow(rDOT,4) 
        + 20*(-771 + 1126*Nu + 3616*pow(Nu,2))*pow(rDOT,5))))/(1584.*sqrt(210)*pow(r,3))
        +(Complex(0,1.125)*sqrt(1.0714285714285714)*(kappa1*(-1 - delta + 2*(2 + delta)*Nu)*pow(S1z,2) 
        + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 + delta)*Nu)*S2z))*pow(x,3.5)));
    }

    else if(vpnorder == 6){
        return (-(delta*pow(mass,2)*Nu*(668*pow(mass,2) 
        + 2*mass*r*(4081*pow(PhiDOT,2)*pow(r,2) + Complex(0,297)*PhiDOT*r*rDOT 
        - 452*pow(rDOT,2)) + 5*pow(r,2)*(1329*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,2926)*pow(PhiDOT,3)*pow(r,3)*rDOT 
        - 384*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,408)*PhiDOT*r*pow(rDOT,3) + 200*pow(rDOT,4))))/(36.*sqrt(210)*pow(r,4))
        +(Complex(0,-0.125)*sqrt(0.04285714285714286)*(-(Nu*((279 + delta)*S1z + (-279 + delta)*S2z)) 
        + 10*(S1z - S2z + delta*(S1z + S2z)) + pow(Nu,2)*(407*(S1z - S2z) + 241*delta*(S1z + S2z)))*pow(x,4)));
    }

    else if(vpnorder == 7){

        return((Complex(0,-0.020833333333333332)*pow(x,4.5)*(-270*(-6*Nu*(-2 + delta*(-2 + Nu) + 6*Nu) 
        + kappa1*(4 - Nu*(13 + 8*Nu) + delta*(4 + Nu*(-5 + 12*Nu))))*pow(S1z,2) + 
       S1z*(810*(-4 + 19*Nu + delta*(-4 + 5*Nu))*M_PI + 5*delta*Nu*(Complex(0,-541) 
       + 216*(-10 + 9*Nu)*S2z) + Complex(0,61560)*Nu*atanh(1/5) +  Complex(0,1)*(2349*delta 
       - 14207*Nu + 81*(29 - 80*log(1.5))) + Complex(0,1620)*delta*(-4 + 5*Nu)*log(1.5)) + 
       S2z*(810*(4 - 19*Nu + delta*(-4 + 5*Nu))*M_PI + Complex(0,1)*(2349*delta + 14207*Nu 
       + 81*(-29 + 80*log(1.5))) - 5*(54*kappa2*(-4 + Nu*(13 + 8*Nu) + delta*(4 + Nu*(-5 + 12*Nu)))*S2z + 
        Nu*(Complex(0,541)*delta - 324*(2 + delta*(-2 + Nu) - 6*Nu)*S2z - Complex(0,648)*(-19 + 5*delta)*atanh(1/5)) 
        - Complex(0,324)*delta*log(1.5802469135802468)) - Complex(0,486)*delta*log(1024))))/sqrt(210));
    }

    else{
        return 0;
    }
    
}


static COMPLEX16 hQC_3_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    double EulerGamma = 0.5772156649015329;

    /* if(vpnorder == 1){
        return(Complex(0,-1.5)*sqrt(1.0714285714285714)*delta*pow(x,1.5));
    } */

    /* else if(vpnorder == 3){
        return((Complex(0,3)*sqrt(4.285714285714286)*delta - 
     Complex(0,3)*sqrt(1.0714285714285714)*delta*Nu)*pow(x,2.5));
    } */

    /* else */ if(vpnorder == 4){
        return((9*sqrt(1.0714285714285714)*delta*pow(x,3)*
     (Complex(0,-1)*M_PI + log(2.25)))/2.);
    }

    /* else if(vpnorder == 5){
        return((Complex(0,-8.386363636363637)*sqrt(0.04285714285714286)*delta + 
     Complex(0,83.54545454545455)*sqrt(0.04285714285714286)*delta*Nu - 
     Complex(0,20.15909090909091)*sqrt(0.04285714285714286)*delta*
      pow(Nu,2))*pow(x,3.5));
    } */

    else if(vpnorder == 6){
        return(M_PI*(Complex(0,9)*sqrt(4.285714285714286)*delta*pow(x,4) - 
        Complex(0,6.75)*sqrt(1.0714285714285714)*delta*Nu*pow(x,4)) + 
        (-18*sqrt(4.285714285714286)*delta*pow(x,4) + 
        (27*sqrt(1.0714285714285714)*delta*Nu*pow(x,4))/2.)*log(1.5));
    }
    else if(vpnorder == 7){
        return(Complex(0,1.1149564720993292e-6)*sqrt(0.04285714285714286)*delta*
        pow(x,4.5)*(-465315528 + 74954880*EulerGamma + 13827800*Nu + 
        124985672*pow(Nu,2) - 19373424*pow(Nu,3) + 
        Complex(0,47279232)*M_PI - 10090080*pow(M_PI,2) - 
        4309305*Nu*pow(M_PI,2) - 94558464*log(1.5) - 
        Complex(0,121080960)*M_PI*log(1.5) + 37477440*log(16*x) + 
        121080960*log(1.5)*log(1.5)));
    }

    else{
        return 0;
    }
}


static COMPLEX16 hl_3_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_3_m_3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R) * (hGO_3_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_3_m_3(Nu,vpnorder,x)) * cpolar(1,-3*Phi);
    }

}

static COMPLEX16 hl_3_m_min3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_3_m_min3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,3) * conj(hGO_3_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) + hQC_3_m_3(Nu,vpnorder,x)) * cpolar(1,3*Phi);
    }
    
}

//H32
static COMPLEX16 hGO_3_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;


    if(vpnorder == 2){
        return (-(sqrt(0.7142857142857143)*mass*(-1 + 3*Nu)*PhiDOT*(4*PhiDOT*r 
        + Complex(0,1)*rDOT))/6.);
    }

    else if(vpnorder == 3){
        return ((sqrt(0.7142857142857143)*pow(mass,2)*Nu*
         (4*PhiDOT*r + Complex(0,1)*rDOT)*(S1z + S2z))/(3.*pow(r,2)) 
         - ((4*sqrt(0.7142857142857143)*Nu*(S1z + S2z)*pow(x,2.5))/3.)
         +((4*sqrt(0.7142857142857143)*Nu*(S1z + S2z)*pow(x,2.5))/3.));
    }


    else if(vpnorder == 4){
        return (-(mass*PhiDOT*(2*mass*((167 - 925*Nu 
        + 1615*pow(Nu,2))*PhiDOT*r + Complex(0,5)*(-82 
        + 239*Nu + 55*pow(Nu,2))*rDOT) - 3*r*(2*(-13 - 25*Nu 
        + 355*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) - Complex(0,60)*(-8 
        + 25*Nu + pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + 12*(-23 + 70*Nu + 10*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,5)*(-13 + 38*Nu 
        + 10*pow(Nu,2))*pow(rDOT,3))))/(108.*sqrt(35)*r));
    }

    else if(vpnorder == 5){
        return ((pow(mass,2)*Nu*PhiDOT*(Complex(0,7)*mass 
        + r*(Complex(0,49)*pow(PhiDOT,2)*pow(r,2) 
        + 90*PhiDOT*r*rDOT - Complex(0,6)*pow(rDOT,2))))/(4.*sqrt(35)*pow(r,2))
        +((sqrt(0.7142857142857143)*(43*delta*Nu*(S1z - S2z) - 3*Nu*(S1z + S2z) 
        - 26*pow(Nu,2)*(S1z + S2z) - 8*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/9.));
    }

    else if(vpnorder == 6){
        return (-(mass*PhiDOT*(4*pow(mass,2)*(2*(5377 + 6438*Nu 
        - 79866*pow(Nu,2) + 37348*pow(Nu,3))*PhiDOT*r - Complex(0,5)*(-4115 
        + 18399*Nu - 20276*pow(Nu,2) + 7*pow(Nu,3))*rDOT) - 4*mass*r*((4599 
        - 15737*Nu + 36259*pow(Nu,2) + 108563*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3) 
        - Complex(0,1)*(-34053 + 59698*Nu + 192949*pow(Nu,2) 
        + 16193*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*rDOT + (-59058 
        + 77983*Nu + 322468*pow(Nu,2) - 4264*pow(Nu,3))*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,5)*(-3387 + 8518*Nu + 8968*pow(Nu,2) + 884*pow(Nu,3))*pow(rDOT,3)) 
        + 3*pow(r,2)*(4*(-710 + 3892*Nu - 10655*pow(Nu,2) 
        + 24000*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5) + Complex(0,11)*(-1484 
        + 11693*Nu - 25006*pow(Nu,2) + 428*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + 4*(4161 - 25618*Nu + 29489*pow(Nu,2) 
        + 22078*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        + Complex(0,44)*(-151 + 1067*Nu - 2419*pow(Nu,2) 
        + 57*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) + 4*(2041 - 11680*Nu 
        + 19334*pow(Nu,2) + 3368*pow(Nu,3))*PhiDOT*r*pow(rDOT,4) 
        + Complex(0,5)*(477 - 2624*Nu 
        + 3862*pow(Nu,2) + 1160*pow(Nu,3))*pow(rDOT,5))))/(4752.*sqrt(35)*pow(r,2))
        +((2*sqrt(0.7142857142857143)*((2*Nu*(-5 - 5*delta + 4*Nu) + 3*kappa1*(1 + delta 
        - 2*(2 + delta)*Nu + 6*pow(Nu,2)))*pow(S1z,2) + S2z*(2*(4 + 9*kappa2)*pow(Nu,2)*S2z 
        - 3*(-1 + delta)*(Complex(0,2) + kappa2*S2z) + 2*Nu*(Complex(0,-15) + 6*M_PI 
        + 5*(-1 + delta)*S2z + 3*(-2 + delta)*kappa2*S2z)) + 2*S1z*(Complex(0,3) 
        + Complex(0,3)*delta + Nu*(Complex(0,-15) + 6*M_PI + 2*S2z - 10*Nu*S2z)))*pow(x,4))/9.));
    }

    else if(vpnorder == 7){

        return(((-26902*pow(Nu,3)*(S1z + S2z) - 4664*(S1z + delta*S1z + S2z - delta*S2z) + 
        Nu*(3960*(1 + delta)*kappa1*pow(S1z,3) + 3960*(1 + delta)*kappa1*pow(S1z,2)*S2z + 
        S1z*(28921 - 18889*delta - 3960*(-1 + delta)*kappa2*pow(S2z,2)) + 
        S2z*(28921 + 18889*delta - 3960*(-1 + delta)*kappa2*pow(S2z,2))) - 
        2*pow(Nu,2)*(1351*delta*(S1z - S2z) + 6*(S1z + S2z)*(6773 + 660*(kappa1*pow(S1z,2) 
       + S2z*(-2*S1z + kappa2*S2z)))))*pow(x,4.5))/(1188.*sqrt(35)));
    }

    else{
        return 0;
    }
    
}


static COMPLEX16 hQC_3_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    /* if(vpnorder == 2){
        return(((2*sqrt(0.7142857142857143))/3. - 2*sqrt(0.7142857142857143)*Nu)*
        pow(x,2));
    }

    else if(vpnorder == 4){
        return((-193/(27.*sqrt(35)) + (145*sqrt(0.7142857142857143)*Nu)/27. - 
        (73*sqrt(0.7142857142857143)*pow(Nu,2))/27.)*pow(x,3));
    } */

    /* else */ if(vpnorder == 5){
        return((4*sqrt(0.7142857142857143)*(1 - 3*Nu)*M_PI*pow(x,3.5))/3.);
    }

    /* else if(vpnorder == 6){
        return((-1451/(1188.*sqrt(35)) - (17387*Nu)/(1188.*sqrt(35)) + 
        (5557*pow(Nu,2))/(66.*sqrt(35)) - 
        (763*sqrt(1.4)*pow(Nu,3))/396.)*pow(x,4));
    } */

    else{
        return 0;
    }

}


static COMPLEX16 hl_3_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_3_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_3_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_3_m_2(Nu,vpnorder,x)) * cpolar(1,-2*Phi);
    }
}

static COMPLEX16 hl_3_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_3_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,3) * conj(hGO_3_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_3_m_2(Nu,vpnorder,x)) * cpolar(1,2*Phi); 
    }
}

//H31
static COMPLEX16 hGO_3_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;
    

    if(vpnorder == 1){
        return (delta*(mass*(Complex(0,7)*PhiDOT*r - 12*rDOT) - Complex(0,6)*r*(PhiDOT*r 
        - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,2)))/(6.*sqrt(14)*r);
    }

    else if(vpnorder == 3){
        return (delta*(Complex(0,6)*(-5 + 19*Nu)*pow(r,2)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,3) 
        + 2*pow(mass,2)*(Complex(0,1)*(-101 + 43*Nu)*PhiDOT*r 
        + (109 - 86*Nu)*rDOT) + 3*mass*r*(Complex(0,-4)*(-9 
        + 14*Nu)*pow(PhiDOT,3)*pow(r,3) + 6*(2 + 9*Nu)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - Complex(0,1)*(33 + 62*Nu)*PhiDOT*r*pow(rDOT,2) 
        + 4*(8 + 17*Nu)*pow(rDOT,3))))/(36.*sqrt(14)*pow(r,2));
    }

    else if(vpnorder == 4){
        return ((Complex(0,0.041666666666666664)*pow(mass,2)*
       (4*mass*(-1 + 5*Nu)*((1 + delta)*S1z + (-1 + delta)*S2z) - 
       r*(2*pow(rDOT,2)*(-6*(1 + delta)*S1z + 
             5*(5 + 3*delta)*Nu*S1z + 6*S2z - 6*delta*S2z + 
             5*(-5 + 3*delta)*Nu*S2z) + 
          pow(PhiDOT,2)*pow(r,2)*
           ((24 + 24*delta - 87*Nu + 31*delta*Nu)*S1z + 
             (-24 + 24*delta + 87*Nu + 31*delta*Nu)*S2z) + 
          Complex(0,2)*PhiDOT*r*rDOT*
           ((6 + 6*delta - 31*Nu + 35*delta*Nu)*S1z + 
             (-6 + 6*delta + 31*Nu + 35*delta*Nu)*S2z))))/
         (sqrt(14)*pow(r,3)) - ((Complex(0,0.041666666666666664)*(-4*S1z + 11*Nu*(S1z - S2z) 
         + 4*S2z - 4*delta*(S1z + S2z) + 13*delta*Nu*(S1z + S2z))*pow(x,3))/sqrt(14)) 
         +((Complex(0,0.041666666666666664)*(-4*S1z 
          + 11*Nu*(S1z - S2z) + 4*S2z - 4*delta*(S1z + S2z) 
          + 13*delta*Nu*(S1z + S2z))*pow(x,3))/sqrt(14)));
    }

    else if(vpnorder == 5){
        return ((delta*(Complex(0,-18)*(183 - 1579*Nu 
        + 3387*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,3)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,4) + 2*pow(mass,3)*(Complex(0,1)*(26473 
        - 27451*Nu + 9921*pow(Nu,2))*PhiDOT*r - 12*(623 - 732*Nu 
        + 1913*pow(Nu,2))*rDOT) + 2*pow(mass,2)*r*(Complex(0,-1)*(-8641 
        - 59189*Nu + 31959*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) + (-32635 
        - 29345*Nu + 29541*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - Complex(0,44)*(-256 + 781*Nu + 840*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) 
        + 6*(-756 + 8238*Nu + 7357*pow(Nu,2))*pow(rDOT,3)) 
        + 3*mass*pow(r,2)*(Complex(0,2)*(-2479 - 4505*Nu 
        + 16785*pow(Nu,2))*pow(PhiDOT,5)*pow(r,5) 
        + 4*(817 + 1220*Nu - 7449*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + Complex(0,6)*(-1679 + 1469*Nu 
        + 12233*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        - 32*(-460 + 421*Nu + 2514*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + Complex(0,1)*(-9851 + 17954*Nu + 40968*pow(Nu,2))*PhiDOT*r*pow(rDOT,4) 
        - 12*(-771 + 1126*Nu 
        + 3616*pow(Nu,2))*pow(rDOT,5))))/(4752.*sqrt(14)*pow(r,3))+((Complex(0,-0.041666666666666664)*(kappa1*(5 
        - 4*Nu + delta*(5 + 6*Nu))*pow(S1z,2) + S2z*(-12*delta*Nu*S1z + kappa2*(-5 + 5*delta + 4*Nu + 6*delta*Nu)*S2z))*
        pow(x,3.5))/sqrt(14)));
    }

    else if(vpnorder == 6){
        return ((delta*pow(mass,2)*Nu*(668*pow(mass,2) 
        - 2*mass*r*(727*pow(PhiDOT,2)*pow(r,2) 
        - Complex(0,99)*PhiDOT*r*rDOT + 452*pow(rDOT,2)) 
        + pow(r,2)*(-499*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,1534)*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 3072*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,680)*PhiDOT*r*pow(rDOT,3) 
        + 1000*pow(rDOT,4))))/(180.*sqrt(14)*pow(r,4))+((Complex(0,0.004629629629629629)*(70*(S1z - S2z 
        + delta*(S1z + S2z)) - pow(Nu,2)*(931*(S1z - S2z) + 45*delta*(S1z + S2z)) - 
         Nu*(59*(-S1z + S2z) + 99*delta*(S1z + S2z)))*pow(x,4))/sqrt(14)));
    }

    else if(vpnorder == 7){

        return((Complex(0,0.001388888888888889)*pow(x,4.5)*(10*(32*(1 + delta)*kappa1 
        - (36*(1 + delta) + (65 + delta)*kappa1)*Nu + 2*(150 + 33*delta + (76 + 22*delta)*kappa1)*pow(Nu,2))*pow(S1z,2) + 
        S2z*(Complex(0,-315)*delta*Nu + 30*(4 - 11*Nu + delta*(-4 + 13*Nu))*M_PI + 60*Nu*(6 - 50*Nu + delta*(-6 + 11*Nu))*S2z + 
          10*kappa2*(-32 + (65 - 152*Nu)*Nu + delta*(32 + Nu*(-1 + 44*Nu)))*S2z + Complex(0,3)*(-117 + 117*delta + 199*Nu - 80*log(2)) - 
          Complex(0,60)*(-11*Nu + delta*(-4 + 13*Nu))*log(2)) + S1z*(30*(-4 + 11*Nu + delta*(-4 + 13*Nu))*M_PI + 5*delta*Nu*(Complex(0,-63) 
          + 8*(-34 + 11*Nu)*S2z) - Complex(0,60)*(11*Nu + delta*(-4 + 13*Nu))*log(2) + Complex(0,3)*(117 + 117*delta - 199*Nu + 80*log(2)))))/sqrt(14));
    }

    else{
        return 0;
    }
    
}


static COMPLEX16 hQC_3_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
     REAL8 delta = sqrt(1-4*Nu);
     double EulerGamma = 0.5772156649015329;

     /* if(vpnorder == 1){
         return((Complex(0,0.16666666666666666)*delta*pow(x,1.5))/sqrt(14));
     } */

     /* else if(vpnorder == 3){
         return((Complex(0,-0.2222222222222222)*sqrt(0.2857142857142857)*delta - 
         (Complex(0,0.1111111111111111)*delta*Nu)/sqrt(14))*pow(x,2.5));
     } */

     /* else */ if(vpnorder == 4){
         return((Complex(0,0.16666666666666666)*delta*M_PI*pow(x,3))/sqrt(14) + 
         (delta*pow(x,3)*log(2))/(3.*sqrt(14)));
     }

     /* else if(vpnorder == 5){
         return(((Complex(0,0.5109427609427609)*delta)/sqrt(14) - 
         Complex(0,0.11447811447811448)*sqrt(0.2857142857142857)*delta*Nu - 
         (Complex(0,0.2079124579124579)*delta*pow(Nu,2))/sqrt(14))*
         pow(x,3.5));
     } */

     else if(vpnorder == 6){
         return((Complex(0,-0.027777777777777776)*delta*(16 + 7*Nu)*pow(x,4)*
         (M_PI - Complex(0,2)*log(2)))/sqrt(14));
     }

     else if(vpnorder == 7){
         return((Complex(0,-2.752978943455134e-9)*delta*pow(x,4.5)*
        (-430135880 + 74954880*EulerGamma + 681626456*Nu - 
        641035640*pow(Nu,2) + 68698000*pow(Nu,3) + 
        Complex(0,47279232)*M_PI - 10090080*pow(M_PI,2) - 
        38783745*Nu*pow(M_PI,2) + 244468224*log(2) + 
        Complex(0,121080960)*M_PI*log(2) + 121080960*pow(log(2),2) + 
        37477440*log(x)))/sqrt(14));
     }
     
     else{
         return 0;
     }

}


static COMPLEX16 hl_3_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_3_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_3_m_1(Nu,vpnorder,x)) * cpolar(1,-1*Phi);
    }

}

static COMPLEX16 hl_3_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,3) * conj(hGO_3_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_3_m_1(Nu,vpnorder,x)) * cpolar(1,1*Phi);
    }
    
}

//H44
static COMPLEX16 hGO_4_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;

    
    if(vpnorder == 2){
        return ((sqrt(0.7142857142857143)*(-1 + 3*Nu)*(7*pow(mass,2) 
        + 6*pow(r,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + 3*mass*r*(17*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,18)*PhiDOT*r*rDOT - 6*pow(rDOT,2))))/(36.*pow(r,2)));
    }
    
    else if(vpnorder == 4){
        return ((40*pow(mass,3)*(314 - 987*Nu + 195*pow(Nu,2)) 
        - 60*(23 - 159*Nu + 291*pow(Nu,2))*pow(r,3)*(PhiDOT*r 
        - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,5) 
        + pow(mass,2)*r*((53143 - 199660*Nu + 127500*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,24)*(967 - 4615*Nu + 5935*pow(Nu,2))*PhiDOT*r*rDOT 
        - 10*(290 - 2033*Nu + 4365*pow(Nu,2))*pow(rDOT,2)) 
        - 3*mass*pow(r,2)*((613 - 920*Nu + 6420*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,8)*(-976 + 1745*Nu + 3150*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 2*(-6141 + 8980*Nu + 31500*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,4)*(-1853 + 1730*Nu + 13230*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        - 20*(-83 + 30*Nu + 762*pow(Nu,2))*pow(rDOT,4)))/(1584.*sqrt(35)*pow(r,3)));
    }

    else if(vpnorder == 5){
        return ((pow(mass,2)*Nu*(6*mass*(Complex(0,-43)*PhiDOT*r + 9*rDOT) 
        + r*(Complex(0,-734)*pow(PhiDOT,3)*pow(r,3) + 129*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,156)*PhiDOT*r*pow(rDOT,2) - 26*pow(rDOT,3))))/(24.*sqrt(35)*pow(r,3))
        +((-32*(39*delta*Nu*(S1z - S2z) + 41*Nu*(S1z + S2z) - 42*pow(Nu,2)*(S1z + S2z) 
        - 10*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/(27.*sqrt(35))));
    }

    else if(vpnorder == 6){
        return ((10*pow(mass,4)*(-4477296 + 12734393*Nu - 6895*pow(Nu,2) 
        + 1043805*pow(Nu,3)) + 3150*(-367 + 4337*Nu - 17462*pow(Nu,2) 
        + 23577*pow(Nu,3))*pow(r,4)*pow(PhiDOT*r - Complex(0,1)*rDOT,2)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,6) + 2*pow(mass,3)*r*((-36967579 + 245501977*Nu 
        - 459916170*pow(Nu,2) + 150200680*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,4)*(7571073 - 10780154*Nu - 56898800*pow(Nu,2) 
        + 43665510*pow(Nu,3))*PhiDOT*r*rDOT - 10*(1283609 - 5800627*Nu + 3725295*pow(Nu,2) 
        + 4771935*pow(Nu,3))*pow(rDOT,2)) - pow(mass,2)*pow(r,2)*((-28258134 
        + 3245207*Nu + 144051250*pow(Nu,2) + 136991820*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,24)*(2371982 - 7733376*Nu - 7948185*pow(Nu,2) 
        + 9074870*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT + 7*(6557973 
        - 50558069*Nu + 59901380*pow(Nu,2) 
        + 104752320*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,168)*(52044 - 1084807*Nu + 1849450*pow(Nu,2) 
        + 4171730*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) - 35*(1083 - 1246819*Nu 
        + 2524240*pow(Nu,2) + 5995845*pow(Nu,3))*pow(rDOT,4)) 
        - 105*mass*pow(r,3)*((116396 - 551405*Nu + 560658*pow(Nu,2) 
        + 293036*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) + Complex(0,2)*(158192 
        - 670661*Nu + 177718*pow(Nu,2) + 2163976*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT 
        + (-393665 + 1322392*Nu + 1589680*pow(Nu,2) 
        - 8622660*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        - Complex(0,8)*(-23048 + 209397*Nu - 487057*pow(Nu,2) 
        + 260396*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) - (630647 
        - 3391000*Nu + 2501958*pow(Nu,2) 
        + 7664096*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        - Complex(0,2)*(218975 - 1037408*Nu + 148970*pow(Nu,2) 
        + 3699480*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) 
        + 10*(10233 - 44864*Nu - 13050*pow(Nu,2) 
        + 203280*pow(Nu,3))*pow(rDOT,6)))/(1.44144e6*sqrt(35)*pow(r,4))
        +((16*sqrt(0.7142857142857143)*(-1 + 3*Nu)*(kappa1*(1 + delta - 2*Nu)*pow(S1z,2) 
        + S2z*(4*Nu*S1z - kappa2*(-1 + delta + 2*Nu)*S2z))*pow(x,4))/9.));
    }

    else if(vpnorder == 7){

        return((8*(6630*pow(Nu,3)*(S1z + S2z) - 786*(S1z + delta*S1z + S2z - delta*S2z) 
        + 21*Nu*(273*delta*(S1z - S2z) + 173*(S1z + S2z)) - 2*pow(Nu,2)*(1557*delta*(S1z - S2z)
         + 2960*(S1z + S2z)))*pow(x,4.5))/(297.*sqrt(35)));
    }

    else{
        return 0;
    }
    
}


static COMPLEX16 hQC_4_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    /* if(vpnorder == 2){
        return(((-16*sqrt(0.7142857142857143))/9. + 
        (16*sqrt(0.7142857142857143)*Nu)/3.)*pow(x,2));
    } */

    /* else if(vpnorder == 4){
        return((4744/(99.*sqrt(35)) - (10184*sqrt(0.7142857142857143)*Nu)/297. + 
        (200*sqrt(35)*pow(Nu,2))/99.)*pow(x,3));
    } */

    /* else */ if(vpnorder == 5){
        return((64*sqrt(0.7142857142857143)*(-1 + 3*Nu)*pow(x,3.5)*
        (M_PI + Complex(0,2)*log(2)))/9.);
    }

    /* else if(vpnorder == 6){
        return((-2137342/(45045.*sqrt(35)) + (2176238*Nu)/(6435.*sqrt(35)) - 
        (587516*pow(Nu,2))/(1053.*sqrt(35)) + 
        (452194*pow(Nu,3))/(3861.*sqrt(35)))*pow(x,4));
    } */

    else{
        return 0;
    }
}

static COMPLEX16 hl_4_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_4_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_4(Nu,vpnorder,x)) * cpolar(1,-4*Phi);
    }
}

static COMPLEX16 hl_4_m_min4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_min4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,4) * conj(hGO_4_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_4(Nu,vpnorder,x)) * cpolar(1,4*Phi); 
    }
}

//H43
static COMPLEX16 hGO_4_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;

    if(vpnorder == 3){
        return (Complex(0,0.16666666666666666)*delta*mass*(-1 + 2*Nu)*PhiDOT*(4*mass 
        + r*(23*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,10)*PhiDOT*r*rDOT - 2*pow(rDOT,2))))/(sqrt(70)*r);
    }   

    else if(vpnorder == 4){
         return ((Complex(0,-0.041666666666666664)*sqrt(0.35714285714285715)*
            pow(mass,2)*Nu*(4*mass + 
          r*(23*pow(PhiDOT,2)*pow(r,2) + Complex(0,10)*PhiDOT*r*rDOT - 
          2*pow(rDOT,2)))*((-1 + delta)*S1z + S2z + delta*S2z))/
           pow(r,3)
           - (Complex(0,-1.125)*sqrt(0.35714285714285715)*Nu*((-1 + delta)*S1z + S2z 
           + delta*S2z)*pow(x,3)) +(Complex(0,-1.125)
           *sqrt(0.35714285714285715)*Nu*(-S1z + S2z + delta*(S1z + S2z))*pow(x,3)));
    } 

    else if(vpnorder == 5){
        return (Complex(0,0.0012626262626262627)*delta*mass*PhiDOT*(2*pow(mass,2)*(972 
        - 2293*Nu + 1398*pow(Nu,2)) + 2*mass*r*((1788 - 9077*Nu 
        + 13416*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) + Complex(0,3)*(-2796 
        + 5299*Nu + 1622*pow(Nu,2))*PhiDOT*r*rDOT - 2*(-1200 + 2545*Nu 
        + 162*pow(Nu,2))*pow(rDOT,2)) - 3*pow(r,2)*((-524 - 489*Nu 
        + 6392*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) + Complex(0,4)*(796 - 1864*Nu 
        + 133*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT + 42*(-51 + 94*Nu 
        + 56*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,4)*(-229 + 366*Nu + 358*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        - 4*(-43 + 62*Nu + 80*pow(Nu,2))*pow(rDOT,4))))/(sqrt(70)*pow(r,2));
    }

    else if(vpnorder == 6){
        return ((delta*pow(mass,2)*Nu*PhiDOT*(6*mass*(181*PhiDOT*r 
        - Complex(0,89)*rDOT) + r*(4847*pow(PhiDOT,3)*pow(r,3) 
        - Complex(0,7338)*pow(PhiDOT,2)*pow(r,2)*rDOT - 408*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,112)*pow(rDOT,3))))/(180.*sqrt(70)*pow(r,2))+((Complex(0,0.03409090909090909)
        *(220*(S1z - S2z + delta*(S1z + S2z)) + Nu*(-2403*(S1z - S2z) + 203*delta*(S1z + S2z)) + 
       pow(Nu,2)*(3359*(S1z - S2z) + 457*delta*(S1z + S2z)))*pow(x,4))/sqrt(70)));
    }

    else if(vpnorder == 7){

        return((Complex(0,-0.020833333333333332)*pow(x,4.5)*(54*(10*(4 + delta)*pow(Nu,2) 
        + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) - 30*Nu + 28*pow(Nu,2)))*pow(S1z,2) + 
        S1z*(Complex(0,297) - Nu*(Complex(0,283) + 810*M_PI) + delta*(Complex(0,297) - 5*Nu*(Complex(0,629) - 162*M_PI + 216*Nu*S2z)) + 
        Complex(0,3240)*(-1 + delta)*Nu*atanh(1/5)) + S2z*(Complex(0,297)*(-1 + delta) + Nu*(Complex(0,283) + 810*M_PI) + 
        54*kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu - 28*pow(Nu,2))*S2z + 
        5*Nu*(Complex(0,-629)*delta + 162*delta*M_PI - 432*Nu*S2z + 108*delta*Nu*S2z + Complex(0,648)*(1 + delta)*atanh(1/5)))))/sqrt(70));
    }

    else{
        return 0;
    }
    
}


static COMPLEX16 hQC_4_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 3){
        return(((Complex(0,-4.5)*delta)/sqrt(70) + (Complex(0,9)*delta*Nu)/sqrt(70))*
         pow(x,2.5));
    } */

   /*  else if(vpnorder == 5){
        return(((Complex(0,15.954545454545455)*delta)/sqrt(70) - 
         Complex(0,6.170454545454546)*sqrt(0.7)*delta*Nu + 
         (Complex(0,17.863636363636363)*delta*pow(Nu,2))/sqrt(70))*
         pow(x,3.5));
    } */

    /* else */ if(vpnorder == 6){
        return((Complex(0,13.5)*delta*(-1 + 2*Nu)*pow(x,4)*
         (M_PI + Complex(0,2)*log(1.5)))/sqrt(70));
    }

    else{
        return 0;
    }
}

static COMPLEX16 hl_4_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_4_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_3(Nu,vpnorder,x)) * cpolar(1,-3*Phi);
    }

}

static COMPLEX16 hl_4_m_min3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_min3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,4) * conj(hGO_4_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_3(Nu,vpnorder,x)) * cpolar(1,3*Phi);
    }
    
}

//H42
static COMPLEX16 hGO_4_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;


    if(vpnorder == 2){
        return (-(sqrt(5)*(-1 + 3*Nu)*(7*pow(mass,2) - 6*pow(r,2)*(PhiDOT*r 
        - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,3) 
        + 3*mass*r*(pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,9)*PhiDOT*r*rDOT - 6*pow(rDOT,2))))/(126.*pow(r,2)));
    }

    else if(vpnorder == 4){
        return (-(40*pow(mass,3)*(314 - 987*Nu + 195*pow(Nu,2)) 
        + 60*(23 - 159*Nu + 291*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + pow(mass,2)*r*((1987 - 11200*Nu + 12960*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,12)*(967 - 4615*Nu + 5935*pow(Nu,2))*PhiDOT*r*rDOT 
        - 10*(290 - 2033*Nu + 4365*pow(Nu,2))*pow(rDOT,2)) 
        - 3*mass*pow(r,2)*((1577 - 7940*Nu + 9920*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,4)*(-454 - 315*Nu + 5980*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        - 2*(549 - 2140*Nu + 2140*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,2)*(-1853 + 1730*Nu + 13230*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        - 20*(-83 + 30*Nu + 762*pow(Nu,2))*pow(rDOT,4)))/(5544.*sqrt(5)*pow(r,3)));
    }

    else if(vpnorder == 5){
        return ((pow(mass,2)*Nu*(mass*(Complex(0,129)*PhiDOT*r - 54*rDOT) 
        + r*(Complex(0,-73)*pow(PhiDOT,3)*pow(r,3) + 21*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - Complex(0,78)*PhiDOT*r*pow(rDOT,2) + 26*pow(rDOT,3))))/(84.*sqrt(5)*pow(r,3))
        +((4*(21*delta*Nu*(S1z - S2z) + 59*Nu*(S1z + S2z) - 78*pow(Nu,2)*(S1z + S2z) 
        - 10*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,3.5))/(189.*sqrt(5))));
    }

    else if(vpnorder == 6){
        return ((-10*pow(mass,4)*(-4477296 + 12734393*Nu - 6895*pow(Nu,2) 
        + 1043805*pow(Nu,3)) + 3150*(-367 + 4337*Nu - 17462*pow(Nu,2) 
        + 23577*pow(Nu,3))*pow(r,4)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,3)*pow(PhiDOT*r + Complex(0,1)*rDOT,5) 
        - 2*pow(mass,3)*r*(7*(-100473 + 3430399*Nu - 9132990*pow(Nu,2) 
        + 2885660*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2) + Complex(0,2)*(7571073 
        - 10780154*Nu - 56898800*pow(Nu,2) + 43665510*pow(Nu,3))*PhiDOT*r*rDOT 
        - 10*(1283609 - 5800627*Nu + 3725295*pow(Nu,2) + 4771935*pow(Nu,3))*pow(rDOT,2)) 
        + 7*pow(mass,2)*pow(r,2)*((1071402 + 3846989*Nu - 27339110*pow(Nu,2) 
        + 17538420*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) + Complex(0,12)*(671714 
        - 1645932*Nu - 1903365*pow(Nu,2) + 3346250*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        - (481563 + 4291961*Nu - 17137220*pow(Nu,2) 
        + 9315720*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,12)*(52044 - 1084807*Nu + 1849450*pow(Nu,2) 
        + 4171730*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) - 5*(1083 - 1246819*Nu 
        + 2524240*pow(Nu,2) + 5995845*pow(Nu,3))*pow(rDOT,4)) 
        - 105*mass*pow(r,3)*((54272 - 58271*Nu - 815454*pow(Nu,2) 
        + 1435572*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) + Complex(0,1)*(73976 - 157355*Nu 
        - 811766*pow(Nu,2) + 2935488*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT 
        + (72637 - 400832*Nu + 282028*pow(Nu,2) 
        + 1063956*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,4)*(90174 - 385167*Nu - 126419*pow(Nu,2) 
        + 1739072*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        + (-55817 + 58920*Nu + 989942*pow(Nu,2) 
        - 2334016*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) + Complex(0,1)*(218975 
        - 1037408*Nu + 148970*pow(Nu,2) + 3699480*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) 
        + 10*(-10233 + 44864*Nu + 13050*pow(Nu,2) 
        - 203280*pow(Nu,3))*pow(rDOT,6)))/(5.04504e6*sqrt(5)*pow(r,4))+((2*sqrt(5)*(kappa1*(1 + delta - 2*Nu 
        + 6*pow(Nu,2))*pow(S1z,2) + S2z*(4*(1 - 3*Nu)*Nu*S1z - kappa2*(-1 + delta + 2*Nu - 6*pow(Nu,2))*S2z))*pow(x,4))/63.));
    }

    else if(vpnorder == 7){

        return(((-5958*delta*pow(Nu,2)*(S1z - S2z) + 45*delta*Nu*(-S1z + S2z) 
        - 6009*Nu*(S1z + S2z) + 7528*pow(Nu,2)*(S1z + S2z) - 1134*pow(Nu,3)*(S1z + S2z) + 
       1104*(S1z + delta*S1z + S2z - delta*S2z))*pow(x,4.5))/(2079.*sqrt(5)));
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hQC_4_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    /* if(vpnorder == 2){
        return(((2*sqrt(5))/63. - (2*sqrt(5)*Nu)/21.)*pow(x,2));
    }

    else if(vpnorder == 4){
        return((-437/(693.*sqrt(5)) + (115*sqrt(5)*Nu)/297. - 
        (19*sqrt(5)*pow(Nu,2))/693.)*pow(x,3));
    } */

    /* else */ if(vpnorder == 5){
        return(M_PI*((4*sqrt(5)*pow(x,3.5))/63. - (4*sqrt(5)*Nu*pow(x,3.5))/21.));
    }
    /* else if(vpnorder == 6){
        return((346013/(420420.*sqrt(5)) - (606751*Nu)/(180180.*sqrt(5)) + 
        (400453*pow(Nu,2))/(162162.*sqrt(5)) + 
        (25783*pow(Nu,3))/(108108.*sqrt(5)))*pow(x,4));
    } */
    else{
        return 0;
    }

}



static COMPLEX16 hl_4_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_4_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_2(Nu,vpnorder,x)) * cpolar(1,-2*Phi);
    }
}

static COMPLEX16 hl_4_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,4) * conj(hGO_4_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_4_m_2(Nu,vpnorder,x)) * cpolar(1,2*Phi); 
    }
}

//H41
static COMPLEX16 hGO_4_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;

    

    if(vpnorder == 3){
        return (delta*mass*(-1 + 2*Nu)*PhiDOT*(Complex(0,-12)*mass 
        + r*(Complex(0,11)*pow(PhiDOT,2)*pow(r,2) + 10*PhiDOT*r*rDOT 
        + Complex(0,6)*pow(rDOT,2))))/(42.*sqrt(10)*r);
    }

    else if(vpnorder == 4){
        return ((Complex(0,0.005952380952380952)*sqrt(2.5)*pow(mass,2)*Nu*
          (12*mass + r*(-11*pow(PhiDOT,2)*pow(r,2) + 
          Complex(0,10)*PhiDOT*r*rDOT - 6*pow(rDOT,2)))*
          ((-1 + delta)*S1z + S2z + delta*S2z))/pow(r,3)
          - (Complex(0,0.005952380952380952)*sqrt(2.5)*Nu*(-S1z + S2z 
          + delta*(S1z + S2z))*pow(x,3)) +(Complex(0,0.005952380952380952)
          *sqrt(2.5)*Nu*(-S1z + S2z + delta*(S1z + S2z))*pow(x,3)));
    }

    else if(vpnorder == 5){
        return (Complex(0,-0.00018037518037518038)*delta*mass*PhiDOT*(6*pow(mass,2)*(972 - 2293*Nu + 1398*pow(Nu,2)) 
        - 2*mass*r*((-340 - 785*Nu + 6504*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2) 
        - Complex(0,3)*(-2796 + 5299*Nu + 1622*pow(Nu,2))*PhiDOT*r*rDOT + 6*(-1200 
        + 2545*Nu + 162*pow(Nu,2))*pow(rDOT,2)) + 3*pow(r,2)*((-540 + 235*Nu + 2648*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,4)*(-1764 + 3536*Nu + 373*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 2*(-723 + 1022*Nu + 1384*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,4)*(-229 + 366*Nu + 358*pow(Nu,2))*PhiDOT*r*pow(rDOT,3) 
        + 12*(-43 + 62*Nu + 80*pow(Nu,2))*pow(rDOT,4))))/(sqrt(10)*pow(r,2));
    }

    else if(vpnorder == 6){
        return (-(delta*pow(mass,2)*Nu*PhiDOT*(mass*(362*PhiDOT*r - Complex(0,534)*rDOT) 
        + r*(149*pow(PhiDOT,3)*pow(r,3) + Complex(0,182)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - 136*PhiDOT*r*pow(rDOT,2) + Complex(0,112)*pow(rDOT,3))))/(420.*sqrt(10)*pow(r,2))
        +((Complex(0,-0.00018037518037518038)*(220*(S1z - S2z + delta*(S1z + S2z)) 
        + Nu*(-2247*(S1z - S2z) + 47*delta*(S1z + S2z)) + 
        pow(Nu,2)*(2891*(S1z - S2z) + 613*delta*(S1z + S2z)))*pow(x,4))/sqrt(10)));
    }

    else if(vpnorder == 7){

        return((Complex(0,0.000992063492063492)*pow(x,4.5)*(6*(10*(4 + delta)*pow(Nu,2)
         + kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) + 6*Nu*(-5 + 2*Nu)))*pow(S1z,2) + 
       S2z*(Complex(0,11)*(-1 + delta) + Nu*(Complex(0,-31) + 30*M_PI) + 6*kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) 
       + 6*(5 - 2*Nu)*Nu)*S2z + 15*Nu*(-4*(4*Nu*S2z + Complex(0,1)*log(2)) + delta*(2*M_PI 
       + 4*Nu*S2z - Complex(0,1)*(5 + log(16))))) + S1z*(Complex(0,11) + Nu*(Complex(0,31) 
       - 30*M_PI + Complex(0,60)*log(2)) + delta*(Complex(0,11) - 15*Nu*(-2*M_PI + 8*Nu*S2z + Complex(0,1)*(5 + log(16)))))))/sqrt(10));
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hQC_4_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 3){
        return(((Complex(0,0.023809523809523808)*delta)/sqrt(10) - 
        (Complex(0,0.047619047619047616)*delta*Nu)/sqrt(10))*pow(x,2.5));
    }
    else if(vpnorder == 5){
        return(((Complex(0,-0.07287157287157287)*delta)/sqrt(10) + 
        (Complex(0,0.18235930735930736)*delta*Nu)/sqrt(10) - 
        (Complex(0,0.05988455988455989)*delta*pow(Nu,2))/sqrt(10))*
        pow(x,3.5));
    } */
    /* else */ if(vpnorder == 6){
        return(M_PI*((Complex(0,0.023809523809523808)*delta*pow(x,4))/sqrt(10) - 
        (Complex(0,0.047619047619047616)*delta*Nu*pow(x,4))/sqrt(10))\
        + ((delta*pow(x,4))/(21.*sqrt(10)) - 
        (sqrt(0.4)*delta*Nu*pow(x,4))/21.)*log(2));
    }
    else{
        return 0;
    }
}


static COMPLEX16 hl_4_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_4_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) + hQC_4_m_1(Nu,vpnorder,x)) * cpolar(1,-1*Phi);
    }

}

static COMPLEX16 hl_4_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_4_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,4) * conj(hGO_4_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+ hQC_4_m_1(Nu,vpnorder,x)) * cpolar(1,1*Phi);
    }
    
}

//H55
static COMPLEX16 hGO_5_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;

    

    if(vpnorder == 3){
        return ((delta*(-1 + 2*Nu)*(24*pow(r,2)*pow(Complex(0,-1)*PhiDOT*r + rDOT,5) 
        + 2*pow(mass,2)*(Complex(0,-86)*PhiDOT*r + 41*rDOT) 
        + 3*mass*r*(Complex(0,-143)*pow(PhiDOT,3)*pow(r,3) + 208*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,132)*PhiDOT*r*pow(rDOT,2) - 32*pow(rDOT,3))))/(48.*sqrt(66)*pow(r,2)));
    }    

    else if(vpnorder == 5){
        return ((delta*(360*(33 - 197*Nu + 294*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,6)*(Complex(0,1)*PhiDOT*r + rDOT) + 2*pow(mass,3)*(Complex(0,-5)*(53311 - 121906*Nu 
        + 42816*pow(Nu,2))*PhiDOT*r + 78*(1141 - 2760*Nu + 1420*pow(Nu,2))*rDOT) 
        + 2*pow(mass,2)*r*(Complex(0,-10)*(40826 - 125981*Nu + 87534*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) 
        + (112966 - 818425*Nu + 1385970*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT + Complex(0,20)*(-2636 - 11335*Nu 
        + 43962*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) - 39*(-639 - 735*Nu + 5690*pow(Nu,2))*pow(rDOT,3)) 
        + 15*mass*pow(r,2)*(Complex(0,4)*(139 - 1687*Nu + 11376*pow(Nu,2))*pow(PhiDOT,5)*pow(r,5) 
        - 2*(11276 - 19559*Nu + 5982*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4)*rDOT + Complex(0,3)*(-14615 + 12440*Nu 
        + 37132*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) - 8*(-4666 + 139*Nu 
        + 22194*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - Complex(0,4)*(-3971 - 3226*Nu + 27360*pow(Nu,2))*PhiDOT*r*pow(rDOT,4) 
        + 48*(-57 - 97*Nu + 518*pow(Nu,2))*pow(rDOT,5))))/(18720.*sqrt(66)*pow(r,3)));
    }

    else if(vpnorder == 6){
        return (-(delta*pow(mass,2)*Nu*(3566*pow(mass,2) + 6*mass*r*(11305*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,3921)*PhiDOT*r*rDOT - 906*pow(rDOT,2)) + pow(r,2)*(104681*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,17192)*pow(PhiDOT,3)*pow(r,3)*rDOT - 27840*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,9968)*PhiDOT*r*pow(rDOT,3) + 1424*pow(rDOT,4))))/(1680.*sqrt(66)*pow(r,4))
        +((Complex(0,-21.70138888888889)*((1 + delta + 3*delta*(-1 + Nu)*Nu + Nu*(-7 + 11*Nu))*S1z 
        + (-1 + delta + (7 - 11*Nu)*Nu + 3*delta*(-1 + Nu)*Nu)*S2z)*pow(x,4))/sqrt(66)));
    }

    else if(vpnorder == 7){

        return((Complex(0,16.276041666666668)*(-1 + 2*Nu)*(kappa1*(-1 - delta + 2*(2 + delta)*Nu)*pow(S1z,2) 
        + S2z*(-4*delta*Nu*S1z + kappa2*(1 - delta + 2*(-2 + delta)*Nu)*S2z))*pow(x,4.5))/sqrt(66));
    }

    else{
        return 0;
    }

}

static COMPLEX16 hQC_5_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 3){
        return(((Complex(0,13.020833333333334)*delta)/sqrt(66) - 
        (Complex(0,26.041666666666668)*delta*Nu)/sqrt(66))*pow(x,2.5));
    }
    else if(vpnorder == 5){
        return(((Complex(0,-87.80715811965813)*delta)/sqrt(66) + 
        (Complex(0,229.7008547008547)*delta*Nu)/sqrt(66) - 
        Complex(0,42.73504273504273)*sqrt(0.06060606060606061)*delta*
        pow(Nu,2))*pow(x,3.5));
    } */
    /* else */ if(vpnorder == 6){
        return((3125*delta*(-1 + 2*Nu)*pow(x,4)*(Complex(0,-1)*M_PI + log(6.25)))/
        (48.*sqrt(66)));
    }
    else{
        return 0;
    }

}

static COMPLEX16 hl_5_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_5_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_5(Nu,vpnorder,x)) * cpolar(1,-5*Phi);
    }

}

static COMPLEX16 hl_5_m_min5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,5) * conj(hGO_5_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_5(Nu,vpnorder,x)) * cpolar(1,5*Phi);
    }
    
}

//H54
static COMPLEX16 hGO_5_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    
    if(vpnorder == 4){
        return (-(mass*(1 - 5*Nu + 5*pow(Nu,2))*PhiDOT*(mass*(82*PhiDOT*r + Complex(0,22)*rDOT) 
        + 3*r*(58*pow(PhiDOT,3)*pow(r,3) + Complex(0,33)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - 12*PhiDOT*r*pow(rDOT,2) - Complex(0,2)*pow(rDOT,3))))/(36.*sqrt(165)*r));
    }

    else if(vpnorder == 5){

        return((32*Nu*((-1 + delta + 2*Nu)*S1z - (1 + delta - 2*Nu)*S2z)*pow(x,3.5))/(3.*sqrt(165)));
    }

    else if(vpnorder == 6){
        return (mass*PhiDOT*(-4*pow(mass,2)*(26*(-5051 + 28623*Nu - 46305*pow(Nu,2) + 29470*pow(Nu,3))*PhiDOT*r 
        + Complex(0,7)*(5684 - 26697*Nu + 14225*pow(Nu,2) + 25355*pow(Nu,3))*rDOT) 
        - 4*mass*r*((-149157 + 1133006*Nu - 2731750*pow(Nu,2) + 2085685*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3) 
        + Complex(0,7)*(101118 - 491779*Nu + 402185*pow(Nu,2) + 172105*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - (363357 - 1825691*Nu + 1714720*pow(Nu,2) + 395255*pow(Nu,3))*PhiDOT*r*pow(rDOT,2) 
        - Complex(0,7)*(9717 - 48896*Nu + 46000*pow(Nu,2) + 10865*pow(Nu,3))*pow(rDOT,3)) 
        + 15*pow(r,2)*(4*(3449 - 6580*Nu - 56728*pow(Nu,2) + 115269*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5) 
        + Complex(0,7)*(-8128 + 45859*Nu - 62702*pow(Nu,2) + 13996*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + 4*(10125 - 47852*Nu + 26635*pow(Nu,2) + 44240*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        + Complex(0,182)*(127 - 548*Nu + 73*pow(Nu,2) + 816*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - 8*(1009 - 4060*Nu - 889*pow(Nu,2) 
        + 7952*pow(Nu,3))*PhiDOT*r*pow(rDOT,4) - Complex(0,28)*(45 - 172*Nu 
        - 85*pow(Nu,2) + 400*pow(Nu,3))*pow(rDOT,5))))/(65520.*sqrt(165)*pow(r,2));
    }

    else if(vpnorder == 7){

        return((16*(-530*pow(Nu,3)*(S1z + S2z) + 104*(S1z + delta*S1z + S2z - delta*S2z) 
        + 2*pow(Nu,2)*(541*delta*(S1z - S2z) - 120*(S1z + S2z)) - 
       Nu*(1139*delta*(S1z - S2z) + 109*(S1z + S2z)))*pow(x,4.5))/(117.*sqrt(165)));
    }

    else{
        return 0;
    }
}

/* static COMPLEX16 hQC_5_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((-64/(9.*sqrt(165)) + (64*sqrt(0.15151515151515152)*Nu)/9. - 
        (64*sqrt(0.15151515151515152)*pow(Nu,2))/9.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((142432/(4095.*sqrt(165)) - (10528*sqrt(0.7333333333333333)*Nu)/585. + 
        (33344*pow(Nu,2))/(117.*sqrt(165)) - 
        (3616*pow(Nu,3))/(39.*sqrt(165)))*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R) * (hGO_5_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_5_m_4(Nu,vpnorder,x) */) * cpolar(1,-4*Phi);
    }
}

static COMPLEX16 hl_5_m_min4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,5) * conj(hGO_5_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_5_m_4(Nu,vpnorder,x) */) * cpolar(1,4*Phi); 
    }
}

//H53
static COMPLEX16 hGO_5_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;


    if(vpnorder == 3){
        return (delta*(-1 + 2*Nu)*(2*pow(mass,2)*(Complex(0,258)*PhiDOT*r - 205*rDOT) 
        - Complex(0,120)*pow(r,2)*(PhiDOT*r - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + 3*mass*r*(Complex(0,-51)*pow(PhiDOT,3)*pow(r,3) - 240*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - Complex(0,396)*PhiDOT*r*pow(rDOT,2) + 160*pow(rDOT,3))))/(144.*sqrt(330)*pow(r,2));
    }

    else if(vpnorder == 5){
        return (delta*(Complex(0,120)*(33 - 197*Nu + 294*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,5) + 2*pow(mass,3)*(Complex(0,1)*(53311 
        - 121906*Nu + 42816*pow(Nu,2))*PhiDOT*r - 26*(1141 - 2760*Nu + 1420*pow(Nu,2))*rDOT) 
        + 2*pow(mass,2)*r*(Complex(0,-2)*(6350 - 12803*Nu + 10314*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) 
        + (-6546 + 64131*Nu - 109702*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT - Complex(0,4)*(-2636 
        - 11335*Nu + 43962*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) + 13*(-639 - 735*Nu + 5690*pow(Nu,2))*pow(rDOT,3)) 
        + 3*mass*pow(r,2)*(Complex(0,-4)*(1223 - 4567*Nu + 3396*pow(Nu,2))*pow(PhiDOT,5)*pow(r,5) + 26*(-412 
        - 437*Nu + 3114*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4)*rDOT + Complex(0,1)*(-2331 - 31496*Nu 
        + 93276*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) + 8*(-1994 + 1541*Nu 
        + 5718*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + Complex(0,4)*(-3971 - 3226*Nu + 27360*pow(Nu,2))*PhiDOT*r*pow(rDOT,4) 
        - 80*(-57 - 97*Nu + 518*pow(Nu,2))*pow(rDOT,5))))/(3744.*sqrt(330)*pow(r,3));
    }

    else if(vpnorder == 6){
        return ((delta*pow(mass,2)*Nu*(17830*pow(mass,2) + 18*mass*r*(5231*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,3921)*PhiDOT*r*rDOT - 1510*pow(rDOT,2)) - pow(r,2)*(48579*pow(PhiDOT,4)*pow(r,4)
        + Complex(0,31304)*pow(PhiDOT,3)*pow(r,3)*rDOT + 33024*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,29904)*PhiDOT*r*pow(rDOT,3) - 7120*pow(rDOT,4))))/(5040.*sqrt(330)*pow(r,4))
        +(Complex(0,0.1875)*sqrt(0.02727272727272727)*((5 + Nu*(-27 + 31*Nu) + delta*(5 + 23*(-1 + Nu)*Nu))*S1z
         + (-5 + (27 - 31*Nu)*Nu + delta*(5 + 23*(-1 + Nu)*Nu))*S2z)*pow(x,4)));
    }

    else if(vpnorder == 7){

        return(Complex(0,-0.140625)*sqrt(0.02727272727272727)*(kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) 
        - 30*Nu + 8*pow(Nu,2))*pow(S1z,2) + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z 
        + kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu - 8*pow(Nu,2))*S2z))*pow(x,4.5));
    }

    else{
        return 0;
    }
}

static COMPLEX16 hQC_5_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 3){
        return((Complex(0,-0.5625)*sqrt(0.02727272727272727)*delta + 
        Complex(0,1.125)*sqrt(0.02727272727272727)*delta*Nu)*pow(x,2.5));
    }
    else if(vpnorder == 5){
        return((Complex(0,2.985576923076923)*sqrt(0.02727272727272727)*delta - 
        Complex(0,6.6923076923076925)*sqrt(0.02727272727272727)*delta*Nu + 
        Complex(0,0.11538461538461539)*sqrt(3.3)*delta*pow(Nu,2))*
        pow(x,3.5));
    } */
    /* else */ if(vpnorder == 6){
        return(Complex(0,1.6875)*sqrt(0.02727272727272727)*delta*(-1 + 2*Nu)*
         pow(x,4)*(M_PI + Complex(0,2)*log(1.5)));
    }
    else{
        return 0;
    }
}


static COMPLEX16 hl_5_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_5_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_3(Nu,vpnorder,x)) * cpolar(1,-3*Phi);
    }

}

static COMPLEX16 hl_5_m_min3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,5) * conj(hGO_5_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_3(Nu,vpnorder,x)) * cpolar(1,3*Phi);
    }
    
}

//H52
static COMPLEX16 hGO_5_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);


    if(vpnorder == 4){
        return (mass*(1 - 5*Nu + 5*pow(Nu,2))*PhiDOT*(mass*(41*PhiDOT*r 
        + Complex(0,22)*rDOT) - 3*r*(11*pow(PhiDOT,3)*pow(r,3) 
        - Complex(0,3)*pow(PhiDOT,2)*pow(r,2)*rDOT + 6*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,2)*pow(rDOT,3))))/(54.*sqrt(55)*r);
    }

    else if(vpnorder == 5){

        return((2*Nu*(S1z - delta*S1z + S2z + delta*S2z - 2*Nu*(S1z + S2z))*pow(x,3.5))/(9.*sqrt(55)));
    }

    else if(vpnorder == 6){
        return (mass*PhiDOT*(4*pow(mass,2)*(13*(-5051 + 28623*Nu - 46305*pow(Nu,2) + 29470*pow(Nu,3))*PhiDOT*r 
        + Complex(0,7)*(5684 - 26697*Nu + 14225*pow(Nu,2) + 25355*pow(Nu,3))*rDOT) 
        - 2*mass*r*((23157 + 154*Nu - 648410*pow(Nu,2) + 1133195*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3) - Complex(0,14)*(6648 - 31729*Nu + 23795*pow(Nu,2) 
        + 13225*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*rDOT + (363357 - 1825691*Nu + 1714720*pow(Nu,2) 
        + 395255*pow(Nu,3))*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,14)*(9717 - 48896*Nu + 46000*pow(Nu,2) + 10865*pow(Nu,3))*pow(rDOT,3)) 
        + 15*pow(r,2)*(2*(2207 - 6076*Nu - 18424*pow(Nu,2) 
        + 38787*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5) - Complex(0,7)*(4744 - 23965*Nu + 23906*pow(Nu,2) 
        + 1892*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + 2*(5667 - 22652*Nu - 5747*pow(Nu,2) + 45416*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        - Complex(0,14)*(97 - 536*Nu 
        + 643*pow(Nu,2) + 36*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + 4*(1009 - 4060*Nu - 889*pow(Nu,2) 
        + 7952*pow(Nu,3))*PhiDOT*r*pow(rDOT,4) + Complex(0,28)*(45 - 172*Nu 
        - 85*pow(Nu,2) + 400*pow(Nu,3))*pow(rDOT,5))))/(98280.*sqrt(55)*pow(r,2));
    }

    else if(vpnorder == 7){

        return(((698*pow(Nu,3)*(S1z + S2z) - 104*(S1z + delta*S1z + S2z - delta*S2z)
         - 2*pow(Nu,2)*(457*delta*(S1z - S2z) + 48*(S1z + S2z)) + 
        Nu*(1055*delta*(S1z - S2z) + 193*(S1z + S2z)))*pow(x,4.5))/(351.*sqrt(55)));
    }
    else{
        return 0;
    }
}

/* static COMPLEX16 hQC_5_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((4/(27.*sqrt(55)) - (4*sqrt(0.45454545454545453)*Nu)/27. + 
        (4*sqrt(0.45454545454545453)*pow(Nu,2))/27.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((-7822/(12285.*sqrt(55)) + (6158*Nu)/(1755.*sqrt(55)) - 
        (1652*pow(Nu,2))/(351.*sqrt(55)) + 
        (14*sqrt(2.2)*pow(Nu,3))/117.)*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_5_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_5_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_5_m_2(Nu,vpnorder,x) */) * cpolar(1,-2*Phi);
    }
}

static COMPLEX16 hl_5_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,5) * conj(hGO_5_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_5_m_2(Nu,vpnorder,x) */) * cpolar(1,2*Phi); 
    }
}

//H51
static COMPLEX16 hGO_5_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    REAL8 kappa1 = 1.0;
    REAL8 kappa2 = 1.0;

    if(vpnorder == 3){
        return (delta*(-1 + 2*Nu)*(120*pow(r,2)*pow(Complex(0,1)*PhiDOT*r 
        - rDOT,3)*pow(PhiDOT*r - Complex(0,1)*rDOT,2) 
        + 2*pow(mass,2)*(Complex(0,-86)*PhiDOT*r + 205*rDOT) 
        + Complex(0,3)*mass*r*(97*pow(PhiDOT,3)*pow(r,3) 
        + Complex(0,160)*pow(PhiDOT,2)*pow(r,2)*rDOT + 132*PhiDOT*r*pow(rDOT,2) 
        + Complex(0,160)*pow(rDOT,3))))/(144.*sqrt(385)*pow(r,2));
    }

    else if(vpnorder == 5){
        return (delta*(-360*(33 - 197*Nu 
        + 294*pow(Nu,2))*pow(r,3)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,4)*pow(Complex(0,1)*PhiDOT*r + rDOT,3) 
        + 2*pow(mass,3)*(Complex(0,-1)*(53311 - 121906*Nu 
        + 42816*pow(Nu,2))*PhiDOT*r + 78*(1141 - 2760*Nu + 1420*pow(Nu,2))*rDOT) 
        + 2*pow(mass,2)*r*(Complex(0,2)*(29938 - 82195*Nu 
        + 59238*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3) + (-27026 + 120623*Nu 
        - 199326*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*rDOT + Complex(0,4)*(-2636 
        - 11335*Nu + 43962*pow(Nu,2))*PhiDOT*r*pow(rDOT,2) - 39*(-639 - 735*Nu 
        + 5690*pow(Nu,2))*pow(rDOT,3)) + 3*mass*pow(r,2)*(Complex(0,-4)*(3115 
        - 15385*Nu + 22098*pow(Nu,2))*pow(PhiDOT,5)*pow(r,5) 
        + 2*(-2108 - 17893*Nu + 56466*pow(Nu,2))*pow(PhiDOT,4)*pow(r,4)*rDOT 
        - Complex(0,3)*(-8473 - 9528*Nu 
        + 65204*pow(Nu,2))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) + 8*(-2692 
        - 6587*Nu + 29754*pow(Nu,2))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - Complex(0,4)*(-3971 - 3226*Nu + 27360*pow(Nu,2))*PhiDOT*r*pow(rDOT,4) 
        + 240*(-57 - 97*Nu 
        + 518*pow(Nu,2))*pow(rDOT,5))))/(11232.*sqrt(385)*pow(r,3));
    }

    else if(vpnorder == 6){
        return (-(delta*pow(mass,2)*Nu*(17830*pow(mass,2) 
        - 6*mass*r*(4723*pow(PhiDOT,2)*pow(r,2) - Complex(0,3921)*PhiDOT*r*rDOT 
        + 4530*pow(rDOT,2)) + pow(r,2)*(12629*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,24248)*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 20064*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,9968)*PhiDOT*r*pow(rDOT,3) 
        + 7120*pow(rDOT,4))))/(5040.*sqrt(385)*pow(r,4))
        +((Complex(0,-0.0023148148148148147)*((5 + Nu*(-23 + 19*Nu) 
        + delta*(5 + 27*(-1 + Nu)*Nu))*S1z + (-5 + (23 - 19*Nu)*Nu 
        + delta*(5 + 27*(-1 + Nu)*Nu))*S2z)*pow(x,4))/sqrt(385)));
    }

    else if(vpnorder == 7){

        return((Complex(0,0.001736111111111111)*(kappa1*(5 + 5*delta*pow(1 - 2*Nu,2) 
        - 2*Nu*(15 + 4*Nu))*pow(S1z,2) + S2z*(20*delta*(1 - 2*Nu)*Nu*S1z 
        + kappa2*(-5 + 5*delta*pow(1 - 2*Nu,2) + 30*Nu + 8*pow(Nu,2))*S2z))*pow(x,4.5))/sqrt(385));
    }

    else{
        return 0;
    }

}


static COMPLEX16 hQC_5_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    /* if(vpnorder == 3){
        return(((Complex(0,0.006944444444444444)*delta)/sqrt(385) - 
        (Complex(0,0.013888888888888888)*delta*Nu)/sqrt(385))*pow(x,2.5));
    }
    else if(vpnorder == 5){
        return(((Complex(0,-0.03187321937321937)*delta)/sqrt(385) + 
         Complex(0,0.005698005698005698)*sqrt(0.3142857142857143)*delta*
         Nu - (Complex(0,0.0007122507122507123)*delta*pow(Nu,2))/
         sqrt(385))*pow(x,3.5));
    } */
    /* else */ if(vpnorder == 6){
        return((Complex(0,-0.006944444444444444)*delta*(-1 + 2*Nu)*pow(x,4)*
        (M_PI - Complex(0,2)*log(2)))/sqrt(385));
    }
    else{
        return 0;
    }
}

static COMPLEX16 hl_5_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_5_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_1(Nu,vpnorder,x)) * cpolar(1,-1*Phi);
    }

}

static COMPLEX16 hl_5_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_5_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,5) * conj(hGO_5_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)+hQC_5_m_1(Nu,vpnorder,x)) * cpolar(1,1*Phi);
    }
    
}

//66
static COMPLEX16 hGO_6_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);
    

    if(vpnorder == 4){
        return (((1 - 5*Nu + 5*pow(Nu,2))*(172*pow(mass,3) 
        + 120*pow(r,3)*pow(PhiDOT*r + Complex(0,1)*rDOT,6) 
        + pow(mass,2)*r*(3269*pow(PhiDOT,2)*pow(r,2) + Complex(0,2920)*PhiDOT*r*rDOT 
        - 806*pow(rDOT,2)) + 15*mass*pow(r,2)*(281*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,494)*pow(PhiDOT,3)*pow(r,3)*rDOT - 444*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,208)*PhiDOT*r*pow(rDOT,3) + 40*pow(rDOT,4))))/(360.*sqrt(143)*pow(r,3)));
    }

    else if(vpnorder == 6){
        return ((14*pow(mass,4)*(-5740 + 29361*Nu - 33348*pow(Nu,2) 
        + 7334*pow(Nu,3)) - 7560*(-1 + 9*Nu - 26*pow(Nu,2) + 23*pow(Nu,3))*pow(r,4)*(PhiDOT*r 
        - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,7) 
        + 2*pow(mass,3)*r*((-539645 + 2950311*Nu - 4086684*pow(Nu,2) 
        + 1644517*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2) + Complex(0,6)*(-47984 + 275121*Nu 
        - 442540*pow(Nu,2) + 255850*pow(Nu,3))*PhiDOT*r*rDOT + 14*(3614 - 21621*Nu + 39684*pow(Nu,2) 
        - 29332*pow(Nu,3))*pow(rDOT,2)) + 3*pow(mass,2)*pow(r,2)*((-311847 + 1966993*Nu 
        - 3502751*pow(Nu,2) + 1752968*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) + Complex(0,4)*(629 
        + 160412*Nu - 846370*pow(Nu,2) + 912975*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        - 3*(65519 - 144403*Nu - 684796*pow(Nu,2) + 1205253*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,32)*(3921 - 11389*Nu - 27265*pow(Nu,2) + 58450*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) 
        + 14*(1867 - 5501*Nu - 12824*pow(Nu,2) + 28137*pow(Nu,3))*pow(rDOT,4)) 
        - 45*mass*pow(r,3)*((195 + 3619*Nu - 36617*pow(Nu,2) 
        + 66836*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) + Complex(0,4)*(-1878 + 10969*Nu 
        - 20741*pow(Nu,2) + 18263*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT + (17169 
        - 75446*Nu + 35497*pow(Nu,2) + 47054*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,2)*(9183 - 30296*Nu - 37835*pow(Nu,2) + 95060*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        - 4*(2781 - 6062*Nu - 28595*pow(Nu,2) + 49070*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        - Complex(0,16)*(228 - 217*Nu - 3871*pow(Nu,2) 
        + 5803*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) + 56*(9 + 4*Nu 
        - 221*pow(Nu,2) + 308*pow(Nu,3))*pow(rDOT,6)))/(15120.*sqrt(143)*pow(r,4)));
    }

    else if(vpnorder == 7){

        return((108*(110*pow(Nu,3)*(S1z + S2z) - 14*(S1z + delta*S1z + S2z - delta*S2z) 
        - 50*pow(Nu,2)*(2*delta*(S1z - S2z) + 3*(S1z + S2z)) +  Nu*(85*delta*(S1z - S2z) 
        + 83*(S1z + S2z)))*pow(x,4.5))/(35.*sqrt(143)));
    }

    else{
        return 0;
    }
   
}

/* static COMPLEX16 hQC_6_m_6(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((108/(5.*sqrt(143)) - (108*Nu)/sqrt(143) + (108*pow(Nu,2))/sqrt(143))*
        pow(x,3));
    }
    else if(vpnorder == 6){
        return((-6102/(35.*sqrt(143)) + (378*sqrt(1.1818181818181819)*Nu)/5. - 
        (6912*pow(Nu,2))/(5.*sqrt(143)) + 
        (162*sqrt(1.1818181818181819)*pow(Nu,3))/5.)*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_6: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_6(Nu,vpnorder,x) */) * cpolar(1,-6*Phi);
    }

}

static COMPLEX16 hl_6_m_min6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min6: Input PN order parameter should be between [0, 7].");
    }

   else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_6(Nu,vpnorder,x) */) * cpolar(1,6*Phi); 
    }
}

//65
static COMPLEX16 hGO_6_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
     REAL8 delta = sqrt(1-4*Nu);

    

    if(vpnorder == 5){
        return (Complex(0,0.003968253968253968)*delta*mass*(1 - 4*Nu
         + 3*pow(Nu,2))*PhiDOT*(82*pow(mass,2) 
        + 2*mass*r*(701*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,343)*PhiDOT*r*rDOT - 62*pow(rDOT,2)) 
        + 3*pow(r,2)*(547*pow(PhiDOT,4)*pow(r,4) + Complex(0,364)*pow(PhiDOT,3)*pow(r,3)*rDOT
         - 180*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) - Complex(0,56)*PhiDOT*r*pow(rDOT,3) 
        + 8*pow(rDOT,4))))/(sqrt(429)*pow(r,2));
    }

    else if(vpnorder == 6){

        return((Complex(0,-21.70138888888889)*Nu*(S1z + delta*(-1 + Nu)*S1z 
        - 3*Nu*S1z - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(429));
    }

    else{
        return 0;
    }
   
}

/* static COMPLEX16 hQC_6_m_5(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return(((Complex(0,12.40079365079365)*delta)/sqrt(429) - 
        (Complex(0,49.6031746031746)*delta*Nu)/sqrt(429) + 
        (Complex(0,37.20238095238095)*delta*pow(Nu,2))/sqrt(429))*
        pow(x,3.5));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_m5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_5(Nu,vpnorder,x) */) * cpolar(1,-5*Phi);
    }

}

static COMPLEX16 hl_6_m_min5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_5(Nu,vpnorder,x) */) * cpolar(1,5*Phi);
    }
    
}

//64
static COMPLEX16 hGO_6_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);


    if(vpnorder == 4){
        return (((1 - 5*Nu + 5*pow(Nu,2))*(-516*pow(mass,3) 
        + 360*pow(r,3)*(PhiDOT*r - Complex(0,1)*rDOT)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,5) + pow(mass,2)*r*(-3587*pow(PhiDOT,2)*pow(r,2) 
        - Complex(0,5840)*PhiDOT*r*rDOT + 2418*pow(rDOT,2)) 
        + 15*mass*pow(r,2)*(113*pow(PhiDOT,4)*pow(r,4) 
        - Complex(0,108)*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 468*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,416)*PhiDOT*r*pow(rDOT,3) - 120*pow(rDOT,4))))/(1980.*sqrt(78)*pow(r,3)));
    }

    else if(vpnorder == 6){
        return (-(14*pow(mass,4)*(-5740 + 29361*Nu - 33348*pow(Nu,2) + 7334*pow(Nu,3)) 
        + 7560*(-1 + 9*Nu - 26*pow(Nu,2) + 23*pow(Nu,3))*pow(r,4)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,6) + 2*pow(mass,3)*r*((-196625 
        + 1082991*Nu - 1522164*pow(Nu,2) + 618457*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,4)*(-47984 + 275121*Nu - 442540*pow(Nu,2) + 255850*pow(Nu,3))*PhiDOT*r*rDOT 
        + 14*(3614 - 21621*Nu + 39684*pow(Nu,2) - 29332*pow(Nu,3))*pow(rDOT,2)) 
        + pow(mass,2)*pow(r,2)*((133599 - 779681*Nu + 1417087*pow(Nu,2) 
        - 1130416*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) + Complex(0,8)*(3849 
        + 4172*Nu - 80290*pow(Nu,2) + 64435*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + (-226971 + 551047*Nu + 2049124*pow(Nu,2) 
        - 3713857*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,64)*(3921 - 11389*Nu - 27265*pow(Nu,2) + 58450*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) 
        + 42*(1867 - 5501*Nu - 12824*pow(Nu,2) + 28137*pow(Nu,3))*pow(rDOT,4)) 
        + 15*mass*pow(r,3)*((2267 - 12733*Nu + 13895*pow(Nu,2) 
        + 6300*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) - Complex(0,8)*(908 - 2597*Nu - 5873*pow(Nu,2) 
        + 11809*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT + (5241 + 10066*Nu - 173159*pow(Nu,2) 
        + 235382*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) + Complex(0,4)*(-1651 + 11312*Nu 
        - 25417*pow(Nu,2) + 20916*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        + 4*(3127 - 8386*Nu - 23569*pow(Nu,2) + 45122*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,32)*(228 - 217*Nu - 3871*pow(Nu,2) 
        + 5803*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) - 168*(9 + 4*Nu 
        - 221*pow(Nu,2) + 308*pow(Nu,3))*pow(rDOT,6)))/(27720.*sqrt(78)*pow(r,4)));
    }

    else if(vpnorder == 7){

        return((256*sqrt(0.05128205128205128)*(-150*pow(Nu,3)*(S1z + S2z) 
        + 14*(S1z + delta*S1z + S2z - delta*S2z) + 10*pow(Nu,2)*(6*delta*(S1z - S2z) + 23*(S1z + S2z)) - 
        Nu*(65*delta*(S1z - S2z) + 103*(S1z + S2z)))*pow(x,4.5))/3465.);
    }

    else{
        return 0;
    }

}

/* static COMPLEX16 hQC_6_m_4(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return(((-256*sqrt(0.05128205128205128))/495. + 
        (256*sqrt(0.05128205128205128)*Nu)/99. - 
        (256*sqrt(0.05128205128205128)*pow(Nu,2))/99.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return(((3968*sqrt(0.05128205128205128))/1155. - 
        (9088*sqrt(0.05128205128205128)*Nu)/495. + 
        (1024*sqrt(0.05128205128205128)*pow(Nu,2))/45. - 
        (2432*sqrt(0.05128205128205128)*pow(Nu,3))/495.)*pow(x,4));
    }
    else{
        return 0;
    }

} */

static COMPLEX16 hl_6_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_4(Nu,vpnorder,x) */) * cpolar(1,-4*Phi);
    }

}

static COMPLEX16 hl_6_m_min4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min4: Input PN order parameter should be between [0, 7].");
    }

   else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_4(Nu,vpnorder,x) */) * cpolar(1,4*Phi); 
    }
}

//63
static COMPLEX16 hGO_6_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);


    if(vpnorder == 5){
        return (Complex(0,-0.00036075036075036075)*delta*mass*(1 - 4*Nu 
        + 3*pow(Nu,2))*PhiDOT*(410*pow(mass,2) + 2*mass*r*(929*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,1029)*PhiDOT*r*rDOT - 310*pow(rDOT,2)) - 3*pow(r,2)*(513*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,28)*pow(PhiDOT,3)*pow(r,3)*rDOT + 228*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,168)*PhiDOT*r*pow(rDOT,3) - 40*pow(rDOT,4))))/(sqrt(65)*pow(r,2));
    }

    else if(vpnorder == 6){

        return((Complex(0,0.4602272727272727)*Nu*(S1z + delta*(-1 + Nu)*S1z - 3*Nu*S1z 
        - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(65));
    }

    else{
        return 0;
    }
   
}

/* static COMPLEX16 hQC_6_m_3(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,-0.262987012987013)*delta*(1 - 4*Nu + 3*pow(Nu,2))*
        pow(x,3.5))/sqrt(65));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_3(Nu,vpnorder,x) */) * cpolar(1,-3*Phi);
    }

}

static COMPLEX16 hl_6_m_min3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
     
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_3(Nu,vpnorder,x) */) * cpolar(1,3*Phi);
    }
    
}

//62
static COMPLEX16 hGO_6_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

   

    if(vpnorder == 4){
        return ((1 - 5*Nu + 5*pow(Nu,2))*(516*pow(mass,3)
        + 360*pow(r,3)*pow(PhiDOT*r - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + pow(mass,2)*r*(-145*pow(PhiDOT,2)*pow(r,2) + Complex(0,2920)*PhiDOT*r*rDOT 
        - 2418*pow(rDOT,2)) - 3*mass*pow(r,2)*(233*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,1050)*pow(PhiDOT,3)*pow(r,3)*rDOT - 252*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,1040)*PhiDOT*r*pow(rDOT,3) - 600*pow(rDOT,4))))/(2376.*sqrt(65)*pow(r,3));
    }

    else if(vpnorder == 6){
        return (14*pow(mass,4)*(-5740 + 29361*Nu - 33348*pow(Nu,2) + 7334*pow(Nu,3)) 
        - 7560*(-1 + 9*Nu - 26*pow(Nu,2) + 23*pow(Nu,3))*pow(r,4)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,3)*pow(PhiDOT*r + Complex(0,1)*rDOT,5) + 2*pow(mass,3)*r*((9187 
        - 37401*Nu + 16548*pow(Nu,2) + 2821*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,2)*(-47984 + 275121*Nu - 442540*pow(Nu,2) + 255850*pow(Nu,3))*PhiDOT*r*rDOT 
        + 14*(3614 - 21621*Nu + 39684*pow(Nu,2) - 29332*pow(Nu,3))*pow(rDOT,2)) 
        + pow(mass,2)*pow(r,2)*((54699 - 336749*Nu + 596995*pow(Nu,2) 
        - 337960*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4) - Complex(0,4)*(-5781 + 89572*Nu 
        - 379358*pow(Nu,2) + 444689*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*rDOT + (-9351 
        + 101899*Nu - 419300*pow(Nu,2) + 566195*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,32)*(3921 - 11389*Nu - 27265*pow(Nu,2) + 58450*pow(Nu,3))*PhiDOT*r*pow(rDOT,3) 
        + 42*(1867 - 5501*Nu - 12824*pow(Nu,2) + 28137*pow(Nu,3))*pow(rDOT,4)) 
        + 3*mass*pow(r,3)*((-7885 + 64211*Nu - 170905*pow(Nu,2) 
        + 146580*pow(Nu,3))*pow(PhiDOT,6)*pow(r,6) + Complex(0,4)*(1438 
        + 9779*Nu - 86023*pow(Nu,2) + 109949*pow(Nu,3))*pow(PhiDOT,5)*pow(r,5)*rDOT 
        + (16353 - 68054*Nu + 10297*pow(Nu,2) + 77294*pow(Nu,3))*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,2)*(14341 - 392*Nu - 316841*pow(Nu,2) 
        + 452508*pow(Nu,3))*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        - 4*(13 + 12530*Nu - 68803*pow(Nu,2) + 80654*pow(Nu,3))*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,80)*(228 - 217*Nu - 3871*pow(Nu,2) + 5803*pow(Nu,3))*PhiDOT*r*pow(rDOT,5) 
        - 840*(9 + 4*Nu - 221*pow(Nu,2) + 308*pow(Nu,3))*pow(rDOT,6)))/(33264.*sqrt(65)*pow(r,4));
    }

    else if(vpnorder == 7){

        return((4*(174*pow(Nu,3)*(S1z + S2z) - 14*(S1z + delta*S1z + S2z - delta*S2z) 
        + Nu*(53*delta*(S1z - S2z) + 115*(S1z + S2z)) - 2*pow(Nu,2)*(18*delta*(S1z - S2z) 
        + 139*(S1z + S2z)))*pow(x,4.5))/(2079.*sqrt(65)));
    }

    else{
        return 0;
    }
    
}

/* static COMPLEX16 hQC_6_m_2(REAL8 Nu, UINT4 vpnorder, REAL8 x){

    if(vpnorder == 4){
        return((4/(297.*sqrt(65)) - (4*sqrt(0.38461538461538464)*Nu)/297. + 
        (4*sqrt(0.38461538461538464)*pow(Nu,2))/297.)*pow(x,3));
    }
    else if(vpnorder == 6){
        return((-6/(77.*sqrt(65)) + (118*Nu)/(297.*sqrt(65)) - 
        (128*pow(Nu,2))/(297.*sqrt(65)) + 
        (14*pow(Nu,3))/(297.*sqrt(65)))*pow(x,4));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_2(Nu,vpnorder,x) */) * cpolar(1,-2*Phi);
    }

}

static COMPLEX16 hl_6_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min2: Input PN order parameter should be between [0, 7].");
    }

   else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_2(Nu,vpnorder,x) */) * cpolar(1,2*Phi); 
    }
}

//61
static COMPLEX16 hGO_6_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);


    if(vpnorder == 5){
        return (Complex(0,0.0002405002405002405)*delta*mass*(1 - 4*Nu + 3*pow(Nu,2))*PhiDOT*(410*pow(mass,2) 
        - 2*mass*r*(359*pow(PhiDOT,2)*pow(r,2) - Complex(0,343)*PhiDOT*r*rDOT + 310*pow(rDOT,2)) 
        + 3*pow(r,2)*(103*pow(PhiDOT,4)*pow(r,4) - Complex(0,196)*pow(PhiDOT,3)*pow(r,3)*rDOT 
        + 108*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,56)*PhiDOT*r*pow(rDOT,3) + 40*pow(rDOT,4))))/(sqrt(26)*pow(r,2));
    }

    else if(vpnorder == 6){

        return((Complex(0,-0.00042087542087542086)*Nu*(S1z + delta*(-1 + Nu)*S1z - 3*Nu*S1z 
        - (1 + delta)*S2z + (3 + delta)*Nu*S2z)*pow(x,4))/sqrt(26));
    }

    else{
        return 0;
    }
    
}

/* static COMPLEX16 hQC_6_m_1(REAL8 Nu, UINT4 vpnorder, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return((Complex(0,0.0002405002405002405)*delta*(1 - 4*Nu + 3*pow(Nu,2))*
        pow(x,3.5))/sqrt(26));
    }
    else{
        return 0;
    }
} */

static COMPLEX16 hl_6_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * (hGO_6_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_1(Nu,vpnorder,x) */) * cpolar(1,-1*Phi);
    }
}

static COMPLEX16 hl_6_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_6_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,6) * conj(hGO_6_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)/* +hQC_6_m_1(Nu,vpnorder,x) */) * cpolar(1,1*Phi);
    }
}

//77
static COMPLEX16 hGO_7_m_7(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return (delta*(1 - 4*Nu + 3*pow(Nu,2))*(720*pow(r,3)*pow(Complex(0,1)*PhiDOT*r - rDOT,7) 
        + 2*pow(mass,3)*(Complex(0,-4559)*PhiDOT*r + 1976*rDOT) 
        + 18*pow(mass,2)*r*(Complex(0,-3317)*pow(PhiDOT,3)*pow(r,3) + 4143*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,2178)*PhiDOT*r*pow(rDOT,2) - 442*pow(rDOT,3)) 
        + 45*mass*pow(r,2)*(Complex(0,-1069)*pow(PhiDOT,5)*pow(r,5) + 2112*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + Complex(0,2370)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) - 1600*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - Complex(0,600)*PhiDOT*r*pow(rDOT,4) + 96*pow(rDOT,5))))/(720.*sqrt(6006)*pow(r,3));
    }

    else{
        return 0;
    }
}

static COMPLEX16 hl_7_m_7(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_7: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_7(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-7*Phi);
    }

}

static COMPLEX16 hl_7_m_min7(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min7: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,7) * conj(hGO_7_m_7(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,7*Phi);
    }

}

//75
static COMPLEX16 hGO_7_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return ((delta*(1 - 4*Nu + 3*pow(Nu,2))*(2*pow(mass,3)*(Complex(0,22795)*PhiDOT*r 
        - 13832*rDOT) - Complex(0,5040)*pow(r,3)*(PhiDOT*r - Complex(0,1)*rDOT)*pow(PhiDOT*r 
        + Complex(0,1)*rDOT,6) + 18*pow(mass,2)*r*(Complex(0,5105)*pow(PhiDOT,3)*pow(r,3) 
        - 13041*pow(PhiDOT,2)*pow(r,2)*rDOT - Complex(0,10890)*PhiDOT*r*pow(rDOT,2) 
        + 3094*pow(rDOT,3)) + 45*mass*pow(r,2)*(Complex(0,-1207)*pow(PhiDOT,5)*pow(r,5) 
        + 336*pow(PhiDOT,4)*pow(r,4)*rDOT - Complex(0,3114)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        + 4928*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + Complex(0,3000)*PhiDOT*r*pow(rDOT,4) - 672*pow(rDOT,5))))/(65520.*sqrt(66)*pow(r,3)));
    }

    else{
        return 0;
    }
}

static COMPLEX16 hl_7_m_5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-5*Phi);
    }
}

static COMPLEX16 hl_7_m_min5(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min5: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,7) * conj(hGO_7_m_5(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,5*Phi);
    }

}

//73
static COMPLEX16 hGO_7_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return (delta*(1 - 4*Nu + 3*pow(Nu,2))*(5040*pow(r,3)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,2)*pow(Complex(0,-1)*PhiDOT*r + rDOT,5) + 2*pow(mass,3)*(Complex(0,-13677)*PhiDOT*r 
        + 13832*rDOT) + 18*pow(mass,2)*r*(Complex(0,1529)*pow(PhiDOT,3)*pow(r,3) + 2401*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + Complex(0,6534)*PhiDOT*r*pow(rDOT,2) - 3094*pow(rDOT,3)) 
        + 15*mass*pow(r,2)*(Complex(0,179)*pow(PhiDOT,5)*pow(r,5) - 4368*pow(PhiDOT,4)*pow(r,4)*rDOT 
        - Complex(0,4878)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) - 2240*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - Complex(0,5400)*PhiDOT*r*pow(rDOT,4) + 2016*pow(rDOT,5))))/(240240.*sqrt(6)*pow(r,3));
    }

    else{
        return 0;
    }
}

static COMPLEX16 hl_7_m_3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-3*Phi);
    }

}

static COMPLEX16 hl_7_m_min3(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min3: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,7) * conj(hGO_7_m_3(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,3*Phi);
    }
}

//71
static COMPLEX16 hGO_7_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){
    REAL8 delta = sqrt(1-4*Nu);

    if(vpnorder == 5){
        return (delta*(1 - 4*Nu + 3*pow(Nu,2))*(2*pow(mass,3)*(Complex(0,4559)*PhiDOT*r - 13832*rDOT) 
        - Complex(0,5040)*pow(r,3)*pow(PhiDOT*r - Complex(0,1)*rDOT,3)*pow(PhiDOT*r + Complex(0,1)*rDOT,4) 
        + 18*pow(mass,2)*r*(Complex(0,-1275)*pow(PhiDOT,3)*pow(r,3) + 2919*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - Complex(0,2178)*PhiDOT*r*pow(rDOT,2) + 3094*pow(rDOT,3)) 
        + Complex(0,27)*mass*pow(r,2)*(699*pow(PhiDOT,5)*pow(r,5) + Complex(0,1120)*pow(PhiDOT,4)*pow(r,4)*rDOT 
        + 1874*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) + Complex(0,2240)*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + 1000*PhiDOT*r*pow(rDOT,4) 
        + Complex(0,1120)*pow(rDOT,5))))/(432432.*sqrt(2)*pow(r,3));}

    else{
        return 0;
    }
}

static COMPLEX16 hl_7_m_1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-1*Phi);
    }

}

static COMPLEX16 hl_7_m_min1(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min1: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,7) * conj(hGO_7_m_1(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,1*Phi);
    }
}

//72
static COMPLEX16 hGO_7_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    

    if(vpnorder == 6){
        return (-(mass*(-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*PhiDOT*(8*pow(mass,2)*(494*PhiDOT*r 
        + Complex(0,411)*rDOT) - 12*mass*r*(530*pow(PhiDOT,3)*pow(r,3) - Complex(0,6)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        + 453*PhiDOT*r*pow(rDOT,2) + Complex(0,197)*pow(rDOT,3)) + 3*pow(r,2)*(824*pow(PhiDOT,5)*pow(r,5) 
        - Complex(0,671)*pow(PhiDOT,4)*pow(r,4)*rDOT + 864*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        + Complex(0,44)*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        + 320*PhiDOT*r*pow(rDOT,4) + Complex(0,120)*pow(rDOT,5))))/(96096.*sqrt(3)*pow(r,2)));
    }

    else if(vpnorder == 7){

        return((4*Nu*(S1z - delta*S1z + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z 
        - 2*(2 + delta)*Nu*S2z + 2*pow(Nu,2)*(S1z + S2z))*pow(x,4.5))/(3003.*sqrt(3)));
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hl_7_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) * cpolar(1,-2*Phi);
    }

}

static COMPLEX16 hl_7_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,7) * conj(hGO_7_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)) * cpolar(1,2*Phi);
    }

}

//74
static COMPLEX16 hGO_7_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
   REAL8 delta = sqrt(1-4*Nu);

    

    if(vpnorder == 6){
        return ((mass*(-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*PhiDOT*(8*pow(mass,2)*(988*PhiDOT*r 
        + Complex(0,411)*rDOT) + 12*mass*r*(844*pow(PhiDOT,3)*pow(r,3) + Complex(0,1518)*pow(PhiDOT,2)*pow(r,2)*rDOT 
        - 906*PhiDOT*r*pow(rDOT,2) - Complex(0,197)*pow(rDOT,3)) - 15*pow(r,2)*(656*pow(PhiDOT,5)*pow(r,5) 
        + Complex(0,179)*pow(PhiDOT,4)*pow(r,4)*rDOT + 192*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        + Complex(0,260)*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) 
        - 128*PhiDOT*r*pow(rDOT,4) - Complex(0,24)*pow(rDOT,5))))/(21840.*sqrt(66)*pow(r,2)));
    }

    else if(vpnorder == 7){

        return((-512*sqrt(0.06060606060606061)*Nu*(S1z - delta*S1z + 2*(-2 + delta)*Nu*S1z 
        + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z + 2*pow(Nu,2)*(S1z + S2z))*pow(x,4.5))/1365.);
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hl_7_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) * cpolar(1,-4*Phi);
    }

}

static COMPLEX16 hl_7_m_min4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,7) * conj(hGO_7_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)) * cpolar(1,4*Phi);
    }
    
}

//76
static COMPLEX16 hGO_7_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    REAL8 delta = sqrt(1-4*Nu);

    

    if(vpnorder == 6){
        return (-(mass*(-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*PhiDOT*(8*pow(mass,2)*(494*PhiDOT*r 
        + Complex(0,137)*rDOT) + 4*mass*r*(6026*pow(PhiDOT,3)*pow(r,3) 
        + Complex(0,4038)*pow(PhiDOT,2)*pow(r,2)*rDOT - 1359*PhiDOT*r*pow(rDOT,2) 
        - Complex(0,197)*pow(rDOT,3)) + 15*pow(r,2)*(1240*pow(PhiDOT,5)*pow(r,5) 
        + Complex(0,911)*pow(PhiDOT,4)*pow(r,4)*rDOT - 544*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,2) 
        - Complex(0,236)*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,3) + 64*PhiDOT*r*pow(rDOT,4) 
        + Complex(0,8)*pow(rDOT,5))))/(3360.*sqrt(429)*pow(r,2)));
    }
    else if(vpnorder == 7){

        return((324*sqrt(0.02097902097902098)*Nu*(S1z - delta*S1z 
        + 2*(-2 + delta)*Nu*S1z + S2z + delta*S2z - 2*(2 + delta)*Nu*S2z 
        + 2*pow(Nu,2)*(S1z + S2z))*pow(x,4.5))/35.);
    }

    else{
        return 0;
    }
    
}

static COMPLEX16 hl_7_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_6: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_7_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x) * cpolar(1,-6*Phi);
    }

}

static COMPLEX16 hl_7_m_min6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder,REAL8 S1z, REAL8 S2z, REAL8 x){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_7_m_min6: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,7) * conj(hGO_7_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder,S1z,S2z,x)) * cpolar(1,6*Phi);
    }
}

//88
static COMPLEX16 hGO_8_m_8(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){

    if(vpnorder == 6){
        return (((-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*(9118*pow(mass,4) 
        + 5040*pow(r,4)*pow(PhiDOT*r + Complex(0,1)*rDOT,8) + 4*pow(mass,3)*r*(82543*pow(PhiDOT,2)*pow(r,2) 
        + Complex(0,67760)*PhiDOT*r*rDOT - 16717*pow(rDOT,2)) + 9*pow(mass,2)*pow(r,2)*(124583*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,192640)*pow(PhiDOT,3)*pow(r,3)*rDOT - 144024*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,56280)*PhiDOT*r*pow(rDOT,3) + 9188*pow(rDOT,4)) + 315*mass*pow(r,3)*(2005*pow(PhiDOT,6)*pow(r,6) 
        + Complex(0,4250)*pow(PhiDOT,5)*pow(r,5)*rDOT - 5538*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        - Complex(0,4760)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) + 2600*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,816)*PhiDOT*r*pow(rDOT,5) 
        - 112*pow(rDOT,6))))/(2016.*sqrt(170170)*pow(r,4)));}

    else{
        return 0;
    }
}

static COMPLEX16 hl_8_m_8(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_8: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_8_m_8(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-8*Phi);
    }
}

static COMPLEX16 hl_8_m_min8(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_min8: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * pow(-1,8) * conj(hGO_8_m_8(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,8*Phi);
    }
  
}

//86
static COMPLEX16 hGO_8_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){

    if(vpnorder == 6){
        return (-((-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*(18236*pow(mass,4) 
        - 10080*pow(r,4)*(PhiDOT*r - Complex(0,1)*rDOT)*pow(PhiDOT*r + Complex(0,1)*rDOT,7) 
        + 8*pow(mass,3)*r*(42923*pow(PhiDOT,2)*pow(r,2) + Complex(0,50820)*PhiDOT*r*rDOT 
        - 16717*pow(rDOT,2)) + 18*pow(mass,2)*pow(r,2)*(15663*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,58240)*pow(PhiDOT,3)*pow(r,3)*rDOT - 74094*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        - Complex(0,42210)*PhiDOT*r*pow(rDOT,3) + 9188*pow(rDOT,4)) 
        - 315*mass*pow(r,3)*(678*pow(PhiDOT,6)*pow(r,6) 
        + Complex(0,681)*pow(PhiDOT,5)*pow(r,5)*rDOT + 852*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,2660)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) 
        - 2640*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        - Complex(0,1224)*PhiDOT*r*pow(rDOT,5) + 224*pow(rDOT,6))))/(10080.*sqrt(51051)*pow(r,4)));}

    else{
        return 0;
    }
}

static COMPLEX16 hl_8_m_6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_6: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_8_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-6*Phi);
    }

}

static COMPLEX16 hl_8_m_min6(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_min6: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,8) * conj(hGO_8_m_6(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,6*Phi);
    }
    
}

//84
static COMPLEX16 hGO_8_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){

    if(vpnorder == 6){
        return (((-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*(9118*pow(mass,4) 
        + 5040*pow(r,4)*pow(PhiDOT*r - Complex(0,1)*rDOT,2)*pow(PhiDOT*r + Complex(0,1)*rDOT,6) 
        + 4*pow(mass,3)*r*(14623*pow(PhiDOT,2)*pow(r,2) + Complex(0,33880)*PhiDOT*r*rDOT 
        - 16717*pow(rDOT,2)) - 9*pow(mass,2)*pow(r,2)*(8377*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,2240)*pow(PhiDOT,3)*pow(r,3)*rDOT + 24144*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,28140)*PhiDOT*r*pow(rDOT,3) - 9188*pow(rDOT,4)) + 45*mass*pow(r,3)*(243*pow(PhiDOT,6)*pow(r,6) 
        - Complex(0,1701)*pow(PhiDOT,5)*pow(r,5)*rDOT + 3762*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,1260)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) + 2840*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,2856)*PhiDOT*r*pow(rDOT,5) - 784*pow(rDOT,6))))/(65520.*sqrt(374)*pow(r,4)));}

    else{
        return 0;
    }
}

static COMPLEX16 hl_8_m_4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
     
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_8_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-4*Phi);
    }

}

static COMPLEX16 hl_8_m_min4(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_min4: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,8) * conj(hGO_8_m_4(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,4*Phi);
    }
    
}

//82
static COMPLEX16 hGO_8_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 PhiDOT,UINT4 vpnorder){

    if(vpnorder == 6){
        return (((-1 + 7*Nu - 14*pow(Nu,2) + 7*pow(Nu,3))*(-18236*pow(mass,4) + 10080*pow(r,4)*pow(PhiDOT*r 
        - Complex(0,1)*rDOT,3)*pow(PhiDOT*r + Complex(0,1)*rDOT,5) + 8*pow(mass,3)*r*(2357*pow(PhiDOT,2)*pow(r,2) 
        - Complex(0,16940)*PhiDOT*r*rDOT + 16717*pow(rDOT,2)) + 18*pow(mass,2)*pow(r,2)*(1297*pow(PhiDOT,4)*pow(r,4) 
        + Complex(0,13440)*pow(PhiDOT,3)*pow(r,3)*rDOT - 5826*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,2) 
        + Complex(0,14070)*PhiDOT*r*pow(rDOT,3) - 9188*pow(rDOT,4)) - 45*mass*pow(r,3)*(758*pow(PhiDOT,6)*pow(r,6) 
        + Complex(0,2891)*pow(PhiDOT,5)*pow(r,5)*rDOT + 564*pow(PhiDOT,4)*pow(r,4)*pow(rDOT,2) 
        + Complex(0,5740)*pow(PhiDOT,3)*pow(r,3)*pow(rDOT,3) - 2000*pow(PhiDOT,2)*pow(r,2)*pow(rDOT,4) 
        + Complex(0,2856)*PhiDOT*r*pow(rDOT,5) - 1568*pow(rDOT,6))))/(288288.*sqrt(85)*pow(r,4)));}

    else{
        return 0;
    }
}

static COMPLEX16 hl_8_m_2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){

    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)  * hGO_8_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder) * cpolar(1,-2*Phi);
    }

}

static COMPLEX16 hl_8_m_min2(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 R,UINT4 vpnorder){
    
    if ((vpnorder < 0) || (vpnorder > 7)) {
    XLAL_ERROR(XLAL_EINVAL,"Error in hl_8_m_min2: Input PN order parameter should be between [0, 7].");
    }

    else{
        return ((4*mass*Nu*sqrt(M_PI/5.))/R)   * pow(-1,8) * conj(hGO_8_m_2(mass,Nu,r,rDOT,PhiDOT,vpnorder)) * cpolar(1,2*Phi);
    }
   
}


static COMPLEX16 h05PNGOresult(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 inc,REAL8 euler_beta,REAL8 R,UINT4 vpnorder,REAL8 S1z,REAL8 S2z, REAL8 x){


    //0 PN

    if(vpnorder==0){
        COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,0,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,0,S1z,S2z,x);
        return hcombined;
    }
    
    
    //0.5PN
    
    else if(vpnorder==1){
        COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,1,S1z,S2z,x);
        return hcombined;
    }
	
    //1PN
    else if(vpnorder==2){
        COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x) +
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 2)*hl_4_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, -2)*hl_4_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 4)*hl_4_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, -4)*hl_4_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,2,S1z,S2z,x);
        return hcombined;
    }

    //1.5PN
    else if(vpnorder==3){
        COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x) +
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5, 5)*hl_5_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5, -5)*hl_5_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 3)*hl_4_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4,-3)*hl_4_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 1)*hl_4_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, -1)*hl_4_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,3,S1z,S2z,x);
        return hcombined;
    }

    //2PN
    else if(vpnorder==4){
        COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2,2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                               XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2,-2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                                XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, 1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 2, -1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3,2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3,-2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                               XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, 1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 3, -1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 3)*hl_4_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, -3)*hl_4_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4,2)*hl_4_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4,-2)*hl_4_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, 1)*hl_4_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                             XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4, -1)*hl_4_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5,2)*hl_5_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5,-2)*hl_5_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6,2)*hl_6_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6,-2)*hl_6_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4,4)*hl_4_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 4,-4)*hl_4_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5,4)*hl_5_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 5,-4)*hl_5_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6,4)*hl_6_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6,-4)*hl_6_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6, 6)*hl_6_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2, 6,-6)*hl_6_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,4,S1z,S2z,x);

       return hcombined;
    }

    //2.5PN
    else if(vpnorder==5){
         COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,5)*hl_6_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-5)*hl_6_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,3)*hl_6_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-3)*hl_6_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,1)*hl_6_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-1)*hl_6_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,5)*hl_5_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-5)*hl_5_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,4)*hl_5_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-4)*hl_5_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,1)*hl_5_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-1)*hl_5_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,3)*hl_5_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-3)*hl_5_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,2)*hl_5_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-2)*hl_5_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,7)*hl_7_m_7(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-7)*hl_7_m_min7(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,5)*hl_7_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-5)*hl_7_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,3)*hl_7_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-3)*hl_7_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,1)*hl_7_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-1)*hl_7_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,1)*hl_4_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-1)*hl_4_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,2)*hl_4_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-2)*hl_4_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,3)*hl_4_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-3)*hl_4_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,4)*hl_4_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-4)*hl_4_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,5,S1z,S2z,x);
        return hcombined;
    }
    //3PN
    else if(vpnorder == 6) {
         COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,2)*hl_6_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-2)*hl_6_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,3)*hl_6_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-3)*hl_6_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,1)*hl_6_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-1)*hl_6_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,5)*hl_6_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-5)*hl_6_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,4)*hl_6_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-4)*hl_6_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,6)*hl_6_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-6)*hl_6_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,5)*hl_5_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-5)*hl_5_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,1)*hl_5_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-1)*hl_5_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,3)*hl_5_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-3)*hl_5_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,2)*hl_5_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-2)*hl_5_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,4)*hl_5_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-4)*hl_5_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,6)*hl_7_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-6)*hl_7_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,4)*hl_7_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-4)*hl_7_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,2)*hl_7_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-2)*hl_7_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,1)*hl_4_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-1)*hl_4_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,2)*hl_4_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-2)*hl_4_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,3)*hl_4_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-3)*hl_4_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,4)*hl_4_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-4)*hl_4_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,8)*hl_8_m_8(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,-8)*hl_8_m_min8(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,6)*hl_8_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,-6)*hl_8_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,4)*hl_8_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,-4)*hl_8_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,2)*hl_8_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,8,-2)*hl_8_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,6);
                             
        return hcombined;
    }
    //3.5PN(only quasi-circular spinning terms)
    else{
         COMPLEX16 hcombined = XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,2)*hl_6_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-2)*hl_6_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,4)*hl_6_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-4)*hl_6_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,6)*hl_6_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,6,-6)*hl_6_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,5)*hl_5_m_5(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-5)*hl_5_m_min5(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,1)*hl_5_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-1)*hl_5_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,3)*hl_5_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-3)*hl_5_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,2)*hl_5_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-2)*hl_5_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,4)*hl_5_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,5,-4)*hl_5_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,6)*hl_7_m_6(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-6)*hl_7_m_min6(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,4)*hl_7_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-4)*hl_7_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,2)*hl_7_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,7,-2)*hl_7_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,1)*hl_4_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-1)*hl_4_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,2)*hl_4_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-2)*hl_4_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,3)*hl_4_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-3)*hl_4_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,4)*hl_4_m_4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,4,-4)*hl_4_m_min4(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,3)*hl_3_m_3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-3)*hl_3_m_min3(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,2)*hl_3_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-2)*hl_3_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,1)*hl_3_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,3,-1)*hl_3_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,1)*hl_2_m_1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-1)*hl_2_m_min1(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,2)*hl_2_m_2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x)+
                              XLALSpinWeightedSphericalHarmonic(inc,euler_beta, -2,2,-2)*hl_2_m_min2(mass,Nu,r,rDOT,Phi,PhiDOT,R,7,S1z,S2z,x);
                             
        return hcombined;
    }

}

static REAL8 hplusGO(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 inc,REAL8 euler_beta,REAL8 R,UINT4 vpnorder,REAL8 S1z,REAL8 S2z, REAL8 x){
// nu = symmetric mass ratio m1m2/M^2
            REAL8 hplus = creal(h05PNGOresult(mass,Nu,r,rDOT,Phi,PhiDOT,inc,euler_beta,R,vpnorder,S1z,S2z,x));
            return hplus;
}

static REAL8 hcrossGO(REAL8 mass,REAL8 Nu,REAL8 r,REAL8 rDOT,REAL8 Phi,REAL8 PhiDOT,REAL8 inc,REAL8 euler_beta,REAL8 R,UINT4 vpnorder,REAL8 S1z,REAL8 S2z, REAL8 x){
// nu = symmetric mass ratio m1m2/M^2
            REAL8 hcross = -1.0*cimag(h05PNGOresult(mass,Nu,r,rDOT,Phi,PhiDOT,inc,euler_beta,R,vpnorder,S1z,S2z,x));
            return hcross;
}
// End of this file