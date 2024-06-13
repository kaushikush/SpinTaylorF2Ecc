/*
*  Copyright (C) 2011 Drew Keppel, 2012 Riccardo Sturani
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>
//#include "LALSimInspiralPNCoefficients.c"

#include <math.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Computes the PN Coefficients for using in the PN energy equation.
 *
 * Terms given in equation 3.1 of: Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 * For the spin terms a good reference are (3.15) and (3.16) of 1303.7412
 *
 * In the latest version coefficients of the terms n.S and L.S are reported
 * "Averaged" spin coefficients refer to the ones obtained by orbital averaging,
 * i.e. by using
 * n_i n_j = 1/2 (\f$\delta_{ij} - \hat LN_i \hat LN_j\f$)
 * However such orbital averaging at 2PN would introduce corrections
 * at 3PN, as LNh is not constant.
 */



/* Non-spin phasing terms - see arXiv:0907.0700, Eq. 3.18 */
//static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_2PNCoeff(
//         REAL8 eta
//     )
// {
//         return 5.*(74.3/8.4 + 11.*eta)/9.;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_3PNCoeff(
//         REAL8 UNUSED eta)
// {
//         return -16.*LAL_PI;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_4PNCoeff(
//         REAL8 eta
//     )
// {
//         return 5.*(3058.673/7.056 + 5429./7.*eta+617.*eta*eta)/72.;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNCoeff(
//         REAL8 eta
//     )
// {
//         return 5./9.*(772.9/8.4-13.*eta)*LAL_PI;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNLogCoeff(
//         REAL8 eta
//     )
// {
//         return 5./3.*(772.9/8.4-13.*eta)*LAL_PI;
// }

// static REAL8 UNUSED XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(
//         REAL8 UNUSED eta
//     )
// {
//   return -684.8/2.1;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNCoeff(
//         REAL8 eta
//     )
// {
//   return 11583.231236531/4.694215680 - 640./3.*LAL_PI*LAL_PI - 684.8/2.1*LAL_GAMMA + eta*(-15737.765635/3.048192 + 225.5/1.2*LAL_PI*LAL_PI) + eta*eta*76.055/1.728 - eta*eta*eta*127.825/1.296 + XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(eta)*log(4.);
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_7PNCoeff(
//         REAL8 eta
//     )
// {
//         return LAL_PI*(770.96675/2.54016 + 378.515/1.512*eta - 740.45/7.56*eta*eta);
// }

// /* Spin-orbit terms - can be derived from arXiv:1303.7412, Eq. 3.15-16 */

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(
//         REAL8 mByM
//     )
// {
//         return  mByM*(25.+38./3.*mByM);
// }

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(
        REAL8 mByM,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (4.3*mByM*(3.4401 + 2.6779*mByM))/2.11950;
                break;
        case 3:
                phase = - (mByM*(10.2 + 5.5*mByM))/5.4;
                break;
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_3PNSOEcc03Coeff(
//         REAL8 mByM
//     )
// {
//         return  - (mByM*(102 + 55*mByM))/54.;
// }

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(
        REAL8 mByM
    )
{
  return -mByM*(1391.5/8.4-mByM*(1.-mByM)*10./3.+ mByM*(1276./8.1+mByM*(1.-mByM)*170./9.));

} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(
        REAL8 mByM,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (-7.31*mByM*(-12.42043275 - 12.73890388*mByM - 2.17526764*mByM*mByM + 4.15391200*mByM*mByM*mByM))/1.136594592;
                break;
        case 2:
                phase = (4.3*mByM*(9.7458033 - 11.3891009*mByM + 4.2042952*mByM*mByM + 14.7712964*mByM*mByM*mByM))/2.136456;
                break;
        case 3:
                phase = (7.31*mByM*(-4.1731566 - 9.5824403*mByM + 3.3785668*mByM*mByM + 3.9536420*mByM*mByM*mByM))/4.426736832;
                break;
        case 5:
                phase = -((mByM*(59.0055 + 6.3902*mByM + 27.9356*mByM*mByM + 34.5760*mByM*mByM*mByM))/3.8880);
                break;
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNSOEcc50Coeff(
//         REAL8 mByM
//     )
// {
//   return (-731*mByM*(-1242043275 - 1273890388*mByM - 217526764*mByM*mByM + 415391200*mByM*mByM*mByM))/1.136594592e10;

// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNSOEcc32Coeff(
//         REAL8 mByM
//     )
// {
//   return (43*mByM*(97458033 - 113891009*mByM + 42042952*mByM*mByM + 147712964*mByM*mByM*mByM))/2.136456e8;

// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNSOEcc23Coeff(
//         REAL8 mByM
//     )
// {
//   return (731*mByM*(-41731566 - 95824403*mByM + 33785668*mByM*mByM + 39536420*mByM*mByM*mByM))/4.426736832e9;

// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_5PNSOEcc05Coeff(
//         REAL8 mByM
//     )
// {
//   return -((mByM*(590055 + 63902*mByM + 279356*mByM*mByM + 345760*mByM*mByM*mByM))/38880.);

// }

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(
        REAL8 mByM)
{
  return LAL_PI*mByM*(1490./3. + mByM*260.);
} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(
        REAL8 mByM,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = LAL_PI*(-7.31*mByM*(54.449547 + 37.080578*mByM))/4.883328;
                break;
        case 3:
                phase = LAL_PI*(4.3*mByM*(19.656399 + 13.701538*mByM))/1.52604;
                break;
        case 6:
                phase = -LAL_PI*((mByM*(8.1717 + 7.5326*mByM))/1.5552);
                break;
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSOEcc60Coeff(
//         REAL8 mByM)
// {
//   return LAL_PI*(-731*mByM*(54449547 + 37080578*mByM))/4.883328e8;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSOEcc33Coeff(
//         REAL8 mByM)
// {
//   return LAL_PI*(43*mByM*(19656399 + 13701538*mByM))/1.52604e7;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSOEcc06Coeff(
//         REAL8 mByM)
// {
//   return -LAL_PI*((mByM*(81717 + 75326*mByM))/15552.);
// }


/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(
        REAL8 mByM
    )
{
  REAL8 eta=mByM*(1.-mByM);
  return mByM*(-17097.8035/4.8384+eta*28764.25/6.72+eta*eta*47.35/1.44 + mByM*(-7189.233785/1.524096+eta*458.555/3.024-eta*eta*534.5/7.2));
} */

/*
 * Spin-squared corrections to TF2 phasing
 * Compute 2.0PN SS, QM, and self-spin
 * See Eq. (6.24) in arXiv:0810.5336
 * 9b,c,d in arXiv:astro-ph/0504538
 * Note that these terms are understood to multiply
 * dimensionless spin magnitudes \chi_i=S_i/m_i^2
 * differently from the convention adopted for SpinTaylorTX
 * whose spinning coefficients multiply \chi_LAL=S_i/M^2
 * where M=m_1+m_2.
 * See also https://dcc.ligo.org/T1800298
 */

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNS1S2Coeff(
        REAL8 eta
    )
{
  return  247./4.8*eta;
} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNS1S2EccCoeff(
        REAL8 eta,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (-119.23341*eta)/8.56592;
                break;
        case 4:
                phase = (5.5*eta)/1.6;
                break;
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_4PNS1S2Ecc40Coeff(
//         REAL8 eta
//     )
// {
//   return  (-11923341*eta)/856592.;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_4PNS1S2Ecc04Coeff(
//         REAL8 eta
//     )
// {
//   return  (55*eta)/16.;
// }

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNS1S2OCoeff(
        REAL8 eta
    )
{
  return  -721./4.8*eta;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(
        REAL8 mByM
    )
{
  return -720./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(
        REAL8 mByM
    )
{
  return 1./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(
        REAL8 mByM
    )
{
  return 240./9.6*mByM*mByM;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(
        REAL8 mByM
    )
{
  return -7./9.6*mByM*mByM;
} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(
        REAL8 mByM,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (-38.350453*mByM*mByM)/5.139552;
                break;
        case 4:
                phase = (19.1*mByM*mByM)/9.6;
                break;
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_4PNSelf2SEcc40Coeff(
//         REAL8 mByM
//     )
// {
//   return (-38350453*mByM*mByM)/5.139552e6;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_4PNSelf2SEcc04Coeff(
//         REAL8 mByM
//     )
// {
//   return (191*mByM*mByM)/96.;
// }

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(
        REAL8 eta
    )
{
  return (326.75/1.12L + 557.5/1.8*eta)*eta;
} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(
        REAL8 eta,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (-7.31*eta*(40.20569403 + 28.73178980*eta))/4.10199552;
                break;
        case 2:
                phase = (3.974447*eta*(-28.33 + 55.16*eta))/2.87814912;
                break;
        case 3:
                phase = (-4.3*eta*(11.641317 + 2.945690*eta))/1.14453;
                break;
        case 4:
                phase = (4.0205*eta*(4.09133 + 7.18844*eta))/1.311625728;
                break;
        case 6:
                phase = -((eta*(-62.21289 + 41.53580*eta))/1.86624);
                break;                
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNS1S2OEcc60Coeff(
//         REAL8 eta
//     )
// {
//   return (-731*eta*(4020569403 + 2873178980*eta))/4.10199552e10;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNS1S2OEcc42Coeff(
//         REAL8 eta
//     )
// {
//   return (3974447*eta*(-2833 + 5516*eta))/2.87814912e8;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNS1S2OEcc33Coeff(
//         REAL8 eta
//     )
// {
//   return (-43*eta*(11641317 + 2945690*eta))/1.14453e7;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNS1S2OEcc24Coeff(
//         REAL8 eta
//     )
// {
//   return (40205*eta*(409133 + 718844*eta))/1.311625728e9;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNS1S2OEcc06Coeff(
//         REAL8 eta
//     )
// {
//   return -((eta*(-6221289 + 4153580*eta))/186624.);
// }

/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(
        REAL8 mByM
    )
{
  return (-4108.25/6.72-108.5/1.2*mByM+125.5/3.6*mByM*mByM)*mByM*mByM;
} */

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(
        REAL8 mByM,
        INT4 v0_order
    )
{
        REAL8 phase = 0.0;
        switch (v0_order) {
        case 0: 
                phase = (7.31*mByM*mByM*(-29.99478969 - 50.20858668*mByM + 43.51928140*mByM*mByM))/8.20399104;
                break;
        case 2:
                phase = -((3.8350453*mByM*mByM*(28.33 - 55.16*mByM + 55.16*mByM*mByM))/5.180668416);
                break;
        case 3:
                phase = -((4.3*mByM*mByM*(3.508902 + 4.623513*mByM + 1.472845*mByM*mByM))/1.14453);
                break;
        case 4:
                phase = -((1.39621*mByM*mByM*(-40.9133 - 71.8844*mByM + 71.8844*mByM*mByM))/7.869754368);
                break;
        case 6:
                phase = (mByM*mByM*(25.711461 + 3.330852*mByM + 44.320444*mByM*mByM))/2.612736;
                break;                
        }
        return  phase;
}

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSelf2SEcc60Coeff(
//         REAL8 mByM
//     )
// {
//   return (731*mByM*mByM*(-2999478969 - 5020858668*mByM + 4351928140*mByM*mByM))/8.20399104e10;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSelf2SEcc42Coeff(
//         REAL8 mByM
//     )
// {
//   return -((38350453*mByM*mByM*(2833 - 5516*mByM + 5516*mByM*mByM))/5.180668416e9);
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSelf2SEcc33Coeff(
//         REAL8 mByM
//     )
// {
//   return -((43*mByM*mByM*(3508902 + 4623513*mByM + 1472845*mByM*mByM))/1.14453e7;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSelf2SEcc24Coeff(
//         REAL8 mByM
//     )
// {
//   return -((139621*mByM*mByM*(-409133 - 718844*m1ByM + 718844*mByM*mByM))/7.869754368e9);
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_6PNSelf2SEcc06Coeff(
//         REAL8 mByM
//     )
// {
//   return (mByM*mByM*(25711461 + 3330852*mByM + 44320444*mByM*mByM))/2.612736e6;
// }


/* static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(
        REAL8 mByM
    )
{
  return (4703.5/8.4+2935./6.*mByM-120.*mByM*mByM)*mByM*mByM;
} */

/*
 * Tidal corrections to F2 phasing
 * See arXiv:1101.1673
 */

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(
//         REAL8 mByM /**< ratio of object mass to total mass */
//     )
// {
//   return (-288. + 264.*mByM)*mByM*mByM*mByM*mByM;

// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(
//         REAL8 mByM /**< ratio of object mass to total mass */
//     )
// {
//   return (-15895./28. + 4595./28.*mByM + 5715./14.*mByM*mByM - 325./7.*mByM*mByM*mByM)*mByM*mByM*mByM*mByM;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(
//       REAL8 mByM /**< ratio of object mass to total mass */
//     ) 
// /*  literature: Agathos et al (arxiv 1503.0545) eq (5)
//  * the coefficient mByM4 conversion & transformation (6.5PN, 7PN, 7.5PN):
//  * mByM=mA/M: mA= mass star A, M is total mass (mA+mB)
//  * Lambda (unitless) = lambda(m) / mA^5 
//  * to call the function: 
//  * Lambda * XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff 
//  * lambda(m)*mByM^4/mA^5= lambda(m)*(mA/M)^4/(mA)^5= lambda/(M^4*mA) 
//  * =lambda/(mByM*M^5) eq (5) 
//  */
// {
//   return mByM*mByM*mByM*mByM * 24.L*(12.L - 11.L*mByM)*LAL_PI;
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(
//       REAL8 mByM /**< ratio of object mass to total mass */
//     )
// /* literature: Agathos et al (arxiv 1503.0545) eq (5)
//  * caveat: these are incomplete terms
//  * conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above
//  * --> completed by the terms given in equation (4) of :
//  * Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
//  * "Gravitational-wave constraints on the GWTC-2 events by measuring
//  * the tidal deformability and the spin-induced quadrupole moment",
//  * Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
//  */
// {
//   REAL8 mByM3 = mByM*mByM*mByM;
//   REAL8 mByM4 = mByM3 * mByM;
//   return - mByM4 * 5.L*(193986935.L/571536.L - 14415613.L/381024.L*mByM - 57859.L/378.L*mByM*mByM - 209495.L/1512.L*mByM3 + 965.L/54.L*mByM4 - 4.L*mByM4*mByM);
// }

// static REAL8 UNUSED
// XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(
//       REAL8 mByM /**< ratio of object mass to total mass */
//     )
// /* literature: Agathos et al (arxiv 1503.0545) eq (5)
//  * conversion see XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff above 
//  * --> corrected by the terms given in equation (4) of :
//  * Tatsuya Narikawa, Nami Uchikata, Takahiro Tanaka,
//  * "Gravitational-wave constraints on the GWTC-2 events by measuring
//  * the tidal deformability and the spin-induced quadrupole moment",
//  * Phys. Rev. D 104, 084056 (2021), arXiv:2106.09193
//  */
// {
//   REAL8 mByM2 = mByM*mByM;
//   REAL8 mByM3 = mByM2*mByM;
//   REAL8 mByM4 = mByM3*mByM;
//   return mByM4 * 1.L/28.L*LAL_PI*(27719.L - 22415.L*mByM + 7598.L*mByM2 - 10520.L*mByM3) ;
// }

/* The phasing function for TaylorF2 frequency-domain waveform.
 * This function is tested in ../test/PNCoefficients.c for consistency
 * with the energy and flux in this file.
 */
// static void UNUSED
// XLALSimInspiralPNPhasing_F2(
// 	PNPhasingSeries *pfa, /**< \todo UNDOCUMENTED */
// 	const REAL8 m1, /**< Mass of body 1, in Msol */
// 	const REAL8 m2, /**< Mass of body 2, in Msol */
// 	const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
// 	const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
// 	const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
// 	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
// 	const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
// 	LALDict *p /**< LAL dictionary containing accessory parameters */
// 	)
// {
//     const REAL8 mtot = m1 + m2;
//     const REAL8 eta = m1*m2/mtot/mtot;
//     const REAL8 m1M = m1/mtot;
//     const REAL8 m2M = m2/mtot;

//     const REAL8 pfaN = 3.L/(128.L * eta);

//     memset(pfa, 0, sizeof(PNPhasingSeries));

//     pfa->v[0] = 1.L;
//     pfa->v[1] = 0.L;
//     pfa->v[2] = XLALSimInspiralTaylorF2Phasing_2PNCoeff(eta);
//     pfa->v[3] = XLALSimInspiralTaylorF2Phasing_3PNCoeff(eta);
//     pfa->v[4] = XLALSimInspiralTaylorF2Phasing_4PNCoeff(eta);
//     pfa->v[5] = XLALSimInspiralTaylorF2Phasing_5PNCoeff(eta);
//     pfa->vlogv[5] = XLALSimInspiralTaylorF2Phasing_5PNLogCoeff(eta);
//     pfa->v[6] =  XLALSimInspiralTaylorF2Phasing_6PNCoeff(eta);
//     pfa->vlogv[6] = XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(eta);
//     pfa->v[7] = XLALSimInspiralTaylorF2Phasing_7PNCoeff(eta);

//     /* modify the PN coefficients if a non null LALSimInspiralTestGRParam structure is passed */
//     /* BEWARE: this is for the non-spinning case only!*/
//     pfa->v[0]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi0(p));
//     pfa->v[1] = XLALSimInspiralWaveformParamsLookupNonGRDChi1(p);
//     pfa->v[2]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi2(p));
//     pfa->v[3]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi3(p));
//     pfa->v[4]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi4(p));
//     pfa->v[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5(p));
//     pfa->vlogv[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5L(p));
//     pfa->v[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6(p));
//     pfa->vlogv[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6L(p));
//     pfa->v[7]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi7(p));

//     const REAL8 qm_def1=1.+XLALSimInspiralWaveformParamsLookupdQuadMon1(p);
//     const REAL8 qm_def2=1.+XLALSimInspiralWaveformParamsLookupdQuadMon2(p);

//     switch( XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) )
//     {
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
// 	    pfa->v[7] += XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m2M)*chi2L;
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
//             pfa->v[6] += XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m1M)*chi1L
// 	      + XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m2M)*chi2L
// 	      + XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(eta)*chi1L*chi2L
// 	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m1M))*chi1sq
// 	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m2M))*chi2sq;
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
//             pfa->v[5] += XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
// 	      + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L;
//             pfa->vlogv[5] += 3.*(XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
// 				 + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L);
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
// 	    /* 2PN SS, QM, and self-spin */
//             pfa->v[4] += XLALSimInspiralTaylorF2Phasing_4PNS1S2Coeff(eta)*chi1dotchi2+XLALSimInspiralTaylorF2Phasing_4PNS1S2OCoeff(eta)*chi1L*chi2L
// 	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m1M))*chi1L*chi1L
// 	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m2M))*chi2L*chi2L
// 	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m1M))*chi1sq
// 	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m2M))*chi2sq;
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
//             pfa->v[3] += XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m2M)*chi2L;
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
//         case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
//             break;
//         default:
//             XLALPrintError("XLAL Error - %s: Invalid spin PN order %i\n",
// 			   __func__, XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) );
//             XLAL_ERROR_VOID(XLAL_EINVAL);
//             break;
//     }

//     REAL8 lambda1=XLALSimInspiralWaveformParamsLookupTidalLambda1(p);
//     REAL8 lambda2=XLALSimInspiralWaveformParamsLookupTidalLambda2(p);
//     switch( XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) )
//     {
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
//             pfa->v[15] = (lambda1*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m2M));
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT:
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
//             pfa->v[14] = (lambda1*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m2M));
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
//             pfa->v[13] = (lambda1*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m2M));
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
//             pfa->v[12] = (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2M) );
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
//             pfa->v[10] = ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2M) );
// #if __GNUC__ >= 7 && !defined __INTEL_COMPILER
//             __attribute__ ((fallthrough));
// #endif
//         case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
//             break;
//         default:
//             XLALPrintError("XLAL Error - %s: Invalid tidal PN order %i\n",
//                            __func__, XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) );
//             XLAL_ERROR_VOID(XLAL_EINVAL);
//     }


//     /* At the very end, multiply everything in the series by pfaN */
//     for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
//     {
//         pfa->v[ii] *= pfaN;
//         pfa->vlogv[ii] *= pfaN;
//         pfa->vlogvsq[ii] *= pfaN;
//     }
// }




/**
 * Computes the PN Coefficients for using in the TaylorF2Ecc equation.
 *
 * 3-dimensional REAL8 array eccPNCoeffs[ORDER][v_power][v0_power] are calculated,
 * where ORDER is relative PN order, v_power is power of v, and v0_power is power of v0.
 * Note that ORDER = v_power + v0_power.
 *
 * Terms given in equation 6.26 of: Blake Moore, Marc Favata,
 * K.G.Arun, and Chandra Kant Mishra, "Gravitational-wave phasing
 * for low-eccentricity inspiralling compact binaries to 3PN order",
 * Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304
 */

/* static INT4 UNUSED
eccentricityPNCoeffs_F2(REAL8 eta, REAL8 eccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1])
{
  INT4 ret = 0;
  memset(eccPNCoeffs, 0x00, (LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*sizeof(REAL8));
  eccPNCoeffs[0][0][0] = 1.0; // lowest order constant term

  eccPNCoeffs[2][2][0] = 29.9076223/8.1976608 + 18.766963/2.927736*eta; //v^2 term
  eccPNCoeffs[2][1][1] = 0.0; //v*v0 term
  eccPNCoeffs[2][0][2] = 2.833/1.008 - 19.7/3.6*eta; //v0^2 term

  eccPNCoeffs[3][3][0] = -28.19123/2.82600*LAL_PI; //v^3 term
  eccPNCoeffs[3][0][3] = 37.7/7.2*LAL_PI; //v0^3 term

  eccPNCoeffs[4][4][0] = 16.237683263/3.330429696 + 241.33060753/9.71375328*eta+156.2608261/6.9383952*eta*eta; //v^4 term
  eccPNCoeffs[4][2][2] = 84.7282939759/8.2632420864-7.18901219/3.68894736*eta-36.97091711/1.05398496*eta*eta; //v^2*v0^2 term
  eccPNCoeffs[4][0][4] = -1.193251/3.048192 - 66.317/9.072*eta +18.155/1.296*eta*eta;  //v0^4 term

  eccPNCoeffs[5][5][0] = -28.31492681/1.18395270*LAL_PI - 115.52066831/2.70617760*LAL_PI*eta; //v^5 term
  eccPNCoeffs[5][3][2] = -79.86575459/2.84860800*LAL_PI + 55.5367231/1.0173600*LAL_PI*eta; //v^3*v0^2 term
  eccPNCoeffs[5][2][3] = 112.751736071/5.902315776*LAL_PI + 70.75145051/2.10796992*LAL_PI*eta; //v^2*v0^3 term
  eccPNCoeffs[5][0][5] = 76.4881/9.0720*LAL_PI - 94.9457/2.2680*LAL_PI*eta;  //v0^5 term

  eccPNCoeffs[6][6][0] = -436.03153867072577087/1.32658535116800000 + 53.6803271/1.9782000*LAL_GAMMA + 157.22503703/3.25555200*LAL_PI*LAL_PI
                         +(2991.72861614477/6.89135247360 - 15.075413/1.446912*LAL_PI*LAL_PI)*eta
                         +345.5209264991/4.1019955200*eta*eta + 506.12671711/8.78999040*eta*eta*eta
                         + 384.3505163/5.9346000*log(2.0) - 112.1397129/1.7584000*log(3.0); //v^6 term except log(16*v^2) term
  eccPNCoeffs[6][4][2] = 46.001356684079/3.357073133568 + 253.471410141755/5.874877983744*eta
                         - 169.3852244423/2.3313007872*eta*eta - 307.833827417/2.497822272*eta*eta*eta; //v^4*v0^2 term
  eccPNCoeffs[6][3][3] = -106.2809371/2.0347200*LAL_PI*LAL_PI; //v^3*v0^3 term
  eccPNCoeffs[6][2][4] = -3.56873002170973/2.49880440692736 - 260.399751935005/8.924301453312*eta
                         + 15.0484695827/3.5413894656*eta*eta + 340.714213265/3.794345856*eta*eta*eta; //v^2*v0^4 term
  eccPNCoeffs[6][0][6] = 265.31900578691/1.68991764480 - 33.17/1.26*LAL_GAMMA + 12.2833/1.0368*LAL_PI*LAL_PI
                         + (91.55185261/5.48674560 - 3.977/1.152*LAL_PI*LAL_PI)*eta - 5.732473/1.306368*eta*eta
                         - 30.90307/1.39968*eta*eta*eta + 87.419/1.890*log(2.0) - 260.01/5.60*log(3.0);  //v0^6 term except log(16*v0^2) term
  //printPNCoeffs_F2(eccPNCoeffs);
  return ret;
} */

// static INT4 UNUSED
// SpinEccentricityPNCoeffs_F2(
//         const REAL8 m1, 
//         const REAL8 m2, 
//         const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
// 	const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
// 	const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
// 	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
//         REAL8 SpinEccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]
//         )
// {
//   const REAL8 mtot = m1 + m2;
//   const REAL8 eta = m1*m2/mtot/mtot;
//   const REAL8 m1M = m1/mtot;
//   const REAL8 m2M = m2/mtot;
//   //const REAL8 chi1sq = chi1L*chi1L;
//   //const REAL8 chi2sq = chi2L*chi2L; 
  
//   INT4 ret = 0;
//   memset(SpinEccPNCoeffs, 0x00, (LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*sizeof(REAL8));
//   SpinEccPNCoeffs[0][0][0] = 0.0; // lowest order constant term

//   SpinEccPNCoeffs[2][2][0] = 0.0; //v^2 term
//   SpinEccPNCoeffs[2][1][1] = 0.0; //v*v0 term
//   SpinEccPNCoeffs[2][0][2] = 0.0; //v0^2 term

//   SpinEccPNCoeffs[3][3][0] = XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m2M, 0)*chi2L; //v^3 term
//   SpinEccPNCoeffs[3][0][3] = XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m2M, 3)*chi2L; //v0^3 term

//   SpinEccPNCoeffs[4][4][0] = XLALSimInspiralTaylorF2Phasing_4PNS1S2EccCoeff(eta, 0)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m1M, 0)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m2M, 0)*chi2sq; //v^4 term
//   SpinEccPNCoeffs[4][2][2] = 0.0; //v^2*v0^2 term
//   SpinEccPNCoeffs[4][0][4] = XLALSimInspiralTaylorF2Phasing_4PNS1S2EccCoeff(eta, 4)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m1M, 4)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m2M, 4)*chi2sq;  //v0^4 term

//   SpinEccPNCoeffs[5][5][0] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 0)*chi2L; //v^5 term
//   SpinEccPNCoeffs[5][3][2] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 2)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 2)*chi2L; //v^3*v0^2 term
//   SpinEccPNCoeffs[5][2][3] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 3)*chi2L; //v^2*v0^3 term
//   SpinEccPNCoeffs[5][0][5] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 5)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 5)*chi2L;  //v0^5 term

//   SpinEccPNCoeffs[6][6][0] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 0)*chi2L 
//                          + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 0)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 0)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 0)*chi2sq; //v^6 term except log(16*v^2) term
//   SpinEccPNCoeffs[6][4][2] = XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 2) + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 2)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 2)*chi2sq; //v^4*v0^2 term
//   SpinEccPNCoeffs[6][3][3] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 3)*chi2L 
//                          + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 3)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 3)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 3)*chi2sq; //v^3*v0^3 term
//   SpinEccPNCoeffs[6][2][4] = XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 4)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 4)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 4)*chi2sq; //v^2*v0^4 term
//   SpinEccPNCoeffs[6][0][6] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 6)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 6)*chi2L 
//                          + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 6)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 6)*chi1sq
//                          + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 6)*chi2sq;  //v0^6 term except log(16*v0^2) term
//   //printSpinPNCoeffs_F2(SpinEccPNCoeffs);
//   return ret;
// }


/* Non-spin phasing terms obtained from Black-hole perturbation (BHP) theory */
static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_8PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        switch (log_order) {
        case 0: 
                phase =  2554404624135128353/830425530654720 - (36812*LAL_GAMMA)/189. - (90490*LAL_PI*LAL_PI)/567. - (1011020*log(2.0))/3969. - (26325*log(3.0))/196.;
                break;
        case 1:
                phase = -2554404624135128353/276808510218240 + (36812*LAL_GAMMA)/63. + (90490*LAL_PI*LAL_PI)/189. + (1011020*log(2.0))/1323. + (78975*log(3.0))/196.;
                break;
        case 2:
                phase =  18406/63;
                break;

        }
        return  phase;
}


static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_9PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        switch (log_order) {
        case 0: 
                phase =  -((640*LAL_PI*LAL_PI*LAL_PI)/3.) + LAL_PI*(105344279473163/18776862720 - (13696*LAL_GAMMA)/21. - (27392*log(2.0))/21.);
                break;
        case 1:
                phase =  -13696*LAL_PI/21;
                break;


        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_10PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        switch (log_order) {
        case 0: 
                phase =  -11385.142032070502 + (6470582647*LAL_GAMMA)/2.750517e7 + (578223115*LAL_PI*LAL_PI)/3.048192e6 + (53992839431*log(2.0))/5.501034e7 - (5512455*log(3.0))/21952.;
                break;
        case 1:
                phase =  6470582647/27505170;
                break;


        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_11PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        switch (log_order) {
        case 0: 
                phase =  -((94390*LAL_PI*LAL_PI*LAL_PI)/567.) + LAL_PI*(1862462447418252011/276808510218240 - (3558011*LAL_GAMMA)/7938. - (862549*log(2.0))/1134. - (26325*log(3.0))/196.);
                break;
        case 1:
                phase =  (-3558011*LAL_PI)/7938.;
                break;


        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_12PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        //REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        /* 
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */

        switch (log_order) {
        case 0: 
                phase =  -2533276330458753684746142171281/350982066886345200815308800 + (2930944*LAL_GAMMA*LAL_GAMMA)/15435. + (1024*LAL_PI4)/21. - (2323495644411647809*log(2.0))/3.700391542848e14 + 
                        (11723776*log(2.0)*log(2.0))/15435. + LAL_PI*LAL_PI*(-1934.9592520858043 + (219136*log(2.0))/441.) + 
                        LAL_GAMMA*(-1031803748014155817/370039154284800 + (109568*LAL_PI*LAL_PI)/441. + (11723776*log(2.0))/15435.) + (4043034784143*log(3.0))/1.08287635456e11 + 
                         (188720703125*log(5.0))/9.13296384e8;
                break;
        case 1:
                phase =  -1031803748014155817/370039154284800 + (5861888*LAL_GAMMA)/15435. + (109568*LAL_PI*LAL_PI)/441. + (11723776*log(2.0))/15435.;
                break;
        case 2:
                phase = -3558011*LAL_PI/7938.;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_13PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        /* REAL8 LAL_PI4 = LAL_PI3*LAL_PI;
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */

        switch (log_order) {
        case 0: 
                phase =  -((255852649*LAL_PI3)/3.048192e6) + LAL_PI*(110116961615984436326569/9094082255703244800 - (214170945619*LAL_GAMMA)/8.8016544e8 - (78840091699*log(2.0))/1.76033088e8 - (6144255*log(3.0))/175616.);
                break;
        case 1:
                phase =  -(214170945619*LAL_PI)/8.8016544e8;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_14PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        //REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        /* REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  543321785858790745520247781639425619/280883928487804337308475326464000 + (365125274*LAL_GAMMA*LAL_GAMMA)/1.250235e6 + (142804*LAL_PI4)/1701. - (1885349108439532357613*log(2.0))/1.252450380413232e17 + 
                (564573274*log(2.0)*log(2.0))/1.250235e6 - (39912634025297167*log(3.0))/4.223217782784e13 + 13455/14*log(2.0)*log(3.0) - (342225*log(3.0)*log(3.0))/1372. + 
                 LAL_PI*LAL_PI*(-184253509229260616617/29895319103569920 + (30012272*log(2.0))/35721. + (2925*log(3.0))/98.) + 
                 LAL_GAMMA*(-121692529309026804633437/14027444260628198400 + (741584*LAL_PI*LAL_PI)/1701. + (1168570484*log(2.0))/1.250235e6 + (158535*log(3.0))/686.) - (1039284912109375*log(5.0))/1.282268123136e12 + 547772*Zeta(3)/567.;
                break;
        case 1:
                phase =  -121353693547838861414237/14027444260628198400 + (730250548*LAL_GAMMA)/1.250235e6 + (741584*LAL_PI*LAL_PI)/1701. + (1168570484*log(2.0))/1.250235e6 + (158535*log(3.0))/686.;
                break;
        case 2:
                phase = 365125274/1.250235e6;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_15PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI3*LAL_PI2;
        /* REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  -8192*LAL_PI5/315. + LAL_PI3*(17430874897721/11496038400 - (438272*LAL_GAMMA)/2205. - (876544*log(2.0))/2205.) + 
                 LAL_PI*(433563913308936569405050851311/30680250601953251819520000 - (23447552*LAL_GAMMA*LAL_GAMMA)/77175. + LAL_GAMMA*(116629249257534007/25876863936000 - (93790208*log(2.0))/77175.) + 
                  (2392348417622159629*log(2.0))/2.84645503296e14 - (93790208*log(2.0)*log(2.0))/77175. + (69120873720189*log(3.0))/2.7071908864e11 + 
                  (37744140625*log(5.0))/2.28324096e8 - (438272*Zeta(3))/735.);
                break;
        case 1:
                phase =  (-438272*LAL_PI3)/2205. + LAL_PI*(116629249257534007/25876863936000 - (46895104*LAL_GAMMA)/77175. - (93790208*log(2.0))/77175.);
                break;
        case 2:
                phase = -23447552*LAL_PI/77175.;
                break;

        }
        return  phase;
}


static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_16PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        //REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        /* REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  7717784949124343229586732254526318153965031/258844651922937257785903118448328704000 + (362659047537859*LAL_GAMMA*LAL_GAMMA)/1.67737528728e12 + (26996009*LAL_PI4)/762048. - 
                 (4097431198874706421969490123*log(2.0))/1.063933734498952e23 + (2565322132571243*log(2.0)*log(2.0))/1.67737528728e12 + 
                 LAL_GAMMA*(-462044918046394592295535811/21278674689979041054720 + (3472481167919*LAL_PI2)/1.1618183808e10 + (876627606865651*log(2.0))/8.3868764364e11 - (137190105*log(3.0))/1.690304e6) - 
                (484717639256500532148913*log(3.0))/7.076085859410248e19 - (2757917565*log(2.0)*log(3.0))/1.690304e6 + (19826120925*log(3.0)*log(3.0))/2.7044864e7 + 
                LAL_PI2*(-43945722625778811724177/2858140137506734080 + (3219899214725*log(2.0))/1.1618183808e10 + (149236425*log(3.0))/965888.) + (109381686589111328125*log(5.0))/9.207710938614989e16 + 
                (2025852318599963*log(7.0))/6.4876660992e12 + (310622273305*Zeta(3))/2.42045496e8;
                break;
        case 1:
                phase =  -458491559217666605970769091/21278674689979041054720 + (362659047537859*LAL_GAMMA)/8.3868764364e11 + (3472481167919*LAL_PI2)/1.1618183808e10 + (876627606865651*log(2.0))/8.3868764364e11 - 
                 (137190105*log(3.0))/1.690304e6;
                break;
        case 2:
                phase = 362659047537859/1.67737528728e12;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_17PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI3*LAL_PI2;
        /* REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  -((272528*LAL_PI5)/5103.) + LAL_PI3*(24337991480444006705/3986042547142656 - (15870500*LAL_GAMMA)/35721. - (31644868*log(2.0))/35721.) + 
                LAL_PI*(300657703907558233213139082820138069/86425824150093642248761638912000 - (1573299001*LAL_GAMMA*LAL_GAMMA)/2.50047e6 - (639497071*log(2.0)*log(2.0))/357210. + log(2.0)*(2032539221128816308113111/56109777042512793600 - (13455*log(3.0))/14.) + 
                 LAL_GAMMA*(25811167832239726073989/1438712231859302400 - (2848666241*log(2.0))/1.250235e6 - (158535*log(3.0))/686.) + (9786240502440407*log(3.0))/5.11905185792e12 + 
                 (342225*log(3.0)*log(3.0))/1372. - (624099365234375*log(5.0))/5.69896943616e11 - (19526807*Zeta(3))/11907.);
                break;
        case 1:
                phase =  -15870500*LAL_PI3/35721. + LAL_PI*(25672158289188262189189/1438712231859302400 - (1573299001*LAL_GAMMA)/1.250235e6 - (2848666241*log(2.0))/1.250235e6 - (158535*log(3.0))/686.);
                break;
        case 2:
                phase = -1573299001*LAL_PI/2.50047e6;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_18PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI3*LAL_PI2;
        REAL8 LAL_PI6 = LAL_PI4*LAL_PI2;
        /* REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  2611287968137573440061272328756049901318885535250051/28200719446330134924944802369483049009152000000 + (10035552256*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/4.5147375e7 + (65536*LAL_PI6)/4095. - (61436010527654028941812522108319616851*log(2.0))/8.232407770952117e32 - 
                (1565361753466011604825507*log(2.0)*log(2.0))/7.738906362923735e19 + (80284418048*log(2.0)*log(2.0)*log(2.0))/4.5147375e7 + 
                LAL_PI4*(-18696174900478561/14874795936000 + (7012352*log(2.0))/20475.) + LAL_GAMMA*LAL_GAMMA*
                (-32762034670645379398153/7035369420839760000 + (187580416*LAL_PI2)/429975. + (20071104512*log(2.0))/1.5049125e7) - (1415161243067749778175775050943*log(3.0))/9.914065669395765e26 + 
                (2684948834816170393*log(2.0)*log(3.0))/6.290834822272e14 - (9441010581979717287*log(3.0)*log(3.0))/5.0326678578176e15 + 
                (320820624461015801123046875*log(5.0))/6.110817064654039e22 + (826906181640625*log(2.0)*log(5.0))/6.36681741696e11 - (145881103515625*log(5.0)*log(5.0))/1.81909069056e11 - 
                (10614597927041548993*log(7.0))/3.910293294336e15 - (1902991615901777*Zeta(3))/2.0331821664e11 + (375160832*log(2.0)*Zeta(3))/143325. + 
                LAL_PI2*(-181992692642605073615070084344909/7521044290421682874613760000 - (29779322706900420631*log(2.0))/2.577058395912e15 + (750321664*log(2.0)*log(2.0))/429975. - 
                (632421821599629*log(3.0))/2.01105608704e12 - (264208984375*log(5.0))/6.36045696e8 + (3506176*Zeta(3))/4095.) + 
                LAL_GAMMA*(-1857079904275735917321329684722559231/56130052983764430108120883200000 + (3506176*LAL_PI4)/20475. - (703733677389283431063491*log(2.0))/3.869453181461868e19 + 
                (40142209024*log(2.0)*log(2.0))/1.5049125e7 + LAL_PI2*(-2070280696258570241/322132299489000 + (750321664*log(2.0))/429975.) - (44111547682284331*log(3.0))/1.2581669644544e14 - 
                (32376923828125*log(5.0))/2.12227247232e11 + (187580416*Zeta(3))/143325.) - (1753088*Zeta(5))/1365.;
                break;
        case 1:
                phase =  -1812020988301897301122044543116159231/56130052983764430108120883200000 + (10035552256*LAL_GAMMA*LAL_GAMMA)/1.5049125e7 + (3506176*LAL_PI4)/20475. - 
                (703733677389283431063491*log(2.0))/3.869453181461868e19 + (40142209024*log(2.0)*log(2.0))/1.5049125e7 + 
                LAL_PI2*(-2070280696258570241/322132299489000 + (375160832*LAL_GAMMA)/429975. + (750321664*log(2.0))/429975.) + 
                LAL_GAMMA*(-32762034670645379398153/3517684710419880000 + (40142209024*log(2.0))/1.5049125e7) - (44111547682284331*log(3.0))/1.2581669644544e14 - (32376923828125*log(5.0))/2.12227247232e11 + 
                (187580416*Zeta(3))/143325.;
                break;
        case 2:
                phase = -32762034670645379398153/7035369420839760000 + (10035552256*LAL_GAMMA)/1.5049125e7 + (187580416*LAL_PI2)/429975. + (20071104512*log(2.0))/1.5049125e7;
                break;
        case 3:
                phase = 10035552256/4.5147375e7;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_19PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI3*LAL_PI;
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        /* REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  -((405069545*LAL_PI5)/1.1002068e7) + LAL_PI3*(3745295331618215890392197/210073300106744954880 - (227231932601*LAL_GAMMA)/7.26136488e8 - (3152920672801*log(2.0))/5.082955416e9 + 
                (5923125*log(3.0))/845152.) + LAL_PI*(-2972949973483323362582395485846960475985281/56622267608142525140666307160571904000 - (695520298806979*LAL_GAMMA*LAL_GAMMA)/1.46770337637e12 - 
                (468133310198633*log(2.0)*log(2.0))/2.93540675274e11 + log(2.0)*(623075439689343589150581375079/6702782527343397932236800 - (143739765*log(3.0))/422576.) + 
                LAL_GAMMA*(346342146222216593179118528671/6702782527343397932236800 - (1315954801615778*log(2.0))/7.33851688185e11 - (214101225*log(3.0))/2.958032e6) + 
                (24373239436490584126577*log(3.0))/1.5478937817459917e19 + (1198129725*log(3.0)*log(3.0))/1.1832128e7 + (495218832177490234375*log(5.0))/1.611349414257623e17 + 
                (289407474085709*log(7.0))/2.703194208e11 - (250153347118*Zeta(3))/2.11789809e8);
                break;
        case 1:
                phase =  -227231932601*LAL_PI3/7.26136488e8 + LAL_PI*(341399109302173248955871783071/6702782527343397932236800 - (695520298806979*LAL_GAMMA)/7.33851688185e11 - (1315954801615778*log(2.0))/7.33851688185e11 - 
                (214101225*log(3.0))/2.958032e6);
                break;
        case 2:
                phase = -695520298806979*LAL_PI/1.46770337637e12;
                break;

        }
        return  phase;
}


static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_20PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI3*LAL_PI2;
        REAL8 LAL_PI6 = LAL_PI4*LAL_PI2;
        /* REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  810993819201019083909129288789745966878050930684209702909867/3529014568515173737835335014441537343108505665536000000 + (326800201792*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/6.56373375e8 + (2438944*LAL_PI6)/59535. + 
                (11302453352831652859762918435368842794590941*log(2.0))/5.676014650653898e38 - (1935901659834424685929931057*log(2.0)*log(2.0))/3.191243569292915e22 + 
                (266845766624*log(2.0)*log(2.0)*log(2.0))/6.56373375e8 - (124898325764776285597768776151510097*log(3.0))/4.24791411310404e31 - 
                (51864882191903105990179*log(2.0)*log(3.0))/1.585290375212544e18 + 6189.3*log(2.0)*log(2.0)*log(3.0) + (2191794881449965617607*log(3.0)*log(3.0))/4.02613428625408e17 - 
                3212.7244897959185*log(2.0)*log(3.0)*log(3.0) + (2669355*log(3.0)*log(3.0)*log(3.0))/4802. + LAL_PI4*(-21861474635316964260311/3986042547142656000 + (46632128*log(2.0))/59535. + (1287*log(3.0))/49.) + 
                LAL_GAMMA*LAL_GAMMA*(-22297093247166523050384303829/1053110377866661994880000 + (6589493084*LAL_PI2)/6.251175e6 + (573052573216*log(2.0))/2.18791125e8 + (8592597*log(3.0))/24010.) - 
                (73242902594746089113509268123046875*log(5.0))/1.6010301600164368e30 - (2213232984150390625*log(2.0)*log(5.0))/2.75046512412672e14 + 
                (188720703125*log(3.0)*log(5.0))/9.7140736e7 + (803367237060546875*log(5.0)*log(5.0))/1.83364341608448e14 + (839337151526409339507875101*log(7.0))/8.781380061671346e22 - 
                (618424409913772295214347*Zeta(3))/1.1689536883856833e19 + (3046519936*log(2.0)*Zeta(3))/297675. - 92313/49*log(3.0)*Zeta(3) + 
                LAL_PI2*(-20554513564716811955601868724904110537/1021396103592015772030819368960000 + (24367229788*log(2.0)*log(2.0))/6.251175e6 - (3246890378020543273*log(3.0))/1.6892871131136e15 - (68445*log(3.0)*log(3.0))/686. + 
                log(2.0)*(-33416935340859630553225597/561097770425127936000 + (2691*log(3.0))/7.) + (5223826806640625*log(5.0))/2.564536246272e12 + (44946392*Zeta(3))/19845.) + 
                LAL_GAMMA*(-2229116754653630564215735452246949072472839/113520293013077952189485209972899840000 + (13448704*LAL_PI4)/33075. + (815664712192*log(2.0)*log(2.0))/2.18791125e8 - 
                (6957474129994211619751*log(3.0))/1.585290375212544e18 - (3709719*log(3.0)*log(3.0))/4802. + 
                LAL_PI2*(-1257293404704810539158171/43161366955779072000 + (25678092728*log(2.0))/6.251175e6 + (31707*log(3.0))/343.) + log(2.0)*(-8762245239397225314045688253/105311037786666199488000 + (729261*log(3.0))/245.) + 
                (2566285304345703125*log(5.0))/1.925325586888704e15 + (8713301312*Zeta(3))/2.083725e6) - (2320448*Zeta(5))/245.;
                break;
        case 1:
                phase =  -1973472053151920484336300004251512620714759/113520293013077952189485209972899840000 + (326800201792*LAL_GAMMA*LAL_GAMMA)/2.18791125e8 + (13448704*LAL_PI4)/33075. + (815664712192*log(2.0)*log(2.0))/2.18791125e8 - 
                (6957474129994211619751*log(3.0))/1.585290375212544e18 - (3709719*log(3.0)*log(3.0))/4802. + 
                LAL_PI2*(-1243948488571870006217371/43161366955779072000 + (13178986168*LAL_GAMMA)/6.251175e6 + (25678092728*log(2.0))/6.251175e6 + (31707*log(3.0))/343.) + 
                LAL_GAMMA*(-22172664623600041789656803029/526555188933330997440000 + (1146105146432*log(2.0))/2.18791125e8 + (8592597*log(3.0))/12005.) + log(2.0)*(-8724916652327280935827438013/105311037786666199488000 + (729261*log(3.0))/245.) + 
                (2566285304345703125*log(5.0))/1.925325586888704e15 + (8713301312*Zeta(3))/2.083725e6;
                break;
        case 2:
                phase = -22172664623600041789656803029/1053110377866661994880000 + (326800201792*LAL_GAMMA)/2.18791125e8 + (6589493084*LAL_PI2)/6.251175e6 + (573052573216*log(2.0))/2.18791125e8 + (8592597*log(3.0))/24010.;
                break;
        case 3:
                phase = 326800201792/6.56373375e8;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_21PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI3*LAL_PI2;
        REAL8 LAL_PI6 = LAL_PI4*LAL_PI2;
        REAL8 LAL_PI7 = LAL_PI4*LAL_PI3;
        //REAL8 LAL_PI8 = LAL_PI7*LAL_PI;
        
        switch (log_order) {
        case 0: 
                phase =  -((8192*LAL_PI7)/819.) + LAL_PI5*(64066827104495183/71399020492800 - (1753088*LAL_GAMMA)/12285. - (3506176*log(2.0))/12285.) + 
                LAL_PI3*(157869342918515513264152507759219/4601109448257970699763712000 - (46895104*LAL_GAMMA*LAL_GAMMA)/85995. + LAL_GAMMA*(1039133415323518217819/131945389870694400 - (187580416*log(2.0))/85995.) + 
                (2231261504151126995563*log(2.0))/1.319453898706944e17 - (187580416*log(2.0)*log(2.0))/85995. + (1475651025*log(3.0))/7.8675968e7 - 
                (188720703125*log(5.0))/3.18022848e8 - (876544*Zeta(3))/819.) + 
                LAL_PI*(-464671899642916359718426239920158253645334880173776869/2183863713923805648587725495492767315268730880000 - (5017776128*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/9.029475e6 + LAL_GAMMA*LAL_GAMMA*(2940869073085264664844101/247645003613559552000 - (10035552256*log(2.0))/3.009825e6) + 
                (964617928589311937689039*log(2.0)*log(2.0))/2.251318214668723e19 - (40142209024*log(2.0)*log(2.0)*log(2.0))/9.029475e6 + 
                (806333478108324665949606919515061*log(3.0))/2.85525091278598e28 - (38008842429969538227*log(3.0)*log(3.0))/1.610453714501632e16 + 
                (3789222134222497757392578125*log(5.0))/9.777307303446463e23 - (729405517578125*log(5.0)*log(5.0))/7.27636276224e11 - 
                (45385748308543142507*log(7.0))/4.1709795139584e15 + log(2.0)*(2643785502925434511518863613073178573269/15806222920228063518446840709120000 + (571504883673423551*log(3.0))/7.189525511168e13 + 
                (4134530908203125*log(5.0))/2.546726966784e12 - (187580416*Zeta(3))/28665.) + 
                LAL_GAMMA*(29380557546812856505609944269585099953/322575977963838030988711034880000 + (5835928855017997165517077*log(2.0))/1.2382250180677978e20 - (20071104512*log(2.0)*log(2.0))/3.009825e6 + 
                (85224987219778941*log(3.0))/1.00653357156352e14 - (161884619140625*log(5.0))/8.48908988928e11 - (93790208*Zeta(3))/28665.) + 
                (69285589503891803879*Zeta(3))/2.7488622889728e15 + (876544*Zeta(5))/273.);
                break;

        case 1:
                phase =  -(-1753088*LAL_PI5)/12285. + 
                LAL_PI3*(1039133415323518217819/131945389870694400 - (93790208*LAL_GAMMA)/85995. - (187580416*log(2.0))/85995.) + 
                LAL_PI*(28058305739213426063343984981244619953/322575977963838030988711034880000 - (5017776128*LAL_GAMMA*LAL_GAMMA)/3.009825e6 + LAL_GAMMA*(2940869073085264664844101/123822501806779776000 - (20071104512*log(2.0))/3.009825e6) + 
                (5835928855017997165517077*log(2.0))/1.2382250180677978e20 - (20071104512*log(2.0)*log(2.0))/3.009825e6 + (85224987219778941*log(3.0))/1.00653357156352e14 - 
                (161884619140625*log(5.0))/8.48908988928e11 - (93790208*Zeta(3))/28665.);
                break;
        case 2:
                phase = (-46895104*LAL_PI3)/85995. + LAL_PI*(2940869073085264664844101/247645003613559552000 - (5017776128*LAL_GAMMA)/3.009825e6 - (10035552256*log(2.0))/3.009825e6);
                break;
        case 3:
                phase = -5017776128*LAL_PI/9.029475e6;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_22PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        /* REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        REAL8 LAL_PI8 = LAL_PI7*LAL_PI; */
        
        switch (log_order) {
        case 0: 
                phase =  220960196188958381666196540330729447949423028998442749599099013871/449547269408002629378702340348272796795988118350321418240000 + 
                (184849556184519742*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/4.71573094827681e14 + (2503246985*LAL_PI6)/2.9755593e7 + 
                (1146793130428215898943331967943744232784906465076087*log(2.0))/7.320832199057229e45 - (15177260009295859280002736624951273*log(2.0)*log(2.0))/6.009425375679682e28 + 
                (118524438841916736226*log(2.0)*log(2.0)*log(2.0))/1.1789327370692024e16 + 
                LAL_PI4*(-400006962822751996113799433/22076794083945197076480 + (1550365266533*log(2.0))/1.40276367e9 - (79419015*log(3.0))/326536.) + 
                LAL_GAMMA*LAL_GAMMA*(-5260923679164553267781679699480017/79324414958971809998953021440 + (20557680798999019*LAL_PI2)/2.721922625268e13 + (376338072967388698*log(2.0))/1.57191031609227e14 - 
                (8648289*log(3.0))/571438.) - (6542775761344716772602377222731421069194673831*log(3.0))/3.005549966589526e40 - 
                (1186925982331109249738411581699*log(2.0)*log(3.0))/3.663725270952405e25 - (68106385737*log(2.0)*log(2.0)*log(3.0))/4.000066e6 + 
                (198684417652302157214797072249*log(3.0)*log(3.0))/5.427741142151712e24 + (953398946565*log(2.0)*log(3.0)*log(3.0))/6.4001056e7 - 
                (78519077325*log(3.0)*log(3.0)*log(3.0))/1.8286016e7 + (71937187114627822471097458122389638671875*log(5.0))/6.95280434943778e35 + 
                (1803825083028134839111328125*log(2.0)*log(5.0))/6.674358323675893e22 - (146565216064453125*log(3.0)*log(5.0))/8.415496241152e12 - 
                (351746766944618555908203125*log(5.0)*log(5.0))/4.449572215783929e22 + (8117250537731844765881589515509*log(7.0))/2.986897891429561e26 + 
                (196673050603673959*log(2.0)*log(7.0))/1.15407961812e14 - (6051220875658089481*log(7.0)*log(7.0))/2.077343312616e15 - 
                (617505778333540317602727400433*Zeta(3))/3.88456714652856e24 - (19822665345315676*log(2.0)*Zeta(3))/1.134134427195e12 + (2956687110*log(3.0)*Zeta(3))/285719. + 
                LAL_PI2*(112018034748925817099887461445553460427158431/2100171380374740932490168483773939712000 + (6186546664726613*log(2.0)*log(2.0))/3.88846089324e12 + (32344266549088514295421*log(3.0))/9.489716948946125e18 - 
                (22867132275*log(3.0)*log(3.0))/1.8286016e7 + log(2.0)*(-92924543016958629184078133479/517942286203808022036480+ (6733929735*log(3.0))/2.285752e6) - 
                (4524184659575439453125*log(5.0))/1.9016586475701535e18 - (2025852318599963*log(7.0))/1.8799486992e12 + (470544021557*Zeta(3))/7.855476552e9) + 
                LAL_GAMMA*(208013473201450521511395504870056139943440842269783/2091666342587779598820124704592306435522560000 + 
                (2480855964989*LAL_PI4)/9.81934569e9 + (1900649391343256338*log(2.0)*log(2.0))/3.57252344566425e14 + 
                log(2.0)*(-266656440538787494535687393416459/1133205927985311571413614592 - (4177048122*log(3.0))/2.000033e6) - (96357151405526000591478951359*log(3.0))/7.32745054190481e24 + 
                (64063082655*log(3.0)*log(3.0))/6.4001056e7 + LAL_PI2*(-150651906306722319727711019/1644261226043834990592 + (34060441335203977*log(2.0))/1.360961312634e13 + (520875225*log(3.0))/2.285752e6) - 
                (314136536274139086669921875*log(5.0))/6.674358323675893e22 - (85623268404500477*log(7.0))/4.154686625232e13 + (866966231851660*Zeta(3))/2.26826885439e11) - 
                (11067016500896*Zeta(5))/3.27311523e8;
                break;

        case 1:
                phase =  211056414577400160578027946620265024804816947946583/2091666342587779598820124704592306435522560000 + 
                (184849556184519742*LAL_GAMMA*LAL_GAMMA)/1.57191031609227e14 + (2480855964989*LAL_PI4)/9.81934569e9 + 
                (1900649391343256338*log(2.0)*log(2.0))/3.57252344566425e14 + log(2.0)*(-1323353319337171830962279106275639/5666029639926557857068072960 - (4177048122*log(3.0))/2.000033e6) + 
                LAL_GAMMA*(-5199946535059299492740268157586897/39662207479485904999476510720 + (752676145934777396*log(2.0))/1.57191031609227e14 - (8648289*log(3.0))/285719.) - 
                (93929953385072884825366695359*log(3.0))/7.32745054190481e24 + (64063082655*log(3.0)*log(3.0))/6.4001056e7 + 
                (LAL_PI2*(-2197789913167280400682672005 + 37255501447457685358174208*LAL_GAMMA + 61725777040303382494552064*log(2.0) + 5620392781507082202624000*log(3.0)))/
                2.4663918390657524e22 - (314136536274139086669921875*log(5.0))/6.674358323675893e22 - (85623268404500477*log(7.0))/4.154686625232e13 + 
                (866966231851660*Zeta(3))/2.26826885439e11;
                break;
        case 2:
                phase = -5217498886999411092368990859447761/79324414958971809998953021440 + (184849556184519742*LAL_GAMMA)/1.57191031609227e14 + (20557680798999019*LAL_PI2)/2.721922625268e13 + 
                (376338072967388698*log(2.0))/1.57191031609227e14 - (8648289*log(3.0))/571438.;
                break;
        case 3:
                phase = 184849556184519742/4.71573094827681e14;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_23PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI3*LAL_PI;
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI5*LAL_PI;
        REAL8 LAL_PI7 = LAL_PI6*LAL_PI;
        //REAL8 LAL_PI8 = LAL_PI7*LAL_PI;
        
        switch (log_order) {
        case 0: 
                phase =  -((5408768*LAL_PI7)/178605.) + LAL_PI5*(2368300727797389721249/498255318392832000 - (922447616*LAL_GAMMA)/2.679075e6 - (1729854208*log(2.0))/2.679075e6 - (1872*log(3.0))/49.) + 
                LAL_PI3*(1308701768722814998919069541415328041531/33706071418536520477017039175680000 - 
                (26304004336*LAL_GAMMA*LAL_GAMMA)/1.8753525e7 + LAL_GAMMA*(869271806706251539139663/22443910817005117440 - (20927981408*log(2.0))/3.750705e6) + 
                (117150134079060691691488337*log(2.0))/1.683293311275384e21 - (20804539184*log(2.0)*log(2.0))/3.750705e6 - (4803045697534840583*log(3.0))/1.6892871131136e15 + 
                (1064799951171875*log(5.0))/2.13711353856e11 - (505548832*Zeta(3))/178605.) + 
                LAL_PI*(-1047439268696043528371334210034399153978237698819404827052699/1764507284257586868917667507220768671554252832768000000 -
                (2696412641228*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/1.969120125e9 - (341415141668*log(2.0)*log(2.0)*log(2.0))/5.6260575e7 + 
                log(2.0)*log(2.0)*(188904811985580567058410139013/789832783399996496160000 - (41262*log(3.0))/5.) + LAL_GAMMA*LAL_GAMMA*
                (2168564482717123206808055407/37611084923809356960000 - (1010530327868*log(2.0))/1.31274675e8 - (5728398*log(3.0))/12005.) + (875492170576714639718946002286610508117*log(3.0))/1.1864424117899583e34 - 
                (2781178256440427416431*log(3.0)*log(3.0))/1.409147000188928e18 - (1779570*log(3.0)*log(3.0)*log(3.0))/2401. - 
                (8515893542439435586808508373046875*log(5.0))/1.000643850010273e29 + (188720703125*log(3.0)*log(5.0))/3.6427776e7 + 
                (482428809326171875*log(5.0)*log(5.0))/6.1121447202816e13 + (3012031912014675879845567569*log(7.0))/6.586035046253509e22 + 
                log(2.0)*(428985521134188225017836717007721598925764499/3405608790392338565684556299186995200000 - 
                (169631042728106894729*log(3.0))/9.7058594400768e15 + (209898*log(3.0)*log(3.0))/49. - 
                (1623992989970703125*log(5.0))/1.06962532604928e14 - (138757016872*Zeta(3))/6.251175e6) + 
                LAL_GAMMA*(261146505082218427275952669906707148949678419/3405608790392338565684556299186995200000 - 
                (349800451436*log(2.0)*log(2.0))/2.6254935e7 + log(2.0)*(8707611690400644990818112433/39491639169999824808000 - (972348*log(3.0))/245.) + 
                (136786574381887610111*log(3.0))/3.7745008933632e16 + (2473146*log(3.0)*log(3.0))/2401. + (103107775068359375*log(5.0))/3.5654177534976e13 - 
                (61686728168*Zeta(3))/6.251175e6) + (19634731061005124618107*Zeta(3))/1.62354678942456e17 + 2511.918367346939*log(3.0)*Zeta(3) + (1008387928*Zeta(5))/59535.);
                break;

        case 1:
                phase =  (-922447616*LAL_PI5)/2.679075e6 + LAL_PI3*(850766856335240666795087/22443910817005117440 - (52608008672*LAL_GAMMA)/1.8753525e7 - (20927981408*log(2.0))/3.750705e6) + 
                LAL_PI*(218629626707380920834224026582302690835425619/3405608790392338565684556299186995200000 - 
                (2696412641228*LAL_GAMMA*LAL_GAMMA)/6.56373375e8 - (349800451436*log(2.0)*log(2.0))/2.6254935e7 + 
                log(2.0)*(8632954516260756234381611953/39491639169999824808000 - (972348*log(3.0))/245.) + LAL_GAMMA*(2144863792513983919050436207/18805542461904678480000 - (2021060655736*log(2.0))/1.31274675e8 - (11456796*log(3.0))/12005.) + 
                (136786574381887610111*log(3.0))/3.7745008933632e16 + (2473146*log(3.0)*log(3.0))/2401. + (103107775068359375*log(5.0))/3.5654177534976e13 - 
                (61686728168*Zeta(3))/6.251175e6);
                break;
        case 2:
                phase = (-26304004336*LAL_PI3)/1.8753525e7 + LAL_PI*(2144863792513983919050436207/37611084923809356960000 - (2696412641228*LAL_GAMMA)/6.56373375e8 - (1010530327868*log(2.0))/1.31274675e8 - 
                  (5728398*log(3.0))/12005.);
                break;
        case 3:
                phase = -2696412641228*LAL_PI/1.969120125e9;
                break;

        }
        return  phase;
}

static REAL8 UNUSED
XLALSimInspiralTaylorF2Phasing_24PNCoeff(INT4 log_order)
{
        REAL8 phase = 0.0;
        REAL8 LAL_PI2 = LAL_PI*LAL_PI; 
        //REAL8 LAL_PI3 = LAL_PI2*LAL_PI;
        REAL8 LAL_PI4 = LAL_PI2*LAL_PI2;
        REAL8 LAL_PI5 = LAL_PI4*LAL_PI;
        REAL8 LAL_PI6 = LAL_PI4*LAL_PI2;
        REAL8 LAL_PI7 = LAL_PI5*LAL_PI2;
        REAL8 LAL_PI8 = LAL_PI4*LAL_PI4;
        
        switch (log_order) {
        case 0: 
                phase =  65647839913816173164493161700101686508707300915591175729603477865160043904011/72948563642860031848691653918712229570552139122585242899484088729600000 + 
                (536902045696*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/1.385677125e9 + (65536*LAL_PI8)/10773. + 
                (32543417811083944120447931463850034809476747172636937195472029*log(2.0))/3.0343806978587067e55 - 
                (1748064471932316636981200302991639077167093221*log(2.0)*log(2.0))/4.422561415301162e39 - (2977919134257534212534529866201*log(2.0)*log(2.0)*log(2.0))/2.9040981322975796e25 + 
                (8590432731136*log(2.0)*log(2.0)*log(2.0)*log(2.0))/1.385677125e9 + LAL_PI6*(-559818838915769777/1067605472102400 + (28049408*log(2.0))/125685.) + 
                LAL_GAMMA*LAL_GAMMA*LAL_GAMMA*(-2516574060639912881119499357/240008110107237992520000 + (40142209024*LAL_PI*LAL_PI)/3.9590775e7 + (4295216365568*log(2.0))/1.385677125e9) + 
                (611505889076804355670597062183574895283779637963053864989*log(3.0))/9.866400683989088e50 - 
                (555616823451183461410496575643364981443*log(2.0)*log(3.0))/2.5298361817239417e34 + (1814733014897445162816209*log(2.0)*log(2.0)*log(3.0))/3.313256884194217e19 - 
                (1897220307891557272518086441430586358427*log(3.0)*log(3.0))/1.4991621817623357e34 - (48858412043337699489279*log(2.0)*log(3.0)*log(3.0))/7.362793075987149e17 + 
                (352575782581627247981607*log(3.0)*log(3.0)*log(3.0))/1.6829241316542054e19 + (3599763503908304539168027921224521272746511328125*log(5.0))/8.190488069738595e43 - 
                (107325532406009880454362363962890625*log(2.0)*log(5.0))/3.638452221287354e30 + (90580130043095703125*log(2.0)*log(2.0)*log(5.0))/1.1177584657214976e16 + 
                (3213439242370148779296875*log(3.0)*log(5.0))/4.829067574258749e19 - (380271898748499099127645122412109375*log(5.0)*log(5.0))/9.702539256766278e30 - 
                (15979961960205078125*log(2.0)*log(5.0)*log(5.0))/1.596797808173568e15 + (2819152325439453125*log(5.0)*log(5.0)*log(5.0))/6.84341917788672e14 - 
                (255327589251374254887597927481693226444269*log(7.0))/5.602697089338102e35 - (37527733708977518967989*log(2.0)*log(7.0))/1.894537101105792e18 + 
                (41343924869387*log(3.0)*log(7.0))/7.56760576e9 + (221940628056511747894637*log(7.0)*log(7.0))/7.578148404423168e18 + 
                (71390211212258470496382517*log(11.0))/4.021936132447509e22 - (18499995898301784053259835637881432629977*Zeta(3))/5.486583148271472e34 + 
                (221231504974634899261381*log(2.0)*Zeta(3))/1.0745173834674879e20 + (80284418048*log(2.0)*log(2.0)*Zeta(3))/4.398975e6 - 
                (697774173785667217*log(3.0)*Zeta(3))/7.3554376383488e13 - (36388408837890625*log(5.0)*Zeta(3))/6.20356568832e11 + (187580416*Power(Zeta(3),2))/41895. + 
                LAL_PI4*(-69652926963465156127275068980602335357/1650778969216074330783224954880000 - 
                (46059996383356816527397*log(2.0))/2.290013737659648e18 + (1500643328*log(2.0)*log(2.0))/628425. + 
                (61856875841407*log(3.0))/5.3440651264e11 + (3359228515625*log(5.0))/1.859210496e9 + (7012352*Zeta(3))/5985.) + 
                LAL_GAMMA*LAL_GAMMA*(-30804053940640623439867345503568384606414037/231483930828435470689417442053324800000 + 
                (375160832*LAL_PI4)/628425. - (601898355141185362998081151169*log(2.0))/9.680327107658598e24 + 
                (4295216365568*log(2.0)*log(2.0))/4.61892375e8 + LAL_PI2*(-178056339765745819397645987/8252293505030307532800 + (80284418048*log(2.0))/1.3196925e7) - 
                (41322005827908123591503*log(3.0))/3.313256884194217e19 + (138864626298828125*log(5.0))/1.241953850801664e15 + (20071104512*Zeta(3))/4.398975e6) + 
                LAL_GAMMA*(15445645271057135917953181946342743021961488190030461400235377/30343806978587066780530907201720116807816074756096000000 + 
                (14024704*LAL_PI6)/125685. - (1152745237447568759293980009613*log(2.0)*log(2.0))/9.680327107658598e24 + 
                (17180865462272*log(2.0)*log(2.0)*log(2.0))/1.385677125e9 + LAL_PI4*(-18792018737252475637861/2290013737659648000 + (1500643328*log(2.0))/628425.) - 
                (74427598315294288500175766565799682299*log(3.0))/1.2046838960590197e33 + (613652094902831301159*log(3.0)*log(3.0))/5.259137911419392e17 - 
                (8716341642547157138442233861328125*log(5.0))/3.638452221287354e30 + (625684052978515625*log(5.0)*log(5.0))/5.32265936057856e14 + 
                (6985987213214917114627*log(7.0))/3.15756183517632e17 - (7213566040081376589462647*Zeta(3))/1.0745173834674879e20 + 
                LAL_PI2*(-492360040122621169127586993970256277289111/2633559911170306583150450536611840000 -
                 (355507166323656162034131299*log(2.0))/4.1261467525151536e21 + (160568836096*log(2.0)*log(2.0))/1.3196925e7 - 
                (3715098855731414131*log(3.0))/3.530610066407424e15 + (161884619140625*log(5.0))/2.65867100928e11 + (750321664*Zeta(3))/125685.) + 
                log(2.0)*(-7294968183929200723600096842362685658314198029/14004777815120345976709755244226150400000 - 
                (142402901279390045029823*log(3.0))/1.6566284420971086e19 - (3546600613056640625*log(5.0))/1.862930776202496e15 + 
                (80284418048*Zeta(3))/4.398975e6) - (375160832*Zeta(5))/41895.) + 
                LAL_PI2*(11850899254412642868832473402164045727848747377540432981/34987047369518957358291519402643002994615910400000 - 
                (642131479357328417739102883*log(2.0)*log(2.0))/8.252293505030307e21 + (321137672192*log(2.0)*log(2.0)*log(2.0))/3.9590775e7 - 
                (30099488951015763225845030440159*log(3.0))/1.1923025789655742e27 + (60860276209445290659*log(3.0)*log(3.0))/1.176870022135808e16 - 
                (30755334759378824984178408203125*log(5.0))/1.140334871806579e27 + (5105838623046875*log(5.0)*log(5.0))/1.595202605568e12 + 
                (226583189018657375989*log(7.0))/2.74322114187264e16 - (20910687278522270007343*Zeta(3))/4.580027475319296e17 + 
                log(2.0)*(-258973189552891961137187812165771265961599/877853303723435527716816845537280000 - 
                (281690099080898123087*log(3.0))/1.765305033203712e16 - (4134530908203125*log(5.0))/7.97601302784e11 + (1500643328*Zeta(3))/125685.) - 
                (7012352*Zeta(5))/1197.) + (9470941644347532641*Zeta(5))/2.3854309767288e14 - (750321664*log(2.0)*Zeta(5))/41895. + (3506176*Zeta(7))/399.;
                 break;

        case 1:
                phase =  14379922693084161928543600440239695177524162323859423949995377/30343806978587066780530907201720116807816074756096000000 +
                 (2147608182784*LAL_GAMMA*LAL_GAMMA*LAL_GAMMA)/1.385677125e9 + (14024704*LAL_PI6)/125685. - 
                (1152745237447568759293980009613*log(2.0)*log(2.0))/9.680327107658598e24 + (17180865462272*log(2.0)*log(2.0)*log(2.0))/1.385677125e9 + 
                LAL_PI4*(-18792018737252475637861/2290013737659648000 + (750321664*LAL_GAMMA)/628425. + (1500643328*log(2.0))/628425.) + 
                LAL_GAMMA*LAL_GAMMA*(-2516574060639912881119499357/80002703369079330840000 + (4295216365568*log(2.0))/4.61892375e8) - (72247773595256810060415457636304642299*log(3.0))/1.2046838960590197e33 + 
                (613652094902831301159*log(3.0)*log(3.0))/5.259137911419392e17 - (8716341642547157138442233861328125*log(5.0))/3.638452221287354e30 + 
                (625684052978515625*log(5.0)*log(5.0))/5.32265936057856e14 + (6985987213214917114627*log(7.0))/3.15756183517632e17 - 
                (7213566040081376589462647*Zeta(3))/1.0745173834674879e20 + LAL_PI*LAL_PI*
                (-452669722926858538985057932999553347369111/2633559911170306583150450536611840000 + 
                (40142209024*LAL_GAMMA*LAL_GAMMA)/1.3196925e7 - (355507166323656162034131299*log(2.0))/4.1261467525151536e21 + 
                (160568836096*log(2.0)*log(2.0))/1.3196925e7 + LAL_GAMMA*(-178056339765745819397645987/4126146752515153766400 + (160568836096*log(2.0))/1.3196925e7) - 
                (3715098855731414131*log(3.0))/3.530610066407424e15 + (161884619140625*log(5.0))/2.65867100928e11 + (750321664*Zeta(3))/125685.) + 
                LAL_GAMMA*(-29712127388445473792380423248420263149774037/115741965414217735344708721026662400000 - 
                (601898355141185362998081151169*log(2.0))/4.840163553829299e24 + (8590432731136*log(2.0)*log(2.0))/4.61892375e8 - 
                (41322005827908123591503*log(3.0))/1.6566284420971086e19 + (138864626298828125*log(5.0))/6.20976925400832e14 + (40142209024*Zeta(3))/4.398975e6) + 
                log(2.0)*(-7124858976678935857257208083643082249418198029/14004777815120345976709755244226150400000 - 
                (142402901279390045029823*log(3.0))/1.6566284420971086e19 - (3546600613056640625*log(5.0))/1.862930776202496e15 + 
                (80284418048*Zeta(3))/4.398975e6) - (375160832*Zeta(5))/41895.;
                break;
        case 2:
                phase = -29884083054862897652172035022704438264014037/231483930828435470689417442053324800000 + 
                (1073804091392*LAL_GAMMA*LAL_GAMMA)/4.61892375e8 + 
                (375160832*LAL_PI4)/628425. - (601898355141185362998081151169*log(2.0))/9.680327107658598e24 + (4295216365568*log(2.0)*log(2.0))/4.61892375e8 + 
                LAL_PI2*(-178056339765745819397645987/8252293505030307532800 + (80284418048*log(2.0))/1.3196925e7) + 
                LAL_GAMMA*(-2516574060639912881119499357/80002703369079330840000 + (40142209024*LAL_PI2)/1.3196925e7 + (4295216365568*log(2.0))/4.61892375e8) - 
                (41322005827908123591503*log(3.0))/3.313256884194217e19 + (138864626298828125*log(5.0))/1.241953850801664e15 + (20071104512*Zeta(3))/4.398975e6;
                break;
        case 3:
                phase = -2516574060639912881119499357/240008110107237992520000 + (2147608182784*LAL_GAMMA)/1.385677125e9 + (40142209024*LAL_PI2)/3.9590775e7 + (4295216365568*log(2.0))/1.385677125e9;
                break;
        case 4:
                phase = 536902045696/1.385677125e9;
        }
        return  phase;
}

static void UNUSED
XLALSimInspiralPNBHPPhasing_F2(
	PNPhasingSeries *pfa, /**< \todo UNDOCUMENTED */
	const REAL8 m1, /**< Mass of body 1, in Msol */
	const REAL8 m2, /**< Mass of body 2, in Msol */
	const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
	const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
	const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
	const REAL8 chi1dotchi2, /**< Dot product of dimensionles spin 1 and spin 2 */
	LALDict *p /**< LAL dictionary containing accessory parameters */
	)
{
    const REAL8 mtot = m1 + m2;
    const REAL8 eta = m1*m2/mtot/mtot;
    const REAL8 m1M = m1/mtot;
    const REAL8 m2M = m2/mtot;

    const REAL8 pfaN = 3.L/(128.L * eta);

    memset(pfa, 0, sizeof(PNPhasingSeries));

    pfa->v[0] = 1.L;
    pfa->v[1] = 0.L;
    pfa->v[2] = XLALSimInspiralTaylorF2Phasing_2PNCoeff(eta);
    pfa->v[3] = XLALSimInspiralTaylorF2Phasing_3PNCoeff(eta);
    pfa->v[4] = XLALSimInspiralTaylorF2Phasing_4PNCoeff(eta);
    pfa->v[5] = XLALSimInspiralTaylorF2Phasing_5PNCoeff(eta);
    pfa->vlogv[5] = XLALSimInspiralTaylorF2Phasing_5PNLogCoeff(eta);
    pfa->v[6] =  XLALSimInspiralTaylorF2Phasing_6PNCoeff(eta);
    pfa->vlogv[6] = XLALSimInspiralTaylorF2Phasing_6PNLogCoeff(eta);
    pfa->v[7] = XLALSimInspiralTaylorF2Phasing_7PNCoeff(eta);
    pfa->v[8] = XLALSimInspiralTaylorF2Phasing_8PNCoeff(0);
    pfa->vlogv[8] = XLALSimInspiralTaylorF2Phasing_8PNCoeff(1);
    pfa->vlogvsq[8] = XLALSimInspiralTaylorF2Phasing_8PNCoeff(2);
    pfa->v[9] = XLALSimInspiralTaylorF2Phasing_9PNCoeff(0);
    pfa->vlogv[9] = XLALSimInspiralTaylorF2Phasing_9PNCoeff(1);
    pfa->v[10] = XLALSimInspiralTaylorF2Phasing_10PNCoeff(0);
    pfa->vlogv[10] = XLALSimInspiralTaylorF2Phasing_10PNCoeff(1);
    pfa->v[11] = XLALSimInspiralTaylorF2Phasing_11PNCoeff(0);
    pfa->vlogv[11] = XLALSimInspiralTaylorF2Phasing_11PNCoeff(1);
    pfa->v[12] = XLALSimInspiralTaylorF2Phasing_12PNCoeff(0);
    pfa->vlogv[12] = XLALSimInspiralTaylorF2Phasing_12PNCoeff(1);
    pfa->vlogvsq[12] = XLALSimInspiralTaylorF2Phasing_12PNCoeff(2);
    pfa->v[13] = XLALSimInspiralTaylorF2Phasing_13PNCoeff(0);
    pfa->vlogv[13] = XLALSimInspiralTaylorF2Phasing_13PNCoeff(1);
    pfa->v[14] = XLALSimInspiralTaylorF2Phasing_14PNCoeff(0);
    pfa->vlogv[14] = XLALSimInspiralTaylorF2Phasing_14PNCoeff(1);
    pfa->vlogvsq[14] = XLALSimInspiralTaylorF2Phasing_14PNCoeff(2);
    pfa->v[15] = XLALSimInspiralTaylorF2Phasing_15PNCoeff(0);
    pfa->vlogv[15] = XLALSimInspiralTaylorF2Phasing_15PNCoeff(1);
    pfa->vlogvsq[15] = XLALSimInspiralTaylorF2Phasing_15PNCoeff(2);
    pfa->v[16] = XLALSimInspiralTaylorF2Phasing_16PNCoeff(0);
    pfa->vlogv[16] = XLALSimInspiralTaylorF2Phasing_16PNCoeff(1);
    pfa->vlogvsq[16] = XLALSimInspiralTaylorF2Phasing_16PNCoeff(2);
    pfa->v[17] = XLALSimInspiralTaylorF2Phasing_17PNCoeff(0);
    pfa->vlogv[17] = XLALSimInspiralTaylorF2Phasing_17PNCoeff(1);
    pfa->vlogvsq[17] = XLALSimInspiralTaylorF2Phasing_17PNCoeff(2);
    pfa->v[18] = XLALSimInspiralTaylorF2Phasing_18PNCoeff(0);
    pfa->vlogv[18] = XLALSimInspiralTaylorF2Phasing_18PNCoeff(1);
    pfa->vlogvsq[18] = XLALSimInspiralTaylorF2Phasing_18PNCoeff(2);
    pfa->vlogvcu[18] = XLALSimInspiralTaylorF2Phasing_18PNCoeff(3);
    pfa->v[19] = XLALSimInspiralTaylorF2Phasing_19PNCoeff(0);
    pfa->vlogv[19] = XLALSimInspiralTaylorF2Phasing_19PNCoeff(1);
    pfa->vlogvsq[19] = XLALSimInspiralTaylorF2Phasing_19PNCoeff(2);
    pfa->v[20] = XLALSimInspiralTaylorF2Phasing_20PNCoeff(0);
    pfa->vlogv[20] = XLALSimInspiralTaylorF2Phasing_20PNCoeff(1);
    pfa->vlogvsq[20] = XLALSimInspiralTaylorF2Phasing_20PNCoeff(2);
    pfa->vlogvcu[20] = XLALSimInspiralTaylorF2Phasing_20PNCoeff(3);
    pfa->v[21] = XLALSimInspiralTaylorF2Phasing_21PNCoeff(0);
    pfa->vlogv[21] = XLALSimInspiralTaylorF2Phasing_21PNCoeff(1);
    pfa->vlogvsq[21] = XLALSimInspiralTaylorF2Phasing_21PNCoeff(2);
    pfa->vlogvcu[21] = XLALSimInspiralTaylorF2Phasing_21PNCoeff(3);
    pfa->v[22] = XLALSimInspiralTaylorF2Phasing_22PNCoeff(0);
    pfa->vlogv[22] = XLALSimInspiralTaylorF2Phasing_22PNCoeff(1);
    pfa->vlogvsq[22] = XLALSimInspiralTaylorF2Phasing_22PNCoeff(2);
    pfa->vlogvcu[22] = XLALSimInspiralTaylorF2Phasing_22PNCoeff(3);
    pfa->v[23] = XLALSimInspiralTaylorF2Phasing_23PNCoeff(0);
    pfa->vlogv[23] = XLALSimInspiralTaylorF2Phasing_23PNCoeff(1);
    pfa->vlogvsq[23] = XLALSimInspiralTaylorF2Phasing_23PNCoeff(2);
    pfa->vlogvcu[23] = XLALSimInspiralTaylorF2Phasing_23PNCoeff(3);
    pfa->v[24] = XLALSimInspiralTaylorF2Phasing_24PNCoeff(0);
    pfa->vlogv[24] = XLALSimInspiralTaylorF2Phasing_24PNCoeff(1);
    pfa->vlogvsq[24] = XLALSimInspiralTaylorF2Phasing_24PNCoeff(2);
    pfa->vlogvcu[24] = XLALSimInspiralTaylorF2Phasing_24PNCoeff(3);
    pfa->vlogvquar[24] = XLALSimInspiralTaylorF2Phasing_24PNCoeff(4);

    /* modify the PN coefficients if a non null LALSimInspiralTestGRParam structure is passed */
    /* BEWARE: this is for the non-spinning case only!*/
    pfa->v[0]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi0(p));
    pfa->v[1] = XLALSimInspiralWaveformParamsLookupNonGRDChi1(p);
    pfa->v[2]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi2(p));
    pfa->v[3]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi3(p));
    pfa->v[4]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi4(p));
    pfa->v[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5(p));
    pfa->vlogv[5]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi5L(p));
    pfa->v[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6(p));
    pfa->vlogv[6]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi6L(p));
    pfa->v[7]*=(1.0+XLALSimInspiralWaveformParamsLookupNonGRDChi7(p));

    const REAL8 qm_def1=1.+XLALSimInspiralWaveformParamsLookupdQuadMon1(p);
    const REAL8 qm_def2=1.+XLALSimInspiralWaveformParamsLookupdQuadMon2(p);

    switch( XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) )
    {
        case LAL_SIM_INSPIRAL_SPIN_ORDER_ALL:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_35PN:
	    pfa->v[7] += XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_7PNSOCoeff(m2M)*chi2L;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_3PN:
            pfa->v[6] += XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m1M)*chi1L
	      + XLALSimInspiralTaylorF2Phasing_6PNSOCoeff(m2M)*chi2L
	      + XLALSimInspiralTaylorF2Phasing_6PNS1S2OCoeff(eta)*chi1L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m1M))*chi1sq
	      + (XLALSimInspiralTaylorF2Phasing_6PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_6PNSelf2SCoeff(m2M))*chi2sq;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_25PN:
            pfa->v[5] += XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
	      + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L;
            pfa->vlogv[5] += 3.*(XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m1M)*chi1L
				 + XLALSimInspiralTaylorF2Phasing_5PNSOCoeff(m2M)*chi2L);
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_2PN:
	    /* 2PN SS, QM, and self-spin */
            pfa->v[4] += XLALSimInspiralTaylorF2Phasing_4PNS1S2Coeff(eta)*chi1dotchi2+XLALSimInspiralTaylorF2Phasing_4PNS1S2OCoeff(eta)*chi1L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m1M))*chi1L*chi1L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SOCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SOCoeff(m2M))*chi2L*chi2L
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m1M)*qm_def1+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m1M))*chi1sq
	      + (XLALSimInspiralTaylorF2Phasing_4PNQM2SCoeff(m2M)*qm_def2+XLALSimInspiralTaylorF2Phasing_4PNSelf2SCoeff(m2M))*chi2sq;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_15PN:
            pfa->v[3] += XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m1M)*chi1L+XLALSimInspiralTaylorF2Phasing_3PNSOCoeff(m2M)*chi2L;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_SPIN_ORDER_1PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_05PN:
        case LAL_SIM_INSPIRAL_SPIN_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid spin PN order %i\n",
			   __func__, XLALSimInspiralWaveformParamsLookupPNSpinOrder(p) );
            XLAL_ERROR_VOID(XLAL_EINVAL);
            break;
    }

    REAL8 lambda1=XLALSimInspiralWaveformParamsLookupTidalLambda1(p);
    REAL8 lambda2=XLALSimInspiralWaveformParamsLookupTidalLambda2(p);
    switch( XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) )
    {
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pfa->v[15] = (lambda1*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_15PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT:
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pfa->v[14] = (lambda1*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_14PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pfa->v[13] = (lambda1*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_13PNTidalCoeff(m2M));
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
            pfa->v[12] = (lambda1*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_12PNTidalCoeff(m2M) );
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pfa->v[10] = ( lambda1*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m1M) + lambda2*XLALSimInspiralTaylorF2Phasing_10PNTidalCoeff(m2M) );
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLALPrintError("XLAL Error - %s: Invalid tidal PN order %i\n",
                           __func__, XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) );
            XLAL_ERROR_VOID(XLAL_EINVAL);
    }


    /* At the very end, multiply everything in the series by pfaN */
    for(int ii = 0; ii <= PN_PHASING_SERIES_MAX_ORDER; ii++)
    {
        pfa->v[ii] *= pfaN;
        pfa->vlogv[ii] *= pfaN;
        pfa->vlogvsq[ii] *= pfaN;
        pfa->vlogvcu[ii] *= pfaN;
        pfa->vlogvquar[ii] *= pfaN;
    }
}

/*PN corrections to eccentric phasing due to BHP theory*/


static INT4 UNUSED
eccentricityPNBHPCoeffs_F2(REAL8 eccPNBHPCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1])
{
  INT4 ret = 0;
  memset(eccPNBHPCoeffs, 0x00, (LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*sizeof(REAL8));
  eccPNBHPCoeffs[7][7][0] = (7374744898223317*LAL_PI)/2.815609724928e13; // lowest order term
  eccPNBHPCoeffs[7][5][2] = (-814045990347257*LAL_PI)/6.78925836288e12; //v^2 term
  eccPNBHPCoeffs[7][4][3] = (-932161278950603*LAL_PI)/9.3998047739904e13; //v*v0 term
  eccPNBHPCoeffs[7][3][4] = (-46517251019233*LAL_PI)/6.8913524736e12; //v0^2 term
  eccPNBHPCoeffs[7][2][5] = (2957938456447*LAL_PI)/8.994004992e10; //v0^2 term
  eccPNBHPCoeffs[7][0][7] = (-41144741521*LAL_PI)/8.77879296e9; //v0^2 term

  eccPNBHPCoeffs[8][8][0] = -10520.658458118267 + (4004693587*LAL_GAMMA)/6.3504e6 + (38662427350307*LAL_PI*LAL_PI)/3.463907328e11 - (98926465489457*log(2.0))/5.68297296e10 + 
                                (3006309464757*log(3.0))/1.8709376e9; //v0^2 term
  eccPNBHPCoeffs[8][6][2] = -79.22408540346447 - (2350835165449*LAL_GAMMA)/3.9880512e10 + (7960580380771*LAL_PI*LAL_PI)/8.750923776e10 + (36095258045573*log(2.0))/1.19641536e11 - 
                         (65811862017*log(3.0))/1.5755264e8; //v0^2 term
  eccPNBHPCoeffs[8][5][3] = (-108328746332833*LAL_PI*LAL_PI)/4.8494702592e11; //v0^2 term
  eccPNBHPCoeffs[8][4][4] = 14.109429464261373; //v0^2 term
  eccPNBHPCoeffs[8][3][5] = (-10839481187249*LAL_PI*LAL_PI)/1.36733184e11; //v0^2 term
  eccPNBHPCoeffs[8][2][6] = 890.7314812237001 - (3062544163169*LAL_GAMMA)/2.0658105216e10 + (2229743359003*LAL_PI*LAL_PI)/4.1972023296e10 - (11163467433487*log(2.0))/3.0987157824e11 - 
                         (2667377174373*log(3.0))/1.020153344e10; //v0^2 term
  eccPNBHPCoeffs[8][0][8] = 226.69952142592206 - (9397061*LAL_GAMMA)/254016. + (684616409*LAL_PI*LAL_PI)/6.967296e7 - (34253803*log(2.0))/3.81024e6 - (8184537*log(3.0))/125440.; //v0^2 term

  eccPNBHPCoeffs[9][9][0] = (61852842716547687233*LAL_PI)/7.780061288595456e16 + (156796692229*LAL_GAMMA*LAL_PI)/1.59522048e8 + (3778827805673*LAL_PI*LAL_PI*LAL_PI)/2.6252771328e10 - 
                        (5170229816917*LAL_PI*log(2.0))/2.39283072e9 + (46201139343*LAL_PI*log(3.0))/1.125376e7; //v^4 term
  eccPNBHPCoeffs[9][7][2] = (34275389625395313173*LAL_PI)/9.933471109545984e16; //v^2*v0^2 term
  eccPNBHPCoeffs[9][6][3] = (-2121859507587722745323*LAL_PI)/5.4579511590912e18 - (312836165681*LAL_GAMMA*LAL_PI)/2.848608e9 + (1059350089499*LAL_PI*LAL_PI*LAL_PI)/6.25065984e9 + 
                         (4803357671437*LAL_PI*log(2.0))/8.545824e9 - (8757879273*LAL_PI*log(3.0))/1.125376e7; //v^2*v0^2 term
  eccPNBHPCoeffs[9][5][4] = (-4659240894806141201*LAL_PI)/5.279327302975488e16; //v^2*v0^2 term
  eccPNBHPCoeffs[9][4][5] = (-1453310632582456687*LAL_PI)/4.7375016060911616e17; //v^2*v0^2 term
  eccPNBHPCoeffs[9][3][6] = (-82066325786060127229*LAL_PI)/9.5514145284096e16 + (10259875723*LAL_GAMMA*LAL_PI)/7.12152e7 - (373870175879*LAL_PI*LAL_PI*LAL_PI)/3.9066624e9 + (37398901829*LAL_PI*log(2.0))/1.068228e9 + 
                        (8936020791*LAL_PI*log(3.0))/3.5168e7; //v^2*v0^2 term
  eccPNBHPCoeffs[9][2][7] = (-634576262425854541*LAL_PI)/4.497847932469248e16; //v^2*v0^2 term
  eccPNBHPCoeffs[9][0][9] = (40733897021722489*LAL_PI)/9.733925634048e13 - (1250509*LAL_GAMMA*LAL_PI)/18144. + (39033449*LAL_PI*LAL_PI*LAL_PI)/2.985984e6 - (4558307*LAL_PI*log(2.0))/272160. - 
                        (1089153*LAL_PI*log(3.0))/8960.; //v^2*v0^2 term

  eccPNBHPCoeffs[10][10][0] = -69490.07946696543 - (1314393019969101407*LAL_GAMMA)/6.8982438260736e14 - (454155027147921331147*LAL_PI*LAL_PI)/5.16024473223168e16 + 
   (8568141654792000809*log(2.0))/1.8813392252928e14 - (6896278250779531*log(3.0))/3.6048044032e11 - (280371337890625*log(5.0))/1.7108739648e10; //v^5 term
  eccPNBHPCoeffs[10][8][2] = -22644.12139010655 + (3472157632121443*LAL_GAMMA)/5.4556540416e12 + (852681317276568653*LAL_PI*LAL_PI)/2.0949711519744e15 - 
   (592143456702332033*log(2.0))/1.145687348736e14 + (1053007765812717*log(3.0))/4.190900224e11; //v^3*v0^2 term
  eccPNBHPCoeffs[10][7][3] = (1990575122973983471*LAL_PI*LAL_PI)/2.18318046363648e15; //v^2*v0^3 term
  eccPNBHPCoeffs[10][6][4] = 713.4284092696258 - (41869030336962953*LAL_GAMMA)/4.82394673152e14 + (63327323239992799*LAL_PI*LAL_PI)/2.2682394427392e15 - 
   (94224174034169723*log(2.0))/1.447184019456e15 - (193009647305323*log(3.0))/1.764589568e12 ; //v^2*v0^3 term
  eccPNBHPCoeffs[10][5][5] = (-7039178702418597241*LAL_PI*LAL_PI)/1.46647980638208e16; //v^2*v0^3 term
  eccPNBHPCoeffs[10][4][6] = 1629.9679707284492 - (202099492979647*LAL_GAMMA)/7.58048772096e11 - (238614585241686307*LAL_PI*LAL_PI)/2.707143774909235e16 - 
   (22837242706700111*log(2.0))/3.5249267902464e14 - (202099492979647*log(3.0))/4.2980360192e11 ; //v^2*v0^3 term
  eccPNBHPCoeffs[10][3][7] = (127265581748693999*LAL_PI*LAL_PI)/4.961773780992e15; //v^2*v0^3 term
  eccPNBHPCoeffs[10][2][8] = 1273.8615495415286 - (8676187614257777*LAL_GAMMA)/4.1646740115456e13 + (632097674824441013*LAL_PI*LAL_PI)/1.142310586023936e16 - 
   (31626103239068671*log(2.0))/6.2470110173184e14 - (839631059444301*log(3.0))/2.28514349056e12; //v^2*v0^3 term
  eccPNBHPCoeffs[10][7][3] = (1990575122973983471*LAL_PI*LAL_PI)/2.18318046363648e15; //v^2*v0^3 term
  eccPNBHPCoeffs[10][0][10] = -189.68296819130052 + (95697675707*LAL_GAMMA)/3.072577536e9 - (71038219571281*LAL_PI*LAL_PI)/1.26414618624e13 + (348833463061*log(2.0))/4.608866304e10 + 
   (3087021797*log(3.0))/5.619712e7; //v^2*v0^3 term

  //printPNCoeffs_F2(eccPNCoeffs);
  return ret;
}



/**
 * Compute eccentric phase correction term due to BHP theory using eccPNBHPCeoffs[k][i][j]
 *
 */

static REAL8 UNUSED
eccentricityPhasingBHP_F2(REAL8 v, REAL8 v0, REAL8 ecc, REAL8 eta, INT4 ecc_order)
{
  static REAL8 v0_power[LAL_MAX_ECC_PN_ORDER+1];
  /* following code is not efficient in memory usage, need to be improved later */
  static REAL8 eccPNBHPCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]; // we want to calculate just one time
  REAL8 v_power[LAL_MAX_ECC_PN_ORDER+1];
  REAL8 phasing = 0.0;
  REAL8 global_factor;
  v0_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v0_power[i] = v0_power[i-1]*v0;
  }
  eccentricityPNBHPCoeffs_F2(eccPNBHPCoeffs);
  //printPNCoeffs_F2(eccPNCoeffs);
  v_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v_power[i] = v_power[i-1]*v;
  }

  global_factor = -2.355/1.462*ecc*ecc*pow(v0/v, 19.0/3.0);
  global_factor *= (3.0/128.0/eta);  // overall factor except v^-5 in phase term, this is Newtonian phase term
  if(ecc_order == -1) {
    ecc_order = LAL_MAX_ECC_PN_ORDER;
  }
  if(ecc_order > LAL_MAX_ECC_PN_ORDER) {
    return XLAL_REAL8_FAIL_NAN;
  }

  REAL8 BHPphaseOrder = 0;
  for(int i=7; i<=ecc_order; i++)
  {
    BHPphaseOrder = 0;
    INT4 k = 0;
    for(int j=i; j>=0; j--)
    {
      k = i - j;

      if (i==8)
      {
        if (j==8) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]+(4004693587*log(v))/6.3504e6)*v_power[j]*v0_power[k];}
        else if (j==6) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (2350835165449*log(v))/3.9880512e10)*v_power[j]*v0_power[k];}
        else if (j==2) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (3062544163169*log(v0))/2.0658105216e10)*v_power[j]*v0_power[k];}
        else if (j==0) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (9397061*log(v0))/254016.)*v_power[j]*v0_power[k];}
   
      }

      else if (i==9)
      {
        if (j==9) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]+ (156796692229*LAL_PI*log(v))/1.59522048e8)*v_power[j]*v0_power[k];}
        else if (j==6) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (312836165681*LAL_PI*log(v))/2.848608e9)*v_power[j]*v0_power[k];}
        else if (j==3) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k] + (10259875723*LAL_PI*log(v0))/7.12152e7)*v_power[j]*v0_power[k];}
        else if (j==0) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (1250509*LAL_PI*log(v0))/18144.)*v_power[j]*v0_power[k];}

      }
      
      else if (i==10)
      {
        if (j==10) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (1314393019969101407*log(v))/6.8982438260736e14)*v_power[j]*v0_power[k];}
        else if (j==8) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]+ (3472157632121443*log(v))/5.4556540416e12)*v_power[j]*v0_power[k];}
        else if (j==6) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (41869030336962953*log(v))/4.82394673152e14)*v_power[j]*v0_power[k];}
        else if (j==4) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (202099492979647*log(v0))/7.58048772096e11)*v_power[j]*v0_power[k];}
        else if (j==2) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]- (8676187614257777*log(v0))/4.1646740115456e13)*v_power[j]*v0_power[k];}
        else if (j==0) {BHPphaseOrder += (eccPNBHPCoeffs[i][j][k]+ (95697675707*log(v0))/3.072577536e9)*v_power[j]*v0_power[k];}
   
      }
      else
      {
        BHPphaseOrder += eccPNBHPCoeffs[i][j][k]*v_power[j]*v0_power[k];
        //phasing += eccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
      }
    }
      phasing += BHPphaseOrder;
      //ecc_phase_order[i] = phaseOrder*global_factor;
  }
  //fprintf(stdout, "======== DEBUG for eccentricity ================\n");
  //fprintf(stdout, "eccentricityPhasing_F2 phasing = %g, global_factor = %g, ecc_order = %d, ecc = %g\n", phasing, global_factor, ecc_order, ecc);
  return phasing*global_factor;
}


static INT4 UNUSED
SpinEccentricityPNCoeffs_F2(
        const REAL8 m1, 
        const REAL8 m2, 
        const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
	const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
	const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
	const REAL8 chi2sq, /**< Magnitude of dimensionless spin 2 */
        REAL8 SpinEccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]
        )
{
  const REAL8 mtot = m1 + m2;
  const REAL8 eta = m1*m2/mtot/mtot;
  const REAL8 m1M = m1/mtot;
  const REAL8 m2M = m2/mtot;
  
  INT4 ret = 0;
  memset(SpinEccPNCoeffs, 0x00, (LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*(LAL_MAX_ECC_PN_ORDER+1)*sizeof(REAL8));
  SpinEccPNCoeffs[0][0][0] = 0.0; // lowest order constant term

  SpinEccPNCoeffs[2][2][0] = 0.0; //v^2 term
  SpinEccPNCoeffs[2][1][1] = 0.0; //v*v0 term
  SpinEccPNCoeffs[2][0][2] = 0.0; //v0^2 term

  SpinEccPNCoeffs[3][3][0] = XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m2M, 0)*chi2L; //v^3 term
  SpinEccPNCoeffs[3][0][3] = XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_3PNSOEccCoeff(m2M, 3)*chi2L; //v0^3 term

  SpinEccPNCoeffs[4][4][0] = XLALSimInspiralTaylorF2Phasing_4PNS1S2EccCoeff(eta, 0)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m1M, 0)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m2M, 0)*chi2sq; //v^4 term
  SpinEccPNCoeffs[4][2][2] = 0.0; //v^2*v0^2 term
  SpinEccPNCoeffs[4][0][4] = XLALSimInspiralTaylorF2Phasing_4PNS1S2EccCoeff(eta, 4)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m1M, 4)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_4PNSelf2SEccCoeff(m2M, 4)*chi2sq;  //v0^4 term

  SpinEccPNCoeffs[5][5][0] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 0)*chi2L; //v^5 term
  SpinEccPNCoeffs[5][3][2] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 2)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 2)*chi2L; //v^3*v0^2 term
  SpinEccPNCoeffs[5][2][3] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 3)*chi2L; //v^2*v0^3 term
  SpinEccPNCoeffs[5][0][5] = XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m1M, 5)*chi1L + XLALSimInspiralTaylorF2Phasing_5PNSOEccCoeff(m2M, 5)*chi2L;  //v0^5 term

  SpinEccPNCoeffs[6][6][0] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 0)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 0)*chi2L 
                         + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 0)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 0)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 0)*chi2sq; //v^6 term except log(16*v^2) term
  SpinEccPNCoeffs[6][4][2] = XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 2)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 2)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 2)*chi2sq; //v^4*v0^2 term
  SpinEccPNCoeffs[6][3][3] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 3)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 3)*chi2L 
                         + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 3)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 3)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 3)*chi2sq; //v^3*v0^3 term
  SpinEccPNCoeffs[6][2][4] = XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 4)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 4)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 4)*chi2sq; //v^2*v0^4 term
  SpinEccPNCoeffs[6][0][6] = XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m1M, 6)*chi1L + XLALSimInspiralTaylorF2Phasing_6PNSOEccCoeff(m2M, 6)*chi2L 
                         + XLALSimInspiralTaylorF2Phasing_6PNS1S2OEccCoeff(eta, 6)*chi1L*chi2L + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m1M, 6)*chi1sq
                         + XLALSimInspiralTaylorF2Phasing_6PNSelf2SEccCoeff(m2M, 6)*chi2sq;  //v0^6 term except log(16*v0^2) term
  //printSpinPNCoeffs_F2(SpinEccPNCoeffs);
  return ret;
}


/**
 * Compute eccentric phase correction term using eccPNCeoffs[k][i][j]
 *
 */

// static REAL8 UNUSED
// eccentricityPhasing_F2(REAL8 v, REAL8 v0, REAL8 ecc, REAL8 eta, INT4 ecc_order)
// {
//   static REAL8 v0_power[LAL_MAX_ECC_PN_ORDER+1];
//   /* following code is not efficient in memory usage, need to be improved later */
//   static REAL8 eccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]; // we want to calculate just one time
//   REAL8 v_power[LAL_MAX_ECC_PN_ORDER+1];
//   REAL8 phasing = 0.0;
//   REAL8 global_factor;
//   v0_power[0] = 1.0;
//   for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
//   {
//     v0_power[i] = v0_power[i-1]*v0;
//   }
//   eccentricityPNCoeffs_F2(eta, eccPNCoeffs);
//   //printPNCoeffs_F2(eccPNCoeffs);
//   v_power[0] = 1.0;
//   for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
//   {
//     v_power[i] = v_power[i-1]*v;
//   }

//   global_factor = -2.355/1.462*ecc*ecc*pow(v0/v, 19.0/3.0);
//   global_factor *= (3.0/128.0/eta);  // overall factor except v^-5 in phase term, this is Newtonian phase term

  
//   if(ecc_order == -1) {
//     ecc_order = LAL_MAX_ECC_PN_ORDER;
//   }
//   if(ecc_order > LAL_MAX_ECC_PN_ORDER) {
//     return XLAL_REAL8_FAIL_NAN;
//   }

//   REAL8 phaseOrder = 0;
//   for(int i=0; i<=ecc_order; i++)
//   {
//     phaseOrder = 0;
//     INT4 k = 0;
//     for(int j=i; j>=0; j--)
//     {
//       k = i - j;
//       if( j==6 )
//       {
//         phaseOrder += (eccPNCoeffs[i][j][k]+53.6803271/3.9564000*log(16.0*v_power[2]))*v_power[j]*v0_power[k];
//         //phasing += (eccPNCoeffs[i][j][k]+53.6803271/3.9564000*log(16.0*v_power[2]))*v_power[j]*v0_power[k];
//       }
//       else if( k == 6 )
//       {
//         phaseOrder += (eccPNCoeffs[i][j][k] - 33.17/2.52*log(16.0*v0_power[2]))*v_power[j]*v0_power[k];
//         //phasing += (eccPNCoeffs[i][j][k] - 33.17/2.52*log(16.0*v0_power[2]))*v_power[j]*v0_power[k];
//       }
//       else
//       {
//         phaseOrder += eccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
//         //phasing += eccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
//       }
//     }
//       phasing += phaseOrder;
//       //ecc_phase_order[i] = phaseOrder*global_factor;
//   }
//   //fprintf(stdout, "======== DEBUG for eccentricity ================\n");
//   //fprintf(stdout, "eccentricityPhasing_F2 phasing = %g, global_factor = %g, ecc_order = %d, ecc = %g\n", phasing, global_factor, ecc_order, ecc);
//   return phasing*global_factor;
// }

static REAL8 UNUSED
SpinEccentricityPhasing_F2(REAL8 v, REAL8 v0, REAL8 ecc, REAL8 m1, REAL8 m2, REAL8 chi1L, REAL8 chi2L, REAL8 chi1sq, REAL8 chi2sq, INT4 ecc_order)
{
  static REAL8 v0_power[LAL_MAX_ECC_PN_ORDER+1];
  /* following code is not efficient in memory usage, need to be improved later */
  static REAL8 SpinEccPNCoeffs[LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1][LAL_MAX_ECC_PN_ORDER+1]; // we want to calculate just one time
  REAL8 v_power[LAL_MAX_ECC_PN_ORDER+1];
  REAL8 phasing = 0.0;
  REAL8 global_factor;
  REAL8 eta= (m1*m2)/((m1+m2)*(m1+m2));
  //REAL8 chi1sq = chi1L*chi1L;
  //REAL8 chi2sq = chi2L*chi2L; 
  v0_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v0_power[i] = v0_power[i-1]*v0;
  }
  SpinEccentricityPNCoeffs_F2(m1, m2, chi1L, chi2L, chi1sq, chi2sq, SpinEccPNCoeffs);
  //printPNCoeffs_F2(SpinEccPNCoeffs);
  v_power[0] = 1.0;
  for(int i=1; i<=LAL_MAX_ECC_PN_ORDER; i++)
  {
    v_power[i] = v_power[i-1]*v;
  }

  global_factor = -2.355/1.462*ecc*ecc*pow(v0/v, 19.0/3.0);
  global_factor *= (3.0/128.0/eta);  // overall factor except v^-5 in phase term, this is Newtonian phase term

  
  if(ecc_order == -1) {
    ecc_order = LAL_MAX_ECC_PN_ORDER;
  }
  if(ecc_order > LAL_MAX_ECC_PN_ORDER) {
    return XLAL_REAL8_FAIL_NAN;
  }

  REAL8 SpinPhaseOrder = 0;
  for(int i=0; i<=ecc_order; i++)
  {
    SpinPhaseOrder = 0;
    INT4 k = 0;
    for(int j=i; j>=0; j--)
    {
      k = i - j;
      SpinPhaseOrder += SpinEccPNCoeffs[i][j][k]*v_power[j]*v0_power[k];
    }
    phasing += SpinPhaseOrder;
      //ecc_phase_order[i] = SpinPhaseOrder*global_factor;
  }
  //fprintf(stdout, "======== DEBUG for eccentricity ================\n");
  //fprintf(stdout, "eccentricityPhasing_F2 phasing = %g, global_factor = %g, ecc_order = %d, ecc = %g\n", phasing, global_factor, ecc_order, ecc);
  return phasing*global_factor;
}