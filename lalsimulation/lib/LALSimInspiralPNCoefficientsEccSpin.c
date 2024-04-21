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