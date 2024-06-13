/*
 *  Copyright (C) 2019 HyunWon Lee, JeongCho Kim, Chunglee Kim, Marc Favata, K.G. Arun
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
 *    - LALSimInspiralTaylorF2.c
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
/*
#include <stdlib.h>
#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/Sequence.h>
#include <lal/LALDatatypes.h>
#include <lal/LALSimInspiralEOS.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALError.h>
#include <lal/AVFactories.h>
#include "LALSimInspiralPNCoefficients.c"
*/

#ifndef _OPENMP
#define omp ignore
#endif

#define INCLUDE_SPIN_ECC_PIECE true
#define INCLUDE_BHP_PIECE true

/**
 * @addtogroup LALSimInspiralTaylorF2Ecc_c
 * @brief Routines for generating eccentric TaylorF2 waveforms
 * @{
 *
 * @review TaylorF2Ecc is reviewed, with review statement and git has details available at https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/home.
 *
 * @name Routines for TaylorF2Ecc Waveforms
 * @sa
 * Section IIIF of Alessandra Buonanno, Bala R Iyer, Evan
 * Ochsner, Yi Pan, and B S Sathyaprakash, "Comparison of post-Newtonian
 * templates for compact binary inspiral signals in gravitational-wave
 * detectors", Phys. Rev. D 80, 084043 (2009), arXiv:0907.0700v1
 *
 * Section IV of Marc, et al paper Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304.
 * review page is https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/Eccentric-phase-PN-coefficient-form.
 *
 * @{
 */

/**
 * \author Jeongcho Kim, Chunglee Kim, Hyung Won Lee, Marc Favata, K.G. Arun
 * \file
 *
 * \brief Module to compute the eccentric TaylorF2 inspiral waveform for small eccentricity.
 * Code is based on Section IV of Marc, et al paper Phys. Rev. D 93, 124061 (2016), arXiv:1605.00304.
 * Code review page is https://git.ligo.org/waveforms/reviews/taylorf2ecc/wikis/Eccentric-phase-PN-coefficient-form.
 */

int XLALSimInspiralSpinTaylorF2CoreEcc(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
	const REAL8Sequence *freqs,            /**< frequency points at which to evaluate the waveform (Hz) */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
	const REAL8 shft,		       /**< time shift to be applied to frequency-domain phase (sec)*/
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 eccentricity,
        const REAL8 chi1L, /**< Component of dimensionless spin 1 along Lhat */
        const REAL8 chi2L, /**< Component of dimensionless spin 2 along Lhat */
        //const REAL8 chi1sq,/**< Magnitude of dimensionless spin 1 */
        //const REAL8 chi2sq,  /**< Magnitude of dimensionless spin 2 */              /**< eccentricity effect control < 0 : no eccentricity effect */
        LALDict *p,                       /**< Linked list containing the extra parameters >**/
        PNPhasingSeries *pfaP /**< Phasing coefficients >**/
        )
{
    REAL8 chi1sq = chi1L*chi1L;
    REAL8 chi2sq = chi2L*chi2L;
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (!freqs) XLAL_ERROR(XLAL_EFAULT);
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    REAL8 amp0;
    size_t i;
    COMPLEX16 *data = NULL;
    LIGOTimeGPS tC = {0, 0};
    INT4 iStart = 0;

    COMPLEX16FrequencySeries *htilde = NULL;

    if (*htilde_out) { //case when htilde_out has been allocated in XLALSimInspiralTaylorF2
	    htilde = *htilde_out;
	    iStart = htilde->data->length - freqs->length; //index shift to fill pre-allocated data
	    if(iStart < 0) XLAL_ERROR(XLAL_EFAULT);
    }
    else { //otherwise allocate memory here
	    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, freqs->data[0], 0., &lalStrainUnit, freqs->length);
	    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
	    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);
    }

    /* phasing coefficients */
    PNPhasingSeries pfa = *pfaP;

    REAL8 pfaN = 0.; REAL8 pfa1 = 0.;
    REAL8 pfa2 = 0.; REAL8 pfa3 = 0.; REAL8 pfa4 = 0.;
    REAL8 pfa5 = 0.; REAL8 pfl5 = 0.;
    REAL8 pfa6 = 0.; REAL8 pfl6 = 0.;
    REAL8 pfa7 = 0.;
    REAL8 pfa8 = 0.; REAL8 pfl8 = 0.; REAL8 pflsq8 = 0;
    REAL8 pfa9 = 0.; REAL8 pfl9 = 0.;
    REAL8 pfa10 = 0.; REAL8 pfl10 = 0.;
    REAL8 pfa11 = 0.; REAL8 pfl11 = 0.; 
    REAL8 pfa12 = 0.; REAL8 pfl12 = 0.; REAL8 pflsq12 = 0;
    REAL8 pfa13 = 0.; REAL8 pfl13 = 0.;
    REAL8 pfa14 = 0.; REAL8 pfl14 = 0.; REAL8 pflsq14 = 0;
    REAL8 pfa15 = 0.; REAL8 pfl15 = 0.; REAL8 pflsq15 = 0;
    REAL8 pfa16 = 0.; REAL8 pfl16 = 0.; REAL8 pflsq16 = 0;
    REAL8 pfa17 = 0.; REAL8 pfl17 = 0.; REAL8 pflsq17 = 0;
    REAL8 pfa18 = 0.; REAL8 pfl18 = 0.; REAL8 pflsq18 = 0; REAL8 pflcu18 = 0;
    REAL8 pfa19 = 0.; REAL8 pfl19 = 0.; REAL8 pflsq19 = 0;
    REAL8 pfa20 = 0.; REAL8 pfl20 = 0.; REAL8 pflsq20 = 0; REAL8 pflcu20 = 0;
    REAL8 pfa21 = 0.; REAL8 pfl21 = 0.; REAL8 pflsq21 = 0; REAL8 pflcu21 = 0;
    REAL8 pfa22 = 0.; REAL8 pfl22 = 0.; REAL8 pflsq22 = 0; REAL8 pflcu22 = 0;
    REAL8 pfa23 = 0.; REAL8 pfl23 = 0.; REAL8 pflsq23 = 0; REAL8 pflcu23 = 0;
    REAL8 pfa24 = 0.; REAL8 pfl24 = 0.; REAL8 pflsq24 = 0; REAL8 pflcu24 = 0; REAL8 pflquar24 = 0;



    INT4 phaseO=XLALSimInspiralWaveformParamsLookupPNPhaseOrder(p);
    switch (phaseO)
    {
        case -1:

        case 24:
            pfa24 = pfa.v[24];
            pfl24 = pfa.vlogv[24];
            pflsq24 = pfa.vlogvsq[24];
            pflcu24 = pfa.vlogvcu[24];
            pflquar24 = pfa.vlogvquar[24];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 23:
            pfa23 = pfa.v[23];
            pfl23 = pfa.vlogv[23];
            pflsq23 = pfa.vlogvsq[23];
            pflcu23 = pfa.vlogvcu[23];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 22:
            pfa22 = pfa.v[22];
            pfl22 = pfa.vlogv[22];
            pflsq22 = pfa.vlogvsq[22];
            pflcu22 = pfa.vlogvcu[22];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 21:
            pfa21 = pfa.v[21];
            pfl21 = pfa.vlogv[21];
            pflsq21 = pfa.vlogvsq[21];
            pflcu21 = pfa.vlogvcu[21];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 20:
            pfa20 = pfa.v[20];
            pfl20 = pfa.vlogv[20];
            pflsq20 = pfa.vlogvsq[20];
            pflcu20 = pfa.vlogvcu[20];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 19:
            pfa19 = pfa.v[19];
            pfl19 = pfa.vlogv[19];
            pflsq19 = pfa.vlogvsq[19];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 18:
            pfa18 = pfa.v[18];
            pfl18 = pfa.vlogv[18];
            pflsq18 = pfa.vlogvsq[18];
            pflcu18 = pfa.vlogvcu[18];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 17:
            pfa17 = pfa.v[17];
            pfl17 = pfa.vlogv[17];
            pflsq17 = pfa.vlogvsq[17];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 16:
            pfa16 = pfa.v[16];
            pfl16 = pfa.vlogv[16];
            pflsq16 = pfa.vlogvsq[16];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 15:
            pfa15 = pfa.v[15];
            pfl15 = pfa.vlogv[15];
            pflsq15 = pfa.vlogvsq[15];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 14:
            pfa14 = pfa.v[14];
            pfl14 = pfa.vlogv[14];
            pflsq14 = pfa.vlogvsq[14];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 13:
            pfa13 = pfa.v[13];
            pfl13 = pfa.vlogv[13];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 12:
            pfa12 = pfa.v[12];
            pfl12 = pfa.vlogv[12];
            pflsq12 = pfa.vlogvsq[12];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 11:
            pfa11 = pfa.v[11];
            pfl11 = pfa.vlogv[11];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 10:
            pfa10 = pfa.v[10];
            pfl10 = pfa.vlogv[10];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 9:
            pfa9 = pfa.v[9];
            pfl9 = pfa.vlogv[9];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 8:
            pfa8 = pfa.v[8];
            pfl8 = pfa.vlogv[8];
            pflsq8 = pfa.vlogvsq[8];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif

        case 7:
            pfa7 = pfa.v[7];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 6:
            pfa6 = pfa.v[6];
            pfl6 = pfa.vlogv[6];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 5:
            pfa5 = pfa.v[5];
            pfl5 = pfa.vlogv[5];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 4:
            pfa4 = pfa.v[4];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 3:
            pfa3 = pfa.v[3];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 2:
            pfa2 = pfa.v[2];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 1:
            pfa1 = pfa.v[1];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case 0:
            pfaN = pfa.v[0];
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid phase PN order %d", phaseO);
    }

    /* Validate expansion order arguments.
     * This must be done here instead of in the OpenMP parallel loop
     * because when OpenMP parallelization is turned on, early exits
     * from loops (via return or break statements) are not permitted.
     */

    /* Validate amplitude PN order. */
    INT4 amplitudeO=XLALSimInspiralWaveformParamsLookupPNAmplitudeOrder(p);
    switch (amplitudeO)
    {
        case -1:
        case 7:
        case 6:
        case 5:
        case 4:
        case 3:
        case 2:
        case 0:
            break;
        default:
            XLAL_ERROR(XLAL_ETYPE, "Invalid amplitude PN order %d", amplitudeO);
    }

    /* Generate tidal terms separately.
     * Enums specifying tidal order are in LALSimInspiralWaveformFlags.h
     */
    REAL8 pft10 = 0.;
    REAL8 pft12 = 0.;
    REAL8 pft13 = 0.;
    REAL8 pft14 = 0.;
    REAL8 pft15 = 0.;
    switch( XLALSimInspiralWaveformParamsLookupPNTidalOrder(p) )
    {
	case LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL:
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_75PN:
            pft15 = pfa.v[15];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_7PN:
            pft14 = pfa.v[14];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_65PN:
            pft13 = pfa.v[13];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
            __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_6PN:
	    pft12 = pfa.v[12];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_5PN:
            pft10 = pfa.v[10];
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
        case LAL_SIM_INSPIRAL_TIDAL_ORDER_0PN:
            break;
        default:
            XLAL_ERROR(XLAL_EINVAL, "Invalid tidal PN order %d", XLALSimInspiralWaveformParamsLookupPNTidalOrder(p));
    }

    /* The flux and energy coefficients below are used to compute SPA amplitude corrections */

    /* flux coefficients */
    const REAL8 FTaN = XLALSimInspiralPNFlux_0PNCoeff(eta);
    const REAL8 FTa2 = XLALSimInspiralPNFlux_2PNCoeff(eta);
    const REAL8 FTa3 = XLALSimInspiralPNFlux_3PNCoeff(eta);
    const REAL8 FTa4 = XLALSimInspiralPNFlux_4PNCoeff(eta);
    const REAL8 FTa5 = XLALSimInspiralPNFlux_5PNCoeff(eta);
    const REAL8 FTl6 = XLALSimInspiralPNFlux_6PNLogCoeff(eta);
    const REAL8 FTa6 = XLALSimInspiralPNFlux_6PNCoeff(eta);
    const REAL8 FTa7 = XLALSimInspiralPNFlux_7PNCoeff(eta);

    /* energy coefficients */
    const REAL8 dETaN = 2. * XLALSimInspiralPNEnergy_0PNCoeff(eta);
    const REAL8 dETa1 = 2. * XLALSimInspiralPNEnergy_2PNCoeff(eta);
    const REAL8 dETa2 = 3. * XLALSimInspiralPNEnergy_4PNCoeff(eta);
    const REAL8 dETa3 = 4. * XLALSimInspiralPNEnergy_6PNCoeff(eta);


    /* Perform some initial checks */
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);

    /* extrinsic parameters */
    amp0 = -4. * m1 * m2 / r * LAL_MRSUN_SI * LAL_MTSUN_SI * sqrt(LAL_PI/12.L);

    data = htilde->data->data;

    REAL8 v_ecc_ref = 0.0;
    REAL8 f_ecc = XLALSimInspiralWaveformParamsLookupEccentricityFreq(p);
    if( eccentricity > 0) {
        v_ecc_ref = cbrt(piM*f_ecc);
    }

    /* Compute the SPA phase at the reference point
     * N.B. f_ref == 0 means we define the reference time/phase at "coalescence"
     * when the frequency approaches infinity. In that case,
     * the integrals Eq. 3.15 of arXiv:0907.0700 vanish when evaluated at
     * f_ref == infinity. If f_ref is finite, we must compute the SPA phase
     * evaluated at f_ref, store it as ref_phasing and subtract it off.
     */
    REAL8 ref_phasing = 0.;
    INT4 ecc_order = XLALSimInspiralWaveformParamsLookupPNEccentricityOrder(p);
    if( f_ref != 0. ) {
        const REAL8 vref = cbrt(piM*f_ref);
        const REAL8 logvref = log(vref);
        const REAL8 logvrefsq = logvref * logvref;
        const REAL8 logvrefcu = logvrefsq * logvref;
        const REAL8 logvrefquar = logvrefcu * logvref;
        const REAL8 v2ref = vref * vref;
        const REAL8 v3ref = vref * v2ref;
        const REAL8 v4ref = vref * v3ref;
        const REAL8 v5ref = vref * v4ref;
        const REAL8 v6ref = vref * v5ref;
        const REAL8 v7ref = vref * v6ref;
        const REAL8 v8ref = vref * v7ref;
        const REAL8 v9ref = vref * v8ref;
        const REAL8 v10ref = vref * v9ref;
        const REAL8 v11ref = vref * v10ref;
        const REAL8 v12ref = v2ref * v10ref;
        const REAL8 v13ref = vref * v12ref;
        const REAL8 v14ref = vref * v13ref;
        const REAL8 v15ref = vref * v14ref;
        /*BHP orders*/
        const REAL8 v16ref = vref * v15ref;
        const REAL8 v17ref = vref * v16ref;
        const REAL8 v18ref = vref * v17ref;
        const REAL8 v19ref = vref * v18ref;
        const REAL8 v20ref = vref * v19ref;
        const REAL8 v21ref = vref * v20ref;
        const REAL8 v22ref = vref * v21ref;
        const REAL8 v23ref = vref * v22ref;
        const REAL8 v24ref = vref * v23ref;        

        ref_phasing += pfa7 * v7ref;
        ref_phasing += (pfa6 + pfl6 * logvref) * v6ref;
        ref_phasing += (pfa5 + pfl5 * logvref) * v5ref;
        ref_phasing += pfa4 * v4ref;
        ref_phasing += pfa3 * v3ref;
        ref_phasing += pfa2 * v2ref;
        ref_phasing += pfa1 * vref;
        ref_phasing += pfaN;

        /* Tidal terms in reference phasing */
        ref_phasing += pft15 * v15ref;
        ref_phasing += pft14 * v14ref;
        ref_phasing += pft13 * v13ref;
        ref_phasing += pft12 * v12ref;
        ref_phasing += pft10 * v10ref;

        /*Circular BHP terms in reference phasing*/
        ref_phasing += (pfa24 + pfl24 * logvref + pflsq24 * logvrefsq + pflcu24 * logvrefcu + pflquar24 * logvrefquar) * v24ref;
        ref_phasing += (pfa23 + pfl23 * logvref + pflsq23 * logvrefsq + pflcu23 * logvrefcu) * v23ref;
        ref_phasing += (pfa22 + pfl22 * logvref + pflsq22 * logvrefsq + pflcu22 * logvrefcu) * v22ref;
        ref_phasing += (pfa21 + pfl21 * logvref + pflsq21 * logvrefsq + pflcu21 * logvrefcu) * v21ref;
        ref_phasing += (pfa20 + pfl20 * logvref + pflsq20 * logvrefsq + pflcu20 * logvrefcu) * v20ref;
        ref_phasing += (pfa19 + pfl19 * logvref + pflsq19 * logvrefsq) * v19ref;
        ref_phasing += (pfa18 + pfl18 * logvref + pflsq18 * logvrefsq + pflcu18 * logvrefcu) * v18ref;
        ref_phasing += (pfa17 + pfl17 * logvref + pflsq17 * logvrefsq) * v17ref;
        ref_phasing += (pfa16 + pfl16 * logvref + pflsq16 * logvrefsq) * v16ref;
        ref_phasing += (pfa15 + pfl15 * logvref + pflsq15 * logvrefsq) * v15ref;
        ref_phasing += (pfa14 + pfl14 * logvref + pflsq14 * logvrefsq) * v14ref;
        ref_phasing += (pfa13 + pfl13 * logvref) * v13ref;
        ref_phasing += (pfa12 + pfl12 * logvref + pflsq12 * logvrefsq) * v12ref;
        ref_phasing += (pfa11 + pfl11 * logvref) * v11ref;
        ref_phasing += (pfa10 + pfl10 * logvref) * v10ref;
        ref_phasing += (pfa9 + pfl9 * logvref) * v9ref;
        ref_phasing += (pfa8 + pfl8 * logvref + pflsq8 * logvrefsq) * v8ref;

        /* Eccentricity terms in phasing */
        if( eccentricity > 0 ) {
          ref_phasing += eccentricityPhasing_F2(vref, v_ecc_ref, eccentricity, eta, ecc_order);
        }

        if (INCLUDE_SPIN_ECC_PIECE){
        if(eccentricity > 0 && (fabs(chi1L) > 0 || fabs(chi2L) > 0)) {
          ref_phasing += SpinEccentricityPhasing_F2(vref, v_ecc_ref, eccentricity, m1, m2, chi1L, chi2L, chi1sq, chi2sq, ecc_order);
        }
        }
        
        /* Eccentric BHP terms in phasing */
        if (INCLUDE_BHP_PIECE){        
        if( eccentricity > 0 ) {
          ref_phasing += eccentricityPhasingBHP_F2(vref, v_ecc_ref, eccentricity, eta, ecc_order);
        }
        }

        ref_phasing /= v5ref;
    } /* End of if(f_ref != 0) block */
    #pragma omp parallel for
    for (i = 0; i < freqs->length; i++) {
        const REAL8 f = freqs->data[i];
        const REAL8 v = cbrt(piM*f);
        const REAL8 logv = log(v);
        const REAL8 logvsq = logv * logv;
        const REAL8 logvcu = logvsq * logv;
        const REAL8 logvquar = logvcu * logv;
        const REAL8 v2 = v * v;
        const REAL8 v3 = v * v2;
        const REAL8 v4 = v * v3;
        const REAL8 v5 = v * v4;
        const REAL8 v6 = v * v5;
        const REAL8 v7 = v * v6;
        const REAL8 v8 = v * v7;
        const REAL8 v9 = v * v8;
        const REAL8 v10 = v * v9;
        const REAL8 v11 = v * v10;
        const REAL8 v12 = v2 * v10;
        const REAL8 v13 = v * v12;
        const REAL8 v14 = v * v13;
        const REAL8 v15 = v * v14;
        /*BHP orders*/
        const REAL8 v16 = v * v15;
        const REAL8 v17 = v * v16;
        const REAL8 v18 = v * v17;
        const REAL8 v19 = v * v18;
        const REAL8 v20 = v * v19;
        const REAL8 v21 = v * v20;
        const REAL8 v22 = v * v21;
        const REAL8 v23 = v * v22;
        const REAL8 v24 = v * v23;

        REAL8 phasing = 0.;
        REAL8 dEnergy = 0.;
        REAL8 flux = 0.;
        REAL8 amp;

        phasing += pfa7 * v7;
        phasing += (pfa6 + pfl6 * logv) * v6;
        phasing += (pfa5 + pfl5 * logv) * v5;
        phasing += pfa4 * v4;
        phasing += pfa3 * v3;
        phasing += pfa2 * v2;
        phasing += pfa1 * v;
        phasing += pfaN;

        /* Tidal terms in phasing */
        phasing += pft15 * v15;
        phasing += pft14 * v14;
        phasing += pft13 * v13;
        phasing += pft12 * v12;
        phasing += pft10 * v10;

        /*Circular BHP terms in phasing*/
        phasing += (pfa24 + pfl24 * logv + pflsq24 * logvsq + pflcu24 * logvcu + pflquar24 * logvquar) * v24;
        phasing += (pfa23 + pfl23 * logv + pflsq23 * logvsq + pflcu23 * logvcu) * v23;
        phasing += (pfa22 + pfl22 * logv + pflsq22 * logvsq + pflcu22 * logvcu) * v22;
        phasing += (pfa21 + pfl21 * logv + pflsq21 * logvsq + pflcu21 * logvcu) * v21;
        phasing += (pfa20 + pfl20 * logv + pflsq20 * logvsq + pflcu20 * logvcu) * v20;
        phasing += (pfa19 + pfl19 * logv + pflsq19 * logvsq) * v19;
        phasing += (pfa18 + pfl18 * logv + pflsq18 * logvsq + pflcu18 * logvcu) * v18;
        phasing += (pfa17 + pfl17 * logv + pflsq17 * logvsq) * v17;
        phasing += (pfa16 + pfl16 * logv + pflsq16 * logvsq) * v16;
        phasing += (pfa15 + pfl15 * logv + pflsq15 * logvsq) * v15;
        phasing += (pfa14 + pfl14 * logv + pflsq14 * logvsq) * v14;
        phasing += (pfa13 + pfl13 * logv) * v13;
        phasing += (pfa12 + pfl12 * logv + pflsq12 * logvsq) * v12;
        phasing += (pfa11 + pfl11 * logv) * v11;
        phasing += (pfa10 + pfl10 * logv) * v10;
        phasing += (pfa9 + pfl9 * logv) * v9;
        phasing += (pfa8 + pfl8 * logv + pflsq8 * logvsq) * v8;

        /* Eccentricity terms in phasing */
        if( eccentricity > 0 ) {
          phasing += eccentricityPhasing_F2(v, v_ecc_ref, eccentricity, eta, ecc_order);
        }
    
        if(INCLUDE_SPIN_ECC_PIECE){
        if(eccentricity > 0 && (fabs(chi1L) > 0 || fabs(chi2L) > 0)) {
          phasing += SpinEccentricityPhasing_F2(v, v_ecc_ref, eccentricity, m1, m2, chi1L, chi2L, chi1sq, chi2sq, ecc_order);
        }
        }

        /* Eccentric BHP terms in phasing */
        if(INCLUDE_BHP_PIECE){    
        if( eccentricity > 0 ) {
          phasing += eccentricityPhasingBHP_F2(v, v_ecc_ref, eccentricity, eta, ecc_order);
        }
        }

        phasing /= v5;

    /* WARNING! Amplitude orders beyond 0 have NOT been reviewed!
     * Use at your own risk. The default is to turn them off.
     * These do not currently include spin corrections.
     * Note that these are not higher PN corrections to the amplitude.
     * They are the corrections to the leading-order amplitude arising
     * from the stationary phase approximation. See for instance
     * Eq 6.9 of arXiv:0810.5336
     */
	switch (amplitudeO)
        {
            case 7:
                flux += FTa7 * v7;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 6:
                flux += (FTa6 + FTl6*logv) * v6;
                dEnergy += dETa3 * v6;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 5:
                flux += FTa5 * v5;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 4:
                flux += FTa4 * v4;
                dEnergy += dETa2 * v4;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 3:
                flux += FTa3 * v3;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case 2:
                flux += FTa2 * v2;
                dEnergy += dETa1 * v2;
#if __GNUC__ >= 7 && !defined __INTEL_COMPILER
                __attribute__ ((fallthrough));
#endif
            case -1: /* Default to no SPA amplitude corrections */
            case 0:
                flux += 1.;
                dEnergy += 1.;
        }

        flux *= FTaN * v10;
        dEnergy *= dETaN * v;
        // Note the factor of 2 b/c phi_ref is orbital phase
        phasing += shft * f - 2.*phi_ref - ref_phasing;
        amp = amp0 * sqrt(-dEnergy/flux) * v;
        data[i+iStart] = amp * cos(phasing - LAL_PI_4)
                - amp * sin(phasing - LAL_PI_4) * 1.0j;
    }

    *htilde_out = htilde;
    return XLAL_SUCCESS;
}

/**
 * Computes the stationary phase approximation to the Fourier transform of
 * a chirp waveform with eccentric correction. The amplitude is given by expanding \f$1/\sqrt{\dot{F}}\f$.
 * If the PN order is set to -1, then the highest implemented order is used.
 *
 * @note f_ref is the GW frequency at which phi_ref is defined. The most common
 * choice in the literature is to choose the reference point as "coalescence",
 * when the frequency becomes infinite. This is the behavior of the code when
 * f_ref==0. If f_ref > 0, phi_ref sets the orbital phase at that GW frequency.
 *
 * See arXiv:0810.5336 and arXiv:astro-ph/0504538 for spin corrections
 * to the phasing.
 * See arXiv:1303.7412 for spin-orbit phasing corrections at 3 and 3.5PN order
 * See Phys. Rev. Lett. 112, 101101(2014) for eccentric phasing corrections upto 3PN order
 *
 * The spin and tidal order enums are defined in LALSimInspiralWaveformFlags.h
 */
int XLALSimInspiralSpinTaylorF2Ecc(
        COMPLEX16FrequencySeries **htilde_out, /**< FD waveform */
        const REAL8 phi_ref,                   /**< reference orbital phase (rad) */
        const REAL8 deltaF,                    /**< frequency resolution */
        const REAL8 m1_SI,                     /**< mass of companion 1 (kg) */
        const REAL8 m2_SI,                     /**< mass of companion 2 (kg) */
        const REAL8 S1z,                       /**<  z component of the spin of companion 1 */
        const REAL8 S2z,                       /**<  z component of the spin of companion 2  */
        const REAL8 fStart,                    /**< start GW frequency (Hz) */
        const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation - if 0, end at Schwarzschild ISCO */
        const REAL8 f_ref,                     /**< Reference GW frequency (Hz) - if 0 reference point is coalescence */
        const REAL8 r,                         /**< distance of source (m) */
        const REAL8 eccentricity,                       /**< eccentricity effect control < 0 : no eccentricity effect */
        LALDict *p                             /**< Linked list containing the extra parameters >**/
        )
{
    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_SI / LAL_MSUN_SI;
    const REAL8 m2 = m2_SI / LAL_MSUN_SI;
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    // const REAL8 eta = m1 * m2 / (m * m);
    const REAL8 piM = LAL_PI * m_sec;
    const REAL8 vISCO = 1. / sqrt(6.);
    const REAL8 fISCO = vISCO * vISCO * vISCO / piM;
    //const REAL8 m1OverM = m1 / m;
    // const REAL8 m2OverM = m2 / m;
    REAL8 shft, f_max;
    size_t i, n;
    INT4 iStart;
    REAL8Sequence *freqs = NULL;
    LIGOTimeGPS tC = {0, 0};
    int ret;
    int retcode;
    REAL8 fCONT;
    INT4 tideO = XLALSimInspiralWaveformParamsLookupPNTidalOrder(p);
    REAL8 lambda1 = XLALSimInspiralWaveformParamsLookupTidalLambda1(p);
    REAL8 lambda2 = XLALSimInspiralWaveformParamsLookupTidalLambda2(p);
    retcode = XLALSimInspiralSetQuadMonParamsFromLambdas(p);
    XLAL_CHECK(retcode == XLAL_SUCCESS, XLAL_EFUNC, "Failed to set quadparams from Universal relation.\n");

    COMPLEX16FrequencySeries *htilde = NULL;

    /* Perform some initial checks */
    if (!htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (*htilde_out) XLAL_ERROR(XLAL_EFAULT);
    if (m1_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (m2_SI <= 0) XLAL_ERROR(XLAL_EDOM);
    if (fStart <= 0) XLAL_ERROR(XLAL_EDOM);
    if (f_ref < 0) XLAL_ERROR(XLAL_EDOM);
    if (r <= 0) XLAL_ERROR(XLAL_EDOM);
    if (eccentricity < 0.0 || eccentricity >= 1.0) XLAL_ERROR(XLAL_EDOM);

    /* allocate htilde */
    if ( (fEnd == 0.) && ( tideO == 0)) // End at ISCO
        f_max = fISCO;
    else if ( (fEnd == 0.) && ( tideO != 0 )) { // End at the minimum of the contact and ISCO frequencies only when tides are enabled
        fCONT = XLALSimInspiralContactFrequency(m1, lambda1, m2, lambda2); /* Contact frequency of two compact objects */
        f_max = (fCONT > fISCO) ? fISCO : fCONT;
    }
    else // End at user-specified freq.
        f_max = fEnd;
    if (f_max <= fStart) XLAL_ERROR(XLAL_EDOM);

    n = (size_t) (f_max / deltaF + 1);
    XLALGPSAdd(&tC, -1 / deltaF);  /* coalesce at t=0 */
    htilde = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, n);
    if (!htilde) XLAL_ERROR(XLAL_EFUNC);
    memset(htilde->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&htilde->sampleUnits, &htilde->sampleUnits, &lalSecondUnit);

    /* Fill with non-zero vals from fStart to f_max */
    iStart = (INT4) ceil(fStart / deltaF);

    /* Sequence of frequencies where waveform model is to be evaluated */
    freqs = XLALCreateREAL8Sequence(n - iStart);

    /* extrinsic parameters */
    shft = LAL_TWOPI * (tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);

    #pragma omp parallel for
    for (i = iStart; i < n; i++) {
        freqs->data[i-iStart] = i * deltaF;
    }

    /* phasing coefficients */
    PNPhasingSeries pfa;
    XLALSimInspiralPNBHPPhasing_F2(&pfa, m1, m2, S1z, S2z, S1z*S1z, S2z*S2z, S1z*S2z, p);
    ret = XLALSimInspiralSpinTaylorF2CoreEcc(&htilde, freqs, phi_ref, m1_SI, m2_SI,
                                      f_ref, shft, r, eccentricity, S1z, S2z, p, &pfa);

    XLALDestroyREAL8Sequence(freqs);

    *htilde_out = htilde;

    return ret;
}

/** @} */
/** @} */
