/* Copyright (C) 2018 Roland Haas, Eliu A. Huerta, Prayush Kumar, Alvin Chua,
 *                    Christopher J. Moore, Daniel George, Eamonn O'Shea,
 *                    Erik Wessel, Joseph Adamo, Sarah Habib, Aditya Vijaykumar
 *
 *  ENIGMA waveform model:
 * - EA Huerta et al., Phys. Rev. D 97, 024031, https://arxiv.org/abs/1711.06276
 * - EA Huerta et al., Phys. Rev. D 100, 064003,
 * https://arxiv.org/abs/1901.07038
 * - EA Huerta et al., Phys. Rev. D 95, 024038, https://arxiv.org/abs/1609.05933
 *
 *
 */

#ifdef __GNUC__
#define UNUSED __attribute__((unused))
#else
#define UNUSED
#endif

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/FileIO.h>
#include <lal/LALConstants.h>
#include <lal/LALSimIMR.h>
#include <lal/LALSimInspiral.h>
#include <lal/Units.h>
#include <lal/XLALGSL.h>

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_vector.h>

/***********************************************************************************/
/****************************** Type definitions
 * ***********************************/
/***********************************************************************************/

#define L_MIN 2
#define L_MAX 8
#define ONLY_LeqM_MODES false

#ifndef ENIGMADEBUG
#define DEBUG 0
#else
#define DEBUG 1
#endif

// Structure to hold plus/cross components of waveform
typedef struct {
  REAL8 *hp;
  REAL8 *hc;
  int length;
} wave;

// Structure to hold training set
typedef struct {
  REAL8 *time; // array to store the time sampling values
  int length;  // length of time

  wave *Waves;       // array of training set waveforms
  REAL8 *xvals;      // array of mass ratios of the training set points
  int N_Dset_points; // length of Waves, and xvals

  int cov_index_amp;   // which covariance function to use for amplitude
                       // interpolation
  REAL8 *inv_cov_amp;  // inverse covariance matrix for amplitude interpolation,
                       // size N_Dset_points by N_Dset_points
  REAL8 Log10DeltaAmp; // Log10 of the lengthscale for amplitude interpolation

  int cov_index_phase;   // which covariance function to use for phase
                         // interpolation
  REAL8 *inv_cov_phase;  // inverse covariance matrix for phase interpolation,
                         // size N_Dset_points by N_Dset_points
  REAL8 Log10DeltaPhase; // Log10 of the lengthscale for phase interpolation

  int Loaded; // Set this to 0 when you first declare a TrainingSet, the
              // function setup() then sets this to 1

  REAL8 *const_phase_array; // This will be a constant array to store the phase
                            // of the equal mass non-spinning binary in the
                            // training set, which is used to zero the phases
                            // for all other waveforms.

} TrainingSet;

struct kepler_vars {
  int pn_order;
  REAL8 eta;
  REAL8 x;
  REAL8 e;
  REAL8 l;
  REAL8 r;
  REAL8 rDOT;
  REAL8 PhiDOT;

  // store commonly used powers of orbital elements
  REAL8 x2, x3, x4, x5, x6, x7, x8, x9;
  REAL8 e2, e3, e4, e5, e6, e7, e8, e9;
  REAL8 r2, r3, r4, r5, r6, r7, r8, r9;
  REAL8 rDOT2, rDOT3, rDOT4, rDOT5, rDOT6, rDOT7, rDOT8, rDOT9;
  REAL8 PhiDOT2, PhiDOT3, PhiDOT4, PhiDOT5, PhiDOT6, PhiDOT7, PhiDOT8, PhiDOT9;
  REAL8 eta2, eta3, eta4, eta5, eta6;
};

struct ode_parameters {
  REAL8 eta;
  REAL8 m1;
  REAL8 m2;
  REAL8 S1z;
  REAL8 S2z;
  int radiation_pn_order;
  struct kepler_vars *orbital_vars;
};

typedef struct {
  REAL8 *tsamples;
  REAL8 *hp_merger, *hc_merger;
  gsl_interp *interpHp, *interpHc;
  gsl_interp_accel *accHp, *accHc;
  REAL8 omega_attach;
} interp_params;

typedef int (*eval_e_func)(const gsl_interp *interp, const REAL8 xa[],
                           const REAL8 ya[], REAL8 x, gsl_interp_accel *acc,
                           REAL8 *y);

/***********************************************************************************/
/************************** Static function declartions
 * ****************************/
/***********************************************************************************/

/* small powers */
static REAL8 pow2(const REAL8 x);

static REAL8 pow3(const REAL8 x);

static REAL8 pow4(const REAL8 x);

static REAL8 pow5(const REAL8 x);

static REAL8 pow6(const REAL8 x);

static REAL8 pow7(const REAL8 x);

/* fractional powers */
static REAL8 pow1_3(const REAL8 x);

// static REAL8 pow3_2(const REAL8 x);

// static REAL8 pow5_2(const REAL8 x);

static REAL8 pow7_2(const REAL8 x);

/* Helpers */
static REAL8 SymMassRatio(REAL8 q);

static REAL8 SmallMassRatio(REAL8 eta);

static REAL8 SymMassRatio(REAL8 q);

static REAL8 SmallMassRatio(REAL8 eta);

static inline void *my_realloc(void *p, size_t sz);

static TrainingSet *CreateEmptyTrainingSet(void);
#if 0 // keep around for multi-threading
static
void ClearTrainingSet(TrainingSet *Dset);
#endif

static int setup(TrainingSet *Dset);

static void WaveToAmpPhase(REAL8 *T, REAL8 *A, REAL8 *B, int length,
                           REAL8 *const_phase_array);

static void AmpPhaseToWave(REAL8 *A, REAL8 *B, int length,
                           REAL8 *const_phase_array);

static int PN_Omega(REAL8 mr, REAL8 tm, int *pn, REAL8 *w1);

/* Keplerian terms */
static void PopulateKeplerParams(struct kepler_vars *params, const REAL8 e,
                                 const REAL8 x, const REAL8 r, const REAL8 rDOT,
                                 const REAL8 PhiDOT);

static REAL8 cosu_factor(REAL8 e, REAL8 u);

static REAL8 pn_kepler_equation(REAL8 eta, REAL8 x, REAL8 e, REAL8 l);
static REAL8 mikkola_finder(REAL8 e, REAL8 l);

static REAL8 separation(REAL8 u, REAL8 eta, REAL8 x, REAL8 e, REAL8 m1,
                        REAL8 m2, REAL8 S1z, REAL8 S2z);

/* PN expansions */
#define c1 (2.0154525)
#define c2 (-39.57186)
#define c3 (-24.30744)
#define c4 (103.93432)
#define c2log (6.4)
#define c3log (-80.254774)
#define d1 (-7.82849889)
#define d2 (20.84938506)
#define d3 (-20.5092515)
#define d4 (5.383192)
#define cz0 (10.22474)
#define mode_pn_order (8)
#define Radiation_PN_Order (8)

static REAL8 x_dot_0pn(REAL8 e, REAL8 eta);
static REAL8 x_dot_1pn(REAL8 e, REAL8 eta);
static REAL8 x_dot_1_5_pn(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z);
static REAL8 x_dot_hereditary_1_5(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 x_dot_2pn(REAL8 e, REAL8 eta);
static REAL8 x_dot_2pn_SS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z);
static REAL8 x_dot_2_5_pn(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                          REAL8 S2z);
static REAL8 x_dot_hereditary_2_5(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 x_dot_hereditary_3(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 x_dot_3pn(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 x_dot_3pnSO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z);
static REAL8 x_dot_3pnSS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z);
static REAL8 x_dot_3_5_pn(REAL8 e, REAL8 eta);
static REAL8 x_dot_3_5pnSO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                           REAL8 S2z);
static REAL8 x_dot_3_5pn_SS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                            REAL8 S2z);
static REAL8 x_dot_3_5pn_cubicSpin(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2,
                                   REAL8 S1z, REAL8 S2z);
static REAL8 x_dot_4pn(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 x_dot_4pnSO(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z);
static REAL8 x_dot_4pnSS(REAL8 e, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                         REAL8 S2z);
static REAL8 x_dot_4_5_pn(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 dxdt_4pn(REAL8 x, REAL8 eta);
static REAL8 dxdt_4_5pn(REAL8 x, REAL8 eta);
static REAL8 dxdt_5pn(REAL8 x, REAL8 eta);
static REAL8 dxdt_5_5pn(REAL8 x, REAL8 eta);
static REAL8 dxdt_6pn(REAL8 x, REAL8 eta);

static REAL8 e_dot_0pn(REAL8 e, REAL8 eta);
static REAL8 e_dot_1pn(REAL8 e, REAL8 eta);
static REAL8 e_rad_hereditary_1_5(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 e_dot_2pn(REAL8 e, REAL8 eta);
static REAL8 e_rad_hereditary_2_5(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 e_rad_hereditary_3(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 e_dot_3pn(REAL8 e, REAL8 eta, REAL8 x);
static REAL8 e_dot_3_5pn(REAL8 e, REAL8 eta);
static REAL8 e_dot_1_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);
static REAL8 e_dot_2pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);

static REAL8 l_dot_1pn(REAL8 e, REAL8 eta);
static REAL8 l_dot_2pn(REAL8 e, REAL8 eta);
static REAL8 l_dot_3pn(REAL8 e, REAL8 eta);
static REAL8 l_dot_1_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);
static REAL8 l_dot_2pn_SS(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);

static REAL8 phi_dot_0pn(REAL8 e, REAL8 eta, REAL8 u);
static REAL8 phi_dot_1pn(REAL8 e, REAL8 eta, REAL8 u);
static REAL8 phi_dot_2pn(REAL8 e, REAL8 eta, REAL8 u);
static REAL8 phi_dot_3pn(REAL8 e, REAL8 eta, REAL8 u);
//static REAL8 phi_dot_3_5pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);
static REAL8 phi_dot_4pn_SO(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z);
static REAL8 phi_dot_1_5_pnSO_ecc(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
                                  REAL8 S2z, REAL8 u);
static REAL8 phi_dot_2_pnSS_ecc(REAL8 e, REAL8 m1, REAL8 m2, REAL8 S1z,
                                REAL8 S2z, REAL8 u);

/* PN evolution equations */
static REAL8 dx_dt(int radiation_pn_order, REAL8 eta, REAL8 m1, REAL8 m2,
                   REAL8 S1z, REAL8 S2z, REAL8 x, REAL8 e);
static REAL8 de_dt(int radiation_pn_order, REAL8 eta, REAL8 m1, REAL8 m2,
                   REAL8 S1z, REAL8 S2z, REAL8 x, REAL8 e);
static REAL8 dl_dt(REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z, REAL8 S2z, REAL8 x,
                   REAL8 e);
static REAL8 dphi_dt(REAL8 u, REAL8 eta, REAL8 m1, REAL8 m2, REAL8 S1z,
                     REAL8 S2z, REAL8 x, REAL8 e);

static int eccentric_x_model_odes(REAL8 t, const REAL8 y[], REAL8 dydt[],
                                  void *params);

/* PN radiation */
// static REAL8 hPlus(REAL8 x, REAL8 x0, REAL8 m1, REAL8 m2, REAL8 i, REAL8 phi,
// UINT4 vpn); static REAL8 hCross(REAL8 x, REAL8 x0, REAL8 m1, REAL8 m2, REAL8
// i, REAL8 phi, UINT4 vpn);
static REAL8 phi_e(REAL8 e);
static REAL8 psi_e(REAL8 e);
static REAL8 zed_e(REAL8 e);
static REAL8 kappa_e(REAL8 e);
static REAL8 phi_e_tilde(REAL8 e);
static REAL8 psi_e_tilde(REAL8 e);
static REAL8 zed_e_tilde(REAL8 e);
static REAL8 kappa_e_tilde(REAL8 e);
static REAL8 phi_e_rad(REAL8 e);
static REAL8 psi_e_rad(REAL8 e);
static REAL8 zed_e_rad(REAL8 e);
static REAL8 kappa_e_rad(REAL8 e);
static REAL8 f_e(REAL8 e);
static REAL8 capital_f_e(REAL8 e);
static REAL8 psi_n(REAL8 e);
static REAL8 zed_n(REAL8 e);

/* Gaussian process emulator */
static REAL8 GPR_CovarianceFunction_SE(REAL8 x1, REAL8 x2, REAL8 *theta,
                                       REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_SE_BC(REAL8 x1, REAL8 x2, REAL8 *theta,
                                          REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland0(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland0_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland1(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland1_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland2(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland2_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland3(REAL8 x1, REAL8 x2, REAL8 *theta,
                                              REAL8 Jitter);

static REAL8 GPR_CovarianceFunction_Wendland3_BC(REAL8 x1, REAL8 x2,
                                                 REAL8 *theta, REAL8 Jitter);

static REAL8 Covariance(REAL8 x1, REAL8 x2, REAL8 log10Delta, REAL8 Jitter,
                        int covindex);

static void GPR_CovarianceMatrix(REAL8 **covMatrix, REAL8 *xvals, int Npoints,
                                 REAL8 log10Delta, int covindex);

static void GPR_InverseCovarianceMatrix(REAL8 **covMatrix, REAL8 *xvals,
                                        int Npoints, REAL8 log10Delta,
                                        int covindex);

static int interp(REAL8 eta, REAL8 **hp, REAL8 **hc, REAL8 **tsamples,
                  int *len_tsamples, TrainingSet *Dset);

static int LoadSurrogateWaveformFile(REAL8 **hp, REAL8 **hc, LALH5File *file,
                                     const char *dsetname);

/* inspiral and ringdown functions */
static int
x_model_eccbbh_imr_waveform(REAL8TimeSeries *h_plus, REAL8TimeSeries *h_cross,
                            REAL8 mass1,     /* mass1 in solar masses  */
                            REAL8 mass2,     /* mass2 in solar masses  */
                            REAL8 S1z,       /* z-ccomponent of spin of mass1*/
                            REAL8 S2z,       /* z-ccomponent of spin of mass2*/
                            REAL8 e_init,    /* initial eccentricity   */
                            REAL8 f_gw_init, /* initial GW frequency   */
                            REAL8 distance,  /* distance of source (m) */
                            REAL8 mean_anom_init, /* initial mean-anomaly   */
                            REAL8 ode_eps,        /* tolerance (relative)   */
                            REAL8 euler_iota,     /* source polar angle     */
                            REAL8 euler_beta,     /* source azimuthal angle */
                            REAL8 sampling_rate,  /* sample rate in Hz      */
                            TrainingSet *Dset);

static int x_model_eccbbh_inspiral_waveform(
    REAL8TimeSeries *h_plus, REAL8TimeSeries *h_cross,
    REAL8 mass1,          /* mass1 in solar masses  */
    REAL8 mass2,          /* mass2 in solar masses  */
    REAL8 S1z,            /* z-component of the spin of companion 1 */
    REAL8 S2z,            /* z-component of the spin of companion 2 */
    REAL8 e_init,         /* initial eccentricity   */
    REAL8 f_gw_init,      /* initial GW frequency   */
    REAL8 distance,       /* distance of source (m) */
    REAL8 mean_anom_init, /* initial mean-anomaly   */
    REAL8 ode_eps,        /* tolerance (relative)   */
    REAL8 euler_iota,     /* source polar angle     */
    REAL8 euler_beta,     /* source azimuthal angle */
    REAL8 sampling_rate   /* sample rate in Hz      */
);

static int Attach_GPE_Merger_Ringdown(REAL8 dt, REAL8 **h_plus, REAL8 **h_cross,
                                      REAL8 inspiral_matching_time,
                                      REAL8 matching_Hp, REAL8 matching_Hc,
                                      int *Length, TrainingSet *Dset, REAL8 m1,
                                      REAL8 m2, REAL8 S1z, REAL8 S2z, REAL8 inc,
                                      REAL8 *time_of_merger);

/***********************************************************************************/
/************************** Static function definitions
 * ****************************/
/***********************************************************************************/
static inline void *my_realloc(void *p, size_t sz) {
  void *retval = LALRealloc(p, sz);
  return retval ? retval : (LALFree(p), NULL);
}

/* small powers */
static REAL8 pow2(const REAL8 x) { return x * x; }

static REAL8 pow3(const REAL8 x) { return x * pow2(x); }

static REAL8 pow4(const REAL8 x) { return pow2(x) * pow2(x); }

static REAL8 pow5(const REAL8 x) { return x * pow4(x); }

static REAL8 pow6(const REAL8 x) { return pow3(x) * pow3(x); }

static REAL8 pow7(const REAL8 x) { return x * pow6(x); }

/* fractional powers */
static REAL8 pow1_3(const REAL8 x) { return cbrt(x); }

// static REAL8 pow3_2(const REAL8 x) { return sqrt(pow3(x)); }

// static REAL8 pow5_2(const REAL8 x) { return sqrt(pow5(x)); }

static REAL8 pow7_2(const REAL8 x) { return sqrt(pow7(x)); }

static void PopulateKeplerParams(struct kepler_vars *params, const REAL8 e,
                                 const REAL8 x, const REAL8 r, const REAL8 rDOT,
                                 const REAL8 PhiDOT) {
  if (x > 0) {
    params->x = x;
    params->x2 = x * x;
    params->x3 = x * params->x2;
    params->x4 = x * params->x3;
    params->x5 = x * params->x4;
    params->x6 = x * params->x5;
    params->x7 = x * params->x6;
    params->x8 = x * params->x7;
    params->x9 = x * params->x8;
  }

  if (e > 0) {
    params->e = e;
    params->e2 = e * e;
    params->e3 = e * params->e2;
    params->e4 = e * params->e3;
    params->e5 = e * params->e4;
    params->e6 = e * params->e5;
    params->e7 = e * params->e6;
    params->e8 = e * params->e7;
    params->e9 = e * params->e8;
  }

  if (r > 0) {
    params->r = r;
    params->r2 = r * r;
    params->r3 = r * params->r2;
    params->r4 = r * params->r3;
    params->r5 = r * params->r4;
    params->r6 = r * params->r5;
    params->r7 = r * params->r6;
    params->r8 = r * params->r7;
    params->r9 = r * params->r8;
  }

  if (rDOT != 0) {
    params->rDOT = rDOT;
    params->rDOT2 = rDOT * rDOT;
    params->rDOT3 = rDOT * params->rDOT2;
    params->rDOT4 = rDOT * params->rDOT3;
    params->rDOT5 = rDOT * params->rDOT4;
    params->rDOT6 = rDOT * params->rDOT5;
    params->rDOT7 = rDOT * params->rDOT6;
    params->rDOT8 = rDOT * params->rDOT7;
    params->rDOT9 = rDOT * params->rDOT8;
  }

  if (PhiDOT != 0) {
    params->PhiDOT = PhiDOT;
    params->PhiDOT2 = PhiDOT * PhiDOT;
    params->PhiDOT3 = PhiDOT * params->PhiDOT2;
    params->PhiDOT4 = PhiDOT * params->PhiDOT3;
    params->PhiDOT5 = PhiDOT * params->PhiDOT4;
    params->PhiDOT6 = PhiDOT * params->PhiDOT5;
    params->PhiDOT7 = PhiDOT * params->PhiDOT6;
    params->PhiDOT8 = PhiDOT * params->PhiDOT7;
    params->PhiDOT9 = PhiDOT * params->PhiDOT8;
  }
}

#include "ENIGMA_GOterms.c"
#include "ENIGMA_GPEMergerRingdown.c"
#include "ENIGMA_PNInspiral.c"
#include "ENIGMA_PNMatch.c"

static REAL8 *interp_to_uniform_and_point(
    const REAL8 *t_vec, const REAL8 *val_vec, const long vec_length,
    const REAL8 t0, const REAL8 dt, const long length, eval_e_func eval_e,
    const REAL8 matching_time, REAL8 *matching_val) {
  int errorcode = XLAL_SUCCESS;

  gsl_interp *dynamics_interp = NULL;
  gsl_interp_accel *dynamics_accel = NULL;

  REAL8 *uniform = LALMalloc(sizeof(*uniform) * length);
  if (uniform == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }

  dynamics_interp = gsl_interp_alloc(gsl_interp_steffen, vec_length);
  dynamics_accel = gsl_interp_accel_alloc();
  if (dynamics_interp == NULL || dynamics_accel == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  gsl_interp_init(dynamics_interp, t_vec, val_vec, vec_length);

  for (long i = 0; i < length; i++) {
    // validate that the target time is within the acceptable range
    const REAL8 t = t0 + dt * (REAL8)i;
    int status =
        eval_e(dynamics_interp, t_vec, val_vec, t, dynamics_accel, uniform + i);
    if (status != GSL_SUCCESS) {
      errorcode = XLAL_EFAILED;
      XLAL_ERROR_FAIL(errorcode);
    }
  }

  int status = eval_e(dynamics_interp, t_vec, val_vec, matching_time,
                      dynamics_accel, matching_val);
  if (status != GSL_SUCCESS) {
    errorcode = XLAL_EFAILED;
    XLAL_ERROR_FAIL(errorcode);
  }

XLAL_FAIL:
  if (dynamics_interp != NULL)
    gsl_interp_free(dynamics_interp);
  if (dynamics_accel != NULL)
    gsl_interp_accel_free(dynamics_accel);

  if (errorcode != XLAL_SUCCESS) {
    LALFree(uniform);
    uniform = NULL;
  }
  return uniform;
}

static void compute_mode_from_dynamics(
    const UINT4 l, const INT4 m, UNUSED REAL8 *t_vec, REAL8 *x_vec,
    REAL8 *phi_vec, REAL8 *phi_dot_vec, REAL8 *r_vec, REAL8 *r_dot_vec,
    const REAL8 mass1, const REAL8 mass2, const REAL8 S1z, const REAL8 S2z,
    const REAL8 R, const long length, const UINT4 vpnorder, COMPLEX16 *h_lm) {

  assert(h_lm != NULL);

  const REAL8 total_mass = mass1 + mass2;
  const REAL8 eta =
      (mass1 * mass2) / (total_mass * total_mass); // Symmetric Mass Ratio

  /* store common powers of elements */
  struct kepler_vars orbital_vars;
  orbital_vars.eta = eta;
  orbital_vars.eta2 = eta * orbital_vars.eta;
  orbital_vars.eta3 = eta * orbital_vars.eta2;
  orbital_vars.eta4 = eta * orbital_vars.eta3;
  orbital_vars.eta5 = eta * orbital_vars.eta4;
  orbital_vars.eta6 = eta * orbital_vars.eta5;

  for (long i = 0; i < length; ++i) {
    PopulateKeplerParams(&orbital_vars, 0., x_vec[i], r_vec[i] * total_mass,
                         r_dot_vec[i], phi_dot_vec[i] / total_mass);
    h_lm[i] = hlmGOresult(l, m, total_mass, eta, r_vec[i] * total_mass,
                          r_dot_vec[i], phi_vec[i], phi_dot_vec[i] / total_mass,
                          R, vpnorder, S1z, S2z, x_vec[i], orbital_vars) *
              LAL_MRSUN_SI;
  }
}

static void compute_strain_from_dynamics(
    REAL8 *t_vec, REAL8 *x_vec, REAL8 *phi_vec, REAL8 *phi_dot_vec,
    REAL8 *r_vec, REAL8 *r_dot_vec, const REAL8 mass1, const REAL8 mass2,
    const REAL8 S1z, const REAL8 S2z,
    /* const REAL8 x0, */ const REAL8 euler_iota, const REAL8 euler_beta,
    const REAL8 R, const long length, const UINT4 vpnorder, REAL8 *h_plus,
    REAL8 *h_cross) {
  if (h_plus == NULL || h_cross == NULL)
    XLAL_ERROR_FAIL(XLAL_EINVAL);
  memset(h_plus, 0, length * sizeof(REAL8));
  memset(h_cross, 0, length * sizeof(REAL8));

  // Initialize and set to zero
  COMPLEX16 *hlm_timeseries =
      (COMPLEX16 *)XLALCalloc(length, sizeof(COMPLEX16));
  if (hlm_timeseries == NULL) {
    XLAL_ERROR_FAIL(XLAL_ENOMEM);
  }

  COMPLEX16 hlm_times_ylm = 0;
  COMPLEX16 ylm = 0;

  for (INT4 ell = L_MIN; ell <= L_MAX; ell++) {
    for (INT4 em = -ell; em <= (INT4)ell; em++) {
      // If set to use only those modes with l == |m|
      if (ONLY_LeqM_MODES && ell != abs(em))
        continue;

      compute_mode_from_dynamics(ell, em, t_vec, x_vec, phi_vec, phi_dot_vec,
                                 r_vec, r_dot_vec, mass1, mass2, S1z, S2z, R,
                                 length, vpnorder, hlm_timeseries);
      ylm = XLALSpinWeightedSphericalHarmonic(euler_iota, euler_beta, -2, ell,
                                              em);

      for (INT4 i = 0; i < length; i++) {
        hlm_times_ylm = hlm_timeseries[i] * ylm;
        h_plus[i] += creal(hlm_times_ylm);
        h_cross[i] -= cimag(hlm_times_ylm);
      }
    }
  }

  // Free temporary memory allocated to store mode timeseries
XLAL_FAIL:
  LALFree(hlm_timeseries);
}

// Note: this modifies t_vec, r_vec and phi_dot_vec by scaling by total_mass
// static void compute_strain_from_dynamics_old(
//     REAL8 *t_vec, REAL8 *x_vec, REAL8 *phi_vec, REAL8 *phi_dot_vec,
//     REAL8 *r_vec, REAL8 *r_dot_vec, const REAL8 mass1, const REAL8 mass2,
//     REAL8 S1z, REAL8 S2z,
//     /* const REAL8 x0, */ const REAL8 euler_iota, const REAL8 euler_beta,
//     const REAL8 R, const long length, REAL8 *h_plus, REAL8 *h_cross) {
//   assert(h_plus != NULL && h_cross != NULL);
//   const REAL8 total_mass = mass1 + mass2;
//   //const REAL8 reduced_mass = mass1 * mass2 / total_mass;
//   const REAL8 eta = (mass1*mass2)/(total_mass*total_mass); //Symmetric Mass
//   Ratio

//   /* overall factor in waveform polarizations */
//   /* REAL8 h_factor = reduced_mass * LAL_MRSUN_SI / R; */

//    /*The desired PN order correction.
//   Note that we are defining this in terms of powers of v = x^2*/
//   UINT4 pn_order_amp = 8; //4PN correction

//   for (long i = 0; i < length; ++i) {
//     t_vec[i] *= total_mass;
//     r_vec[i] *= total_mass;
//     phi_dot_vec[i] /= total_mass;

//     /* some math used below:
//      * sin(a-b) = sin(a) * cos(b) - cos(a) * sin(b)
//      * cos(a-b) = cos(a) * cos(b) + sin(a) * sin(b)
//      */

//     /* The leading order (quadrupolar) post-Newtonian GW polarizations *
//      * Reference: Damour, Gopakumar, and Iyer (PRD 70 064028).         */

//     //Defining the hplus and hcross general orbit term corrections
//     REAL8 hplusGOtotal=0;
//     REAL8 hcrossGOtotal=0;

//     /* printf("Show me the value of total mass:%e\n", total_mass);
//     fflush(NULL); */

//     for(long j=0;j<=pn_order_amp;j++){
//       hplusGOtotal = hplusGOtotal +   LAL_MRSUN_SI *
//       hplusGO(total_mass,eta,r_vec[i],r_dot_vec[i],phi_vec[i],phi_dot_vec[i],euler_iota,euler_beta,R,j,S1z,S2z,x_vec[i]);
//       hcrossGOtotal = hcrossGOtotal + LAL_MRSUN_SI *
//       hcrossGO(total_mass,eta,r_vec[i],r_dot_vec[i],phi_vec[i],phi_dot_vec[i],euler_iota,euler_beta,R,j,S1z,S2z,x_vec[i]);
//     }

//     h_plus[i] =
//         (REAL8)(/* h_factor *
//                 (-((cos(euler_iota) * cos(euler_iota) + 1.0) *
//                        ((total_mass / r_vec[i] +
//                          r_vec[i] * r_vec[i] * phi_dot_vec[i] *
//                          phi_dot_vec[i] - r_dot_vec[i] * r_dot_vec[i]) *
//                             (cos(2.0 * phi_vec[i]) * cos(2.0 * euler_beta)
//                             +
//                              sin(2.0 * phi_vec[i]) * sin(2.0 * euler_beta))
//                              +
//                         2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] *
//                             sin(2.0 * phi_vec[i])) +
//                    (total_mass / r_vec[i] -
//                     r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] -
//                     r_dot_vec[i] * r_dot_vec[i]) *
//                        sin(euler_iota) * sin(euler_iota))

//                  + hPlus(x_vec[i], x0, mass1, mass2, euler_iota,
//                  phi_vec[i], pn_order_amp ))
//                  + */ hplusGOtotal);

//     h_cross[i] =
//         (REAL8)(/* h_factor *
//                 ((-(2.0 * cos(euler_iota)) *
//                   ((total_mass / r_vec[i] +
//                     r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] -
//                     r_dot_vec[i] * r_dot_vec[i]) *
//                        (sin(2.0 * phi_vec[i]) * cos(2.0 * euler_beta) -
//                         cos(2.0 * phi_vec[i]) * sin(2.0 * euler_beta)) -
//                    2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] *
//                        (cos(2.0 * phi_vec[i]) * cos(2.0 * euler_beta) +
//                         sin(2.0 * phi_vec[i]) * sin(2.0 * euler_beta))))

//                  + hCross(x_vec[i], x0, mass1, mass2, euler_iota,
//                  phi_vec[i], pn_order_amp ))
//                  + */ hcrossGOtotal);
//   }

//   /* printf("GOING TO OPEN file and write %ld lines to it...\n", length);
//   fflush(NULL);
//   {
//     FILE *fout;
//     fout =
//     fopen("/Users/kaushikpaul/Dropbox/Enigma_spins_2023/fork_data_19Apr23.txt",
//     "w");

//     for(int i=0; i < length; ++i){
//       printf("Show me the value of hplus and hcross: (%e, %e)\n",h_plus[i],
//       h_cross[i]); fflush(NULL); fprintf(fout, "%e,%e,%e,%e,%e,%e,%e,%e\n",
//       t_vec[i], x_vec[i], phi_vec[i], phi_dot_vec[i], r_vec[i],
//       r_dot_vec[i], h_plus[i], h_cross[i]); fflush(fout);
//     }
//     fclose(fout);
//   } */
// }

static int x_model_eccbbh_inspiral_waveform(
    REAL8TimeSeries *h_plus, REAL8TimeSeries *h_cross,
    REAL8 mass1,          /* mass1 in solar mass    */
    REAL8 mass2,          /* mass2 in solar mass    */
    REAL8 S1z,            /* z-component of the spin of companion 1 */
    REAL8 S2z,            /* z-component of the spin of companion 2 */
    REAL8 e_init,         /* initial eccentricity   */
    REAL8 f_gw_init,      /* initial GW frequency   */
    REAL8 distance,       /* distance of source (m) */
    REAL8 mean_anom_init, /* initial mean-anomaly   */
    REAL8 ode_eps,        /* tolerance (relative)   */
    REAL8 euler_iota,     /* source polar angle     */
    REAL8 euler_beta,     /* source azimuthal angle */
    REAL8 sampling_rate   /* sample rate in Hz      */
) {

  int errorcode = XLAL_SUCCESS; /* the current error state value */

  /* Declare memory variables for storing dynamics */
  REAL8TimeSeries *time_evol = NULL;         /* Time steps */
  REAL8TimeSeries *x_evol = NULL;            /* Orbital x parameter */
  REAL8TimeSeries *eccentricity_evol = NULL; /* Orbital eccentricity */
  REAL8TimeSeries *mean_ano_evol = NULL;     /* Orbital mean anomaly */
  REAL8TimeSeries *phi_evol = NULL;          /* Orbital phase */
  REAL8TimeSeries *phi_dot_evol = NULL; /* Time derivative of orbital phase */
  REAL8TimeSeries *r_evol = NULL;       /* Orbital radius */
  REAL8TimeSeries *r_dot_evol = NULL;   /* Time derivative of orbital radius */
  REAL8Vector *Hp = NULL;               /* temporary store for strain */
  REAL8Vector *Hc = NULL;               /* temporary store for strain */

  if (!(h_plus && h_cross)) { // need to be passed existing object
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  /* computed variables */
  REAL8 total_mass = mass1 + mass2;
  int Length = -1;

  /* sampling interval [sec] */
  REAL8 dt = 1.0 / sampling_rate;

  /* scale the step times by t_sun */
  dt /= total_mass * LAL_MTSUN_SI;

  /* Declare memory variables for storing info on inspiral-merger attachment
   */
  REAL8 imr_matching_time;         /* attachment time */
  REAL8 imr_matching_x;            /* attachment x */
  REAL8 imr_matching_eccentricity; /* attachment e */
  REAL8 imr_matching_mean_ano;     /* attachment mean anomaly */
  REAL8 imr_matching_phi;          /* attachment phase */
  REAL8 imr_matching_phi_dot;      /* attachment phase time derivative */
  REAL8 imr_matching_r;            /* attachment radius */
  REAL8 imr_matching_r_dot;        /* attachment radial time derivative */

  // get inspiral dynamics
  errorcode = XLALSimInspiralENIGMADynamics(
      &time_evol, &x_evol, &eccentricity_evol, &mean_ano_evol, &phi_evol,
      &phi_dot_evol, &r_evol, &r_dot_evol, &imr_matching_time, &imr_matching_x,
      &imr_matching_eccentricity, &imr_matching_mean_ano, &imr_matching_phi,
      &imr_matching_phi_dot, &imr_matching_r, &imr_matching_r_dot, mass1, mass2,
      S1z, S2z, e_init, f_gw_init, mean_anom_init, ode_eps, sampling_rate,
      false);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  Length = time_evol->data->length;

  // get inspiral waveform
  errorcode = XLALSimInspiralENIGMAStrainFromDynamics(
      &Hp, &Hc, time_evol->data, x_evol->data, phi_evol->data,
      phi_dot_evol->data, r_evol->data, r_dot_evol->data, mass1, mass2, S1z,
      S2z, euler_iota, euler_beta, distance);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  /* store the sample interval in the output vectors */
  h_plus->data->length = h_cross->data->length = Length;
  h_plus->data->data = Hp->data;
  h_cross->data->data = Hc->data;
  Hc->data = Hp->data = NULL; // ownership has passed to h_plus and h_cross

XLAL_FAIL:
  XLALDestroyREAL8Sequence(Hp);
  XLALDestroyREAL8Sequence(Hc);
  XLALDestroyREAL8TimeSeries(time_evol);
  XLALDestroyREAL8TimeSeries(x_evol);
  XLALDestroyREAL8TimeSeries(eccentricity_evol);
  XLALDestroyREAL8TimeSeries(mean_ano_evol);
  XLALDestroyREAL8TimeSeries(phi_evol);
  XLALDestroyREAL8TimeSeries(phi_dot_evol);
  XLALDestroyREAL8TimeSeries(r_evol);
  XLALDestroyREAL8TimeSeries(r_dot_evol);
  time_evol = x_evol = eccentricity_evol = mean_ano_evol = phi_evol =
      phi_dot_evol = r_evol = r_dot_evol = NULL;

  // We do not free h_plus or h_minus as this is the callers duty

  return errorcode;
}

static int x_model_eccbbh_imr_waveform(
    REAL8TimeSeries *h_plus, REAL8TimeSeries *h_cross,
    REAL8 mass1,          /* mass1 in solar mass    */
    REAL8 mass2,          /* mass2 in solar mass    */
    REAL8 S1z,            /* z-component of the spin of companion 1 */
    REAL8 S2z,            /* z-component of the spin of companion 2 */
    REAL8 e_init,         /* initial eccentricity   */
    REAL8 f_gw_init,      /* initial GW frequency   */
    REAL8 distance,       /* distance of source (m) */
    REAL8 mean_anom_init, /* initial mean-anomaly   */
    REAL8 ode_eps,        /* tolerance (relative)   */
    REAL8 euler_iota,     /* source polar angle     */
    REAL8 euler_beta,     /* source azimuthal angle */
    REAL8 sampling_rate,  /* sample rate in Hz      */
    TrainingSet *Dset) {
  int errorcode = XLAL_SUCCESS; /* the current error state value */

  /* Declare memory variables for storing dynamics */
  REAL8TimeSeries *time_evol = NULL;         /* Time steps */
  REAL8TimeSeries *x_evol = NULL;            /* Orbital x parameter */
  REAL8TimeSeries *eccentricity_evol = NULL; /* Orbital eccentricity */
  REAL8TimeSeries *mean_ano_evol = NULL;     /* Orbital mean anomaly */
  REAL8TimeSeries *phi_evol = NULL;          /* Orbital phase */
  REAL8TimeSeries *phi_dot_evol = NULL; /* Time derivative of orbital phase */
  REAL8TimeSeries *r_evol = NULL;       /* Orbital radius */
  REAL8TimeSeries *r_dot_evol = NULL;   /* Time derivative of orbital radius */
  REAL8Vector *Hp = NULL;               /* temporary store for strain */
  REAL8Vector *Hc = NULL;               /* temporary store for strain */

  if (!(h_plus && h_cross)) { // need to be passed existing object
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  /* computed variables */
  REAL8 total_mass = mass1 + mass2;
  int Length = -1;

  /* sampling interval [sec] */
  REAL8 dt = 1.0 / sampling_rate;

  /* scale the step times by t_sun */
  dt /= total_mass * LAL_MTSUN_SI;

  /* Declare memory variables for storing info on inspiral-merger attachment
   */
  REAL8 imr_matching_time;         /* attachment time */
  REAL8 imr_matching_x;            /* attachment x */
  REAL8 imr_matching_eccentricity; /* attachment e */
  REAL8 imr_matching_mean_ano;     /* attachment mean anomaly */
  REAL8 imr_matching_phi;          /* attachment phase */
  REAL8 imr_matching_phi_dot;      /* attachment phase time derivative */
  REAL8 imr_matching_r;            /* attachment radius */
  REAL8 imr_matching_r_dot;        /* attachment radial time derivative */

  // get inspiral dynamics
  errorcode = XLALSimInspiralENIGMADynamics(
      &time_evol, &x_evol, &eccentricity_evol, &mean_ano_evol, &phi_evol,
      &phi_dot_evol, &r_evol, &r_dot_evol, &imr_matching_time, &imr_matching_x,
      &imr_matching_eccentricity, &imr_matching_mean_ano, &imr_matching_phi,
      &imr_matching_phi_dot, &imr_matching_r, &imr_matching_r_dot, mass1, mass2,
      S1z, S2z, e_init, f_gw_init, mean_anom_init, ode_eps, sampling_rate,
      true);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  Length = (int)ceil(imr_matching_time / dt) + 1;

  // get inspiral waveform
  errorcode = XLALSimInspiralENIGMAStrainFromDynamics(
      &Hp, &Hc, time_evol->data, x_evol->data, phi_evol->data,
      phi_dot_evol->data, r_evol->data, r_dot_evol->data, mass1, mass2, S1z,
      S2z, euler_iota, euler_beta, distance);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  // get strain at matching time
  REAL8 matching_Hp, matching_Hc;
  REAL8 t_val =
      imr_matching_time; // is modified by compute_strain_from_dynamics
  compute_strain_from_dynamics(
      &t_val, &imr_matching_x, &imr_matching_phi, &imr_matching_phi_dot,
      &imr_matching_r, &imr_matching_r_dot, mass1, mass2, S1z,
      S2z, /* x_evol->data->data[0], */
      euler_iota, euler_beta, distance, 1, (UINT4)mode_pn_order, &matching_Hp,
      &matching_Hc);

  // get merger and ringdown and attach at the end of the existing time series
  REAL8 time_of_merger = 0.0;
  errorcode = Attach_GPE_Merger_Ringdown(
      dt, &(Hp->data), &(Hc->data), imr_matching_time, matching_Hp, matching_Hc,
      &Length, Dset, mass1, mass2, S1z, S2z, euler_iota, &time_of_merger);
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(errorcode);

  /* store the sample interval in the output vectors */
  h_plus->data->length = h_cross->data->length = Length;
  h_plus->data->data = Hp->data;
  h_cross->data->data = Hc->data;
  Hc->data = Hp->data = NULL; // ownership has passed to h_plus and h_cross

  /* Set epoch to the length of the inspiral waveform, as is usual convention
   */
  XLALGPSSetREAL8(&(h_plus->epoch),
                  -time_of_merger * total_mass * LAL_MTSUN_SI);
  XLALGPSSetREAL8(&(h_cross->epoch),
                  -time_of_merger * total_mass * LAL_MTSUN_SI);

XLAL_FAIL:
  XLALDestroyREAL8Sequence(Hp);
  XLALDestroyREAL8Sequence(Hc);
  XLALDestroyREAL8TimeSeries(time_evol);
  XLALDestroyREAL8TimeSeries(x_evol);
  XLALDestroyREAL8TimeSeries(eccentricity_evol);
  XLALDestroyREAL8TimeSeries(mean_ano_evol);
  XLALDestroyREAL8TimeSeries(phi_evol);
  XLALDestroyREAL8TimeSeries(phi_dot_evol);
  XLALDestroyREAL8TimeSeries(r_evol);
  XLALDestroyREAL8TimeSeries(r_dot_evol);
  time_evol = x_evol = eccentricity_evol = mean_ano_evol = phi_evol =
      phi_dot_evol = r_evol = r_dot_evol = NULL;

  // We do not free h_plus or h_minus as this is the callers duty

  return errorcode;
}

/***********************************************************************************/
/*************************** Global function Definitions
 * ***************************/
/***********************************************************************************/

// Note: this modifies t_vec, r_vec and phi_dot_vec by scaling by total_mass
int XLALSimInspiralENIGMAModeFromDynamics(
    COMPLEX16Vector **h_lm, const UINT4 l, const INT4 m, REAL8Vector *t_vector,
    REAL8Vector *x_vector, REAL8Vector *phi_vector, REAL8Vector *phi_dot_vector,
    REAL8Vector *r_vector, REAL8Vector *r_dot_vector, const REAL8 mass1,
    const REAL8 mass2, const REAL8 S1z, const REAL8 S2z, const REAL8 R) {
  int errorcode = XLAL_SUCCESS;

  UINT4 length = t_vector->length;

  // check input parameters
  if (!h_lm) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }

  *h_lm = XLALCreateCOMPLEX16Sequence(length);

  if (!h_lm) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  compute_mode_from_dynamics(
      l, m, t_vector->data, x_vector->data, phi_vector->data,
      phi_dot_vector->data, r_vector->data, r_dot_vector->data, mass1, mass2,
      S1z, S2z, R, length, (UINT4)mode_pn_order, (*h_lm)->data);

XLAL_FAIL:
  // We do not free h_lm
  return errorcode;
}

// Note: this modifies t_vec, r_vec and phi_dot_vec by scaling by total_mass
int XLALSimInspiralENIGMAStrainFromDynamics(
    REAL8Vector **h_plus, REAL8Vector **h_cross, REAL8Vector *t_vector,
    REAL8Vector *x_vector, REAL8Vector *phi_vector, REAL8Vector *phi_dot_vector,
    REAL8Vector *r_vector, REAL8Vector *r_dot_vector, const REAL8 mass1,
    const REAL8 mass2, const REAL8 S1z, const REAL8 S2z, const REAL8 euler_iota,
    const REAL8 euler_beta, const REAL8 R) {

  int errorcode = XLAL_SUCCESS;

  UINT4 length = t_vector->length;

  // check input parameters
  if (!h_plus || !h_cross) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }

  *h_plus = XLALCreateREAL8Vector(length);
  *h_cross = XLALCreateREAL8Vector(length);

  if (!h_plus || !h_cross) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  compute_strain_from_dynamics(
      t_vector->data, x_vector->data, phi_vector->data, phi_dot_vector->data,
      r_vector->data, r_dot_vector->data, mass1, mass2, S1z,
      S2z, /* x_vector->data[0], */
      euler_iota, euler_beta, R, length, (UINT4)mode_pn_order, (*h_plus)->data,
      (*h_cross)->data);

XLAL_FAIL:
  // We do not free h_plus or h_cross as this is the callers duty
  return errorcode;
}

int XLALSimInspiralENIGMADynamics(
    REAL8TimeSeries **time_evol,         /* Time steps */
    REAL8TimeSeries **x_evol,            /* Orbital x parameter */
    REAL8TimeSeries **eccentricity_evol, /* Orbital eccentricity */
    REAL8TimeSeries **mean_ano_evol,     /* Orbital mean anomaly */
    REAL8TimeSeries **phi_evol,          /* Orbital phase */
    REAL8TimeSeries **phi_dot_evol,      /* Time derivative of orbital phase */
    REAL8TimeSeries **r_evol,            /* Orbital radius */
    REAL8TimeSeries **r_dot_evol,        /* Time derivative of orbital radius */
    REAL8 *imr_matching_time,            /* I + MR attachment time */
    REAL8 *imr_matching_x,               /* I + MR attachment x */
    REAL8 *imr_matching_eccentricity,    /* I + MR attachment e */
    REAL8 *imr_matching_mean_ano,        /* I + MR attachment mean anomaly */
    REAL8 *imr_matching_phi,             /* I + MR attachment phase */
    REAL8 *imr_matching_phi_dot, /* I + MR attachment phase time derivative */
    REAL8 *imr_matching_r,       /* I + MR attachment radius */
    REAL8 *imr_matching_r_dot,   /* I + MR attachment radial time derivative */
    REAL8 mass1,                 /* mass1 in solar mass    */
    REAL8 mass2,                 /* mass2 in solar mass   dynamics */
    REAL8 S1z,                   /* z-component of the spin of companion 1 */
    REAL8 S2z,                   /* z-component of the spin of companion 2*/
    REAL8 e_init,                /* initial eccentricity   */
    REAL8 f_gw_init,             /* initial GW frequency   */
    REAL8 mean_anom_init,        /* initial mean-anomaly   */
    REAL8 ode_eps,               /* tolerance (relative)   */
    REAL8 sampling_rate,         /* sample rate in Hz      */
    bool will_attach_mr /* Enables integration to attachment frequency */
) {

  int errorcode = XLAL_SUCCESS; /* the current error state value */

  // check input parameters
  if (!time_evol || !x_evol || !eccentricity_evol || !mean_ano_evol ||
      !phi_evol || !phi_dot_evol || !r_evol || !r_dot_evol) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }
  if (*time_evol || *x_evol || *eccentricity_evol || *mean_ano_evol ||
      *phi_evol || *phi_dot_evol || *r_evol || *r_dot_evol) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }

  /* parameters for the ODE system */
  struct kepler_vars orbital_vars;
  struct ode_parameters ecc_params;
  ecc_params.orbital_vars = &orbital_vars;

  int rad_pn_order = -1;

  long int i, j, final_i_reached = 0;
  int badnumber = 0;
  int status = 0;

  /* dynamical quantities */
  REAL8 y[4];     /* array containing x, e, l, phi                  */
  REAL8 y_dot[4]; /* array containing x_dot, e_dot, l_dot, phi_dot  */
  REAL8 y_temp[4];
  REAL8 yerr[4];
  REAL8 y_dot_in[4];
  REAL8 y_dot_out[4];
  REAL8 y_dot_temp[4];

  /* vectors to store the data */
  REAL8 *uniform_t_vec = NULL;   /* uniformly sampled time */
  REAL8 *uniform_x_vec = NULL;   /* uniformly sampled PN expansion variable */
  REAL8 *uniform_phi_vec = NULL; /* uniformly sampled orbital phase */
  REAL8 *uniform_phi_dot_vec = NULL; /* uniformly sampled orbital frequency */
  REAL8 *uniform_r_vec = NULL;       /* uniformly sampled orbital separation */
  REAL8 *uniform_r_dot_vec =
      NULL; /* uniformly sampled orbital separation change rate */
  REAL8 *uniform_e_vec = NULL; /* uniformly sampled eccentricity */
  REAL8 *uniform_l_vec = NULL; /* uniformly sampled mean anomaly */
  REAL8 *t_vec = NULL;         /* time                  */
  REAL8 *x_vec = NULL;         /* PN expansion variable */
  REAL8 *e_vec = NULL;         /* eccentricity          */
  REAL8 *l_vec = NULL;         /* mean anomaly          */
  REAL8 *u_vec = NULL;         /* eccentric anomaly     */
  REAL8 *phi_vec = NULL;       /* orbital phase         */
  REAL8 *phi_dot_vec = NULL;   /* orbital frequency     */
  REAL8 *r_vec = NULL;         /* orbital separation    */
  REAL8 *p_vec = NULL;         /* semi-latus rectum     */

  /* computed variables */
  REAL8 omega_init;
  REAL8 x_init;
  REAL8 x_final;
  REAL8 f_gw_isco; /* GW frequency and isco     */
  REAL8 total_mass;
  REAL8 mass_ratio;
  REAL8 reduced_mass, sym_mass_ratio;

  /* time stepping variables */
  REAL8 t, dt, t_next;
  /* step size returned by the ODE stepping routine */
  REAL8 step_size;

  /* interpolation helpers */
  gsl_interp *phi_dot_interp = NULL;
  gsl_interp_accel *phi_dot_accel = NULL;

  /* ODE solverr objects */
  gsl_odeiv_step *solver_step = NULL;
  gsl_odeiv_control *solver_control = NULL;
  gsl_odeiv_evolve *solver_evolve = NULL;

  /* Allocate memory for dynamics */
  LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
  *time_evol = XLALCreateREAL8TimeSeries("time_evol", &epoch, f_gw_init,
                                         1 / sampling_rate, &lalStrainUnit, 0);
  *x_evol = XLALCreateREAL8TimeSeries("x_evol", &epoch, f_gw_init,
                                      1 / sampling_rate, &lalStrainUnit, 0);
  *eccentricity_evol =
      XLALCreateREAL8TimeSeries("eccentricity_evol", &epoch, f_gw_init,
                                1 / sampling_rate, &lalStrainUnit, 0);
  *mean_ano_evol = XLALCreateREAL8TimeSeries(
      "mean_ano_evol", &epoch, f_gw_init, 1 / sampling_rate, &lalStrainUnit, 0);
  *phi_evol = XLALCreateREAL8TimeSeries("phi_evol", &epoch, f_gw_init,
                                        1 / sampling_rate, &lalStrainUnit, 0);
  *phi_dot_evol = XLALCreateREAL8TimeSeries(
      "phi_dot_evol", &epoch, f_gw_init, 1 / sampling_rate, &lalStrainUnit, 0);
  *r_evol = XLALCreateREAL8TimeSeries("r_evol", &epoch, f_gw_init,
                                      1 / sampling_rate, &lalStrainUnit, 0);
  *r_dot_evol = XLALCreateREAL8TimeSeries("r_dot_evol", &epoch, f_gw_init,
                                          1 / sampling_rate, &lalStrainUnit, 0);
  if (!time_evol || !x_evol || !eccentricity_evol || !mean_ano_evol ||
      !phi_evol || !phi_dot_evol || !r_evol || !r_dot_evol) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  /* check the command line arguments */
  if (e_init < 0.0 || e_init >= 1.0) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode,
                    "ERROR: Invalid eccentricity, must be in range [0,1)");
  }

  total_mass = mass1 + mass2;
  reduced_mass = mass1 * mass2 / total_mass;
  sym_mass_ratio = reduced_mass / total_mass;
  mass_ratio = mass1 / mass2;

  /* initial orbital angular frequency in units of solar mass */
  omega_init = LAL_PI * f_gw_init * LAL_MTSUN_SI;
  x_init = pow(total_mass * omega_init, 2. / 3.);

  /* Set con and pn orders to cover (q,M) space*/

  REAL8 omega_attach;

  if (will_attach_mr) {
    errorcode = PN_Omega(mass_ratio, total_mass, &rad_pn_order, &omega_attach);
    if (errorcode != XLAL_SUCCESS)
      XLAL_ERROR_FAIL(XLAL_EFUNC);
    // omega_attach*= 0.5;
  } else {
    rad_pn_order = (UINT4) Radiation_PN_Order;
    omega_attach = 1.0; // deliberately use unphysical value.
  }
  /* store the mass and pn params in the param structure */
  ecc_params.eta = sym_mass_ratio;
  ecc_params.radiation_pn_order = rad_pn_order;
  ecc_params.m1 = mass1;
  ecc_params.m2 = mass2;
  ecc_params.S1z = S1z;
  ecc_params.S2z = S2z;

  /* GSL ode solver error tolerances */
  REAL8 absolute_step_error = 0.0;
  REAL8 relative_step_error = ode_eps;

  /* GSL Runge-Kutta Fehlberg 4-5 ode stepping method */
  const gsl_odeiv_step_type *solver_type = gsl_odeiv_step_rkf45;

  solver_step = gsl_odeiv_step_alloc(solver_type, 4);
  solver_control =
      gsl_odeiv_control_y_new(absolute_step_error, relative_step_error);
  solver_evolve = gsl_odeiv_evolve_alloc(4);
  gsl_odeiv_system solver_system = {eccentric_x_model_odes, NULL, 4,
                                    &ecc_params};
  if (solver_step == NULL || solver_control == NULL || solver_evolve == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }

  /* Number of allowed retries at each step of the adaptive-stepping rk45
   * integrator */
  int num_retries = 500;
  int tempretries = num_retries;

  /* sampling interval [sec] */
  dt = 1.0 / sampling_rate;

  /* Transition to merger */
  REAL8 TRANS = 4.;

  /* grab memory for vector variables */
  int statevec_allocated;
#define allocate_statevec(n)                                                   \
  do {                                                                         \
    t_vec = (REAL8 *)my_realloc(t_vec, (n) * sizeof(REAL8));                   \
    x_vec = (REAL8 *)my_realloc(x_vec, (n) * sizeof(REAL8));                   \
    e_vec = (REAL8 *)my_realloc(e_vec, (n) * sizeof(REAL8));                   \
    phi_vec = (REAL8 *)my_realloc(phi_vec, (n) * sizeof(REAL8));               \
    l_vec = (REAL8 *)my_realloc(l_vec, (n) * sizeof(REAL8));                   \
    p_vec = (REAL8 *)my_realloc(p_vec, (n) * sizeof(REAL8));                   \
    u_vec = (REAL8 *)my_realloc(u_vec, (n) * sizeof(REAL8));                   \
    phi_dot_vec = (REAL8 *)my_realloc(phi_dot_vec, (n) * sizeof(REAL8));       \
    r_vec = (REAL8 *)my_realloc(r_vec, (n) * sizeof(REAL8));                   \
    if (t_vec == NULL || x_vec == NULL || e_vec == NULL || phi_vec == NULL ||  \
        l_vec == NULL || p_vec == NULL || u_vec == NULL ||                     \
        phi_dot_vec == NULL || r_vec == NULL) {                                \
      errorcode = XLAL_ENOMEM;                                                 \
      XLAL_ERROR_FAIL(errorcode);                                              \
    }                                                                          \
    statevec_allocated = (n);                                                  \
  } while (0)
  allocate_statevec(1);

  /* scale the start, end and step times by t_sun */
  dt /= total_mass * LAL_MTSUN_SI;

  /* termination condition:                         *
   * Schwarzschild innermost stable circular orbit */

  /* gw frequency (geometrized) at isco */
  f_gw_isco = 1.0 / (TRANS * sqrt(TRANS) * LAL_PI * total_mass);

  /* final value of x is at x(f_isco) = 1/6 */
  x_final = pow(LAL_PI * total_mass * f_gw_isco, 2. / 3.);

  /* initial conditions on dynamical variables */
  t_vec[0] = t = 0.;
  x_vec[0] = y[0] = x_init;
  e_vec[0] = y[1] = e_init;
  l_vec[0] = y[2] = mean_anom_init; /* initial mean anomaly */
  phi_vec[0] = y[3] = 0.0;          /* inital phase */
  p_vec[0] = (1.0 - e_vec[0] * e_vec[0]) /
             pow(LAL_PI * total_mass * LAL_MTSUN_SI * f_gw_init, 2. / 3.);

  /* initial eccentric anomaly */
  u_vec[0] = pn_kepler_equation(sym_mass_ratio, x_vec[0], e_vec[0], l_vec[0]);

  /* compute the intital value of r using Eq. (5) */
  r_vec[0] = separation(u_vec[0], sym_mass_ratio, x_vec[0], e_vec[0], mass1,
                        mass2, S1z, S2z);

  /* compute derivatives at the current value of t */
  eccentric_x_model_odes(t, y, y_dot_in, (void *)&ecc_params);

  /* initial orbital frequency */
  phi_dot_vec[0] = y_dot_in[3];

  /* store common powers of elements */
  orbital_vars.eta = sym_mass_ratio;
  orbital_vars.eta2 = sym_mass_ratio * orbital_vars.eta;
  orbital_vars.eta3 = sym_mass_ratio * orbital_vars.eta2;
  orbital_vars.eta4 = sym_mass_ratio * orbital_vars.eta3;
  orbital_vars.eta5 = sym_mass_ratio * orbital_vars.eta4;
  orbital_vars.eta6 = sym_mass_ratio * orbital_vars.eta5;
  PopulateKeplerParams(&orbital_vars, e_init, x_init, r_vec[0], 0,
                       phi_dot_vec[0]);

  /* evolve the dynamical variables forward until we reach  */
  /* termination condition: phi_dot > omega_attach or x > x_final */
  t = t_next = 0.;
  step_size = 2 * LAL_PI / phi_dot_vec[0] /
              100; /* start with ~100 steps for one orbit */

  // sanity check on omega_attach
  if (will_attach_mr && (phi_dot_vec[0] > omega_attach)) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode,
                    "initial phi_dot > omega_attach: %g > %g, lower f_min to "
                    "produce this waveform",
                    phi_dot_vec[0], omega_attach);
  }

  REAL8 time_omega_attach_reached = -1.;
  int i_omega_attach_reached = -1;

  for (i = 1; /* no end */; ++i) { /*{{{*/

    if (i >= statevec_allocated) {
      // Limit the maximum size a waveform can have to 2048 seconds at 16kHz.
      if (statevec_allocated >= 2048 * 16384) {
#if DEBUG
        FILE *fp = fopen("large_mem_waves.txt", "a");
        fprintf(fp, "%lf, %lf, %e, %e\n", mass1, mass2, e_init, mean_anom_init);
        fflush(fp);
        fclose(fp);
#endif
        break;
      }
      allocate_statevec(2 * statevec_allocated);
    }
    /* integrate to the next time step */
    do {
      memcpy(y_dot_temp, y_dot_in, 4 * sizeof(REAL8));
      memcpy(y_temp, y, 4 * sizeof(REAL8));

      do {
        status = gsl_odeiv_step_apply(solver_step, t, step_size, y, yerr,
                                      y_dot_in, y_dot_out, &solver_system);

        if (status != GSL_SUCCESS) {
          if (tempretries > 0) {
            tempretries--;
            step_size /= 10.;
            continue;
          } else {
            errorcode = XLAL_EFAILED;
            XLAL_ERROR_FAIL(errorcode);
          }
        }
      } while (status != GSL_SUCCESS);

      t_next = t + step_size;

      status = gsl_odeiv_control_hadjust(solver_control, solver_step, y, yerr,
                                         y_dot_out, &step_size);

      if (status == GSL_ODEIV_HADJ_DEC) {
        memcpy(y, y_temp, 4 * sizeof(REAL8));
        memcpy(y_dot_in, y_dot_temp, 4 * sizeof(REAL8));
      }
    } while (status == GSL_ODEIV_HADJ_DEC);

    t = t_next;
    memcpy(y_dot_in, y_dot_out, 4 * sizeof(REAL8));

    /* check for nan or inf in dynamical variables */
    for (j = 0; j < 4; ++j) {
      if (isnan(y[j]) || isinf(y[j])) {
        badnumber = 1;
      }
    }
    /* stop integrating if we get a nan or an inf */
    if (badnumber) {
      break;
    }

    /* compute derivatives at the current value of t */
    eccentric_x_model_odes(t, y, y_dot, (void *)&ecc_params);

    /* store the computed time, dynamical variables and their derivatives */
    t_vec[i] = t;
    x_vec[i] = y[0];
    e_vec[i] = y[1];
    l_vec[i] = y[2];
    phi_vec[i] = y[3];
    phi_dot_vec[i] = y_dot[3];

    /* semi-latus rectum: correct only for 0pN R.R. 0pN Con. *
     * p = (1-e^2)/(M n)^{2/3}, where M n = x^{3/2}          */
    p_vec[i] = (1.0 - e_vec[i] * e_vec[i]) / (x_vec[i]);

    /* compute the current value of the eccentric anomaly */
    u_vec[i] = pn_kepler_equation(sym_mass_ratio, x_vec[i], e_vec[i], l_vec[i]);

    /* compute the values of r using Eq. (5) */
    r_vec[i] = separation(u_vec[i], sym_mass_ratio, x_vec[i], e_vec[i], mass1,
                          mass2, S1z, S2z);

    /* store common powers of elements */
    PopulateKeplerParams(&orbital_vars, e_vec[i], x_vec[i], r_vec[i], 0,
                         phi_dot_vec[i]);

    /* check if we reached ISCO (we should never get this far) */
    if (x_vec[i] >= x_final) {
      break;
    }

    if (will_attach_mr && (phi_dot_vec[i] >= omega_attach)) {
      if (time_omega_attach_reached < 0.) {
        if (i - 1 < 0) { /* reached omega_attach (way) too early */
          errorcode = XLAL_EMAXITER;
          XLAL_ERROR_FAIL(errorcode);
        }
        time_omega_attach_reached =
            t_vec[i - 1]; /* time just before omega_attach is reached */
        i_omega_attach_reached = i - 1;
      } else if (t > time_omega_attach_reached + 10 * dt) {
        /* we need at least one point i*dt past time_omega_attach_reached and
         * add a couple more for spline interpolation to work nicely */
        break;
      }
    }

  } /*}}}*/ /* end integration for loop */

  if (will_attach_mr &&
      (time_omega_attach_reached <
       0.)) { /* never reached omega_attach before encountering ISCO */
    errorcode = XLAL_EMAXITER;
    XLAL_ERROR_FAIL(errorcode);
  }

  final_i_reached = i;
  if (badnumber == 0) {
    final_i_reached = i + 1;
  }

  // determine when attachment point was reached
  /*==========================================================================================================*/
  int Length = -1;
  REAL8 matching_time = -1;
  if (will_attach_mr) {
    phi_dot_interp = gsl_interp_alloc(gsl_interp_steffen, final_i_reached);
    phi_dot_accel = gsl_interp_accel_alloc();
    if (phi_dot_interp == NULL || phi_dot_accel == NULL) {
      errorcode = XLAL_ENOMEM;
      XLAL_ERROR_FAIL(errorcode);
    }
    gsl_interp_init(phi_dot_interp, t_vec, phi_dot_vec, final_i_reached);

    assert(i_omega_attach_reached + 1 < final_i_reached);
    REAL8 t_lo = t_vec[i_omega_attach_reached];
    REAL8 t_hi = t_vec[i_omega_attach_reached + 1];

    // sanity check on initial bracket for bisection
    {
      REAL8 phi_dot_val = -1;
      // Sets phi_dot_val by interpolating the data
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec, t_lo,
                                 phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val > omega_attach) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(
            errorcode,
            "Unexpected ordering of phi_dot[] array. Expected phi_dot(t_lo) <= "
            "omega_attach, but found phi_dot(%g) = %g = %g > omega_attach = %g",
            t_lo, phi_dot_val, phi_dot_vec[i_omega_attach_reached],
            omega_attach);
      }
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec, t_hi,
                                 phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(
            errorcode,
            "Unexpected ordering of phi_dot[] array. Expected phi_dot(t_hi) >= "
            "omega_attach, but found phi_dot(%g) = %g < omega_attach = %g",
            t_hi, phi_dot_val, omega_attach);
      }
      if (phi_dot_val < omega_attach) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(errorcode);
      }
    }

    // locate the matching time to within ~1e-3 of the sampling rate
    const double epsabs = 1e-3 * dt;
    const double epsrel = 0;
    const int max_iter = 100;
    for (int niter = 0;
         (status = gsl_root_test_interval(t_lo, t_hi, epsabs, epsrel)) !=
         GSL_SUCCESS;
         niter++) {
      if (status != GSL_CONTINUE) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (niter > max_iter) {
        errorcode = XLAL_EMAXITER;
        XLAL_ERROR_FAIL(errorcode);
      }

      matching_time = 0.5 * (t_lo + t_hi);
      REAL8 phi_dot_val = -1;
      // Sets phi_dot_val by interpolating the data
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec,
                                 matching_time, phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val < omega_attach) {
        t_lo = matching_time;
      } else if (phi_dot_val > omega_attach) {
        t_hi = matching_time;
      } else {
        break;
      }
    }

    // sanity check on final bracket after bisection
    {
      assert(dt > epsabs);
      if (t_hi - t_lo >
          1.05 * epsabs) // give some sleeve for roundoff and the like
      {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(errorcode,
                        "bisection did not reduce bracket width [%.15e,%.18e] "
                        "below threshold of %.18e for sampling rate %.18e",
                        t_lo, t_hi, epsabs, dt);
      }
      REAL8 phi_dot_val = -1;
      // Sets phi_dot_val by interpolating the data
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec, t_lo,
                                 phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val > omega_attach) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(errorcode);
      }
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec, t_hi,
                                 phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val < omega_attach) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(errorcode);
      }
    }
    Length = (int)ceil(t_lo / dt) + 1;
    if (Length < 2) {
      errorcode = XLAL_EBADLEN;
      XLAL_ERROR_FAIL(errorcode);
    }
    // if the bracket [t_lo,t_hi] straddles a grid point then we need to move
    // the bracket fully to one side of the grid point using one more bisection
    // step
    if (t_lo < (Length - 1) * dt && (Length - 1) * dt < t_hi) {
      REAL8 phi_dot_val = -1;
      matching_time = (Length - 1) * dt;
      status = gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec,
                                 matching_time, phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val < omega_attach) {
        t_lo = matching_time;
        Length += 1;
      } else if (phi_dot_val > omega_attach) {
        t_hi = matching_time;
      }
    }
    matching_time = 0.5 * (t_lo + t_hi);

    // sanity check on length
    {
      REAL8 phi_dot_val = -1;
      // Sets phi_dot_val by interpolating the data
      status =
          gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec,
                            (Length - 2) * dt, phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val > omega_attach || (Length - 2) * dt > matching_time) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(
            errorcode,
            "phi_dot_val: %.18e, omega_attach: %.18e, (Length - 2) * "
            "dt: %.18e, t_lo: %.18e, matching_time: %.18e, t_hi: "
            "%.18e, (Length - 1) * dt: %.18e",
            phi_dot_val, omega_attach, (Length - 2) * dt, t_lo, matching_time,
            t_hi, (Length - 1) * dt);
      }
      status =
          gsl_interp_eval_e(phi_dot_interp, t_vec, phi_dot_vec,
                            (Length - 1) * dt, phi_dot_accel, &phi_dot_val);
      if (status != GSL_SUCCESS) {
        errorcode = XLAL_EFAILED;
        XLAL_ERROR_FAIL(errorcode);
      }
      if (phi_dot_val < omega_attach || (Length - 1) * dt < matching_time) {
        errorcode = XLAL_ETOL;
        XLAL_ERROR_FAIL(
            errorcode,
            "phi_dot_val: %.18e, omega_attach: %.18e, (Length - 2) * "
            "dt: %.18e, t_lo: %.18e, matching_time: %.18e, t_hi: "
            "%.18e, (Length - 1) * dt: %.18e",
            phi_dot_val, omega_attach, (Length - 2) * dt, t_lo, matching_time,
            t_hi, (Length - 1) * dt);
      }
    }
  } else {
    Length = (int)ceil(t_vec[final_i_reached - 2] / dt);
    matching_time = t_vec[final_i_reached - 2];
  }
  // get uniform data for inspiral
  uniform_t_vec = (REAL8 *)LALMalloc(Length * sizeof(REAL8));
  if (uniform_t_vec == NULL) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(errorcode);
  }
  for (i = 0; i < Length; i++)
    uniform_t_vec[i] = dt * (REAL8)i;

  REAL8 matching_x, matching_phi, matching_phi_dot, matching_r, matching_r_dot,
      matching_e, matching_l;
  struct interp_args_t {
    REAL8 *vals, **interp_vals;
    REAL8 *matching_val;
    eval_e_func eval_e;
  } uniform_args[] = {
      {x_vec, &uniform_x_vec, &matching_x, gsl_interp_eval_e},
      {phi_vec, &uniform_phi_vec, &matching_phi, gsl_interp_eval_e},
      {phi_dot_vec, &uniform_phi_dot_vec, &matching_phi_dot, gsl_interp_eval_e},
      {r_vec, &uniform_r_vec, &matching_r, gsl_interp_eval_e},
      {r_vec, &uniform_r_dot_vec, &matching_r_dot, gsl_interp_eval_deriv_e},
      {e_vec, &uniform_e_vec, &matching_e, gsl_interp_eval_e},
      {l_vec, &uniform_l_vec, &matching_l, gsl_interp_eval_e},
  };
  for (i = 0; i < (long)(sizeof(uniform_args) / sizeof(uniform_args[0])); ++i) {
    REAL8 *uniform_vals = interp_to_uniform_and_point(
        t_vec, uniform_args[i].vals, final_i_reached, 0., dt, Length,
        uniform_args[i].eval_e, matching_time, uniform_args[i].matching_val);
    if (uniform_vals == NULL) {
      XLAL_ERROR_FAIL(XLAL_EFUNC);
    }
    *uniform_args[i].interp_vals = uniform_vals;
  }

  /* store the sample interval in the output vectors */
  *imr_matching_time = matching_time;
  *imr_matching_x = matching_x;
  *imr_matching_eccentricity = matching_e;
  *imr_matching_mean_ano = matching_l;
  *imr_matching_phi = matching_phi;
  *imr_matching_phi_dot = matching_phi_dot;
  *imr_matching_r = matching_r;
  *imr_matching_r_dot = matching_r_dot;

  (*time_evol)->data->length = (*x_evol)->data->length =
      (*eccentricity_evol)->data->length = (*mean_ano_evol)->data->length =
          (*phi_evol)->data->length = (*phi_dot_evol)->data->length =
              (*r_evol)->data->length = (*r_dot_evol)->data->length = Length;

  (*time_evol)->data->data = uniform_t_vec;
  (*x_evol)->data->data = uniform_x_vec;
  (*eccentricity_evol)->data->data = uniform_e_vec;
  (*mean_ano_evol)->data->data = uniform_l_vec;
  (*phi_evol)->data->data = uniform_phi_vec;
  (*phi_dot_evol)->data->data = uniform_phi_dot_vec;
  (*r_evol)->data->data = uniform_r_vec;
  (*r_dot_evol)->data->data = uniform_r_dot_vec;

  // Ownership of these has passed onto (*x)->data->data
  uniform_t_vec = uniform_x_vec = uniform_e_vec = uniform_l_vec =
      uniform_phi_vec = uniform_phi_dot_vec = uniform_r_vec =
          uniform_r_dot_vec = NULL;

XLAL_FAIL:
  if (phi_dot_interp)
    gsl_interp_free(phi_dot_interp);
  if (phi_dot_accel)
    gsl_interp_accel_free(phi_dot_accel);

  if (solver_evolve)
    gsl_odeiv_evolve_free(solver_evolve);
  if (solver_control)
    gsl_odeiv_control_free(solver_control);
  if (solver_step)
    gsl_odeiv_step_free(solver_step);

  LALFree(t_vec);
  LALFree(x_vec);
  LALFree(e_vec);
  LALFree(phi_vec);
  LALFree(l_vec);
  LALFree(u_vec);
  LALFree(phi_dot_vec);
  LALFree(p_vec);
  LALFree(r_vec);

  LALFree(uniform_t_vec);
  LALFree(uniform_x_vec);
  LALFree(uniform_e_vec);
  LALFree(uniform_phi_vec);
  LALFree(uniform_l_vec);
  LALFree(uniform_phi_dot_vec);
  LALFree(uniform_r_vec);
  LALFree(uniform_r_dot_vec);

  return errorcode;
}

/** This function evaluates the ENIGMA model to obtain the + and x
 * polarizations.*/
int XLALSimInspiralENIGMA(
    REAL8TimeSeries **hplus,  /**< OUTPUT h_+ vector */
    REAL8TimeSeries **hcross, /**< OUTPUT h_x vector */
    REAL8 phiRef,             /**< orbital phase at reference pt. */
    REAL8 inclination,        /**< inclination angle */
    REAL8 eccentricity,       /**< eccentricity at reference pt. */
    REAL8 meanPerAno,         /**< mean anomaly of periastron */
    REAL8 deltaT,             /**< sampling interval (s) */
    REAL8 m1,                 /**< mass of companion 1 (kg) */
    REAL8 m2,                 /**< mass of companion 2 (kg) */
    REAL8 S1z,                /**< z-component of the spin of companion 1 */
    REAL8 S2z,                /**< z-component of the spin of companion 2 */
    REAL8 distance,           /**< distance of source (m) */
    REAL8 fMin,               /**< start GW frequency (Hz) */
    REAL8 fRef                /**< reference GW frequency (Hz) */
) {
  (void)phiRef;
  (void)fRef;
  double fsamp = 1 / deltaT;
  double beta = 0.; // for now

  int errorcode = XLAL_SUCCESS;

  m1 /= LAL_MSUN_SI;
  m2 /= LAL_MSUN_SI;

  // check input parameters
  if (m1 <= 0.0 || m2 <= 0.0 || eccentricity < 0.0 || eccentricity >= 1.0) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }

  if (!hplus || !hcross) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  if (*hplus || *hcross) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  // the actually input tolerance (reset to 1.0e-16 for eccentricity=0.0)
  const REAL8 tol_in = eccentricity == 0.0 ? 1.0e-16 : 1e-12;

  LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
  *hplus = XLALCreateREAL8TimeSeries("h_plus", &epoch, fMin, 1 / fsamp,
                                     &lalStrainUnit, 0);
  *hcross = XLALCreateREAL8TimeSeries("h_cross", &epoch, fMin, 1 / fsamp,
                                      &lalStrainUnit, 0);
  if (!hplus || !hcross) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  XLAL_CALLGSL(errorcode = x_model_eccbbh_inspiral_waveform(
                   *hplus, *hcross, m1, m2, S1z, S2z, eccentricity, fMin,
                   distance, meanPerAno, tol_in, inclination, beta, fsamp));
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(XLAL_EFUNC);

  /* Set epoch to the length of the inspiral waveform, as is usual convention */
  XLALGPSSetREAL8(&(*hplus)->epoch, -deltaT * (REAL8)(*hplus)->data->length);
  XLALGPSSetREAL8(&(*hcross)->epoch, -deltaT * (REAL8)(*hplus)->data->length);

XLAL_FAIL:
  if (errorcode != XLAL_SUCCESS) {
    XLALDestroyREAL8TimeSeries(*hplus);
    XLALDestroyREAL8TimeSeries(*hcross);
    *hplus = *hcross = NULL;
  }

  return errorcode;
}

int XLALSimIMRENIGMA(REAL8TimeSeries **hplus,  /**< OUTPUT h_+ vector */
                     REAL8TimeSeries **hcross, /**< OUTPUT h_x vector */
                     REAL8 phiRef,       /**< orbital phase at reference pt. */
                     REAL8 inclination,  /**< inclination angle */
                     REAL8 eccentricity, /**< eccentricity at reference pt. */
                     REAL8 meanPerAno,   /**< mean anomaly of periastron */
                     REAL8 deltaT,       /**< sampling interval (s) */
                     REAL8 m1,           /**< mass of companion 1 (kg) */
                     REAL8 m2,           /**< mass of companion 2 (kg) */
                     REAL8 S1z, /**< z-component of the spin of companion 1 */
                     REAL8 S2z, /**< z-component of the spin of companion 2 */
                     REAL8 distance, /**< distance of source (m) */
                     REAL8 fMin,     /**< start GW frequency (Hz) */
                     REAL8 fRef      /**< reference GW frequency (Hz) */
) {
  (void)phiRef;
  (void)fRef;
  double fsamp = 1 / deltaT;
  double beta = 0.; // for now

  int errorcode = XLAL_SUCCESS;

  m1 /= LAL_MSUN_SI;
  m2 /= LAL_MSUN_SI;

  // check input parameters
  if (m1 <= 0.0 || m2 <= 0.0 || eccentricity < 0.0 || eccentricity >= 1.0) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR(errorcode);
  }

  if (!hplus || !hcross) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  if (*hplus || *hcross) {
    errorcode = XLAL_EINVAL;
    XLAL_ERROR_FAIL(errorcode);
  }

  // Initialise an empty training set
  static TrainingSet *Dset = NULL;
  if (Dset == NULL) { // TODO: make thread safe
    Dset = CreateEmptyTrainingSet();
    if (!Dset)
      XLAL_ERROR_FAIL(XLAL_EFAILED);
  }

  // the actually input tolerance (reset to 1.0e-16 for eccentricity=0.0)
  const REAL8 tol_in = eccentricity == 0.0 ? 1.0e-16 : 1e-12;

  LIGOTimeGPS epoch = LIGOTIMEGPSZERO;
  *hplus = XLALCreateREAL8TimeSeries("h_plus", &epoch, fMin, 1 / fsamp,
                                     &lalStrainUnit, 0);
  *hcross = XLALCreateREAL8TimeSeries("h_cross", &epoch, fMin, 1 / fsamp,
                                      &lalStrainUnit, 0);
  if (!hplus || !hcross) {
    errorcode = XLAL_ENOMEM;
    XLAL_ERROR_FAIL(XLAL_EFUNC);
  }

  XLAL_CALLGSL(errorcode = x_model_eccbbh_imr_waveform(
                   *hplus, *hcross, m1, m2, S1z, S2z, eccentricity, fMin,
                   distance, meanPerAno, tol_in, inclination, beta, fsamp,
                   Dset));
  if (errorcode != XLAL_SUCCESS)
    XLAL_ERROR_FAIL(XLAL_EFUNC);

XLAL_FAIL:
  if (errorcode != XLAL_SUCCESS) {
    XLALDestroyREAL8TimeSeries(*hplus);
    XLALDestroyREAL8TimeSeries(*hcross);
    *hplus = *hcross = NULL;
  }

  return errorcode;
}
