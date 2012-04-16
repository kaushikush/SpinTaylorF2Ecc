/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Kipp Cannon
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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#ifndef _REALFFT_H
#define _REALFFT_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#ifdef  __cplusplus
#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif
#else
#define RESTRICT restrict
#endif

/** Plan to perform FFT of REAL4 data */
typedef struct tagREAL4FFTPlan REAL4FFTPlan;
/** Plan to perform FFT of REAL8 data */
typedef struct tagREAL8FFTPlan REAL8FFTPlan;
#define tagRealFFTPlan tagREAL4FFTPlan
#define RealFFTPlan REAL4FFTPlan

/**
 *
 * \addtogroup RealFFT_h
 * \brief Performs real-to-complex and complex-to-real FFTs.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/RealFFT.h>
 * \endcode
 *
 * Perform real-to-complex and complex-to-real fast Fourier
 * transforms of vectors, and sequences of vectors using the package
 * FFTW [\ref fj_1998].
 *
 *
 * \section sec_RealFFT_XLAL XLAL Functions
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/RealFFT.h>
 *
 * REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl );
 * REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl );
 * REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl );
 * void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan );
 *
 * int XLALREAL4ForwardFFT( COMPLEX8Vector *output, REAL4Vector *input, REAL4FFTPlan *plan );
 * int XLALREAL4ReverseFFT( REAL4Vector *output, COMPLEX8Vector *input, REAL4FFTPlan *plan );
 * int XLALREAL4VectorFFT( REAL4Vector *output, REAL4Vector *input, REAL4FFTPlan *plan );
 * int XLALREAL4PowerSpectrum( REAL4Vector *spec, REAL4Vector *data, REAL4FFTPlan *plan );
 *
 * REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl );
 * REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl );
 * REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl );
 * void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan );
 *
 * int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input, REAL8FFTPlan *plan );
 * int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input, REAL8FFTPlan *plan );
 * int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input, REAL8FFTPlan *plan );
 * int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data, REAL8FFTPlan *plan );
 * \endcode
 *
 * \heading{Description}
 *
 * The \c REAL4 routines are described below.  These use single-precision
 * FFTs, i.e., they convert \c REAL4Vectors into \c COMPLEX8Vectors
 * and vice-versa.  The \c REAL8 versions of the routines are the same
 * but they are double-precision versions, i.e., they convert
 * \c REAL8Vectors into \c COMPLEX16Vectors.
 *
 * The routine XLALCreateREAL4FFTPlan() creates a REAL4FFTPlan
 * structure to perform FFTs of vectors of length \c size.  If
 * \c fwdflg is non-zero then the plan is created to perform forward
 * (real-to-complex) FFTs with a negative exponential sign.  Otherwise
 * the plan is created to perform reverse (complex-to-real) FFTs with a
 * positive exponential sign.  The value of \c measurelvl determines
 * how much optimization of the plan FFTW will do with the most optimization
 * taking the most amount of time.  Reasonable values for \c measurelvl
 * would be 0 for the fasted plan creation (FFTW does not measure the speed
 * of any transform with this level but rather estimates which plan will
 * be the fastet) or 1 to measure a few likely plans to determine the fastest.
 *
 * XLALCreateForwardREAL4FFTPlan() is equivalent to
 * XLALCreateREAL4FFTPlan() with \c fwdflg set to 1.
 * XLALCreateReverseREAL4FFTPlan() is equivalent to
 * XLALCreateREAL4FFTPlan() with \c fwdflg set to 0.
 *
 * XLALDestroyREAL4FFTPlan() is used to destroy the plan, freeing all
 * memory that was allocated in the structure as well as the structure
 * itself.  It can be used on either forward or reverse plans.
 *
 * XLALREAL4ForwardFFT() and
 * XLALREAL4ReverseFFT() perform forward (real to complex) and
 * reverse (complex to real) transforms respectively.  The plan supplied
 * to these routines must be correctly generated for the direction of the
 * transform.  I.e., XLALREAL4ForwardFFT() cannot be supplied with
 * a plan generated by XLALCreateReverseREAL4FFTPlan().
 *
 * XLALREAL4VectorFFT() is a low-level routine that transforms
 * a real vector to a half-complex real vector (with a forward plan) or
 * a half-complex real vector to a real vector (with a reverse plan).
 * If you're not sure what this means, don't use this routine.
 * The input and output vectors (and their data) must be distinct pointers.
 *
 * XLALREAL4PowerSpectrum() computes a real power spectrum of the
 * input real vector and a forward FFT plan.
 *
 * \heading{Return Values}
 *
 * Upon success,
 * XLALCreateREAL4FFTPlan(),
 * XLALCreateForwardREAL4FFTPlan(), and
 * XLALCreateReverseREAL4FFTPlan() return a pointer to a newly-allocated
 * FFT plan.  Upon failure, they return a \c NULL pointer and set
 * \c xlalErrno to one of the following values:
 * #XLAL_EBADLEN if \c size is not greater than zero,
 * #XLAL_ENOMEM if a memory allocation failed, or
 * #XLAL_EFAILED if the FFTW plan creation routine failed.
 *
 * XLALDestroyREAL4FFTPlan() does not return any value but, upon
 * failure, it will set \c xlalErrno to one of the following values:
 * #XLAL_EFAULT if the routine is provided a \c NULL pointer, or
 * #XLAL_EINVAL if the contents of the plan are invalid (e.g., if the
 * routine is provided a plan that had been previously destroyed).
 *
 * XLALREAL4ForwardFFT(),
 * XLALREAL4ReverseFFT(),
 * XLALREAL4VectorFFT(), and
 * XLALREAL4PowerSpectrum() return the value 0 upon succes; upon
 * failure they return #XLAL_FAILURE and set xlalErrno to
 * one of the following values:
 * #XLAL_EFAULT if one of the input pointers is \c NULL,
 * #XLAL_EINVAL if the input, output, or plan structures appears
 * invalid or if the routine is passed a plan for the wrong transform
 * directions or if the input and output data pointers are not distinct
 * for XLALREAL4VectorFFT(),
 * #XLAL_EBADLEN if the input and output vectors and the plan have
 * incompatible lengths,
 * #XLAL_ENOMEM if a memory allocation of temporary internal memory
 * fails.
 *
 * As before, the \c REAL8 versions of these routines behave the
 * same way but for double-precision transforms.
 *
*/
/*@{ */
/** \name Error Codes */ /*@{*/
#define REALFFTH_ENULL 1	/**< Null pointer */
#define REALFFTH_ENNUL 2	/**< Non-null pointer */
#define REALFFTH_ESIZE 4	/**< Invalid input size */
#define REALFFTH_ESZMM 8	/**< Size mismatch */
#define REALFFTH_ESLEN 16	/**< Invalid/mismatched sequence lengths */
#define REALFFTH_ESAME 32	/**< Input/Output data vectors are the same */
#define REALFFTH_ESIGN 64	/**< Incorrect plan sign */
#define REALFFTH_EDATA 128	/**< Bad input data: DC/Nyquist should be real */
#define REALFFTH_EALOC 256	/**< Memory allocation failed */
#define REALFFTH_EFFTW 512	/**< Error in FFTW */
#define REALFFTH_ESNGL 1024	/**< FFTW library is not single-precision */
#define REALFFTH_EINTL 2048	/**< Error in Intel FFT library */
/*@}*/

/** \cond DONT_DOXYGEN */
#define REALFFTH_MSGENULL "Null pointer"
#define REALFFTH_MSGENNUL "Non-null pointer"
#define REALFFTH_MSGESIZE "Invalid input size"
#define REALFFTH_MSGESZMM "Size mismatch"
#define REALFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define REALFFTH_MSGESAME "Input/Output data vectors are the same"
#define REALFFTH_MSGESIGN "Incorrect plan sign"
#define REALFFTH_MSGEDATA "Bad input data: DC/Nyquist should be real"
#define REALFFTH_MSGEALOC "Memory allocation failed"
#define REALFFTH_MSGEFFTW "Error in FFTW"
#define REALFFTH_MSGESNGL "FFTW library is not single-precision"
#define REALFFTH_MSGEINTL "Error in Intel FFT library"
/** \endcond */

/*
 *
 * XLAL REAL4 functions
 *
 */

/** Returns a new REAL4FFTPlan
 *
 * A REAL4FFTPlan is required to perform a FFT that involves real data.
 * A different plan is required for each size of the real data vector
 * and for each direction of transform (forward or reverse).
 * A forward transform performs
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 * A reverse transform performs
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a \"reverse\"
 * transform rather than an \"inverse\" transform).
 *
 * @param[in] size The number of points in the real data.
 * @param[in] fwdflg Set non-zero for a forward FFT plan;
 *                    otherwise create a reverse plan
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL4FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateREAL4Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl );

/** Returns a new REAL4FFTPlan for a forward transform
 *
 * A REAL4FFTPlan is required to perform a FFT that involves real data.
 * A different plan is required for each size of the real data vector.
 * A forward transform performs
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 *
 * @param[in] size The number of points in the real data.
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL4FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateForwardREAL4Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl );

/** Returns a new REAL4FFTPlan for a reverse transform
 *
 * A REAL4FFTPlan is required to perform a FFT that involves real data.
 * A reverse transform performs
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a \"reverse\"
 * transform rather than an \"inverse\" transform).
 *
 * @param[in] size The number of points in the real data.
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL4FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateReverseREAL4Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl );

/** Destroys a REAL4FFTPlan
 * @param[in] plan A pointer to the REAL4FFTPlan to be destroyed.
 * @return None.
 */
void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan );

/** Performs a forward FFT of REAL4 data
 *
 * This routine performs the transformation:
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 *
 * @note
 * Due to the reality of the input data x, the following identity
 * holds for the complex FFT data: \f$z[N-k]=z^\ast[k]\f$.  Therefore,
 * the length of the output vector is equal to \f$\lfloor N/2\rfloor + 1\f$
 * since the remaining \"negative\" frequency components can be obtained from the
 * \"positive\" frequency components.
 *
 * @param[out] output The complex data vector z of length [N/2] + 1
 * that results from the transform
 * @param[in] input The real data vector x of length to be transformed
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL4ForwardFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL4ForwardFFT( COMPLEX8Vector *output, const REAL4Vector *input,
    const REAL4FFTPlan *plan );

/** Performs a reverse FFT of REAL4 data
 *
 * This routine performs the transformation:
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * Due to the reality of the output data y, the following identity
 * holds for the complex data: \f$z[N-k]=z^\ast[k]\f$.  Therefore,
 * the length of the input vector is equal to \f$\lfloor N/2\rfloor + 1\f$
 * since the remaining \"negative\" frequency components can be obtained from the
 * \"positive\" frequency components.
 *
 * @param[out] output The real data vector y of length N
 * that results from the transform
 * @param[in] input The complex data vector z of length [N/2] + 1
 * to be transformed
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL4ForwardFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EDOM] Domain error if the DC component of the input data, z[0],
 *    is not purely real
 *    or if the length of the output vector N is even and the Nyquist
 *    component of the input data, z[N/2], is not purely real.
 *  .
 */
int XLALREAL4ReverseFFT( REAL4Vector *output, const COMPLEX8Vector *input,
    const REAL4FFTPlan *plan );

/** Perform a REAL4Vector to REAL4Vector FFT
 *
 * This routine computes
 * \f[y[k]=\left\{\begin{array}{ll}\Re z[k]&0\le k\le\lfloor N/2\rfloor\\\Im z[N-k]&\lfloor N/2\rfloor<k<N\end{array}\right.\f]
 * where \f[z[k] = \sum_{j=0}^{N-1} e^{\mp2\pi ijk/N}\,x[j],\f]
 * and where the minus sign is used if a forward plan is provided as the argument
 * and the plus sign is used if a reverse plan is provided as the argument;
 * here N is the length of the input vector x.
 *
 * @param[out] output The real output data vector y of length N
 * @param[in] input The input real data vector x of length N
 * @param[in] plan The FFT plan to use for the transform
 * @note
 * The input and output vectors must be distinct.
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL4VectorFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL4VectorFFT( REAL4Vector * RESTRICT output, const REAL4Vector * RESTRICT input,
    const REAL4FFTPlan *plan );

/** Computes the power spectrum of REAL4 data
 *
 * This routine computes
 * \f[P[k]=\left\{\begin{array}{ll}|z[0]|^2&k=0\\2|z[k]|^2&1\leq \lfloor (N+1)/2\rfloor\\|z[N/2]|^2&k=N/2,\;\mbox{$N$ even}\end{array}\right.\f]
 * where \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j],\f]
 * and N is the length of the input vector x.
 *
 * @param[out] spec The real power spectrum P of length [N/2] + 1 of the data x
 * @param[in] data The input real data vector x of length N
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL4PowerSpectrum() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the input and output
 *    data vectors are the same.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL4PowerSpectrum( REAL4Vector *spec, const REAL4Vector *data,
    const REAL4FFTPlan *plan );

/*
 *
 * XLAL REAL8 functions
 *
 */

/** Returns a new REAL8FFTPlan
 *
 * A REAL8FFTPlan is required to perform a FFT that involves real data.
 * A different plan is required for each size of the real data vector
 * and for each direction of transform (forward or reverse).
 * A forward transform performs
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 * A reverse transform performs
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a \"reverse\"
 * transform rather than an \"inverse\" transform).
 *
 * @param[in] size The number of points in the real data.
 * @param[in] fwdflg Set non-zero for a forward FFT plan;
 *                    otherwise create a reverse plan
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateREAL8Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl );

/** Returns a new REAL8FFTPlan for a forward transform
 *
 * A REAL8FFTPlan is required to perform a FFT that involves real data.
 * A different plan is required for each size of the real data vector.
 * A forward transform performs
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 *
 * @param[in] size The number of points in the real data.
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateForwardREAL8Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl );

/** Returns a new REAL8FFTPlan for a reverse transform
 *
 * A REAL8FFTPlan is required to perform a FFT that involves real data.
 * A reverse transform performs
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * The reverse transform of the forward transform of some data is
 * equal to N times the original data (we therefore call it a \"reverse\"
 * transform rather than an \"inverse\" transform).
 *
 * @param[in] size The number of points in the real data.
 * @param[in] measurelvl Measurement level for plan creation:
 *                      - 0: no measurement, just estimate the plan;
 *                      - 1: measure the best plan;
 *                      - 2: perform a lengthy measurement of the best plan;
 *                      - 3: perform an exhasutive measurement of the best plan.
 * @return A pointer to an allocated \c REAL8FFTPlan structure is returned
 * upon successful completion.  Otherwise, a \c NULL pointer is returned
 * and \c xlalErrno is set to indicate the error.
 * @par Errors:
 * The \c XLALCreateReverseREAL8Plan() function shall fail if:
 *  - [\c XLAL_EBADLEN] The size of the requested plan is 0.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EFAILED] The call to the underlying FFTW routine failed.
 *  .
 */
REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl );

/** Destroys a REAL8FFTPlan
 * @param[in] plan A pointer to the REAL8FFTPlan to be destroyed.
 * @return None.
 */
void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan );

/** Performs a forward FFT of REAL8 data
 *
 * This routine performs the transformation:
 * \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j]\f]
 * where N, the size of the transform, is the length of the vector x.
 *
 * @note
 * Due to the reality of the input data x, the following identity
 * holds for the complex FFT data: \f$z[N-k]=z^\ast[k]\f$.  Therefore,
 * the length of the output vector is equal to \f$\lfloor N/2\rfloor + 1\f$
 * since the remaining \"negative\" frequency components can be obtained from the
 * \"positive\" frequency components.
 *
 * @param[out] output The complex data vector z of length [N/2] + 1
 * that results from the transform
 * @param[in] input The real data vector x of length to be transformed
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL8ForwardFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan );

/** Performs a reverse FFT of REAL8 data
 *
 * This routine performs the transformation:
 * \f[y[j] = \sum_{k=0}^{N-1} e^{+2\pi ijk/N}\,z[k]\f]
 * where N, the size of the transform, is the length of the vector y.
 *
 * @note
 * Due to the reality of the output data y, the following identity
 * holds for the complex data: \f$z[N-k]=z^\ast[k]\f$.  Therefore,
 * the length of the input vector is equal to \f$\lfloor N/2\rfloor + 1\f$
 * since the remaining \"negative\" frequency components can be obtained from the
 * \"positive\" frequency components.
 *
 * @param[out] output The real data vector y of length N
 * that results from the transform
 * @param[in] input The complex data vector z of length [N/2] + 1
 * to be transformed
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL8ForwardFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  - [\c XLAL_EDOM] Domain error if the DC component of the input data, z[0],
 *    is not purely real
 *    or if the length of the output vector N is even and the Nyquist
 *    component of the input data, z[N/2], is not purely real.
 *  .
 */
int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input,
    const REAL8FFTPlan *plan );

/** Perform a REAL8Vector to REAL8Vector FFT
 *
 * This routine computes
 * \f[y[k]=\left\{\begin{array}{ll}\Re z[k]&0\le k\le\lfloor N/2\rfloor\\\Im z[N-k]&\lfloor N/2\rfloor<k<N\end{array}\right.\f]
 * where \f[z[k] = \sum_{j=0}^{N-1} e^{\mp2\pi ijk/N}\,x[j],\f]
 * and where the minus sign is used if a forward plan is provided as the argument
 * and the plus sign is used if a reverse plan is provided as the argument;
 * here N is the length of the input vector x.
 *
 * @param[out] output The real output data vector y of length N
 * @param[in] input The input real data vector x of length N
 * @param[in] plan The FFT plan to use for the transform
 * @note
 * The input and output vectors must be distinct.
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL8VectorFFT() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the input and output data
 *    vectors are the same.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input,
    const REAL8FFTPlan *plan );

/** Computes the power spectrum of REAL8 data
 *
 * This routine computes
 * \f[P[k]=\left\{\begin{array}{ll}|z[0]|^2 & k=0\\2|z[k]|^2 & 1\leq \lfloor (N+1)/2\rfloor\\ |z[N/2]|^2 & k=N/2,\;\mbox{$N$ even}\end{array}\right.\f]
 * where \f[z[k] = \sum_{j=0}^{N-1} e^{-2\pi ijk/N}\,x[j],\f]
 * and N is the length of the input vector x.
 *
 * @param[out] spec The real power spectrum P of length [N/2] + 1 of the data x
 * @param[in] data The input real data vector x of length N
 * @param[in] plan The FFT plan to use for the transform
 * @return 0 upon successful completion or non-zero upon failure.
 * @par Errors:
 * The \c XLALREAL8PowerSpectrum() function shall fail if:
 *  - [\c XLAL_EFAULT] A \c NULL pointer is provided as one of the arguments.
 *  - [\c XLAL_EINVAL] A argument is invalid or the plan is for a
 *    reverse transform.
 *  - [\c XLAL_EBADLEN] The input vector, output vector, and plan size are
 *    incompatible.
 *  - [\c XLAL_ENOMEM] Insufficient storage space is available.
 *  .
 */
int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data,
    const REAL8FFTPlan *plan );

/*
 *
 * LAL REAL4 functions
 *
 */

/** \b DEPRECATED
 * @deprecated Use XLALCreateForwardREAL4FFTPlan() instead.
 */
void
LALCreateForwardREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );
/** \b DEPRECATED
 * @deprecated Use XLALCreateForwardREAL4FFTPlan() instead.
 */
#define LALCreateForwardRealFFTPlan LALCreateForwardREAL4FFTPlan

/** \b DEPRECATED
 * @deprecated Use XLALCreateReverseREAL4FFTPlan() instead.
 */
void
LALCreateReverseREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );
/** \b DEPRECATED
 * @deprecated Use XLALCreateReverseREAL4FFTPlan() instead.
 */
#define LALCreateReverseRealFFTPlan LALCreateReverseREAL4FFTPlan

/** \b DEPRECATED
 * @deprecated Use XLALDestroyREAL4FFTPlan() instead.
 */
void
LALDestroyREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan
    );
/** \b DEPRECATED
 * @deprecated Use XLALDestroyREAL4FFTPlan() instead.
 */
#define LALDestroyRealFFTPlan LALDestroyREAL4FFTPlan

/** \b DEPRECATED
 * @deprecated Use XLALREAL4ForwardFFT() instead.
 */
void
LALForwardREAL4FFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    REAL4FFTPlan    *plan
    );
/** \b DEPRECATED
 * @deprecated Use XLALREAL4ForwardFFT() instead.
 */
#define LALForwardRealFFT LALForwardREAL4FFT

/** \b DEPRECATED
 * @deprecated Use XLALREAL4ReverseFFT() instead.
 */
void
LALReverseREAL4FFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    REAL4FFTPlan    *plan
    );
/** \b DEPRECATED
 * @deprecated Use XLALREAL4ReverseFFT() instead.
 */
#define LALReverseRealFFT LALReverseREAL4FFT

/** \b DEPRECATED
 * @deprecated Use XLALREAL4PowerSpectrum() instead.
 */
void
LALREAL4PowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    REAL4FFTPlan *plan
    );
/** \b DEPRECATED
 * @deprecated Use XLALREAL4PowerSpectrum() instead.
 */
#define LALRealPowerSpectrum LALREAL4PowerSpectrum

/** \b DEPRECATED
 * @deprecated Use XLALREAL4VectorFFT() instead.
 */
void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    REAL4FFTPlan *plan
    );

/*
 *
 * LAL REAL8 functions
 *
 */

/** \b DEPRECATED
 * @deprecated Use XLALCreateForwardREAL8FFTPlan() instead.
 */
void
LALCreateForwardREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

/** \b DEPRECATED
 * @deprecated Use XLALCreateReverseREAL8FFTPlan() instead.
 */
void
LALCreateReverseREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

/** \b DEPRECATED
 * @deprecated Use XLALDestroyREAL8FFTPlan() instead.
 */
void
LALDestroyREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan
    );

/** \b DEPRECATED
 * @deprecated Use XLALREAL8ForwardFFT() instead.
 */
void
LALForwardREAL8FFT(
    LALStatus      *status,
    COMPLEX16Vector *output,
    REAL8Vector    *input,
    REAL8FFTPlan    *plan
    );

/** \b DEPRECATED
 * @deprecated Use XLALREAL8ReverseFFT() instead.
 */
void
LALReverseREAL8FFT(
    LALStatus      *status,
    REAL8Vector    *output,
    COMPLEX16Vector *input,
    REAL8FFTPlan    *plan
    );

/** \b DEPRECATED
 * @deprecated Use XLALREAL8PowerSpectrum() instead.
 */
void
LALREAL8PowerSpectrum (
    LALStatus   *status,
    REAL8Vector *spec,
    REAL8Vector *data,
    REAL8FFTPlan *plan
    );

/** \b DEPRECATED
 * @deprecated Use XLALREAL8VectorFFT() instead.
 */
void
LALREAL8VectorFFT(
    LALStatus   *status,
    REAL8Vector *output,
    REAL8Vector *input,
    REAL8FFTPlan *plan
    );

/*@}*/

#undef RESTRICT

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _REALFFT_H */
