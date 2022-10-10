/*
 * Copyright (C) 2020 Rodrigo Tenorio
 * Copyright (C) 2020 David Keitel
 * Copyright (C) 2019 Pep Covas
 * Copyright (C) 2010--2016, 2019--2022 Karl Wette
 * Copyright (C) 2009 Chris Messenger
 * Copyright (C) 2004--2006, 2008--2011, 2013, 2014, 2018 Reinhard Prix
 * Copyright (C) 2004, 2005 Bernd Machenschalk
 * Copyright (C) 2004, 2005 Alicia Sintes
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

#ifndef _SFTFILEIO_H  	/* Double-include protection. */
#define _SFTFILEIO_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/* includes */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/StringVector.h>
#include <lal/Date.h>
#include <lal/Segments.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALRunningMedian.h>
#include <lal/RngMedBias.h>
#include <lal/UserInputParse.h>
#include <lal/PulsarDataTypes.h>

/**
 * \defgroup SFTfileIO_h Header SFTfileIO.h
 * \ingroup lalpulsar_sft
 * \author R. Prix, B. Machenschalk, A.M. Sintes, B. Krishnan, I. Gholami,
 *         C. Messenger, K. Wette, P. Covas, D. Keitel, R. Tenorio
 *
 * \brief Module for reading/writing/manipulating SFTs (Short Fourier Transforms)
 *
 * This implements the SFT (Short Fourier Transforms) version 2 standard defined
 * in \cite SFT-spec . It contains various helper functions to create, handle,
 * combine, and destroy SFTs and related data structures, including SFTtype,
 * SFTVector, SFTCatalog (and their multi-detector generalizations) as well as
 * tools for dealing with timestamps, segments and ASD/PSD (Amplitude/Power
 * Spectral Density) estimates.
 *
 * <p>
 * <h2>Overview</h2>
 * - \ref Time-freq-convention-func "Time range/frequency bin convention functions".
 * - SFT types:
 *   \ref SFT-type-cdtor-func "create/destroy functions",
 *   \ref SFT-type-mod-func "modify functions",
 *   \ref SFT-type-prop-func "property functions".
 * - SFT timestamp types:
 *   \ref SFT-timestamp-cdtor-func "create/destroy functions",
 *   \ref SFT-timestamp-prop-func "property functions",
 *   \ref SFT-timestamp-gen-func "generation functions".
 * - SFT catalog types:
 *   \ref SFT-catalog-cdtor-func "create/destroy functions",
 *   \ref SFT-catalog-prop-func "property functions",
 *   \ref SFT-catalog-gen-func "generation functions".
 * - SFT files:
 *   \ref SFT-file-naming-func "file naming convention functions",
 *   \ref SFT-file-read-func "file reading functions",
 *   \ref SFT-file-write-func "file writing functions".
 * - PSD types:
 *   \ref PSD-type-cdtor-func "create/destroy functions",
 *   \ref PSD-type-mod-func "modify functions",
 *   \ref PSD-type-prop-func "property functions",
 *   \ref PSD-type-gen-func "generation functions",
 *   \ref PSD-file-write-func "file writing functions".
 *
 * <p>
 * <h2>Usage: Reading of SFT-files</h2>
 *
 * The basic operation of <b>reading SFTs</b> from files proceeds in two simple steps:
 *
 * - XLALSFTdataFind(): get an SFTCatalog of SFTs matching certain requirements (free with XLALDestroySFTCatalog())
 * - XLALLoadSFTs(): load a frequency-band into a single-IFO SFTVector defined by the catalogue, OR <br>
 *   XLALLoadMultiSFTs(): load a frequency-band into a MultiSFTVector defined by the catalogue
 *
 * <b>Note 1:</b> currently supported SFT file-formats are (merged or single) SFT-v2 files.
 * This might be extended in the future to support further file-formats (frames?).
 * None of the following API depends on the details of the underlying file-format. This will ensure that
 * codes using the following functions will NOT have to be changed irrespective of SFT file-format used.
 *
 * <b>Note 2:</b> irrespective of the underlying SFT file-format, the returned SFTs (SFTVector) will
 * <em>ALWAYS</em> be normalized according the the LAL-specification for frequency-series
 * (<tt>LIGO-T010095-00</tt>), that is the pure DFT of the time-series \f$x_j\f$ is <em>multiplied</em>
 * by the time-step \f$\Delta t\f$:
 * \f[
 * \mathrm{data}[k] = X^\mathrm{d}_k = \Delta t \,\sum_{j=0}^{N-1} x_j \,e^{-i2\pi \,k \,j / N}
 * \f]
 *
 * <h4>Details to 1: find matching SFTs and get the SFTCatalog:</h4>
 *
 * Thes function XLALSFTdataFind() returns an SFTCatalog of matching SFTs for a given file-pattern
 * (e.g. "SFT.*", "SFT.000", "/some/path/some_files_[0-9]?.sft", etc ) and additional, optional SFTConstraints.
 *
 * The optional constraints are:
 * - detector-prefix (e.g. "H1", "H2", "L1", "G1", "V1", etc..)
 * - GPS start-time + end-time
 * - a list of GPS-timestamps
 *
 * <b>Note 1:</b> Any constraint can be specified as \c NULL, all given constraints will be
 * combined by logical \c AND.
 *
 * <b>Note 2:</b> if a timestamps-list is given, *ALL* timestamps within
 * <tt>[minStartTime, maxStartTime)</tt> MUST be found!]
 *
 * <b>Note 3:</b> XLALSFTdataFind() will refuse to return any SFTs without their detector-name
 * properly set.
 *
 * The returned SFTCatalog is a vector of SFTDescriptor describing one SFT, with the fields
 * - \c locator:  an opaque data-type describing where to read this SFT from.
 * - \c header:	the SFts header
 * - \c comment: the comment-string found in the SFT, if any
 * - \c numBins: the number of frequency-bins in the SFT
 * - \c version: version-number of SFT file-format
 * - \c crc64: the crc64 checksum reported by this SFT
 *
 * One can use the following catalog-handling API functions:
 * - XLALDestroySFTCatalog(): free up a complete SFTCatalog
 * - XLALSFTtimestampsFromCatalog(): extract the list of SFT timestamps found in the SFTCatalog
 * - XLALDestroyTimestampVector(): free up a timestamps-vector (\c LIGOTimeGPSVector)
 * - XLALshowSFTLocator(): [*debugging only*] show a static string describing the 'locator'
 *
 * <b>NOTE:</b> The SFTs in the returned catalogue are \em guaranteed to
 * - be sorted in order of increasing GPS-epoch
 * - contain a valid detector-name
 *
 * <h4>Details to 2: load frequency-band from SFTs described in an SFTCatalog</h4>
 *
 * The function XLALLoadSFTs() takes an SFTCatalog and reads the smallest frequency-band containing <tt>[fMin, fMax]</tt>
 * from the SFTs, returning the resulting SFTVector. Note that this function will return an error if the
 * SFTCatalog contains SFTs from different detectors, for which XLALLoadMultiSFTs() must be used.
 *
 * The frequency-bounds are optional and \c -1 can be used to specify an 'open bound', i.e.<br>
 * <tt>[-1, fMax]</tt>: read from first frequency-bin in the SFT up to \c fMax.<br>
 * <tt>[fMin, -1]</tt>: read from \c fMin up to last frequency-bin in the SFTS<br>
 * <tt>[-1, -1]</tt>: read ALL frequency-bins from SFT.
 *
 * The function XLALLoadMultiSFTs() is similar to the above, except that it accepts an SFTCatalog with different detectors,
 * and returns corresponding multi-IFO vector of SFTVectors.
 *
 * <p><h2>Usage: Writing of SFT-files</h2>
 *
 * For <b>writing SFTs</b>:
 * - XLALWriteSFT2file(): write a single SFT (SFTtype) into an SFT-file following the specification \cite SFT-spec .
 */

/** @{ */

/*---------- exported types ----------*/

/** A so-called 'SFT' (short-Fourier-transform) will be stored in a COMPLEX8FrequencySeries */
typedef COMPLEX8FrequencySeries 	SFTtype;

/** The corresponding vector-type to hold a vector of 'SFTs' */
typedef COMPLEX8FrequencySeriesVector 	SFTVector;

/** A collection of SFT vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiSFTVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiSFTVector, SFTVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  SFTVector  **data; 	/**< sftvector for each ifo */
} MultiSFTVector;

/** A vector of 'timestamps' of type LIGOTimeGPS */
typedef struct tagLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(LIGOTimeGPSVector, LIGOTimeGPS, data, UINT4, length));
#endif /* SWIG */
  UINT4 	length;		/**< number of timestamps */
  LIGOTimeGPS 	*data;		/**< array of timestamps */
  REAL8		deltaT;		/**< 'length' of each timestamp (e.g. typically Tsft) */
} LIGOTimeGPSVector;

/** A collection of (multi-IFO) LIGOTimeGPSVector time-stamps vectors */
typedef struct tagMultiLIGOTimeGPSVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiLIGOTimeGPSVector, LIGOTimeGPSVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 	        length;	   /**< number of timestamps vectors or ifos */
  LIGOTimeGPSVector 	**data;    /**< timestamps vector for each ifo */
} MultiLIGOTimeGPSVector;

/**
 * 'Constraints' for SFT-matching: which detector, within which time-stretch and which
 * timestamps exactly should be loaded ?
 * Any of the entries is optional, and they will be combined by logical AND.
 * Note however, that *ALL* timestamps within [minStartTime, maxStartTime) MUST be found if specified.
 */
typedef struct tagSFTConstraints
{
  CHAR *detector;			/**< 2-char channel-prefix describing the detector (eg 'H1', 'H2', 'L1', 'G1' etc) */
  LIGOTimeGPS *minStartTime;		/**< only include SFTs whose epoch is >= minStartTime */
  LIGOTimeGPS *maxStartTime;		/**< only include SFTs whose epoch is <  maxStartTime */
  LIGOTimeGPSVector *timestamps;	/**< list of timestamps  */
} SFTConstraints;

/**
 * A 'descriptor' of an SFT: basically containing the header-info plus an opaque description
 * of where exactly to load this SFT from.
 */
typedef struct tagSFTDescriptor
{
  struct tagSFTLocator *locator; 	/**< *internal* description of where to find this SFT [opaque!] */
  SFTtype header;			/**< SFT-header info */
  CHAR *comment;			/**< comment-entry in SFT-header (v2 only) */
  UINT4 numBins;			/**< number of frequency-bins in this SFT */
  UINT4 version;			/**< SFT-specification version */
  UINT8 crc64;				/**< crc64 checksum */
} SFTDescriptor;

/** An "SFT-catalogue": a vector of SFTdescriptors, as returned by XLALSFTdataFind() */
typedef struct tagSFTCatalog
{
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(SFTCatalog, SFTDescriptor, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of SFTs in catalog */
  SFTDescriptor *data;		/**< array of data-entries describing matched SFTs */
} SFTCatalog;

/**
 * A multi-SFT-catalogue "view": a multi-IFO vector of SFT-catalogs
 *
 * Note: this is only a multi-IFO "view" of an existing SFTCatalog,
 * various allocated memory of the original catalog is only
 * pointed to, not duplicated!
 * This means one must not free the original catalog
 * while this multi-view is still in use!
 */
typedef struct tagMultiSFTCatalogView
{
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiSFTCatalogView, SFTCatalog, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;			/**< number of detectors */
  SFTCatalog *data;		/**< array of SFT-catalog pointers */
} MultiSFTCatalogView;

/** Special type for holding a PSD vector (over several SFTs) */
typedef REAL8FrequencySeriesVector PSDVector;

/** A collection of PSD vectors -- one for each IFO in a multi-IFO search */
typedef struct tagMultiPSDVector {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiPSDVector, PSDVector*, data, UINT4, length));
#endif /* SWIG */
  UINT4      length;  	/**< number of ifos */
  PSDVector  **data; 	/**< sftvector for each ifo */
} MultiPSDVector;

/** One noise-weight (number) per SFT (therefore indexed over IFOs and SFTs */
typedef struct tagMultiNoiseWeights {
#ifdef SWIG /* SWIG interface directives */
  SWIGLAL(ARRAY_1D(MultiNoiseWeights, REAL8Vector*, data, UINT4, length));
#endif /* SWIG */
  UINT4 length;		/**< number of detectors */
  REAL8Vector **data;	/**< weights-vector for each detector */
  REAL8 Sinv_Tsft;	/**< normalization factor used: \f$\mathcal{S}^{-1}\,T_\mathrm{SFT}\f$ (using single-sided PSD!) */
  BOOLEAN isNotNormalized;  /**< if true: weights are saved unnormalized (divide by Sinv_Tsft to get normalized version). */
} MultiNoiseWeights;

/** common types of mathematical operations over an array */
typedef enum tagMathOpType {
  MATH_OP_ARITHMETIC_SUM = 0,   /**< \f$\sum_k x_k\f$      */
  MATH_OP_ARITHMETIC_MEAN,      /**< \f$\sum_k x_k / N\f$ */
  MATH_OP_ARITHMETIC_MEDIAN,    /**< \f$x_1 \leq \dots \leq  x_{N/2} \leq \dots \leq x_n\f$ */
  MATH_OP_HARMONIC_SUM,         /**< \f$1 / \sum_k (1/x_k)\f$ */
  MATH_OP_HARMONIC_MEAN,        /**< \f$N / \sum_k (1/x_k)\f$ */
  MATH_OP_POWERMINUS2_SUM,      /**< \f$1 / \sqrt{ \sum_k (1/x_k^2) }\f$ */
  MATH_OP_POWERMINUS2_MEAN,     /**< \f$1 / \sqrt{ \sum_k (1/x_k^2) / N }\f$ */
  MATH_OP_MINIMUM,              /**< \f$\min_k(x_k)\f$ */
  MATH_OP_MAXIMUM,              /**< \f$\max_k(x_k)\f$ */
  MATH_OP_LAST
} MathOpType;

/*---------- global variables ----------*/

extern const UserChoices MathOpTypeChoices;

/*---------- exported prototypes [API] ----------*/

/**
 * \name Time range/frequency bin convention functions
 * \anchor Time-freq-convention-func
 */
/** @{ */

int XLALCWGPSinRange( const LIGOTimeGPS gps, const LIGOTimeGPS* minGPS, const LIGOTimeGPS* maxGPS );

UINT4 XLALRoundFrequencyDownToSFTBin( const REAL8 freq, const REAL8 df );
UINT4 XLALRoundFrequencyUpToSFTBin( const REAL8 freq, const REAL8 df );

int XLALFindCoveringSFTBins ( UINT4 *firstBin, UINT4 *numBins, REAL8 fMinIn, REAL8 BandIn, REAL8 Tsft );

/** @} */

/**
 * \name SFT type create/destroy functions
 * \anchor SFT-type-cdtor-func
 */
/** @{ */

SFTtype* XLALCreateSFT ( UINT4 numBins );
void XLALDestroySFT (SFTtype *sft);
int XLALCopySFT ( SFTtype *dest, const SFTtype *src );

SFTVector* XLALCreateSFTVector (UINT4 numSFTs, UINT4 numBins );
SFTVector* XLALCreateEmptySFTVector (UINT4 numSFTs );
void XLALDestroySFTVector (SFTVector *vect);
SFTVector *XLALDuplicateSFTVector ( const SFTVector *sftsIn );

MultiSFTVector *XLALCreateMultiSFTVector ( UINT4 length, UINT4Vector *numsft );
MultiSFTVector *XLALCreateEmptyMultiSFTVector ( UINT4Vector *numsft );
void XLALDestroyMultiSFTVector ( MultiSFTVector *multvect );

int XLALExtractBandFromSFT ( SFTtype **outSFT, const SFTtype *inSFT, REAL8 fMin, REAL8 Band );
SFTVector *XLALExtractBandFromSFTVector ( const SFTVector *inSFTs, REAL8 fMin, REAL8 Band );
MultiSFTVector *XLALExtractBandFromMultiSFTVector ( const MultiSFTVector *inSFTs, REAL8 fMin, REAL8 Band );

int XLALExtractStrictBandFromSFT ( SFTtype **outSFT, const SFTtype *inSFT, REAL8 fMin, REAL8 Band );
SFTVector *XLALExtractStrictBandFromSFTVector ( const SFTVector *inSFTs, REAL8 fMin, REAL8 Band );
MultiSFTVector *XLALExtractStrictBandFromMultiSFTVector ( const MultiSFTVector *inSFTs, REAL8 fMin, REAL8 Band );

/** @} */

/**
 * \name SFT type modify functions
 * \anchor SFT-type-mod-func
 */
/** @{ */

int XLALAppendSFT2Vector (SFTVector *vect, const SFTtype *sft );

int XLALReorderMultiSFTVector( MultiSFTVector *multiSFTs, const LALStringVector *IFOs);

int XLALSFTAdd ( SFTtype *a, const SFTtype *b );
int XLALSFTVectorAdd ( SFTVector *a, const SFTVector *b );
int XLALMultiSFTVectorAdd ( MultiSFTVector *a, const MultiSFTVector *b );

int XLALSFTResizeBand ( SFTtype *SFT, REAL8 f0, REAL8 Band );
int XLALSFTVectorResizeBand ( SFTVector *SFTs, REAL8 f0, REAL8 Band );
int XLALMultiSFTVectorResizeBand ( MultiSFTVector *multiSFTs, REAL8 f0, REAL8 Band );

/** @} */

/**
 * \name SFT type property functions
 * \anchor SFT-type-prop-func
 */
/** @{ */

int XLALEarliestMultiSFTsample ( LIGOTimeGPS *out, const MultiSFTVector *multisfts );
int XLALLatestMultiSFTsample ( LIGOTimeGPS *out, const MultiSFTVector *multisfts );

SFTVector *XLALExtractSFTVectorWithTimestamps ( const SFTVector *sfts, const LIGOTimeGPSVector *timestamps );
MultiSFTVector *XLALExtractMultiSFTVectorWithMultiTimestamps ( const MultiSFTVector *multiSFTs, const MultiLIGOTimeGPSVector *multiTimestamps );

/** @} */

/**
 * \name SFT timestamp type create/destroy functions
 * \anchor SFT-timestamp-cdtor-func
 */
/** @{ */

LIGOTimeGPSVector *XLALCreateTimestampVector (UINT4 len);
void XLALDestroyTimestampVector (LIGOTimeGPSVector *vect);
LIGOTimeGPSVector *XLALResizeTimestampVector ( LIGOTimeGPSVector *vector, UINT4 length );

MultiLIGOTimeGPSVector *XLALCreateMultiLIGOTimeGPSVector ( UINT4 numDetectors );
void XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS );

/** @} */

/**
 * \name SFT timestamp type property functions
 * \anchor SFT-timestamp-prop-func
 */
/** @{ */

int XLALFindTimesliceBounds ( UINT4 *iStart, UINT4 *iEnd, const LIGOTimeGPSVector *timestamps, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS );

/** @} */

/**
 * \name SFT timestamp generation functions
 * \anchor SFT-timestamp-gen-func
 */
/** @{ */

LIGOTimeGPSVector *XLALMakeTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap );
MultiLIGOTimeGPSVector *XLALMakeMultiTimestamps ( LIGOTimeGPS tStart, REAL8 Tspan, REAL8 Tsft, REAL8 Toverlap, UINT4 numDet );

LIGOTimeGPSVector *XLALReadTimestampsFile ( const CHAR *fname );
MultiLIGOTimeGPSVector *XLALReadMultiTimestampsFiles ( const LALStringVector *fnames );

LIGOTimeGPSVector *XLALReadTimestampsFileConstrained ( const CHAR *fname, const LIGOTimeGPS *minGPS, const LIGOTimeGPS *maxGPS );
MultiLIGOTimeGPSVector *XLALReadMultiTimestampsFilesConstrained ( const LALStringVector *fnames, const LIGOTimeGPS *minGPS, const LIGOTimeGPS *maxGPS );

LIGOTimeGPSVector *XLALExtractTimestampsFromSFTs ( const SFTVector *sfts );
MultiLIGOTimeGPSVector *XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs );

LIGOTimeGPSVector *XLALTimestampsFromSFTCatalog ( const SFTCatalog *catalog );
MultiLIGOTimeGPSVector *XLALTimestampsFromMultiSFTCatalogView ( const MultiSFTCatalogView *multiView );

LALSegList *XLALReadSegmentsFromFile ( const char *fname );

LIGOTimeGPSVector *XLALTimestampsFromSegmentFile( const char *filename, REAL8 Tsft, REAL8 Toverlap, BOOLEAN adjustSegExtraTime, BOOLEAN synchronize );

/** @} */

/**
 * \name SFT catalog create/destroy functions
 * \anchor SFT-catalog-cdtor-func
 */
/** @{ */

SFTCatalog *XLALSFTdataFind ( const CHAR *file_pattern, const SFTConstraints *constraints );
void XLALDestroySFTCatalog ( SFTCatalog *catalog );

MultiSFTCatalogView *XLALGetMultiSFTCatalogView ( const SFTCatalog *catalog );
void XLALDestroyMultiSFTCatalogView ( MultiSFTCatalogView *multiView );

/** @} */

/**
 * \name SFT catalog property functions
 * \anchor SFT-catalog-prop-func
 */
/** @{ */

int XLALCheckCRCSFTCatalog( BOOLEAN *crc_check, SFTCatalog *catalog );

LALStringVector *XLALListIFOsInCatalog( const SFTCatalog *catalog );
INT4 XLALCountIFOsInCatalog( const SFTCatalog *catalog );

const CHAR * XLALshowSFTLocator ( const struct tagSFTLocator *locator );

#ifndef SWIG // exclude from SWIG interface
int XLALSFTCatalogTimeslice( SFTCatalog *slice, const SFTCatalog *catalog, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS );
#endif
#ifdef SWIG // SWIG interface directives
SWIGLAL( RETURN_OWNED_BY_1ST_ARG( SFTCatalog*, XLALReturnSFTCatalogTimeslice) );
#endif
SFTCatalog *XLALReturnSFTCatalogTimeslice(const SFTCatalog *catalog, const LIGOTimeGPS *minStartGPS, const LIGOTimeGPS *maxStartGPS );

/** @} */

/**
 * \name SFT catalog generation functions
 * \anchor SFT-catalog-gen-func
 */
/** @{ */

SFTCatalog *XLALAddToFakeSFTCatalog( SFTCatalog *catalog, const CHAR *detector, const LIGOTimeGPSVector *timestamps );
SFTCatalog *XLALMultiAddToFakeSFTCatalog( SFTCatalog *catalog, const LALStringVector *detectors, const MultiLIGOTimeGPSVector *timestamps );

/** @} */

/**
 * \name SFT file naming convention functions
 * \anchor SFT-file-naming-func
 */
/** @{ */

int XLALFindCWDetector ( CHAR** prefix, INT4 *lalCachedIndex, const CHAR *name, const BOOLEAN exactMatch );
BOOLEAN XLALIsValidCWDetector ( const CHAR *name );
CHAR *XLALGetChannelPrefix ( const CHAR *name );
LALDetector *XLALGetSiteInfo ( const CHAR *name );

char *XLALOfficialSFTFilename ( char site, char channel, UINT4 numSFTs, UINT4 Tsft, UINT4 GPS_start, UINT4 Tspan, const char *Misc );
char *XLALGetOfficialName4SFT ( const SFTtype *sft, const char *Misc );
char *XLALGetOfficialName4MergedSFTs ( const SFTVector *sfts, const char *Misc );

int XLALCheckValidDescriptionField ( const char *desc );

/** @} */

/**
 * \name SFT file reading functions
 * \anchor SFT-file-read-func
 */
/** @{ */

LALStringVector *XLALFindFiles (const CHAR *globstring);

SFTVector* XLALLoadSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);

MultiSFTVector* XLALLoadMultiSFTs (const SFTCatalog *catalog, REAL8 fMin, REAL8 fMax);
MultiSFTVector *XLALLoadMultiSFTsFromView ( const MultiSFTCatalogView *multiCatalogView, REAL8 fMin, REAL8 fMax );

MultiSFTVector* XLALReadSFDB(REAL8 f_min, REAL8 f_max, const CHAR *file_pattern, const CHAR *timeStampsStarting, const CHAR *timeStampsFinishing);

/** @} */

/**
 * \name SFT file writing functions
 * \anchor SFT-file-write-func
 */
/** @{ */

int XLALWriteSFT2fp   ( const SFTtype *sft, FILE *fp, const CHAR *SFTcomment );
int XLALWriteSFT2file ( const SFTtype *sft, const CHAR *fname, const CHAR *SFTcomment );

int XLALWriteSFTVector2Dir  ( const SFTVector *sftVect, const CHAR *dirname, const CHAR *SFTcomment, const CHAR *Misc );
int XLALWriteSFTVector2File ( const SFTVector *sftVect, const CHAR *dirname, const CHAR *SFTcomment, const CHAR *Misc );
int XLALWriteSFTVector2NamedFile ( const SFTVector *sftVect, const CHAR *filename, const CHAR *SFTcomment );

int XLALValidateSFTFile ( const char *fname );

/** @} */

/**
 * \name PSD type create/destroy functions
 * \anchor PSD-type-cdtor-func
 */
/** @{ */

void XLALDestroyPSDVector ( PSDVector *vect );
void XLALDestroyMultiPSDVector ( MultiPSDVector *multvect );

void XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights );

/** @} */

/**
 * \name PSD type modify functions
 * \anchor PSD-type-mod-func
 */
/** @{ */

int XLALCropMultiPSDandSFTVectors ( MultiPSDVector *multiPSDVect, MultiSFTVector *multiSFTVect, UINT4 firstBin, UINT4 lastBin );

/** @} */

/**
 * \name PSD type property functions
 * \anchor PSD-type-prop-func
 */
/** @{ */

MultiNoiseWeights *XLALComputeMultiNoiseWeights ( const MultiPSDVector *rngmed, UINT4 blocksRngMed, UINT4 excludePercentile);

/** @} */

/**
 * \name PSD generation functions
 * \anchor PSD-type-gen-func
 */
/** @{ */

REAL8FrequencySeries *XLALComputeSegmentDataQ ( const MultiPSDVector *multiPSDVect, LALSeg segment );
REAL8 XLALMathOpOverArray(const REAL8* data, const size_t length, const MathOpType optype);
REAL8 XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs (  const UINT4 totalNumSFTs, const MathOpType optypeSFTs );

int XLALComputePSDandNormSFTPower ( REAL8Vector **finalPSD, MultiPSDVector **multiPSDVector, REAL8Vector **normSFT, MultiSFTVector *inputSFTs, const BOOLEAN returnMultiPSDVector, const BOOLEAN returnNormSFT,
                                    const UINT4 blocksRngMed, const MathOpType PSDmthopSFTs, const MathOpType PSDmthopIFOs, const MathOpType nSFTmthopSFTs, const MathOpType nSFTmthopIFOs,
                                    const BOOLEAN normalizeByTotalNumSFTs, const REAL8 FreqMin, const REAL8 FreqBand, const BOOLEAN normalizeSFTsInPlace );
int XLALComputePSDfromSFTs ( REAL8Vector **finalPSD, MultiSFTVector *inputSFTs, const UINT4 blocksRngMed, const MathOpType PSDmthopSFTs, const MathOpType PSDmthopIFOs,
                             const BOOLEAN normalizeByTotalNumSFTs, const REAL8 FreqMin, const REAL8 FreqBand );

/** @} */

/**
 * \name PSD file writing functions
 * \anchor PSD-file-write-func
 */
/** @{ */

int XLALDumpMultiPSDVector ( const CHAR *outbname, const MultiPSDVector *multiPSDVect );

int XLALWritePSDtoFilePointer ( FILE *fpOut, REAL8Vector *PSDVect, REAL8Vector *normSFTVect, BOOLEAN outputNormSFT, BOOLEAN outFreqBinEnd, INT4 PSDmthopBins, INT4 nSFTmthopBins,
                                INT4 binSize, INT4 binStep, REAL8 Freq0, REAL8 dFreq );

/** @} */

/** @} */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif

#endif     /* Close double-include protection */
