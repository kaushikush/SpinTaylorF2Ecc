/*
 *  Copyright (C) 2007 Badri Krishnan
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

 
/*
 *   Protection against double inclusion (include-loop protection)
 *     Note the naming convention!
 */

#ifndef _RADIOMETER_H
#define _RADIOMETER_H

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glob.h>
#include <time.h>
#include <errno.h> 

#include <lal/Date.h>
#include <lal/DetectorSite.h>
#include <lal/LALDatatypes.h>
#include <lal/LALHough.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/Velocity.h>
#include <lal/Statistics.h>
#include <lal/ComputeFstat.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LALInitBarycenter.h>
#include <lal/SFTClean.h>
#include <lalapps.h>
#include <gsl/gsl_cdf.h>


/******************************************************
 *   Protection against C++ name mangling
 */

#ifdef  __cplusplus
extern "C" {
#endif


/******************************************************
 *  Assignment of Id string using NRCSID()
 */

NRCSID (RADIOMETERH, "$Id$");

/******************************************************
 *  Error codes and messages.
 */
 
#define RADIOMETER_ENORM 0
#define RADIOMETER_ESUB  1
#define RADIOMETER_EARG  2
#define RADIOMETER_EBAD  3
#define RADIOMETER_EFILE 4
#define RADIOMETER_EDIR 4
#define RADIOMETER_ENULL 5
#define RADIOMETER_ENONULL 6

#define RADIOMETER_MSGENORM "Normal exit"
#define RADIOMETER_MSGESUB  "Subroutine failed"
#define RADIOMETER_MSGEARG  "Error parsing arguments"
#define RADIOMETER_MSGEBAD  "Bad argument values"
#define RADIOMETER_MSGEFILE "Could not create output file"
#define RADIOMETER_MSGEDIR  "Could not create directory"
#define RADIOMETER_MSGENULL "Null pointer"
#define RADIOMETER_MSGENONULL "Non-null pointer"

/* A global pointer for debugging. */
#ifndef NDEBUG
char *lalWatch;
#endif

#define PIXELFACTOR  2 


/* ******************************************************************
 *  Structure, enum, union, etc., typdefs.
 */

  /** struct holding info about skypoints */
  typedef struct tagSkyPatchesInfo{
    UINT4 numSkyPatches;
    REAL8 *alpha;
    REAL8 *delta;
    REAL8 *alphaSize;
    REAL8 *deltaSize;
  } SkyPatchesInfo;

  /** A pair of SFTs which will be correlated 
      -- this struct contains all information needed 
      to calculate the cross correlation for some given 
      parameter space point   */
  typedef struct tagSingleSFTpair{
    COMPLEX8FrequencySeries  *sft1; 
    COMPLEX8FrequencySeries  *sft2; 
    REAL8 vel1[3];
    REAL8 vel2[3];
    REAL8 pos1[3];
    REAL8 pos2[3];
    REAL8 a1, a2, b1, b2;
  } SingleSFTpair;

  /** vector of SFT pairs */
  typedef struct tagSFTPairVec{
    UINT4 length;
    SingleSFTpair *data;
  } SFTPairVec;  

  typedef struct tagSFTPairparams{
    REAL8 lag;
  } SFTPairParams;

  typedef struct tagSFTDetectorInfo{
    COMPLEX8FrequencySeries *sft;
    REAL8 vDetector[3];
    REAL8 rDetector[3];
    REAL8 a;
    REAL8 b;
  } SFTDetectorInfo;

/*
 *  Functions Declarations (i.e., prototypes).
 */


void SetUpRadiometerSkyPatches(LALStatus *status, 
			       SkyPatchesInfo *out,  
			       CHAR *skyFileName, 
			       CHAR *skyRegion, 
			       REAL8 dAlpha, 
			       REAL8 dDelta);


void CreateSFTPairs(LALStatus                *status,
		    SFTPairVec               *out,
		    MultiSFTVector           *inputSFTs,
		    MultiDetectorStateSeries *mdetStates,
		    MultiAMCoeffs            *multiAMcoef,
		    SFTPairParams            *par);


void CreateSFTPairsFrom2SFTvectors(LALStatus                 *status,
				   SFTPairVec                *out,
				   const SFTVector           *in1,
				   const SFTVector           *in2,
				   const DetectorStateSeries *det1,
				   const DetectorStateSeries *det2,
				   const AMCoeffs            *amc1,
				   const AMCoeffs            *amc2,
				   SFTPairParams             *par);


void FillSFTPair(LALStatus                 *status,
		 SingleSFTpair             *out,
		 COMPLEX8FrequencySeries   *sft1, 
		 COMPLEX8FrequencySeries   *sft2, 
		 DetectorState             *det1,
		 DetectorState             *det2,
		 REAL8                     a1,
		 REAL8                     a2,
		 REAL8                     b1,
		 REAL8                     b2);


/* ****************************************************** */

#ifdef  __cplusplus
}                /* Close C++ protection */
#endif


#endif     /* Close double-include protection _RADIOMETER_H */
