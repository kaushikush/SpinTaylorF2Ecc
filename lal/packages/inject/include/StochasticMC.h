
/*
<lalVerbatim file="StochasticMCHV">
Author: Tania Regimbau, Sukanta Bose, Jeff Noel
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{StochasticMC.h}}
\label{s:StochasticMC.h}
Routine used by the stochastic DSO to do software injection.
\subsection*{Synopsis}
\begin{verbatim}
#include <lal/StochasticMC.h>
\end{verbatim}
\noindent 

\subsection*{Error conditions}
\input{StochasticMCHE}
\subsection*{Structures}
</lalLaTeX>
*/
#ifndef _STOCHASTICMC_H
#define _STOCHASTICMC_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/Calibration.h>
#include <lal/SimulateSB.h>

#ifdef  __cplusplus
extern "C" {
#endif

  NRCSID( STOCHASTICMCH, "$Id$" );
  
/***************** <lalErrTable file="StochasticMCHE"> */

#define STOCHASTICMCH_ENULLP        1
#define STOCHASTICMCH_ENULLLEN      2
#define STOCHASTICMCH_ENULLSEG      3
#define STOCHASTICMCH_ENULLSRATE    4
#define STOCHASTICMCH_ENEGFMIN      5
#define STOCHASTICMCH_EMMEPOCH      6
#define STOCHASTICMCH_EMMUNITS      7


#define STOCHASTICMCH_MSGENULLP       "Null pointer"
#define STOCHASTICMCH_MSGENULLLEN     "Negative or zero length for data member of time series" 
#define STOCHASTICMCH_MSGENULLSEG      "Negative or zero number of segment" 
#define STOCHASTICMCH_MSGENULLSRATE    "Negative or zero sampling rate" 
#define STOCHASTICMCH_MSGENEGFMIN      "negative start frequency"
#define STOCHASTICMCH_MSGEMMEPOCH      "Mismatch in epochs"
#define STOCHASTICMCH_MSGEMMUNITS      "Mismatch in units"
 
/************************************ </lalErrTable> */

  typedef struct tagStochasticMCInput {
    CalibrationFunctions         calfuncs1;
    CalibrationFunctions         calfuncs2 ;
    CalibrationUpdateParams      calfacts1;
    CalibrationUpdateParams      calfacts2 ;
  } StochasticMCInput;

   
  typedef struct tagStochasticMCParams {
    UINT4 lengthseg;
    UINT4 numseg;
    REAL8 sRate;
    INT4  starttime;
    INT4  seed;
    REAL8 fRef;
    REAL8 f0;
    REAL8 omegaRef;
    REAL8 alpha;
    INT4 site1;
    INT4 site2; 
    
  } StochasticMCParams;



 void
 LALStochasticMCDso( LALStatus                  *status,
		  SSSimStochBGOutput         *output, 
		  StochasticMCInput          *input,
		  StochasticMCParams         *params );

 void
 LALStochasticMCDsoSplice( LALStatus                  *status,
		  SSSimStochBGOutput         *output, 
		  StochasticMCInput          *input,

#ifdef  __cplusplus
}
#endif

#endif /* _STOCHASTICMC_H */
  
