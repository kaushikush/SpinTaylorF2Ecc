/*
 * Copyright (C) 2010 Reinhard Prix (xlalified)
 * Copyright (C) 2004, 2005 Reinhard Prix
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

#ifndef _USERINPUT_H  /* Double-include protection. */
#define _USERINPUT_H

#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

#include <lal/ConfigFile.h>
#include <lal/ParseStringValue.h>
#if 0
#include <lal/LIGOMetadataTables.h>
#endif

/**
 * \defgroup UserInput_h Header UserInput.h
 * \ingroup lal_support
 * \author Reinhard Prix
 * \brief Module for simple unified handling of user-input from config-file and/or command-line.
 *
 *
 * ### Description ###
 *
 * This module provides functions and macros to 'register' a set off C-variables as 'User Variables',
 * which can then be read in from the commandline and/or an input config file, as parsed by \ref ConfigFile_h.
 *
 * The module also handles generating and outputting a help-string on the available input options when requested, and
 * can deal with enforcing input of required options and using default values.
 *
 * ### Usage ###
 *
 * The general approach consists of the following steps:
 * <ol>
 * <li> set default-values for optional user-variables as appropriate</li>
 * <li> \c register all user-variables using calls to \c XLALRegister<TYPE>UserVar(), or more conveniently, using the shortcut-macros
 * XLALreg<TYPE>UserStruct() that assume a pointer named 'uvar' to a struct containing all user-variables in the form 'uvar->UserVariable'.</li>
 * <li> parse user-input using XLALUserVarReadAllInput()</li>
 * </ol>
 *
 * One can use XLALUserVarWasSet() to determine whether a given user-input option has been set by the user.
 *
 * The function XLALUserVarGetLog() can be used to obtain a log-string containing the full user-input, either in \c commandline- or \c ConfigFile format.
 *
 * Here is a worked simple example of its recommended use:
 * \code
 * #include <stdio.h>
 * #include <lal/XLALError.h>
 * #include <lal/LALDatatypes.h>
 *
 * #include <lal/UserInput.h>
 *
 * // these are the C-variables we want to read in from user-input
 * typedef struct {
 *   BOOLEAN help;                // did user request help-output?
 *   INT4 anInteger;
 *   REAL8 aDoubleVar;
 *   CHAR *andAString;
 *   REAL8 specialGeekSwitch;
 *   LIGOTimeGPS someEpoch;
 *   REAL8 RA;
 *   REAL8 DEC;
 * } UserInput_t;
 *
 *
 * int main(int argc,char *argv[])
 * {
 *   UserInput_t XLAL_INIT_DECL(UserVariables); // initializes this struct to {0}
 *   UserInput_t *uvar = &UserVariables;          // struct-pointer allows us to use the XLALreg<TYPE>UserStruct() macros...
 *
 *   // 1. step: set default-values for optional user-input variables
 *   uvar->anInteger = 0;
 *   uvar->andAString = NULL;     // Note: need to assign allocated strings here as default!!
 *
 *   // 2. step: Register all user-variables using the shortcut macros:
 *   XLALregBOOLUserStruct  ( help,               'h',  UVAR_HELP,     "Output this help-message");
 *   XLALregINTUserStruct   ( anInteger,          'i',  UVAR_OPTIONAL, "An example user-variable of an optional integer");
 *   XLALregREALUserStruct  ( aDoubleVar,         'r',  UVAR_REQUIRED, "This REAL8 user-variable is required");
 *   XLALregSTRINGUserStruct( andAString,          0,   UVAR_OPTIONAL, "Optional string-input, has no short-option");
 *   XLALregEPOCHUserStruct ( someEpoch,           0,   UVAR_OPTIONAL, "Reference epoch (format 'xx.yy[GPS|MJD]')");
 *   XLALregLONGITUDEUserStruct(RA,                0,   UVAR_OPTIONAL, "Sky location: right ascension in [0,2pi] (in radians or hours:minutes:seconds)");
 *   XLALregLATITUDEUserStruct(DEC,                0,   UVAR_OPTIONAL, "Sky location: declination [-pi/2,pi/2] (in radians or degrees:minutes:seconds)");
 *
 *   XLALregREALUserStruct  ( specialGeekSwitch,   'g', UVAR_DEVELOPER, "This REAL8 user-variable may not be relevant for standard usage");
 *
 *   // 3. step: parse all user-input, from either config-file if given, or commandline (overloads config-file values)
 *   XLAL_CHECK ( XLALUserVarReadAllInput ( argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC);
 *
 *   if (uvar->help){      // if user had requested help, then we're already done here
 *     return 0;
 *   }
 *
 *   printf ("User-input was: anInteger = %d, aDoubleVar = %f, andAString = %s\n", uvar->anInteger, uvar->aDoubleVar, uvar->andAString );
 *   printf ("someEpoch = {%d s, %d ns}, RA = %f rad, DEC = %f rad\n", uvar->someEpoch.gpsSeconds, uvar->someEpoch.gpsNanoSeconds, uvar->RA, uvar->DEC );
 *
 *   // 4. step: free user-input module memory
 *   XLALDestroyUserVars();
 *
 *   LALCheckMemoryLeaks();
 *   return 0;
 * } // main()
 * \endcode
 *
 * \note This code can be compiled <b>as is</b> within lalapps, and yields
 *
 * \verbatim
$ LAL_DEBUG_LEVEL=1 ./tmp --help

Usage: ./tmp [@ConfigFile] [options], where options are:

-h, --help                BOOLEAN    Output this help-message []
-i, --anInteger           INT4       An example user-variable of an optional integer [0]
-r, --aDoubleVar          REAL8      This REAL8 user-variable is required [REQUIRED]
    --andAString          STRING     Optional string-input, has no short-option [NULL]
    --someEpoch           EPOCH      Reference epoch (format 'xx.yy[GPS|MJD]') [0.000000000GPS]
    --RA                  LONGITUDE  Sky location: right ascension in [0,2pi] (in radians or hours:minutes:seconds) [0.0]
    --DEC                 LATITUDE   Sky location: declination [-pi/2,pi/2] (in radians or degrees:minutes:seconds) [0.0]

---------- The following are 'Developer'-options not useful for most users:----------

-g, --specialGeekSwitch   REAL8      This REAL8 user-variable may not be relevant for standard usage [0.0]

\endverbatim
 *
 * And if called
 * \verbatim
$ ./tmp -r 3.1415 --andAString="stupid example" --someEpoch=5555MJD --RA=10:25:10.123 --DEC=-30:0:0
User-input was: anInteger = 0, aDoubleVar = 3.141500, andAString = stupid example
someEpoch = {2147483596 s, 816000000 ns}, RA = 2.727813 rad, DEC = -0.523599 rad
\endverbatim
 *
 * \note For a real-world example of usage, see various codes under lalapps/src/pulsar/{Injections,Fstatistic}
 *
 */
/*@{*/

/**
 * \name Shortcut Macros
 * With this family of short-cut macros one can conveniently register User-variables
 * that are accessible via a \e struct-pointer &quot;uvar->&quot;
 */
/*@{*/
#define XLALregREALUserStruct(name,option,flag,help) \
  XLALRegisterREALUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregINTUserStruct(name,option,flag,help) \
  XLALRegisterINTUserVar(#name, option,flag, help, &(uvar-> name))

#define XLALregBOOLUserStruct(name,option,flag,help) \
  XLALRegisterBOOLUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregSTRINGUserStruct(name,option,flag,help) \
  XLALRegisterSTRINGUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregLISTUserStruct(name,option,flag,help)                    \
  XLALRegisterLISTUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregEPOCHUserStruct(name,option,flag,help)                    \
  XLALRegisterEPOCHUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregLONGITUDEUserStruct(name,option,flag,help)                    \
  XLALRegisterLONGITUDEUserVar(#name, option, flag, help, &(uvar-> name))

#define XLALregLATITUDEUserStruct(name,option,flag,help)               \
  XLALRegisterLATITUDEUserVar(#name, option, flag, help, &(uvar-> name))

/*@}*/

/** State-flags: variable is optional, required, help, developer or was_set */
typedef enum {
  UVAR_OPTIONAL		= 0,	/**< not required, and hasn't been set */
  UVAR_REQUIRED 	= 1<<0,	/**< we require the user to set this variable */
  UVAR_HELP		= 1<<1,	/**< special variable: trigger output of help-string */
  UVAR_DEVELOPER	= 1<<2,	/**< OPTIONAL and hidden in help-output at lalDebugLevel==0 */
  UVAR_SPECIAL		= 1<<3,	/**< OPTIONAL and *turns off* checking of required variables (LALUserVarCheckRequired) */
  UVAR_WAS_SET 		= 1<<7	/**< flag that this user-var has been set by user */
} UserVarFlag;

/**
 * Format for logging User-input: configFile- or cmdLine-style.
 * This determines the format of the string returned from XLALUserVarGetLog().
 */
typedef enum {
  UVAR_LOGFMT_CFGFILE,		/**< return UserVars as a config-file */
  UVAR_LOGFMT_CMDLINE,		/**< return UserVars as a command-line */
  UVAR_LOGFMT_PROCPARAMS, 	/**< return UserVars suitable for filling in process-params struct */
  UVAR_LOGFMT_LAST
} UserVarLogFormat;

/* Function prototypes */
void XLALDestroyUserVars( void );
int XLALUserVarReadCmdline (int argc, char *argv[]);
int XLALUserVarReadCfgfile ( const CHAR *cfgfile );
CHAR *XLALUserVarHelpString ( const CHAR *progname );
int XLALUserVarReadAllInput ( int argc, char *argv[] );
int XLALUserVarCheckRequired( void );
int XLALUserVarWasSet (const void *cvar);
CHAR * XLALUserVarGetLog ( UserVarLogFormat format );

/* type-specific wrappers to XLALRegisterUserVar() to allow type-checking! */
int XLALRegisterREALUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar );
int XLALRegisterINTUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, INT4 *cvar );
int XLALRegisterBOOLUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, BOOLEAN *cvar );
int XLALRegisterSTRINGUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, CHAR **cvar );
int XLALRegisterLISTUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, LALStringVector **cvar);
int XLALRegisterEPOCHUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, LIGOTimeGPS *cvar);
int XLALRegisterLONGITUDEUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar);
int XLALRegisterLATITUDEUserVar ( const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar);


/* ========== Deprecated LAL interface wrappers ========== */

/**
 * \name Deprecated LAL-interface
 * These functions and macros are deprecated, and you should user their XLAL-equivalents instead.
 */
/*@{*/

void LALRegisterREALUserVar(LALStatus *, const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, REAL8 *cvar);
void LALRegisterINTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, INT4 *cvar);
void LALRegisterBOOLUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, BOOLEAN *cvar);
void LALRegisterSTRINGUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, CHAR **cvar);
void LALRegisterLISTUserVar (LALStatus *, const CHAR *name, CHAR optchar, UserVarFlag flag, const CHAR *helpstr, LALStringVector **cvar);

void LALDestroyUserVars (LALStatus *);

void LALUserVarReadAllInput(LALStatus *, int argc, char *argv[]);
void LALUserVarReadCmdline (LALStatus *, int argc, char *argv[]);
void LALUserVarReadCfgfile (LALStatus *, const CHAR *cfgfile);

void LALUserVarHelpString (LALStatus *, CHAR **helpstring, const CHAR *progname);
void LALUserVarCheckRequired (LALStatus *);
INT4 LALUserVarWasSet (const void *cvar);
void LALUserVarGetLog (LALStatus *, CHAR **logstr,  UserVarLogFormat format);

/*@}*/

/*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif  /* Double-include protection. */
