 /*
 * Copyright (C) 2020 Pep Covas, David Keitel, Rodrigo Tenorio
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

 /**
 * \file
 * \ingroup lalpulsar_WriteSFTsfromSFDBs
 * \author Pep Covas, David Keitel
 * \brief Read data in the form of Virgo Collaboration's SFDBs, convert to SFTs in the usual LALSuite format, and write those to disk.
 */

 #include <sys/stat.h>
 #include <sys/types.h>
 #include <lal/LALPulsarVCSInfo.h>
 #include <lal/UserInput.h>
 #include <lal/PulsarDataTypes.h>
 #include <lal/SFTfileIO.h>
 #include <lal/LogPrintf.h>

typedef struct
{
    REAL8 fmin;
    REAL8 fmax;
    CHAR *file_pattern;
    CHAR *timeStampsStarting;
    CHAR *timeStampsFinishing;
    CHAR *outSFTnames;
    CHAR *outSFTdir;
    BOOLEAN outSingleSFT;
} UserInput_t;

int initUserVars ( UserInput_t *uvar );

BOOLEAN is_directory ( const CHAR *fname );  // Function taken from makefakedata_v5

int main(int argc, char *argv[]) {

    UserInput_t XLAL_INIT_DECL(uvar);

    /* register all user-variables */
    XLAL_CHECK_MAIN ( initUserVars ( &uvar ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* read cmdline & cfgfile  */
    BOOLEAN should_exit = 0;
    XLAL_CHECK_MAIN( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( should_exit ) {
      exit (1);
    }

    XLAL_CHECK_MAIN( ( ( XLALUserVarWasSet ( &uvar.outSFTnames ) == 1 && XLALUserVarWasSet ( &uvar.outSFTdir ) == 0 ) || ( XLALUserVarWasSet(&uvar.outSFTnames ) == 0 && XLALUserVarWasSet ( &uvar.outSFTdir ) == 1 ) ), XLAL_EINVAL, "Only one of the two inputs (outSFTnames or outSFTdir) has to be used\n");

    MultiSFTVector *inputSFTs = NULL;
    XLAL_CHECK_MAIN ( (inputSFTs = XLALReadSFDB(uvar.fmin, uvar.fmax, uvar.file_pattern, uvar.timeStampsStarting, uvar.timeStampsFinishing)) != NULL, XLAL_EFUNC );

    LALStringVector *fnames = NULL;
    if (uvar.outSFTnames) {  // Parsing the input filenames into a string vector
        XLAL_CHECK_MAIN ( (fnames = XLALFindFiles (uvar.outSFTnames)) != NULL, XLAL_EFUNC, "Failed to parse file(s) with pattern '%s'.\n\n", uvar.outSFTnames );
    }
    else {
        XLAL_CHECK_MAIN ( is_directory ( uvar.outSFTdir ), XLAL_EFUNC, "Directory '%s' is not valid\n", uvar.outSFTdir );
    }

    /* generate comment string */
    size_t len;
    char *VCSInfoString = XLALVCSInfoString(lalPulsarVCSInfoList, 0, "%% ");
    XLAL_CHECK ( VCSInfoString != NULL, XLAL_EFUNC, "XLALVCSInfoString failed.\n" );
    char *logstr;
    XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
    char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(VCSInfoString) + 512 );
    XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zu) failed.\n", len );
    sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, VCSInfoString );

    SFTFilenameSpec XLAL_INIT_DECL(spec);
    XLAL_CHECK_MAIN( XLALFillSFTFilenameSpecStrings( &spec, uvar.outSFTdir, NULL, NULL, "unknown", "fromSFDBs", NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );

    for ( UINT4 X = 0; X < inputSFTs->length; X++) {

        SFTVector *sfts = inputSFTs->data[X];
        
        /* either write whole SFT-vector to single concatenated file or as individual SFT-files */
        if (uvar.outSFTnames) {
          XLAL_CHECK_MAIN ( XLALWriteSFTVector2NamedFile(sfts, fnames->data[X], spec.window_type, spec.window_param, comment) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        else {
          XLAL_CHECK_MAIN ( XLALWriteSFTVector2StandardFile( sfts, &spec, comment, uvar.outSingleSFT ) == XLAL_SUCCESS, XLAL_EFUNC );
        }

    }

    XLALDestroyMultiSFTVector (inputSFTs); 
    if (uvar.outSFTnames) XLALDestroyStringVector(fnames);
    XLALFree ( logstr );
    XLALFree ( comment );
    XLALFree ( VCSInfoString );

    return 0;

}


int initUserVars ( UserInput_t *uvar ) {
   
   XLAL_CHECK ( uvar != NULL, XLAL_EINVAL );

   uvar->fmin = 0;
   uvar->fmax = 0;
   uvar->file_pattern = NULL;
   uvar->timeStampsStarting = NULL;
   uvar->timeStampsFinishing = NULL;
   uvar->outSFTnames = NULL;
   uvar->outSFTdir = NULL;
   uvar->outSingleSFT = 1;

   /* now register all our user-variable */
   XLALRegisterUvarMember( file_pattern,            STRING, 'i', REQUIRED, "String of SFDB files (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember( timeStampsStarting,      STRING, 's', OPTIONAL, "File(s) containing the starting timestamps of science segments (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember( timeStampsFinishing,     STRING, 'f', OPTIONAL, "File(s) containing the starting timestamps of science segments (possibly from more than one detector, separated by a ;)");
   XLALRegisterUvarMember(   fmin,               REAL8, 0, REQUIRED, "Lowest frequency to extract from SFTs");
   XLALRegisterUvarMember(   fmax,               REAL8, 0, REQUIRED, "Highest frequency to extract from SFTs");
   XLALRegisterUvarMember( outSFTnames,         STRING, 'd', OPTIONAL, "Single SFT output custom path(s) (possibly from more than one detector, separated by a ;). Ignored if outSingleSFT is FALSE. Only use either outSFTnames or outSFTdir");
   XLALRegisterUvarMember( outSFTdir,          STRING, 'n', OPTIONAL,  "Directory for officially-named output SFTs. Only use either outSFTnames or outSFTdir.");
   XLALRegisterUvarMember(   outSingleSFT,       BOOLEAN, 0, OPTIONAL, "Write a single concatenated SFT file instead of individual files (default: TRUE)" );


   return XLAL_SUCCESS;
 
} /* initUserVars() */

/* determine if the given filepath is an existing directory or not */
BOOLEAN
is_directory ( const CHAR *fname )
{
  struct stat stat_buf;

  if ( stat (fname, &stat_buf) )
    return 0;

  if ( ! S_ISDIR (stat_buf.st_mode) )
    return 0;
  else
    return 1;

} /* is_directory() */