/*--------------------------------------------------------------------
 * $Id: lex_starter.l 3054 2014-07-21 09:43:44Z Claudia.Emde $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/
/* the "incl" state is used for picking up the name
 * of an include file
 */
%x include

%{

#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sun.h>
#include <string.h>
#include <float.h>
       
#include "uvspec.h"
#include "ascii.h"
#include "atmosphere.h"
#include "ancillary.h"
#include "netCDF_functions.h"
#include "solver.h"
#include "cdisort.h"
#include "uvspecrandom.h"

#define Q(x) #x
#define QUOTE(x) Q(x)

#define THEVERSION QUOTE(LRTVERSION)

/* #define COPYPROTECTION 1 */

/* following line needed to suppress compiler warning about unused function "input" */
#define YY_NO_INPUT

#define QSORT_CAST (int (*)(const void *, const void *))

#if HAVE_OPENGL
  extern int GLmystic_hiddenline; /* hidden line removal */
#else
  int GLmystic_hiddenline; /* dummy if OPENGL not available */
#endif

static float *add_cmu_to_umu(int nstr, int numu, int maxumu, float* oldumu, int* cmuind, int* umuind);

static float  yytext2float(char* yytext, int which_token, int max_token);
static float* yytext2floats(char* yytext, int nary);
static double yytext2double(char* yytext, int which_token, int max_token);
static int    yytext2int(char* yytext, int which_token, int max_token);
static int*   yytext2integers(char* yytext, int nary);
static char*  yytext2string(char* yytext, int which_token, int max_token);
static int    yytext2ntokens(char* yytext);
static int    yytextCstring(char *dest, char* yytext, int which_token, int max_token);
static void free_output(input_struct input, output_struct *output);
#ifndef UVSPEC_FUNCTION
static int print_output(input_struct input, output_struct *output);
#else
static int print_output(input_struct input, output_struct *output, char *outfilename);
#endif
void write_flexstor_output     (input_struct Input, output_struct *output, FILE *f);
void write_polradtran_output   (input_struct Input, output_struct *output, FILE *f);
int write_header              (input_struct Input, output_struct *output, FILE *f);
void write_spectral_heat       (input_struct Input, output_struct *output, FILE *f);
int write_user_defined_output (input_struct Input, output_struct *output, FILE *f);
int check_if_result_exists     (int *output_user, int n_output_user, int solver, int **valid_result, int quiet);
int write_netCDF_output        (input_struct Input, output_struct *output);

char *s, *t, *u;
int ntokens;
int mol_id;
int id=0, id1=0, id2=0;
float tmpfloat;
int   tmpint;
float* tmp;
int*   tmpi;
int ierror=0;
int line_number;
input_struct Input;
output_struct Output;
int i=0, j=0, isp=0, iscr=0, ispoff=0;
int iscloudcover=0, isipafiles=0;
int position=0;
float epsilon=1E-7;
  
#define EPSILON 1.E-6        /*relative rounding error for floats */
#define MAX_INCLUDE_DEPTH 10
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0;

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* This is a fix needed at MIM, need to test whether this also works on DLR */
#if COMPILING_CONDOR_AT_MIM
#define fscanf __fscanf
#endif

#define UVSPEC_TOKEN " |\t "

%}

%option noyywrap
%option nounput

DIGIT    [0-9]
SPACE    [ \t]
CCHAR    [#]
COMMENT  ({CCHAR}.*)
INTEGER  {DIGIT}+
SPACE_INTEGER {SPACE}+{INTEGER}
INTEGERS {SPACE_INTEGER}+ 
EXPON   ([DdEe][+-]?{DIGIT}+)
FLOAT    ((({DIGIT}*+"."{DIGIT}*)|({DIGIT}*+"."{DIGIT}*{EXPON}))|({INTEGER}))
SPACE_FLOAT {SPACE}+{FLOAT}
FLOATS   {SPACE_FLOAT}+ 
SIGN     [-+]
SIGNED_FLOAT {SIGN}?{FLOAT}
SPACE_SIGNED_FLOAT {SPACE}+{SIGNED_FLOAT}
SIGNED_FLOATS {SPACE_SIGNED_FLOAT}+
SIGNED_INTEGER {SIGN}?{INTEGER}
WORD     [_\.\/A-Za-z][_\.\/A-Za-z0-9\-]*
DIM      [1,3][d,D]
SPACE_PROFILETYPE {SPACE}+({DIM}|{WORD})
SPACE_WORD {SPACE}+{WORD}
WORDS {SPACE_WORD}*
FLOAT_OR_WORD ({SPACE_FLOAT}|{SPACE_WORD})
FLOATS_OR_WORDS {FLOAT_OR_WORD}*
SIGNED_FLOAT_OR_WORD ({SPACE_SIGNED_FLOAT}|{SPACE_WORD})
SIGNED_FLOATS_OR_WORDS {SIGNED_FLOAT_OR_WORD}*

%%

  /* These are the allowed input lines in the uvspec input file 
     The input line specifications must come before the eating of
     blanks, otherwise things get messed up 
     Obviously this whole lot could have been handled more properly
     using bison etc., but then rather few would know what was going 
     on??? 
   */


 /* PYTHON_OPTIONS */
aerosol_species_file{SPACE_WORD}{SPACE_WORD}+ {
				     ntokens = yytext2ntokens(yytext);
                                     yytextCstring (Input.aer.filename[FN_AER_SPECIES], yytext, 1+1, ntokens+1);
				     Input.aer.n_species = ntokens - 1; /* one token for type */
                                     Input.aer.species_names = (char **) calloc (Input.aer.n_species, sizeof(char *));
                                     for (i=3;i<=ntokens+1; i++) {
				         Input.aer.species_names[i-3] = yytext2string(yytext, i, ntokens+1);
				     }
                                     Input.aer.spec =1; } 

ECHAM_atmosphere_file{SPACE}+{WORD} {Input.atm.rs_source=RS_FROM_ECHAM;
                                      yytextCstring(Input.filename[FN_ECHAM],yytext,2,2);
                                      Input.atm.n_rs_gas=0;
                                    Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_H2O; Input.atm.n_rs_gas++; mol_id = MOL_H2O;
                                      Input.atm.rs_unit[mol_id] = MMR;
                                     }

ECHAM_clouds_file{SPACE}+{WORD}     {isp=get_caoth_index(&Input.caoth,&Input.n_caoth,"wc",1);
                                     Input.caoth[isp].source = CAOTH_FROM_ECHAM;
                                     Input.caoth[isp].layer=TRUE;
                                     yytextCstring (Input.caoth[isp].filename,yytext,2,2);
                                     isp=get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
                                     Input.caoth[isp].source = CAOTH_FROM_ECHAM;
                                     Input.caoth[isp].layer=TRUE;
                                     yytextCstring (Input.caoth[isp].filename,yytext,2,2);
				    }
ECMWF_atmosphere_file{SPACE}+{WORD} {Input.atm.rs_source=RS_FROM_ECMWF;
                                     yytextCstring(Input.filename[FN_ECMWF],yytext,2,2);
                                     Input.atm.n_rs_gas=0;
                                     Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_O3;
				     Input.atm.n_rs_gas++;
				     mol_id = MOL_O3;
                                     Input.atm.rs_unit[mol_id] = MMR;
                                     Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_H2O;
				     Input.atm.n_rs_gas++;
				     mol_id = MOL_H2O;
                                     Input.atm.rs_unit[mol_id] = MMR;
                                    }
ECMWF_ic_reff{SPACE_WORD}{FLOAT_OR_WORD}* {ntokens = yytext2ntokens(yytext);
                                     if (ntokens < 1 || ntokens > 2) {
                                       fprintf (stderr, "Error, wrong number of arguments for option 'ECMWF_ic_reff' on line %d \n",line_number+1); 
                                       ierror++;
                                     }
                                     isp=get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
                                     t = yytext2string(yytext, 2,ntokens+1);
				     if      (strncasecmp("fixed",t,5)==0)  {Input.caoth[isp].reff_prop = REFF_FIXED;
                                       if (ntokens != 2) {
                                         fprintf (stderr, "Error, wrong number of arguments for option 'ECMWF_ic_reff fixed' on line %d \n",line_number+1); 
                                         ierror++;
                                       }
                                       else Input.caoth[isp].reff =  yytext2float(yytext,2+1,ntokens+1);    
                                     }
				     else if (strncasecmp("file",t,5)==0)  {Input.caoth[isp].reff_prop = REFF_FILE;
                                       if (ntokens != 2) {
                                         fprintf (stderr, "Error, wrong number of arguments for option 'ECMWF_ic_reff file' on line %d \n",line_number+1); 
                                         ierror++;
                                       }
                                       yytextCstring(Input.filename[FN_IC_REFF],yytext,2+1,ntokens+1);    
                                     }
				     else if (strncasecmp("Ou"   ,t,2)==0)  {Input.caoth[isp].reff_prop = REFF_OU;}
				     else {fprintf(stderr,"Unknown 'ECMWF_ic_reff %s' on line %d\n", s, line_number+1);
				           ierror++;}
                                     free(t);
                                    }
include                             BEGIN(include);
mc_circ{SPACE_FLOAT}*{SPACE_INTEGER}* {
                                     s = (char *) calloc (strlen(yytext)+1,sizeof(char));
                                     strcpy (s, yytext);
                                     ntokens = yytext2ntokens(s);
				     Input.rte.mc.dcirc=1000.;
				     Input.rte.mc.ncirc=10;
				     if (ntokens>0)
				       Input.rte.mc.dcirc=yytext2double(yytext,2,ntokens+1);
				     if (ntokens>1)
				       Input.rte.mc.ncirc=yytext2int(yytext,3,ntokens+1);
				     free(s); }
mc_escape{SPACE_WORD}*             {
                                    s = (char *) calloc (strlen(yytext)+1,sizeof(char));
				    strcpy (s, yytext);
                                    ntokens = yytext2ntokens(s);
                                    free(s);            
                                    if (ntokens == 0) { Input.rte.mc.escape = 1;
                                                        Input.rte.mc.locest = 0;}
                                    if (ntokens == 1) {
				      s = yytext2string(yytext,2,ntokens+1);
				      if (strncasecmp(s,"on",2)==0) {
				        Input.rte.mc.escape = 1;
				        Input.rte.mc.locest = 0;
				      } else {
				        Input.rte.mc.escape = 0;
                                      }				      
                         	      free(s);
				    }
                                   }
mc_spectral_is_wvl{FLOATS}          {
                                     Input.rte.mc.spectral_is = TRUE;
                                     Input.rte.mc.spectral_is_nwvl = yytext2ntokens(yytext);
                                     Input.rte.mc.spectral_is_wvl = calloc (Input.rte.mc.spectral_is_nwvl,sizeof(float));
                                     for (i=2;i<=Input.rte.mc.spectral_is_nwvl+1; i++)
                                       Input.rte.mc.spectral_is_wvl[i-2] = yytext2float(yytext,i,Input.rte.mc.spectral_is_nwvl+1);
                                    }
output_user{WORDS}                  {s = (char *) calloc (strlen(yytext)+1,sizeof(char));
				     strcpy (s, yytext);
                                     Input.output_user_flag = 1; 
                                     Input.n_output_user = yytext2ntokens(s); 
				     free(s);
				     Input.output_user     = calloc (Input.n_output_user, sizeof (int));
				     Input.output_user_gas = calloc (Input.n_output_user, sizeof (int));
				     for (i=2;i<=Input.n_output_user+1; i++) {
				       t = yytext2string(yytext, i, Input.n_output_user+1);  
				       if      (strcmp("lambda",t)==0 ||  
                                                strcmp("wavelength",t)==0)    {Input.output_user[i-2] = OUTPUT_USER_WAVE;}
				       else if (strcmp("lambda_min",t)==0)    {Input.output_user[i-2] = OUTPUT_USER_WAVE_MIN;}
				       else if (strcmp("lambda_max",t)==0)    {Input.output_user[i-2] = OUTPUT_USER_WAVE_MAX;}
				       else if (strcmp("wavenumber",t)==0)    {Input.output_user[i-2] = OUTPUT_USER_WAVENUMBER;}
				       else if (strcmp("wavenumber_min",t)==0){Input.output_user[i-2] = OUTPUT_USER_WAVENUMBER_MIN;}
				       else if (strcmp("wavenumber_max",t)==0){Input.output_user[i-2] = OUTPUT_USER_WAVENUMBER_MAX;}
				       else if (strcmp("zout",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_ZOUT_SUR;}
				       else if (strcmp("zout_sur",t)==0)      {Input.output_user[i-2] = OUTPUT_USER_ZOUT_SUR;}
				       else if (strcmp("zout_sea",t)==0)      {Input.output_user[i-2] = OUTPUT_USER_ZOUT_SEA;}
				       else if (strcmp("z_sur",t)   ==0)      {Input.output_user[i-2] = OUTPUT_USER_Z_SUR;}
				       else if (strcmp("p",t)       ==0)      {Input.output_user[i-2] = OUTPUT_USER_P;}
				       else if (strcmp("T",t)       ==0)      {Input.output_user[i-2] = OUTPUT_USER_T;}
				       else if (strcmp("T_d",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_T_D;}
				       else if (strcmp("T_sur",t)   ==0)      {Input.output_user[i-2] = OUTPUT_USER_T_SUR;}
				       else if (strcmp("theta",t)   ==0)      {Input.output_user[i-2] = OUTPUT_USER_THETA;}
				       else if (strcmp("theta_e",t) ==0)      {Input.output_user[i-2] = OUTPUT_USER_THETA_E;}
				       else if (strcmp("edir",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_EDIR;}
				       else if (strcmp("eglo",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_EGLO;}
				       else if (strcmp("edn",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_EDN;}
				       else if (strcmp("eup",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_EUP;}
				       else if (strcmp("enet",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_ENET;}
				       else if (strcmp("esum",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_ESUM;}
				       else if (strcmp("fdir",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_FDIR;}
				       else if (strcmp("fglo",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_FGLO;}
				       else if (strcmp("fdn",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_FDN;}
				       else if (strcmp("fup",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_FUP;}
				       else if (strcmp("f",t)       ==0)      {Input.output_user[i-2] = OUTPUT_USER_F;}
				       else if (strcmp("uavgdir",t) ==0)      {Input.output_user[i-2] = OUTPUT_USER_UDIR;}
				       else if (strcmp("uavgglo",t) ==0)      {Input.output_user[i-2] = OUTPUT_USER_UGLO;}
				       else if (strcmp("uavgdn",t)  ==0)      {Input.output_user[i-2] = OUTPUT_USER_UDN;}
				       else if (strcmp("uavgup",t)  ==0)      {Input.output_user[i-2] = OUTPUT_USER_UUP;}
				       else if (strcmp("uavg",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_U;}
                                       else if (strcmp("uu",t)      ==0)      {Input.output_user[i-2] = OUTPUT_USER_UU;}
				       else if (strcmp("albedo",t)  ==0)      {Input.output_user[i-2] = OUTPUT_USER_ALB;}
				       else if (strcmp("sza",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_SZA;}
				       else if (strcmp("c_p",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_C_P;}
				       else if (strcmp("heat",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_HEAT;}
				       else if (strcmp("emis",t)    ==0)      {Input.output_user[i-2] = OUTPUT_USER_EMIS;}
				       else if (strcmp("abs",t)     ==0)      {Input.output_user[i-2] = OUTPUT_USER_ABS;}
				       else if (strcmp("w_rad",t)   ==0)      {Input.output_user[i-2] = OUTPUT_USER_W_RAD;}
				       else if (strcmp("m_rad",t)   ==0)      {Input.output_user[i-2] = OUTPUT_USER_M_RAD;}
                                       else if (strncasecmp("n_",  t,2)==0)   {Input.output_user[i-2] = OUTPUT_USER_N;}
                                       else if (strncasecmp("rho_",t,4)==0)   {Input.output_user[i-2] = OUTPUT_USER_RHO;}
                                       else if (strncasecmp("vmr_",t,4)==0)   {Input.output_user[i-2] = OUTPUT_USER_VMR;}
                                       else if (strncasecmp("mmr_",t,4)==0)   {Input.output_user[i-2] = OUTPUT_USER_MMR;}
                                       else if (strcmp("rh",t)        ==0)    {Input.output_user[i-2] = OUTPUT_USER_RH;}
                                       else if (strcmp("rh_ice",t)    ==0)    {Input.output_user[i-2] = OUTPUT_USER_RH_ICE;}
                                       else if (strcasecmp("clwc",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_CLWC;}
                                       else if (strcasecmp("clwd",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_CLWD;}
                                       else if (strcasecmp("tclw",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_TCLW;}
                                       else if (strcasecmp("reff_water",t)==0){Input.output_user[i-2] = OUTPUT_USER_REFF_WAT;}
                                       else if (strcasecmp("ciwc",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_CIWC;}
                                       else if (strcasecmp("ciwd",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_CIWD;}
                                       else if (strcasecmp("tciw",t)  ==0)    {Input.output_user[i-2] = OUTPUT_USER_TCIW;}
                                       else if (strcasecmp("reff_ice",t) ==0) {Input.output_user[i-2] = OUTPUT_USER_REFF_ICE;}
                                       else if (strcasecmp("tcw",t)  ==0)     {Input.output_user[i-2] = OUTPUT_USER_TCW;}
                                       else if (strcasecmp("cc",t)   ==0)     {Input.output_user[i-2] = OUTPUT_USER_CC;}
                                       else if (strcasecmp("tcc",t)   ==0)    {Input.output_user[i-2] = OUTPUT_USER_TCC;}
                                       else if (strcasecmp("clouds",t)==0)    {Input.output_user[i-2] = OUTPUT_USER_CLOUDS;}
                                       else if (strcmp("u",t)         ==0)    {Input.output_user[i-2] = OUTPUT_USER_WIND_U;}
                                       else if (strcmp("v",t)         ==0)    {Input.output_user[i-2] = OUTPUT_USER_WIND_V;}
                                       else if (strcmp("w",t)         ==0)    {Input.output_user[i-2] = OUTPUT_USER_WIND_W;}
                                       else if (strcmp("dtheta_dx",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_DTDX;}
                                       else if (strcmp("dtheta_dy",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_DTDY;}
                                       else if (strcmp("dtheta_dz",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_DTDZ;}
                                       else if (strcmp("heat_ad_x",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_HEAT_AD_X;}
                                       else if (strcmp("heat_ad_y",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_HEAT_AD_Y;}
                                       else if (strcmp("heat_ad_z",t) ==0)    {Input.output_user[i-2] = OUTPUT_USER_HEAT_AD_Z;}
                                       else if (strcmp("heat_ad",t)   ==0)    {Input.output_user[i-2] = OUTPUT_USER_HEAT_AD;}
                                       else if (strcmp("spher_alb",t)   ==0)  {Input.output_user[i-2] = OUTPUT_USER_SPHER_ALB;;}
				       else {fprintf(stderr,"Unknown 'output_user %s' on line %d\n", t, line_number+1);
				       ierror++;}
                                       if (Input.output_user[i-2] == OUTPUT_USER_N || 
                                           Input.output_user[i-2] == OUTPUT_USER_RHO || 
                                           Input.output_user[i-2] == OUTPUT_USER_VMR || 
                                           Input.output_user[i-2] == OUTPUT_USER_MMR) {
                                         if      ((s = strstr(t,"air")) != NULL || (s = strstr(t,"AIR")) != NULL)
                                            Input.output_user_gas[i-2] = MOL_AIR;
                                         else if ((s = strstr(t,"o3"))  != NULL || (s = strstr(t,"O3"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_O3;
                                         else if ((s = strstr(t,"h2o")) != NULL || (s = strstr(t,"H2O")) != NULL)
                                            Input.output_user_gas[i-2] = MOL_H2O;
                                         else if ((s = strstr(t,"co2")) != NULL || (s = strstr(t,"CO2")) != NULL)
                                            Input.output_user_gas[i-2] = MOL_CO2;
                                         else if ((s = strstr(t,"no2")) != NULL || (s = strstr(t,"NO2")) != NULL)
                                            Input.output_user_gas[i-2] = MOL_NO2;
                                         else if ((s = strstr(t,"bro")) != NULL || (s = strstr(t,"BRO")) != NULL)
                                            Input.output_user_gas[i-2] = MOL_BRO;
                                         else if ((s = strstr(t,"oclo"))!= NULL || (s = strstr(t,"OCLO"))!= NULL)
                                            Input.output_user_gas[i-2] = MOL_OCLO;
                                         else if ((s = strstr(t,"hcho"))!= NULL || (s = strstr(t,"HCHO"))!= NULL)
                                            Input.output_user_gas[i-2] = MOL_HCHO;
                                         else if ((s = strstr(t,"o4"))  != NULL || (s = strstr(t,"O4"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_O4;
                                         else if ((s = strstr(t,"so2"))  != NULL || (s = strstr(t,"SO2"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_SO2;
                                         else if ((s = strstr(t,"o2"))  != NULL || (s = strstr(t,"O2"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_O2;
                                         else if ((s = strstr(t,"ch4"))  != NULL || (s = strstr(t,"CH4"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_CH4;
                                         else if ((s = strstr(t,"n2o"))  != NULL || (s = strstr(t,"N2O"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_N2O;
                                         else if ((s = strstr(t,"co"))  != NULL || (s = strstr(t,"CO"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_CO;
                                         else if ((s = strstr(t,"n2"))  != NULL || (s = strstr(t,"N2"))  != NULL)
                                            Input.output_user_gas[i-2] = MOL_N2;
                                         else {
                                           fprintf(stderr,"Unknown 'output_user %s' on line %d\n", t, line_number+1);
				           ierror++;
                                         }
                                       }                                         
                                       free(t);
				     }
                                    }

pressure_out{SIGNED_FLOATS_OR_WORDS}|p_out{SIGNED_FLOATS_OR_WORDS} {
                                      Input.atm.nzout = yytext2ntokens(yytext);
                                      Input.atm.press_zout = calloc (Input.atm.nzout,sizeof(float));
                                      for (i=2;i<=Input.atm.nzout+1; i++) {
                                       Input.atm.zout_source = OUTLEVEL_PRESS;
				       t = yytext2string(yytext, i, Input.atm.nzout+1);
				       if      (strncasecmp("toa",t,3)==0) {Input.atm.press_zout[i-2] = ZOUT_TOA;}
				       else if (strncasecmp("sur",t,3)==0) {Input.atm.press_zout[i-2] = ZOUT_SURFACE;}
				       else if (strncasecmp("cpt",t,3)==0) {Input.atm.press_zout[i-2] = ZOUT_CPT;}
                                       else                                {Input.atm.press_zout[i-2] = yytext2float(yytext,i,Input.atm.nzout+1);}
                                       free(t);
				      }
				     }
radiosonde{SPACE}+{WORD}{WORDS}     {Input.atm.rs_source=RS_FROM_FILE;
                                     s = (char *) calloc (strlen(yytext)+1,sizeof(char));
                                     strcpy (s, yytext);
                                     ntokens = yytext2ntokens(s);
                                     free(s);
				     s =  yytext2string(yytext,2,ntokens+1);
				     strcpy (Input.filename[FN_RADIOSONDE], s);
                                     Input.atm.n_rs_gas=0;
                                     free(s);
                                     for (i=3;i<=ntokens+1; i++) {
				       t = yytext2string(yytext, i, ntokens+1);  
				       if      (strcmp("O3",t)  ==0||strcmp("o3",t)  ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_O3;  Input.atm.n_rs_gas++; mol_id = MOL_O3;  }
				       else if (strcmp("O2",t)  ==0||strcmp("o2",t)  ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_O2;  Input.atm.n_rs_gas++; mol_id = MOL_O2;  }
				       else if (strcmp("H2O",t) ==0||strcmp("h2o",t) ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_H2O; Input.atm.n_rs_gas++; mol_id = MOL_H2O; }
				       else if (strcmp("CO2",t) ==0||strcmp("co2",t) ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_CO2; Input.atm.n_rs_gas++; mol_id = MOL_CO2; }
				       else if (strcmp("NO2",t) ==0||strcmp("no2",t) ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_NO2; Input.atm.n_rs_gas++; mol_id = MOL_NO2; }
				       else if (strcmp("BRO",t) ==0||strcmp("bro",t) ==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_BRO; Input.atm.n_rs_gas++; mol_id = MOL_BRO; }
				       else if (strcmp("OCLO",t)==0||strcmp("oclo",t)==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_OCLO;Input.atm.n_rs_gas++; mol_id = MOL_OCLO;}
				       else if (strcmp("HCHO",t)==0||strcmp("hcho",t)==0) {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_HCHO;Input.atm.n_rs_gas++; mol_id = MOL_HCHO;}
				       else if (strcmp("O4",t)==0  ||strcmp("o4",t)==0)   {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_O4;  Input.atm.n_rs_gas++; mol_id = MOL_O4;}
				       else if (strcmp("SO2",t)==0 ||strcmp("so2",t)==0)  {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_SO2; Input.atm.n_rs_gas++; mol_id = MOL_SO2;}
				       else if (strcmp("CH4",t)==0 ||strcmp("ch4",t)==0)  {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_CH4; Input.atm.n_rs_gas++; mol_id = MOL_CH4;}
				       else if (strcmp("N2O",t)==0 ||strcmp("n2o",t)==0)  {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_N2O; Input.atm.n_rs_gas++; mol_id = MOL_N2O;}
				       else if (strcmp("CO",t)==0 ||strcmp("co",t)==0)    {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_CO;  Input.atm.n_rs_gas++; mol_id = MOL_CO;}
				       else if (strcmp("N2",t)==0 ||strcmp("n2",t)==0)    {Input.atm.rs_gas[Input.atm.n_rs_gas]=MOL_N2;  Input.atm.n_rs_gas++; mol_id = MOL_N2;}
				       else if (strcmp("CM_3",t)==0||strcmp("cm-3",t)==0) { if (mol_id==0 || mol_id>=MOL_NN) ierror++; 
                                                                                            else Input.atm.rs_unit[mol_id] = CM_3; }
				       else if (strcmp("M_3",t) ==0||strcmp("m-3",t) ==0) { if (mol_id==0 || mol_id>=MOL_NN) ierror++;
                                                                                            else Input.atm.rs_unit[mol_id] = M_3;  }
				       else if (strcmp("MMR",t) ==0||strcmp("mmr",t) ==0) { if (mol_id==0 || mol_id>=MOL_NN) ierror++;
                                                                                            else Input.atm.rs_unit[mol_id] = MMR;  }
				       else if (strcmp("VMR",t) ==0||strcmp("vmr",t) ==0) { if (mol_id==0 || mol_id>=MOL_NN) ierror++;
                                                                                            else Input.atm.rs_unit[mol_id] = VMR;  }
                                       else if (strcmp("RH" ,t) ==0||strcmp("rh",t) == 0) {
                                         if (mol_id==0 || mol_id>=MOL_NN) 
                                           ierror++;
                                         else if (mol_id != MOL_H2O) {
                                           fprintf(stderr,"Error, unit \"relative humidity\" only possible for H2O mol_files.\n");
                                           fprintf(stderr,"See 'mol_modify ... %s' in line %d\n", s, line_number+1);
                                           ierror++;
                                         }
                                         else 
                                           Input.atm.rs_unit[mol_id]= RH; 
                                        }
				       else {fprintf(stderr,"Unknown 'radiosonde ... %s' on line %d\n", t, line_number+1);
				             ierror++;}
                                       free(t);
				     }
                                     if (Input.atm.n_rs_gas>2)
                                     for (i=0;i<Input.atm.n_rs_gas-1;i++) {
                                       for (j=i+1;j<Input.atm.n_rs_gas;j++)
                                         if (Input.atm.rs_gas[i] == Input.atm.rs_gas[j]) {
                                           fprintf(stderr, "Error in input file, 'radiosonde', on line %d\n", line_number+1);
                                           fprintf(stderr, "At least one species is specified more than once!\n");
                                           ierror++;
                                         }
                                     }  
                                     if (ierror!=0) {
                                       fprintf(stderr,"Error in input file, 'radiosonde', on line %d\n",line_number+1);  
                                       fprintf(stderr,"Please check the arguments and their order (first species, than unit)\n"); 
                                     }
                                    }
time{WORDS}{FLOATS}                 {ntokens = yytext2ntokens(yytext);
                                     tmpint=0;
                                     s =  yytext2string(yytext,2,ntokens+1);
                                     if      (strncasecmp(s,"UTC",3) == 0) { Input.time_zone = TIME_UTC; tmpint=1; }
				     else if (strncasecmp(s,"LAT",3) == 0) { Input.time_zone = TIME_LAT; tmpint=1; }
                                     else {Input.time_zone = TIME_UTC; tmpint=0;}
                                     if (Input.time_zone == TIME_UTC) {
                                      if (ntokens == 6+tmpint) {
                                        Input.UTC.tm_year  = yytext2int  (yytext, 2+tmpint, ntokens+1) - 1900; /* !!! years after 1900             !!! */
                                        Input.UTC.tm_mon   = yytext2int  (yytext, 3+tmpint, ntokens+1) -    1; /* !!! month after January 0 ... 11 !!! */
                                        Input.UTC.tm_mday  = yytext2int  (yytext, 4+tmpint, ntokens+1);   
                                        Input.UTC.tm_hour  = yytext2int  (yytext, 5+tmpint, ntokens+1); 
                                        Input.UTC.tm_min   = yytext2int  (yytext, 6+tmpint, ntokens+1);
                                        Input.UTC.tm_sec   = yytext2float(yytext, 7+tmpint, ntokens+1);                                     
                                      }
                                      else {fprintf(stderr,"Error, wrong number of arguments %d for %s on line %d\n",ntokens, s, line_number+1);
                                            fprintf(stderr,"  Please specify the time with the format: YYYY MM DD HH MM SS.\n");
				            ierror++;
                                      }
                                     }
                                     else if (Input.time_zone == TIME_LAT) {
                                      if (ntokens == 6+tmpint) {
                                        Input.LAT.tm_year  = yytext2int  (yytext, 2+tmpint, ntokens+1) - 1900; /* !!! years after 1900             !!! */
                                        Input.LAT.tm_mon   = yytext2int  (yytext, 3+tmpint, ntokens+1) -    1; /* !!! month after January 0 ... 11 !!! */
                                        Input.LAT.tm_mday  = yytext2int  (yytext, 4+tmpint, ntokens+1);   
                                        Input.LAT.tm_hour  = yytext2int  (yytext, 5+tmpint, ntokens+1); 
                                        Input.LAT.tm_min   = yytext2int  (yytext, 6+tmpint, ntokens+1);
                                        Input.LAT.tm_sec   = yytext2float(yytext, 7+tmpint, ntokens+1);                                     
                                      }
                                      else {fprintf(stderr,"Error, wrong number of arguments %d for %s on line %d\n",ntokens, s, line_number+1);
                                            fprintf(stderr,"  Please specify the time with the format: YYYY MM DD HH MM SS.\n");
				            ierror++;
                                      }
                                     }
                                     else {
                                      fprintf(stderr,"Error, unknown time zone %d, please report this program bug!\n", Input.time_zone);
				      ierror++;
                                     }
                                     free(s);
                                    } 
z_interpolate{WORDS}                {//s = (char *) calloc (strlen(yytext)+1,sizeof(char));
                                     //strcpy (s, yytext);
                                     ntokens = yytext2ntokens(yytext); 
                                     //free(s);
                                     switch (ntokens) {
                                     case 0:
                                       Input.atm.z_interpolate = Z_INTERPOLATE_LINMIX;
                                       break;
                                     case 1:
                                       t = yytext2string(yytext, 2, ntokens+1);
                                       if      (strncasecmp("default",t, 7)==0)  { Input.atm.z_interpolate = Z_INTERPOLATE_LINMIX; }
                                       else if (strncasecmp("lin_mix",t, 7)==0 || strncasecmp("linmix",t, 6)==0)
                                                                                 { Input.atm.z_interpolate = Z_INTERPOLATE_LINMIX; }
                                       else if (strncasecmp("log_lin",t, 7)==0 || strncasecmp("loglin" ,t, 6)==0)  
                                                                                 { Input.atm.z_interpolate = Z_INTERPOLATE_LOGLIN; }
                                       else if (strncasecmp("spline" ,t, 6)==0)  { Input.atm.z_interpolate = Z_INTERPOLATE_SPLINE; }
                                       else {fprintf (stderr, "Error, unknown 'z_interpolate %s' on line %d\n", t, line_number+1);
				             ierror++;}
                                       free(t);
                                       break;
                                     case 2:
				       t =  yytext2string(yytext,2,ntokens+1);
				       if      (strcmp("O3",t)  ==0 || strcmp("o3",t)  ==0)  mol_id = MOL_O3; 
                                       else if (strcmp("O2",t)  ==0 || strcmp("o2",t)  ==0)  mol_id = MOL_O2; 
                                       else if (strcmp("H2O",t) ==0 || strcmp("h2o",t) ==0)  mol_id = MOL_H2O; 
                                       else if (strcmp("CO2",t) ==0 || strcmp("co2",t) ==0)  mol_id = MOL_CO2; 
                                       else if (strcmp("NO2",t) ==0 || strcmp("no2",t) ==0)  mol_id = MOL_NO2; 
                                       else if (strcmp("BRO",t) ==0 || strcmp("bro",t) ==0)  mol_id = MOL_BRO; 
                                       else if (strcmp("OCLO",t)==0 || strcmp("oclo",t)==0)  mol_id = MOL_OCLO; 
                                       else if (strcmp("HCHO",t)==0 || strcmp("hcho",t)==0)  mol_id = MOL_HCHO; 
                                       else if (strcmp("O4",t)==0   || strcmp("o4",t)==0)    mol_id = MOL_O4; 
                                       else if (strcmp("SO2",t)==0  || strcmp("so2",t)==0)   mol_id = MOL_SO2; 
                                       else if (strcmp("CH4",t)==0  || strcmp("ch4",t)==0)   mol_id = MOL_CH4; 
                                       else if (strcmp("N2O",t)==0  || strcmp("n2o",t)==0)   mol_id = MOL_N2O; 
                                       else if (strcmp("CO",t)==0   || strcmp("co",t)==0)    mol_id = MOL_CO; 
                                       else if (strcmp("N2",t)==0   || strcmp("n2",t)==0)    mol_id = MOL_N2; 
                                       else if (strcmp("T" ,t)==0 || strcmp("t" ,t)==0)  mol_id = -1;
                                       else {fprintf (stderr,"Error, unknown z_interpolate option: %s on line %d\n", t, line_number+1);
				             ierror++;}
                                       free(t);  
                                       t =  yytext2string(yytext,3,ntokens+1);
                                       if (mol_id > 0) {
                                         if      (strcmp("SPLINE",t)==0 || strcmp("spline",t)==0) 
                                           Input.atm.interpol_method_gas[mol_id] = INTERP_METHOD_SPLINE;  
                                         else if (strcmp("LINEAR",t)==0 || strcmp("linear",t)==0) 
                                           Input.atm.interpol_method_gas[mol_id] = INTERP_METHOD_LINEAR; 
                                         else if (strcmp("LOG",t)   ==0 || strcmp("log"   ,t)==0) 
                                           Input.atm.interpol_method_gas[mol_id] = INTERP_METHOD_LOG;
                                         else if (strcmp("LINMIX",t)==0 || strcmp("linmix",t)==0) 
                                           Input.atm.interpol_method_gas[mol_id] = INTERP_METHOD_LINMIX; 
                                         else if (strcmp("LOG_SPLINE",t)==0 || strcmp("log_spline",t)==0) 
                                           Input.atm.interpol_method_gas[mol_id] = INTERP_METHOD_LOG_SPLINE;
                                         else {fprintf (stderr,"Error, unknown z_interpolate option: %s on line %d\n", t, line_number+1);
				               ierror++;}
                                       }
                                       else if (mol_id == -1) {
                                         if      (strcmp("SPLINE",t)==0 || strcmp("spline",t)==0) 
                                           Input.atm.interpol_method_temper = INTERP_METHOD_SPLINE;  
                                         else if (strcmp("LINEAR",t)==0 || strcmp("linear",t)==0) 
                                           Input.atm.interpol_method_temper = INTERP_METHOD_LINEAR; 
                                         else if (strcmp("LOG",t)   ==0 || strcmp("log"   ,t)==0) 
                                           Input.atm.interpol_method_temper = INTERP_METHOD_LOG;
                                         else if (strcmp("LINMIX",t)==0 || strcmp("linmix",t)==0) 
                                           Input.atm.interpol_method_temper = INTERP_METHOD_LINMIX; 
                                         else if (strcmp("LOG_SPLINE",t)==0 || strcmp("log_spline",t)==0) 
                                           Input.atm.interpol_method_temper = INTERP_METHOD_LOG_SPLINE;
                                         else {fprintf (stderr,"Error, unknown z_interpolate option: %s on line %d\n", t, line_number+1);
				               ierror++;}
                                       }
                                       free(t);
                                       break;
                                     default:
                                       fprintf (stderr,"Error, z_interpolate has wrong number of arguments on line %d.\n", line_number+1);
                                       ierror++;
                                     }
				    }
zout{SIGNED_FLOATS_OR_WORDS}|zout_sur{SIGNED_FLOATS_OR_WORDS}|zout_sea{SIGNED_FLOATS_OR_WORDS}  {
                                     /* determine input option */
                                     ntokens = yytext2ntokens(yytext);
                                     s = yytext2string(yytext, 1, ntokens+1);
                                     if      ( strncasecmp("zout"    ,s,4) == 0 && strlen(s) == 4 ) {Input.atm.zout_source = OUTLEVEL_ZOUT_ABOVE_SUR;}
				     else if ( strncasecmp("zout_sur",s,8) == 0 && strlen(s) == 8 ) {Input.atm.zout_source = OUTLEVEL_ZOUT_ABOVE_SUR;}
                                     else if ( strncasecmp("zout_sea",s,8) == 0 && strlen(s) == 8 ) {Input.atm.zout_source = OUTLEVEL_ZOUT_ABOVE_SEA;}
                                     else {fprintf (stderr,"Error, program bug during reading input file: %s on line %d\n", yytext, line_number+1);
				           ierror++;}
                                     free(s);

                                     /* check for keywords */
                                     t = yytext2string(yytext, 2, ntokens+1);
				     if      (strncasecmp("atm_levels",  t,10)==0)            {Input.atm.zout_source = OUTLEVEL_ATM_LEVELS;}
				     else if (strncasecmp("all_levels",  t,10)==0)            {Input.atm.zout_source = OUTLEVEL_ALL_LEVELS;}
				     else if (strncasecmp("model_levels",t,13)==0)            {Input.atm.zout_source = OUTLEVEL_MODEL_LEVELS;}
				     else if (strncasecmp("model_layers",t,13)==0)            {Input.atm.zout_source = OUTLEVEL_MODEL_LAYERS;}
				     else if (strncasecmp("model_levels_and_layers",t,24)==0) {Input.atm.zout_source = OUTLEVEL_MODEL_LEVELS_AND_LAYERS;}
                                     else {
                                       Input.atm.nzout = yytext2ntokens(yytext);
                                       if ( Input.atm.zout_source == OUTLEVEL_ZOUT_ABOVE_SUR) {
                                         Input.atm.zout_sur = calloc (Input.atm.nzout,sizeof(float));
                                         for (i=2;i<=Input.atm.nzout+1; i++) {
				           u = yytext2string(yytext, i, Input.atm.nzout+1);
				           if      (strncasecmp("toa"         ,u, 3)==0) {Input.atm.zout_sur[i-2] = ZOUT_TOA;}
				           else if (strncasecmp("sur"         ,u, 3)==0) {Input.atm.zout_sur[i-2] = ZOUT_SURFACE;}
                                           else if (strncasecmp("cpt"         ,u, 3)==0) {Input.atm.zout_sur[i-2] = ZOUT_CPT;}
                                           else {Input.atm.zout_sur[i-2] = yytext2float(yytext,i,ntokens+1);}
                                           free(u);
                                         }
                                       }
                                       else if ( Input.atm.zout_source == OUTLEVEL_ZOUT_ABOVE_SEA ) {
                                         Input.atm.zout_sea = calloc (Input.atm.nzout,sizeof(float));
                                         for (i=2;i<=Input.atm.nzout+1; i++) {
				           u = yytext2string(yytext, i, Input.atm.nzout+1);
				           if      (strncasecmp("toa"         ,u, 3)==0) {Input.atm.zout_sea[i-2] = ZOUT_TOA;}
				           else if (strncasecmp("sur"         ,u, 3)==0) {Input.atm.zout_sea[i-2] = ZOUT_SURFACE;}
                                           else if (strncasecmp("cpt"         ,u, 3)==0) {Input.atm.zout_sea[i-2] = ZOUT_CPT;}
                                           else {Input.atm.zout_sea[i-2] = yytext2float(yytext,i,ntokens+1);}
                                           free(u);
                                         }
                                       }
                                     }
                                     free(t);

                                     /* check: these keyword do not allow further levels */
                                     if (Input.atm.zout_source == OUTLEVEL_ATM_LEVELS ||
                                         Input.atm.zout_source == OUTLEVEL_ALL_LEVELS ||
                                         Input.atm.zout_source == OUTLEVEL_MODEL_LEVELS ||
                                         Input.atm.zout_source == OUTLEVEL_MODEL_LAYERS ||
                                         Input.atm.zout_source == OUTLEVEL_MODEL_LEVELS_AND_LAYERS) {
                                       if (ntokens != 1 && ntokens != 2) {
                                         fprintf (stderr,"Error, using '%s %s' in line %d accept no further levels!!!\n", s, t, line_number+1);
				         ierror++;
                                       }
                                       if (ntokens == 2)
                                         Input.atm.nzout = yytext2int(yytext,2+1,ntokens+1);
                                     }
				    }

  /* OBSOLETE OPTIONS */
 /*aerosol_sizedist_file{SPACE_WORD}   {fprintf(stderr,"aerosol_sizedist_file is obsolete!"); ierror++;}*/
aerosol_sizedist_file{SPACE_WORD} {
	Input.aer.spec = 1;
	yytextCstring ( Input.aer.filename[FN_AER_SIZ], yytext, 2, 2);
}
 /*aerosol_refrac_index{SPACE_FLOAT}{SPACE_FLOAT}   {fprintf(stderr,"aerosol_refrac_index is obsolete!"); ierror++;} */
aerosol_refrac_index{SPACE_FLOAT}{SPACE_FLOAT} {
	Input.aer.spec = 1;
	Input.aer.re = yytext2float(yytext, 2, 3);
	Input.aer.im = yytext2float(yytext, 3, 3);
}



mc_cloud_grid{SPACE_INTEGER}{SPACE_INTEGER}{SPACE_INTEGER} { fprintf (stderr, " ... option mc_cloud_grid is obsolete!\n"); ierror++;}
profile_layer{SPACE}+{WORD}         {fprintf(stderr,"profile_layer is obsolete!"); ierror++;}
wc_ipa{SPACE_WORD} {fprintf(stderr,"wc_ipa is obsolete !"); ierror++;}
wc_layer{SPACE}+{WORD}		    {fprintf(stderr,"wc_layer is obsolete!"); ierror++;}
wc_saturate{SPACE}+{WORD}+{SPACE_FLOAT}* {fprintf(stderr,"wc_saturate obsolete!"); ierror++;}
ic_layer{SPACE}+{WORD}              {fprintf(stderr,"ic_layer is obsolete!"); ierror++;}
ic_saturate{SPACE}+{WORD}+{SPACE_FLOAT}* {fprintf(stderr,"ic_saturate obsolete!"); ierror++;}
mc_momcld_file{SPACE}+{WORD}        {fprintf (stderr, "Error, option mc_momcld_file does not exist anymore.\n");
                                     fprintf (stderr, "Use wc_properties instead!\n"); exit(-1);}
mc_panorama{SIGNED_FLOATS}          {fprintf (stderr, "Error, option mc_panorama does not exist anymore.\n");
                                     fprintf (stderr, "Use mc_panorama_view instead!\n"); exit(-1);}
mc_pcld_file{SPACE}+{WORD}          {fprintf(stderr,"Error, option 'mc_pcld_file' on line %d is obsolete!\n", line_number+1); ierror++; }
mc_panorama_vr_sc0 		    {fprintf(stderr,"Error, option 'mc_panorama_vr_sc0' on line %d is obsolete!\n", line_number+1); ierror++; }
mc_lidar_abserr{SPACE}+{FLOAT}      {fprintf(stderr,"'mc_lidar_abserr on line %d is obsolete option!\n", line_number+1); ierror++;}
mc_new_delta_scaling {fprintf(stderr,"mc_new_delta_scaling obsolete, use mc_delta_scaling! read documentation for correct usage"); ierror++;}
mc_relerr{SPACE}+{FLOAT}            {fprintf(stderr,"'mc_relerr on line %d is obsolete option!\n", line_number+1); ierror++;}
aerosol_refrac_file{SPACE_WORD}   {fprintf(stderr,"aerosol_refrac_file is obsolete!"); ierror++;}
mc_readrandomseed                   {fprintf(stderr,"mc_readrandomseed is obsolete!"); ierror++;}
header                              {fprintf(stderr,"header is obsolete!"); ierror++;}
optimize_fortran{SPACE_INTEGER}*    {fprintf(stderr,"optimize_fortran is obsolete!"); ierror++;}
time_interpolate{WORDS}             {fprintf(stderr,"time_interpolate is obsolete!"); ierror++;}
time_interval{SIGNED_FLOATS_OR_WORDS} {fprintf(stderr,"time_interval is obsolete!"); ierror++;}



{SPACE}+       ;         /* eat up blanks */ 

{COMMENT}*     ;         /* eat up comments */

\n                        { ++line_number;};

{WORD}                    {fprintf(stderr,"Unknown command on line %d\n", line_number+1);
                           fprintf(stderr,"%s\n", yytext);
                           ierror++;}


<include>[ \t]*      /* eat the whitespace */
<include>[^ \t\n]+   { /* got the include file name */
        if (include_stack_ptr >= MAX_INCLUDE_DEPTH) {
            fprintf (stderr, "Includes nested too deeply");
            exit (1);
        }

        include_stack[include_stack_ptr++] = YY_CURRENT_BUFFER;

        yyin = fopen(yytext, "r");

        if (!yyin) {
            fprintf (stderr, "Error including %s\n", yytext);
            exit(1);
        }

        yy_switch_to_buffer (yy_create_buffer (yyin, YY_BUF_SIZE));

        BEGIN(INITIAL);
        }

<<EOF>> {
        if (--include_stack_ptr < 0) {
            yyterminate();
        }
        else {
            yy_delete_buffer (YY_CURRENT_BUFFER);
            yy_switch_to_buffer (include_stack[include_stack_ptr]);
        }
        }

%%

/*************************************************************/
/* print usage information                                   */
/*************************************************************/

#ifndef UVSPEC_FUNCTION
static void version()
{
  fprintf (stderr, "uvspec, version %s\n", THEVERSION);
}

static void usage()
{

  fprintf (stderr, "\n");  
  version();
  fprintf (stderr, "\nCalculates the radiation field of the Earth\'s atmosphere.\n");
  fprintf (stderr, "See libRadtran documentation for more information.\n\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "Radiative transfer is a complex problem. Make sure\n");
  fprintf (stderr, "you understand both the problem you try to solve and the\n");
  fprintf (stderr, "tools you use to solve it.\n\n");
  fprintf (stderr, "\n");

  fprintf (stderr, "Usage: uvspec [-h] < input_file > Output_file\n");
  fprintf (stderr, " -h             Print this message.\n");  
  fprintf (stderr, " -v             Print version number.\n");  
  fprintf (stderr, " -c             Checks if the syntax of the input file is correct,\n");  
  fprintf (stderr, "                then exit. NO radiation calculated.\n");  
  fprintf (stderr, " -i input_file  Alternative way to specifiy the input file.\n");  
  fprintf (stderr, " -f input_file  Same as -i.\n");  

  fprintf (stderr, "\n");  
}
#endif

#ifndef UVSPEC_FUNCTION
int main(int argc, char **argv) 
#else
int calluvspec (char *outfilename)
#endif
{
  

/* @11c@

% @include mystic.sdoc

<lpdoc>
\section{Radiative transfer tool - \code{uvspec}}
\label{sec:uvspec_options}

The \code{uvspec} input file consists of single line entries, each making
up a complete input to the \code{uvspec} program. First on the line comes
the parameter name, followed by one or more parameter values. The 
parameter name and the parameter values are separated by white space.
Filenames are entered without any surrounding single or double quotes.
Comments are introduced by a \code{\#}. Blank lines are ignored. 
The order of the lines is not important, with one exception: if the
same input option is used more than once, the second one will usually
over-write the first one. Be aware that also options in another \code{included} 
input file will overwrite options specified before.

Since libRadtran version 1.8 input option names have changed.
See section \ref{sec:translate} on how to translate old style input files to new style input files.

The various input parameters are described in detail below.
\begin{description}

\undocumented{
\option{ECHAM_atmosphere_file}
Undocumented option.
}
\undocumented{
\option{ECHAM_clouds_file}
Undocumented option.
}

\undocumented{
\option{ECMWF_atmosphere_file}
Reads in combination with the options \code{latitude},
\code{longitude}, and \code{time} (all mandatory) the pressure,
temperature, ozone, and water vapour from an ECMWF \emph{netCDF} data
file and will combine it with the data given by the
\code{atmosphere_file}.
\fcode{
ECMWF_atmosphere_file file
}
No spatial interpolation of the values is done.  The atmospheric data
nearest to the specified \code{time} will be selected (or linearly
interpolated if \code{time_interpolate} is switched on).  Atmospheric
profiles, which are not provided by the ECMWF file (O2, CO2, NO2) are
taken from the \code{atmosphere_file}. Per default, also the
atmosphere above the ECMWF data is taken from the
\code{atmosphere_file}. In order to avoid this, please have a look at
the option: \code{ECMWF_levels_only}.
}

\undocumented{
\option{ECMWF_ic_reff}
\emph{This option is preliminary and still subject to change!}
The ECMWF data only contains cloud water content, but no effective radius.
With this option, the effective radius can be specified. 
There are two possibilities:
For a fixed effective radius use the keyword \code{fixed} 
and specify the \code{reff} in micrometer.
\fcode{
ECMWF_ic_reff fixed reff
}
In order to use the parametrisation by \citet{Ou1995} use the keyword \code{Ou}.
\fcode{
ECMWF_ic_reff Ou
}
The default option is Ou.
}

\option{no2_column_du}
Obsolete, use \code{mol_modify} instead.
%Set the NO2 column to a given value. The column is integrated between the 
%user-defined \code{altitude} and TOA (top of atmosphere). The value must be in 
%Dobson units. If value is negative or \code{no2_column_du}
%is not specified, the NO2 column is not scaled.

\option{no2_column_moleccm-2}
Obsolete, use \code{mol_modify} instead.
%Set the NO2 column to a given value. The column is integrated between the 
%user-defined \code{altitude} and TOA (top of atmosphere). The value must be in 
%molecules / cm2. If value is negative or \code{no2_column_moleccm-2}
%is not specified, the NO2 column is not scaled.

\option{ozone_column}
Obsolete, use \code{mol_modify} instead.
%Set the ozone column to a given value. The column is integrated between the 
%user-defined \code{altitude} and TOA (top of atmosphere). The value must be in 
%Dobson units. If value is negative or \code{ozone_column}
%is not specified, the ozone column is not scaled.

%CE Commented for ESA version
\undocumented{
\option{refractive_index_file}
Undocumented option to read a refractive index profile. 
\fcode{
refractive_index_file file
}
Actually not the refractive index profile 
but refractive index - 1 must be provided for better accuracy. The file must have two columns: 
\code{altitude [km]   refractive_index-1}
}

\option{time}
Specifies the time to simulate.
\fcode{
time YYYY MM DD hh mm ss
}
where \code{YYYY} is the year, \code{MM} the month, \code{DD} the day, 
\code{hh} the hour, \code{mm} the minute, \code{ss} the second in UTC.
The time information will be used for a couple of things: 

%\code{time} in combination with \code{latitude}, \code{longitude}, and any \code{map}-option
%is used to select the location where to read the input data.

\code{time} is used to correct extraterrestrial 
irradiance for the Sun-Earth distance with the day of year. If not given, 
the Earth-Sun distance is 1 AU (i.e. equinox distance).

\code{time} in combination with \code{latitude} and \code{longitude} is used to calculate 
the solar zenith angle if no \code{sza} is specified.

\code{time} in combination with \code{latitude} and \code{longitude} is used to choose a
suitable default atmosphere file, if no \code{atmosphere_file} is specified.

%\code{time} in combination with an \code{ECMWF_atmosphere_file} is used to choose
%a date in the ECMWF input file. 

%\option{time_interpolate}
%If a map option is used in combination with \code{time}, the data, which is nearest to the 
%specified \code{time} is used for the simulation. This means \code{time_interpolate} is switched 
%off per default.
%
%If this option is switched on, the data fields stored in the \emph{netCDF} files are interpolated 
%to the specified \code{time}. (Be aware, that this might cause strange effects for data field 
%of moving properties. E.g. an interpolated cloud field might have double horizontal extent, 
%but only half the optical depth.)

%\option{time_interval}
%This option can be used in order to calculate an effective solar zenith angle 
%for a time interval, instead of a distinct point in time. 
%The cosine of the solar zenith angle is here replaced by its time average. 
%The azimuth of the sun is replaced by an average of the azimuth position weighted 
%with the cosine of the solar zenith angle.
%\fcode{time_interval dtime_start dtime_end [unit]\\
%time_interval    -180        180     min
%}
%the time interval reaches from \code{time} + \code{dtime_start} to 
%\code{time} + \code{dtime_end}, in the example from 180 minutes before \code{time}
%to 180 minutes after \code{time}. The \code{unit} argument is optional, and can be one of the 
%following: \code{s} (seconds), \code{min} (minutes), or \code{h} (hour). The default is \code{s}.
%This option makes only an effect in combination with \code{time}, \code{latitude}, \code{longitude}, 
%and only has an effect for solar simulations (and of course if no \code{sza} defined).

%DOCUMENTATION END
\end{description} % @c11\_03@ */


/*
</lpdoc>
 */


  /***************** End of documentation *********************/

  int i=0;
  int brdftype=0;
  int status=0;
  int io=0;
  int check_input_file_only = 0;

#ifndef UVSPEC_FUNCTION
#ifdef COPYPROTECTION
  /* copy protection :-) */
  struct tm *tmstruct = calloc (1, sizeof(struct tm));
  time_t t=0;
#endif
#endif

#if HAVE_LIBNETCDF
  int ncid=0;
  int id_var_test=0;
#endif

  char infilename[FILENAME_MAX]="";


  char randomseedfilename [FILENAME_MAX] = "randomseed";
  FILE *randomseedfile=NULL;
#ifdef WRITERANDOMSTATUS
  FILE *randomstatusfile=NULL;
#endif
  int rseed=0;

  char tmp_filename [FILENAME_MAX] = "";
  FILE *f=NULL;

  /* char base[255]=""; */

  int time_UTC_in_s = NOT_DEFINED_INTEGER;
  int time_LAT_in_s = NOT_DEFINED_INTEGER;

  
#ifndef UVSPEC_FUNCTION
  int c=0;
#endif

  char function_name[]="main";
  char file_name[]="uvspec_lex.l";

  /* Initialization */
  
  /* allocate memory for filenames */

  Input.filename         = calloc (FN_NN,     sizeof(char *));
  Input.aer.filename     = calloc (FN_AER_NN, sizeof(char *));
  Input.atm.filename     = calloc (MOL_NN,    sizeof(char *));
  Input.rte.mc.filename  = calloc (FN_MC_NN,  sizeof(char *));
  Output.crs.filename    = calloc (MOL_NN,    sizeof(char *));

  for (i=0; i<FN_NN; i++)
    Input.filename[i]        = calloc (FILENAME_MAX+1, sizeof(char));
  
  for (i=0; i<FN_AER_NN; i++) 
    Input.aer.filename[i]    = calloc (FILENAME_MAX+1, sizeof(char));

  for (i=0; i<MOL_NN; i++) 
    Input.atm.filename[i]    = calloc (FILENAME_MAX+1, sizeof(char));
 
  for (i=0; i<FN_MC_NN; i++)
    Input.rte.mc.filename[i] = calloc (FILENAME_MAX+1, sizeof(char));

  for (i=0; i<MOL_NN; i++) 
    Output.crs.filename[i]   = calloc (FILENAME_MAX+1, sizeof(char));

  Input.ck_scheme_filename   = calloc (FILENAME_MAX+1, sizeof(char));

  s = NULL; 
  s = getenv("LIBRADTRAN_DATA_FILES");
  if ( s == NULL ) { s = "../data/"; }
  strcpy (Input.filename[FN_PATH], s);
  
  if (Input.filename[FN_PATH][strlen(Input.filename[FN_PATH])-1]!='/') {
    strcat(Input.filename[FN_PATH],"/");
  }

  Input.aer.modify = (float **) calloc (MODIFY_VAR_NN, sizeof (float *));
  for (int i=0; i< MODIFY_VAR_NN; i++) {
      Input.aer.modify[i] = (float*) calloc (MODIFY_TYPE_NN, sizeof (float));
      for (int j=0; j< MODIFY_TYPE_NN; j++)	Input.aer.modify[i][j] = NOT_DEFINED_FLOAT;
  }


  Input.quiet             = 0;
  Input.verbose           = 0;
  Input.output_format     = OUTPUT_FORMAT_ASCII;
  Input.header            = 0;
  Input.spectrum_unit     = UNIT_NOT_DEFINED;        /* default: unit of solar spectrum is unknown */
  Input.output_unit       = UNIT_NOT_DEFINED;        /* default: output units are unknown, same as unit of solar spectrum */
  Input.bandwidth         = 1.0;                     /* default: bandwidth in cm-1 (thermal) or 1nm (solar) */
  Input.bandwidth_unit    = UNIT_NOT_DEFINED;        /* cm-1 for thermal, nm for solar */
  Output.print_phi        = 1;                       /* default: print radiances   */
  Input.source		  = SRC_SOLAR;
  Input.ck_scheme         = CK_NONE;
  Input.ck_h2ocont        = 1;                       /* default is h2o-continuum is on, version 2.1 */

  Input.ck_abs 		  = calloc (CK_ABS_NN, sizeof(int));
  for ( i=0; i<CK_ABS_NN; i++ ) {
	Input.ck_abs[i]   = SWITCH_ON;
  }

  Input.n_caoth           = 0;
  Input.i_ic              = -1;
  Input.i_wc              = -1;
  
  Input.cloud_overlap     = CLOUD_OVERLAP_OFF;
  Output.cf.cf_total      = NOT_DEFINED_FLOAT;
  Output.cf.nlev          = 0; 

  Input.rte.mc.spectral_is_wvl = calloc (1, sizeof(double));

  Input.crs_model = (int *) calloc ( CRS_MOL_NN, sizeof (int ));
  Input.crs_model[CRS_MOL_RAYLEIGH] = CRS_MODEL_BODHAINE;
  Input.crs_model[CRS_MOL_O3]       = CRS_MODEL_MOLINA;
  Input.crs_model[CRS_MOL_NO2]      = CRS_MODEL_BURROWS;

  Input.so2_spectrum      = SO2_BREMEN;
  Input.co2_spectrum      = CO2_YOSHINO;
  Input.o2_spectrum       = O2_DEFAULT;
  Input.bro_spectrum      = BRO_IASB;
  Input.oclo_spectrum     = OCLO_WAHNER;
  Input.hcho_spectrum     = HCHO_IASB;
  Input.o4_spectrum       = O4_GREENBLATT;

  Input.UTC.tm_year       = NOT_DEFINED_INTEGER;
  Input.UTC.tm_mon        = NOT_DEFINED_INTEGER;
  Input.UTC.tm_mday       = NOT_DEFINED_INTEGER;
  Input.UTC.tm_hour       = NOT_DEFINED_INTEGER;
  Input.UTC.tm_min        = NOT_DEFINED_INTEGER;
  Input.UTC.tm_sec        = NOT_DEFINED_FLOAT;

  Input.LAT.tm_yday       = NOT_DEFINED_INTEGER;
  Input.LAT.tm_year       = NOT_DEFINED_INTEGER;
  Input.LAT.tm_mon        = NOT_DEFINED_INTEGER;
  Input.LAT.tm_mday       = NOT_DEFINED_INTEGER;
  Input.LAT.tm_hour       = NOT_DEFINED_INTEGER;
  Input.LAT.tm_min        = NOT_DEFINED_INTEGER;
  Input.LAT.tm_sec        = NOT_DEFINED_FLOAT;

  Input.raman               = RAMAN_NONE;
  Input.raman_fast          = 0;
  Input.raman_original      = 0;
  Input.delta_wvl_raman  = -1;

  Input.aer.standard      = 0;
  strcpy(Input.aer.filename[FN_AER_SPECIES_LIB],"OPAC");  /* default aerosol library   */


  Input.latitude          = NOT_DEFINED_FLOAT; 
  Input.longitude         = NOT_DEFINED_FLOAT; 
  Input.sat_pixel_x       = NOT_DEFINED_FLOAT;
  Input.sat_pixel_y       = NOT_DEFINED_FLOAT;

  Input.alb.surface       = NOT_DEFINED_INTEGER;    /* there is no default */
  Input.alb.scale_factor  = 1.0;
  Input.alb.netCDF_alb_name  = calloc (FILENAME_MAX+1, sizeof(char));
  Input.alb.netCDF_surf_name = calloc (FILENAME_MAX+1, sizeof(char));

  Input.flu.source        = NOT_DEFINED_INTEGER;    /* first of all, no defined albedo, default setting see below */
  Input.flu.fluorescence  = 0.0;

  Input.alt.source        = ALT_NOT_DEFINED;
  Input.alt.scale_factor  = 1.0;
  Input.alt.netCDF_alt_name  = calloc (FILENAME_MAX+1, sizeof(char));

  Input.netCDF_name_surf_T = calloc (FILENAME_MAX+1, sizeof(char));

  Input.atm.rs_source     = RS_NO_DATA;                /* default no radiosonde-file */
  Input.atm.nz_atm_forced_sea = NOT_DEFINED_INTEGER;   /* no forced atmosphere level */

  Input.delta_time_start  = 0.0;  /* no time period for averaging */
  Input.delta_time_end    = 0.0;  /* no time period for averaging */
  Input.atm.sza_spher     = 0.0;
  Input.atm.phi0_spher    = 0.0;

  Input.atm.nzout         = NOT_DEFINED_INTEGER;
  Input.atm.zout_source   = OUTLEVEL_ZOUT_ABOVE_SUR;
  Output.atm.nz_model_level = NOT_DEFINED_INTEGER;
  Output.atm.nz_model_layer = NOT_DEFINED_INTEGER;
  Input.atm.z_interpolate = Z_INTERPOLATE_LINMIX; /* alternative default: Z_INTERPOLATE_LOGLIN */
  Input.atm.interpol_method_press  = NOT_DEFINED_INTEGER;
  Input.atm.interpol_method_temper = NOT_DEFINED_INTEGER;
  Input.atm.interpol_method_refind = NOT_DEFINED_INTEGER;   
  Input.atm.interpol_method_gas = calloc(MOL_NN, sizeof(int));
  Input.atm.well_mixed_gas      = calloc(MOL_NN, sizeof(int));
  for (i=0; i<MOL_NN; i++) {
    Input.atm.interpol_method_gas[i] = NOT_DEFINED_INTEGER;
    Input.atm.well_mixed_gas     [i] = NOT_DEFINED_INTEGER;
  }
  Input.atm.ECMWF_new_format = NOT_DEFINED_INTEGER;
  Input.atm.ECMWF_ozone_climatology = FALSE;
  Input.atm.time_interpolate = TIME_NEAREST_DATE;  
  Output.atm.microphys.z_cpt_sea  = NOT_DEFINED_FLOAT; /* cold point tropopause height      */ 
  Output.atm.microphys.temper_cpt = NOT_DEFINED_FLOAT; /* cold point tropopause temperature */ 
  Output.atm.microphys.press_cpt  = NOT_DEFINED_FLOAT; /* cold point tropopause pressure    */ 

  Input.cm.pcl_netCDF_name  = calloc (FILENAME_MAX+1, sizeof(char));
  Input.cm.pcl_scale_factor = 1.0;
  Input.cm.sal_netCDF_name  = calloc (FILENAME_MAX+1, sizeof(char));
  Input.cm.sal_scale_factor = 1.0;

  Input.ambrals.ambrals = calloc( BRDF_AMBRALS_NN, sizeof(float) );
  for ( i=0; i<BRDF_AMBRALS_NN; i++) {
	Input.ambrals.ambrals[i] = NOT_DEFINED_FLOAT;
  }

  Input.cm.param = calloc( BRDF_CAM_NN, sizeof(float) );
  for ( i=0; i<BRDF_CAM_NN; i++) {
	Input.cm.param[i] = NOT_DEFINED_FLOAT;
  }
  Input.cm.param[BRDF_CAM_UPHI] = 0.0;

  Input.rpv.rpv = calloc( BRDF_RPV_NN, sizeof(float) );
  for ( i=0; i<BRDF_RPV_NN; i++) {
	Input.rpv.rpv[i] = 0.0;
  }
  Input.rpv.rpv[BRDF_RPV_SCALE] = 1.0;

  Input.bpdf.type           = BPDF_NONE; 
  Input.bpdf.u10            = NOT_DEFINED_FLOAT;

  Input.rte.disort_icm    = DISORT_ICM_PHASE;
  Input.rte.pseudospherical = 0;
  Input.rte.ibcnd         = 0; 
  Input.rte.maxphi        = 0;
  Input.rte.solver        = SOLVER_DISORT;
  Input.rte.pol_stokes[0] = 'I';
  Input.rte.pol_stokes[1] = 'Q';
  Input.rte.pol_stokes[2] = 'U';
  Input.rte.pol_stokes[3] = 'V';

  Input.rte.mc.spherical = calloc ( DIM_NN, sizeof(int) );
  Input.rte.mc.ris = calloc ( MC_RIS_NN, sizeof(double) );
  for ( i=0; i<MC_RIS_NN;i++) {
	Input.rte.mc.ris[i] = 1.0;
  }
  Input.rte.mc.abs_unit      = MCABS_UNIT_W_PER_M2_AND_DZ;
  Input.rte.mc.absorption    = MCFORWARD_ABS_NONE;
  Input.rte.mc.photons       = 100000;
  Input.rte.mc.delta_scaling_mucut = NOT_DEFINED_FLOAT;
  Input.rte.mc.delta_scaling_start = -1;
  Input.rte.mc.dcirc = NOT_DEFINED_FLOAT;
  Input.rte.mc.ncirc = 0;
  Input.rte.mc.truncate            = NOT_DEFINED_INTEGER;
  Input.rte.mc.bcond         = MCBCOND_PERIODIC;
  Input.rte.mc.randomseed    = NOT_DEFINED_INTEGER;
  Input.rte.mc.escape = 1;
  Input.rte.mc.vroom  = -1;
  Input.rte.mc.sensorposition     = MC_SENSORPOSITION_NONE;
  
  Input.rte.mc.pan = calloc ( PAN_FLAG_NN, sizeof(int) );
  Input.rte.mc.panorama = PAN_MODE_NONE;
  Input.rte.mc.pan_alignment = PAN_ALIGNMENT_NONE;
  Input.rte.mc.nstokes = 1;

  Input.rte.mc.blitz_position = NULL;

  Input.rte.mc.backward.output = MCBACKWARD_EDN;
  
  Input.optimize_fortran  = FALSE;  /* no verbose output for FORTRAN optimization */
  Input.optimize_delta    = 0;      /* no security levels */

  /* Define polradtran defaults for polarized radiance simulations, this is 
  the standard application for the polradtran solver */ 
  Input.rte.pol_quad_type = "G";
  Input.mixing_ratio = calloc(MX_NN, sizeof(float));
  for (i=0; i<MX_NN; i++)
    Input.mixing_ratio[i] = NOT_DEFINED_INTEGER;
  Input.atm.column       = calloc(MOL_NN, sizeof(float));
  for (i=0; i<MOL_NN; i++) 
    Input.atm.column[i] = NOT_DEFINED_INTEGER;
  Input.atm.unit_column           = calloc(MOL_NN, sizeof(int));
  Input.atm.unit_column[MOL_AIR]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_O3]   = MOL_UNIT_DU;
  Input.atm.unit_column[MOL_O2]   = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_H2O]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_CO2]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_NO2]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_BRO]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_OCLO] = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_HCHO] = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_O4]   = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_SO2]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_CH4]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_N2O]  = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_CO]   = MOL_UNIT_CM_2;
  Input.atm.unit_column[MOL_N2]   = MOL_UNIT_CM_2;
  Input.atm.unit_profile           = calloc(MOL_NN, sizeof(int));
  Input.atm.unit_profile[MOL_AIR]  = CM_3;
  Input.atm.unit_profile[MOL_O3]   = CM_3;
  Input.atm.unit_profile[MOL_O2]   = CM_3;
  Input.atm.unit_profile[MOL_CO2]  = CM_3;
  Input.atm.unit_profile[MOL_NO2]  = CM_3;
  Input.atm.unit_profile[MOL_BRO]  = CM_3;
  Input.atm.unit_profile[MOL_OCLO] = CM_3;
  Input.atm.unit_profile[MOL_HCHO] = CM_3;
  Input.atm.unit_profile[MOL_O4]   = CM_3;
  Input.atm.unit_profile[MOL_SO2]  = CM_3;
  Input.atm.unit_profile[MOL_CH4]  = CM_3;
  Input.atm.unit_profile[MOL_N2O]  = CM_3;
  Input.atm.unit_profile[MOL_CO]   = CM_3;
  Input.atm.unit_profile[MOL_N2]   = CM_3;
  Input.atm.well_mixed_gas[MOL_AIR]  = YES;
  Input.atm.well_mixed_gas[MOL_O3]   = NO;
  Input.atm.well_mixed_gas[MOL_O2]   = YES;
  Input.atm.well_mixed_gas[MOL_H2O]  = NO;
  Input.atm.well_mixed_gas[MOL_CO2]  = YES;
  Input.atm.well_mixed_gas[MOL_NO2]  = NO;
  Input.atm.well_mixed_gas[MOL_BRO]  = NO;
  Input.atm.well_mixed_gas[MOL_OCLO] = NO;
  Input.atm.well_mixed_gas[MOL_HCHO] = NO;
  Input.atm.well_mixed_gas[MOL_O4]   = NO;
  Input.atm.well_mixed_gas[MOL_SO2]  = NO;
  Input.atm.well_mixed_gas[MOL_CH4]  = NO;
  Input.atm.well_mixed_gas[MOL_N2O]  = NO;
  Input.atm.well_mixed_gas[MOL_CO]   = NO;
  Input.atm.well_mixed_gas[MOL_N2]   = YES;
  Input.atm.n_rs_gas = 0;
  Input.atm.rs_gas  = calloc(MOL_NN, sizeof(int));
  Input.atm.rs_unit = calloc(MOL_NN, sizeof(int));
  Input.atm.rs_unit[MOL_AIR]  = CM_3;
  Input.atm.rs_unit[MOL_O3]   = MMR;
  Input.atm.rs_unit[MOL_O2]   = CM_3;
  Input.atm.rs_unit[MOL_H2O]  = RH;
  Input.atm.rs_unit[MOL_CO2]  = CM_3;
  Input.atm.rs_unit[MOL_NO2]  = CM_3;
  Input.atm.rs_unit[MOL_BRO]  = CM_3;
  Input.atm.rs_unit[MOL_OCLO] = CM_3;
  Input.atm.rs_unit[MOL_HCHO] = CM_3;
  Input.atm.rs_unit[MOL_O4]   = CM_3;
  Input.atm.rs_unit[MOL_SO2]  = CM_3;
  Input.atm.rs_unit[MOL_CH4]  = CM_3;
  Input.atm.rs_unit[MOL_N2O]  = CM_3;
  Input.atm.rs_unit[MOL_CO]   = CM_3;
  Input.atm.rs_unit[MOL_N2]   = CM_3;
  Input.atm.rs_add_upper_levels = TRUE;
  Input.atm.mol_mass = calloc(MOL_NN, sizeof(float));
  Input.atm.mol_mass[MOL_AIR]  = mol_mass_air; /* 28.977  from function ck.c                    in g/MOL */
  Input.atm.mol_mass[MOL_O3]   = 48.000;       /*         from function ck.c                    in g/MOL */
  Input.atm.mol_mass[MOL_O2]   = 31.998609;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_H2O]  = mol_mass_wv;  /* 18.015  from function ck.c                    in g/MOL */
  Input.atm.mol_mass[MOL_CO2]  = 44.009646;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_NO2]  = 46.005333;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_BRO]  = 95.902831;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_OCLO] = 67.451347;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_HCHO] = 30.026294;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_O4]   = 63.979110;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_SO2]  = 64.062998;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_CH4]  = 16.042941;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_N2O]  = 44.012751;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_CO]   = 28.010342;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */
  Input.atm.mol_mass[MOL_N2]   = 28.013446;    /* http://www.ch.cam.ac.uk/magnus/MolWeight.html in g/MOL */

  Input.r_earth = 6370.0;

  Input.processing         = PROCESS_NONE;
  Input.n_output_user      = NOT_DEFINED_INTEGER;
  Input.output_user_flag   = 0;

  Input.rte.sdisort = calloc( SDISORT_NN, sizeof(int) );
  Input.rte.sdisort[SDISORT_ICHAPMAN] = 1;

  Input.rte.polradtran = calloc( POLRADTRAN_NN, sizeof(int) );
  Input.rte.polradtran[POLRADTRAN_AZIORDER] = 0;
  Input.rte.polradtran[POLRADTRAN_NSTOKES]  = 1;
  Input.rte.polradtran[POLRADTRAN_SRC_CODE] = 1;

  Input.sslidar = calloc( SSLIDAR_NN, sizeof(double) ); 
  Input.sslidar[SSLIDAR_E0]=0.1; /* Joule */
  Input.sslidar[SSLIDAR_POSITION]=0.0; /* km */
  Input.sslidar_nranges=100;
  Input.sslidar[SSLIDAR_RANGE]=0.1; /* km */
  Input.sslidar[SSLIDAR_EFF]=0.5;
  Input.sslidar[SSLIDAR_AREA]=1.0; /* m^2 */
  Input.sslidar_polarisation=0;
 
  Input.ck_reptran_option = REPTRAN_OPTION_NONE;
  Input.ck_reptran_channel  = calloc (FILENAME_MAX+1, sizeof(char));

  Input.absorption = 1;
  Input.molabs = MOLABS_CALC;
  Input.scattering = 1;
  Input.rayleigh = RAYLEIGH_CALC;

  /*PYTHON INITIALISATIONS*/




  
#ifndef UVSPEC_FUNCTION
  while ((c=getopt (argc, argv, "hvcf:i:")) != EOF)  {
    switch(c)  {
    case 'h': 
      usage();
      return (0);
      break;
    case 'c': 
      check_input_file_only = TRUE;
      break;
    case 'f':
    case 'i':
      strcpy (infilename, optarg);
      break;
    case 'v': 
      version();
      return (0);
      break;
    default:
      usage();
      return (0);
    }
  }
#endif

  /* copy protection :-) */
#ifndef UVSPEC_FUNCTION
#ifdef COPYPROTECTION
  tmstruct->tm_mday=1;   /* day of the month 1-31      */
  tmstruct->tm_mon=3;    /* months since January 0-11  */
  tmstruct->tm_year=109; /* years since 1900           */
  
  t=mktime (tmstruct);
 
  fprintf (stderr, "%s\n", asctime (tmstruct));

  if ((int) time(NULL) > t) {
    remove(argv[0]);
    return -1;
  }
#endif  
#endif  

  /***********************/
  /* Read the input file */
  /***********************/

  if (strlen(infilename)>0)
    yyin = fopen (infilename, "r");

  yylex();

  if (ierror!=0)  {
    fprintf (stderr, "Error in input file! Please check your input file and retry!\n");
    return -1;
  }
  /* Everything is fine ininput file. Thus exit if only checking file */
  if (check_input_file_only)  { return 0; }

  if (strlen(infilename)>0)
    fclose (yyin);


  /*****************************************/
  /* check combinition of parameters       */
  /* and do input dependent initialisation */
  /*****************************************/

  /* append slash to data_files_path*/
  if (Input.filename[FN_PATH] != NULL) {
	  if (Input.filename[FN_PATH][strlen(Input.filename[FN_PATH])-1]!= '/')      strcat (Input.filename[FN_PATH], "/");
  }

  /* initialize wavelength end if not specified */
  if ( Input.wl.end < 0 && Input.wl.start > 0 )
	Input.wl.end = Input.wl.start;
  if ( Input.wl.end_index < 0 && Input.wl.start_index > 0 )
	Input.wl.end_index = Input.wl.start_index;

 /* Reverse order of levels of atm-z-grid */
  if ( Input.atm.nz_atm_forced_sea > 0 ) {
    for (i=0;i<Input.atm.nz_atm_forced_sea/2;i++){
      tmpfloat = Input.atm.z_atm_forced_sea[i];
      Input.atm.z_atm_forced_sea[i] = Input.atm.z_atm_forced_sea[Input.atm.nz_atm_forced_sea-1-i];
      Input.atm.z_atm_forced_sea[Input.atm.nz_atm_forced_sea-1-i] = tmpfloat;
    }
  }

  /* initialization of correlated-k scheme */
  if ( Input.ck_scheme == CK_REPTRAN && Input.ck_reptran_arg != NULL ) {
    if (strncasecmp("fine", Input.ck_reptran_arg,4)==0)         {Input.ck_reptran_option = REPTRAN_OPTION_FINE;}
    else if (strncasecmp("medium", Input.ck_reptran_arg,6)==0)  {Input.ck_reptran_option = REPTRAN_OPTION_MEDIUM;}
    else if (strncasecmp("coarse", Input.ck_reptran_arg,6)==0)  {Input.ck_reptran_option = REPTRAN_OPTION_COARSE;}
    else {Input.ck_reptran_option = REPTRAN_OPTION_NONE;}
  }
  
  if ( Input.ck_scheme == CK_REPTRAN_CHANNEL )
     strcpy (Input.ck_reptran_channel, Input.ck_reptran_arg);
  
  free(Input.ck_reptran_arg);

  /* add standard profiles of CH4, N2O, CO, and N2 if representative wavelengths are used */
  if ( Input.ck_scheme == CK_REPTRAN || Input.ck_scheme == CK_REPTRAN_CHANNEL ) {

    if ( strlen(Input.atm.filename[MOL_CH4]) == 0 ) {
      if (strlen(Input.filename[FN_PATH])>0)
        strcpy (Input.atm.filename[MOL_CH4],  Input.filename[FN_PATH]);
      else 
        strcpy (Input.atm.filename[MOL_CH4],  "../data/");
      strcat (Input.atm.filename[MOL_CH4], "atmmod/afglus_ch4_vmr.dat");
      Input.atm.unit_profile[MOL_CH4]=VMR;
      if(!Input.quiet) fprintf(stderr, " ... using CH4 volume mixing ratio profile from US standard atmosphere.\n");
    }

    if ( strlen(Input.atm.filename[MOL_N2O]) == 0 ) {
      if (strlen(Input.filename[FN_PATH])>0)
        strcpy (Input.atm.filename[MOL_N2O],  Input.filename[FN_PATH]);
      else 
        strcpy (Input.atm.filename[MOL_N2O],  "../data/");
      strcat (Input.atm.filename[MOL_N2O], "atmmod/afglus_n2o_vmr.dat");
      Input.atm.unit_profile[MOL_N2O]=VMR;
      if(!Input.quiet) fprintf(stderr, " ... using N2O volume mixing ratio profile from US standard atmosphere.\n");
    }

    if ( strlen(Input.atm.filename[MOL_CO]) == 0 ) {
      if (strlen(Input.filename[FN_PATH])>0)
        strcpy (Input.atm.filename[MOL_CO],  Input.filename[FN_PATH]);
      else 
        strcpy (Input.atm.filename[MOL_CO],  "../data/");
      strcat (Input.atm.filename[MOL_CO], "atmmod/afglus_co_vmr.dat");
      Input.atm.unit_profile[MOL_CO]=VMR;
      if(!Input.quiet) fprintf(stderr, " ... using CO volume mixing ratio profile from US standard atmosphere.\n");
    }

    if ( strlen(Input.atm.filename[MOL_N2]) == 0 ) {
      if (strlen(Input.filename[FN_PATH])>0)
        strcpy (Input.atm.filename[MOL_N2],  Input.filename[FN_PATH]);
      else 
        strcpy (Input.atm.filename[MOL_N2],  "../data/");
      strcat (Input.atm.filename[MOL_N2], "atmmod/afglus_n2_vmr.dat");
      Input.atm.unit_profile[MOL_N2]=VMR;
      if(!Input.quiet) fprintf(stderr, " ... using N2 volume mixing ratio profile from US standard atmosphere.\n");
    }

  }
  
  
 /* check crs_model */
  switch (Input.crs_model[CRS_MOL_RAYLEIGH]) {
	case CRS_MODEL_BODHAINE:
	case CRS_MODEL_NICOLET:
	case CRS_MODEL_PENNDORF:
	case CRS_MODEL_BODHAINE29:
		break;
	default:
		fprintf(stderr, "Unknown 'crs_model' for 'rayleigh'.\nPlease read the documentation for valid crs_models.\n");
		ierror++;
  }
  switch (Input.crs_model[CRS_MOL_O3]) {
	case CRS_MODEL_BASS_AND_PAUR:
	case CRS_MODEL_MOLINA:
	case CRS_MODEL_DAUMONT:
	case CRS_MODEL_BOGUMIL:
		break;
	default:
		fprintf(stderr, "Unknown 'crs_model' for 'o3'.\nPlease read the documentation for valid crs_models.\n");
		ierror++;
  }
  switch (Input.crs_model[CRS_MOL_NO2]) {
	case CRS_MODEL_BOGUMIL:
	case CRS_MODEL_BURROWS:
		break;
	default:
		fprintf(stderr, "Unknown 'crs_model' for 'no2'.\nPlease read the documentation for valid crs_models.\n");
		ierror++;
  }

  /*Calculate latitude and longitude*/
  if (Input.lat_signum)		Input.latitude  = Input.lat_signum * ( Input.lat_degrees + Input.lat_minutes / 60.0 + Input.lat_seconds / 3600.0 );
  if (Input.lon_signum)		Input.longitude = Input.lon_signum * ( Input.lon_degrees + Input.lon_minutes / 60.0 + Input.lon_seconds / 3600.0 );

  /* initialization for albedo_library option */
  if (Input.alb.library != NULL) {
    if (strncasecmp("IGBP", Input.alb.library, 4)==0) {Input.alb.source = ALBEDO_IGBP_LIBRARY;}
    else {
         Input.alb.source = ALBEDO_USER_LIBRARY;
         yytextCstring (Input.filename[FN_ALBEDO_LIB_PATH],Input.alb.library, 2, 2);
         /* append slash */
         if (Input.filename[FN_ALBEDO_LIB_PATH][strlen(Input.filename[FN_ALBEDO_LIB_PATH])-1]!= '/')
	     strcat (Input.filename[FN_ALBEDO_LIB_PATH], "/");
    }
    free(Input.alb.library);
  }
  /* Initialize Input.processing or Input.Output_unit from option output{WORD} */
  if ( Input.output_unit_processing != NULL) {
     if      (strncasecmp(Input.output_unit_processing,"integrate",9) == 0) Input.processing = PROCESS_INT;
     else if (strncasecmp(Input.output_unit_processing,"sum",3)       == 0) Input.processing = PROCESS_SUM;
     else if (strncasecmp(Input.output_unit_processing,"rgbnorm",7)   == 0) Input.processing = PROCESS_RGBNORM;
     else if (strncasecmp(Input.output_unit_processing,"rgb_norm",8)  == 0) Input.processing = PROCESS_RGBNORM;
     else if (strncasecmp(Input.output_unit_processing,"rgb",3)       == 0) Input.processing = PROCESS_RGB;
     else if (strncasecmp(Input.output_unit_processing,"per_nm",6)    == 0) Input.output_unit= UNIT_PER_NM;
     else if (strncasecmp(Input.output_unit_processing,"per_cm-1",8)  == 0) Input.output_unit= UNIT_PER_CM_1;
     else if (strncasecmp(Input.output_unit_processing,"per_band",11) == 0) Input.output_unit= UNIT_PER_BAND;
     else {fprintf(stderr,"Unknown 'output %s'\n", Input.output_unit_processing);
         ierror++;}
     if ((Input.processing==PROCESS_SUM) && (Input.ck_scheme==CK_REPTRAN)) Input.output_unit= UNIT_PER_BAND; /* to make summation over band-integrated quantities */
     free(Input.output_unit_processing);
  }

  /* initialization for albedo_library option */
  if (Input.rpv.library != NULL) {
    if (strncasecmp("IGBP", Input.rpv.library, 4)==0) {Input.rpv.source = RPV_IGBP_LIBRARY;}
    else {
         Input.rpv.source = RPV_USER_LIBRARY;
         yytextCstring (Input.filename[FN_RPV_LIB_PATH],yytext, 1+1, ntokens+1);
         /* append slash */
         if (Input.filename[FN_RPV_LIB_PATH][strlen(Input.filename[FN_RPV_LIB_PATH])-1]!= '/')
	       strcat (Input.filename[FN_RPV_LIB_PATH], "/");
    }
    free(Input.rpv.library);
  }

  /*Initializing Monte Carlo options */

 #if !HAVE_OPENGL
  if (Input.rte.mc.visualize) {
	Input.rte.mc.visualize = 0; 
	fprintf (stderr, "*** Warning, you selected mc_visualize but OpenGL is not available.\n");
	fprintf (stderr, "*** Ignoring mc_visualize - please install OpenGL if you need it.\n");
  }
 #endif

  /* mc_lidar: check if lidar is available */
  if (Input.source == SRC_LIDAR) {
    #if !HAVE_LIDAR
	fprintf(stderr,"Error, unknown 'mc_lidar' or 'mc_radar'.\nPlease check your input file and retry!\n");
	Input.rte.mc.locest = -1;
	ierror++;
    #endif
  }

  /* check and initialize mc_backward */

  if (Input.rte.mc.backward.yes) {
    if (Input.rte.mc.backward.islower != NOT_DEFINED_INTEGER) {
	if (Input.rte.mc.backward.jslower == NOT_DEFINED_INTEGER) {
	  fprintf(stderr,"Error, input option 'mc_backward' need either 0, 2, or 4 integers as arguments!\n");
          ierror++; }
	else {
	  if (Input.rte.mc.backward.isupper != NOT_DEFINED_INTEGER && Input.rte.mc.backward.jsupper == NOT_DEFINED_INTEGER ) {
	    fprintf(stderr,"Error, input option 'mc_backward' need either 0, 2, or 4 integers as arguments!\n");
            ierror++; }
	  else if ( Input.rte.mc.backward.isupper == NOT_DEFINED_INTEGER && Input.rte.mc.backward.jsupper == NOT_DEFINED_INTEGER ) {
		Input.rte.mc.backward.isupper = Input.rte.mc.backward.islower;
		Input.rte.mc.backward.jsupper = Input.rte.mc.backward.jslower;
	  }
	}
    }
  }

  /* check mc_sample_grid */
  
  if ( Input.rte.mc.dx_sample != 0 && Input.rte.mc.dy_sample == 0 ) {
	 fprintf(stderr,"Error, input option 'mc_sample_grid' needs either 2 or 4 arguments!\n");
         ierror++; 
  }

  /* Set mc backward absorption for output */

  switch ( Input.rte.mc.backward.output ) {
    case MCBACKWARD_ACT:
    case MCBACKWARD_ABS:
    case MCBACKWARD_EMIS:
    case MCBACKWARD_HEAT:
	Input.rte.mc.backward.absorption = 1;
	break;
    default:
	break;
  }

  /* if uvspec is told to be verbose, it cannot be quiet at the same time */
  if (Input.verbose)
    Input.quiet = 0;

  /*Check if NetCDF filename for output is provided*/
  switch ( Input.output_format ) {
    case OUTPUT_FORMAT_NETCDF:
    case OUTPUT_FORMAT_SAT_PICTURE:
      if ( strlen(Input.filename[FN_OUTFILE]) == 0 ) {
        fprintf (stderr, "NetCDF-Filename not provided! Please define 'output_file' in your input file!\n");
        ierror++;
      }
      break;
    default:
      break;
  }

  /* Quick fix for mc_bw_umu_file, mc_panorama_view, mc_blitz_position */
  /* Initialize Input.rte.umu and Input.rte.phi ************************/
  /* TODO: This needs to be cleaned up somewhere in the code ***********/

  if ( Input.rte.mc.allocate_umu_and_phi == TRUE ) {
    switch ( Input.rte.numu ) {
      case 0:
	Input.rte.numu   = 1;
	Input.rte.maxumu = 1;
	Input.rte.umu    = calloc(1, sizeof(float));
	break;
      case 1:
	break;
      default:
	fprintf (stderr, "It doesn't make sense to define more than one 'umu' with options 'mc_bw_umu_file', 'mc_panorama_view' or 'mc_blitz_position'. Please check you input file again.\n");
	ierror++;
    }
    switch ( Input.rte.nphi ) {
      case 0:
	Input.rte.nphi   = 1;
	Input.rte.maxphi = 1;
	Input.rte.phi    = calloc(1, sizeof(float));
	break;
      case 1:
	break;
      default:
	fprintf (stderr, "It doesn't make sense to define more than one 'phi' with options 'mc_bw_umu_file', 'mc_panorama_view' or 'mc_blitz_position'. Please check you input file again.\n");
	ierror++;
    }
  }

  /* To calculate radiances, we want to increase nstr to at least 16 */
  if (Input.rte.numu > 0){
    if(Input.rte.solver==SOLVER_POLRADTRAN){ 
      if(!Input.quiet){  
             fprintf (stderr, "*** Warning, you want to calculate radiances with POLRADTRAN. \n");
	     fprintf (stderr, "*** POLRADTRAN uses only the first NSTR moments of the phase \n");
	     fprintf (stderr, "*** function which might result in bad accuracy, especially if \n");
	     fprintf (stderr, "*** the phase function has a strong forward scattering peak. \n");
	     fprintf (stderr, "*** If you are only interested in the scalar intensity, it is \n");
	     fprintf (stderr, "*** strongly recommended to use disort. For calculations with \n");
	     fprintf (stderr, "*** polarization, there is currently no other option, but please \n");
	     fprintf (stderr, "*** be aware, that the radiances might not be very accurate. \n\n");
	     fprintf (stderr, " ... setting polradtran parameters:\n");
	     fprintf (stderr, " ... polradtran_aziorder 4\n");
	     fprintf (stderr, " ... polradtran_quad_type E\n");
            } 
       Input.rte.polradtran[POLRADTRAN_AZIORDER]  = 4;
       Input.rte.pol_quad_type = "E";
       if (Input.rte.nstr<8) {
          if(!Input.quiet)
              fprintf (stderr, " ... setting number_of_streams to 8 for radiance calculation!\n");
	  Input.rte.nstr=8;
        }
       }
     else{
       if (Input.rte.nstr<16) {
          if(!Input.quiet)
             fprintf (stderr, " ... setting number_of_streams to 16 for radiance calculation!\n");
	  Input.rte.nstr=16;
	  }
       }
    }
   
  /* nstr is implemented differently for polradtran, 2*nstr (polradtran)*/
  /* corresponds to nstr (disort) */
  if(Input.rte.solver==SOLVER_POLRADTRAN) 
     Input.rte.nstr*=2;
     
  /* MYSTIC with polarisation requires all 4 Stokes components */  
  if( Input.rte.mc.polarisation==1 )
      Input.rte.mc.nstokes = 4;

  /* if tica is specified, turn on ipa_3d (ulrike) */  
  if( Input.tipa>0 )
      Input.ipa3d = 1;

  if (Input.rte.numu>0 && Input.rte.solver==SOLVER_FDISORT1) {
    fprintf (stderr, "\n");
    fprintf (stderr, "*** Warning, you want to calculate radiances with disort 1.3.\n");
    fprintf (stderr, "*** It is strongly recommended to use disort in this case because\n");
    fprintf (stderr, "*** the older disort 1.3 version uses only the first NSTR moments of the\n");
    fprintf (stderr, "*** phase function which might result in really bad distortions.\n");
    fprintf (stderr, "*** disort is set up to consider at least the first 1024 moments.\n");
    fprintf (stderr, "\n");
  }


  /* defaut zout-level is 0 or -999 for MYSTIC (if no keyword is set) */
  if (Input.atm.nzout == NOT_DEFINED_INTEGER) {
    if (Input.rte.mc.spectral_is || Input.rte.mc.concentration_is){
      fprintf(stderr, "Error, calculations with importance sampling ('mc_spectral_is' and 'mc_aerosol_is') require a  \n");
      fprintf(stderr, "       user defined output altitude. Please specify option 'zout'.\n");
      return -1;
    }  
    switch ( Input.atm.zout_source ) {
    case OUTLEVEL_ZOUT_ABOVE_SUR: 
    case OUTLEVEL_ZOUT_ABOVE_SEA:
    case OUTLEVEL_PRESS:
      Input.atm.nzout    = 1;
      Input.atm.zout_sur = (float *) calloc (Input.atm.nzout, sizeof(float)); 
      if (Input.rte.solver == SOLVER_MONTECARLO )
        Input.atm.zout_sur[0] = -999.0;
      else 
        Input.atm.zout_sur[0] = 0.0;
      break;
    case OUTLEVEL_ATM_LEVELS:
    case OUTLEVEL_ALL_LEVELS:
    case OUTLEVEL_MODEL_LEVELS:
    case OUTLEVEL_MODEL_LAYERS:
    case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
      break;
    default:
      fprintf (stderr, "Error, determining output levels zout.\n");
      fprintf (stderr, "       unknown option %d for zout source! (line %d, function '%s' in '%s') \n", 
                               Input.atm.zout_source, __LINE__, __func__, __FILE__);
      return -1;
    }
  }

  if (Input.rte.solver== SOLVER_MONTECARLO && Input.rte.mc.surfaceparallel) {
    if (Input.atm.zout_sur[0] != -999.0) {
      fprintf (stderr, "*** Warning, mc_surfaceparallel does only make sense for radiation\n");
      fprintf (stderr, "*** calculations at the surface. Ignoring mc_surfaceparallel!\n");
      Input.rte.mc.surfaceparallel=0;
    }

    if (Input.atm.zout_sur[0] == -999.0 && Input.atm.nzout>1) {
      fprintf (stderr, "*** Warning, mc_surfaceparallel affects only radiation at the surface.\n");
      fprintf (stderr, "*** Switching off all higher output levels for safety reasons.\n");
      Input.atm.nzout=1;
    }
  }      

  if (Input.rte.solver== SOLVER_MONTECARLO && Input.alt.altitude != NOT_DEFINED_FLOAT) {
    fprintf (stderr, "Error, option altitude does not work with\n");
    fprintf (stderr, "       solver montecarlo! Use mc_elevation_file!\n");
    return -1;
  }
  
  if (Input.rte.solver == SOLVER_MONTECARLO && Input.source == SRC_THERMAL) {
       if (!Input.quiet)
            fprintf(stderr, " ... mc_surface_reflectalways disabled for thermal radiative transfer.\n");    
       Input.rte.mc.reflectalways = 0;
  }
  
  /* surface temperature, this is == NOT_DEFINED_FLOAT, if not specified in the input file */
  Output.surface_temperature = Input.surface_temperature;

  /* define default output as a user defined output */
  if (Input.n_output_user < 0) {
    if (Input.heating == HEAT_NONE) {
      switch (Input.rte.solver) {
      case SOLVER_FTWOSTR:
      case SOLVER_TWOSTR:
      case SOLVER_RODENTS:
      case SOLVER_TWOSTREBE:
        Input.n_output_user = 5;
        Input.output_user   = calloc (Input.n_output_user, sizeof (int));
        Input.output_user[0] = OUTPUT_USER_WAVE;
        Input.output_user[1] = OUTPUT_USER_EDIR;
        Input.output_user[2] = OUTPUT_USER_EDN;
        Input.output_user[3] = OUTPUT_USER_EUP;
        Input.output_user[4] = OUTPUT_USER_U;
	break;
      case SOLVER_FDISORT1:
      case SOLVER_FDISORT2:
      case SOLVER_DISORT:
      case SOLVER_SDISORT:
      case SOLVER_SPSDISORT:
      case SOLVER_SOS:
      case SOLVER_MONTECARLO:
      case SOLVER_SSS:
      case SOLVER_TZS:
      case SOLVER_SSSI:
      case SOLVER_NULL:
        Input.n_output_user = 7;
        Input.output_user = calloc (Input.n_output_user, sizeof (int));
        Input.output_user[0] = OUTPUT_USER_WAVE;
        Input.output_user[1] = OUTPUT_USER_EDIR;
        Input.output_user[2] = OUTPUT_USER_EDN;
        Input.output_user[3] = OUTPUT_USER_EUP;
        Input.output_user[4] = OUTPUT_USER_UDIR;
        Input.output_user[5] = OUTPUT_USER_UDN;
        Input.output_user[6] = OUTPUT_USER_UUP;
        break;
      case SOLVER_SSLIDAR:
        Input.n_output_user  = 3 + Input.sslidar_polarisation;
        Input.output_user    = calloc (Input.n_output_user, sizeof (int));
        Input.output_user[0] = OUTPUT_USER_ZOUT_SUR;
        Input.output_user[1] = OUTPUT_USER_SSLIDAR_NPHOT;
        Input.output_user[2] = OUTPUT_USER_SSLIDAR_RATIO;
        if (Input.sslidar_polarisation)
          Input.output_user[3] = OUTPUT_USER_SSLIDAR_NPHOT_Q;
	break;
      case SOLVER_POLRADTRAN:
        /* polradtran has its special output function, do nothing here */
        break;
      default:
        fprintf (stderr, "Error, unknown rte_solver %d (in uvspec_lex.l)\n", Input.rte.solver);
        return -2;
      }
    }
    else { 

      /* default output for spectral heating rates, as a 2D matrix */
      if (Input.processing != PROCESS_INT && Input.processing != PROCESS_SUM) {

        if( Input.output_format == OUTPUT_FORMAT_ASCII)
          Input.output_format = OUTPUT_FORMAT_MATRIX;

        Input.n_output_user = 3;
        Input.output_user = calloc (Input.n_output_user, sizeof (int));
        Input.output_user[0] = OUTPUT_USER_WAVE;
        switch (Input.atm.zout_source) {
        case OUTLEVEL_ZOUT_ABOVE_SUR:
          Input.output_user[1] = OUTPUT_USER_ZOUT_SUR;
          break;
        case OUTLEVEL_ZOUT_ABOVE_SEA:
        case OUTLEVEL_ATM_LEVELS:
        case OUTLEVEL_ALL_LEVELS:
        case OUTLEVEL_MODEL_LEVELS:
	case OUTLEVEL_MODEL_LAYERS:
	case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
          Input.output_user[1] = OUTPUT_USER_ZOUT_SEA;
          break;
        case OUTLEVEL_PRESS:
          Input.output_user[1] = OUTPUT_USER_P;
          break;
        default:
          fprintf (stderr, "Error, unknown zout marker %d when defining default output (uvspec_lex.l)\n", Input.atm.zout_source);
          return -2;
        }
        Input.output_user[2] = OUTPUT_USER_HEAT;
      }
      else { /* integrated output */
        Input.n_output_user = 2;
        Input.output_user = calloc (Input.n_output_user, sizeof (int));
        switch (Input.atm.zout_source) {
        case OUTLEVEL_ZOUT_ABOVE_SUR:
          Input.output_user[0] = OUTPUT_USER_ZOUT_SUR;
          break;
        case OUTLEVEL_ZOUT_ABOVE_SEA:
        case OUTLEVEL_ATM_LEVELS:
        case OUTLEVEL_ALL_LEVELS:
        case OUTLEVEL_MODEL_LEVELS:
	case OUTLEVEL_MODEL_LAYERS:
	case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
          Input.output_user[0] = OUTPUT_USER_ZOUT_SEA;
          break;
        case OUTLEVEL_PRESS:
          Input.output_user[0] = OUTPUT_USER_P;
          break;
        default:
          fprintf (stderr, "Error, unknown zout marker %d when defining default output (uvspec_lex.l)\n", Input.atm.zout_source);
          return -2;
        }
        Input.output_user[1] = OUTPUT_USER_HEAT;
      }
    }
  }

  /* if output_user requires heating_rate -> set heating_rate default calculation method */
  if (Input.heating == HEAT_NONE) {
    for (io=0; io<Input.n_output_user; io++) {
      if (Input.output_user[io] == OUTPUT_USER_HEAT) { 
        Input.heating = HEAT_LAYER_CD;
        if (!Input.quiet)
          fprintf (stderr, " ... heating rates calculated with central differences\n");
      }
    }
  }

  /* if output_user any horizontal derivations -> atmosphere maps and wind file are mandatory */
  for (io=0; io<Input.n_output_user; io++) {
    if ( Input.output_user[io] == OUTPUT_USER_HEAT_AD_X || Input.output_user[io] == OUTPUT_USER_DTDX ||
         Input.output_user[io] == OUTPUT_USER_HEAT_AD_Y || Input.output_user[io] == OUTPUT_USER_DTDY ||
         Input.output_user[io] == OUTPUT_USER_HEAT_AD ) {
      if ( (strlen(Input.filename[FN_ECMWF]) == 0 && strlen(Input.filename[FN_ECHAM]) == 0) ) {
        fprintf (stderr, "Error. for output_user dTdx, dTdy, or any heat advection  \n");
        fprintf (stderr, "       An 3D atmosphere file in netCDF format are nessesary \n");
        return -1;
      }
    }
  }

  /* if output_user any horizontal derivations -> atmosphere maps and wind file are mandatory */
  for (io=0; io<Input.n_output_user; io++) {
    if ( Input.output_user[io] == OUTPUT_USER_HEAT_AD_X ||
         Input.output_user[io] == OUTPUT_USER_HEAT_AD_Y ||
         Input.output_user[io] == OUTPUT_USER_HEAT_AD ) {
      if ( strlen(Input.filename[FN_ECMWF_WIND_MAP]) == 0 ) {
        fprintf (stderr, "Error. for output_user heat_ad_x, heat_ad_y, or heat_ad  \n");
        fprintf (stderr, "       a wind file in netCDF format are nessesary \n");
        return -1;
      }
    }
  }

  /* activate z_interpolate for heating rate calculations */

  /* COMMENTED THIS OUT! We don't really want that because adding a zout would affect */
  /* all results - this is against the redistribution philosophy. It was introduced   */
  /* in revision 858 back in 2005 in order to prevent funny wiggles in the heating    */
  /* rate; these occur because the optical properties are assumed constant in the     */
  /* layers defined in the atmosphere_file; these wiggles occur only if the zout's    */
  /* are defined at higher resolution than the atmosphere_file and they vanish if     */
  /* zout_interpolate is specified. It should be left to the user to define an        */
  /* appropriate grid with atm_z_grid.                                                */

  /*
  if (Input.heating != HEAT_NONE && Input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) {
    Input.atm.zout_interpolate = ZOUT_INTERPOLATE;
    if (!Input.quiet)
      fprintf (stderr, " ... zout_interpolate is turned on by default for heating rates calculations\n");
  }
  */
  
  /* activate zout_interpolate for calculations with rte_solver rodents (i.e. Robert Buras' twostream model;
                                                                         ulrike, June 2010) */
  
  /* tipa dir works only with rte_solver rodents (ulrike, June 2010) */
  if (Input.tipa == TIPA_DIR && Input.rte.solver != SOLVER_RODENTS) {
    fprintf (stderr, "Error: tipa dir works only with rte_solver rodents. Thus, choose rodents or mc_tipa dir!\n");
    return -1;
  }
  
  /* activate zout_interpolate for calculations with rte_solver twostrebe (i.e. Bernhard Mayers twostream model;
                                                           ulrike, June 2010) */
  if (Input.rte.solver == SOLVER_TWOSTREBE && Input.atm.zout_interpolate == NO_ZOUT_INTERPOLATE) {
    Input.atm.zout_interpolate = ZOUT_INTERPOLATE;
    if (!Input.quiet)
      fprintf (stderr, " ... zout_interpolate is turned on by default for calculations with rte_solver twostrebe\n");
  }

  /* turn on all zout levels to get range bins for sslidar        */
  /* this is a trick to allocate the memory needed for the output */
  if (Input.rte.solver == SOLVER_SSLIDAR)
    Input.atm.zout_source = OUTLEVEL_ALL_LEVELS;

  if (Input.sslidar_polarisation)
    Input.rte.polradtran[POLRADTRAN_NSTOKES] = 2;


  /* check, if only time or day_of_year are used */
  if (Input.UTC.tm_yday > 0 && Input.UTC.tm_mday > 0) {
    fprintf (stderr, "Error. Either 'day_of_year' or 'time' can be specified  \n");
    fprintf (stderr, "       in the input file, but not both. Please remove one of these lines. \n");
    return -1; 
  }

  /* copy date to uninitialised UTC structure, if LAT input */
  if (Input.time_zone == TIME_LAT) {
    Input.UTC.tm_year = Input.LAT.tm_year;
    Input.UTC.tm_mon  = Input.LAT.tm_mon;
    Input.UTC.tm_mday = Input.LAT.tm_mday;
  }

  /* calculate day of year */
  if ( Input.UTC.tm_mday > 0. && (Input.UTC.tm_mon + 1) > 0. ) { /* +1 == C-standard Jan=0, Feb=1, ...*/
    if ((Input.UTC.tm_yday = day_of_year(Input.UTC.tm_mday, (Input.UTC.tm_mon+1))) < 0)  { /* +1 == C-standard Jan=0, Feb=1, ...*/
      fprintf (stderr, "\n*** Error, calculating day_of_year, (in uvspec_lex.l)");
      fprintf (stderr, "\n    Date out of range! day = %d, month = %d\n", Input.UTC.tm_mday, Input.UTC.tm_mon+1);
      return -1;
    }
  }
  Input.LAT.tm_yday = Input.UTC.tm_yday;
  

  /* convert local apperent time to UTC, needs correct day_of_year */
  switch (Input.time_zone) {
  case TIME_UTC:
    /* do nothing already UTC */
    /* fprintf (stderr, " *** UTC! ***  year=%4d, month=%2d, day=%2d, hour=%2d, min=%2d, sec=%2d\n", */
    /*                  Input.UTC.tm_year+1900, Input.UTC.tm_mon+1, Input.UTC.tm_mday, Input.UTC.tm_hour, Input.UTC.tm_min, Input.UTC.tm_sec ); */
    break;
  case TIME_LAT:

    if (Input.verbose) {
      fprintf (stderr, " ... Convert local time to UTC\n"); 
      fprintf (stderr, "     local apparent time:        year=%4d, month=%2d, day=%2d, hour=%2d, min=%2d, sec=%2d, longitude = %7.2f \n",
                     Input.LAT.tm_year+1900, Input.LAT.tm_mon+1, Input.LAT.tm_mday, Input.LAT.tm_hour, Input.LAT.tm_min, Input.LAT.tm_sec, Input.longitude );
    }

    time_LAT_in_s = (int) ( 3600*Input.LAT.tm_hour + 60*Input.LAT.tm_min + Input.LAT.tm_sec );

    /* convert LAT to UTC, last argument standard longitude of UTC, ATTENTION: for function LAT() POSITIVE longitude is WEST */
    time_UTC_in_s = standard_time ( time_LAT_in_s, Input.LAT.tm_yday, -Input.longitude, 0.0); 

    /* correct negative UTCs and UTCs larger than 24:00 */
    /* to get human-readable output                     */ 
    while (time_UTC_in_s<0)
      time_UTC_in_s += 86400; 
    while (time_UTC_in_s>86400)
      time_UTC_in_s -= 86400;

    Input.UTC.tm_hour =  time_UTC_in_s / 3600;
    Input.UTC.tm_min  = (time_UTC_in_s - Input.UTC.tm_hour*3600) / 60;
    Input.UTC.tm_sec  = (time_UTC_in_s - Input.UTC.tm_hour*3600 - Input.UTC.tm_min*60);
           
    if (Input.verbose) {
      fprintf (stderr, "     universal time correlated:  year=%4d, month=%2d, day=%2d, hour=%2d, min=%2d, sec=%2d\n",
                     Input.UTC.tm_year+1900, Input.UTC.tm_mon+1, Input.UTC.tm_mday, Input.UTC.tm_hour, Input.UTC.tm_min, Input.UTC.tm_sec );
    }

    break;
  default:
    fprintf (stderr, "Error, unknown time_zone %d!\n", Input.time_zone);
    fprintf (stderr, "This is a program bug, please contact the programmers!\n");
    return -1;
  }

  /* check that both satellite_pixel and satellite_geometry are activated */
  if (Input.sat_pixel_x > 0 && strlen(Input.filename[FN_SATGEOMETRY]) == 0) {
    fprintf (stderr, "Error, satellite_pixel is specified, but no satellite_geometry is given\n");
    return -1; 
  }
  if (Input.sat_pixel_x < 0 && strlen(Input.filename[FN_SATGEOMETRY]) > 0) {
    fprintf (stderr, "Error, satellite_geometry is given, but no satellite_pixel is specified\n");
    return -1; 
  }


  /* set default SATELLITE GEOMETRY map for MSG */
  if (strncasecmp("MSG",Input.filename[FN_SATGEOMETRY],3)==0 && strlen(Input.filename[FN_SATGEOMETRY]) == 3) {
    strcpy (Input.filename[FN_SATGEOMETRY], Input.filename[FN_PATH]);
    strcat (Input.filename[FN_SATGEOMETRY], "satellites/MSG_seviri/MSG_seviri_geometry.nc");
    if (Input.verbose)
      fprintf (stderr, " ... using MSG satellite geometry %s\n", Input.filename[FN_SATGEOMETRY] );
  }
  if (strlen(Input.filename[FN_SATGEOMETRY]) > 0 && Input.output_format == OUTPUT_FORMAT_NETCDF && !Input.quiet)
    fprintf (stderr, "*** Warning, using satellite_geometry, but not 'output_format sat_picture'\n");

  /* initialise satellite geometry */
  if (strlen(Input.filename[FN_SATGEOMETRY]) > 0) {
    /* read satellite geometry: lat/lon/theta_sat/phi_sat */
    status = read_sat_geometry (Input.sat_pixel_x, Input.sat_pixel_y, Input.UTC, Input.filename[FN_SATGEOMETRY], 
                                &(Input.latitude), &(Input.longitude), 
                                &(Input.rte.numu),&(Input.rte.maxumu), &(Input.rte.umu),
                                &(Input.rte.nphi),&(Input.rte.maxphi), &(Input.rte.phi),
                                Input.rte.solver, Input.verbose, Input.quiet);    /* in ancillary.c */
    if (status!=0) {
      fprintf (stderr, "Error %d reading satellite geometry from %s in uvspec_lex.l\n", status, Input.filename[FN_SATGEOMETRY]);
      return status;
    }
    if (Input.latitude == 999) {
      /* if (!Input.quiet) */
      /*  fprintf (stderr, "Pixel outside Earth plate, skip this calculation\n"); */
      return 0;
    }
  }

  /* if no atmosphere file is given, check if we can find a appropriate atmosphere */
   if (Input.atmosphere_filename != NULL) { 
	strcpy (Input.filename[FN_ATMOSPHERE], Input.atmosphere_filename); 
	free( Input.atmosphere_filename); } 
  if (strlen(Input.filename[FN_ATMOSPHERE]) == 0) {
    /* if lat, lon AND time is defined than search suitable atmosphere file */
    if ( ( fabs(Input.latitude  - NOT_DEFINED_FLOAT) > EPSILON ) &&
         ( fabs(Input.longitude - NOT_DEFINED_FLOAT) > EPSILON ) &&
         ( Input.UTC.tm_mday > 0 ) ) {

      if ( fabs (Input.latitude) < 23.0 )
        strcpy (Input.filename[FN_ATMOSPHERE], "tropics");
      else if ( fabs (Input.latitude) < 67.0 )
        strcpy (Input.filename[FN_ATMOSPHERE],"midlatitude");
      else 
        strcpy (Input.filename[FN_ATMOSPHERE],"subarctic");

      if ( fabs (Input.latitude) >= 23.0 ) { 
        if ( ( Input.latitude > 0 &&  (Input.UTC.tm_mon+1) >=  4 && (Input.UTC.tm_mon+1) <=  9  ) || 
             ( Input.latitude < 0 && ((Input.UTC.tm_mon+1) >= 10 || (Input.UTC.tm_mon+1) <=  3 )) )
          strcat (Input.filename[FN_ATMOSPHERE], "_summer");
        else 
          strcat (Input.filename[FN_ATMOSPHERE], "_winter");
      }
    }
    else /* not all nessesary information are given, choose standard atmosphere */
      strcpy (Input.filename[FN_ATMOSPHERE],"US-standard");

    if (Input.verbose)
      fprintf (stderr, " ... choosing atmosphere file %s\n",Input.filename[FN_ATMOSPHERE]);
  }
  


  /* replace keywords by atmosphere path + files */
  if      (strncasecmp("subarctic_winter",  Input.filename[FN_ATMOSPHERE],16)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglsw.dat");
  }
  else if (strncasecmp("subarctic_summer",  Input.filename[FN_ATMOSPHERE],16)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglss.dat");
  }
  else if (strncasecmp("midlatitude_summer",Input.filename[FN_ATMOSPHERE],18)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglms.dat");
  }
  else if (strncasecmp("midlatitude_winter",Input.filename[FN_ATMOSPHERE],18)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglmw.dat");
  }
  else if (strncasecmp("tropics",           Input.filename[FN_ATMOSPHERE], 7)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglt.dat");
  }
  else if (strncasecmp("US-standard",       Input.filename[FN_ATMOSPHERE],11)==0) {
    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.filename[FN_ATMOSPHERE],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.filename[FN_ATMOSPHERE],  "../data/");
    strcat (Input.filename[FN_ATMOSPHERE], "atmmod/afglus.dat");
  }

  /* try to attach default path for solar_flux file, if file does not exit */
  if (strlen(Input.filename[FN_EXTRATERRESTRIAL]) > 0) {
    if ((f = fopen(Input.filename[FN_EXTRATERRESTRIAL], "r")) == NULL) {

      /* file does not exit, try to add an default path */
      if (strlen(Input.filename[FN_PATH])>0)
        strcpy (tmp_filename,  Input.filename[FN_PATH]);
      else
        strcpy (tmp_filename,  "../data/");
      strcat (tmp_filename, "solar_flux/");
      strcat (tmp_filename, Input.filename[FN_EXTRATERRESTRIAL]);

      /* try to open file with default path */
      if ((f = fopen(tmp_filename, "r")) == NULL) {
        /* (default path + file) does not exits, also */
        /* final check if the solar spectrum is needed is done in */
        /* setup_extraterrestrial() (in extraterrestrial.c) */
      }
      else {
        /* (default path + file) does exit => copy filename */
        fclose(f);
        strcpy (Input.filename[FN_EXTRATERRESTRIAL],tmp_filename);
      }
    }
    else
      fclose(f);
  }

  /* CE: Generate aerosol default filenames and set aerosol types if standard OPAC aerosol is specified */
  /* May be this should moved somewhere else ??? */
  if(Input.aer.n_species==-1){

    if (strlen(Input.filename[FN_PATH])>0)
      strcpy (Input.aer.filename[FN_AER_SPECIES],  Input.filename[FN_PATH]);
    else 
      strcpy (Input.aer.filename[FN_AER_SPECIES],  "../data/");

   if (strncasecmp(Input.aer.mixture_name,"continental_clean",17)==0){ 
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/continental_clean.dat");
      Input.aer.n_species=4;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="suso";
    }
    else if (strncasecmp(Input.aer.mixture_name,"continental_average",19)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/continental_average.dat");
      Input.aer.n_species=4;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="suso";
    }
    else if (strncasecmp(Input.aer.mixture_name,"continental_polluted",20)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/continental_polluted.dat");
      Input.aer.n_species=4;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"maritime_clean",14)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/maritime_clean.dat");
      Input.aer.n_species=6;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";   
      Input.aer.species_names[3]="ssam";
      Input.aer.species_names[4]="sscm";
      Input.aer.species_names[5]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"maritime_tropical",17)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/maritime_tropical.dat");
      Input.aer.n_species=6;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";   
      Input.aer.species_names[3]="ssam";
      Input.aer.species_names[4]="sscm";
      Input.aer.species_names[5]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"maritime_polluted",17)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/maritime_polluted.dat");
      Input.aer.n_species=6;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";   
      Input.aer.species_names[3]="ssam";
      Input.aer.species_names[4]="sscm";
      Input.aer.species_names[5]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"urban",5)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/urban.dat");
      Input.aer.n_species=4;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"desert_spheroids",16)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/desert_spheroids.dat");
      Input.aer.n_species=7;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="minm_spheroids";
      Input.aer.species_names[4]="miam_spheroids";
      Input.aer.species_names[5]="micm_spheroids";
      Input.aer.species_names[6]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"desert",6)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/desert.dat");
      Input.aer.n_species=7;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="minm";
      Input.aer.species_names[4]="miam";
      Input.aer.species_names[5]="micm";
      Input.aer.species_names[6]="suso";
    }  
    else if (strncasecmp(Input.aer.mixture_name,"antarctic",9)==0){
      strcat(Input.aer.filename[FN_AER_SPECIES],  "aerosol/OPAC/standard_aerosol_files/antarctic.dat");
      Input.aer.n_species=6;
      status = ASCII_calloc_char(&Input.aer.species_names,Input.aer.n_species,FILENAME_MAX+1);
      Input.aer.species_names[0]="inso";
      Input.aer.species_names[1]="waso";
      Input.aer.species_names[2]="soot";
      Input.aer.species_names[3]="ssam";
      Input.aer.species_names[4]="mitr";
      Input.aer.species_names[5]="suso";
    }
    else{
      fprintf(stderr, "Error, standard aerosol mixture %s not definded.\n", Input.aer.filename[FN_AER_SPECIES]);
      return -1; 
    }
  } 
  
  /* specify the individual interpolation methods for all gases*/
  /* !!! ATTENTION !!!: INTERP_METHOD_LINMIX not possible for press, temper, refind, rh, MOL_AIR !!! */ 
  switch (Input.atm.z_interpolate) {
  case Z_INTERPOLATE_LOGLIN:    
    if (Input.atm.interpol_method_press  < 0) 
        Input.atm.interpol_method_press  = INTERP_METHOD_LOG;    /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_temper < 0)
        Input.atm.interpol_method_temper = INTERP_METHOD_LINEAR; 
    if (Input.atm.interpol_method_refind < 0)
        Input.atm.interpol_method_refind = INTERP_METHOD_LINEAR;

    if (Input.atm.interpol_method_gas[MOL_AIR]  < 0)
        Input.atm.interpol_method_gas[MOL_AIR]  = INTERP_METHOD_LOG;    /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_gas[MOL_O3]   < 0)
        Input.atm.interpol_method_gas[MOL_O3]   = INTERP_METHOD_LINEAR; /* ozone is assumed to be not well-mixed which */
    if (Input.atm.interpol_method_gas[MOL_O2]   < 0)    
        Input.atm.interpol_method_gas[MOL_O2]   = INTERP_METHOD_LOG;    /* is certainly true for the stratosphere      */
    if (Input.atm.interpol_method_gas[MOL_H2O]  < 0)
        Input.atm.interpol_method_gas[MOL_H2O]  = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_CO2]  < 0)
        Input.atm.interpol_method_gas[MOL_CO2]  = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_NO2]  < 0)
        Input.atm.interpol_method_gas[MOL_NO2]  = INTERP_METHOD_LINEAR; /* NO2 is assumed to be not well-mixed which   */
    if (Input.atm.interpol_method_gas[MOL_BRO]  < 0)
        Input.atm.interpol_method_gas[MOL_BRO]  = INTERP_METHOD_LOG;    /* is certainly true for the stratosphere      */ 
    if (Input.atm.interpol_method_gas[MOL_OCLO] < 0)
        Input.atm.interpol_method_gas[MOL_OCLO] = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_HCHO] < 0)
        Input.atm.interpol_method_gas[MOL_HCHO] = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_O4]   < 0)
        Input.atm.interpol_method_gas[MOL_O4]   = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_SO2]  < 0)
        Input.atm.interpol_method_gas[MOL_SO2]  = INTERP_METHOD_LINEAR; /* SO2 is assumed to be not well-mixed which   */
    if (Input.atm.interpol_method_gas[MOL_CH4]  < 0)
        Input.atm.interpol_method_gas[MOL_CH4]  = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_N2O]  < 0)
        Input.atm.interpol_method_gas[MOL_N2O]  = INTERP_METHOD_LOG;    
    if (Input.atm.interpol_method_gas[MOL_CO]   < 0)
        Input.atm.interpol_method_gas[MOL_CO]   = INTERP_METHOD_LINEAR; /* CO is assumed to be not well-mixed which   */
    if (Input.atm.interpol_method_gas[MOL_N2]   < 0)
        Input.atm.interpol_method_gas[MOL_N2]   = INTERP_METHOD_LOG;    
    break;
  case Z_INTERPOLATE_LINMIX:
    /* LINMIX means interpolation assuming linear mixing ratio for gases O3, ..., HCLO while       */
    /* air density is assumed to vary logarithically                                               */ 
    /* densities profiles interpolated             with LINMIX  */
    /* columns are summed up                       with LINMIX  */ 
    /* optical properties per layer are calculated with LINMIX  */
    /* in general: INTERP_METHOD_LINMIX not possible for press, temper, refind, rh, MOL_AIR        */  

    if (Input.atm.interpol_method_press  < 0) 
        Input.atm.interpol_method_press  = INTERP_METHOD_LOG;    /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_temper < 0)
        Input.atm.interpol_method_temper = INTERP_METHOD_LINEAR; 
    if (Input.atm.interpol_method_refind < 0)
        Input.atm.interpol_method_refind = INTERP_METHOD_LINEAR;

    if (Input.atm.interpol_method_gas[MOL_AIR]  < 0)
        Input.atm.interpol_method_gas[MOL_AIR]  = INTERP_METHOD_LOG;    /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_gas[MOL_O3]   < 0)
        Input.atm.interpol_method_gas[MOL_O3]   = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_O2]   < 0)    
        Input.atm.interpol_method_gas[MOL_O2]   = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_H2O]  < 0)
        Input.atm.interpol_method_gas[MOL_H2O]  = INTERP_METHOD_LINMIX;    
    if (Input.atm.interpol_method_gas[MOL_CO2]  < 0)
        Input.atm.interpol_method_gas[MOL_CO2]  = INTERP_METHOD_LINMIX;    
    if (Input.atm.interpol_method_gas[MOL_NO2]  < 0)
        Input.atm.interpol_method_gas[MOL_NO2]  = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_BRO]  < 0)
        Input.atm.interpol_method_gas[MOL_BRO]  = INTERP_METHOD_LINMIX;  
    if (Input.atm.interpol_method_gas[MOL_OCLO] < 0)
        Input.atm.interpol_method_gas[MOL_OCLO] = INTERP_METHOD_LINMIX;    
    if (Input.atm.interpol_method_gas[MOL_HCHO] < 0)
        Input.atm.interpol_method_gas[MOL_HCHO] = INTERP_METHOD_LINMIX;  
    if (Input.atm.interpol_method_gas[MOL_O4]   < 0)
        Input.atm.interpol_method_gas[MOL_O4]   = INTERP_METHOD_LINMIX;  
    if (Input.atm.interpol_method_gas[MOL_SO2]  < 0)
        Input.atm.interpol_method_gas[MOL_SO2]  = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_CH4]  < 0)
        Input.atm.interpol_method_gas[MOL_CH4]  = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_N2O]  < 0)
        Input.atm.interpol_method_gas[MOL_N2O]  = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_CO]   < 0)
        Input.atm.interpol_method_gas[MOL_CO]   = INTERP_METHOD_LINMIX; 
    if (Input.atm.interpol_method_gas[MOL_N2]   < 0)
        Input.atm.interpol_method_gas[MOL_N2]   = INTERP_METHOD_LINMIX; 
    break;
  case Z_INTERPOLATE_SPLINE:
    /*        CONTINUOUS FIRST DERIVATIVES, BUT COLUMNS AND OPTICAL PROPERTIES NOT CONSISTENT              */
    /* SPLINE means spline interpolation for additional levels zout and altitude levels as specified below */
    /* columns summed up  with LINMIX assumption (this might be improved, if nessesary)                    */   
    /* optical properties with LINMIX assumption (this might be improved, if nessesary)                    */
    /* in general: INTERP_METHOD_LINMIX not possible for press, temper, refind, rh, MOL_AIR                */

    if (Input.atm.interpol_method_press  < 0) 
        Input.atm.interpol_method_press  = INTERP_METHOD_LOG_SPLINE; /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_temper < 0)
        Input.atm.interpol_method_temper = INTERP_METHOD_SPLINE; 
    if (Input.atm.interpol_method_refind < 0)
        Input.atm.interpol_method_refind = INTERP_METHOD_SPLINE;

    if (Input.atm.interpol_method_gas[MOL_AIR]  < 0)
        Input.atm.interpol_method_gas[MOL_AIR]  = INTERP_METHOD_LOG_SPLINE;    /* do not make changes here, it will not affect all important parts*/
    if (Input.atm.interpol_method_gas[MOL_O3]   < 0)
        Input.atm.interpol_method_gas[MOL_O3]   = INTERP_METHOD_SPLINE; 
    if (Input.atm.interpol_method_gas[MOL_O2]   < 0)
        Input.atm.interpol_method_gas[MOL_O2]   = INTERP_METHOD_LOG_SPLINE; 
    if (Input.atm.interpol_method_gas[MOL_H2O]  < 0)
        Input.atm.interpol_method_gas[MOL_H2O]  = INTERP_METHOD_LOG_SPLINE;    
    if (Input.atm.interpol_method_gas[MOL_CO2]  < 0)
        Input.atm.interpol_method_gas[MOL_CO2]  = INTERP_METHOD_LOG_SPLINE;    
    if (Input.atm.interpol_method_gas[MOL_NO2]  < 0)
        Input.atm.interpol_method_gas[MOL_NO2]  = INTERP_METHOD_SPLINE;
    if (Input.atm.interpol_method_gas[MOL_BRO]  < 0)
        Input.atm.interpol_method_gas[MOL_BRO]  = INTERP_METHOD_LOG_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_OCLO] < 0)
        Input.atm.interpol_method_gas[MOL_OCLO] = INTERP_METHOD_LOG_SPLINE;    
    if (Input.atm.interpol_method_gas[MOL_HCHO] < 0)
        Input.atm.interpol_method_gas[MOL_HCHO] = INTERP_METHOD_LOG_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_O4]   < 0)
        Input.atm.interpol_method_gas[MOL_O4]   = INTERP_METHOD_LOG_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_SO2]  < 0)
        Input.atm.interpol_method_gas[MOL_SO2]  = INTERP_METHOD_SPLINE;
    if (Input.atm.interpol_method_gas[MOL_CH4]  < 0)
        Input.atm.interpol_method_gas[MOL_CH4]  = INTERP_METHOD_LOG_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_N2O]  < 0)
        Input.atm.interpol_method_gas[MOL_N2O]  = INTERP_METHOD_LOG_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_CO]   < 0)
        Input.atm.interpol_method_gas[MOL_CO]   = INTERP_METHOD_SPLINE;  
    if (Input.atm.interpol_method_gas[MOL_N2]   < 0)
        Input.atm.interpol_method_gas[MOL_N2]   = INTERP_METHOD_LOG_SPLINE;  
    break;
  default:
    fprintf (stderr, "Error, unknown interpolation method z_interpolate %d!\n", Input.atm.z_interpolate);
    fprintf (stderr, "This is a program bug, please contact the programmers!\n");
    return -1;
    break;
  }

  if (Input.atm.interpol_method_temper != INTERP_METHOD_LINEAR &&
      Input.atm.interpol_method_temper != INTERP_METHOD_LOG && 
      Input.atm.interpol_method_temper != INTERP_METHOD_SPLINE) {
    fprintf (stderr, "Error, interpolation method %d not possible\n", Input.atm.interpol_method_temper);
    fprintf (stderr, "for vertical interpolation of temperature.\n");
    return -1;
  }

  /* no user-defined transmittance grid if correlated-k */
  if (Input.ck_scheme!=CK_NONE && Input.ck_scheme!=CK_LOWTRAN && Input.ck_scheme!=CK_REPTRAN && strlen(Input.filename[FN_WLTRANS])>0) {
    strcpy (Input.filename[FN_WLTRANS], "");
    fprintf (stderr, " ... option wavelength_grid_file ignored\n");
    fprintf (stderr, " ... (does not make sense with mol_abs_param)\n");
  }

  if (Input.alb.source == NOT_DEFINED_INTEGER) {
    /* default setting */
    Input.alb.source = ALBEDO_CONSTANT;
  }

  /* set default SURFACE TYPE map */
  if (strncasecmp("IGBP",Input.filename[FN_SURFACE_TYPE_MAP],4)==0) {
    strcpy (Input.filename[FN_SURFACE_TYPE_MAP], Input.filename[FN_PATH]);
    strcat (Input.filename[FN_SURFACE_TYPE_MAP], "albedo/IGBP_map/SURFACE_TYPE_IGBP_10min.cdf");

    if (Input.alb.source != ALBEDO_IGBP_LIBRARY && !Input.quiet)
      fprintf (stderr, "*** Warning, using IGBP map, but not 'albedo_library IGBP'\n");
  }

  /* message to the user because of name change */
  if (strncasecmp("IGBP",Input.filename[FN_ALBEDO_MAP],4)==0) {
    fprintf (stderr, "\n*** ERROR, starting from libRadtran version 1.2 \n");
    fprintf (stderr,   "           'albedo_map' contains broadband albedo values \n");
    fprintf (stderr,   "           and not surface_type information any more. \n");
    fprintf (stderr,   "           Sorry for the inconvenience! \n");
    fprintf (stderr,   "           please use instead: \n");
    fprintf (stderr,   "           surface_type_map IGBP \n");
    return -1;
  }

  /* set default SURFACE TYPE netCDF variable name */
  if ( strlen(Input.alb.netCDF_surf_name) == 0 )
    strcpy (Input.alb.netCDF_surf_name, "surface_type");

  /* set default RPV Polder library path  */
  if (Input.rpv.source == RPV_IGBP_LIBRARY) {

    strcpy (Input.filename[FN_RPV_LIB_PATH], Input.filename[FN_PATH]);
    strcat (Input.filename[FN_RPV_LIB_PATH], "albedo/IGBP_rpv_polder/");

  }

  /* set default ALTITUDE map */
  if (strncasecmp("GTOPO",Input.filename[FN_ALTITUDE_MAP],5)==0) {
    strcpy (Input.filename[FN_ALTITUDE_MAP], Input.filename[FN_PATH]);
    strcat (Input.filename[FN_ALTITUDE_MAP], "altitude/ELEVATION_GTOPO_10min.cdf");
  }


  /* check format of the ECMWF file */
  if (strlen(Input.filename[FN_ECMWF])>0) {

#if HAVE_LIBNETCDF

    /* open netcdf file */
    status = nc_open (Input.filename[FN_ECMWF], NC_NOWRITE, &ncid);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d opening netCDF file %s\n", status, Input.filename[FN_ECMWF]);
      return status;
    }
    status = nc_inq_varid (ncid, "z", &id_var_test);
    if (status==NC_NOERR) {
      Input.atm.ECMWF_new_format=TRUE;
    }
    else {
      status = nc_inq_varid (ncid, "Z", &id_var_test);
      if (status==NC_NOERR) {
        Input.atm.ECMWF_new_format=FALSE;
      }
      else {
        fprintf (stderr, "Error, unknown format of the ECMWF atmosphere file %s\n", Input.filename[FN_ECMWF]);
        fprintf (stderr, "       neither 'z' nor 'Z' was found in data file\n");
        return status;
      }
    }
    nc_close (ncid);

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use the ECMWF_atmosphere_file option. Please get netcdf and rebuild.*\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

  }
  

  /* set default ALTITUDE netCDF variable name */
  if ( strlen(Input.alt.netCDF_alt_name) == 0 ) {
    if ( Input.atm.ECMWF_new_format == TRUE && strlen(Input.filename[FN_ECMWF])>0 )
      strcpy (Input.alt.netCDF_alt_name, "z");
    else
      strcpy (Input.alt.netCDF_alt_name, "Z");
  }

  /* set default surface temperature netCDF variable name */
  if ( strlen(Input.filename[FN_SURFACE_TEMP_MAP]) > 0 ) {
#if HAVE_LIBNETCDF

    status = nc_open (Input.filename[FN_SURFACE_TEMP_MAP], NC_NOWRITE, &ncid);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d opening netCDF file %s\n", status, Input.filename[FN_SURFACE_TEMP_MAP]);
      return status;
    }
    
    /* if no name is given in the input file, try 'skt' and 'SKT' */
    if ( strlen(Input.netCDF_name_surf_T) == 0 ) {
      if      ( (status = nc_inq_varid (ncid, "skt", &id_var_test)) == NC_NOERR )
        strcpy (Input.netCDF_name_surf_T, "skt");
      else if ( (status = nc_inq_varid (ncid, "SKT", &id_var_test)) == NC_NOERR )
        strcpy (Input.netCDF_name_surf_T, "SKT");
      else {
        fprintf (stderr, "Error '%s' while getting id for surface temperature from '%s' \n", 
                                nc_strerror(status), Input.filename[FN_SURFACE_TEMP_MAP]);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return -1;
      }
    }
    nc_close (ncid);

#else
    fprintf (stderr, " ***************************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot                *\n");
    fprintf (stderr, " * use the 'surface_temperature_map' option. Please get netcdf and rebuild.*\n");
    fprintf (stderr, " ***************************************************************************\n");
    return -1;
#endif
  }

  int caoth_error=caoth_set(&Input,&Output);
  if (caoth_error != 0) {
	fprintf (stderr, "%d error(s) in caoth settings.\nPlease check uvspec input-file again.\n", abs(caoth_error));
	ierror++;
  }

  /* Check if modify for aerosol. ic, wc, and profile are set correctly */

  if ( Input.aer.modify[MODIFY_VAR_GG][MODIFY_TYPE_SCALE] != NOT_DEFINED_FLOAT || Input.aer.modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SCALE] != NOT_DEFINED_FLOAT ) {
	fprintf ( stderr, "Error 'aerosol_modify gg scale' and 'aerosol_modify tau550 scale' are not yet implemented.\n" ) ;
	ierror++;
  }
  for ( isp = 0; isp < Input.n_caoth; isp ++ ) {
	if ( Input.caoth[isp].modify[MODIFY_VAR_TAU550][MODIFY_TYPE_SCALE] != NOT_DEFINED_FLOAT ) {
		fprintf ( stderr, "Error '%s_modify tau550 scale' and 'profile_modify %s tau550 scale' are not yet implemented.\n", Input.caoth[isp].name, Input.caoth[isp].name ) ;
		ierror++;
	}
/*	if ( Input.caoth[isp].modify[MODIFY_VAR_TAU][MODIFY_TYPE_SCALE] != NOT_DEFINED_FLOAT ) {
	   if ( strcasecmp ( Input.caoth[isp].name, "wc" ) != 0 || strcasecmp ( Input.caoth[isp].name, "ic" ) != 0 ) {
		fprintf ( stderr, "Error 'profile_modify %s tau scale' is not yet implemented.\n", Input.caoth[isp].name ) ;
		ierror++;
	   }
	}*/
  }

  /* Check if no_absorption and no_scattering is set correctly */
  for (ispoff = 0; ispoff < Input.n_caothoff; ispoff ++ ) {
	if  (Input.caothoff[ispoff].no_absorption == 1) {
		if (strcasecmp ( Input.caothoff[ispoff].name, "all") == 0 ) { 
			Input.absorption = 0;
			Input.molabs = MOLABS_NONE;
		}
		else if (strcasecmp ( Input.caothoff[ispoff].name, "mol") == 0 ) Input.molabs = MOLABS_NONE;
		else {
        		fprintf (stderr, "Error 'no_absorption %s' not yet implemented.' \n", Input.caothoff[ispoff].name );
			ierror++;
		}
	}
	if (Input.caothoff[ispoff].no_scattering == 1) {
		if (strcasecmp ( Input.caothoff[ispoff].name, "all") == 0 ) {
		        Input.scattering = 0;
		        Input.rayleigh = RAYLEIGH_NONE;
                }
		else if (strcasecmp ( Input.caothoff[ispoff].name, "mol") == 0 ) Input.rayleigh = RAYLEIGH_NONE;
		else if (strcasecmp ( Input.caothoff[ispoff].name, "aer") == 0 ) Input.aer.no_scattering = TRUE;
		else {
			int found_caoth=0;
			for (isp=0; isp<Input.n_caoth; isp++) {
				if ( strcasecmp ( Input.caothoff[ispoff].name, Input.caoth[isp].name) == 0 ) {
					Input.caoth[isp].no_scattering = TRUE;
					found_caoth=1;
					break;
				}
			}
			if (!found_caoth) {
		        	fprintf (stderr, "Error 'no_scattering %s' not yet implemented.' \n", Input.caothoff[ispoff].name );
				ierror++;
			}
		}
	}
  }

  if (ierror!=0)  {
    fprintf (stderr, "Error in input file! Please check your input file and retry!\n");
    return -1;
  }

  /* Default MYSTIC output filename */
  if (strlen(Input.rte.mc.filename[FN_MC_BASENAME]) == 0)
    strcpy (Input.rte.mc.filename[FN_MC_BASENAME], "mc");
  
  /* big caoth loop */
  for (isp=0; isp<Input.n_caoth;isp++) {

    /* caoth properties */
    switch (Input.caoth[isp].properties) {
    case PROP_MIE:
      /* Mie water cloud properties */ 

      switch (Input.ck_scheme) {
      case CK_KATO:
      case CK_KATO2:
      case CK_KATO2_96:
      case CK_FU:
      case CK_AVHRR_KRATZ:

        if (Input.rte.disort_icm == DISORT_ICM_MOMENTS) {
          /* need to use old files, this will be soon be obsolete */
          if (Input.caoth[isp].interpolate)  {
            fprintf (stderr, "* It does not make sense to use \"wc_properties mie interpolate\" \n");
            fprintf (stderr, "* for \"mol_abs_param KATO/FU/AVHRR_KRATZ\" because in the latter\n");
            fprintf (stderr, "* case the cloud properties are pre-calculated for the required\n");
            fprintf (stderr, "* wavelength and interpolation is not useful. Please remove the\n");
            fprintf (stderr, "* optional argument 'interpolate' from the 'wc_properties' uvspec input.\n");
            fprintf (stderr, "* Note that this only happens if you are using disort_intcor moments.\n");
            fprintf (stderr, "* In such case, the old Mie files specifically create for the kato/fu/kratz bands\n");
  	  fprintf (stderr, "* is being used. In case of using disort_intcor phases (default), the\n");
  	  fprintf (stderr, "* new Mie files are being used, and interpolation is performed.\n");
            return -1;
          }
     
          strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
	  switch (Input.ck_scheme) {
	      case CK_KATO:
	      case CK_KATO2:
	      case CK_KATO2_96:
		strcat (Input.caoth[isp].properties_filename, "/correlated_k/kato/");
		if (Input.verbose)	fprintf(stderr, "... using old mie files specifically created for kato bands\n");
		break;
	      case CK_FU:
		strcat (Input.caoth[isp].properties_filename, "/correlated_k/fu/");
		if (Input.verbose)	fprintf(stderr, "... using old mie files specifically created for Fu bands\n");
		break;
	      case CK_AVHRR_KRATZ:
		strcat (Input.caoth[isp].properties_filename, "/correlated_k/kratz/");
          	if (Input.verbose)	fprintf(stderr, "... using old mie files specifically created for kratz bands\n");
		break;
	  }
          strcat (Input.caoth[isp].properties_filename, "wc.");
        }
        else {
          strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
          strcat (Input.caoth[isp].properties_filename, "/wc/mie/");
          strcat (Input.caoth[isp].properties_filename, "wc.");
          Input.caoth[isp].interpolate=1; 
          if (Input.verbose) {
	    switch (Input.ck_scheme) {
	      case CK_KATO:
              case CK_KATO2:
              case CK_KATO2_96:
		fprintf(stderr, "... interpolate optical properties to central wavelengths of kato bands\n");
		break;
	      case CK_FU:
		fprintf(stderr, "... interpolate optical properties to central wavelengths of fu bands\n");
		break;
	      case CK_AVHRR_KRATZ:
		fprintf(stderr, "... interpolate optical properties to central wavelengths of kratz bands\n");
		break;
	    }
	  }
	}
        break;	    
  
      case CK_NONE:
      case CK_RAMAN:
      case CK_LOWTRAN:
      case CK_REPTRAN:
      case CK_REPTRAN_CHANNEL:
        if (!Input.caoth[isp].interpolate)  {
          fprintf (stderr, "* 'wc_properties mie' for a spectral calculation requires that the\n");
          fprintf (stderr, "* cloud optical properties are interpolated to the internal wavelength\n");
  	fprintf (stderr, "* grid. Depending on the number of grid points this may require huge\n");
  	fprintf (stderr, "* amounts of memory and some additional computational time. To use this\n");
  	fprintf (stderr, "* option, please specify explicitely 'wc_properties mie interpolate'\n");
  	fprintf (stderr, "* in the uvspec input!\n");
  	return -1;
        }
        else {
          strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
          strcat (Input.caoth[isp].properties_filename, "/wc/mie/");
          strcat (Input.caoth[isp].properties_filename, "wc.");
        }
        break;
  
      case CK_FILE:
        fprintf (stderr, "Error, 'wc_properties mie' not implemented for a generic\n");
        fprintf (stderr, "mol_abs_param calculation\n");
        return -1;
  
      default:
        fprintf (stderr, "Error, unsupported mol_abs_param scheme %d\n", Input.ck_scheme);
        return -1;
      }
      break;

    case PROP_IC_MIE:
      /* Mie ice cloud properties */ 

      switch (Input.ck_scheme) {
      case CK_KATO:
      case CK_KATO2:
      case CK_KATO2_96:
      case CK_FU:	 
      case CK_AVHRR_KRATZ:
        if (Input.rte.disort_icm == DISORT_ICM_PHASE) {
          fprintf(stderr, "* ic_properties mie for mol_abs_param kato/fu/avhrr_kratz are not implemented for the disort_intcor phase.\n");
          fprintf(stderr, "* Please use disort_intcor moments. \n");
          return -1;
        }
        if (Input.caoth[isp].interpolate)  {
          fprintf (stderr, "* It does not make sense to use \"ic_properties mie interpolate\"\n");
          fprintf (stderr, "* for \"mol_abs_param KATO/FU/AVHRR_KRATZ\" because\n");
  	  fprintf (stderr, "* in the latter case the cloud properties are pre-calculated\n");
          fprintf (stderr, "* for the required wavelength and interpolation is not useful.\n");
  	  fprintf (stderr, "* Please remove 'interpolate' from the 'ic_properties' uvspec input.\n");
  	  return -1;
        }
  
        strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
	switch (Input.ck_scheme) {
	  case CK_KATO:
	  case CK_KATO2:
	  case CK_KATO2_96:
		strcat (Input.caoth[isp].properties_filename, "/correlated_k/kato/");
		break;
          case CK_FU:
        	strcat (Input.caoth[isp].properties_filename, "/correlated_k/fu/");
		break;
	  case CK_AVHRR_KRATZ:
        	strcat (Input.caoth[isp].properties_filename, "/correlated_k/kratz/");
		break;
	}
        strcat (Input.caoth[isp].properties_filename, "ic.");
        break;	    
  		    
      case CK_NONE:
      case CK_RAMAN:
      case CK_LOWTRAN:
      case CK_REPTRAN:
      case CK_REPTRAN_CHANNEL:
        if (!Input.caoth[isp].interpolate)  {
          fprintf (stderr, "* 'ic_properties mie' for a spectral calculation requires that the\n");
          fprintf (stderr, "* cloud optical properties are interpolated to the internal wavelength\n");
	  fprintf (stderr, "* grid. Depending on the number of grid points this may require huge\n");
  	  fprintf (stderr, "* amounts of memory and some additional computational time. To use this\n");
  	  fprintf (stderr, "* option, please specify explicitely 'ic_properties mie interpolate'\n");
  	  fprintf (stderr, "* in the uvspec input!\n");
  	  return -1;
        }
        else {
          strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
          strcat (Input.caoth[isp].properties_filename, "/ic/mie/");
          strcat (Input.caoth[isp].properties_filename, "ic.");
        }
        break;
  
      case CK_FILE:
        fprintf (stderr, "Error, 'ic_properties mie' not implemented for a generic\n");
        fprintf (stderr, "correlated-k calculation\n");
        return -1;
  
      default:
        fprintf (stderr, "Error, unsupported mol_abs_param scheme %d\n", Input.ck_scheme);
        return -1;
      }
  
      if (!Input.quiet)
         fprintf (stderr, " ... reading ice cloud optical properties from %s*\n", 
  	        Input.caoth[isp].properties_filename);
      break;
  
    case PROP_BAUM:
    /* Baum ice cloud properties */ 
  
      if (!Input.caoth[isp].interpolate)  {
        fprintf (stderr, "* 'ic_properties baum' for a spectral calculation requires that the\n");
        fprintf (stderr, "* cloud optical properties are interpolated to the internal wavelength\n");
        fprintf (stderr, "* grid. Depending on the number of grid points this may require huge\n");
        fprintf (stderr, "* amounts of memory and some additional computational time. To use this\n");
        fprintf (stderr, "* option, please specify explicitely 'ic_properties baum interpolate'\n");
        fprintf (stderr, "* in the uvspec input!\n");
        
        return -1;
      }
      else if(Input.wl.start < 3076.92 &&  Input.wl.end > 2200){
        fprintf(stderr,"* 'ic_properties baum' is not available in the\n");
        fprintf(stderr,"* wavelength range from 2200 to 3076.92 nm. \n");
        return -1;
        }
      else {
        strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
  
        if (Input.caoth[isp].properties == PROP_BAUM)
          strcat (Input.caoth[isp].properties_filename, "ic/baum/");
    
        strcat (Input.caoth[isp].properties_filename, "ic.");
      }
  
      if (!Input.quiet)
         fprintf (stderr, " ... reading ice cloud optical properties from %s*\n", 
  	        Input.caoth[isp].properties_filename);

      break;
      case PROP_BAUM_V36:
      /* Baum V3.6 ice cloud properties */ 
      if (!Input.caoth[isp].interpolate)  {
        fprintf (stderr, "* 'ic_properties baum_v36' for a spectral calculation requires that the\n");
        fprintf (stderr, "* cloud optical properties are interpolated to the internal wavelength\n");
        fprintf (stderr, "* grid. Depending on the number of grid points this may require huge\n");
        fprintf (stderr, "* amounts of memory and some additional computational time. To use this\n");
        fprintf (stderr, "* option, please specify explicitely 'ic_properties baum interpolate'\n");
        fprintf (stderr, "* in the uvspec input!\n");
        
        return -1;
      }
      else if(Input.wl.start < 200 &&  Input.wl.end > 100000){
        fprintf(stderr,"* 'ic_properties baum_v36' is available in the\n");
        fprintf(stderr,"* wavelength range from 200 nm to 100 microns. \n");
        return -1;
      }
      else {
        strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
        strcat (Input.caoth[isp].properties_filename, "ic/baum_v3.6/");
        strcat (Input.caoth[isp].properties_filename, "ic.");
      }
  
      switch(Input.caoth[isp].habit){
      case IC_HABIT_GHM:
        strcat(Input.caoth[isp].properties_filename, "ghm");
        break;
      case IC_HABIT_SOLID_COLUMN:
        strcat(Input.caoth[isp].properties_filename, "col");
        break;
      case IC_HABIT_ROUGH_AGGREGATE:
        strcat(Input.caoth[isp].properties_filename, "agg");
        break;
       case IC_HABIT_ROSETTE_6:
       case IC_HABIT_PLATE:
       case IC_HABIT_DROXTAL:
       case IC_HABIT_ROSETTE_4:
        fprintf (stderr, "Error, ic_habit is not implemented for \n"); 
        fprintf (stderr, "ic_properties baum_v36.\n");
        return -1;
        break; 
      default: 
          fprintf (stderr, "Error, unknown ice particle habit. \n"); 
          return -1; 
      }
  
      if (!Input.quiet)
        fprintf (stderr, " ... reading ice cloud optical properties from %s*\n", 
			 Input.caoth[isp].properties_filename);
      break;
      case PROP_HEY:
      /* HEY ice cloud properties */ 
  
      if (!Input.caoth[isp].interpolate)  {
        fprintf (stderr, "* 'ic_properties hey' for a spectral calculation requires that the\n");
        fprintf (stderr, "* cloud optical properties are interpolated to the internal wavelength\n");
        fprintf (stderr, "* grid. Depending on the number of grid points this may require huge\n");
        fprintf (stderr, "* amounts of memory and some additional computational time. To use this\n");
        fprintf (stderr, "* option, please specify explicitely 'ic_properties hey interpolate'\n");
        fprintf (stderr, "* in the uvspec input!\n");
        
        return -1;
      }
      else if(Input.wl.start < 200 &&  Input.wl.end > 5000){
        fprintf(stderr,"* 'ic_properties hey' is only available in the\n");
        fprintf(stderr,"* wavelength range from 200 to 5000 nm. \n");
        return -1;
      }
      else {
        strcpy (Input.caoth[isp].properties_filename, Input.filename[FN_PATH]);
        strcat (Input.caoth[isp].properties_filename, "ic/hey/");
        strcat (Input.caoth[isp].properties_filename, "ic.");
      }
  
      switch(Input.caoth[isp].habit){
      case IC_HABIT_SOLID_COLUMN:
        strcat(Input.caoth[isp].properties_filename, "col");
        break;
      case IC_HABIT_HOLLOW_COLUMN:
        strcat(Input.caoth[isp].properties_filename, "hol");
        break;
      case IC_HABIT_ROUGH_AGGREGATE:
        strcat(Input.caoth[isp].properties_filename, "agg");
        break;
      case IC_HABIT_ROSETTE_6:
        strcat(Input.caoth[isp].properties_filename, "ros");
        break;	
      case IC_HABIT_PLATE:
        strcat(Input.caoth[isp].properties_filename, "pla");
        break;
      case IC_HABIT_DROXTAL:
        strcat(Input.caoth[isp].properties_filename, "dro");
        break;
      case IC_HABIT_GHM:
        strcat(Input.caoth[isp].properties_filename, "ghm");
        break;
      case IC_HABIT_ROSETTE_4:
        fprintf (stderr, "Error, ic_habit rosette-4 is not implemented for \n"); 
        fprintf (stderr, "ic_properties hey.\n");
        return -1;
        break; 
      default: 
          fprintf (stderr, "Error, unknown ice particle habit. \n"); 
          return -1; 
      }
  
      if (!Input.quiet)
        fprintf (stderr, " ... reading ice cloud optical properties from %s*\n", 
			 Input.caoth[isp].properties_filename);
      break;
      
   default:
   break; 
   }
 } /* end loop isp */

  /* aerosol properties */
  if (Input.aer.standard == 1) {
    if (Input.aer.seasn==0)
      Input.aer.seasn=1;
    if (Input.aer.vulcan==0)
      Input.aer.vulcan=1;
    if (Input.aer.haze==0)
      Input.aer.haze=1;
    if (Input.aer.visibility==NOT_DEFINED_FLOAT)
      Input.aer.visibility=50.0;
  }

  if (Input.aer.seasn>0 && Input.aer.vulcan>0 && Input.aer.haze>0 && Input.aer.visibility>0)
    Input.aer.standard=1;


  /* BRDR stuff */
  Input.disort2_brdf = BRDF_NONE;
  brdftype=0;

  /******************************************/
  /* First, RPV (Rahman, Pinty, Verstraete) */
  /******************************************/

  /* homogeneous RPV */
  /* homogeneous, wavelength-dependent RPV */
  if (Input.rpv.rpv[BRDF_RPV_K]!=0 || Input.rpv.rpv[BRDF_RPV_RHO0]!=0 || Input.rpv.rpv[BRDF_RPV_THETA]!=0 || strlen(Input.filename[FN_RPV])>0) {
    brdftype++;
    Input.disort2_brdf = BRDF_RPV;

    if (!Input.quiet)
      fprintf (stderr, " ... using RPV parameterization of the BRDF (disort, fdisort2 and MYSTIC only!)\n");
  }


  /* 2D RPV */
  if (strlen(Input.rte.mc.filename[FN_MC_RPV])>0)  {
    brdftype++;

    if (!Input.quiet)
      fprintf (stderr, " ... using 2D RPV parameterization of the BRDF (MYSTIC only!)\n");
  }


  /* check if both rpv_file and rpv_type were defined */
  if ((strlen(Input.rte.mc.filename[FN_MC_RPV])>0  && strlen(Input.rte.mc.filename[FN_MC_RPVTYPE])==0) ||
      (strlen(Input.rte.mc.filename[FN_MC_RPV])==0 && strlen(Input.rte.mc.filename[FN_MC_RPVTYPE])>0)) {
    fprintf (stderr, "Error, you need to define both mc_rpv_file and mc_rpv_type if you want 2D RPV\n");
    status--;
  }

  
  /******************************************/
  /* Second, Cox and Munk                   */
  /******************************************/

  if (Input.cm.param[BRDF_CAM_U10] >= 0 || Input.cm.param[BRDF_CAM_PCL] >= 0 || Input.cm.param[BRDF_CAM_SAL] >= 0 || Input.disort2_brdf == BRDF_CAM) {
    brdftype++;

    Input.disort2_brdf = BRDF_CAM;

    if (!Input.quiet)
      fprintf (stderr, " ... using Cox and Munk parameterization of the BRDF (disort, fdisort2 and MYSTIC only!)\n");
  }

  /* set default values for windspeed, pigment concentration, and salinity */
  if (Input.disort2_brdf == BRDF_CAM)  {

    if (Input.cm.param[BRDF_CAM_U10]<0)  {
      if (strlen(Input.filename[FN_U10_MAP])==0) {
        fprintf (stderr, "Error setting up Cox and Munk ocean BRDF: At least a wind speed brdf_cam u10\n");
	fprintf (stderr, "is required. Please define or switch off ocean BRDF!\n");
	return -1;
      }
      else {
        /* read wind speed from netCDF wind map */
        status = read_u10_from_map (Input.latitude, Input.longitude, Input.UTC, Input.atm.time_interpolate, 
                                    Input.filename[FN_U10_MAP],
                                    &Input.cm.param[BRDF_CAM_U10], Input.verbose, Input.quiet);
        if (status!=0) {
          fprintf (stderr, "Error %d during read_u10_from_map in uvspec_lex.l \n", status);
          return status;
        }
      }
      if (Input.verbose)
        fprintf (stderr, " ... setting wind speed for Cox and Munk to %f m/s\n", Input.cm.param[BRDF_CAM_U10]);
    }

    if (Input.cm.param[BRDF_CAM_PCL]<0)  {

      if (strlen(Input.filename[FN_SALINITY_MAP])==0) {
        Input.cm.param[BRDF_CAM_PCL] = 0.01;
      }
      else {

        if ( strlen(Input.cm.pcl_netCDF_name) == 0 )
          strcpy (Input.cm.pcl_netCDF_name, "chlorophyll");

        /* read ocean water pigment concentration from map */
        status = get_number_from_netCDF_map (Input.latitude, Input.longitude, Input.UTC, Input.atm.time_interpolate,
                                             Input.filename[FN_PIGMENTS_MAP], &(Input.cm.param[BRDF_CAM_PCL]), TYPE_FLOAT, Input.cm.pcl_netCDF_name,
                                             Input.verbose, Input.quiet);
        if (status!=0) {
          if ( status == ERROR_READ_MISSING_VALUE )
            fprintf (stderr, "     Error '%s' reading '%s' from %s in %s (%s) \n", 
                                "read missing value", "chlorophyll", Input.filename[FN_PIGMENTS_MAP], function_name, file_name ); 
          else { 
#if HAVE_LIBNETCDF
            fprintf (stderr, "     Error '%s' reading '%s' from %s in %s (%s) \n", 
                                nc_strerror(status), "chlorophyll", Input.filename[FN_PIGMENTS_MAP], function_name, file_name );
#else
            fprintf (stderr, "     Error %d reading '%s' from %s in %s (%s) \n", 
                                            status , "chlorophyll", Input.filename[FN_PIGMENTS_MAP], function_name, file_name );
#endif
          }
          return status;
        }
      }

      /* take care of user scale value */
      Input.cm.param[BRDF_CAM_PCL] *= Input.cm.pcl_scale_factor;

      /* automatic check for plausible pigment concentration */
      if ( (Input.cm.param[BRDF_CAM_PCL] < 0.0 || 100.0 < Input.cm.param[BRDF_CAM_PCL]) && !Input.quiet )
        fprintf (stderr, "*** Warning: ocean pigment concentration is not plausible: %f mg/m3\n", Input.cm.param[BRDF_CAM_PCL]  );
      if ( Input.cm.param[BRDF_CAM_PCL] < 0.0 ) 
        return -1;

      if (Input.verbose)
        fprintf (stderr, " ... setting pigment concentration for Cox and Munk to %f mg/m3\n", Input.cm.param[BRDF_CAM_PCL]);
    }

    if (Input.cm.param[BRDF_CAM_SAL]<0)  {

      if (strlen(Input.filename[FN_SALINITY_MAP])==0) {
        Input.cm.param[BRDF_CAM_SAL] = 34.3;
      }
      else {

        if ( strlen(Input.cm.sal_netCDF_name) == 0 )
          strcpy (Input.cm.sal_netCDF_name, "salinity");

        /* read salinity from map */
        status = get_number_from_netCDF_map (Input.latitude, Input.longitude, Input.UTC, Input.atm.time_interpolate,
                                             Input.filename[FN_SALINITY_MAP], &(Input.cm.param[BRDF_CAM_SAL]), TYPE_FLOAT, Input.cm.sal_netCDF_name,
                                             Input.verbose, Input.quiet);
        if ( status == ERROR_READ_MISSING_VALUE ) {
          Input.cm.param[BRDF_CAM_SAL] = 34.3;
          if ( !Input.quiet ) {
            fprintf (stderr, "*** Warning: read missing value from map. \n" );
            fprintf (stderr, "***          set salinity to default value of %f \n", Input.cm.param[BRDF_CAM_SAL] );
          }
          status=0;
        }

        if (status!=0) {
#if HAVE_LIBNETCDF
          fprintf (stderr, "     Error '%s' reading '%s' from %s in %s (%s) \n", 
                                nc_strerror(status), "SALT", Input.filename[FN_SALINITY_MAP], function_name, file_name );
#else
          fprintf (stderr, "     Error %d reading '%s' from %s in %s (%s) \n", 
                                            status , "SALT", Input.filename[FN_SALINITY_MAP], function_name, file_name );
#endif
          return status;
        }
      }

      /* take care of user scale value */
      Input.cm.param[BRDF_CAM_SAL] *= Input.cm.sal_scale_factor;

      /* automatic check for plausible salinity range */
      if ( (Input.cm.param[BRDF_CAM_SAL] < 0.0 || 55.0 < Input.cm.param[BRDF_CAM_SAL]) && !Input.quiet )
        fprintf (stderr, "*** Warning: salinity value is not plausible: %f ppt\n", Input.cm.param[BRDF_CAM_SAL]  );

      if ( Input.cm.param[BRDF_CAM_SAL] < 0.0 ) 
        return -1;

      if (Input.verbose)
        fprintf (stderr, " ... setting salinity for Cox and Munk to %f ppm\n", Input.cm.param[BRDF_CAM_SAL]);
    }
    
    /* we don't allow wind speeds smaller than 1 m/s because the  */
    /* BRDF is expressed as a series of nstr Legendre polynomials */
    /* which causes unwanted structures not only in the BRDF but  */
    /* also in the albedo                                         */
    /*     ??? Is this really required for MYSTIC - check ???     */
    if (Input.cm.param[BRDF_CAM_U10]<1)  {
      Input.cm.param[BRDF_CAM_U10] = 1.0;
      fprintf (stderr, " ... increasing wind speed for Cox and Munk to %.2f m/s\n", Input.cm.param[BRDF_CAM_U10]);
    }
    
    /* we also don't want nstr < 16 for the same reason, even for */
    /* the calculation of fluxes                                  */

    if (Input.rte.solver!=SOLVER_MONTECARLO && Input.rte.nstr<16)  {
      fprintf (stderr, " ... increasing number_of_streams to 16 because of Cox and Munk!\n");
      Input.rte.nstr=16;
    }
  }


  /******************************************/
  /* Third, AMBRALS                         */
  /******************************************/

  if (Input.ambrals.ambrals[BRDF_AMBRALS_ISO] >= 0 || Input.ambrals.ambrals[BRDF_AMBRALS_VOL] >= 0 || Input.ambrals.ambrals[BRDF_AMBRALS_GEO] >= 0) {
    brdftype++;

    Input.disort2_brdf = BRDF_AMB;
    
    fprintf (stderr, " ... AMBRALS BRDF (disort, fdisort2 and MYSTIC only!)\n");
    fprintf (stderr, " ... iso = %f, vol = %f, geo = %f\n", Input.ambrals.ambrals[BRDF_AMBRALS_ISO], Input.ambrals.ambrals[BRDF_AMBRALS_VOL], Input.ambrals.ambrals[BRDF_AMBRALS_GEO]);
  }

  
  /* Check if only one or more BRDF types have been defined */
  /* but allow 2D RPV in combination with Cox and Munk      */
  if (brdftype > 1 && !(Input.disort2_brdf==BRDF_CAM && strlen(Input.rte.mc.filename[FN_MC_RPV])))  {
    fprintf (stderr, "Error, it does not make sense to define more than one BRDF type!\n");
    return -1;
  }




  /**********************************************/
  /* Tsang BPDF                                 */
  /**********************************************/
  
  if (Input.bpdf.u10 >=0){
     Input.bpdf.type=BPDF_TSANG;
     
     if (!Input.quiet) 
        fprintf(stderr, " ... TSANG BPDF - newly implemented, please check results carefully !!!! \n");

     if(Input.rte.solver!=SOLVER_MONTECARLO){
       fprintf(stderr, "Error: Polarized reflectance matrix can currently only be specified with\n");
       fprintf(stderr, "rte_solver montecarlo \n");
       return -1;
     }
  }   
       
  
  /* Number of phase function moments equals the number of streams  */
  /* per default, but may be changed later in the code;             */
  /* in case of disort2 we increase nmom to at least 1024 to assure */
  /* correct representation of the phase function with only little  */
  /* computational cost                                             */
  Output.atm.nmom = Input.rte.nstr;      

  switch (Input.rte.solver) {
  case SOLVER_FDISORT2:
  case SOLVER_DISORT:
  case SOLVER_SSS:
  case SOLVER_SSSI:
    Output.atm.nmom = (Input.rte.nstr>1024 ? Input.rte.nstr : 1024);
    break;
  case SOLVER_SDISORT:
    if (Input.rte.sdisort[SDISORT_NSCAT] != 1 && Input.rte.sdisort[SDISORT_NSCAT] != 2)
       Input.rte.sdisort[SDISORT_NSCAT] = 2;
    break;
  case SOLVER_FDISORT1:
  case SOLVER_FTWOSTR:
  case SOLVER_SOS:
  case SOLVER_MONTECARLO:
  case SOLVER_POLRADTRAN:
  case SOLVER_SPSDISORT:
  case SOLVER_TZS:
  case SOLVER_NULL:
  case SOLVER_RODENTS:
  case SOLVER_TWOSTREBE:
  case SOLVER_TWOSTR:
  case SOLVER_SSLIDAR:
    break;
  default:
    fprintf (stderr, "Error, unknown rte_solver %d (line %d, function '%s' in '%s')\n", Input.rte.solver, __LINE__, __func__, __FILE__);
    return -1;
  }

  Output.atm.molabs   = Input.molabs;
  Output.atm.rayleigh = Input.rayleigh;
  
  if (Input.rte.nphi <= 0 && Input.rte.ibcnd==0)  {
    switch (Input.rte.solver) {
    case SOLVER_SDISORT:
    case SOLVER_SPSDISORT:
    case SOLVER_FDISORT1:
    case SOLVER_FDISORT2:
    case SOLVER_DISORT:
      Input.rte.maxphi = 1;  /* Minimum number required by disort. */
      Input.rte.nphi   = 1;
      Input.rte.phi    = (float *) calloc (Input.rte.maxphi, sizeof(float));
      Output.print_phi = 0;  /* Do not print radiances if user does not want them */
      break;
    case SOLVER_POLRADTRAN:
      if (Input.rte.polradtran[POLRADTRAN_AZIORDER] == 0) 
	Input.rte.nphi = 0;
      break; 
    case SOLVER_FTWOSTR:
    case SOLVER_SOS:
    case SOLVER_MONTECARLO:
    case SOLVER_TZS:
    case SOLVER_SSS:
    case SOLVER_SSSI:
    case SOLVER_NULL:
    case SOLVER_RODENTS:
    case SOLVER_TWOSTREBE:
    case SOLVER_TWOSTR:
    case SOLVER_SSLIDAR:
      break;
    default:
      fprintf (stderr, "Error, unknown rte_solver %d (line %d, function '%s' in '%s')\n", Input.rte.solver, __LINE__, __func__, __FILE__);
      return -1;
    }
  }

/*aky  if ( Input.rte.solver == SOLVER_SOS ) {
        fprintf (stderr, "Successive order of scattering code not yet in production.\n");
        fprintf (stderr, "Please change to another solver.\n");
	return -1;
  }aky*/
  if ( Input.flu.source >= 0 ) {
    Input.calibration = OUTCAL_TRANSMITTANCE; /* For fluorescence fbeam is its real value and not 1.0 */
                                              /* Hence, extraterrestrial is included in the rte*/
                                              /* solution and should not be multiplied with    */
                                              /* afterwards                                    */
  }
  if ( Input.raman ) {
    
    if (Input.rte.solver != SOLVER_DISORT) {
      fprintf (stderr, "Only the 'rte_solver disort' is allowed with the 'raman' option.\n");
      return -1;
    }
    if (Input.ck_scheme != CK_RAMAN) {
      fprintf (stderr, "Only spectral calculations are allowed with the 'raman' option.\n");
      return -1;
    }

    if ( !Input.raman_original ) Input.raman_fast = 1;

    Input.n_raman_transitions_N2 =  48;           /* Same number as Spurr et al 2008 */
    Input.n_raman_transitions_O2 = 185;           /* Same number as Spurr et al 2008 */

    Output.crs.number_of_ramanshifts = Input.n_raman_transitions_N2 + Input.n_raman_transitions_O2;
    Output.crs.number_of_ramanwavelengths = Output.crs.number_of_ramanshifts + 1 ;

    Input.calibration = OUTCAL_TRANSMITTANCE; /* For Raman fbeam is its real value and not 1.0 */
                                              /* Hence, extraterrestrial is included in the rte*/
                                              /* solution and should not be multiplied with    */
                                              /* afterwards                                    */
    /* For Raman scattering only include wavelengths that the user asked. */
    /* Internally we have to include more wavelengths to account for      */
    /* Raman scattered radiation                                          */
    if ( Input.delta_wvl_raman <= 0 )  Output.wl.delta_wvl_raman = 8.0; 
    else                                  Output.wl.delta_wvl_raman = Input.delta_wvl_raman;
    
    /* Internally the lower and upper wavelengths for Raman scattering are different from */
    /* what the user specifies. This is setup in setup_wlgrid (ancillary.c).              */

    if ( Input.rte.numu > 0 ) {
      Input.rte.maxumu = Input.rte.nstr + Input.rte.numu;
    }
    else {
      Input.rte.maxumu = Input.rte.nstr;
    }

    Input.rte.cmuind = (int *) calloc (Input.rte.nstr, sizeof(int));
    Input.rte.umuind = (int *) calloc (Input.rte.numu, sizeof(int));
    Input.rte.umu    = add_cmu_to_umu(Input.rte.nstr, Input.rte.numu, Input.rte.maxumu, 
				      Input.rte.umu, Input.rte.cmuind, Input.rte.umuind);

    Input.rte.numu   = Input.rte.nstr + Input.rte.numu;
  }
  else if (Input.rte.numu <= 0) {
    Input.rte.maxumu = Input.rte.nstr;
    Input.rte.umu    = (float *) calloc (Input.rte.maxumu, sizeof(float));
    }
  
  /* Output processing: thermal source for brightness temperature */
  if (Input.calibration == OUTCAL_BRIGHTNESS && Input.source != SRC_THERMAL) {
    fprintf (stderr, "Error, brightness temperatures does only make sense for thermal sources!\n");
    return -1;
  }

  /* Output processing: thermal source for brightness temperature */
  if (Input.calibration == OUTCAL_BRIGHTNESS && ( Input.output_unit == UNIT_PER_NM || Input.output_unit == UNIT_PER_BAND ) ) {
    fprintf (stderr, "Error, unsupported combination of brightness temperatures with output unit!\n");
    return -1;
  }

  /* Output processing: no transmittance for thermal calculations */
  if (Input.calibration == OUTCAL_TRANSMITTANCE && Input.source == SRC_THERMAL) {
    fprintf (stderr, "Error, transmittance does not make sense for thermal sources!\n");
    return -1;
  }

  /* if thermal source, and no brightness temperature requested, use thermal calibration */
  if (Input.source == SRC_THERMAL && Input.calibration == OUTCAL_ABSOLUTE) {
    Input.calibration = OUTCAL_THERMAL;
  }


  if (Input.wl.start > Input.wl.end) {
    if (!Input.quiet) 
      fprintf(stderr," ... swap first and second wavelength.\n");
    tmpfloat = Input.wl.start;
    Input.wl.start = Input.wl.end;
    Input.wl.end = tmpfloat;
  }

  if (Input.wl.start <= 0 && Input.wl.start != NOT_DEFINED_INTEGER  ) {
    fprintf (stderr, "Error, start wavelength %f is negative.\n", Input.wl.start);
    return -1;
  }

  if (Input.wl.start_index > Input.wl.end_index) {
    if (!Input.quiet) 
      fprintf(stderr," ... swap first and second wavelength index.\n");
    tmpint = Input.wl.start_index;
    Input.wl.start_index = Input.wl.end_index;
    Input.wl.end_index = tmpint;
  }

  if (Input.wl.start_index <= 0 && Input.wl.start_index != -999) {
    fprintf (stderr, "Error, start wavelength_index %d is negative.\n", Input.wl.start_index);
    return -1;
  }

  /* check the wavelength range */
  if (Input.wl.end>850 && (Input.ck_scheme==CK_NONE && strlen(Input.filename[FN_MOL_TAU_ABS])==0))  {
    status=-1;
    fprintf (stderr, "Error %d, you want to do a spectral calculation for wavelengths larger than 850 nm. While uvspec\n", status);
    fprintf (stderr, "    treats ozone absorption correctly, molecular absorption is NOT considered in monochromatic\n");
    fprintf (stderr, "    uvspec calculations, as absorption cross-section are highly variable with wavelength.\n");
    fprintf (stderr, "    To consider molecular absorption other than ozone you have two choices \n");
    fprintf (stderr, "    with uvspec:\n");
    fprintf (stderr, "     (1) Do a line-by-line calculation using 'mol_tau_file abs' to specify\n");
    fprintf (stderr, "         the wavelength-dependent absorption profile; to calculate the\n");
    fprintf (stderr, "         latter, you need something like David Edwards' genln2.\n");
    fprintf (stderr, "         ATTENTION: line-by-line calculations are very time-consuming!\n");
    fprintf (stderr, "     (2) Use the correlated-k approximation which is the most accurate\n");
    fprintf (stderr, "         solution after the line-by-line calculation; use either the\n");
    fprintf (stderr, "         pre-defined parameterization that come with libRadtran or provide\n");
    fprintf (stderr, "         your own; both options are selected with 'mol_abs_param ...'\n");
    fprintf (stderr, "\n");
    return -1;
  }

  if (Input.wl.end>500 && (Input.ck_scheme==CK_NONE && strlen(Input.filename[FN_MOL_TAU_ABS])==0))  {
    if (!Input.quiet)  {
      fprintf (stderr, "\n");
      fprintf (stderr, "*** WARNING, you want to do a spectral calculation for wavelengths larger\n");
      fprintf (stderr, "*** than 500 nm. Be warned that from about 500 nm upwards there are highly\n");
      fprintf (stderr, "*** structured molecular absorption bands which are NOT considered in the\n");
      fprintf (stderr, "*** spectral uvspec calculation. While uvspec treats ozone absorption\n");
      fprintf (stderr, "*** correctly, the water vapour absorption bands around 510nm, 545nm,\n");
      fprintf (stderr, "*** 595nm, 655nm, 725nm, 825nm, 935nm, and higher, as well as the\n");
      fprintf (stderr, "*** oxygen absorption at 630nm, 690nm, and 765nm are not implemented.\n");
      fprintf (stderr, "*** Below around 600nm these may usually be neglected; above 600nm it\n");
      fprintf (stderr, "*** starts getting more serious, and above 800nm it does not make any\n");
      fprintf (stderr, "*** sense to do the RTE calculation without considering molecular\n");
      fprintf (stderr, "*** absorption.\n");
      fprintf (stderr, "*** \n");
      fprintf (stderr, "*** To consider molecular absorption other than ozone you have two choices\n");
      fprintf (stderr, "*** with uvspec:\n");
      fprintf (stderr, "***  (1) Do a line-by-line calculation using 'mol_tau_file abs' to specify\n");
      fprintf (stderr, "***      the wavelength-dependent absorption profile; to calculate the\n");
      fprintf (stderr, "***      latter, you need something like David Edwards' genln2.\n");
      fprintf (stderr, "***      ATTENTION: line-by-line calculations are very time-consuming!\n");
      fprintf (stderr, "***  (2) Use the correlated-k approximation which is the most accurate\n");
      fprintf (stderr, "***      solution after the line-by-line calculation; use either the\n");
      fprintf (stderr, "***      pre-defined parameterization that come with libRadtran or provide\n");
      fprintf (stderr, "***      your own; both options are selected with 'mol_abs_param ...'\n");
      fprintf (stderr, "\n");
    }
  }


  /**************************************/
  /* initialize random number generator */
  /**************************************/
  /* initialize random number generator */

  if (Input.rte.mc.readrandomseed) {
    /* read random seed for debugging purposes */
    if ((randomseedfile = fopen(randomseedfilename, "r")) != NULL) {
      fscanf(randomseedfile,"%d",&rseed);
      (void) fclose (randomseedfile);
    }
  }

  rseed = init_uvspec_random (&Input.rte.mc.randomseed, Input.rte.solver == SOLVER_MONTECARLO, Input.quiet);
  if (rseed<0) {
    fprintf (stderr, "Error %d returned from init_libRadtran_random\n", rseed);
    return rseed;
  }

  /* store random seed for debugging purposes */
  if ((randomseedfile = fopen(randomseedfilename, "w")) != NULL) {

    /* tmpnam generates nonexistent filename, but by the time the file is opened this may no longer be true */
    /* mkstemp should fix this, but does not work properly. instead, I create the file as quickly as possible */

    /* strcpy("/tmp/rand-XXXXXX", Input.filename[FN_RANDOMSTATUS]); */
    /* ensure exclusive access on buggy implementations of mkstemp */
    /* umask(0077); */
    /* int fd = mkstemp(Input.filename[FN_RANDOMSTATUS]); */
    /* close(fd); */

    fprintf (randomseedfile, "%d\n", rseed);

#ifdef WRITERANDOMSTATUS
    if (Input.rte.mc.readrandomstatus) {
      Input.filename[FN_RANDOMSTATUS] = "randomstatus";
    } else {
      tmpnam(Input.filename[FN_RANDOMSTATUS]);
      if ( ( randomstatusfile = fopen(Input.filename[FN_RANDOMSTATUS], "w")) != NULL)
        (void) fclose (randomstatusfile);
      else {
        fprintf (stderr, "Error during creation of uvspec randomstatusfile %s\n", Input.filename[FN_RANDOMSTATUS] );
        return -1;
      }
    }

    fprintf (randomseedfile, "Randomstatusfile %s\n", Input.filename[FN_RANDOMSTATUS]);
#endif

    (void) fclose (randomseedfile);
  }
    

  /********************************/
  /* Do the radiation calculation */
  /********************************/
  
  status = uvspec(Input,&Output);
  if (status!=0) {
    fprintf (stderr, "Error %d during execution of uvspec\n", status);
    return status;
  }

  
#ifndef UVSPEC_FUNCTION
  print_output(Input, &Output); 
#else
  print_output(Input, &Output, outfilename); 
#endif
  if (status!=0) {
    fprintf (stderr, "Error %d during print_output \n", status);
    return status;
  }

  free_output(Input, &Output);
  return status;

}

  /*************************/
  /* Output of the results */
  /*************************/

#ifndef UVSPEC_FUNCTION
static int print_output(input_struct Input, output_struct *output) 
#else 
static int print_output(input_struct Input, output_struct *output, char *outfilename) 
#endif
{
int status=0;

#ifndef UVSPEC_FUNCTION
  FILE *f;
  switch(Input.output_format) {
  case OUTPUT_FORMAT_FLEXSTOR:
  case OUTPUT_FORMAT_ASCII:
  case OUTPUT_FORMAT_MATRIX:
          /*If OUTPUT_FORMAT_ASCII is used with a filename then open this file so that the output*/
          /*goes there instead of stdout*/
          if (strlen(Input.filename[FN_OUTFILE]) != 0) {
                  f = fopen (Input.filename[FN_OUTFILE], "w");
          }
          else{
                  f = stdout;
          }
          break;
  default:
          f = stdout;
  }
#else
  FILE *f = fopen (outfilename, "w");
#endif


  switch(Input.output_format) {

  case OUTPUT_FORMAT_FLEXSTOR:
    write_flexstor_output(Input, output, f);
    break;
  case OUTPUT_FORMAT_ASCII:
    
    if (Input.header)
      status=write_header(Input, output, f);
      if (status) {
        fprintf (stderr, "Error, returned from write_header (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	return -1;
      }

    if (Input.rte.solver != SOLVER_POLRADTRAN) {
      /* write user defined output to file with pointer f */
      status=write_user_defined_output(Input, output, f);
      if (status) {
        fprintf (stderr, "Error, returned from write_user_defined_output (line %d, function '%s' in '%s')\n", __LINE__, __func__, __FILE__);
	return -1;
      }
    }
    else { /* -> Input.rte.solver == SOLVER_POLRADTRAN) */
      /* write output in a special Polradtran format */
      write_polradtran_output(Input, output, f);
    }
    break;
  case OUTPUT_FORMAT_MATRIX:

    /* special output for spectral heating rates as 2D matrix */
    write_spectral_heat (Input, output, f);

    break; 
  case OUTPUT_FORMAT_NETCDF:
  case OUTPUT_FORMAT_SAT_PICTURE:

#if HAVE_LIBNETCDF

    status = write_netCDF_output(Input, output);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d during write_netCDF_output() to %s (uvspec_lex.l) \n", status, Input.filename[FN_OUTFILE]);
      return -1;
    }
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use the netCDF output option. Please get netcdf and rebuild.        *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

    break;
  default:
    fprintf (stderr, "Error, unknown output format %d (line %d, function '%s' in '%s')\n", Input.output_format, __LINE__, __func__, __FILE__);
    return -1; 
  }


  fflush(f);

  if (f!=stdout)
    fclose(f);

  return 0;
}  


void write_flexstor_output(input_struct Input, output_struct *output, FILE *f)
{

  int lev=0;
  int iv=0;

  fprintf (f, "flexstor 1 0 suvdama 1 1\n");
  fprintf (f, "# This flexstor output has not been checked by flxcheck.\n");
  fprintf (f, "# Hence it may not follow flexstor/suvdama rules. This will\n");
  fprintf (f, "# be improved soon.\n");
  fprintf (f, "model_name() uvspec, version %s\n", THEVERSION);
  if (fabs(Input.aer.beta-0.0) > 0.0) {
    fprintf (f, "angstrom_alpha()  %f\n", Input.aer.alpha);
    fprintf (f, "angstrom_beta()  %f\n",  Input.aer.beta);
  }
  fprintf (f, "mol_modify O3 (DU) %f\n", output->column[MOL_O3]);
  fprintf (f, "pressure(mb) %f\n", Input.pressure);
  fprintf (f, "self_day_of_year() %d\n", Input.UTC.tm_yday);
  
  /*
  fprintf (f, "self_zout(km)");
  for (lev=0;lev<output->atm.nzoutatm.nzout;lev++)
    fprintf (f, "  %f", Input.atm.zout[lev]);
  fprintf (f, "\n");
  */

  fprintf (f, "self_temper(K)");
  for (lev=0;lev<output->atm.nzout;lev++)
    fprintf (f, "  %f", output->atm.temper_out[lev]);
  
  fprintf (f, "\n");
  fprintf (f, "self_day_of_year() %d\n", Input.UTC.tm_yday);
  fprintf (f, "self_eccentricity_factor() %f\n", output->eccentricity_factor);
  fprintf (f, "begin table data\n");
  for (lev=0;lev<output->atm.nzout;lev++) {
    fprintf (f, "WL(nm) SZA(deg.dec) self_rfldir(mW m-2 nm-2) self_rfldn(mW m-2 nm-2) self_rflup(mW m-2 nm-2)");
    fprintf (f, " self_uavgso(mW m-2 nm-2) self_uavgdn(mW m-2 nm-2) self_uavgup(mW m-2 nm-2) ");
  }
  fprintf (f, "\n");
  for (iv=0;iv<output->wl.nlambda_h;iv++) {
    for (lev=0;lev<output->atm.nzout;lev++) {
      if (Input.rte.solver == SOLVER_FTWOSTR) {
        fprintf (f, "%7.3f %13.6e %13.6e %13.6e %13.6e ", 
        	 output->wl.lambda_h[iv], 
		 output->rfldir[lev][iv], 
	         output->rfldn[lev][iv], 
	         output->flup[lev][iv],
	         output->uavg[lev][iv]);
      }
      else {
        fprintf (f, "%7.3f %7.3f %8.5f %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e",
      	           output->wl.lambda_h[iv], 
		   output->sza_h[iv], 
		   output->atm.zout_sur[lev],
		   output->rfldir[lev][iv], 
		   output->rfldn[lev][iv], 
		   output->flup[lev][iv],
		   output->uavgso[lev][iv], 
		   output->uavgdn[lev][iv], 
		   output->uavgup[lev][iv]);
      }
      fprintf (f, "\n");
    }
  }
  fprintf (f, "end table data\n");
}


int write_header (input_struct Input, output_struct *output, FILE *f)
{
  int lev=0;

  fprintf (f, "#aerosol_haze                 : %16d\n",   Input.aer.haze);
  fprintf (f, "#aerosol_season               : %16d\n",   Input.aer.seasn);
  fprintf (f, "#aerosol_modify gg  set       : %16.3f\n", Input.aer.modify[MODIFY_VAR_GG][MODIFY_TYPE_SET]);
  fprintf (f, "#aerosol_modify ssa set       : %16.3f\n", Input.aer.modify[MODIFY_VAR_SSA][MODIFY_TYPE_SET]);
  fprintf (f, "#aerosol_file tau             : %s\n",     Input.aer.filename[FN_AER_TAU]);
  fprintf (f, "#aerosol_file ssa             : %s\n",     Input.aer.filename[FN_AER_SSA]);
  fprintf (f, "#aerosol_file gg              : %s\n",     Input.aer.filename[FN_AER_GG]);
  fprintf (f, "#aerosol_file moments         : %s\n",     Input.aer.filename[FN_AER_MOMENTS]);
  fprintf (f, "#aerosol_visibility           : %16.3f\n", Input.aer.visibility);
  fprintf (f, "#albedo                       : %16.3f\n", Input.alb.albedo);
  fprintf (f, "#aerosol_angstrom             : %16.3f  %16.3f\n", Input.aer.alpha, Input.aer.beta);
  fprintf (f, "#atmosphere_file              : %s\n",     Input.filename[FN_ATMOSPHERE]);
  fprintf (f, "#day_of_year                  : %16d\n",   Input.UTC.tm_yday);

  switch (Input.rte.deltam) {
  case 0:
    printf("#deltam                       : off\n");
    break;
  case 1:
    printf("#deltam                       : on\n");
    break;
  default:
    printf("#deltam                       : unknown\n");
    break;
  }

  switch (Input.crs_model[CRS_MOL_O3]) {
  case CRS_MODEL_BASS_AND_PAUR:
    printf("#crs_model o3                 : Bass_and_Paur\n");
    break;
  case CRS_MODEL_MOLINA:
    printf("#crs_model o3                 : Molina\n");
    break;
  case CRS_MODEL_DAUMONT:
    printf("#crs_model o3                 : Daumont\n");
    break;
  default:
    printf("#crs_model o3                 : unknown\n");
    break;
  }

  fprintf (f, "#mol_modify O3             : %16.3f\n", Output.column[MOL_O3]);

  switch (Input.rte.solver) {
  case SOLVER_SDISORT:
    printf("#rte_solver                   : sdisort\n");
    break;
  case SOLVER_SPSDISORT:
    printf("#rte_solver                   : spsdisort\n");
    break;
  case SOLVER_FDISORT1:
    printf("#rte_solver                   : fdisort1\n");
    break;
  case SOLVER_FDISORT2:
    printf("#rte_solver                   : fdisort2\n");
    break;
  case SOLVER_DISORT:
    printf("#rte_solver                   : disort\n");
    break;
  case SOLVER_MONTECARLO:
    printf("#rte_solver                   : mystic\n");
    break;
  case SOLVER_POLRADTRAN:
    printf("#rte_solver                   : polradtran\n");
    break;
  case SOLVER_SOS:
    printf("#rte_solver                   : sos\n");
    break;
  case SOLVER_TZS: 
    printf("#rte_solver                   : tzs\n");
    break;
  case SOLVER_SSS: 
    printf("#rte_solver                   : sss\n");
    break;
  case SOLVER_SSSI: 
    printf("#rte_solver                   : sssi\n");
    break;
  case SOLVER_FTWOSTR: 
    printf("#rte_solver                   : ftwostr\n");
    break;
  case SOLVER_TWOSTR: 
    printf("#rte_solver                   : twostr\n");
    break;
  case SOLVER_RODENTS:
    printf("#rte_solver                   : rodents\n");
    break;
  case SOLVER_TWOSTREBE:
    printf("#rte_solver                   : twostrebe\n");
    break;
  case SOLVER_SSLIDAR:
    printf("#rte_solver                   : sslidar\n");
    break;
  default:
    printf("#rte_solver                   : unknown (%d)\n",Input.rte.solver);
    return -1;
  }

  fprintf (f, "#slit_function_file           : %s\n", 
           Input.filename[FN_SLITFUNCTION]);
  fprintf (f, "#source solar file                   : %s\n", 
           Input.filename[FN_EXTRATERRESTRIAL]);
  fprintf (f, "#wavelength                          : %16.3f  %16.3f\n", 
           Input.wl.start, Input.wl.end);
  fprintf (f, "#zout                         : ");
  for (lev=0;lev<output->atm.nzout;lev++)
    fprintf (f, " %16.3f", Input.atm.zout_sur[lev]);
  fprintf (f, "\n");

  return 0;
}


void write_polradtran_output (input_struct Input, output_struct *output, FILE *f)
{

  int iv=0;
  int lev=0;
  int is=0;  /* counter for stokes element (1...4) */
  int iu=0;  /* counter for cos(SZA) */

  for (iv=0;iv<output->wl.nlambda_h;iv++) {
    for (lev=0;lev<output->atm.nzout;lev++) {

      fprintf (f, "%7.3f ", output->wl.lambda_h[iv]);

      for (is=0;is<Input.rte.polradtran[POLRADTRAN_NSTOKES];is++)
        fprintf (f, "%13.6e %13.6e ",
                 output->down_flux[lev][is][iv], 
	         output->up_flux[lev][is][iv]);

      fprintf (f, "\n");

      if (output->print_phi && Input.rte.nphi > 0) {
        fprintf (f, "                      ");
        for (j=0; j<Input.rte.nphi; j++)
          fprintf (f, " %13.3f ", Input.rte.phi[j]);
        fprintf (f, "\n");
        for (is=0; is<Input.rte.polradtran[POLRADTRAN_NSTOKES]; is++) {
          fprintf (f, "Stokes vector %c\n", Input.rte.pol_stokes[is]);
          for (iu=0; iu<Input.rte.numu; iu++) {
            /* printf("%7.4f ", output->mu_values[iu+Input.rte.nstr/2]);*/
	    printf("%7.4f ", Input.rte.umu[iu]);
	    printf("%13.6e ", 0.0);
	    for (j=0; j<Input.rte.nphi; j++) {
	      if (Input.rte.umu[iu] <= 0)
	        fprintf (f, " %13.6e ", output->down_rad[lev][j][iu][is][iv]);
              else
	        fprintf(f, " %13.6e ", output->up_rad[lev][j][iu][is][iv]);}
	    printf("\n");
	  }
        }
      }
    }
  }
}

void write_spectral_heat (input_struct Input, output_struct *output, FILE *f)
{
  int iv=0;
  int lev=0;
  int lev2=0;

  for (iv=0;iv<output->wl.nlambda_h;iv++) {
    for (lev=0;lev<output->atm.nzout;lev++) {

      /* spectral output */
      if (iv == 0 && lev == 0) {
        fprintf (f, "%9.3e ", 0.0);
        for (lev2=0;lev2<output->atm.nzout;lev2++)
          fprintf (f, "%13.6f ", output->atm.zout_sur[lev2]);
        fprintf (f, "\n");
      }
      if (lev == 0)
        fprintf (f, "%9.3f ", output->wl.lambda_h[iv]);
      fprintf (f, "%13.6e ", output->heat[lev][iv]);
      if (lev == output->atm.nzout-1)
        fprintf (f, "\n");
    }
  }
}


int write_user_defined_output(input_struct Input, output_struct *output, FILE *f)
{
  int io    = NOT_DEFINED_INTEGER;
  int i_gas = NOT_DEFINED_INTEGER;
  int iv    = NOT_DEFINED_INTEGER;
  int lev   = NOT_DEFINED_INTEGER;
  int llev  = NOT_DEFINED_INTEGER;
  int iu    = NOT_DEFINED_INTEGER;
  int iiu   = NOT_DEFINED_INTEGER;
  int numu  = NOT_DEFINED_INTEGER;
  int write_rad = TRUE; 

  for (iv=0;iv<output->wl.nlambda_h;iv++) {
    for (llev=0;llev<output->atm.nzout_user;llev++) {
      if ( Input.raman ) lev = output->atm.zout_user_index[llev];
      else lev =llev;

      for (io=0; io<Input.n_output_user; io++) {
        switch (Input.output_user[io]) {
        case OUTPUT_USER_WAVE:
          fprintf (f, "%9.3f ", output->wl.lambda_h[iv]);
          break;
        case OUTPUT_USER_WAVE_MAX:
          fprintf (f, "%9.3f ", 1.0e7/output->wl.wvnmlo_h[iv]); 
          break;
        case OUTPUT_USER_WAVE_MIN:
          fprintf (f, "%9.3f ", 1.0e7/output->wl.wvnmhi_h[iv]); 
          break;
        case OUTPUT_USER_WAVENUMBER:
          fprintf (f, "%9.3f ", 1E7/output->wl.lambda_h[iv]);
          break;
        case OUTPUT_USER_WAVENUMBER_MAX:
          fprintf (f, "%9.3f ", output->wl.wvnmhi_h[iv]);        
          break;
        case OUTPUT_USER_WAVENUMBER_MIN:
          fprintf (f, "%9.3f ", output->wl.wvnmlo_h[iv]);        
          break;
        case OUTPUT_USER_ZOUT_SUR:
          fprintf (f, "%11.6f ", output->atm.zout_sur[lev]);
          break;
        case OUTPUT_USER_ZOUT_SEA:
          fprintf (f, "%11.6f ", output->atm.zout_sea[lev]);
          break;
        case OUTPUT_USER_Z_SUR:
          fprintf (f, "%11.6f ", output->alt.altitude);
          break;
        case OUTPUT_USER_EDIR:
          fprintf (f, "%13.6e ", output->rfldir[lev][iv]);
          break;
        case OUTPUT_USER_EGLO:
          fprintf (f, "%13.6e ", output->rfldir[lev][iv]+output->rfldn[lev][iv]);
          break;
        case OUTPUT_USER_EDN:
          fprintf (f, "%13.6e ", output->rfldn[lev][iv]);
          break;
        case OUTPUT_USER_EUP:
          fprintf (f, "%13.6e ", output->flup[lev][iv]);
          break;
        case OUTPUT_USER_ENET:
          fprintf (f, "%13.6e ", output->rfldir[lev][iv]+output->rfldn[lev][iv]-output->flup[lev][iv]);
          break;
        case OUTPUT_USER_ESUM:
          fprintf (f, "%13.6e ", output->rfldir[lev][iv]+output->rfldn[lev][iv]+output->flup[lev][iv]);
          break;
        case OUTPUT_USER_FDIR:
          fprintf (f, "%13.6e ", 4.0*PI*output->uavgso[lev][iv]);
          break;
        case OUTPUT_USER_FGLO:
          fprintf (f, "%13.6e ", 4.0*PI*(output->uavgso[lev][iv]+output->uavgdn[lev][iv]));
          break;
        case OUTPUT_USER_FDN:
          fprintf (f, "%13.6e ", 4.0*PI*output->uavgdn[lev][iv]);
          break;
        case OUTPUT_USER_FUP:
          fprintf (f, "%13.6e ", 4.0*PI*output->uavgup[lev][iv]);
          break;
        case OUTPUT_USER_SPHER_ALB:
	  fprintf (f, "%13.6e ", output->albmed[iv]);
	  break;
        case OUTPUT_USER_UU:
   	  if ( Input.raman ) numu = Input.rte.maxumu-Input.rte.nstr; 
	  else numu = Input.rte.numu;
          for (iu=0;iu<numu;iu++){
	    if ( Input.raman ) iiu = Input.rte.umuind[iu];
	    else iiu = iu;
            for (j=0;j<Input.rte.nphi;j++) {
              fprintf (f, " %16.9e ", output->uu[lev][j][iiu][iv]);
            }
          }
          write_rad = FALSE;
          break;
        case OUTPUT_USER_F:
          switch(Input.rte.solver) {
          case SOLVER_MONTECARLO:
            fprintf (f, "%13.6e ", 4.0*PI*(output->uavgso[lev][iv]+output->uavgdn[lev][iv]+output->uavgup[lev][iv]));
            break;
          default:
            fprintf (f, "%13.6e ", 4.0*PI*output->uavg[lev][iv]);
          }
	  break;
        case OUTPUT_USER_UDIR:
          fprintf (f, "%13.6e ", output->uavgso[lev][iv]);
          break;
        case OUTPUT_USER_UGLO:
          fprintf (f, "%13.6e ", output->uavgso[lev][iv]+output->uavgdn[lev][iv]);
          break;
        case OUTPUT_USER_UDN:
          fprintf (f, "%13.6e ", output->uavgdn[lev][iv]);
          break;
        case OUTPUT_USER_UUP:
          fprintf (f, "%13.6e ", output->uavgup[lev][iv]);
          break;
        case OUTPUT_USER_U:
          fprintf (f, "%13.6e ", output->uavg[lev][iv]);
          break;
        case OUTPUT_USER_ALB:
	  fprintf (f, "%13.6e ", output->flup[lev][iv]/(output->rfldir[lev][iv]+output->rfldn[lev][iv]));
          break;
        case OUTPUT_USER_SZA:
          fprintf (f, "%9.3f ", output->sza_h[iv]);
          break;
        case OUTPUT_USER_C_P:
          fprintf (f, "%13.6e ", output->atm.microphys.c_p[lev]);
          break;
        case OUTPUT_USER_HEAT:
          fprintf (f, "%13.6e ", output->heat[lev][iv]);
          break;
        case OUTPUT_USER_EMIS:
          fprintf (f, "%13.6e ", output->emis[lev][iv]);
          break;
        case OUTPUT_USER_ABS:
          fprintf (f, "%13.6e ", output->heat[lev][iv] - output->emis[lev][iv]);
          break;
        case OUTPUT_USER_W_RAD:
          fprintf (f, "%13.6e ", output->w_zout[lev][iv]/3600.0); /* m/day -> m/s */
          break;
        case OUTPUT_USER_M_RAD:
          fprintf (f, "%13.6e ", output->w_zout[lev][iv] * output->atm.microphys.dens_zout[MOL_AIR][lev] 
                                 * 1.e+6 * 1.e-3 * Input.atm.mol_mass[MOL_AIR] / avogadro );   /* kg/ (m2 day) */
          break;
        case OUTPUT_USER_P:
          fprintf (f, "%14.9f ", output->atm.microphys.press_zout[lev]);
          break;
        case OUTPUT_USER_T:
          fprintf (f, "%11.5f ", output->atm.microphys.temper_zout[lev]);
          break;
        case OUTPUT_USER_T_D:
          fprintf (f, "%11.5f ", output->atm.microphys.temper_d_zout[lev]);
          break;
        case OUTPUT_USER_T_SUR:
          fprintf (f, "%11.5f ", output->surface_temperature);
          break;
        case OUTPUT_USER_THETA:
          fprintf (f, "%11.5f ", output->atm.microphys.theta_zout[lev]);
          break;
        case OUTPUT_USER_THETA_E:
          fprintf (f, "%11.5f ", output->atm.microphys.theta_e_zout[lev]);
          break;
        case OUTPUT_USER_CLWC: /* (cwd[g/m3]/1000.)/rho_air[kg/m3] == cwc[kg/kg] */
          if (Input.i_wc!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.lwc_zout[lev] * 1.e-3 / 
			(output->atm.microphys.dens_zout[MOL_AIR][lev]*1.e+6*1.e-3*Input.atm.mol_mass[MOL_AIR]/avogadro) );
          else
            fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_CLWD:
          if (Input.i_wc!=-1)
	    fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.lwc_zout[lev]);
	  else
	    fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_TCLW:
          if (Input.i_wc!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.tcw);
	  else
	    fprintf (f, "%13.6e ", 0.0 );
	  break;
        case OUTPUT_USER_REFF_WAT:
          if (Input.i_wc!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.effr_zout[lev]);
	  else
	    fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_CIWC: /* (cwd[g/m3]/1000.)/rho_air[kg/m3] == cwc[kg/kg] */
          if (Input.i_ic!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.lwc_zout[lev] * 1.e-3 / 
			(output->atm.microphys.dens_zout[MOL_AIR][lev]*1.e+6*1.e-3*Input.atm.mol_mass[MOL_AIR]/avogadro) );
	  else
	    fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_CIWD:
          if (Input.i_ic!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.lwc_zout[lev]);
	  else
	    fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_TCIW:
          if (Input.i_ic!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.tcw);
	  else
	    fprintf (f, "%13.6e ", 0.0 );
	  break;
        case OUTPUT_USER_REFF_ICE:
          if (Input.i_ic!=-1) {
            if (Input.caoth[isp].properties==PROP_FU && Input.caoth[isp].fu2yang) 
              fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.effr_zout[lev]*(3.0*sqrt(3.0)/4.0));
            else
              fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.effr_zout[lev]);
	  }
	  else
	    fprintf (f, "%13.6e ", 0.0 );
          break;
        case OUTPUT_USER_TCW:
          if (Input.i_wc!=-1 && Input.i_ic!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.tcw + output->caoth[Input.i_ic].microphys.tcw);
	  else if (Input.i_wc!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_wc].microphys.tcw);
	  else if (Input.i_ic!=-1)
            fprintf (f, "%13.6e ", output->caoth[Input.i_ic].microphys.tcw);
          else
	    fprintf (f, "%13.6e ", 0.0 );
	  break;
        case OUTPUT_USER_CC:
          fprintf (f, "%13.6e ", output->cf.cf_zout[lev]);
          break;
        case OUTPUT_USER_TCC:
          fprintf (f, "%13.6e ", output->cf.cf_total);
          break;
        case OUTPUT_USER_CLOUDS:
          if (Input.i_wc!=-1 && Input.i_ic!=-1)
            fprintf (f, "%13d ", 3);
	  else if (Input.i_wc!=-1)
            fprintf (f, "%13d ", 2); 
	  else if (Input.i_ic!=-1)
            fprintf (f, "%13d ", 1);
          else
	    fprintf (f, "%13d ", 0 );
          break;
        case OUTPUT_USER_N:
          i_gas = Input.output_user_gas[io];
          fprintf (f, "%13.6e ", output->atm.microphys.dens_zout[i_gas][lev]);
          break;
        case OUTPUT_USER_RHO:   /* 1.e+6: convert from cm-3 to m-3; 1.e-3: convert g -> kg */
          i_gas = Input.output_user_gas[io];
          fprintf (f, "%13.6e ", output->atm.microphys.dens_zout[i_gas][lev] 
                                 * 1.e+6 * 1.e-3 * Input.atm.mol_mass[i_gas] / avogadro);                   
          break;
        case OUTPUT_USER_VMR:
          i_gas = Input.output_user_gas[io];
          fprintf (f, "%13.6e ", output->atm.microphys.dens_zout[i_gas][lev]/
                                 output->atm.microphys.dens_zout[MOL_AIR][lev]);
          break;
        case OUTPUT_USER_MMR:
          i_gas = Input.output_user_gas[io];
          fprintf (f, "%13.6e ", (output->atm.microphys.dens_zout[i_gas][lev]*
                                  Input.atm.mol_mass[i_gas]) /
                                  (output->atm.microphys.dens_zout[MOL_AIR][lev]*
                                  Input.atm.mol_mass[MOL_AIR]));
          break;
        case OUTPUT_USER_RH:
          fprintf (f, "%13.6e ", output->atm.microphys.dens_zout[MOL_H2O][lev]/
                                 vapor_pressure(output->atm.microphys.temper_zout[lev]) * 100.0 );
          break;
        case OUTPUT_USER_RH_ICE:
          fprintf (f, "%13.6e ", output->atm.microphys.dens_zout[MOL_H2O][lev]/
                                 vapor_pressure_over_ice(output->atm.microphys.temper_zout[lev]) * 100.0 );
          break;
        case OUTPUT_USER_WIND_U:
          fprintf (f, "%13.6e ", output->wind.u_zout[lev] );  /* in m/s */
          break;
        case OUTPUT_USER_WIND_V:
          fprintf (f, "%13.6e ", output->wind.v_zout[lev] );  /* in m/s */
          break;
        case OUTPUT_USER_WIND_W:
          fprintf (f, "%13.6e ", output->wind.w_zout[lev] );  /* in Pa/s !!! in m/s */
          break;
        case OUTPUT_USER_DTDX:
          fprintf (f, "%13.6e ", output->atm.microphys.dtheta_dx[lev] );  /* in K/m */
          break;
        case OUTPUT_USER_DTDY:
          fprintf (f, "%13.6e ", output->atm.microphys.dtheta_dy[lev] );  /* in K/m */
          break;
        case OUTPUT_USER_DTDZ:
          fprintf (f, "%13.6e ", output->atm.microphys.dtheta_dz[lev] );  /* in K/m */
          break;
        case OUTPUT_USER_HEAT_AD_X:
          fprintf (f, "%13.6e ", output->wind.u_zout[lev] * output->atm.microphys.dtheta_dx[lev] *24*60*60 );     /* in K/day */
          break;
        case OUTPUT_USER_HEAT_AD_Y:
          fprintf (f, "%13.6e ", output->wind.v_zout[lev] * output->atm.microphys.dtheta_dy[lev] *24*60*60 );     /* in K/day */
          break;
        case OUTPUT_USER_HEAT_AD_Z:
          fprintf (f, "%13.6e ", output->wind.w_zout[lev] * output->atm.microphys.dtheta_dz[lev] *24*60*60 );     /* in K/day */
          break;
        case OUTPUT_USER_HEAT_AD:
          fprintf (f, "%13.6e ", (output->wind.u_zout[lev] * output->atm.microphys.dtheta_dx[lev] + 
                                  output->wind.v_zout[lev] * output->atm.microphys.dtheta_dy[lev] +
                                  output->wind.w_zout[lev] * output->atm.microphys.dtheta_dz[lev]  ) *24*60*60 ); /* in K/day */
          break;
	case OUTPUT_USER_SSLIDAR_NPHOT:
          fprintf (f, "%13.6e ", output->sslidar_nphot[lev][iv]);
          break;
	case OUTPUT_USER_SSLIDAR_NPHOT_Q:
          fprintf (f, "%13.6e ", output->sslidar_nphot_q[lev][iv]);
          break;
	case OUTPUT_USER_SSLIDAR_RATIO:
          fprintf (f, "%13.6e ", output->sslidar_ratio[lev][iv]);
          break;
        default:
          fprintf (stderr, "Error, unknown user output %d\n", Input.output_user[io]);
	  return -1;
        }
      }
      if (Input.n_output_user>0)
        fprintf (f, "\n");


      /* radiance output */
      if ( Input.output_user_flag == 0) {  /* User has not specified any output. We supply what we think the 
                                              user wants. At least this worked for us as users for many years. 
                                              AK20091116  */ 
      
         if (Input.rte.numu > 0 && !output->print_phi) {
	    if ( Input.raman ) {
	      numu = Input.rte.maxumu-Input.rte.nstr; 
   	    }
	    else {
	      numu = Input.rte.numu;
	    }
            for (iu=0;iu<numu;iu++) {
	      if ( Input.raman ) iiu = Input.rte.umuind[iu];
	      else iiu = iu;
              fprintf (f, "%7.4f %16.6e \n", Input.rte.umu[iiu], output->u0u[lev][iiu][iv]);
	    }
         }
         if (output->print_phi && Input.rte.numu> 0 && write_rad == TRUE) {
	   switch (Input.rte.solver) {
	   case SOLVER_FDISORT1:
           case SOLVER_FDISORT2:
           case SOLVER_DISORT:
           case SOLVER_SDISORT:
	   case SOLVER_SPSDISORT:
	   case SOLVER_SSS:
	   case SOLVER_SSSI:
	   case SOLVER_TZS:
	     if ( Input.raman ) {
	       numu = Input.rte.maxumu-Input.rte.nstr; 
  	     }
  	     else {
  	       numu = Input.rte.numu;
  	     }
             fprintf (f, "                      ");
             for (j=0;j<Input.rte.nphi;j++)
               fprintf (f, " %16.3f ", Input.rte.phi[j]);
             fprintf (f, "\n");
             for (iu=0;iu<numu;iu++) {
  	     if ( Input.raman ) iiu = Input.rte.umuind[iu];
  	     else iiu = iu;
               fprintf (f, "%7.4f %16.9e ", Input.rte.umu[iiu], output->u0u[lev][iiu][iv]);
               for (j=0;j<Input.rte.nphi;j++)
                  fprintf (f, " %16.9e ", output->uu[lev][j][iiu][iv]);
               fprintf (f, "\n");
             }
	     break;
	   case SOLVER_FTWOSTR:
	   case SOLVER_SOS:
	   case SOLVER_MONTECARLO:
	   case SOLVER_POLRADTRAN:
	   case SOLVER_NULL:
	   case SOLVER_RODENTS:
	   case SOLVER_TWOSTREBE:
	   case SOLVER_TWOSTR:
	   case SOLVER_SSLIDAR:
	     break;
	   default:
	     fprintf (stderr, "Error, unknown rte_solver %d (line %d, function '%s' in '%s')\n", Input.rte.solver, __LINE__, __func__, __FILE__);
	     return -1;
	   }
         }
       } /* radiance output */

    } /* for (lev=0;lev<output->atm.nzout;lev++) */
  } /* for (iv=0;iv<output->wl.nlambda_h;iv++) */ 

  return 0;
}



int write_netCDF_output (input_struct Input, output_struct *output)
{

#if HAVE_LIBNETCDF

  int status=0;
  int io=0, iv=0, lev=0, i_gas=0;

  int ncid=0;

  int *id_dim  = NULL;
  int id_nx      = 0;
  int id_ny      = 0;
  int id_ntime   = 0;

  char dim_x_name[10];
  char dim_y_name[10];
  char x_name[10];
  char y_name[10];

  int periodic_x = NOT_DEFINED_INTEGER;
  int periodic_y = NOT_DEFINED_INTEGER;

  double x = NOT_DEFINED_FLOAT;
  double y = NOT_DEFINED_FLOAT;

  int id_nzout   = 0;
  int id_nlambda = 0;

  int *id_data  = NULL;
  int id_levels = 0;
  int id_lambda = 0;

  float *model_levels=NULL;

  int id_phi  = 0;
  int id_nphi = 0;
  int id_umu  = 0;
  int id_numu = 0;

  size_t nx   = NOT_DEFINED_INTEGER;
  size_t ny   = NOT_DEFINED_INTEGER;
  float dt = 0.0;

  double *x_grid         = NULL;
  double *y_grid         = NULL;

  int itime1=0, itime2=0, nt=0;

  int ix  = NOT_DEFINED_INTEGER;
  int iy  = NOT_DEFINED_INTEGER;
  int itime = NOT_DEFINED_INTEGER;

  size_t index_3D[3] = {0,0,0};
  size_t index_4D[4] = {0,0,0,0};
  size_t index_5D[5] = {0,0,0,0,0};
  size_t index_6D[6] = {0,0,0,0,0,0};
  size_t index_7D[7] = {0,0,0,0,0,0,0};

  size_t n_4D[4] = {0,0,0,0};
  size_t n_5D[5] = {0,0,0,0,0};
  size_t n_6D[6] = {0,0,0,0,0,0};
  size_t n_7D[7] = {0,0,0,0,0,0,0};

  int id_dim_array3[3];
  int id_dim_array4[4]; 
  int id_dim_array5[5];
  int id_dim_array6[6];
  int id_dim_array7[7];
  const int dim3=3, dim4=4, dim5=5;

  int *valid_result=NULL;    /* marker, if the RTE solver can calculate the disired quantity */
  int simulation_done=TRUE;  /* marker for valid (not NaN and not Infinity) results          */

  float tmp_float = NOT_DEFINED_FLOAT;
  float  *tmp_data_1D = NULL;
  float **tmp_data_3D = NULL;

  char *gas      = NULL;
  char var_name[20];
  char long_name[40];
  char irradiance_units[40];
  char radiance_units[40];

  char function_name[]="write_netCDF_output";
  char file_name[]="uvspec_lex.l";

  if (Input.verbose)
    fprintf (stderr, " ... start writing result to netCDF file %s\n", Input.filename[FN_OUTFILE]);

  /* open netcdf file */
  status = nc_open (Input.filename[FN_OUTFILE], NC_WRITE, &ncid); /* read and write access */
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, Input.filename[FN_OUTFILE]);
    return status;
  }

  /* ----------------------------------------------------------------------- */

  /* get identification number and read data */

  /* get id for dimension "time" (the dimension, not the variable) */
  if      (Input.time_zone==TIME_UTC)
    status = nc_inq_dimid (ncid, "time", &(id_ntime));
  else if (Input.time_zone==TIME_LAT)
    status = nc_inq_dimid (ncid, "local_apparent_time", &(id_ntime)); 
  else {
    fprintf (stderr, "Error, unknown time_zone %d (line %d, function '%s' in '%s')\n", Input.time_zone, __LINE__, __func__, __FILE__);
    return status;
  }
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d locating 'ntime' in %s \n", status, Input.filename[FN_OUTFILE]);
    return status;
  }


  switch (Input.output_format)  {
  case OUTPUT_FORMAT_NETCDF:
    strcpy (dim_x_name, "lon");
    strcpy (dim_y_name, "lat");
    strcpy (x_name, "lon");
    strcpy (y_name, "lat");
    x = Input.longitude;
    y = Input.latitude;
    periodic_x = FALSE;
    periodic_y = TRUE;
    break;
  case OUTPUT_FORMAT_SAT_PICTURE: 
    strcpy (dim_x_name, "pixel_x");
    strcpy (dim_y_name, "pixel_y");
    strcpy (x_name, "pixel_x");
    strcpy (y_name, "pixel_y");
    x = Input.sat_pixel_x;
    y = Input.sat_pixel_y;
    periodic_x = FALSE;
    periodic_y = FALSE;
    break;
  default:
    fprintf (stderr, "Error, unknown output_format %d (line %d, function '%s' in '%s')\n", Input.output_format, __LINE__, __func__, __FILE__);
    return (-1);
  }

  /* get id for dimension "lat/pixel_x" (the dimension, not the variable) */
  status = nc_inq_dimid (ncid, dim_x_name, &(id_nx));
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d locating '%s' in %s \n", status, dim_x_name, Input.filename[FN_OUTFILE]);
    return status;
  }

  /* read latitude array (must already be stored in the file) */
  status = alloc_and_read_netCDF_1D_double(ncid, dim_x_name, &nx, x_name, &(x_grid));
  if (status != 0) {
    fprintf (stderr, "Error %d reading %s from %s\n", status, x_name, Input.filename[FN_OUTFILE]);
    return status;
  }

  /* get id for dimension "lon/pixel_y" (the dimension, not the variable) */
  status = nc_inq_dimid (ncid, dim_y_name, &(id_ny));
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d locating '%s' in %s \n", status, dim_y_name, Input.filename[FN_OUTFILE]);
    return status;
  }

  /* read longitude array (must already be stored in the file) */
  status = alloc_and_read_netCDF_1D_double( ncid, dim_y_name, &ny, y_name, &(y_grid) );
  if (status != 0) {
    fprintf (stderr, "Error %d reading %s from %s\n", status, y_name, Input.filename[FN_OUTFILE]);
    return status;
  }


  /* ----------------------------------------------------------------------- */

  /* search index */

  if (Input.time_zone==TIME_UTC)
    status = get_time_index (ncid, Input.UTC, TIME_NEAREST_DATE, 
                             &(nt), &(itime1), &(itime2), &(dt),
                             Input.verbose, Input.quiet);
  else if (Input.time_zone==TIME_LAT)
    status = get_local_apparent_time_index (ncid, Input.LAT, TIME_NEAREST_DATE, 
                                            &(nt), &(itime1), &(itime2), &(dt),
                                            Input.verbose, Input.quiet);
  else {
    fprintf (stderr, "Error, unknown time_zone %d (line %d, function '%s' in '%s')\n", Input.time_zone, __LINE__, __func__, __FILE__);
    return status;
  }


  if (status != 0) {
    fprintf (stderr, "Error %d, during get_time_index in write_netCDF_output (uvspec_lex.l)\n", status);
    return status;
  }
  itime = itime1;


  /* get lat/pixel_x-index in the result grid, where to write data */
  ix = get_grid_index(x, x_grid, nx, periodic_x);
  if (ix < 0) {
    fprintf (stderr, "Error %d, while searching index in %s_grid \n", ix, x_name);
    return ix;
  }

  /* get lon/pixel_y-index in the result grid, where to write data */
  iy = get_grid_index(y, y_grid, ny, periodic_y);
  if (iy < 0) {
    fprintf (stderr, "Error %d, while searching index in %s_grid \n", iy, y_name);
    return iy;
  }



  /* ----------------------------------------------------------------------- */

  /* check if lambda exists, (write lambda if not), and get dimension id for nlambda */
  write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                         "lambda", &id_lambda, "wavelength", "nm", output->wl.lambda_h, 
                         "lambda", &id_nlambda, output->wl.nlambda_h);

  /* check if zout exists, (write zout if not), and get dimension id for zout */
  switch (Input.atm.zout_source) {
  case OUTLEVEL_ZOUT_ABOVE_SUR:
    write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "zout_sur", &id_levels, "altitude above surface", "km", output->atm.zout_sur_org, 
                           "zout_sur", &id_nzout, output->atm.nzout_org );
    break;
  case OUTLEVEL_ZOUT_ABOVE_SEA:
  case OUTLEVEL_ATM_LEVELS:
  case OUTLEVEL_ALL_LEVELS:
    write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "zout_sea", &id_levels, "altitude above sea level", "km", output->atm.zout_sea_org, 
                           "zout_sea", &id_nzout, output->atm.nzout_org );
    break;
  case OUTLEVEL_PRESS:
    write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "p", &id_levels, "pressure", "hPa", Input.atm.press_zout, 
                           "p", &id_nzout, Input.atm.nzout);
    break;
  case OUTLEVEL_MODEL_LEVELS:
  case OUTLEVEL_MODEL_LAYERS:
  case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:

    /* write simple the numbers from 1 to mlev, instead of p or z */
    if ((model_levels = calloc (output->atm.nzout_org, sizeof(float) )) == NULL) {
      fprintf (stderr,"Error, allocating memory for 'model_levels' in %s (in %s)\n", function_name, file_name);
      return -1;
    }

    for (lev=0;lev<output->atm.nzout;lev++)
      model_levels[lev]=lev+1;

    write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "mlev", &id_levels, "model_levels", "", model_levels, 
                           "mlev", &id_nzout, output->atm.nzout_org);

    break;
  default:
    fprintf (stderr, "Error, unknown zout_source %d (line %d, function '%s' in '%s')\n", Input.atm.zout_source, __LINE__, __func__, __FILE__);
    return -2;
  }

  /* ----------------------------------------------------------------------- */


  /* initialisation of the index, where to write */
  id_dim_array3[POSITION_TIME_3D] = id_ntime;  
  id_dim_array3[POSITION_LAT_3D ] = id_ny;     /* one before last dimension is default y-axis in ncview */
  id_dim_array3[POSITION_LON_3D ] = id_nx;     /*            last dimension is default x-axis in ncview */

  /* initialisation of the index, where to write */
  index_3D[POSITION_TIME_3D] = itime;
  index_3D[POSITION_LAT_3D]  = iy;               /* latitude or x-pixel */
  index_3D[POSITION_LON_3D]  = ix;               /* longitude or y-pixel */
 
  /* order of dimensions */
  id_dim_array4[POSITION_TIME_4D] = id_ntime;
  id_dim_array4[POSITION_ZOUT_4D] = id_nzout;
  id_dim_array4[POSITION_LAT_4D ] = id_ny;     /* one before last dimension is default y-axis in ncview */
  id_dim_array4[POSITION_LON_4D ] = id_nx;     /*            last dimension is default x-axis in ncview */

  /* initialisation of the index, where to write */
  index_4D[POSITION_TIME_4D] = itime;
  n_4D[POSITION_ZOUT_4D]     = output->atm.nzout; /* loop over zout */
  index_4D[POSITION_LAT_4D]  = iy;                /* latitude or x-pixel */
  index_4D[POSITION_LON_4D]  = ix;                /* longitude or y-pixel */

  /* order of dimensions */
  id_dim_array5[POSITION_TIME_5D  ] = id_ntime;
  id_dim_array5[POSITION_LAMBDA_5D] = id_nlambda;
  id_dim_array5[POSITION_ZOUT_5D  ] = id_nzout;
  id_dim_array5[POSITION_LAT_5D   ] = id_ny;  /* one before last dimension is default y-axis in ncview */
  id_dim_array5[POSITION_LON_5D   ] = id_nx;  /*            last dimension is default x-axis in ncview */

  /* initialisation of the index, where to write */
  index_5D[POSITION_TIME_5D] = itime;
  n_5D[POSITION_LAMBDA_5D]   = output->wl.nlambda_h; /* loop over lambda */    
  n_5D[POSITION_ZOUT_5D]     = output->atm.nzout;    /* loop over zout */ 
  index_5D[POSITION_LAT_5D]  = iy;                   /* latitude or x-pixel */
  index_5D[POSITION_LON_5D]  = ix;                   /* longitude or y-pixel */

  /* ----------------------------------------------------------------------- */

  /* allocate memory for id_dimensions */
  if ((id_dim  = (int *) calloc (Input.n_output_user, sizeof (int *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for id_dim  in write_netCDF_output (in uvspec_lex.l)\n");
    return -1;
  }

  /* allocate memory for id_variables */
  if ((id_data = (int *) calloc (Input.n_output_user, sizeof (int *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for id_data in write_netCDF_output (in uvspec_lex.l)\n");
    return -2;
  }

  /* allocate memory for marker, if result can be calculated by the solver */
  if ((valid_result = (int *) calloc (Input.n_output_user, sizeof (int *))) == NULL) {
    fprintf (stderr,"Error, allocating memory for valid_result in write_netCDF_output (in uvspec_lex.l)\n");
    return -2;
  }

  /* check if result exist for the specific solver */
  for (io=0; io<Input.n_output_user; io++)
    status += check_if_result_exists(Input.output_user, Input.n_output_user, Input.rte.solver, 
              &(valid_result), Input.quiet);
  if (status != 0) {
    fprintf (stderr,"Error, during check_if_result_exists in write_netCDF_output (in uvspec_lex.l)\n");
    return -2;
  }

  /* ----------------------------------------------------------------------- */

  switch(Input.processing) {
  case PROCESS_NONE:
    switch(Input.output_unit) {
    case UNIT_NOT_DEFINED:
      strcpy (irradiance_units, "unknown units");
      strcpy (radiance_units,  "unknown units");
      break;
    case UNIT_PER_NM:
      strcpy (irradiance_units, "(W or mW)/(m^2 nm)");
      strcpy (radiance_units,  "(W or mW)/(m^2 nm sr)");
      break;
    case UNIT_PER_CM_1:
      strcpy (irradiance_units, "(W or mW)/(m^2 cm^-1)");
      strcpy (radiance_units,  "(W or mW)/(m^2 cm^-1 sr)");
      break;
    case UNIT_PER_BAND:
      strcpy (irradiance_units, "(W or mW)/(m^2 band)");
      strcpy (radiance_units,  "(W or mW)/(m^2 band sr)");
      break;
      break;
    default:
      fprintf (stderr, "Error, unknown output units %d (line %d, function '%s' in '%s')\n", Input.output_unit, __LINE__, __func__, __FILE__);
      return -1;
    }
    break;
  case PROCESS_INT:
  case PROCESS_SUM:
    strcpy (irradiance_units, "(W or mW)/m^2");
    strcpy (radiance_units,  "(W or mW)/(m^2 sr)");
    break;

  case PROCESS_RGBNORM:
    strcpy (irradiance_units, "(red, green, blue) (normalised)");
    strcpy (radiance_units,   "(red, green, blue) (normalised)");
    break;

  case PROCESS_RGB:
    strcpy (irradiance_units, "(red, green, blue)");
    strcpy (radiance_units,   "(red, green, blue)");
    break;

  default:
    fprintf (stderr, "Error, unknown output processing %d (line %d, function '%s' in '%s')\n", Input.processing, __LINE__, __func__, __FILE__);
    return -1;
  }


  /* ----------------------------------------------------------------------- */
  /* write the results into the netCDF file                                  */

  for (io=0; io<Input.n_output_user; io++) {

    if (valid_result[io] == TRUE) {

      switch (Input.output_user[io]) {

      case OUTPUT_USER_WAVE:
        /* is done in each case, already done above, do nothing here */
      break;

      case OUTPUT_USER_ZOUT_SUR:
        switch (Input.atm.zout_source) {
        case OUTLEVEL_ZOUT_ABOVE_SUR:
        case OUTLEVEL_ZOUT_ABOVE_SEA:
        case OUTLEVEL_ATM_LEVELS:
        case OUTLEVEL_ALL_LEVELS:
          /* is done in each case, already done above, do nothing here  */

          /* if cold point tropopause is output level, than we like to save the height in the output */
          if (output->atm.microphys.z_cpt_sea != NOT_DEFINED_INTEGER) {
            tmp_float = output->atm.microphys.z_cpt_sea - output->alt.altitude;
            status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                            "z_cpt_sur", "cold point tropopause altitude above surface", "km", tmp_float, 
                                            id_dim_array3, index_3D, dim3);
            if (status!=0) {
              fprintf (stderr, " *** Error, while writing 'z_cpt_sur' into %s, %s (%s)\n", 
                               Input.filename[FN_OUTFILE], function_name, file_name );
            }
          }
          break;
        case OUTLEVEL_PRESS:
        case OUTLEVEL_MODEL_LEVELS:
	case OUTLEVEL_MODEL_LAYERS:
	case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "zout_sur", "altitude above surface", "km", output->atm.zout_sur, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'zout_sur' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          } 
          break;
        default:
	  fprintf (stderr, "Error, unknown zout_source %d (line %d, function '%s' in '%s')\n", Input.atm.zout_source, __LINE__, __func__, __FILE__);
          return -2;
        }
        break;

      case OUTPUT_USER_ZOUT_SEA:
        switch (Input.atm.zout_source) {
        case OUTLEVEL_ZOUT_ABOVE_SUR:
        case OUTLEVEL_ZOUT_ABOVE_SEA:
        case OUTLEVEL_ATM_LEVELS:
        case OUTLEVEL_ALL_LEVELS:
          /* is done in each case, already done above, do nothing here  */

          /* if cold point tropopause is output level, than we like to save the height in the output */
          if (output->atm.microphys.z_cpt_sea != NOT_DEFINED_INTEGER) {
            status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                            "z_cpt_sea", "cold point tropopause altitude above sealevel", "km", output->atm.microphys.z_cpt_sea, 
                                            id_dim_array3, index_3D, dim3);
            if (status!=0) {
              fprintf (stderr, " *** Error, while writing 'z_cpt_sea' into %s, %s (%s)\n", 
                               Input.filename[FN_OUTFILE], function_name, file_name );
            }
          }
          break;
        case OUTLEVEL_PRESS:
        case OUTLEVEL_MODEL_LEVELS:
	case OUTLEVEL_MODEL_LAYERS:
	case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "zout_sea", "altitude above sealevel", "km", output->atm.zout_sea, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'zout_sea' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          } 
          break;
        default:
	  fprintf (stderr, "Error, unknown zout_source %d (line %d, function '%s' in '%s')\n", Input.atm.zout_source, __LINE__, __func__, __FILE__);
          return -2;
        }
        break;

      case OUTPUT_USER_Z_SUR:
        status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "z_sur", "surface altitude", "km", output->alt.altitude, 
                                        id_dim_array3, index_3D, dim3);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'z_sur' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        } 
        break;

      case OUTPUT_USER_WAVENUMBER:
       /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->wl.nlambda_h, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          tmp_data_1D[iv] = 1.0E7/output->wl.lambda_h[iv];  /* convert from nm to cm-1 */

        status = write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "wavenumber", &(id_data[io]), "wavenumber", "cm-1", output->wl.lambda_h, 
                                        "wavenumber", &(id_dim[io]),  output->wl.nlambda_h);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'wavenumber' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }  
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_EDIR:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "edir", "direct irradiance", irradiance_units, output->rfldir, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'edir' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }  
        break;

      case OUTPUT_USER_EGLO:
        if ((status = ASCII_calloc_float (&tmp_data_3D, Input.n_output_user, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++)
            tmp_data_3D[lev][iv] = output->rfldir[lev][iv]+output->rfldn[lev][iv];
  
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "eglo", "global ( = direct + diffuse downward ) irradiance", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'eglo' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }  

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;
  
      case OUTPUT_USER_EDN:

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "edn", "diffuse downward irradiance", irradiance_units, output->rfldn, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'edn' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;
  
      case OUTPUT_USER_EUP:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "eup", "diffuse upward irradiance", irradiance_units, output->flup, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'eup' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_ENET:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->rfldir[lev][iv]+output->rfldn[lev][iv]-output->flup[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "enet", "net ( = global - upward ) irradiance", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'enet' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_ESUM:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->rfldir[lev][iv]+output->rfldn[lev][iv]+output->flup[lev][iv];
  
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "esum", "total ( = global + upward ) irradiance", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'esum' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_FDIR:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = 4.0*PI*output->uavgso[lev][iv];
  
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "fdir", "direct actinic flux", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'fdir' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_FGLO:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = 4.0*PI*(output->uavgso[lev][iv]+output->uavgdn[lev][iv]);

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "fglo", "global ( = direct + diffuse downward ) actinic flux", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'fglo' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_FDN:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = 4.0*PI*output->uavgdn[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "fdn", "diffuse downward actinic flux", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'fdn' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_FUP:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = 4.0*PI*output->uavgup[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "fup", "diffuse upward actinic flux", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'fup' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_F:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = 4.0*PI*output->uavg[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "f", "total actinic flux", irradiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'f' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_UDIR:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "uavgdir", "direct mean radiance ( = actinic flux/ (4 pi) )", 
                                        radiance_units, output->uavgso, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'uavgdir' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_UGLO:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->uavgso[lev][iv]+output->uavgdn[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "uavgglo", "global mean radiance ( = actinic flux/ (4 pi) )", 
                                        radiance_units, tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'uavgglo' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_UDN:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "uavgdn", "diffusive downward mean radiance ( = actinic flux/ (4 pi) )", 
                                        radiance_units, output->uavgdn, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'uavgdn' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_UUP:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "uavgup", "diffusive upward mean radiance ( = actinic flux/ (4 pi) )", 
                                        radiance_units, output->uavgup, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'uavgup' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_U:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "uavg", "diffusive upward mean radiance ( = actinic flux/ (4 pi) )", 
                                        radiance_units, output->uavg,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'uavg' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_ALB:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for eglo in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->flup[lev][iv]/(output->rfldir[lev][iv]+output->rfldn[lev][iv]);

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "albedo", "albedo", " ", tmp_data_3D, 
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'albedo' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_SZA:
        status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "SZA", "solar zenith angle", "degrees", output->sza_h[0], 
                                        id_dim_array3, index_3D, dim3);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'SZA' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_C_P:

       /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "c_p", "c_p", "J/(kg K)", output->atm.microphys.c_p, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_HEAT:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "heat", "heating rate", "K/day", output->heat,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'heat' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_EMIS:
        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "emis", "emission rate", "K/day", output->emis,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'emis' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_ABS:

        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for tmp_data in %s (%s)\n", function_name, file_name);
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->heat[lev][iv] - output->emis[lev][iv];

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "abs", "absorption rate", "K/day", tmp_data_3D,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'abs' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_W_RAD:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for tmp_data in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->w_zout[lev][iv]/3600.0; /* m/day -> m/s */

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "w_rad", "radiative vertical wind", "m/s", tmp_data_3D,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'w_rad' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_M_RAD:
        if ((status = ASCII_calloc_float (&tmp_data_3D, output->atm.nzout, output->wl.nlambda_h)) !=0 ) {
          fprintf (stderr,"Error, allocating memory for tmp_data in write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        for (iv=0;iv<output->wl.nlambda_h;iv++)
          for (lev=0;lev<output->atm.nzout;lev++) 
            tmp_data_3D[lev][iv] = output->w_zout[lev][iv] * output->atm.microphys.dens_zout[MOL_AIR][lev] 
                                   * 1.0e+6 * 1.0e-3 * Input.atm.mol_mass[MOL_AIR] / avogadro;   /* kg/ (m2 day) */

        status = write_netCDF_5D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "m_rad", "radiative vertical mass flux", "kg/(m2 day)", tmp_data_3D,
                                        id_dim_array5, index_5D, dim5, n_5D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'm_rad' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        ASCII_free_float(tmp_data_3D, output->wl.nlambda_h);
        break;

      case OUTPUT_USER_P:
        switch (Input.atm.zout_source) {
        case OUTLEVEL_ZOUT_ABOVE_SUR:
        case OUTLEVEL_ZOUT_ABOVE_SEA:
        case OUTLEVEL_ATM_LEVELS:
        case OUTLEVEL_ALL_LEVELS:
        case OUTLEVEL_MODEL_LEVELS:
	case OUTLEVEL_MODEL_LAYERS:
	case OUTLEVEL_MODEL_LEVELS_AND_LAYERS:
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "p", "pressure", "hPa", output->atm.microphys.press_zout, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'p' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          }
          break;
        case OUTLEVEL_PRESS:
          /* is done in each case, already done above in 1D, do nothing here  */
          /* if cold point tropopause is output level, than we like to save the height in the output */
          if (output->atm.microphys.press_cpt != NOT_DEFINED_INTEGER) {
            status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                            "p_cpt", "cold point tropopause pressure", "hPa", output->atm.microphys.press_cpt, 
                                            id_dim_array3, index_3D, dim3);
            if (status!=0) {
              fprintf (stderr, " *** Error, while writing 'p_cpt' into %s, %s (%s)\n", 
                               Input.filename[FN_OUTFILE], function_name, file_name );
            }
          }
          break;
        default:
	  fprintf (stderr, "Error, unknown zout_source %d (line %d, function '%s' in '%s')\n", Input.atm.zout_source, __LINE__, __func__, __FILE__);
          return -2;
        }
        break;

      case OUTPUT_USER_T:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "T", "temperature", "K", output->atm.microphys.temper_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'T' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_T_D:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "T_d", "dew point temperature", "K", output->atm.microphys.temper_d_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'T_d' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_T_SUR:
        status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "T_sur", "surface skin temperature", "K", output->surface_temperature, 
                                        id_dim_array3, index_3D, dim3);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'T_sur' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_THETA:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "theta", "potential temperature", "K", 
                                        output->atm.microphys.theta_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'theta' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_THETA_E:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "theta_e", "equivalent potential temperature", "K", 
                                        output->atm.microphys.theta_e_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'theta_e' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_CLWC: /* lwc[g/m3]*1000./rho_air[kg/m3] == clwc[kg/kg] */

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in %s (%s)\n", var_name, function_name, file_name);
          return -1;
        }

        /* calculate data */  /* 1.e-3 == convert g -> kg, devided by rho air */
        for (lev=0;lev<output->atm.nzout;lev++) 
          if (Input.i_wc!=-1)
            tmp_data_1D[lev] = output->caoth[Input.i_wc].microphys.lwc_zout[lev] * 1.e-3 / 
                             (output->atm.microphys.dens_zout[MOL_AIR][lev]*1.e+6*1.e-3*Input.atm.mol_mass[MOL_AIR]/avogadro) ;
	  else
	    tmp_data_1D[lev] = 0.0;

        /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "CLWC", "cloud liquid water content", "kg/kg", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_CLWD:
        if (Input.i_wc!=-1) {
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "CLWD", "cloud ice water density",    "g/m3", output->caoth[Input.i_wc].microphys.lwc_zout, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'CLWD' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
	  }
        }
        break;

      case OUTPUT_USER_TCLW:
        if (Input.i_wc!=-1) {
          status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "TCLW", "Total Column Liquid Water", "kg/m2", output->caoth[Input.i_wc].microphys.tcw, 
                                          id_dim_array3, index_3D, dim3);
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'TCLW' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          }
	}
        break;

      case OUTPUT_USER_REFF_WAT:
        if (Input.i_wc!=-1) {
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "reff_water", "effective radius liquid water cloud", "micro m", output->caoth[Input.i_wc].microphys.effr_zout, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'reff_wat' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          }
	}
        break;

      case OUTPUT_USER_CIWC: /* lwc[g/m3]*1000./rho_air[kg/m3] == clwc[kg/kg] */

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */  /* 1.e-3 == convert g -> kg, devided by rho air */
        for (lev=0;lev<output->atm.nzout;lev++) 
          if (Input.i_ic!=-1)
            tmp_data_1D[lev] = output->caoth[Input.i_ic].microphys.lwc_zout[lev] * 1.e-3 / 
	      (output->atm.microphys.dens_zout[MOL_AIR][lev]*1.e+6*1.e-3*Input.atm.mol_mass[MOL_AIR]/avogadro) ;
          else
            tmp_data_1D[lev] = 0.0;
        /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "CIWC", "cloud ice water content", "kg/kg", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_CIWD:
        if (Input.i_ic!=-1) {
          status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "CIWD", "cloud ice water density",    "g/m3", output->caoth[Input.i_ic].microphys.lwc_zout, 
                                          id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'CIWD' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          }
	}
        break;

      case OUTPUT_USER_TCIW:
        if (Input.i_ic!=-1) {
          status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                          "TCIW", "Total Column Ice Water", "kg/m2", output->caoth[Input.i_ic].microphys.tcw, 
                                          id_dim_array3, index_3D, dim3);
          if (status!=0) {
            fprintf (stderr, " *** Error, while writing 'TCIW' into %s, %s (%s)\n", 
                             Input.filename[FN_OUTFILE], function_name, file_name );
          }
	}
        break;

      case OUTPUT_USER_REFF_ICE:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        if (Input.i_ic!=-1) {
          if (Input.caoth[isp].properties==PROP_FU && Input.caoth[isp].fu2yang) 
	    for (lev=0; lev<output->atm.nzout; lev++)
	      tmp_data_1D[lev] = output->caoth[Input.i_ic].microphys.effr_zout[lev] * (3.0*sqrt(3.0)/4.0);
          else
            for (lev=0; lev<output->atm.nzout; lev++)
	      tmp_data_1D[lev] = output->caoth[Input.i_ic].microphys.effr_zout[lev];
        }
	else
	  for (lev=0; lev<output->atm.nzout; lev++)
	    tmp_data_1D[lev] = 0.0;

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "reff_ice", "effective radius ice water cloud",    "micro m", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'reff_ice' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        /* free memory */
        free(tmp_data_1D);

        break;

      case OUTPUT_USER_TCW:
        if (Input.i_wc!=-1 && Input.i_ic!=-1)
          status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE],
                                          "TCW", "Total Column Water", "kg/m2", output->caoth[Input.i_wc].microphys.tcw + output->caoth[Input.i_ic].microphys.tcw, 
                                          id_dim_array3, index_3D, dim3);
	else if (Input.i_wc!=-1)
          status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE],
                                          "TCW", "Total Column Water", "kg/m2", output->caoth[Input.i_wc].microphys.tcw,
                                          id_dim_array3, index_3D, dim3);
	else if (Input.i_ic!=-1)
          status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE],
                                          "TCW", "Total Column Water", "kg/m2", output->caoth[Input.i_ic].microphys.tcw, 
                                          id_dim_array3, index_3D, dim3);

        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'TCW' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_CC:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "CC", "cloud cover",    "0-1", output->cf.cf_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'CC' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_TCC:
        status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "TCC", "Total Cloud Cover", "0-1", output->cf.cf_total, 
                                        id_dim_array3, index_3D, dim3);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'TCC' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_CLOUDS:
        iscr=0;
	if (Input.i_wc!=-1)
	  iscr+=2;
	if (Input.i_ic!=-1)
	  iscr+=1;
	  
        status = write_netCDF_3D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "clouds", "Clouds: IC=1, WC=2, IC+WC=3", "0-3", iscr, 
                                        id_dim_array3, index_3D, dim3);
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'clouds' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;

      case OUTPUT_USER_N:
        i_gas = Input.output_user_gas[io];
        /* create variable name */
        gas = gas_number2string(i_gas);
        if (strcmp(" -1 ",gas)==0) {
          fprintf (stderr, "Error wrong gas identifier during write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        strcpy (var_name, "n_");
        strcat (var_name, gas);
        strcpy (strtrim(var_name," "),var_name); /* get rid of spaces before and behind text */

        strcpy (long_name, "number density of ");
        strcat (long_name, gas);
        strcpy (strtrim(long_name," "), long_name); /* get rid of spaces before and behind text */

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        var_name, long_name, "1/cm^3", output->atm.microphys.dens_zout[i_gas], 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }
        break;
  
      case OUTPUT_USER_RHO:   
        i_gas = Input.output_user_gas[io];
        /* create variable name */
        gas = gas_number2string(i_gas);
        if (strcmp(" -1 ",gas)==0) {
          fprintf (stderr, "Error wrong gas identifier during write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        strcpy (var_name, "rho_");
        strcat (var_name, gas);
        strcpy (strtrim(var_name," "),var_name); /* get rid of spaces before and behind text */

        strcpy (long_name, "mass density of ");
        strcat (long_name, gas);
        strcpy (strtrim(long_name," "), long_name); /* get rid of spaces before and behind text */

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */  /* 1.e+6: convert from cm-3 to m-3; 1.e-3: convert g -> kg */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->atm.microphys.dens_zout[i_gas][lev] 
                               * 1.e+6 * 1.e-3 * Input.atm.mol_mass[i_gas] / avogadro ;  

        /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        var_name, long_name, "kg/m^3", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_VMR:
        i_gas = Input.output_user_gas[io];
        /* create variable name */
        gas = gas_number2string(i_gas);
        if (strcmp(" -1 ",gas)==0) {
          fprintf (stderr, "Error wrong gas identifier during write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        strcpy (var_name, "vmr_");
        strcat (var_name, gas);
        strcpy (strtrim(var_name," "),var_name); /* get rid of spaces before and behind text */

        strcpy (long_name, "volume mixing ratio of ");
        strcat (long_name, gas);
        strcpy (strtrim(long_name," "), long_name); /* get rid of spaces before and behind text */

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->atm.microphys.dens_zout[i_gas]  [lev]/
                             output->atm.microphys.dens_zout[MOL_AIR][lev];  

        /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        var_name, long_name, " ", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_MMR:
        i_gas = Input.output_user_gas[io];
        /* create variable name */
        gas = gas_number2string(i_gas);
        if (strcmp(" -1 ",gas)==0) {
          fprintf (stderr, "Error wrong gas identifier during write_netCDF_output (in uvspec_lex.l)\n");
          return -1;
        }
        strcpy (var_name, "mmr_");
        strcat (var_name, gas);
        strcpy (strtrim(var_name," "),var_name); /* get rid of spaces */

        strcpy (long_name, "mass mixing ratio of ");
        strcat (long_name, gas);
        strcpy (strtrim(long_name," "), long_name); /* get rid of spaces before and behind text */

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = (output->atm.microphys.dens_zout[i_gas]  [lev] * Input.atm.mol_mass[i_gas]  ) /
                             (output->atm.microphys.dens_zout[MOL_AIR][lev] * Input.atm.mol_mass[MOL_AIR]);

        /* write data */
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        var_name, long_name, "kg/kg", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing '%s' into %s, %s (%s)\n", 
                           var_name, Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_RH:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->atm.microphys.dens_zout[MOL_H2O][lev]/
                             vapor_pressure(output->atm.microphys.temper_zout[lev]) * 100.0;

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "rh", "relative humidity over water", "per cent", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'rh' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_RH_ICE:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->atm.microphys.dens_zout[MOL_H2O][lev]/
                             vapor_pressure_over_ice(output->atm.microphys.temper_zout[lev]) * 100.0;

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "rh_ice", "relative humidity over ice", "per cent", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'rh_ice' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);
        break;

      case OUTPUT_USER_WIND_U:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "u", "zonal wind", "m/s", output->wind.u_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'u' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_WIND_V:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "v", "meridional wind", "m/s", output->wind.v_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'v' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_WIND_W:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "w", "vertical wind", "Pa/s", output->wind.w_zout, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'w' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_DTDX:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "dtheta_dx", "zonal gradiant of potential temperature", "K/m",    
                                        output->atm.microphys.dtheta_dx, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'dtheta_dx' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_DTDY:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "dtheta_dy", "merional gradiant of potential temperature", "K/m", 
                                        output->atm.microphys.dtheta_dy, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'dtheta_dy' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_DTDZ:
        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                       "dtheta_dz", "vertical gradiant of potential temperature", "K/m", 
                                        output->atm.microphys.dtheta_dz, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'dtheta_dz' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }
       break;

      case OUTPUT_USER_HEAT_AD_X:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->wind.u_zout[lev] * output->atm.microphys.dtheta_dx[lev] *24*60*60; /* K/day */

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "heat_ad_x", "zonal heat advection", "K/day", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'heat_ad_x' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);

       break;

     case OUTPUT_USER_HEAT_AD_Y:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->wind.v_zout[lev] * output->atm.microphys.dtheta_dy[lev] *24*60*60; /* K/day */

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "heat_ad_y", "meridional heat advection", "K/day", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'heat_ad_y' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);

       break;

     case OUTPUT_USER_HEAT_AD_Z:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = output->wind.w_zout[lev] * output->atm.microphys.dtheta_dz[lev] *24*60*60;

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "heat_ad_z", "vertical heat advection", "K/day", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'heat_ad_z' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);

       break;

     case OUTPUT_USER_HEAT_AD:

        /* callocate data */
        if ((tmp_data_1D  = (float *) calloc (output->atm.nzout, sizeof (float))) == NULL) {
          fprintf (stderr,"Error, allocating memory for %s in write_netCDF_output (in uvspec_lex.l)\n", var_name);
          return -1;
        }

        /* calculate data */
        for (lev=0;lev<output->atm.nzout;lev++) 
          tmp_data_1D[lev] = (output->wind.u_zout[lev] * output->atm.microphys.dtheta_dx[lev] + 
                              output->wind.v_zout[lev] * output->atm.microphys.dtheta_dy[lev] +
                              output->wind.w_zout[lev] * output->atm.microphys.dtheta_dz[lev]  ) *24*60*60;

        status = write_netCDF_4D_float (ncid, Input.filename[FN_OUTFILE], 
                                        "heat_ad", "heat advection", "K/day", tmp_data_1D, 
                                        id_dim_array4, index_4D, dim4, n_4D, output->atm.zout_index, &(simulation_done));
        if (status!=0) {
          fprintf (stderr, " *** Error, while writing 'heat_ad' into %s, %s (%s)\n", 
                           Input.filename[FN_OUTFILE], function_name, file_name );
        }

        /* free memory */
        free(tmp_data_1D);

       break;

      default:
        fprintf (stderr, "Error, unknown user output %d (line %d, function '%s' in '%s')\n", 
                         Input.output_user[io], __LINE__, __func__, __FILE__);
        return -1;
      }
    }
  }



  if (Input.rte.numu > 0 && Input.rte.nphi > 0) {
    /* write phi array */
    status = write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "phi", &id_phi, "radiance azimuth angle", "degrees", Input.rte.phi, 
                           "phi", &id_nphi, Input.rte.nphi );
    if (status!=0) {
      fprintf (stderr, " *** Error, while writing 'phi' into %s, %s (%s)\n", 
                       Input.filename[FN_OUTFILE], function_name, file_name );
    }

    /* write umu array */
    status = write_netCDF_1D_float (ncid, Input.filename[FN_OUTFILE], 
                           "umu", &id_umu, "radiance zenith angle", "degrees", Input.rte.umu, 
                           "umu", &id_numu, Input.rte.numu );
    if (status!=0) {
      fprintf (stderr, " *** Error, while writing 'umu' into %s, %s (%s)\n", 
                       Input.filename[FN_OUTFILE], function_name, file_name );
    }

    /* order of dimensions */
    id_dim_array6[POSITION_TIME_6D  ] = id_ntime;
    id_dim_array6[POSITION_LAMBDA_6D] = id_nlambda;
    id_dim_array6[POSITION_ZOUT_6D  ] = id_nzout;
    id_dim_array6[POSITION_UMU_6D   ] = id_numu;
    id_dim_array6[POSITION_LAT_6D   ] = id_ny;  /* one before last dimension is default y-axis in ncview */
    id_dim_array6[POSITION_LON_6D   ] = id_nx;  /*            last dimension is default x-axis in ncview */

    /* initialisation of the index, where to write */
    index_6D[POSITION_TIME_6D]   = itime;
    n_6D    [POSITION_LAMBDA_6D] = output->wl.nlambda_h; /* loop over lambda */    
    n_6D    [POSITION_ZOUT_6D]   = output->atm.nzout;      /* loop over zout   */ 
    n_6D    [POSITION_UMU_6D]    = Input.rte.numu;       /* loop over umu    */ 
    index_6D[POSITION_LAT_6D]    = iy;                   /* latitude or x-pixel */
    index_6D[POSITION_LON_6D]    = ix;                   /* longitude or y-pixel */

    /* order of dimensions */
    id_dim_array7[POSITION_TIME_7D  ] = id_ntime;
    id_dim_array7[POSITION_LAMBDA_7D] = id_nlambda;
    id_dim_array7[POSITION_ZOUT_7D  ] = id_nzout;
    id_dim_array7[POSITION_UMU_7D   ] = id_nphi;
    id_dim_array7[POSITION_PHI_7D   ] = id_numu;
    id_dim_array7[POSITION_LAT_7D   ] = id_ny;  /* one before last dimension is default y-axis in ncview */
    id_dim_array7[POSITION_LON_7D   ] = id_nx;  /*            last dimension is default x-axis in ncview */

    /* initialisation of the index, where to write */
    index_7D[POSITION_TIME_7D]   = itime;
    n_7D    [POSITION_LAMBDA_7D] = output->wl.nlambda_h; /* loop over lambda */    
    n_7D    [POSITION_ZOUT_7D]   = output->atm.nzout;      /* loop over zout   */ 
    n_7D    [POSITION_UMU_7D]    = Input.rte.numu;       /* loop over umu    */ 
    n_7D    [POSITION_PHI_7D]    = Input.rte.nphi;       /* loop over umu    */ 
    index_7D[POSITION_LAT_7D]    = iy;                   /* latitude  or x-pixel */
    index_7D[POSITION_LON_7D]    = ix;                   /* longitude or y-pixel */

    /* write horizontal averaged radiance u0u */
    status = write_netCDF_6D_float (ncid, Input.filename[FN_OUTFILE], 
                                    "u0u", "horizontally averaged radiance", radiance_units, output->u0u, 
                                    id_dim_array6, index_6D, 6, n_6D, output->atm.zout_index, &(simulation_done));
    if (status!=0) {
      fprintf (stderr, " *** Error, while writing 'u0u' into %s, %s (%s)\n", 
                       Input.filename[FN_OUTFILE], function_name, file_name );
    }

    /* write radiance */
    status = write_netCDF_7D_float (ncid, Input.filename[FN_OUTFILE], 
                                    "uu", "radiance", radiance_units, output->uu, 
                                    id_dim_array7, index_7D, 7, n_7D, output->atm.zout_index, &(simulation_done));
    if (status!=0) {
      fprintf (stderr, " *** Error, while writing 'uu' into %s, %s (%s)\n", 
                       Input.filename[FN_OUTFILE], function_name, file_name );
    }
  }

  if (Input.verbose)
    fprintf (stderr, " ... write simulation_done = %d (0=not_done,1=done) in netCDF file %s\n", simulation_done, Input.filename[FN_OUTFILE]);

  status += write_netCDF_3D_byte (ncid, Input.filename[FN_OUTFILE], 
                                 "done", "0=not jet calculated, 1=calculation done", " ", simulation_done, 
                                 id_dim_array3, index_3D, dim3);
  if (status!=0) {
    fprintf (stderr, " *** Error, while writing 'done' into %s, %s (%s)\n", 
                     Input.filename[FN_OUTFILE], function_name, file_name );
  }

  nc_close (ncid);
  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use the netCDF output option. Please get netcdf and rebuild.        *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}


int check_if_result_exists (int *output_user, int n_output_user, int solver, int **valid_result, int quiet)
{
  int status=0;
  int io=0;

  for ( io=0; io< n_output_user; io++ ) {
    switch (solver) {
    case SOLVER_TWOSTR:
    case SOLVER_FTWOSTR:
      /* no dir, dn up parts of the actinic flux and radiance */

      switch (Input.output_user[io]) {
      case OUTPUT_USER_FDIR:
      case OUTPUT_USER_FGLO:
      case OUTPUT_USER_FDN:
      case OUTPUT_USER_FUP:
      case OUTPUT_USER_UDIR:
      case OUTPUT_USER_UGLO:
      case OUTPUT_USER_UDN:
        (*valid_result)[io] = FALSE;
        if (!quiet) {
          fprintf (stderr, "*** Warning: twostream solver does not separate between upward and \n");
          fprintf (stderr, "***          downward actinic flux or radiances. THIS OUTPUT IS SKIPPED. \n");
        }
        break;
      case OUTPUT_USER_WAVE:
      case OUTPUT_USER_WAVENUMBER:
      case OUTPUT_USER_ZOUT_SUR:
      case OUTPUT_USER_ZOUT_SEA:
      case OUTPUT_USER_Z_SUR:
      case OUTPUT_USER_EDIR:
      case OUTPUT_USER_EGLO:
      case OUTPUT_USER_EDN:
      case OUTPUT_USER_EUP:
      case OUTPUT_USER_ENET:
      case OUTPUT_USER_ESUM:
      case OUTPUT_USER_F:
      case OUTPUT_USER_U:
      case OUTPUT_USER_ALB:
      case OUTPUT_USER_SZA:
      case OUTPUT_USER_C_P:
      case OUTPUT_USER_HEAT:
      case OUTPUT_USER_W_RAD:
      case OUTPUT_USER_M_RAD:
      case OUTPUT_USER_P:
      case OUTPUT_USER_T:
      case OUTPUT_USER_T_D:
      case OUTPUT_USER_THETA:
      case OUTPUT_USER_THETA_E:
      case OUTPUT_USER_N:
      case OUTPUT_USER_RHO:   
      case OUTPUT_USER_VMR:
      case OUTPUT_USER_MMR:
      case OUTPUT_USER_CLWC:
      case OUTPUT_USER_CLWD:
      case OUTPUT_USER_TCLW:
      case OUTPUT_USER_REFF_WAT: 
      case OUTPUT_USER_CIWC:
      case OUTPUT_USER_CIWD:
      case OUTPUT_USER_TCIW:
      case OUTPUT_USER_REFF_ICE:
      case OUTPUT_USER_TCW:
      case OUTPUT_USER_CC:
      case OUTPUT_USER_TCC:
      case OUTPUT_USER_CLOUDS:
      case OUTPUT_USER_RH:
      case OUTPUT_USER_RH_ICE:
      case OUTPUT_USER_WIND_U:
      case OUTPUT_USER_WIND_V:
      case OUTPUT_USER_WIND_W:
      case OUTPUT_USER_HEAT_AD_X:
      case OUTPUT_USER_HEAT_AD_Y:
      case OUTPUT_USER_HEAT_AD_Z:
      case OUTPUT_USER_HEAT_AD:
        (*valid_result)[io] = TRUE;
        break;
      default:
        fprintf (stderr, "Error, unknown user output %d (line %d, function '%s' in '%s')\n", 
                         Input.output_user[io], __LINE__, __func__, __FILE__);
        return -1;
      }
      break;      

    case SOLVER_FDISORT1:
    case SOLVER_FDISORT2:
    case SOLVER_DISORT:
    case SOLVER_SDISORT:
    case SOLVER_SPSDISORT:
    case SOLVER_MONTECARLO:
    case SOLVER_NULL:
    case SOLVER_SOS:
      /* everything possible */
      (*valid_result)[io] = TRUE;
      break;

    case SOLVER_POLRADTRAN:
      fprintf (stderr, "Error, netCDF output for POLRADTRAN not yet implemented (line %d, function '%s' in '%s')\n",
                        __LINE__, __func__, __FILE__);
      return -1;
      break;

    case SOLVER_SSS:
    case SOLVER_SSSI:
    case SOLVER_TZS:
      /* only radiance in a direction (no fluxes, no actinic fluxes, no averaged fluxes) */
      fprintf (stderr, "Error, netCDF for 'solar single scattering' and 'thermal_zero_scattering'\n");
      fprintf (stderr, "       not yet implemented (in uvspec_lex.l)\n");
      return -2;
      break;

    default:
      fprintf (stderr, "Error, unknown rte_solver %d (line %d, function '%s' in '%s')\n", solver, __LINE__, __func__, __FILE__);
      return -3;
    }
  }

  return status;
}




static void free_output(input_struct Input, output_struct *output) 
{

  if (output->down_flux != NULL)
    ASCII_free_float_3D (output->down_flux, 
	                 output->atm.nzout, Input.rte.polradtran[POLRADTRAN_NSTOKES]);

  if (output->down_flux_t != NULL)
    ASCII_free_float_3D (output->down_flux_t, 
			 output->atm.nzout, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->up_flux != NULL)
    ASCII_free_float_3D (output->up_flux, 
			 output->atm.nzout, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->up_flux_t != NULL)
    ASCII_free_float_3D (output->up_flux_t, 
			 output->atm.nzout, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->down_rad != NULL)
    ASCII_free_float_5D (output->down_rad, 
			 output->atm.nzout, Input.rte.nphi, 
			 Input.rte.numu, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->down_rad_t != NULL)
    ASCII_free_float_5D (output->down_rad_t, 
			 output->atm.nzout, Input.rte.nphi,
			 Input.rte.numu, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->up_rad != NULL)
    ASCII_free_float_5D (output->up_rad, 
			 output->atm.nzout, Input.rte.nphi,
			 Input.rte.numu, Input.rte.polradtran[POLRADTRAN_NSTOKES]);
  
  if (output->up_rad_t != NULL)
    ASCII_free_float_5D (output->up_rad_t, 
			 output->atm.nzout, Input.rte.nphi,
			 Input.rte.numu, Input.rte.polradtran[POLRADTRAN_NSTOKES]);

  for (isp=0; isp<Input.n_caoth; isp++) {
    if (output->mc.dt != NULL) {
      if (output->mc.dt[isp] != NULL)
        free(output->mc.dt[isp]);
      if (output->mc.om[isp] != NULL)
        free(output->mc.om[isp]);
      if (output->mc.g1[isp] != NULL)
        free(output->mc.g1[isp]);
      if (output->mc.g2[isp] != NULL)
        free(output->mc.g2[isp]);
      if (output->mc.ff[isp] != NULL)
        free(output->mc.ff[isp]);
      if (output->mc.ds[isp] != NULL)
        free(output->mc.ds[isp]);
      if (output->mc.re[isp] != NULL)
        free(output->mc.re[isp]);
    }
  }
  if (output->mc.dt != NULL)
    free(output->mc.dt);
  if (output->mc.om != NULL)
    free(output->mc.om);
  if (output->mc.g1 != NULL)
    free(output->mc.g1);
  if (output->mc.g2 != NULL)
    free(output->mc.g2);
  if (output->mc.ff != NULL)
    free(output->mc.ff);
  if (output->mc.ds != NULL)
    free(output->mc.ds);
  if (output->mc.re != NULL)
    free(output->mc.re);

  if (output->dtauc     != NULL) free(output->dtauc);
  if (output->dtauc_md  != NULL) free(output->dtauc_md);
  if (output->sza_h     != NULL) free(output->sza_h);
  if (output->sza_s     != NULL) free(output->sza_s);
  if (output->mu_values != NULL) free(output->mu_values);

}


static float yytext2float(char* yytext, int which_token, int max_token) {
  double dnumber=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL, *s=NULL;
  int count=0;
  char *tmpstring=NULL;
  
  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }
  
  if ((count != which_token) && (count>max_token))
    fprintf(stderr,"Wrong number of arguments on line %d\n", line_number+1);
  
  dnumber = atof(s);
  free(tmpstring);
  return (float) dnumber;
}


static double yytext2double(char* yytext, int which_token, int max_token) {
  double dnumber=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL, *s=NULL;
  int count=0;
  char *tmpstring=NULL;
  
  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);
  
  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }
  
  if ((count != which_token) && (count>max_token))
    fprintf(stderr,"Wrong number of arguments on line %d\n", line_number+1);
  
  dnumber = atof(s);

  free(tmpstring);
  return dnumber;
}


static float* yytext2floats(char* yytext, int nary) {
  int i=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL;
  char *tmpstring=NULL;
  
  /* Fill up ary */
  float *ary = (float *) calloc (nary, sizeof(float));
  
  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  i = 0;

  while ((p = strtok(NULL,token)) != NULL)
    ary[i++] = atof(p);

  free(tmpstring);
  return ary;
}


static int* yytext2integers(char* yytext, int nary) {
  int i=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL;
  char *tmpstring=NULL;
  
  /* Fill up ary */
  int *ary    = (int *) calloc (nary, sizeof(int));

  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  i = 0;

  while ((p = strtok(NULL,token)) != NULL)
    ary[i++] = atoi(p);

  free(tmpstring);
  return ary;
}


static int yytext2ntokens(char* yytext) {
  int nary=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL;
  char *tmpstring=NULL;
  
  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  /* Get number of ary elements */
  nary = 0;
  p = strtok(tmpstring,token);

  while ((p = strtok(NULL,token)) != NULL) 
    nary++;

  free(tmpstring);
  return nary;
}


static int yytext2int(char* yytext, int which_token, int max_token) {
  int inumber=0;
  char token[]=UVSPEC_TOKEN;
  char *p=NULL, *s=NULL;
  int count=0;
  char *tmpstring=NULL;
  
  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }

  if ((count != which_token) && (count>max_token))
    fprintf(stderr,"Wrong number of arguments on line %d\n", line_number+1);

  inumber = atoi(s);

  free(tmpstring);
  return inumber;
}


static char* yytext2string(char* yytext, int which_token, int max_token) {
  char token[]=UVSPEC_TOKEN;
  char *p=NULL, *s=NULL, *string=NULL, *tmpstring=NULL;
  int count=0;

  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  if (which_token == 1) {
    s = p;
  }
  else { 
    while ((p = strtok(NULL,token)) != NULL) {
      count++;
      if (count == which_token)
        s = p;
    }
  }

  if ((count != which_token) && (count > max_token))
    fprintf(stderr,"Wrong number of arguments on line %d (yytext2string) \n", line_number+1);

  string = (char *) calloc (strlen(s)+1, sizeof (char));
  strcpy (string, s);

  free(tmpstring);

  return string;
}


/* like yytext2string, but string is copied to a pre-allocated string dest */
static int yytextCstring(char *dest, char* yytext, int which_token, int max_token) {
  char token[]=UVSPEC_TOKEN;
  char *p=NULL, *s=NULL, *tmpstring=NULL;
  int count=0;

  tmpstring = (char *) calloc (strlen(yytext)+1, sizeof (char));
  strcpy (tmpstring, yytext);

  p = strtok(tmpstring,token);
  count = 1;
  while ((p = strtok(NULL,token)) != NULL) {
    count++;
    if (count == which_token)
      s = p;
  }

  if ((count != which_token) && (count>max_token))
    fprintf(stderr,"Wrong number of arguments on line %d\n", line_number+1);

  strcpy (dest, s);

  free(tmpstring);
  return 0;
}

int compare_ascending_float (void *ap, void *bp) { /* Sort such that umu is in ascending order */
  float *a = (float *) ap;
  float *b = (float *) bp;
  
  if (*a == *b) 
    return 0;

  if (*a > *b) 
    return 1;
  else 
    return -1;
}

static float* add_cmu_to_umu(int nstr, int numu, int maxumu, float* oldumu, int* cmuind, int* umuind) {
  
  int nn, i, k=0;
  float *cmu = (float *) calloc (nstr/2, sizeof(float));
  float *cwt = (float *) calloc (nstr/2, sizeof(float));
  float *umu = NULL;
  int not = 0;

  umu = (float *) calloc (nstr+numu, sizeof(float));

  nn = nstr/2;

  F77_FUNC (qgausn, QGAUSN)(&nn, cmu, cwt); 

  for (i=0;i<nn;i++){
    umu[i]=-cmu[i];
    umu[i+nn]=cmu[i];
  }
  for (i=0;i<numu;i++)
    umu[nstr+i]=oldumu[i];

  free(cwt);
  free(cmu);

  /* Sort umu in increasing order */
  qsort ( umu, maxumu, sizeof(float),  QSORT_CAST compare_ascending_float);

  k = 0;
  for (j=0;j<numu;j++) {
    for (i=0;i<maxumu;i++) {
      if (fabs(umu[i]-oldumu[j]) < FLT_MIN ) {
	umuind[k++] = i;   /* These are the indices for the user angles */
      }
    }
  }

  free(oldumu);

  /* Get the indics for the computational angles */
  k=0;
  i=0;
  for (j=0;j<maxumu;j++) {
    not = 0;
    for (i=0;i<numu;i++) {
      if ( j == umuind[i] ) {	not = 1; }
      

    }
    if ( ! not ) 
      cmuind[k++] = j; 
  }
    
  return umu;
}


