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


ECHAM_sza{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_SZA], yytext, 2, 2);
	Input.atm.sza = NOT_DEFINED_FLOAT;
	Input.atm.sza_source = SZA_ECHAM;
}
ECMWF_ic_file{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	Input.caoth[isp].source = CAOTH_FROM_ECMWF;
	yytextCstring ( Input.caoth[isp].filename, yytext, 2, 2);
}
ECMWF_levels_only {
	Input.atm.rs_add_upper_levels = FALSE;
}
ECMWF_ozone_climatology {
	Input.atm.ECMWF_ozone_climatology = TRUE;
}
ECMWF_wc_file{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"wc",1);
	Input.caoth[isp].source = CAOTH_FROM_ECMWF;
	yytextCstring ( Input.caoth[isp].filename, yytext, 2, 2);
}
ECMWF_wind_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_ECMWF_WIND_MAP], yytext, 2, 2);
}
aerosol_angstrom{SPACE_SIGNED_FLOAT}{SPACE_FLOAT} {
	Input.aer.alpha = yytext2float(yytext, 2, 3);
	Input.aer.beta = yytext2float(yytext, 3, 3);
	Input.aer.spec = 1;
}
aerosol_default {
	Input.aer.standard = 1;
}
aerosol_file{SPACE_WORD}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("explicit", s)==0 )	{ id = FN_AER_EXPLICIT; }
	else if ( strcasecmp("gg", s)==0 )	{ id = FN_AER_GG; }
	else if ( strcasecmp("moments", s)==0 )	{ id = FN_AER_MOMENTS; }
	else if ( strcasecmp("ssa", s)==0 )	{ id = FN_AER_SSA; }
	else if ( strcasecmp("tau", s)==0 )	{ id = FN_AER_TAU; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.aer.filename[id], yytext, 3, 3);
	Input.aer.spec = 1;
}
aerosol_haze{SPACE_INTEGER} {
	Input.aer.haze = yytext2int(yytext, 2, 2);
}
aerosol_king_byrne{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT} {
	Input.aer.alpha_0 = yytext2float(yytext, 2, 4);
	Input.aer.alpha_1 = yytext2float(yytext, 3, 4);
	Input.aer.alpha_2 = yytext2float(yytext, 4, 4);
	Input.aer.spec = 1;
}
aerosol_modify{SPACE_WORD}{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 4);
	if ( strcasecmp("gg", s)==0 )	{ id1 = MODIFY_VAR_GG; }
	else if ( strcasecmp("ssa", s)==0 )	{ id1 = MODIFY_VAR_SSA; }
	else if ( strcasecmp("tau", s)==0 )	{ id1 = MODIFY_VAR_TAU; }
	else if ( strcasecmp("tau550", s)==0 )	{ id1 = MODIFY_VAR_TAU550; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 4);
	if ( strcasecmp("scale", s)==0 )	{ id2 = MODIFY_TYPE_SCALE; }
	else if ( strcasecmp("set", s)==0 )	{ id2 = MODIFY_TYPE_SET; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.aer.modify[id1][id2] = yytext2float(yytext, 4, 4);
	Input.aer.spec = 1;
}
aerosol_profile_modtran {
	Input.aer.profile_modtran = TRUE;
}
aerosol_refrac_index{SPACE_SIGNED_FLOAT}{SPACE_FLOAT} {
	Input.aer.nreal = yytext2float(yytext, 2, 3);
	Input.aer.nimag = yytext2float(yytext, 3, 3);
	Input.aer.spec = 1;
}
aerosol_season{SPACE_INTEGER} {
	Input.aer.seasn = yytext2int(yytext, 2, 2);
	Input.aer.spec = 1;
}
aerosol_set_tau_at_wvl{SPACE_FLOAT}{SPACE_FLOAT} {
	Input.aer.tau_wvl_lambda = yytext2float(yytext, 2, 3);
	Input.aer.tau_wvl_tau = yytext2float(yytext, 3, 3);
}
aerosol_sizedist_file{SPACE_WORD} {
	yytextCstring ( Input.aer.filename[Id], yytext, 2, 2);
}
aerosol_species_file{SPACE_WORD} {
	Input.aer.mixture_name = yytext2string(yytext, 2, 2);
	Input.aer.n_species = -1;
	Input.aer.spec = 1;
}
aerosol_species_library{SPACE_WORD} {
	yytextCstring ( Input.aer.filename[FN_AER_SPECIES_LIB], yytext, 2, 2);
}
aerosol_visibility{SPACE_FLOAT} {
	Input.aer.visibility = yytext2float(yytext, 2, 2);
}
aerosol_vulcan{SPACE_INTEGER} {
	Input.aer.vulcan = yytext2int(yytext, 2, 2);
}
albedo{SPACE_FLOAT} {
	Input.alb.albedo = yytext2float(yytext, 2, 2);
	Input.alb.source = ALBEDO_CONSTANT;
}
albedo_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_ALBEDO], yytext, 2, 2);
	Input.alb.source = ALBEDO_FROM_ALBEDO_FILE;
	Input.alb.surf_type_map = FALSE;
	Input.alb.albedo = -1.0;
}
albedo_library{SPACE_WORD} {
	Input.alb.library = yytext2string(yytext, 2, 2);
}
albedo_map{SPACE_WORD}{SPACE_WORD}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.alb.source = ALBEDO_FROM_ALBEDO_MAP;
	yytextCstring ( Input.filename[FN_ALBEDO_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.alb.netCDF_alb_name, yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.alb.scale_factor = yytext2float(yytext, 4, ntokens+1);
	}
}
altitude{SPACE_SIGNED_FLOAT}{SPACE_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.alt.source = ALT_DIRECT_INPUT;
	Input.alt.altitude = yytext2float(yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		Input.alt.altitude_dz = yytext2float(yytext, 3, ntokens+1);
	}
}
altitude_map{SPACE_WORD}{SPACE_WORD}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.alt.source = ALT_FROM_MAP;
	yytextCstring ( Input.filename[FN_ALTITUDE_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.alt.netCDF_alt_name, yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.alt.scale_factor = yytext2float(yytext, 4, ntokens+1);
	}
}
atm_z_grid{SIGNED_FLOATS} {
	ntokens = yytext2ntokens(yytext);
	Input.atm.z_atm_forced_sea =  yytext2floats(yytext, ntokens);
	Input.atm.nz_atm_forced_sea = ntokens;
}
atmos_region{SPACE_INTEGER}{SPACE_INTEGER}{SPACE_INTEGER}{SPACE_INTEGER} {
	Input.rte.mc.backward.yes = 1;
	Input.rte.mc.ixmin = yytext2int(yytext, 2, 5);
	Input.rte.mc.iymin = yytext2int(yytext, 3, 5);
	Input.rte.mc.ixmax = yytext2int(yytext, 4, 5);
	Input.rte.mc.iymax = yytext2int(yytext, 5, 5);
}
atmosphere_file{SPACE_WORD} {
	Input.atmosphere_filename = yytext2string(yytext, 2, 2);
}
bpdf_tsang_u10{SPACE_SIGNED_FLOAT} {
	Input.bpdf.u10 = yytext2float(yytext, 2, 2);
	Input.bpdf.type = BPDF_TSANG;
}
brdf_ambrals{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	Input.rossli.source = ROSSLI_AMBRALS_CONSTANT;
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("geo", s)==0 )	{ id = BRDF_ROSSLI_GEO; }
	else if ( strcasecmp("iso", s)==0 )	{ id = BRDF_ROSSLI_ISO; }
	else if ( strcasecmp("vol", s)==0 )	{ id = BRDF_ROSSLI_VOL; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rossli.rossli[id] = yytext2float(yytext, 3, 3);
}
brdf_ambrals_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_ROSSLI], yytext, 2, 2);
	Input.rossli.source = ROSSLI_FROM_AMBRALS_FILE;
}
brdf_ambrals_hotspot {
	Input.rossli.hotspot = BRDF_ROSSLI_HOTSPOT_ON;
}
brdf_cam{SPACE_WORD}{SPACE_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("pcl", s)==0 )	{ id = BRDF_CAM_PCL; }
	else if ( strcasecmp("sal", s)==0 )	{ id = BRDF_CAM_SAL; }
	else if ( strcasecmp("u10", s)==0 )	{ id = BRDF_CAM_U10; }
	else if ( strcasecmp("uphi", s)==0 )	{ id = BRDF_CAM_UPHI; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.cm.param[id] = yytext2float(yytext, 3, 3);
	Input.disort2_brdf = BRDF_CAM;
}
brdf_cam_solar_wind {
	Input.disort2_brdf = BRDF_CAM;
	Input.cm.solar_wind = 1;
}
brdf_hapke{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("b0", s)==0 )	{ id = BRDF_HAPKE_B0; }
	else if ( strcasecmp("h", s)==0 )	{ id = BRDF_HAPKE_H; }
	else if ( strcasecmp("w", s)==0 )	{ id = BRDF_HAPKE_W; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.hapke.hapke[id] = yytext2float(yytext, 3, 3);
}
brdf_hapke_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_HAPKE], yytext, 2, 2);
	Input.hapke.source = HAPKE_FROM_HAPKE_FILE;
}
brdf_rossli{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("geo", s)==0 )	{ id = BRDF_ROSSLI_GEO; }
	else if ( strcasecmp("iso", s)==0 )	{ id = BRDF_ROSSLI_ISO; }
	else if ( strcasecmp("vol", s)==0 )	{ id = BRDF_ROSSLI_VOL; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rossli.rossli[id] = yytext2float(yytext, 3, 3);
}
brdf_rossli_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_ROSSLI], yytext, 2, 2);
	Input.rossli.source = ROSSLI_FROM_ROSSLI_FILE;
}
brdf_rossli_hotspot {
	Input.rossli.hotspot = BRDF_ROSSLI_HOTSPOT_ON;
}
brdf_rpv{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("k", s)==0 )	{ id = BRDF_RPV_K; }
	else if ( strcasecmp("rho0", s)==0 )	{ id = BRDF_RPV_RHO0; }
	else if ( strcasecmp("scale", s)==0 )	{ id = BRDF_RPV_SCALE; }
	else if ( strcasecmp("sigma", s)==0 )	{ id = BRDF_RPV_SIGMA; }
	else if ( strcasecmp("t1", s)==0 )	{ id = BRDF_RPV_T1; }
	else if ( strcasecmp("t2", s)==0 )	{ id = BRDF_RPV_T2; }
	else if ( strcasecmp("theta", s)==0 )	{ id = BRDF_RPV_THETA; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rpv.rpv[id] = yytext2float(yytext, 3, 3);
}
brdf_rpv_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_RPV], yytext, 2, 2);
	Input.rpv.source = RPV_FROM_RPV_FILE;
}
brdf_rpv_library{SPACE_WORD} {
	Input.rpv.library = yytext2string(yytext, 2, 2);
}
brdf_rpv_type{SPACE_INTEGER} {
	Input.alb.surface = yytext2int(yytext, 2, 2);
	Input.alb.surf_type_map = FALSE;
}
ck_fu_h2o_continuum{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("off", s)==0 )	{ Input.ck_h2ocont = CK_H2OCONT_OFF; }
	else if ( strcasecmp("on", s)==0 )	{ Input.ck_h2ocont = CK_H2OCONT_ON; }
	else if ( strcasecmp("v2.1", s)==0 )	{ Input.ck_h2ocont = CK_H2OCONT_V2_1; }
	else if ( strcasecmp("v2.4", s)==0 )	{ Input.ck_h2ocont = CK_H2OCONT_V2_4; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
ck_lowtran_absorption{SPACE_WORD}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("CO", s)==0 )	{ mol_id = CK_ABS_CO; }
	else if ( strcasecmp("HNO3", s)==0 )	{ mol_id = CK_ABS_HNO3; }
	else if ( strcasecmp("N2", s)==0 )	{ mol_id = CK_ABS_N2; }
	else if ( strcasecmp("NH3", s)==0 )	{ mol_id = CK_ABS_NH3; }
	else if ( strcasecmp("NO", s)==0 )	{ mol_id = CK_ABS_NO; }
	else if ( strcasecmp("O4", s)==0 )	{ mol_id = CK_ABS_O4; }
	else if ( strcasecmp("SO2", s)==0 )	{ mol_id = CK_ABS_SO2; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 3);
	if ( strcasecmp("off", s)==0 )	{ Input.ck_abs[mol_id] = SWITCH_OFF; }
	else if ( strcasecmp("on", s)==0 )	{ Input.ck_abs[mol_id] = SWITCH_ON; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
cloud_fraction_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_CLOUD_FRACTION], yytext, 2, 2);
}
cloud_fraction_map{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_CLOUD_FRACTION_MAP], yytext, 2, 2);
}
cloud_overlap{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("max", s)==0 )	{ Input.cloud_overlap = CLOUD_OVERLAP_MAX; }
	else if ( strcasecmp("maxrand", s)==0 )	{ Input.cloud_overlap = CLOUD_OVERLAP_MAXRAND; }
	else if ( strcasecmp("off", s)==0 )	{ Input.cloud_overlap = CLOUD_OVERLAP_OFF; }
	else if ( strcasecmp("rand", s)==0 )	{ Input.cloud_overlap = CLOUD_OVERLAP_RAND; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
cloudcover{SPACE_WORD}{SPACE_FLOAT} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, 2, 3),0);
	Input.caoth[isp].cloudcover = yytext2float(yytext, 3, 3);
	Input.caoth[isp].ipa = 1;
	Input.ipa = 1;
}
cox_and_munk_pcl_map{SPACE_WORD}{SPACE_WORD}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.disort2_brdf = BRDF_CAM;
	yytextCstring ( Input.filename[FN_PIGMENTS_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.cm.pcl_netCDF_name, yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.cm.pcl_scale_factor = yytext2float(yytext, 4, ntokens+1);
	}
}
cox_and_munk_sal_map{SPACE_WORD}{SPACE_WORD}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.disort2_brdf = BRDF_CAM;
	yytextCstring ( Input.filename[FN_SALINITY_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.cm.sal_netCDF_name, yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.cm.sal_scale_factor = yytext2float(yytext, 4, ntokens+1);
	}
}
cox_and_munk_u10_map{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_U10_MAP], yytext, 2, 2);
	Input.disort2_brdf = BRDF_CAM;
}
crs_file{SPACE_WORD}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("BRO", s)==0 )	{ mol_id = MOL_BRO; }
	else if ( strcasecmp("CH4", s)==0 )	{ mol_id = MOL_CH4; }
	else if ( strcasecmp("CO", s)==0 )	{ mol_id = MOL_CO; }
	else if ( strcasecmp("CO2", s)==0 )	{ mol_id = MOL_CO2; }
	else if ( strcasecmp("H2O", s)==0 )	{ mol_id = MOL_H2O; }
	else if ( strcasecmp("HCHO", s)==0 )	{ mol_id = MOL_HCHO; }
	else if ( strcasecmp("N2", s)==0 )	{ mol_id = MOL_N2; }
	else if ( strcasecmp("N2O", s)==0 )	{ mol_id = MOL_N2O; }
	else if ( strcasecmp("NO2", s)==0 )	{ mol_id = MOL_NO2; }
	else if ( strcasecmp("O2", s)==0 )	{ mol_id = MOL_O2; }
	else if ( strcasecmp("O3", s)==0 )	{ mol_id = MOL_O3; }
	else if ( strcasecmp("O4", s)==0 )	{ mol_id = MOL_O4; }
	else if ( strcasecmp("OCLO", s)==0 )	{ mol_id = MOL_OCLO; }
	else if ( strcasecmp("SO2", s)==0 )	{ mol_id = MOL_SO2; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Output.crs.filename[mol_id], yytext, 3, 3);
}
crs_model{SPACE_WORD}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("no2", s)==0 )	{ mol_id = CRS_MOL_NO2; }
	else if ( strcasecmp("o3", s)==0 )	{ mol_id = CRS_MOL_O3; }
	else if ( strcasecmp("rayleigh", s)==0 )	{ mol_id = CRS_MOL_RAYLEIGH; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 3);
	if ( strcasecmp("Bass_and_Paur", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_BASS_AND_PAUR; }
	else if ( strcasecmp("Bodhaine", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_BODHAINE; }
	else if ( strcasecmp("Bodhaine29", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_BODHAINE29; }
	else if ( strcasecmp("Bogumil", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_BOGUMIL; }
	else if ( strcasecmp("Burrows", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_BURROWS; }
	else if ( strcasecmp("Daumont", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_DAUMONT; }
	else if ( strcasecmp("Molina", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_MOLINA; }
	else if ( strcasecmp("Nicolet", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_NICOLET; }
	else if ( strcasecmp("Penndorf", s)==0 )	{ Input.crs_model[mol_id] = CRS_MODEL_PENNDORF; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
data_files_path{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_PATH], yytext, 2, 2);
}
day_of_year{SPACE_INTEGER} {
	Input.UTC.tm_yday = yytext2int(yytext, 2, 2);
}
deltam{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("off", s)==0 )	{ Input.rte.deltam = SWITCH_OFF; }
	else if ( strcasecmp("on", s)==0 )	{ Input.rte.deltam = SWITCH_ON; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
disort_intcor{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("moments", s)==0 )	{ Input.rte.disort_icm = DISORT_ICM_MOMENTS; }
	else if ( strcasecmp("off", s)==0 )	{ Input.rte.disort_icm = DISORT_ICM_OFF; }
	else if ( strcasecmp("phase", s)==0 )	{ Input.rte.disort_icm = DISORT_ICM_PHASE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
disort_spherical_albedo {
	Input.rte.ibcnd = 1;
	Input.rte.nphi = 0;
}
earth_radius{SPACE_SIGNED_FLOAT} {
	Input.r_earth = yytext2float(yytext, 2, 2);
}
emissivity_map{SPACE_WORD}{SPACE_WORD}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.alb.source = ALBEDO_FROM_EMISSIVITY_MAP;
	yytextCstring ( Input.filename[FN_ALBEDO_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.alb.netCDF_alb_name, yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.alb.scale_factor = yytext2float(yytext, 4, ntokens+1);
	}
}
filter_function_file{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	yytextCstring ( Input.filename[FN_FILTERFUNCTION], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("normalize", s)==0 )	{ Input.filter_function_normalize = 1; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
fluorescence{SPACE_SIGNED_FLOAT} {
	Input.flu.fluorescence = yytext2float(yytext, 2, 2);
	Input.flu.source = FLUORESCENCE_CONSTANT;
}
fluorescence_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_FLUORESCENCE], yytext, 2, 2);
	Input.flu.source = FLUORESCENCE_FROM_FLUORESCENCE_FILE;
}
heating_rate{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	Input.heating = HEAT_LAYER_CD;
	if ( ntokens >= 1 )	 {
		s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("layer_cd", s)==0 )	{ Input.heating = HEAT_LAYER_CD; }
	else if ( strcasecmp("layer_fd", s)==0 )	{ Input.heating = HEAT_LAYER_FD; }
	else if ( strcasecmp("local", s)==0 )	{ Input.heating = HEAT_LOCAL; }
	else if ( strcasecmp("none", s)==0 )	{ Input.heating = HEAT_NONE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
ic_file{SPACE_option_definition.ProfileType}{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("1d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_1D; }
	else if ( strcasecmp("3d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_3D; }
	else if ( strcasecmp("ipa_files", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_IPA_FILES; }
	else if ( strcasecmp("moments", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_MOMENTS; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.caoth[isp].filename, yytext, 3, 3);
}
ic_fu{SPACE_WORD}{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("deltascaling", s)==0 )	{ id1 = IC_FU_DELTASCALING; }
	else if ( strcasecmp("reff_def", s)==0 )	{ id1 = IC_FU_REFF_DEF; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 3);
	if ( strcasecmp("off", s)==0 )	{ id2 = SWITCH_OFF; }
	else if ( strcasecmp("on", s)==0 )	{ id2 = SWITCH_ON; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.caoth[isp].ic_fu[id1] = id2;
}
ic_habit{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("dendrite", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_DENDRITE; }
	else if ( strcasecmp("droxtal", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_DROXTAL; }
	else if ( strcasecmp("ghm", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_GHM; }
	else if ( strcasecmp("hollow-column", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_HOLLOW_COLUMN; }
	else if ( strcasecmp("plate", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_PLATE; }
	else if ( strcasecmp("rosette-4", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_ROSETTE_4; }
	else if ( strcasecmp("rosette-6", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_ROSETTE_6; }
	else if ( strcasecmp("rough-aggregate", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_ROUGH_AGGREGATE; }
	else if ( strcasecmp("solid-column", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_SOLID_COLUMN; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
ic_habit_yang2013{SPACE_WORD}{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("column_8elements", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_COLUMN_8ELEMENTS; }
	else if ( strcasecmp("droxtal", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_DROXTAL; }
	else if ( strcasecmp("hollow_bullet_rosette", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_HOLLOW_BULLET_ROSETTE; }
	else if ( strcasecmp("hollow_column", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_HOLLOW_COLUMN; }
	else if ( strcasecmp("plate", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_PLATE; }
	else if ( strcasecmp("plate_10elements", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_PLATE_10ELEMENTS; }
	else if ( strcasecmp("plate_5elements", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_PLATE_5ELEMENTS; }
	else if ( strcasecmp("solid_bullet_rosette", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_SOLID_BULLET_ROSETTE; }
	else if ( strcasecmp("solid_column", s)==0 )	{ Input.caoth[isp].habit = IC_HABIT_SOLID_COLUMN; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 3);
	if ( strcasecmp("moderate", s)==0 )	{ Input.caoth[isp].roughness = IC_ROUGHNESS_MODERATE; }
	else if ( strcasecmp("severe", s)==0 )	{ Input.caoth[isp].roughness = IC_ROUGHNESS_SEVERE; }
	else if ( strcasecmp("smooth", s)==0 )	{ Input.caoth[isp].roughness = IC_ROUGHNESS_SMOOTH; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
ic_modify{SPACE_WORD}{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, 4);
	if ( strcasecmp("gg", s)==0 )	{ id1 = MODIFY_VAR_GG; }
	else if ( strcasecmp("ssa", s)==0 )	{ id1 = MODIFY_VAR_SSA; }
	else if ( strcasecmp("tau", s)==0 )	{ id1 = MODIFY_VAR_TAU; }
	else if ( strcasecmp("tau550", s)==0 )	{ id1 = MODIFY_VAR_TAU550; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 4);
	if ( strcasecmp("scale", s)==0 )	{ id2 = MODIFY_TYPE_SCALE; }
	else if ( strcasecmp("set", s)==0 )	{ id2 = MODIFY_TYPE_SET; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.caoth[isp].modify[id1][id2] = yytext2float(yytext, 4, 4);
}
ic_properties{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"ic",1);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("baum", s)==0 )	{ Input.caoth[isp].properties = PROP_BAUM; }
	else if ( strcasecmp("baum_v36", s)==0 )	{ Input.caoth[isp].properties = PROP_BAUM_V36; }
	else if ( strcasecmp("echam4", s)==0 )	{ Input.caoth[isp].properties = PROP_ECHAM4; }
	else if ( strcasecmp("fu", s)==0 )	{ Input.caoth[isp].properties = PROP_FU; }
	else if ( strcasecmp("hey", s)==0 )	{ Input.caoth[isp].properties = PROP_HEY; }
	else if ( strcasecmp("key", s)==0 )	{ Input.caoth[isp].properties = PROP_KEY; }
	else if ( strcasecmp("yang", s)==0 )	{ Input.caoth[isp].properties = PROP_YANG; }
	else if ( strcasecmp("yang2013", s)==0 )	{ Input.caoth[isp].properties = PROP_YANG2013; }
	else	{ Input.caoth[isp].properties = PROP_FILE;
		 strcpy ( Input.caoth[isp].properties_filename, s ); }
	free(s); 
	if ( ntokens >= 2 )	 {
		s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("interpolate", s)==0 )	{ Input.caoth[isp].interpolate = 1; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
interpret_as_level{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, 2, 2),0);
	Input.caoth[isp].layer = FALSE;
}
isotropic_source_toa {
	Input.rte.fisot = 1;
}
latitude{SPACE_WORD}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("n", s)==0 )	{ Input.lat_signum = LATITUDE_SIGNUM_N; }
	else if ( strcasecmp("s", s)==0 )	{ Input.lat_signum = LATITUDE_SIGNUM_S; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.lat_degrees = yytext2float(yytext, 3, ntokens+1);
	if ( ntokens >= 3 )	 {
		Input.lat_minutes = yytext2float(yytext, 4, ntokens+1);
	}
	if ( ntokens >= 4 )	 {
		Input.lat_seconds = yytext2float(yytext, 5, ntokens+1);
	}
}
longitude{SPACE_WORD}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}?{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("e", s)==0 )	{ Input.lon_signum = LONGITUDE_SIGNUM_E; }
	else if ( strcasecmp("w", s)==0 )	{ Input.lon_signum = LONGITUDE_SIGNUM_W; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.lon_degrees = yytext2float(yytext, 3, ntokens+1);
	if ( ntokens >= 3 )	 {
		Input.lon_minutes = yytext2float(yytext, 4, ntokens+1);
	}
	if ( ntokens >= 4 )	 {
		Input.lon_seconds = yytext2float(yytext, 5, ntokens+1);
	}
}
mc_aerosol_is{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_AERIS], yytext, 2, 2);
	Input.rte.mc.concentration_is = TRUE;
}
mc_albedo_spectral {
	Input.rte.mc.spectral = 1;
}
mc_albedo_spectral_file{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_ALBEDO_SPECTRAL], yytext, 2, 2);
}
mc_azimuth_old {
	Input.rte.mc.azimuth = MCAZIMUTH_OLD;
}
mc_backward{SPACE_INTEGER}?{SPACE_INTEGER}?{SPACE_INTEGER}?{SPACE_INTEGER}? {
	ntokens = yytext2ntokens(yytext);
	Input.rte.mc.backward.yes = 1;
	if ( ntokens >= 1 )	 {
		Input.rte.mc.backward.islower = yytext2int(yytext, 2, ntokens+1);
	}
	if ( ntokens >= 2 )	 {
		Input.rte.mc.backward.jslower = yytext2int(yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		Input.rte.mc.backward.isupper = yytext2int(yytext, 4, ntokens+1);
	}
	if ( ntokens >= 4 )	 {
		Input.rte.mc.backward.jsupper = yytext2int(yytext, 5, ntokens+1);
	}
}
mc_backward_output{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("abs", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_ABS; }
	else if ( strcasecmp("act", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_ACT; }
	else if ( strcasecmp("edir", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_EDIR; }
	else if ( strcasecmp("edn", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_EDN; }
	else if ( strcasecmp("emis", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_EMIS; }
	else if ( strcasecmp("eup", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_EUP; }
	else if ( strcasecmp("fdir", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_FDIR; }
	else if ( strcasecmp("fdn", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_FDN; }
	else if ( strcasecmp("fup", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_FUP; }
	else if ( strcasecmp("heat", s)==0 )	{ Input.rte.mc.backward.output = MCBACKWARD_HEAT; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	if ( ntokens >= 2 )	 {
		s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("K_per_day", s)==0 )	{ Input.rte.mc.abs_unit = MCABS_UNIT_K_PER_DAY; }
	else if ( strcasecmp("W_per_m2_and_dz", s)==0 )	{ Input.rte.mc.abs_unit = MCABS_UNIT_W_PER_M2_AND_DZ; }
	else if ( strcasecmp("W_per_m3", s)==0 )	{ Input.rte.mc.abs_unit = MCABS_UNIT_W_PER_M3; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
mc_backward_sunshape_file{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_SUNSHAPE_FILE], yytext, 2, 2);
}
mc_backward_writeback {
	Input.rte.mc.backward.writeback = 1;
}
mc_basename{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_BASENAME], yytext, 2, 2);
}
mc_bcond{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("absorb", s)==0 )	{ Input.rte.mc.bcond = MCBCOND_ABSORB; }
	else if ( strcasecmp("mirror", s)==0 )	{ Input.rte.mc.bcond = MCBCOND_MIRROR; }
	else if ( strcasecmp("periodic", s)==0 )	{ Input.rte.mc.bcond = MCBCOND_PERIODIC; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
mc_coherent_backscatter {
	Input.rte.mc.coherent_backscatter = 1;
}
mc_ipa {
	Input.rte.mc.ipa = 1;
}
mc_jacobian {
	Input.rte.mc.jacobian = 1;
}
mc_lidar_file{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_LIDAR], yytext, 2, 2);
}
mc_maxscatters{SPACE_INTEGER} {
	Input.rte.mc.maxscatters = yytext2int(yytext, 2, 2);
}
mc_panorama{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("circumsolar_var_red", s)==0 )	{ id = PAN_FLAG_CIRCUMSOLAR_VAR_RED; }
	else if ( strcasecmp("distr_photons_over_pixel", s)==0 )	{ id = PAN_FLAG_DISTR_PHOTONS_OVER_PIXEL; }
	else if ( strcasecmp("no_pixel", s)==0 )	{ id = PAN_FLAG_NO_PIXEL; }
	else if ( strcasecmp("quicklook", s)==0 )	{ id = PAN_FLAG_QUICKLOOK; }
	else if ( strcasecmp("weight_with_cos", s)==0 )	{ id = PAN_FLAG_WEIGHT_WITH_COS; }
	else if ( strcasecmp("with_direct_rad", s)==0 )	{ id = PAN_FLAG_WITH_DIRECT_RAD; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rte.mc.pan[id] = 1;
}
mc_panorama_view{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT} {
	Input.rte.mc.pan_angles = calloc (4, sizeof(double));
	Input.rte.mc.pan_angles[0] = yytext2float(yytext, 2, 5);
	Input.rte.mc.pan_angles[1] = yytext2float(yytext, 3, 5);
	Input.rte.mc.pan_angles[2] = yytext2float(yytext, 4, 5);
	Input.rte.mc.pan_angles[3] = yytext2float(yytext, 5, 5);
	Input.rte.mc.panorama = PAN_MODE_CAMERA;
	Input.rte.mc.allocate_umu_and_phi = TRUE;
}
mc_photons{SPACE_INTEGER} {
	Input.rte.mc.photons = yytext2int(yytext, 2, 2);
}
mc_photons_file{SPACE_WORD} {
	yytextCstring ( Input.rte.mc.filename[FN_MC_PHOTONS], yytext, 2, 2);
}
mc_polarisation{SPACE_INTEGER}? {
	ntokens = yytext2ntokens(yytext);
	Input.rte.mc.polarisation = 1;
	Input.rte.mc.reflectalways = 1;
	if ( ntokens >= 1 )	 {
		Input.rte.mc.polarisation_state = yytext2int(yytext, 2, ntokens+1);
	}
}
mc_rad_alpha{SPACE_FLOAT} {
	Input.rte.mc.alpha = yytext2float(yytext, 2, 2);
}
mc_radial_pathlength_dt{SPACE_FLOAT} {
	Input.rte.mc.dt = yytext2float(yytext, 2, 2);
}
mc_randomseed{SPACE_INTEGER} {
	Input.rte.mc.randomseed = yytext2int(yytext, 2, 2);
}
mc_refraction {
	Input.rte.mc.refraction = 1;
}
mc_ris{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("factor", s)==0 )	{ id = MC_RIS_FACTOR; }
	else if ( strcasecmp("optical_depth", s)==0 )	{ id = MC_RIS_OPTICAL_DEPTH; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rte.mc.ris[id] = yytext2float(yytext, 3, 3);
}
mc_satellite_view{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT} {
	Input.rte.mc.pan_angles = calloc (4, sizeof(double));
	Input.rte.mc.pan_angles[0] = yytext2float(yytext, 2, 5);
	Input.rte.mc.pan_angles[1] = yytext2float(yytext, 3, 5);
	Input.rte.mc.pan_angles[2] = yytext2float(yytext, 4, 5);
	Input.rte.mc.pan_angles[3] = yytext2float(yytext, 5, 5);
	Input.rte.mc.panorama = PAN_MODE_SATELLITE;
	Input.rte.mc.allocate_umu_and_phi = TRUE;
}
mc_sensorposition{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT} {
	Input.rte.mc.sensorposition_x = yytext2float(yytext, 2, 4);
	Input.rte.mc.sensorposition_y = yytext2float(yytext, 3, 4);
	Input.rte.mc.sensorposition_z = yytext2float(yytext, 4, 4);
	Input.rte.mc.sensorposition = MC_SENSORPOSITION_CARTESIAN;
}
mc_spherical{SPACE}+{DIM} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("1D", s)==0 )	{ id = DIM_1D; }
	else if ( strcasecmp("3D", s)==0 )	{ id = DIM_3D; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rte.mc.spherical[id] = 1;
}
mc_std {
	Input.rte.mc.std = 1;
}
mc_sun_angular_size{SPACE_SIGNED_FLOAT} {
	Input.rte.mc.sun_radius = yytext2float(yytext, 2, 2);
}
mc_sun_position{SPACE_SIGNED_FLOAT}{SPACE_SIGNED_FLOAT} {
	Input.atm.sza_source = SZA_DIRECT_INPUT;
	Input.atm.phi0_spher = yytext2float(yytext, 2, 3);
	Input.atm.sza_spher = yytext2float(yytext, 3, 3);
}
mc_surface_reflectalways {
	Input.rte.mc.reflectalways = 1;
}
mc_tipa{SPACE_WORD} {
	Input.rte.mc.ipa = 1;
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("dir3d", s)==0 )	{ Input.rte.mc.tipa = TIPA_DIR3D; }
	else if ( strcasecmp("dirdiff", s)==0 )	{ Input.rte.mc.tipa = TIPA_DIRDIFF; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
mc_truncate{SPACE_SIGNED_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.rte.mc.truncate = 0.99;
	if ( ntokens >= 1 )	 {
		Input.rte.mc.truncate = yytext2float(yytext, 2, ntokens+1);
	}
}
mc_vroom{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	Input.rte.mc.vroom = 1;
	if ( ntokens >= 1 )	 {
		s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("off", s)==0 )	{ Input.rte.mc.vroom = SWITCH_OFF; }
	else if ( strcasecmp("on", s)==0 )	{ Input.rte.mc.vroom = SWITCH_ON; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
mixing_ratio{SPACE_WORD}{SPACE_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("CH4", s)==0 )	{ id = MX_CH4; }
	else if ( strcasecmp("CO2", s)==0 )	{ id = MX_CO2; }
	else if ( strcasecmp("F11", s)==0 )	{ id = MX_F11; }
	else if ( strcasecmp("F12", s)==0 )	{ id = MX_F12; }
	else if ( strcasecmp("F22", s)==0 )	{ id = MX_F22; }
	else if ( strcasecmp("H2O", s)==0 )	{ id = MX_H2O; }
	else if ( strcasecmp("N2O", s)==0 )	{ id = MX_N2O; }
	else if ( strcasecmp("NO2", s)==0 )	{ id = MX_NO2; }
	else if ( strcasecmp("O2", s)==0 )	{ id = MX_O2; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.mixing_ratio[id] = yytext2float(yytext, 3, 3);
}
mol_abs_param{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("avhrr_kratz", s)==0 )	{ Input.ck_scheme = CK_AVHRR_KRATZ; }
	else if ( strcasecmp("crs", s)==0 )	{ Input.ck_scheme = CK_CRS; }
	else if ( strcasecmp("fu", s)==0 )	{ Input.ck_scheme = CK_FU; }
	else if ( strcasecmp("kato", s)==0 )	{ Input.ck_scheme = CK_KATO; }
	else if ( strcasecmp("kato2", s)==0 )	{ Input.ck_scheme = CK_KATO2; }
	else if ( strcasecmp("kato2.96", s)==0 )	{ Input.ck_scheme = CK_KATO2_96; }
	else if ( strcasecmp("lowtran", s)==0 )	{ Input.ck_scheme = CK_LOWTRAN; }
	else if ( strcasecmp("reptran", s)==0 )	{ Input.ck_scheme = CK_REPTRAN; }
	else if ( strcasecmp("reptran_channel", s)==0 )	{ Input.ck_scheme = CK_REPTRAN_CHANNEL; }
	else if ( strcasecmp("sbdart", s)==0 )	{ Input.ck_scheme = CK_SBDART; }
	else	{ Input.ck_scheme = CK_FILE;
		 strcpy ( Input.ck_scheme_filename, s ); }
	free(s); 
	if ( ntokens >= 2 )	 {
		Input.ck_reptran_arg = yytext2string(yytext, 3, ntokens+1);
	}
}
mol_file{SPACE_WORD}{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("BRO", s)==0 )	{ mol_id = MOL_BRO; }
	else if ( strcasecmp("CH4", s)==0 )	{ mol_id = MOL_CH4; }
	else if ( strcasecmp("CO", s)==0 )	{ mol_id = MOL_CO; }
	else if ( strcasecmp("CO2", s)==0 )	{ mol_id = MOL_CO2; }
	else if ( strcasecmp("H2O", s)==0 )	{ mol_id = MOL_H2O; }
	else if ( strcasecmp("HCHO", s)==0 )	{ mol_id = MOL_HCHO; }
	else if ( strcasecmp("N2", s)==0 )	{ mol_id = MOL_N2; }
	else if ( strcasecmp("N2O", s)==0 )	{ mol_id = MOL_N2O; }
	else if ( strcasecmp("NO2", s)==0 )	{ mol_id = MOL_NO2; }
	else if ( strcasecmp("O2", s)==0 )	{ mol_id = MOL_O2; }
	else if ( strcasecmp("O3", s)==0 )	{ mol_id = MOL_O3; }
	else if ( strcasecmp("O4", s)==0 )	{ mol_id = MOL_O4; }
	else if ( strcasecmp("OCLO", s)==0 )	{ mol_id = MOL_OCLO; }
	else if ( strcasecmp("SO2", s)==0 )	{ mol_id = MOL_SO2; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.atm.filename[mol_id], yytext, 3, ntokens+1);
	if ( ntokens >= 3 )	 {
		s = yytext2string(yytext, 4, ntokens+1);
	if ( strcasecmp("MMR", s)==0 )	{ Input.atm.unit_profile[mol_id] = MMR; }
	else if ( strcasecmp("RH", s)==0 )	{ Input.atm.unit_profile[mol_id] = RH; }
	else if ( strcasecmp("VMR", s)==0 )	{ Input.atm.unit_profile[mol_id] = VMR; }
	else if ( strcasecmp("cm_3", s)==0 )	{ Input.atm.unit_profile[mol_id] = CM_3; }
	else if ( strcasecmp("m_3", s)==0 )	{ Input.atm.unit_profile[mol_id] = M_3; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 3 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
mol_modify{SPACE_WORD}{SPACE_SIGNED_FLOAT}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 4);
	if ( strcasecmp("BRO", s)==0 )	{ id = MOL_BRO; }
	else if ( strcasecmp("CH4", s)==0 )	{ id = MOL_CH4; }
	else if ( strcasecmp("CO", s)==0 )	{ id = MOL_CO; }
	else if ( strcasecmp("CO2", s)==0 )	{ id = MOL_CO2; }
	else if ( strcasecmp("H2O", s)==0 )	{ id = MOL_H2O; }
	else if ( strcasecmp("HCHO", s)==0 )	{ id = MOL_HCHO; }
	else if ( strcasecmp("N2", s)==0 )	{ id = MOL_N2; }
	else if ( strcasecmp("N2O", s)==0 )	{ id = MOL_N2O; }
	else if ( strcasecmp("NO2", s)==0 )	{ id = MOL_NO2; }
	else if ( strcasecmp("O2", s)==0 )	{ id = MOL_O2; }
	else if ( strcasecmp("O3", s)==0 )	{ id = MOL_O3; }
	else if ( strcasecmp("O4", s)==0 )	{ id = MOL_O4; }
	else if ( strcasecmp("OCLO", s)==0 )	{ id = MOL_OCLO; }
	else if ( strcasecmp("SO2", s)==0 )	{ id = MOL_SO2; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.atm.column[id] = yytext2float(yytext, 3, 4);
	s = yytext2string(yytext, 4, 4);
	if ( strcasecmp("CM_2", s)==0 )	{ Input.atm.unit_column[id] = MOL_UNIT_CM_2; }
	else if ( strcasecmp("DU", s)==0 )	{ Input.atm.unit_column[id] = MOL_UNIT_DU; }
	else if ( strcasecmp("MM", s)==0 )	{ Input.atm.unit_column[id] = MOL_UNIT_MM; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 3 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
mol_tau_file{SPACE_WORD}{SPACE_WORD} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("abs", s)==0 )	{ id = FN_MOL_TAU_ABS; }
	else if ( strcasecmp("sca", s)==0 )	{ id = FN_MOL_TAU_SCA; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.filename[id], yytext, 3, 3);
}
no_absorption{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	ispoff = get_caothoff_index(&Input.caothoff,&Input.n_caothoff, "all");
	if ( ntokens >= 1 )	 {
		ispoff = get_caothoff_index(&Input.caothoff,&Input.n_caothoff, yytext2string(yytext, 2, ntokens+1));
	}
	Input.caothoff[ispoff].no_absorption = 1;
}
no_scattering{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	ispoff = get_caothoff_index(&Input.caothoff,&Input.n_caothoff, "all");
	if ( ntokens >= 1 )	 {
		ispoff = get_caothoff_index(&Input.caothoff,&Input.n_caothoff, yytext2string(yytext, 2, ntokens+1));
	}
	Input.caothoff[ispoff].no_scattering = 1;
}
number_of_streams{SPACE_INTEGER} {
	Input.rte.nstr = yytext2int(yytext, 2, 2);
}
output_file{SPACE_WORD} {
	Input.filename[FN_OUTFILE] = yytext2string(yytext, 2, 2);
}
output_format{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("ascii", s)==0 )	{ Input.output_format = OUTPUT_FORMAT_ASCII; }
	else if ( strcasecmp("flexstor", s)==0 )	{ Input.output_format = OUTPUT_FORMAT_FLEXSTOR; }
	else if ( strcasecmp("netCDF", s)==0 )	{ Input.output_format = OUTPUT_FORMAT_NETCDF; }
	else if ( strcasecmp("sat_picture", s)==0 )	{ Input.output_format = OUTPUT_FORMAT_SAT_PICTURE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
output_process{SPACE_WORD} {
	Input.output_unit_processing = yytext2string(yytext, 2, 2);
}
output_quantity{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("brightness", s)==0 )	{ Input.calibration = OUTCAL_BRIGHTNESS; }
	else if ( strcasecmp("reflectivity", s)==0 )	{ Input.calibration = OUTCAL_REFLECTIVITY; }
	else if ( strcasecmp("transmittance", s)==0 )	{ Input.calibration = OUTCAL_TRANSMITTANCE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
phi{SIGNED_FLOATS} {
	ntokens = yytext2ntokens(yytext);
	Input.rte.phi =  yytext2floats(yytext, ntokens);
	Input.rte.nphi = ntokens;
	Input.rte.maxphi = ntokens;
}
phi0{SPACE_SIGNED_FLOAT} {
	Input.atm.phi0 = yytext2float(yytext, 2, 2);
}
polradtran{SPACE_WORD}{SPACE_INTEGER} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("aziorder", s)==0 )	{ id = POLRADTRAN_AZIORDER; }
	else if ( strcasecmp("nstokes", s)==0 )	{ id = POLRADTRAN_NSTOKES; }
	else if ( strcasecmp("src_code", s)==0 )	{ id = POLRADTRAN_SRC_CODE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rte.polradtran[id] = yytext2int(yytext, 3, 3);
}
polradtran_max_delta_tau{SPACE_SIGNED_FLOAT} {
	Input.rte.pol_max_delta_tau = yytext2float(yytext, 2, 2);
}
polradtran_quad_type{SPACE_WORD} {
	Input.rte.pol_quad_type = yytext2string(yytext, 2, 2);
}
pressure{SPACE_FLOAT} {
	Input.pressure = yytext2float(yytext, 2, 2);
}
print_disort_info{INTEGERS} {
	ntokens = yytext2ntokens(yytext);
	Input.rte.prndis =  yytext2integers(yytext, ntokens);
	Input.rte.nprndis = ntokens;
}
profile_file{SPACE_WORD}{SPACE_option_definition.ProfileType}{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, 2, 4),0);
	s = yytext2string(yytext, 3, 4);
	if ( strcasecmp("1d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_1D; }
	else if ( strcasecmp("3d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_3D; }
	else if ( strcasecmp("ipa_files", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_IPA_FILES; }
	else if ( strcasecmp("moments", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_MOMENTS; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.caoth[isp].filename, yytext, 4, 4);
}
profile_modify{SPACE_WORD}{SPACE_WORD}{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, 2, 5),0);
	s = yytext2string(yytext, 3, 5);
	if ( strcasecmp("gg", s)==0 )	{ id1 = MODIFY_VAR_GG; }
	else if ( strcasecmp("ssa", s)==0 )	{ id1 = MODIFY_VAR_SSA; }
	else if ( strcasecmp("tau", s)==0 )	{ id1 = MODIFY_VAR_TAU; }
	else if ( strcasecmp("tau550", s)==0 )	{ id1 = MODIFY_VAR_TAU550; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 4, 5);
	if ( strcasecmp("scale", s)==0 )	{ id2 = MODIFY_TYPE_SCALE; }
	else if ( strcasecmp("set", s)==0 )	{ id2 = MODIFY_TYPE_SET; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 3 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.caoth[isp].modify[id1][id2] = yytext2float(yytext, 5, 5);
}
profile_properties{SPACE_WORD}{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth, yytext2string(yytext, 2, ntokens+1),0);
	s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("baum", s)==0 )	{ Input.caoth[isp].properties = PROP_BAUM; }
	else if ( strcasecmp("baum_v36", s)==0 )	{ Input.caoth[isp].properties = PROP_BAUM_V36; }
	else if ( strcasecmp("echam4", s)==0 )	{ Input.caoth[isp].properties = PROP_ECHAM4; }
	else if ( strcasecmp("fu", s)==0 )	{ Input.caoth[isp].properties = PROP_FU; }
	else if ( strcasecmp("hey", s)==0 )	{ Input.caoth[isp].properties = PROP_HEY; }
	else if ( strcasecmp("hu", s)==0 )	{ Input.caoth[isp].properties = PROP_HU; }
	else if ( strcasecmp("key", s)==0 )	{ Input.caoth[isp].properties = PROP_KEY; }
	else if ( strcasecmp("mie", s)==0 )	{ Input.caoth[isp].properties = PROP_MIE; }
	else if ( strcasecmp("yang", s)==0 )	{ Input.caoth[isp].properties = PROP_YANG; }
	else	{ Input.caoth[isp].properties = PROP_FILE;
		 strcpy ( Input.caoth[isp].properties_filename, s ); }
	free(s); 
	if ( ntokens >= 3 )	 {
		s = yytext2string(yytext, 4, ntokens+1);
	if ( strcasecmp("interpolate", s)==0 )	{ Input.caoth[isp].interpolate = 1; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 3 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
pseudospherical {
	Input.rte.pseudospherical = 1;
}
quiet {
	Input.quiet = 1;
}
radiosonde_levels_only {
	Input.atm.rs_add_upper_levels = FALSE;
}
raman{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	Input.raman = RAMAN_CALC;
	Input.ck_scheme = CK_RAMAN;
	Input.rte.solver = SOLVER_DISORT;
	Input.processing = PROCESS_RAMAN;
	if ( ntokens >= 1 )	 {
		s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("original", s)==0 )	{ Input.raman_original = 1; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
rayleigh_depol{SPACE_SIGNED_FLOAT} {
	Input.rayleigh_depol = yytext2float(yytext, 2, 2);
}
refractive_index_file{SPACE_WORD} {
	Input.filename[FN_REFIND] = yytext2string(yytext, 2, 2);
}
reptran_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_REPTRAN], yytext, 2, 2);
}
reverse_atmosphere {
	Input.rte.reverse = 1;
}
rte_solver{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("disort", s)==0 )	{ Input.rte.solver = SOLVER_DISORT; }
	else if ( strcasecmp("fdisort1", s)==0 )	{ Input.rte.solver = SOLVER_FDISORT1; }
	else if ( strcasecmp("fdisort2", s)==0 )	{ Input.rte.solver = SOLVER_FDISORT2; }
	else if ( strcasecmp("ftwostr", s)==0 )	{ Input.rte.solver = SOLVER_FTWOSTR; }
	else if ( strcasecmp("montecarlo", s)==0 )	{ Input.rte.solver = SOLVER_MONTECARLO; }
	else if ( strcasecmp("mystic", s)==0 )	{ Input.rte.solver = SOLVER_MYSTIC; }
	else if ( strcasecmp("null", s)==0 )	{ Input.rte.solver = SOLVER_NULL; }
	else if ( strcasecmp("polradtran", s)==0 )	{ Input.rte.solver = SOLVER_POLRADTRAN; }
	else if ( strcasecmp("rodents", s)==0 )	{ Input.rte.solver = SOLVER_RODENTS; }
	else if ( strcasecmp("sdisort", s)==0 )	{ Input.rte.solver = SOLVER_SDISORT; }
	else if ( strcasecmp("sos", s)==0 )	{ Input.rte.solver = SOLVER_SOS; }
	else if ( strcasecmp("spsdisort", s)==0 )	{ Input.rte.solver = SOLVER_SPSDISORT; }
	else if ( strcasecmp("sslidar", s)==0 )	{ Input.rte.solver = SOLVER_SSLIDAR; }
	else if ( strcasecmp("sss", s)==0 )	{ Input.rte.solver = SOLVER_SSS; }
	else if ( strcasecmp("sssi", s)==0 )	{ Input.rte.solver = SOLVER_SSSI; }
	else if ( strcasecmp("twostr", s)==0 )	{ Input.rte.solver = SOLVER_TWOSTR; }
	else if ( strcasecmp("twostrebe", s)==0 )	{ Input.rte.solver = SOLVER_TWOSTREBE; }
	else if ( strcasecmp("tzs", s)==0 )	{ Input.rte.solver = SOLVER_TZS; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
satellite_geometry{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_SATGEOMETRY], yytext, 2, 2);
}
satellite_pixel{SPACE_INTEGER}{SPACE_INTEGER} {
	Input.sat_pixel_x = yytext2int(yytext, 2, 3);
	Input.sat_pixel_y = yytext2int(yytext, 3, 3);
}
sdisort{SPACE_WORD}{SPACE_INTEGER} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("ichapman", s)==0 )	{ id = SDISORT_ICHAPMAN; }
	else if ( strcasecmp("nrefrac", s)==0 )	{ id = SDISORT_NREFRAC; }
	else if ( strcasecmp("nscat", s)==0 )	{ id = SDISORT_NSCAT; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.rte.sdisort[id] = yytext2int(yytext, 3, 3);
}
slit_function_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_SLITFUNCTION], yytext, 2, 2);
	Input.convolve = 1;
}
sos_nscat{SPACE_INTEGER} {
	Input.rte.sos_nscat = yytext2int(yytext, 2, 2);
}
source{SPACE_WORD}{SPACE_WORD}?{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("solar", s)==0 )	{ Input.source = SRC_SOLAR; }
	else if ( strcasecmp("thermal", s)==0 )	{ Input.source = SRC_THERMAL; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.filename[FN_EXTRATERRESTRIAL], yytext, 3, ntokens+1);
	}
	if ( ntokens >= 3 )	 {
		s = yytext2string(yytext, 4, ntokens+1);
	if ( strcasecmp("per_band", s)==0 )	{ Input.spectrum_unit = UNIT_PER_BAND; }
	else if ( strcasecmp("per_cm-1", s)==0 )	{ Input.spectrum_unit = UNIT_PER_CM_1; }
	else if ( strcasecmp("per_nm", s)==0 )	{ Input.spectrum_unit = UNIT_PER_NM; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 3 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
spline{SPACE_FLOAT}{SPACE_FLOAT}{SPACE_FLOAT} {
	Input.spline_lambda_0 = yytext2float(yytext, 2, 4);
	Input.spline_lambda_1 = yytext2float(yytext, 3, 4);
	Input.spline_lambda_step = yytext2float(yytext, 4, 4);
	Input.spline = 1;
}
spline_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_SPLINE], yytext, 2, 2);
	Input.spline = 1;
}
sslidar{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("E0", s)==0 )	{ id = SSLIDAR_E0; }
	else if ( strcasecmp("area", s)==0 )	{ id = SSLIDAR_AREA; }
	else if ( strcasecmp("eff", s)==0 )	{ id = SSLIDAR_EFF; }
	else if ( strcasecmp("position", s)==0 )	{ id = SSLIDAR_POSITION; }
	else if ( strcasecmp("range", s)==0 )	{ id = SSLIDAR_RANGE; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.sslidar[id] = yytext2float(yytext, 3, 3);
}
sslidar_nranges{SPACE_INTEGER} {
	Input.sslidar_nranges = yytext2int(yytext, 2, 2);
}
sslidar_polarisation {
	Input.sslidar_polarisation = 1;
}
sur_temperature{SPACE_SIGNED_FLOAT} {
	Input.surface_temperature = yytext2float(yytext, 2, 2);
}
surface_temperature_map{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	yytextCstring ( Input.filename[FN_SURFACE_TEMP_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.netCDF_name_surf_T, yytext, 3, ntokens+1);
	}
}
surface_type_map{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	Input.alb.surf_type_map = TRUE;
	yytextCstring ( Input.filename[FN_SURFACE_TYPE_MAP], yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		yytextCstring ( Input.alb.netCDF_surf_name, yytext, 3, ntokens+1);
	}
}
sza{SPACE_FLOAT} {
	Input.atm.sza = yytext2float(yytext, 2, 2);
	Input.atm.sza_source = SZA_DIRECT_INPUT;
}
sza_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_SZA], yytext, 2, 2);
	Input.atm.sza_source = SZA_FROM_SZA_FILE;
}
test_optical_properties {
	Input.test_optical_properties = 1;
}
thermal_bands_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_WLBANDS], yytext, 2, 2);
}
thermal_bandwidth{SPACE_FLOAT}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	Input.bandwidth = yytext2float(yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("cm-1", s)==0 )	{ Input.bandwidth_unit = UNIT_PER_CM_1; }
	else if ( strcasecmp("nm", s)==0 )	{ Input.bandwidth_unit = UNIT_PER_NM; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
tipa{SPACE_WORD} {
	s = yytext2string(yytext, 2, 2);
	if ( strcasecmp("dir", s)==0 )	{ Input.tipa = TIPA_DIR; }
	else if ( strcasecmp("dirdiff", s)==0 )	{ Input.tipa = TIPA_DIRDIFF; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
}
umu{SIGNED_FLOATS} {
	ntokens = yytext2ntokens(yytext);
	Input.rte.umu =  yytext2floats(yytext, ntokens);
	Input.rte.numu = ntokens;
	Input.rte.maxumu = ntokens;
}
verbose {
	Input.verbose = 1;
}
wavelength{SPACE_FLOAT}{SPACE_FLOAT}? {
	ntokens = yytext2ntokens(yytext);
	Input.wl.start = yytext2float(yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		Input.wl.end = yytext2float(yytext, 3, ntokens+1);
	}
}
wavelength_grid_file{SPACE_WORD} {
	yytextCstring ( Input.filename[FN_WLTRANS], yytext, 2, 2);
}
wavelength_index{SPACE_INTEGER}{SPACE_INTEGER}? {
	ntokens = yytext2ntokens(yytext);
	Input.wl.start_index = yytext2int(yytext, 2, ntokens+1);
	if ( ntokens >= 2 )	 {
		Input.wl.end_index = yytext2int(yytext, 3, ntokens+1);
	}
}
wc_file{SPACE_option_definition.ProfileType}{SPACE_WORD} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"wc",1);
	s = yytext2string(yytext, 2, 3);
	if ( strcasecmp("1d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_1D; }
	else if ( strcasecmp("3d", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_3D; }
	else if ( strcasecmp("ipa_files", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_IPA_FILES; }
	else if ( strcasecmp("moments", s)==0 )	{ Input.caoth[isp].source = CAOTH_FROM_MOMENTS; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	yytextCstring ( Input.caoth[isp].filename, yytext, 3, 3);
}
wc_modify{SPACE_WORD}{SPACE_WORD}{SPACE_SIGNED_FLOAT} {
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"wc",1);
	s = yytext2string(yytext, 2, 4);
	if ( strcasecmp("gg", s)==0 )	{ id1 = MODIFY_VAR_GG; }
	else if ( strcasecmp("ssa", s)==0 )	{ id1 = MODIFY_VAR_SSA; }
	else if ( strcasecmp("tau", s)==0 )	{ id1 = MODIFY_VAR_TAU; }
	else if ( strcasecmp("tau550", s)==0 )	{ id1 = MODIFY_VAR_TAU550; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 1 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	s = yytext2string(yytext, 3, 4);
	if ( strcasecmp("scale", s)==0 )	{ id2 = MODIFY_TYPE_SCALE; }
	else if ( strcasecmp("set", s)==0 )	{ id2 = MODIFY_TYPE_SET; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	Input.caoth[isp].modify[id1][id2] = yytext2float(yytext, 4, 4);
}
wc_properties{SPACE_WORD}{SPACE_WORD}? {
	ntokens = yytext2ntokens(yytext);
	isp = get_caoth_index(&Input.caoth,&Input.n_caoth,"wc",1);
	s = yytext2string(yytext, 2, ntokens+1);
	if ( strcasecmp("echam4", s)==0 )	{ Input.caoth[isp].properties = PROP_ECHAM4; }
	else if ( strcasecmp("hu", s)==0 )	{ Input.caoth[isp].properties = PROP_HU; }
	else if ( strcasecmp("mie", s)==0 )	{ Input.caoth[isp].properties = PROP_MIE; }
	else	{ Input.caoth[isp].properties = PROP_FILE;
		 strcpy ( Input.caoth[isp].properties_filename, s ); }
	free(s); 
	if ( ntokens >= 2 )	 {
		s = yytext2string(yytext, 3, ntokens+1);
	if ( strcasecmp("interpolate", s)==0 )	{ Input.caoth[isp].interpolate = 1; }
	else	{fprintf(stderr,"Option %s  has an invalid argument 2 on line %d: %s\n", yytext, line_number+1, s); ierror++;} 
	free(s); 
	}
}
zout_interpolate {
	Input.atm.zout_interpolate = ZOUT_INTERPOLATE;
}

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

\option{aerosol_angstrom}
	Scale the aerosol optical depth using the {\AA}ngstr{\"o}m formula:
	\begin{equation}
	\tau = \beta \lambda^{-\alpha}
	\end{equation}
	where $\lambda$ is in units of micrometer \citep{Angstrom1929}. Specify 
	the {\AA}ngstr{\"o}m alpha and beta coefficients by
	\fcode{
	   aerosol\_angstrom alpha beta
	}
	The optical thickness defined here  is the integral from the user-defined 
	\code{altitude} to TOA (top of atmosphere). 
		
\option{aerosol_default}
	Set up a default aerosol according to \citet{shettle89}. 
	The default properties are 
	a rural type aerosol in the boundary layer, background aerosol above 2km,
	spring-summer conditions and a visibility of 50km. These settings may
	be modified with \code{aerosol\_haze}, \code{aerosol\_vulcan}, \code{aerosol\_season},
	and \code{aerosol\_visibility}.
		
\option{aerosol_file}
	Location of file defining aerosol optical properties.
	\fcode{
	   aerosol\_file type file
	}
	\code{type} defines the file type, which can be one of the following:
	\begin{description}
	\parameter{gg} Location of aerosol asymmetry parameter file.

	The file must have two columns.
	Column 1 is the altitude in km. 
	Column 2 is the asymmetry parameter of each layer. 
	The asymmetry parameter defined with this option is constant with wavelength.
	If you require spectral dependence please use \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{ssa} Location of aerosol single scattering albedo file. 

	The file must have two columns.
	Column 1 is the altitude in km. The altitude grid must be exactly equal to 
	the altitude grid specified in the file \code{atmosphere\_file}.
	Column 2 is the single scattering albedo of each layer. 
	The single scattering albedo defined with this option is constant with wavelength. If you require spectral dependence please use \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{tau} Location of aerosol optical depth file. 

	The file must have two columns.
	Column 1 is the altitude in km. The altitude grid must be exactly equal to 
	the altitude grid specified in the file \code{atmosphere\_file}.
	Column 2 is the aerosol optical depth of each layer. 
	To allow wavelength-dependent aerosol optical thickness please
	use either \code{aerosol\_angstrom} or \code{aerosol\_file explicit}.
	Comments start with \code{\#}. Empty lines are ignored.
	\parameter{moments} Set the aerosol phase function moments to the values specified in 
	the aerosol moments file.

	The file contains one column with arbitrary number of 
	Legendre terms of the phase function. The phase function 
	 $p(\mu)$ is
	\begin{equation}
	   p (\mu) = \sum_{m=0}^{\infty} (2m+1) \cdot k_m \cdot P_m (\mu) 
	\end{equation}
	where $k_m$ is the m'th moment and $P_m (\mu)$ is the m'th Legendre
	polynomial.  If not specified, a Henyey-Greenstein phase function is
	assumed where the asymmetry parameter g is either a default value
	depending on the aerosol type or it may be specified using
	\code{aerosol\_modify set gg}.  The phase function will be the same for all
	altitudes and wavelengths.  See \code{aerosol\_file explicit} if more
	flexibility is wanted. May only be used together with the
	\code{disort} or \code{fdisort2} solver in combination with the option
	\code{disort_intcor moments}.
	\parameter{explicit} A way to specify aerosol optical depth, single scattering albedo, 
	and phase function moments for each layer. 

	The file must have two columns where column 1 is the altitude in km. The second
	column is a the name of a file which defines the optical properties of the layer
	starting at the given altitude. The files specified in the second column must have the 
	following format:
	\begin{description}
	\item[Column 1:]
	The wavelength in nm. These wavelengths may be different from those in \code{source solar file}. 
	Optical properties are interpolated to the requested wavelengths.
	\item[Column 2:]
	The extinction coefficient of the layer in units km-1. 
	\item[Column 3:]
	The aerosol single scattering albedo of the layer.
	\item[Column 4-(nmom+4):]
	The moments of the aerosol phase function. 
	\end{description}
	For some simple examples see the files \file{examples/AERO\_*.LAYER}. Note that 
	if using the \code{rte_solver disort} it makes good sense to make the 
	number of moments larger than \code{number_of_streams}. For \code{rte_solver fdisort1} and 
	\code{rte_solver polradtran} the number of moments included in the calculations 
	will be \code{number_of_streams}+1. Higher order moments will be ignored for these solvers.
	Please note that the uppermost line of the \code{aerosol_file explicit}
	denotes simply the top altitude of the uppermost layer. The optical
	properties of this line are consequently ignored. There are two
	options for this line: either an optical property file with zero
	optical thickness is specified or "NULL" us used.
	\end{description}
		
\option{aerosol_haze}
	Specify the aerosol type in the lower 2 km of the atmosphere as
	\fcode{
	   aerosol\_haze type
	}
	where \code{type} is an integer identifying the following aerosol types:
	\begin{description}
	\item[1] Rural type aerosols.
	\item[4] Maritime type aerosols.
	\item[5] Urban type aerosols.
	\item[6] Tropospheric type aerosols.
	\end{description}
	For a description of the different aerosol types see \citet{shettle89}.
		
\option{aerosol_king_byrne}
	Scale the aerosol optical depth using the King Byrne formula \citep{king1976}:
	\begin{equation}
	\tau = e^{\alpha_0} \cdot \lambda^{\alpha_1} \lambda^{-\alpha}
	\end{equation}
	where $\lambda$ is in units of micrometer \citep{Angstrom1929}. Specify 
	the King Byrne $alpha_0$ $alpha_1$ $alpha_2$ coefficients by
	\fcode{
	   aerosol\_king\_byrne alpha_0 alpha_1 alpha_2
	}
	The optical thickness defined here  is the integral from the user-defined 
	\code{altitude} to TOA (top of atmosphere). 
		
\option{aerosol_modify}
	Modify aerosol optical properties.
	\fcode{
	aerosol\_modify variable scale/set value
	}
	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
		
\option{aerosol_profile_modtran}
	Squeeze aerosol profile up to 6 km when altitude is non-zero as in MODTRAN.
	Per default the aerosol profile is shifted upwards
	and remains unchanged. 
		
\option{aerosol_refrac_index}
        Calculate optical properties from size distribution and index of refraction using Mie
        theory. Here is an exception from the rule that ALL values defined above are overwritten
        because the optical thickness profile is re-scaled so that the optical thickness
        at the first internal wavelength is unchanged. It is done that way to give the user an
        easy means of specifying the optical thickness at a given wavelength.
        
\option{aerosol_season}
	Specify season to get appropriate aerosol profile.
	\fcode{
	   aerosol\_season season
	}
	where \code{season} is either \code{1} or \code{2}:
	\begin{description}
	\item[1] Spring-summer profile.
	\item[2] Fall-winter profile.
	\end{description}
		
\option{aerosol_set_tau_at_wvl}
	Set the aerosol optical thickness at wavelength lambda (nm). Other wavelengths are scaled accordingly.
	Note that this option requires for technical reasons that the
	wavelength interval defined by \code{wavelength} does contain \code{lambda}.
	The optical thickness defined here is the integral from the user-definded 
	\code{altitude} to TOA (top of atmosphere). 
	\fcode{
	   aerosol\_set\_tau\_at\_wvl lambda tau
	}
		
\option{aerosol_sizedist_file}
        Calculate optical properties from size distribution and index of refraction using Mie
        theory. Here is an exception from the rule that ALL values defined above are overwritten
        because the optical thickness profile is re-scaled so that the optical thickness
        at the first internal wavelength is unchanged. It is done that way to give the user an
        easy means of specifying the optical thickness at a given wavelength.
        
\option{aerosol_species_file}
	Specify mass density profiles of a mixture of aerosol types. 
	\fcode{
	aerosol\_species\_file profile [aero\_1 aero\_2 ... aero\_n]
	}
	where \code{aero\_1} to \code{aero\_n} are the aerosol species to be included. 
	For each of these species, the optical properties are read from the 
	\code{aerosol\_species\_library}, e.g. the OPAC data set provided with 
	libRadtran. The profile file needs to include vertical profiles for each
	of these species. This file can be either in \emph{netCDF}-format 
	(automatically recognized filename extension \code{.nc} or \code{.cdf}) 
	or in ASCII format. The format of the ASCII file is: 
	\fcode{
	z1    dens(aero\_1, z1) dens(aero\_2, z1)  ... dens(aero\_n, z1)\\
	z2    dens(aero\_1, z2) dens(aero\_2, z2)  ...\\
	 .        .           .\\
	 .        .           .\\
	}
	where \code{z} is the height in km, and \code{dens} are the aerosol mass
	densities in g/m3. Please make sure to include one column for each of 
	the species \code{aero\_1} to \code{aero\_n} listed after 
	\code{aerosol\_species\_file}. 
	For netCDF input it is also possible to specify the unit 'kg kg$^{-1}$'; the data 
	are then automatically converted to g/m$^3$.
	
	Some default aerosol mixtures are provided, corresponding to the definitions in 
	\citet{hess98:_optic_proper_aeros_cloud}.
	They can simply be invoked by
	\fcode{
	aerosol\_species\_file mixture\_name
	}
	where \code{mixture\_name} can be one of the following: 
	\fcode{\\
	continental\_clean\\
	continental\_average\\
	continental\_polluted\\
	urban\\
	maritime\_clean\\
	maritime\_polluted\\
	maritime\_tropical\\
	desert\\
	antarctic
	}
	A variation of the desert mixture containing nonspherical particles is 
	\fcode{
	desert\_spheroids
	}
		
\option{aerosol_species_library}
	With this option the \emph{directory} is specified where the optical property
	files for all aerosols species used in the \code{aerosol\_species\_file} are
	expected: For each species defined in \code{aerosol\_species\_file},
	\emph{netCDF}-file \emph{species\_name}\code{.nc},
	(e.g. \code{INSO.nc}), which contains the optical properties of the aerosol
	species, has to be provided. The netcdf format is the one produced by 
	the {\sl libRadtran} \code{mie} tool.
	 
	At the libRadtran webpage we provide the OPAC data set 
	\citep{hess98:_optic_proper_aeros_cloud}
	which can be directly used with \code{uvspec}: 
	\fcode{
	aerosol\_species\_library OPAC
	}
	OPAC contains following aerosol species: 
	\fcode{INSO   insoluble \\
	WASO   water\_soluble\\
	SOOT   soot\\
	SSAM   sea\_salt\_accumulation\_mode\\
	SSCM   sea\_salt\_coarse\_mode\\
	MINM   mineral\_nucleation\_mode\\
	MIAM   mineral\_accumulation\_mode\\
	MICM   mineral\_coarse_mode\\
	MITR   mineral\_transported\\
	SUSO   sulfate\_droplets\\
	}
	Variations of the mineral aerosol species containing nonspherical particles are:
	\fcode{MINM_SPHEROIDS   mineral\_nucleation\_mode\\
	MIAM_SPHEROIDS   mineral\_accumulation\_mode\\
	MICM_SPHEROIDS   mineral\_coarse\_mode
	}
	The nonspherical mineral components are described by \citet{koepke2015}.
		
\option{aerosol_visibility}
	Horizontal visibility in km. Affects the profile according to \citet{shettle89}
	and the optical thickness. 
	\fcode{
	   aerosol\_visibility value
	}
		
\option{aerosol_vulcan}
	Aerosol situation above 2 km as defined in \citet{shettle89}.
	\fcode{
	   aerosol\_vulcan value
	}
	where \code{value} is an integer choosing between the following models
	\begin{description}
	\item[1] Background aerosols.
	\item[2] Moderate volcanic aerosols.
	\item[3] High volcanic aerosols.
	\item[4] Extreme volcanic aerosols.
	\end{description}
		
\option{albedo}
	The Lambertian surface albedo
	\fcode{
	   albedo value
	}
	where \code{value} is a number between 0.0 and 1.0,
	constant for all wavelengths. For wavelength dependent surface albedo 
	use \code{albedo\_file}. The default albedo is 0.0.
		
\option{albedo_file}
	Location of surface albedo file for wavelength dependent surface albedo. 
	\fcode{
	albedo\_file file
	}
	The file must have two columns.
	Column 1 is the wavelength in nm, and column 2 the corresponding 
	Lambertian surface albedo. An arbitrary wavelength grid may be chosen as the 
	albedo will be interpolated linearely to the wavelength grid used for the 
	radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
	A large collection of spectral albedos are available e.g. at http://speclib.jpl.nasa.gov/ 
	\citep{Baldridge2009}.
		
\undocumented{
\option{albedo_library}
	Albedo libraries are a collection of spectral albedos of different surface types. 
	This option must be used either with \code{bdrf_rpv_type} 
	or \code{surface_type_map}, in order to select the specific surface type.
	There are two possibilities for libraries: 
	the built-in IGBP library or a user defined albedo library.

	The built-in library of the International Geosphere Biosphere Programme is 
	selected with 
	\fcode{
	     albedo_library IGBP
	}
	The IGBP library contains 20 surface types which are set by \code{brdf_rpv_type}:
	\fcode{
	 1 evergreen_needle_forest\\
	 2 evergreen_broad_forest\\
	 3 deciduous_needle_forest\\
	 4 deciduous_broad_forest\\
	 5 mixed_forest\\
	 6 closed_shrub\\
	 7 open_shrubs\\
	 8 woody_savanna\\
	 9 savanna\\
	10 grassland\\
	11 wetland\\
	12 cropland\\
	13 urban\\
	14 crop_mosaic\\
	15 antarctic_snow\\
	16 desert\\
	17 ocean_water\\
	18 tundra\\
	19 fresh_snow\\
	20 sea_ice\\
	}
	Surface types 1 - 17 are defined by the International Geosphere Biosphere 
	Programme (IGBP); additionally there are tundra, fresh\_snow, and sea\_ice surface types.
	The spectral albedo of the ground is determined as a function of solar 
	zenith angle, precitable water, and clouds. The spectral resolution 
	equals the grid of the correlated-k Fu/Liou parameterisation.
	This library originates from the NASA CERES/SARB Surface Properties Project, 
	see \citet{Belward1996}.

	For creating your own albedo library use \code{albedo_library path}, where \code{path} is 
	the path of the directory where the albedo data is stored. The files are expected to have the 
	names \code{albedo_01.dat, albedo_02.dat, ...} If \code{brdf_rpv_type 1} is specified
	the albedo from \code{albedo_01.dat} will be used, and so on.  
	Each file is required to have two columns:
	Column 1 is the wavelength in nm, and column 2 the corresponding 
	Lambertian surface albedo. The wavelength grid may be freely set. The 
	albedo will be interpolated linearely to the wavelength grid used for the 
	radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
	This option is similar to \code{albedo_file}, except that it offers an
	easy way to use the option \code{surface_type_map} in combinition with albedo files.
		}
\undocumented{
\option{albedo_map}
	\emph{This option is preliminary and still subject to change (no wavelength dependency yet)!}
	It gives the possibility to specify a wavelength independent albedo with the help of 
	a \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	albedo_map file [variable_name]
	}
	Here \code{file} is the location of the \emph{netCDF} file,
	The optional argument allows the name of the albedo variable in 
	the \emph{netCDF} file to be specified (the default name is \code{AL}).
	The albedo must be provided as function of \code{latitude} and \code{longitude} \code{AL(lat, lon)}, 
	and may also depend on time \code{AL(time, lat, lon)}. 
	The latitude, longitude, and time grids must be provided as doubles 
	\code{double lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
	\code{uvspec} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values are performed.
	If a time-dependent albedo is provided, the albedo data nearest to the specified time 
	will be selected (or linear interpolated if \code{time_interpolate} is switched on). 
		}
\option{altitude}
	Set the bottom level in the model atmosphere provided in 
	\code{atmosphere\_file} to be at the given altitude above sea level (km).
	\fcode{
	altitude 0.73   \# Altitude of IFU, Garmisch-Partenkirchen 
	}
	The profiles of pressure, temperature, molecular absorbers, 
	ice and water clouds are cut at the specified altitude. 
	The aerosol profile is not affected by \code{altitude} but starts
	right from the model surface. This is a convenient way for the user to calculate the 
	radiation at other altitudes than sealevel. Note that \code{altitude} is very different
	from \code{zout} where the radiation is calculated at an altitude of zout
	above the surface. E.g. to calculate the radiation field 1 km above the surface 
	at a location at 0.73 km above sealevel, one would specify '\code{altitude} 0.73' 
	and '\code{zout} 1.0'.
	If an altitude is specified which is below the lowest level in the 
	\code{atmosphere\_file}, the atmospheric profiles are extrapolated assuming a constant 
	gradient for temperature and mixing ratios.
	A second optional argument may be given to \code{altitude} as e.g.
	\fcode{
	altitude 0.73 0.5
	}
	Here the bottom level will be at 0.73 km and the vertical resolution
	of the model atmosphere will be redistributed to have a spacing
	between levels specified by the second number, here 0.5 km, starting
	however from 0km. (Levels 0.73, 1., 1.5 ... will be added to the
	original atmosphere grid and optical properties are devided into the
	new layers. In order to use interpolated properties use
	\code{zout\_interpolate}.  See verbose output for details.)  Be aware
	that specifying a fine vertical spacing will produce many layers thus
	increasing the computing time. Also the radiative transfer equation
	solvers implemented in Fortran 77 might need to have some array sizes
	increased (see \file{src\_f/DISORT.MXD}).
		
\undocumented{
\option{altitude_map}
	Specifies an altitude map which is used in combinition with  
	\code{latitude}, \code{longitude} in order to  
	select the altitude for the simulation. 
	No interpolation is done between the pixels of the map. 
	The format of the call is: 
	\fcode{ 
	altitude_map file [variable_name] 
	} 
	where \code{file} is the location of the altitude map file. The map is expected to  
	be in \emph{netCDF} format. The file must contain \code{double lat(lat)},  
	\code{double lon(lon)}, and the altitude variable, where \code{variable_name} is the  
	name of the surface elevation variable in the \emph{netCDF} file. The default name is Z. 
	The altitude variable must be \code{altitude(lat, lon)}.  
	For format discribtion see also the example map included in {\sl libRadtran},  
	\file{data/altitude/ELEVATION_GTOPO_10min.cdf}. 
	To use this map in \code{uvspec}, you may also use \code{altitude_map GTOPO}. This map has a resolution  
	of 10 arc minutes and the unit of the altitude is meter. Please note that this resolution 
	might not ne adequate for your application.  
	If an altitude in the map is below the lowest level of the  
	\code{atmosphere_file}, the atmospheric profiles are extrapolated assuming a constant  
	gradient for temperature and mixing ratios. 
		}
\option{atm_z_grid}
	With this option the vertical resolution of the \code{atmosphere\_file} data is changed to 
	the levels (in km above sea surface) given as argument. This might be useful in oder to reduce the number 
	of levels (save computational time)\ifmystic{ or in order to easily adjust the atmosphere profile to the resolution 
	of a Monte Carlo cloud file \code{wc\_file 3D} or \code{ic\_file 3D}}.
	\fcode{
	atm\_z\_grid 0 2 4 6 8 10 20 30 ...
	}
		
\undocumented{
\option{atmos_region}
                Very experimental option, ask Robert Buras or Luca Bugliaro for usage.
		}
\option{atmosphere_file}
	Location of the atmospheric data file. 
	\fcode{
	atmosphere\_file file
	}
	The file must have at least three columns containing the altitude,
	pressure, and temperature. Missing profiles are filled with 0 (e.g., if you did not specify 
	the ozone profile, there will be no ozone absorption!), with exception of the air density which 
	is calculated from pressure and temperature. Other trace gases may be set by \code{mol\_file}.  
	The columns are interpreted as follows:
	\begin{description}
	\item[1] Altitude above sea level in km
	\item[2] Pressure in hPa
	\item[3] Temperature in K
	\item[4] air density in cm$^{-3}$
	\item[5] Ozone density in cm$^{-3}$
	\item[6] Oxygen density in cm$^{-3}$
	\item[7] Water vapour density in cm$^{-3}$
	\item[8] CO2 density in cm$^{-3}$
	\item[9] NO2 density in cm$^{-3}$
	\end{description}
	The atmosphere is specified
	top-down, that is, the top level is the first line in the file, the bottom
	(surface) level the last line. All properties refer to model \emph{level} z,
	not to model \emph{layer}. It is important that the correct units are 
	used, otherwise unpredictable results are guaranteed.
	Comments start with \code{\#}. Empty lines are ignored. Please note that there 
	is some redundancy: For air as an ideal gas the density $\rho$, can be
	calculated from pressure and temperature, $\rho = p / kT$. \code{uvspec} will check 
	if this relation is fulfilled and will stop if it is not. 
	{\sl libRadtran} provides the six standard atmospheres by \citet{Anderson1986}:
	\begin{description}
	\item[afglt]   Tropical            (\code{tropics})                                           
       	\item[afglms]  Midlatitude Summer  (\code{midlatitude\_summer})
	\item[afglmw]  Midlatitude Winter  (\code{midlatitude\_winter})
	\item[afglss]  Subarctic Summer    (\code{subarctic\_summer})
	\item[afglsw]  Subarctic Winter    (\code{subarctic\_winter})
	\item[afglus]  U.S. Standard       (\code{US-standard})                   
	\end{description}
	which may be chosen by for example
	\fcode{
	atmosphere\_file tropics
	}
	or by specifying the full file name. These atmosphere files are found in 
	\file{data/atmmod}.
	If no \code{atmosphere\_file} is defined, {\sl uvspec} will automatically select one. 
	If the information  \code{time}, \code{latitude} and \code{longitude} are provided in 
	the input file {\sl uvspec} will choose from the first 5 files, otherwise it takes 
	the U.S. Standard atmosphere.
		
\ifmystic{
\option{bpdf_tsang_u10}
	Wind speed for ocean BPDF (in m/s) at present only available with 
	\code{rte\_solver mystic}.
	\fcode{
	bpdf\_tsang\_u10 value
	}
	The BPDF model has been developed by \citet{Tsang1985}. 
	The wind speed is the most important parameter affecting the ocean reflectance
	matrix. The BPDF model also takes into account shadowing by surface waves. 
	The model has been implemented in a FORTRAN routine by Mishchenko 
	(http://www.giss.nasa.gov/staff/mmishchenko/brf/)
	which has been included into {\sl libRadtran}.
		}
\option{brdf_ambrals}
	AMBRALS (Ross-Li) BRDF, a three-parameter BRDF fit for vegetated and 
	non-vegetated surfaces \citep{roujean1992,wanner97,lucht2000,
        schaaf2002}, based on the RossThickLiSparseReciprocal
        model. The implementation is according to
        \cite{lucht2000}. See also http://www-modis.bu.edu/brdf/models.html.
	\fcode{
	brdf\_ambrals variable value 
	}
	\code{variable} can be one of the following:
	\begin{description}
        \item[iso] f\_iso, Eq. (37) \cite{lucht2000}.
	\item[vol] f\_vol, Eq. (37) \cite{lucht2000}.
	\item[geo] f\_geo, Eq. (37) \cite{lucht2000}.
	\end{description}
        The hot spot correction factor introduced by \cite{maignan2004}
        can be turned on by the option \code{brdf\_ambrals\_hotspot}.
        May be combined with \ifmystic{\code{mystic}},
	\code{disort}, and \code{fdisort2}.
		
\option{brdf_ambrals_file}
	        4 column file, containing the AMBRALS (Ross-Li) BDRF
	parameterization, see \code{brdf\_ambrals}.
	\fcode{
	brdf\_ambrals\_file file
	}
	This option is only supported with
	solvers: \code{disort}\ifmystic{ and \code{rte\_solver
	mystic}}. The columns of the input file are wavelength
	[nm], iso, vol, and geo.  The parameters are interpolated linearly to
	the internal wavelength grid.  To make sure that the results are
	reasonable, specify the Ambrals data on a wavelength grid similar or equal
	to that used internally for the radiative transfer calculation!
		
\option{brdf_ambrals_hotspot}
	Turn on hot spot correction factor by \cite{maignan2004} when using the AMBRALS (Ross-Li) BRDF.
		
\option{brdf_cam}
	Set \citet{cox54a,cox54b} ocean BRDF properties.
	\fcode{
	brdf\_cam variable value
	}
	\code{variable} can be one of the following:
	\begin{description}
	\parameter{pcl} 
	Pigment concentration for \citet{cox54a,cox54b} ocean BRDF (in mg/m$^{-3}$).
	The default value is 0.01 mg/m$^{-3}$. 
	\parameter{sal}
	Salinity for \citet{cox54a,cox54b} ocean BRDF (in "per mille", 0.1\%; this 
        unit is equivalent to the other common units for salinity, 
        ppt - parts per thousand, psu - practical salinity unit).
	The default value is 34.3. 
	\parameter{u10}
	Wind speed for \citet{cox54a,cox54b} ocean BRDF (in m/s).
	The wind speed is the most important parameter affecting ocean BRDF.  The
        minimum allowed wind speed is 1 m/s because otherwise the strong
        specular reflection causes numerical problems. If a lower value is
        specified, the wind speed is automatically set to 1m/s.
	\parameter{uphi}
	Wind direction for \citet{cox54a,cox54b} ocean BRDF.
	Default value is 0 degrees,
        which is wind from the South. 90 degrees corresponds to wind from the West, etc. (Honestly, this was
        never truly validated. It could possibly be that 0 is wind from the North, 90 is wind from the East, etc.)
	The option is only implemented for the \code{mystic} solver. 
        \end{description}

	At present only available with \code{rte\_solver disort}, 
        \ifmystic{\code{rte\_solver mystic}}
        and \code{rte\_solver fdisort2}.  The
        number of streams (\code{number\_of\_streams}) is automatically increased to 16 if
        \code{cox\_and\_munk} BRDF is switched on, to avoid numerical
        problems.
	To switch on Cox and Munk BRDF,
        specify any of the \code{brdf\_cam} options and define at least
        \code{brdf\_cam u10}.
		
\option{brdf_cam_solar_wind}
	Use old definition of wind direction for Monte Carlo simulations. If this switch is set,
	the wind azimuth is identical to the incoming photon azimuth. Else, the wind azimuth is set by
	\code{brdf\_cam uphi} or is 0 by default.
		
\option{brdf_hapke}

	Wavelength-independent Hapke BRDF values.
        Hapke is a parameterization of extraterrestial solid bodies such as
        the moon, asteroids or the inner planets by \cite{hapke1993}.
	\fcode{
	brdf\_hapke variable value
	}
	\code{variable} can be one of the following:
	\begin{description}
	\item[w] Single scattering albedo of surface.
	\item[B0] Hot spot factor.
	\item[h] Surface roughness parameter.
	\end{description}
		
\option{brdf_hapke_file}
	        4 column file, containing the Hapke BDRF
	        parameterization, see \code{brdf\_hapke}.
	\fcode{
	brdf\_hapke\_file file
	}
	This option is only supported with
	solvers: \code{disort}\ifmystic{ and \code{rte\_solver
	mystic}}. The columns of the input file are wavelength
	[nm], w, B0, and h.  The parameters are interpolated linearly to
	the internal wavelength grid.  To make sure that the results are
	reasonable, specify the Hapke data on a wavelength grid similar or equal
	to that used internally for the radiative transfer calculation!
		
\undocumented{
\option{brdf_rossli}
        This is the \cite{Lin2015} implementation of the Ross-Li (AMBRALS) BRDF, see also 
        \code{brdf\_ambrals} BRDF. For unknown reasons, these formulas
        differ in their overall factors of the three terms from the factors in
        \cite{lucht2000}. Identical results can be obtained by using: f\_iso=$\pi$ k\_iso, f\_geo=$\pi$ k\_geo, f\_vol = 4 k\_vol / 3.
	\fcode{
	brdf\_rossli variable value 
	}
	\code{variable} can be one of the following:
	\begin{description}
	\item[iso] k\_iso, Eq. (31) \cite{Lin2015}.
	\item[vol] k\_vol, Eq. (31) \cite{Lin2015}.
	\item[geo] k\_geo, Eq. (31) \cite{Lin2015}.
	\end{description}
        The hot spot correction factor introduced by \cite{maignan2004}
        can be turned on by the option \code{brdf\_rossli\_hotspot}.
        May be combined with \ifmystic{\code{mystic}},
	\code{disort}, and \code{fdisort2}.
		}
\undocumented{
\option{brdf_rossli_file}
        4 column file, containing the Ross-Li (AMBRALS) BDRF
        parameterization, see \code{brdf\_rossli} and \code{brdf\_ambrals\_file}.
		}
\undocumented{
\option{brdf_rossli_hotspot}
	Turn on hot spot correction factor by \cite{maignan2004} when using the Ross-Li (AMBRALS) BRDF.
		}
\option{brdf_rpv}

	Wavelength independent RPV values, see \code{brdf\_rpv\_file}.
	\fcode{
	brdf\_rpv variable value
	}
	\code{variable} can be one of the following:
	\begin{description}
	\item[k] Overwrite the wavelength-dependent
	k value defined in \code{brdf\_rpv\_file}.
	\item[rho0] Overwrite the wavelength-dependent
	rpv0 value defined in \code{brdf\_rpv\_file}.
	\item[theta] Overwrites the wavelength-dependent
	theta value defined in \code{brdf\_rpv\_file}.
	\item[sigma] Constant RPV sigma, to be used for snow \citep{Degunther2000b}.
	\item[t1] Constant RPV t1, to be used for snow \citep{Degunther2000b}.
	\item[t2] Constant RPV t2, to be used for snow \citep{Degunther2000b}.
	\item[scale] Apply a constant scaling factor for the RPV BRDF. 
	Required e.g. if the the albedo should be set to a certain value. This
	factor is only used by \code{rte\_solver disort}, \code{rte\_solver
	fdisort2}\ifmystic{ and \code{rte\_solver mystic}}.
	\end{description}
		
\option{brdf_rpv_file}
	4 to 7 column file, containing the Rahman, Pinty, and Verstraete (RPV) BDRF
	parameterization \citep{Rahman1993} and the snow extension by \citep{Degunther2000b}.
	\fcode{
	brdf\_rpv\_file file
	}
	Bidirectional reflectance distribution functions for a variety of
	surfaces are given in the paper. This option is only supported with
	solvers: \code{disort}, \code{fdisort2}\ifmystic{ and \code{rte\_solver
	mystic}}. The columns of the input file are wavelength
	[nm], rho0, k, and theta.  The parameters are interpolated linearly to
	the internal wavelength grid.  To make sure that the results are
	reasonable, specify the RPV data on a wavelength grid similar or equal
	to that used internally for the radiative transfer calculation!
	Optionally, a fifth column with a constant scaling factor may be
	defined.  If it has seven columns, the fifth to seventh are sigma, t1,
	t2, and if it has eight, the eighth is scale again.
		
\option{brdf_rpv_library}
	The rpv libraries are collections of spectral BRDFs of different surface types, 
	This option must be used either with \code{brdf\_rpv\_type} 
	or \code{surface\_type\_map}, in order to select the specific surface type.
				
	For using a \code{brdf\_rpv\_library} write 
	\fcode{
	brdf\_rpv\_library library\_path 
	}
	where \code{library\_path} is the path of the directory, where the BRDF data is stored. The files are 
	expected to have the names \code{IGBP.01.rpv, IGBP.02.rpv, ...} If \code{brdf\_rpv\_type 1} is specified
	the BRDF from \code{IGBP.01.rpv} will be used, and so on.
	Each file must have the structure like an \code{brdf\_rpv\_file}.
	(This option is quite the same as \code{brdf\_rpv\_file}, except that it offers you an
	easy way to use the option \code{surface\_type\_map} in combinition with your \code{brdf\_rpv\_files}.)
	
	%ak \emph{For our group there is also a built-in library for IGBP surface types.}
	\fcode{
	brdf\_rpv\_library IGBP
	}
	The built-in library contains the first 17 surface types see
	\code{albedo\_library}.  The data is given for the wavelengths 443nm,
	565nm, 670nm, and 865nm. Stay near this wavelength in order to get
	reasonable results. In future this the rpv-library will be NDVI
	dependent, but until now the most common NDVI class is selected
	automatically.
		
\option{brdf_rpv_type}
	With this option the (RPV) BRDF surface type is selected.
	This option can be used with \code{albedo\_library} in order to select a spectral albedo 
	or with \code{brdf\_rpv\_library} in order to select a BRDF function.
	\fcode{
	brdf\_rpv\_type surface\_type\_number
	}
	where \code{surface\_type\_number} is an integer starting from 0,
	where 0 refers to a black surface and the following numbers to 
	the entries in the specified \code{library}.
		
\undocumented{
\option{ck_fu_h2o_continuum}
	Undocumented option to switch the H2O continuum on or off or select a specific 
	version of the continuum.
		}
\option{ck_lowtran_absorption}
	Switch off absorption by individual minor trace gases which are currently only
	included when \code{mol\_abs\_param lowtran} is chosen. The syntax is
	\fcode{
	  ck\_lowtran\_absorption species on/off
	}
	where species may be one of O4, N2, CO, SO2, NH3, NO, HNO3. By default all
	are switched on. 
	
	This option may also be used to turn on/off absorption by O4 in spectral 
	resolution. It is on by default.
		
\option{cloud_fraction_file}
        File containing a cloud fraction profile. 
        \fcode{
        cloud\_fraction\_file file
        }
        Two columns are expected:
        altitude [km] and cloud fraction, including ice and water clouds. If \code{cloud\_fraction\_file}
        is defined, effective cloud properties are calculated assuming either random overlap or maximum random overlap
        of the cloud layers (see also \code{cloud\_overlap}). An example is provided in \file{examples/CF.DAT}.
                
\undocumented{
\option{cloud_fraction_map}
	Undocumented option.
		}
\option{cloud_overlap}
        Cloud overlap assumption. 
        \fcode{
        cloud\_overlap type
        }
        Following types are implemented:
        \begin{description}
        \item[rand] Random overlap of cloud layers 
        \item[maxrand]  Maximum random overlap scheme
        \item[max]      Maximum overlap scheme
        \item[off]      Turn off cloud overlap for ECMWF clouds                                        
        \end{description}
        Per default the \code{cloud\_overlap} scheme is switched off.
                
\option{cloudcover}
        Set the fraction of the horizontal sky area which is covered by clouds.
        \fcode{
        cloudcover typename value
        }
        \code{typename} describes the name of the cloud type, which can be "wc" and "ic":

        \begin{description}
        \item[ic] Cloud cover of ice cloud, where the cloud properties are taken
        from \code{ic\_file}.
        \item[wc] Cloud cover of water cloud, where the cloud properties are taken 
        from \code{wc\_file}.
        \end{description}

        When a cloud cover is specified, the result will be calculated by
        the independent pixel approximation (IPA), that is, as weighted average 
        of cloudless sky and overcast sky. 
        Please note that, if both \code{cloudcover ic} and
        \code{cloudcover wc} are set, both must be equal.
        
        This option is ignored, if the option \code{cloud\_fraction\_file} is used.
                
\undocumented{
\option{cox_and_munk_pcl_map}
	A possibility to specify pigment concentration (in mg/m3) for the Cox and Munk 
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_pcl_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
	
	The default name of the pigment concentration variable is
	\code{chlorophyll}, but can be changed with the optional argument
	\code{variable_name}.  The pigment concentration must be provided as
	function of \code{latitude} and \code{longitude},
	\code{chlorophyll(lat, lon)}, or additionally may also depend on time
	\code{chlorophyll(time, lat, lon)}.  If a time-dependent pigment
	concentration is specified, the pigment concentration will be
	interpolated according to the option \code{time_interpolate}. All
	grids must be provided in the file as \code{double lat(lat)},
	\code{double lon(lon)}, and \code{double time(time)}.
		}
\undocumented{
\option{cox_and_munk_sal_map}
	Specify ocean salinity (in ppt) for the \citet{cox54a,cox54b}
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_pcl_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
	
	The expected name of the pigment concentration variable is per default
	\code{salinity}, but can be changed with the optional argument
	\code{variable_name}.  The pigment concentration must be provided as
	function of \code{latitude} and \code{longitude}, \code{salinity(lat,
	lon)}, or additionally may also depent on time \code{salinity(time,
	lat, lon)}.  If a time-dependent salinity is specified, the salinity
	will be interpolated according to the option
	\code{time_interpolate}. All grids must be provided as \code{double
	lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
		}
\undocumented{
\option{cox_and_munk_u10_map}
	Specify wind speed (in m/s) for the \citet{cox54a,cox54b}
	ocean BRDF with the help of an \emph{netCDF} file, which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	cox_and_munk_u10_map file
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.

	The file must contain the elements of the wind vector \code{U10} and
	\code{V10}.  These must be specified as functions of \code{latitude}
	and \code{longitude} \code{U10(lat, lon)}, \code{V10(lat, lon)}, or
	additionally may also depent on time \code{U10(time, lat, lon)},
	\code{V10(time, lat, lon)}.  If the variable time is present in the
	file, the wind speed will be interpolated according to the option
	\code{time_interpolate}. All grids must be provided as \code{double
	lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}.
		}
\option{crs_file}
	May be used to specify cross sections of O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, 
        O4, SO2, CH4, N2O, CO, or N2 to be used instead of those supplied with 
	{\sl libRadtran}. No temperature dependence may be specified. Use as follows:
	\fcode{
	crs\_file NO2 ../examples/no2\_crs.dat
	}
	The species, e.g. \code{NO2}, must be specified to identify the
	species for which the cross section applies.
	The cross section file has two columns:
	\begin{description}
	\item[1] wavelength (nm)
	\item[2] cross section (cm$^2$)
	\end{description}
		
\option{crs_model}
	Choose between various cross sections.
	\fcode{
	crs\_model species crs
	}
	Following \code{species} are included:
	\begin{description}

	\parameter{rayleigh} Specify the Rayleigh cross section.
	Choose between the following Rayleigh scattering cross sections (\code{crs}):
	\begin{description}
	\item[Bodhaine]		\citet{Bodhaine1999} Rayleigh scattering cross section using their Eqs. 22-23.
	\item[Bodhaine29]	\citet{Bodhaine1999} Rayleigh scattering cross section using their Eq. 29.
	\item[Nicolet] 		\citet{Nicolet1984} Rayleigh scattering cross section.
	\item[Penndorf]		\citet{Penndorf1957} Rayleigh scattering cross section.
	\end{description}
	\citet{Bodhaine1999} is default.

	\parameter{o3} Choose ozone cross section. 
	\code{crs} can be one of the following:
	\begin{description}
	\item[Bass\_and\_Paur]	\citet{Bass1985} ozone cross section.
	\item[Molina] 		\citet{Molina1986} ozone cross section.
	\item[Daumont]		Ozone cross section by \citet{Daumont1992}, \citet{Malicet1995}.
	\item[Bogumil]		Ozone cross section from \citet{Bogumil2003}.
	\end{description}
	\citet{Molina1986} is default.

	\parameter{no2} Choose between the various no2 cross sections.
	\code{crs} is one of:
	\begin{description}
	\item[Burrows]		\citet{Burrows1998} no2 cross section.
	\item[Bogumil]		no2 cross section from \citet{Bogumil2003}.
	\end{description}
	\citet{Burrows1998} is default.

	\end{description}
		
\option{data_files_path}
	The path to the directory where all \code{uvspec} internal data files live, e.g.
	the files that are in the subdirectories of the \file{data} directory
	that comes with the \code{uvspec} distribution. 
	\fcode{
	data\_files\_path path
	}
	The default for \code{path} is \code{../data/}.
		
\option{day_of_year}
	Integer, to correct the calculated radiation quantities for the 
	Sun-Earth distance for the specified Julian day (1-365). 
	\fcode{
	   day_of_year value
	}
	If not specified, the Earth-Sun distance is 1 AU (i.e. equinox
	distance), that is, no correction is applied to the extraterrestrial
	irradiance \code{source solar file}. Alternatively \code{time} may be used
	for that purpose.
		
\option{deltam}
	Turn delta-M scaling on/off. Set to either \code{on} or
	\code{off}. Note that for the \code{rte_solver disort} and
	\code{rte_solver fdisort2} delta-M scaling is hardcoded to be always
	on.
		
\option{disort_intcor}
	Intensity correction method for \code{rte_solver disort} or
	\code{rte_solver fdisort2}. Valid options are \code{phase}, i.e.~the
	phase function is used for the Nakajima intensity correction, and
	\code{moments}, i.e.~the Legendre moments are used for the
	correction. Optionally, the option \code{off} turns off the intensity
	correction. Default is \code{phase}.
		
\undocumented{
\option{disort_spherical_albedo}
	Calculate spherical albedo using \code{disort}. When this option is enabled,
	only the spherical albedo is calculated. The output is enabled by 
	\fcode{output_user spher_alb}
		}
\option{earth_radius}
	Specify the earth radius in km. 
	\fcode{
	    earth_radius value
	}
	This is needed by all solvers in spherical geometry\ifmystic{, e.g. \code{mystic}
	in combination with option \code{mc_spherical}}. The default value is 6370 km.
		
\undocumented{
\option{ECHAM_atmosphere_file}
Undocumented option.
}
\undocumented{
\option{ECHAM_clouds_file}
Undocumented option.
}

\undocumented{
\option{ECHAM_sza}
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
\option{ECMWF_ic_file}
	\fcode{ 
	ECMWF_ic_file file 
	} 
	For further information see \code{ECMWF_wc_file}.
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

\undocumented{
\option{ECMWF_levels_only}
	The atmosphere considered in the simulation has the same height range 
	as the data in the 
	\code{ECMWF_atmosphere_file}/\code{radiosonde}-file. No further levels 
	are added above those.  This option has only an effect in combination 
	with \code{ECMWF_atmosphere_file} or \code{radiosonde} 
	(this option is identical to \code{radiosonde_levels_only}). 
		}
\undocumented{
\option{ECMWF_ozone_climatology}
	The Integrated Forecast System (IFS) of the ECMWF uses a ozone 
	climatology for radiative transfer instead of the ozone simulated by 
	the IFS. If this option is activated the ozone profile of the 
	\code{atmosphere_file} or \code{ECMWF_atmosphere_file} is replaced by 
	the ozone climatology by Fortuin and Langematz (1995). 
	(If there is also a \code{mol_file} for ozone, it modifies the ozone 
	climatology profile.) 
		}
\undocumented{
\option{ECMWF_wc_file}
	Reads in combination with the options \code{latitude}, 
	\code{longitude}, and \code{time} (all mandatory) the pressure, 
	temperature, and cloud liqid water content (CLWC) and cloud cover (CC) 
	from an ECMWF \emph{netCDF} data file. 
	\fcode{ 
	ECMWF_wc_file file 
	} 
	No spatial interpolation of the values is done.  The data nearest to 
	the specified \code{time} will be selected (or linearly interpolated 
	if \code{time_interpolate} is switched on).  In order to use the ECMWF 
	data without cloud overlap assumption, use \code{cloud_overlap off}. 
		}
\undocumented{
\option{ECMWF_wind_file}
	Reads in combination with the options \code{latitude}, 
	\code{longitude}, and \code{time} (all mandatory) the wind components 
	U, V, and W from an ECMWF \emph{netCDF} data file. 
	\fcode{ 
	ECMWF_wind_file file 
	} 
	The data nearest to the specified \code{time} will be selected (or 
	linearly interpolated, if \code{time_interpolate} is switched on). 
		}
\undocumented{
\option{emissivity_map}
	\emph{This option is preliminary and still subject to change (no wavelength dependency yet)!} 
	Specify a wavelength independent emissivity with the help of  
	an \emph{netCDF} file, which is used in combinition  
	with the options \code{latitude}, \code{longitude}, and \code{time}. 
	\fcode{ 
	emissivity_map file [variable_name] 
	} 
	where \code{file} is the location of the \emph{netCDF} file.  
	With the optional argument \code{variable_name} the name of the emissivity variable in  
	the \emph{netCDF} file can be specified. (By default the expected name is \code{EMIS}.) 
	The emissivity must be specified as function of \code{latitude} and \code{longitude} \code{EMIS(lat, lon)},  
	or additionally may also depent on time \code{EMIS(time, lat, lon)}.  
	All grids must be provided as  
	\code{double lat(lat)}, \code{double lon(lon)}, and \code{double time(time)}. 
	{\sl libRadtran} reads the value at the nearest pixel to the given \code{latitude} and  
	\code{longitude}. No spatial interpolation or averaging of the values is done. 
	If the variable time is present in the file, the emissivity data nearest to the specified time  
	will be selected (or interpolated if \code{time_interpolate} is switched on).  
		}
\option{filter_function_file}
	If specified, the calculated spectrum is multiplied with a filter function 
	defined in file.
	\fcode{
	filter\_function\_file file [normalize]
	}
	The file must contain two columns. 
	Column 1 is the wavelength, in nm. Column 2 is the corresponding filter 
	function value. Comments start with \code{\#}. Empty lines are ignored.
	In combination with \code{output\_process sum} or \code{output\_process integrate} this option 
	is useful e.g. to calculate weighted irradiances or actinic fluxes or 
	to simulate broadband or satellite observations.

	If the optional second argument \code{normalize} is specified, 
	the integral of the filter function over wavelength is normalized such
	that \code{output\_process integrate} gives radiative properties per wavelength, averaged over
	the filter function.
			
\option{fluorescence}
	Specifies the magnitude of a bottom surface isotropic fluorescence source.
	\fcode{
	   fluorescence value
	}
	where \code{value} is a number greater or equal to 0.0, 
	constant for all wavelengths.
	Must be used together with  \code{source solar file}. The units of the fluorescence should 
	obviously be the same as for the solar source in \code{source solar file}.
	For wavelength dependent fluorescence use \code{fluorescence\_file}. The default fluorescence is 0.0. 
	Currently only works with the \code{disort} solver.
		
\option{fluorescence_file}
	Location of fluorescence file for wavelength dependent fluorescence emission 
	from the bottom surface. 
	\fcode{
	   fluorescence\_file file
	}
	The file must have two columns.
	Column 1 is the wavelength in nm, and column 2 the corresponding 
	fluorescence. An arbitrary wavelength grid may be chosen as the 
	fluorescence will be interpolated linearely to the wavelength grid used for the 
	radiation calculation. Comments start with \code{\#}. Empty lines are ignored.
	Currently only works with the \code{disort} solver. Furthermore. if \code{raman} is 
	not set, \code{wavelength\_grid\_file} 
	must be specified with the same resolution as the \code{source solar file}, and the first value
	must be the value specified by \code{wavelength}. The units of the fluorescence should 
	obviously be the same as for the \code{source solar file}.
		
\option{heating_rate}
	Calculation of heating rates. Output is only provided at altitudes specified by \code{zout}.
	To get heating rate profiles a number of altitudes must thus be specified by \code{zout}.
	Heating rates is the change of temperature with time in units of K/day.
	For spectral calculations the default output is a matrix: 
	\begin{Verbatim}
		0.0        zout1          zout2 ...
		lambda1    heating\_rates  ...
		lambda2      .
		   .         .
		   .         .
	\end{Verbatim}
	For integrated calculations (\code{output\_process sum} or \code{output\_process integrate}) the default output 
	is in two columns with column 1 being the altitude and column 2 the heating rates.
	The output of \code{heating\_rate} can also be specified with the \code{output\_user} option.
	Note that heating rates are only well-behaved up to altitudes
	for which the respective correlated-k options are valid. E.g. about 60 km for 
	\code{fu} and about 80 km for \code{kato}, \code{kato2}, \code{kato2.96}, and \code{lowtran}.
	Attention: For spectral calculations, the extraterrestrial spectrum is assumed to be in 
	mW/(m2 nm).

	Three different methods are implemented to calculate the heating rate, which can be selected
	with an optional keyword:
	\fcode{
	heating\_rate [method]
	}
	where \code{method} may be either \code{layer\_cd} (heating rates are derived from centered 
	differences of the flux (the default method), \code{local} (heating rates are derived 
	from the actinic flux), or \code{layer\_fd} (heating rates are derived from forward differences of the flux 
        over one layer. Attention: \code{heating\_rate local} introduces new levels into the profile which  
	slightly affects the model output.
	With \code{layer\_fd}, the output is 
	not representative for a \emph{level}, but for the \emph{layer} from the z-level of the line in the output file, 
	where it is written, up to next \emph{output level} above!
		
\option{ic_file}
	Location of file defining ice cloud properties. 
	\fcode{
	ic\_file type file
	}
	
	\code{type} defines the file type, which is identical to \code{wc\_file type}.
	See \code{wc\_file} for choices of \code{type} and a description on the file structures.
		
\option{ic_fu}
	\fcode{
	ic\_fu reff\_def on/off
	}
	Specify wich definition of the effective radius is used.
	
	If \code{on} the parameterization uses the original definition of the effective radius
	as specified in \citet{Fu1996,Fu1998}. 

	Default is \code{off}. The same definition of the effective radius is used as the \citet{Key2002},
	\citet{Yang2000} and \citet{baum05a:_bulk,baum05b:_bulk}
	parameterizations; see discussion of \code{ic\_properties}.

	\fcode{
	ic\_fu deltascaling on/off
	}
	Specify if the \citet{Fu1996} optical properties are delta-scaled or not. 

	If \code{on} delta-scaling is switched on.

	If  \code{off} delta-scaling is switched off. The default is without delta-scaling.
	Please note that this was changed on July 22, 2008: 
	Before, delta-scaling was switched on by default which might have
	caused some confusion, because irradiance calculations were not
	consistent with the other ice cloud parameterizations implemented in
	\code{uvspec}. 
	Using the \citet{Fu1996} parameterization in combination with one of 
	\code{ic\_modify} you now get
	consistent results with all other ice cloud parameterizations.
	
	%THIS IS THE OLD DOCUMENTATION (BEFORE CHANGE OF DELTA-SCALING
	%				      DEFAULT).
	%
	%It has been confirmed that the difference with and 
	%without delta-scaling is typically on the order of less than 1 - 2\% 
	%and therefore it was decided to switch delta-scaling off by default.
	%Please note that this has nothing to do with the internal delta-scaling
	%by the solvers: e.g. disort2 always applies it's own internal delta-scaling 
	%anyway. (THE FOLLOWING NEEDS TO BE MODIFIED).
	%If you define a cloud only by its microphysical properties (ice water
	%content, effective radius), delta-scaling should certainly be switched
	%on and you do not need to read further. 
	%If, however, you want to
	%use the Fu (1996) parameterization in combination with one of 
	%\code{ic_set_tau/tau550/gg/ssa} or \code{ic_scale_gg/ssa} it might be
	%reasonable to switch delta-scaling off and you should make sure that
	%you understand the following. Citing from Fu (1996): "For nonspherical
	%particles in cirrus clouds, it is found that a simple representation
	%of the scattering phase function through the asymmetry factor is
	%inadequate (Fu and Takano 1994). As demonstrated in appendix A, the
	%fraction of scattered energy residing in the forward peak, f, needs to
	%be removed from the scattering parameters to incorporate the strong
	%forward peak contribution in multiple scattering." Or in other words,
	%the sharp forward peak is truncated and added to the unscattered
	%direct radiation. The remaining phase function (excluding the sharp
	%forward peak) can be safely approximated by a Henyey-Greenstein
	%function. The scaling implies a reduction of the optical thickness,
	%the asymmetry parameter, and the single scattering albedo. This
	%reduction can be rather severe, e.g. a factor of about 3 for the
	%optical thickness in the visible spectral range. This implies
	%seemingly inconsistent optical properties: For idential IWC content
	%and effective radius, \code{ic_properties key/yang} would give an
	%(unscaled) optical thickness about three times higher than
	%\code{ic_properties fu}. The effect on the radiation field, however,
	%will be comparable, due the consistent scaling of optical thickness,
	%asymmetry parameter, and single scattering albedo. If you, however,
	%adjust the optical thickness using e.g. \code{ic_set_tau}, the effect
	%on the radiation field will be completely different because the
	%(unscaled) optical thickness by Key (2002) has a completely different
	%meaning as the (scaled) optical thickness by Fu (1996). In such cases
	%it might be reasonable to switch scaling off. This is a complicated
	%and confusing topic and it is suggested that you play around a bit
	%with the options, read the Fu (1996) paper, and make heavy use of the
	%\code{verbose} feature.
		
\option{ic_habit}
	Ice crystal habit for the \citet{Yang2000}, \citet{Key2002} and
	\code{hey} parameterizations, see also  \code{ic\_properties key/yang/hey}.
	\fcode{
	ic\_habit type
	} 
	For Key/Yang \code{type} may be one of \code{solid-column}, 
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-4}, 
	\code{rosette-6}, \code{plate}, \code{droxtal}, \code{dendrite}
	and \code{spheroid}. Please note that this parameterization is only valid for
	a restricted size range, depending on the habit (see table 1 in 
	\citet{Key2002}. Also, some of the habits are only available for
	wavelengths below 5 micrometer (\code{rosette-4}) while others are only available
	for wavelengths larger than 3 micrometer (\code{droxtal}, \code{spheroid}).

	For \code{hey} the following habits can be chosen: \code{solid-column}, 
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-6}, 
	\code{plate}, \code{droxtal}, and the general habit mixture \code{ghm} 
        which follows the ``recipe'' by \citet{baum05a:_bulk}. All
	habits and the habit mixture are available for effective radii from 5 to 90 micrometers in
	the wavelength region from 0.2 to 5 micrometers.

        The parameterization \code{baum\_v36} includes the general habit mixture
        \code{ghm} and the habits \code{solid-column} and \code{aggregate}, 
        all crystals modeled with severe roughness.          
		
\option{ic_habit_yang2013}
	Ice crystal habit for the \citet{yang2013} parameterization,
        \code{ic\_properties yang2013}.
	\fcode{
	ic\_habit type roughness
	} 
	The following habits are available: \code{column\_8elements}, 
        \code{droxtal}, \code{hollow\_bullet\_rosette}, \code{hollow\_column}, 
        \code{plate}, \code{plate\_10elements}, \code{plate\_5elements}, \code{solid\_bullet\_rosette}, 
        and \code{solid\_column}.

        For each habit three degrees of roughness needs to be specified, options are \code{smooth}, 
        \code{moderate} and \code{severe}.         
        
		
\option{ic_modify}
	Modify ice cloud optical properties.
	\fcode{
	ic\_modify variable scale/set value
	}
	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
	
	If you use this option in combination with the ice
	cloud properties by \citet{Fu1996}, please make sure that you understand
	the explanation of \code{ic\_fu}.
		
\option{ic_properties}
	Defines how ice water content and effective particle radius are translated 
	to optical properties. 
	\fcode{
	ic\_properties property [interpolate]
	}
	Possible choices for \code{property} are
	\begin{description}
	\parameter{fu}
	Parameterization by \citet{Fu1996,Fu1998}, see \code{ic\_file}; 
	this is the default setting. Note that this is a parameterization
	which has been created to calculate fluxes but not radiances. 
	Note also that the optical properties in the solar range provided by 
	\citet{Fu1996} are delta-scaled properties (that is, the forward peak of
	the phase function is truncated and optical thickness, asymmetry
	parameter, and single scattering albedo are reduced accordingly),
	whereas \code{uvspec} uses non delta-scaled properties unless the option
	\code{ic\_fu deltascaling on} is specified. By default the parameterization
	by \citet{Fu1996} is treated consistently with all other ice cloud
	parameterizations. 
	For wavelengths up to 4 micrometer \citet{Fu1996} is used while for wavelengths 
	larger than 4 micrometer \citet{Fu1998} is chosen. Please note that 
	\citet{Fu1996} is based on ray-tracing calculations while \citet{Fu1998}
	is a mixture of ray-tracing and Mie calculations (which is required for 
	the infrared wavelengths where the geometrical assumption does not hold).
	Hence, both parameterizations are not fully consistent. Rather, differences
	of some \% are to be expected in the wavelength region where both 
	parameterizations overlap. Also, the wavelength dependence in the solar 
	and infrared parts is treated differently: In the solar part \citep{Fu1996}
	the optical properties are defined for wavelength bands - hence they 
	are assumed constant within each band. In the infrared \citep{Fu1998} 
	they are defined at certain wavelengths and linearely interpolated 
	in between. If you use this option, please see also the
	discussion of \code{ic\_fu deltascaling} and \code{ic\_fu reff\_def}.
	The allowed range for the effective radius is from 9.315 -  65.120 micrometer. 
	\parameter{echam4}
	Use the simple two-band parameterization of the ECHAM4 climate model, described 
	in \cite{Roeckner1996}; this is probably only meaningful if you want to compare
	your results with ECHAM4, the two bands are 0.2 - 0.68 micrometer and 0.68 - 4.0 micrometer.
	Within the two ECHAM4 bands, the optical properties are assumed constant.
	\parameter{key}
	Parameterization by \citet{Key2002}. This parameterization can also
	be used to calculate radiances because it uses a
	double-Henyey-Greenstein phase function which better represents both
	forward and backward peaks. This parameterization covers the wavelength region  
	from 0.2 to 5.0 micrometer and is available for the following \code{habit}:
	\code{solid-column}, \code{hollow-column}, \code{rough-aggregate}, \code{rosette-4}, 
	\code{rosette-6}, and \code{plate}.
	\parameter{yang}
	Parameterization similar to \citet{Key2002} but based on more recent
	single scattering calculations. Below 3.4 micrometer it actually equals
	the \citet{Key2002} parameterization while from 3.4 - 100 micrometer new 
	coefficients have been calculated with much higher wavelength
	resolution and better accuracy. Hence, \code{yang} should give a reasonably
	consistent approximation from 0.2 - 100 micrometer, suitable for spectrally
	resolved calculations of radiance and irradiance.
	The covered range for the effective radius depends on the \code{ic\_habit}.
	(In micrometer: \code{solid-column} [5.96, 84.22], \code{hollow-column} [4.97, 70.24],
	 \code{rough-aggregate} [3.55, 108.10], \code{rosettes-4} [2.77, 45.30], 
	\code{rosettes-6} [2.85, 46.01], \code{plate} [4.87, 48.18], \code{dendrites}
	 [0.45, 1.88], \code{droxtal} [9.48, 293.32], 
	\code{spheroid}  [6.58, 203.39]). 
	\parameter{baum}
	Use ice cloud parameterization from \citet{baum05a:_bulk,baum05b:_bulk},
	\url{http://www.ssec.wisc.edu/\~baum/Cirrus/IceCloudModels.html}.
	In combination with the radiative transfer solvers \code{disort}, \ifmystic{\code{montecarlo},}
	 and \code{fdisort2}, accurate phase functions are used.
        \parameter{baum_v36} Use cloud parameterization from
	\citet{heymsfield2013,yang2013,baum2014} covering the 
	spectral range from 0.2 to 100~$\mu$m and effective radii from
	5 to 60 $\mu$m. The parameterization
	assumes severly roughened and randomly oriented ice particles
	and includes the full phase matrix. Three set of models are
	available, they can be selected using \code{ic\_habit}:
	\code{ghm} is based on a general habit mixture involving 9
	habits, \code{solid-column} assumes severely roughened solid
	columns, and \code{rough-aggregate} is based on severly
        roughened aggregates. The default is \code{solid-column}. 
	%\parameter{baum\_hufit }
	%Similar to the option \code{baum} but here the phase function
	%is parameterized by 128 Legendre coefficients, calculated with the
	%delta-fit method from 
	%\citet{Hu2000}. This parameterization covers the region  
	%from 0.4 to 2.2 micrometer. If high accuracy is needed e.g. in the vicinity of the halo, 
	%the forward peak, or the backscatter peak, \code{ic\_properties baum}
	%is recommended.
	\parameter{hey}
	Use pre-calculated ice cloud optical properties including full phase
	matrices.
	The parameterization
	is currently only available for the spectral region from 0.2 to 5
	micrometers. The single scattering
	properties have been been generated by Hong Gang using the
	models by \citet{Yang2000}. The parameterization is based on
	simple gamma distributions 
	\begin{equation}
	n(r) = n_0 r^{\alpha} \exp\left(-\frac{(\alpha+3)r}{r_e}\right),
	\end{equation}
	where $n_0$ is found by normalization and $\alpha$ is set to 1. In case of
	spherical particles the parameter $r_e$ would be the effective
	radius. For aspherical particles, the parameter $r_e$ is found
	iteratively so that the size distribution yields the required
	effective radius. The
	parameterization is availabe for the following habits: \code{solid-column},
	\code{hollow-column}, \code{rough-aggregate}, \code{rosette-6},
	\code{plate}, and \code{droxtal}. The
	default habit is \code{solid-column}. Furthermore a general
        habit mixture \code{ghm} similar to the one defined by
	\citet{baum05b:_bulk} may be selected.  For the HEY
	parameterization the ice crystals are assumed to be smooth, in
	contrast to the severely roughened particles assumed by
	\code{baum_v36}. The habit can be specified using the option
	\code{ic\_habit}.
        \parameter{yang2013} Pre-calculated ice optical properties
	including full phase matrices for the spectral region from 0.2
	to 99 $\mu$m. The parameterization uses single scattering data by
	\citet{yang2013} and assumes gamma size distributions as the
	\code{hey} parameterization. It is available for 9 habits and
	3 degrees of roughness, which can be selected using the option
	\code{ic\_habit\_yang2013}. 
	%\parameter{ic-mie}
	%Use pre-calculated Mie tables; useful for \code{mol\_abs\_param}; 
	%the tables are expected in \code{data\_files\_path}\file{/correlated\_k/.../}.
	%For spectral or pseudo-spectral (\code{mol\_abs\_param sbdart}) calculations, 
	%a set of pre-calculated tables is also available; the wavelength grid points of these 
	%data has been carefully selected such that the extinction cross section, 
	%single scattering albedo, and the asymmetry parameter are accurate to 1\% 
	%(compared to the fully-resolved Mie calculation) for all wavelengths
	%between 250nm and 100 micrometer. 
	%For spectral or pseudo-spectral
	%calculations the optional argument \code{interpolate} has to be defined explicitely to
	%initiate the interpolation of the optical properties to the internal wavelength grid.
	%Please note that this option may be extremely memory-consuming because for each 
	%internal wavelength a full set of Legendre moments of the phase function is 
	%stored (up to several thousands). The Mie tables are not part of the standard distribution 
	%(because of their large size) but they are freely available from http://www.libradtran.org. 
	%Note that a Mie calculation assumes spherical ice particles, the scattering function of 
	%which differs systematically from non-spherical particles. Hence, \code{ic\_properties mie}
	%is usually not representative of natural ice clouds.
	\parameter{filename}
	Read optical properties from specified filename; file format is as produced 
	by the \code{mie} tool of {\sl libRadtran} (see \code{output\_user netcdf}). 
	% or by Frank Evans' \code{cloudprp}.
	\end{description}
	The default property is \code{fu}.

	Please note also that, in contrast to spherical particles, there is no unique 
	definition of effective size for non-spherical particles. In particular, the
	above parameterizations use different definitions which, however, differ only by 
	a constant factor. 
	\citet{Yang2000}, \cite{Key2002}, and \citet{baum05a:_bulk,baum05b:_bulk}
	 use the general definition 
	\begin{equation}
	   r_{\rm eff} = {{3}\over{4}}{{\int V(h) n(h) dh}\over{\int A(h) n(h) dh}}
	\end{equation}
	  where $h$ is the maximum dimension of an ice crystal, $n(h)$ is the
	  number of particles with maximum dimension $h$ in the size distribution,
	  and $V$ and $A$ are the volume and mean
	  projected area of the particles, respectively. The volume and area are
	  based on the spherical diameter with equivalent volume and the
	  spherical diameter with equivalent projected area as defined by 
	  \citet{Yang2000}. On the other hand, \citet{Fu1996,Fu1998} use
	  hexagonal columns and use the following definition
	\begin{equation}
	   r_{\rm eff} =  {{\int D^2 L n(L) dL}\over{2 \int (D L + {\sqrt{3}\over{4}} D^2) n(L) dL}}
	\end{equation}
	  where $D$ is the width of the ice crystal (that is, the maximum diameter of the 
	  hexagonal area) and $L$ is the length. The integrand in the numerator is proportional 
	  to the volume while that in the denominator is proportional to the projected area. 
	  Evaluating these formulas one finds that, for the same hexagonal particle, the effective
	  radius would be $3 \sqrt{3} / 4 = 1.299$ times larger following the 
	  \citet{Yang2000}, \citet{Key2002} definition rather than the \citet{Fu1996,Fu1998} definition.
	  As an example, an effective radius of 20$\mu m$ with 
	  \code{ic\_properties fu} and \code{ic\_fu reff\_def on} and
	  1.299 $\cdot$ 20$\mu m$ = 26$\mu m$ with \code{ic\_properties yang} would give comparable results
	  for hexagonal columns. 
	To use the original definition of the effective radius by \citet{Fu1996,Fu1998} use
	\code{ic\_fu reff\_def on}!

	With the optional argument \code{interpolate} the ice cloud optical properties are interpolated 
	over wavelength; useful for precalculated optical property files defined with \code{ic\_properties}. 
	Please note that this option may be extremely memory-consuming because for each internal wavelength 
	a full set of Legendre moments of the phase function is stored (up to several thousands). 
		
\option{include}
	Include a file into the \code{uvspec} input.  
	\fcode{ 
	include file 
	} 
	Works exactly like the C \code{\#include} or  
	the Fortran \code{INCLUDE} statements. 
		
\option{interpret_as_level}
	Interpret profile properties as level properties (this was the default
	behaviour before version 1.4). 
	\fcode{
	interpret\_as\_level profile
	}
	profile can be either \code{wc}, \code{ic} or any profile type specified in \code{profile\_file}.

	If \code{interpret\_as\_level wc} is defined, a \code{wc\_file} would be interpreted as
	follows:
	\begin{verbatim}
	#      z     LWC    R\_eff
	#     (km)  (g/m^3) (um)  
	       5.000    0      0   
	       4.000   0.2   12.0 
	       3.000   0.1   10.0 
	       2.000   0.1    8.0 
	\end{verbatim}
	The value 0.2 g/m$^3$ refers to altitude 4.0km, as e.g. in a
	radiosonde profile. The properties of each layer are calculated as
	average over the adjacent levels. E.g. the single scattering
	properties for the model layer between 3 and 4km are obtained by
	averaging over the two levels 3km and 4km. To allow easy definition of
	sharp cloud boundaries, clouds are only formed if both liquid water
	contents above and below the respective layer are larger than
	0. Hence, in the above example, the layers between 2 and 3 as well as
	between 3 and 4km are cloudy while those between 1 and 2km and between
	4 and 5km are not.

	Note that since version 1.4 the default is to interpret profile properties as 
	layer properties. For example wc properties are assumed to be constant over the layer. 
	The layer reaches from the level, where the properties are defined in the
	\code{wc\_file} to the level above that one.  The following lines
	\begin{verbatim}
	#      z     LWC    R_eff
	#     (km)  (g/m^3) (um) 
	       4.000   0.0   0.0
	       3.000   1.0  10.0
	\end{verbatim}
	define a cloud in the layer between 3 and 4 km with sharp boundaries.
		
\option{isotropic_source_toa}
	Specifies that isotropic illumination is used at top-boundary instead
	of beam source. Useful for those who want to calculate the reflectance 
	for a homogeneous or inhomogeneous atmosphere. The intensity is still set by
	\code{source solar file}. Only works with \code{disort} and \code{twostr}.
		
\option{latitude}
	This option can be used to specify the latitude of the location to simulate.
	(This option only has an effects, if \code{longitude} is specified, too.)
	\fcode{
	latitude N/S deg [min] [sec]
	}
	where \code{N} and \code{S} stands for the northern and southern hemisphere, respectively.
	\code{deg min sec} is the position in degrees, arc minutes, and arc seconds.
	\code{deg} might also be a float number. \code{min} and \code{sec} may be obmitted.
	The \code{latitude} information will be used for the following: 
	
	%\code{latitude} in combination with \code{longitude}, \code{time}, and any \code{map}-option
	%is used to select the location where to read the input data.
	
	\code{latitude} in combination with \code{longitude} and \code{time} is used to calculate 
	the solar zenith angle, if no \code{sza} is specified (see also \code{time_interval}).
	
	\code{latitude} in combination with \code{longitude} and \code{time} is used to choose a
	suitable default atmosphere file, if no \code{atmosphere_file} is specified. 
		
\option{longitude}
	This option can be used to specify the longitude of the location to simulate.
	(This option only has an effects, if \code{latitude} is specified, too.)
	\fcode{
	longitude E/W deg [min] [sec]
	}
	where \code{E} and \code{W} stand for the eastern and western hemisphere, respectively.
	\code{deg min sec} is the position in degrees, arc minutes, and arc seconds.
	\code{deg} might also be a float number. \code{min} and \code{sec} may be obmitted.
	For possible usage of the \code{longitude} information, see \code{latitude}. 
		
\iflidar{
\option{mc_aerosol_is}
	Enable aerosol concentration importance sampling. The input file
	includes one column containing the desired scaling factors: 
	\fcode{
	mc_aerosol_is file
	}
	USE ONLY IF YOU REALLY KNOW WHAT YOU ARE DOING!
  		}
\iflidar{
\option{mc_albedo_spectral}
	To be done!
		}
\ifthreedmystic{
\option{mc_albedo_spectral_file}
	Define a MYSTIC 2D spectral albedo input file.
	\fcode{
	mc\_albedo\_spectral\_file file
	}
		}
\ifmystic{
\option{mc_azimuth_old}
	Use old MYSTIC azimuth convention (0 degree = looking from the
	direction of the sun; 180 degree = looking into the direction of the
	sun; that is, exactly opposite to the disort convention). The MYSTIC
	azimuth was changed March 1, 2004 - hence this option was introduced
	for compatibility reasons.
		}
\ifmystic{
\option{mc_backward}
	Backward tracing of photons. \code{mc_backward} takes either zero, two or four coordinates:
	\fcode{
	  mc_backward [ix_start iy_start] [ix_end iy_end]
	}
	where \code{ix_start}, \code{iy_start} is the index of the sample pixel to be calculated
	or the pixel area from \code{ix_start} to \code{ix_end} and \code{iy_start} to \code{iy_end}.
	All x-indices must be in the range of 0 ... (\code{Nx}-1) and y-indices the range of 0 ... (\code{Ny}-1).
	If no coordinates are specified, all sample pixels will be calculated.
	\code{mc_backward} computes radiances and downward diffuse irradiances. If 
	a different quantity is required, please use \code{mc_backward_output}.
       		}
\ifmystic{
\option{mc_backward_output}
	Specify quantity to be calculated using backward Monte Carlo. 
	\fcode{
	mc\_backward\_output output [unit]
	}
	So far the
	following \code{output} options have been implemented: 
	\begin{description}
	\parameter{edir} 
	  direct horizontal irradiance
	\parameter{edn} 
	  diffuse downward irradiance (default)
	\parameter{eup} 
	  diffuse upward irradiance
	\parameter{act} 
	  actinic flux
	\parameter{abs} 
	  absorption
	\parameter{emis} 
	  emission
	\parameter{heat} 
	  heating rates, that is absorption + emission
	\end{description}
	For \code{abs}, \code{emis}, \code{heat} an optional argument \code{W\_per\_m2\_and\_dz} (default), 
	\code{W\_per\_m3}, or \code{K\_per\_day} may be specified which converts the result e.g. to heating rates. 
		}
\iflidar{
\option{mc_backward_sunshape_file}
	Path to a 2-coloumned text file containing the extraterrestrial sun shape. 1. column:
	A probability density or if you want to a radiance distribution. Does not need
	to be normalized 2. column: Relative angular distance to center from sun inside sun disc,
	values from 0 to 1 (0: center, 1: limb).  Use \code{mc_sun_angular_size} to manually 
	set the extent of the sun. If filename is set to "default", or if only sun radius is given
	 a spectrally resolved sun shape as in \citet{koepke2001} is used.
		}
\ifmystic{
\option{mc_backward_writeback}
	If set, the distribution of photons contributing to the result is written 
	to a file with extension .bac which may be useful for some interpretations
	(it basically tells you where the photons come from whichn contribute to the result). 
		}
\ifmystic{
\option{mc_basename}
	Filename for MYSTIC 3D output (default: mc).
	\fcode{
	mc\_basename basename
	}
		}
\undocumented{
\option{mc_bcond}
	Undocumented option to define MYSTIC boundary conditions: periodic (default),
	absorbing (\code{absorb}), mirror (not really supported - don't expect meaningful results). 
		}
\ifmystic{
\option{mc_coherent_backscatter}
	Switches on coherent backscattering, use only with \code{mc_polarisation}.
        \begin{itemize}
        \item If used without \code{mc_lidar} or \code{mc_radar} it will calculate enhancement values directly 
        from averaged scattering matrices (\citet{mishchenko2006}) and output them 
        into a .mish.cb file. This method is only correct in the exact backscattering
        direction, its accuracy is around $5\%$ compared to \citet{Mishchenko1992}.
        \item If used with \code{mc_lidar} or \code{mc_radar} it will use the Stokes Vector method which explicitely
        simulates backward photon paths, calculate the enhancement using the Stokes vector-method
        (\citet{Muinonen2004}) and enhance the result accordingly.
        To get enhancement values repeat the calculation without coherent backscattering and
        calculate the quotients from the results. This method is not limited to the exact
        backscattering direction, but the accuracy especially in backscattering direction can be
        problematic, so it should be viewed as an
        experimental option. Additional output like the value of $0.75/(kl^*)$ which gives the cone 
        halfwidth in radian (\citet{fiebig2010_phd}) is written into a .lid.cb file.
        \end{itemize}
		}
\ifmystic{
\option{mc_escape}
	Calculate MYSTIC radiances via escape probabilities; slows down the tracing 
	but usually speeds up the computation considerably since it reduces
	noise. Switched on per default since it should basically be used always 
	when calculating radiances. Only meaningful with \code{rte_solver montecarlo}. The syntax is
	\fcode{
	  mc_escape on/off
	} 
  		}
\ifthreedmystic{
\option{mc_ipa}
	Run MYSTIC in independent pixel mode. Only meaningful with \code{rte_solver montecarlo}.
		}
\iflidar{
\option{mc_jacobian}
	Calculate jacobi matrix. Only meaningful with \code{rte\_solver 
	montecarlo}. 
 
	This option works with \code{mc\_lidar}, or with \code{mc\_backward}. In 
	the first case, it delivers the jacobi matrix for an elastic or HSRL 
	lidar. In the latter case, it delivers the jacobi matrix for height 
	resolved changes in molecular absorption and scattering coefficients. 
		}
\iflidar{
\option{mc_lidar_file}
	File containing positions, looking directions, and opening angles of 
	lasers and detectors for lidar simulations in MYSTIC.  Only meaningful 
	with \code{mc_lidar}. 
	\fcode{ 
	mc_lidar_file file 
	}
		}
\ifthreedmystic{
\option{mc_maxscatters}
	If set, photons are destroyed after n scatters. Please note that this
	is only for testing, debugging and process understanding! The result
	will certainly be wrong because photons are lost.
		}
\ifthreedmystic{
\option{mc_panorama}
	Use this option only together with \code{mc_panorama_view}.

    \fcode{
    mc_panorama circumsolar_var_red
    }
    Variance reduction method for the calculation of circumsolar radiation.
    Re-distributes the initial backward photon directions according to DDIS phase\_max.
    Use only if you calculate the average radiance for one large FOV (at least two times of the angular radius of the sun disk, e.g. 
    if your sensor convers 0$^\circ$ to 1$^\circ$ from the sun center or 0$^\circ$ to 5$^\circ$ from the sun
    center). For smaller FOVs (e.g 0$^\circ$ to 0.35$^\circ$) it may slow down convergence as this method does not yet consider
    the extent of the sun disk. Results should not be wrong but you may need more photons to achieve the same 
    standard deviation. Consideration of the sun's extent in this method should be easy to implement, however.
    For the simulations of radiance profiles (so called sunshapes) this option should probably not be used as it will
    only increase computing time without any benefit. For further information see Reinhardt (2013) Sect. 3.4.2 
    (Dissertation, LMU, urn:nbn:de:bvb:19-164380, http://edoc.ub.uni-muenchen.de/16438/).

	\fcode{
	mc_panorama distr_photons_over_pixel
	}
	With only \code{mc_panorama_view} all photons for the same pixel
	are emitted into exactly the same direction. The finite solid angle that
	is covered by one pixel is not accounted for. With 
	\code{mc_panorama distr_photons_over_pixel} the photons are distributed
	over the FOV covered by the pixel. 
	
	\fcode{
	mc_panorama no_pixel
	}
	Useful for calculation of radiance distributions which shall be exactly
	charted. Originally \code{mc_panorama_view} was meant to calculate images
	like produced by a CCD camera. You set the edges of your desired image
	and the resolution of the sensor. MYSTIC will then calculate the
	radiances at the center of each pixel. This implies that values for
	the edges of the image as given by the input file are never
	computed. The border angles for which radiances are computed are
	shifted half a pixel to the inside of the domain. With
	\code{mc_panorama no_pixel} the edging values as given in the input
	file are directly hit. A short example: Assume you only have three
	pixels. You choose your edges to be at $\phi=0$ and $\phi=180$. With
	\code{mc_panorama no_pixel} MYSTIC will calculate radiances at
	$\phi=0$, $\phi=90$ and $\phi=180$. Without this option you will get
	values for $\phi=30$, $\phi=90$ and $\phi=150$.
	
	\fcode{
	mc_panorama quicklook
	}
	This option allows to do a fast calculation in order to see whether you
	have set your parameters right. All panorama pixels are run in a fast
	mode: only one photon is started per pixel, it is moved till to a
	fixed optical depth (0.1), and then a local estimate is performed. In
	case the local estimate is done scattering (i.e.~the photon has not
	hit a surface before reaching tau=0.1), isotropic scattering is
	assumed. Then the photon is killed. We suggest you turn off Rayleigh
	scattering and aerosols while using this option. The result lets you
	see roughly where your panorama is looking at surface, sunlit clouds,
	shaded clouds, or slanted clouds.
	
	\fcode{
	mc_panorama with_direct_rad
	} 
	Output is direct plus diffuse radiance.
	Use together with \code{mc_panorama distr_photons_over_pixel}.

	\fcode{
	mc_panorama weight_with_cos
	} 
	Output is weighted with the cosine of the angle between photon and zenith of camera.
	Use together with \code{mc_panorama distr_photons_over_pixel}.
		}
\ifthreedmystic{
\option{mc_panorama_view}
	Simulate a panorama. Only works with \code{mc_sensorposition},
	\code{mc_sample_grid} and \code{mc_backward}. Can best be explained by
	following example:

	\fcode{\\
	mc_panorama_view PHI_1 PHI_2 THETA_1 THETA_2\\
	mc_sensorposition X Y Z \\
	mc_sample_grid NPHI NTHETA [] []\\
	mc_backward NPHI1 NTHETA1 NPHI2 NTHETA2\\
	}

	or, with numbers:

	\fcode{\\
	mc_panorama_view 45 135 80 180\\
	mc_sensorposition 3000 3000 26 \\
	mc_sample_grid 90 100 [] []\\
	mc_backward 10 0 10 99\\
	}

	In this example, \code{mc_panorama_view} defines the camera to cover
	the area between 10 degrees below the horizon (80 degrees from nadir,
	looking into the ground) and the zenith (180 degrees from nadir), and
	looking into the directions between south-west (45) and
	north-west (135). The camera is positioned
	at x=3km, y=3km and z=26m, see \code{mc_sensorposition}. The number of
	pixels that the camera has is defined by \code{mc_sample_grid}, hence
	the camera horizontal resolution is 1 degree in the
	horizon: (135-45)/90. The camera vertical resolution is also 1
	degree: (180-80)/100. This run only calculates one vertical line at
	phi=55.5 (see \code{mc_backward}). The square brackets in the line
	\code{mc_sample_grid} are usually dx and dy, they are set
	automatically in the presence of 3d clouds. If no clouds are present,
	you need to set them, if you use \code{mc_elevation_file}, you need to
	set them consistent with the elevation file. NOTE: If you have
	formerly used \code{mc_panorama}, this is what has changed: The order
	of the numbers has changed: Formerly, we wrote \code{mc_panorama MU_1
	MU_2 PHI_1 PHI_2}, where MU\_i=cos(THETA\_i). In the options
	\code{mc_sample_grid} and \code{mc_backward}, X and Y have been
	swapped. Although these changes are far from being nice, they are
	necessary to introduce a conform nomenclatur throughout
	libRadtran. Sorry for that. (CAUTION! The definition of phi for
	panoramas was formerly different than that of \code{phi}!  This has
	been changed on 19.11.2010!)
		}
\ifmystic{
\option{mc_photons}
	Total number of photons to be traced by the Monte Carlo solver, MYSTIC.
	\fcode{
	mc_photons value
	}
	Only meaningful with \code{rte_solver montecarlo}.
		}
\ifmystic{
\option{mc_photons_file}
	Distribution of photons over wavelength bands; to be used with \code{mol_abs_param}.
	\fcode{
	mc_photons_file file
	}
	For an example see \file{data/correlated_k/kato2/x_solar.dat}. No error checking!
	Do only use if you are absolutely sure what you are doing.
	Only meaningful with \code{rte_solver montecarlo}.
		}
\ifmystic{
\option{mc_polarisation}
	Switch on polarisation for \code{rte_solver montecarlo}. You can provide one of the 
        following optional numbers to specify the initial Stokes Vector:
        \begin{itemize}
        \item 0 : (1,0,0,0) (default)
        \item 1 : (1,1,0,0)
        \item 2 : (1,0,1,0)
        \item 3 : (1,0,0,1)
        \item -1 : (1,-1,0,0)
        \item -2 : (1,0,-1,0)
        \item -3 : (1,0,0,-1)
        \item 4 : Each Stokes Vector is randomly determined fulfilling the condition $I^2=Q^2+U^2+V^2$,
        use this e.g. if you need unpolarized radiation with \code{mc_coherent_backscatter} in combination
        with \code{mc_lidar}.
        \end{itemize}
	Details about the implementation of polarisation are described in \citet{emde2010}.
		}
\ifmystic{
\option{mc_rad_alpha}
	Define opening angle for radiance calculations without local
	estimate. Ths oprion is useful for all-sky simulations.	
                }
\ifthreedmystic{
\option{mc_radial_pathlength_dt}
	Specify time increment for \code{mc_radial_pathlength}. Time is converted to pathlength 
	assuming a speed of light of 3$\cdot$10$^8$~m/s
		}
\option{mc_randomseed}
       	Provide your own random seed (positive integer) for the random number generator. 
       	\fcode{
       	mc_randomseed value
       	}
       	Usually a random seed is determined from current time plus process id. 
       	This option is useful to re-run a simulation for debugging.
       		
\undocumented{
\option{mc_refraction}
	Enable refraction for \code{rte_solver montecarlo}. Works only 
	in 1D spherical geometry (with option \code{mc_spherical}).
		}
\iflidar{
\option{mc_ris}
        Use with \code{mc_radar} or \code{mc_lidar} to get good statistics. Artificially enhances
        the scattering optical depth by a constant factor.
        \fcode{
        mc_ris factor value
        }
        Set the ris-factor manually.
        \fcode{
        mc_ris optical_depth value
        }        Automatically determine a ris-factor to reach the provided scattering optical
        depth value along a straight line from the transmitter to the border of the domain.
        In case of Lidar/Radar four additional paths at the rim of the transmitter cone are
        simulated and if needed an average over all paths is calculated.
		}
\iflidar{
\option{mc_satellite_view}
	Simulate a satellite image. Only works with
	\code{mc_satellite_position}/\code{mc_sensorposition},
	\code{mc_sample_grid}, \code{mc_backward} and
	\code{mc_spherical 3D}. Can best be explained by following example:
	
	\fcode{\\
	mc_spherical 3D \\
	mc_satellite_view PHI_W PHI_E THETA_S THETA_N \\
	mc_sample_grid N_PHI N_THETA \\
	mc_satellite_position SAT_LON SAT_LAT SAT_ALT \\
	mc_sun_position SUN_LON SUN_LAT \\
	mc_backward \\
	}
	
	In this example, the earth is a sphere (not plane parallel). The
	satellite is in zenith above longitude SAT\_LON degrees, latitude
	SAT\_LAT degrees, and at an altitude of SAT\_ALT
	meters. (Alternatively, you can use \code{mc_sensorposition S_X S_Y
	S_Z}, with [S\_X,0,0] being toward MSG, [0,S\_Y,0] being above the
	Maledives, and [0,0,S\_Z] being above the north pole.) The image taken
	has a resolution of N\_THETA times N\_PHI, and the picture (taken with
	a camera with spherical picture geometry) spans from PHI\_W to PHI\_E
	degrees (positive being east of nadir) and from THETA\_S to THETA\_N
	degrees (positive being north of nadir). In this example, the sun is
	in zenith at the latitude SUN\_LAT degrees (positive meaning north) and
	longitude SUN\_LON degrees (positive meaning east). See also
	\code{mc_panorama_view}, which is a similar option.  Example for MSG:
	
	\fcode{\\
	mc_spherical 3D \\
	mc_satellite_view -10 10 -10 10 \\
	mc_sample_grid 4000 4000 \\
	mc_satellite_position 0 0 3.6e7 \\
	mc_sun_position -10 23\# about 1PM UTC on mid summer \\
	mc_backward \\
	}
	
	Only meaningful with \code{rte_solver montecarlo}.
		}
\ifthreedmystic{
\option{mc_sensorposition}
	Define the position of a sensor. 
	\fcode{
	mc_sensorposition x-value y-value z-value
	}
	Has been introduced for irradiance calculations in topography
	and might not properly work with all options. 
		}
\ifmystic{
\option{mc_spherical}
	Spherical geometry in MYSTIC. 
	\fcode{
	mc_spherical 1D
	}
	Works only in "1D" - \code{wc_file 3D} 
	and \code{ic_file 3D} are not yet considered. If \code{mc_spherical} is selected
	\code{mc_backward} is switched on automatically. Viewing
	direction (\code{umu}, \code{phi}) and sun position
	(\code{sza}, \code{phi0}) are defined with respect to the
	sensor position specified by \code{zout}. For details about
	the implementation of spherical geometry please refer to
	\citet{emde2007}. 

	\iflidar{
	\fcode{
	mc_spherical 3D
	}
	3D spherical geometry in MYSTIC. Only works with
	\code{mc_satellite_position}/\code{mc_sensorposition}, and
	\code{mc_satellite_view}, and only for solar backward. You have to set
	the sun position using \code{mc_sun_position}. Still very
	experimental. Does not work with \code{mc_elevation_file}.
	}
		}
\ifthreedmystic{
\option{mc_std}
	Calculate standard deviation of the average.
		}
\ifthreedmystic{
\option{mc_sun_angular_size}
	At the moment only useful together with \code{mc_panorama_view}.
	Set the angular radius of the sun in degrees. If omitted the radius is calculated 
	(more or less correctly) via the earth-sun distance (not well tested).
	If no \code{mc_backward_sun_shape_file} is given a spectral sun shape according to
	Koepke2001 is used.
		}
\iflidar{
\option{mc_sun_position}
	Sun position when using \code{mc_spherical 3D}. See \code{mc_satellite_view}
	for an example and explanation.
	\fcode{
	mc_sun_position longitude latitude
	}
		}
\ifmystic{
\option{mc_surface_reflectalways}
	Usually, a photon is either absorbed or reflected at the surface, 
	with a probability defined by the surface albedo. If 
	\code{mc_surface_reflectalways} is specified, each photon is reflected and 
	the albedo is considered by reducing the photon weight. In case of BRDF,
	\code{mc_surface_reflectalways} is switched on automatically because the 
	other method is no longer implemented for non-Lambertian BRDFs, due to 
	implementation and numerical problems. For small albedos, the computational 
	time is increased if \code{mc_surface_reflectalways} is used; however, 
	the accuracy of the upward radiance (reflected by the surface) is increased 
	considerably. In case of clouds, however, computational time might be 
	increased considerably without gaining accuracy. 
       		}
\undocumented{
\option{mc_tipa}
	Run MYSTIC in tilted independent pixel mode. Only meaningful with \code{rte_solver
	montecarlo}. Choose either \code{mc_tipa dir3d} or \code{mc_tipa
	dirdiff}. The latter is the analogon to \code{tipa dirdiff}, where the
	3d-cloud-matrix is tilted according to the solar zenith angle and is
	then used as input for all the following calculations. 
	In \code{mc_tipa dir3d}, photons are propagated in 3d until the first
        scattering occurs, after that the photon is captured within the current
        column.
		}
\iflidar{
\option{mc_truncate}
	Truncate phase function at the specified polar angle mu. 
	USE ONLY IF YOU REALLY KNOW WHAT YOU ARE DOING!
       	 	}
\ifmystic{
\option{mc_vroom}
	Variance Reduction Optimal Options Method (VROOM). Options are "on"
	and "off". Needs to be specified if you are calculating radiances and
	spiky phase functions are present in your atmosphere. If you are using
	VROOM, please cite: \cite{buras2011a}.
       		}
\option{mixing_ratio}
	Mixing ratio in ppm.
	\fcode{
	mixing\_ratio species value
	}
	\code{species} can be one of the following:
	\begin{description}
	\item[O2] The mixing ratio of O2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
	\item[H2O] The mixing ratio of H2O in ppm. Scale the profile so that the mixing
        ratio at the user-define \code{altitude} assumes the specified value.
        \item[CO2] The mixing ratio of CO2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
	\item[NO2] The mixing ratio of NO2 in ppm. Scale the profile so that the mixing
        ratio at the user-defined \code{altitude} assumes the specified value.
        \item[CH4] The mixing ratio of CH4 in ppm (default: 1.6 ppm). 
	\item[N2O] The mixing ratio of N2O in ppm (default: 0.28 ppm).
	\item[F11] The mixing ratio of F11 in ppm (default: 0.000268 ppm).
	\item[F12] The mixing ratio of F12 in ppm (default: 0.000503 ppm).
	\item[F22] The mixing ratio of F22 in ppm (default: 0.000105 ppm).
	\end{description}
        The \code{mixing_ratio} of F11, F12, and F22 and the default values for CH4 and N2O are ignored in case of \code{mol_abs_param reptran}.

		
\option{mol_abs_param}
	To calculate integrated shortwave or longwave irradiance, or to simulate 
	satellite instrument channels, use
	\fcode{
	   mol\_abs\_param type
	}
	to choose between the following types of schemes:
	\begin{description}
        \item[reptran]
        Representative wavelengths parameterization adapted for spectral bands.
        This parameterization is used by default if no \code{mol\_abs\_param} option is given 
        in the {\sl uvspec} input file.
        Different band widths may be selected by
        \fcode{
          mol\_abs\_param reptran [fine|medium|coarse]
        } 
        (fine: 1cm$^{-1}$; medium: 5cm$^{-1}$; coarse: 15cm$^{-1}$; coarse is default).
        The data files for coarse resolution are included in the libRadtran package.
        The files required for fine and medium resolution can be downloaded from the libRadtran homepage.
        Absorption data is mainly based on HITRAN 2004. Absorption by H2O, CO2, O3, N2O, CO, CH4, O2, N2, and NO2 
        is considered, and absorption by all other gases is zero. 
        By default volume mixing ratios of N2O, CO, CH4, and N2 (those are not in the 
        \code{atmosphere\_file}) from the US standard atmosphere are applied. 
        Use \code{mol\_file} or \code{mol\_modify} to change the gas profiles.
        In case of radiative transfer problems with solar source, the extraterrestrial spectrum from
        Kurudz is applied by default. This parameterization is described in detail by \citet{gasteiger2014}.
        \item[reptran\_channel]
        Representative wavelengths parameterization for satellite channels. Usage 
        \fcode{
          mol\_abs\_param reptran\_channel channel\_name
        }
        Channel-integrated quantities are obtained using \code{output\_process per\_band}.
        The file \file{data/correlated\_k/reptran/channel\_list.txt} provides a list of available channels; 
        more information on the channels is provided in \file{data/filter/}. 
        See \citet{gasteiger2014} for details about the approach.
        \item[crs]
        Switch off spectral parameterizations. Only molecular absorption cross sections from 
        \code{crs_file} (including the default ones) are considered.
	\item[kato] 
	\citet{Kato1999b} correlated-k distribution, shortwave; based on HITRAN 96. Please note that the 
	bands above 2.5 micrometer are not very reliable which, however, this has only little impact
	on integrated shortwave radiation.
	\item[kato2] 
	\citet{Kato1999b}, shortwave; optimized version (Seiji Kato, personal communication, 2003);
	please note that \code{kato2} only has 148 subbands (that is, calls to the \code{rte\_solver}) 
	compared to 575 for \code{kato} which translates to a increase in computational speed by 
	up to a factor of 4 with only little increase in uncertainty. The absorption data are
	based on HITRAN 2000. Please note that the bands above 2.5 micrometer are not very reliable which, 
	however, this has only little impact on integrated shortwave radiation.
	\item[kato2.96]
	\citet{Kato1999b}, shortwave; optimized version (Seiji Kato, personal communication, 2003);
	similar to \code{kato2} but based on HITRAN96. Please note that the bands above 2.5 micrometer 
	are not very reliable which, however, has only little impact on integrated shortwave radiation.
	\item[fu]
	\citet{fu92,fu93}, shortwave and longwave; fast parameterization, developed for climate models.
	\item[avhrr\_kratz] 
	\citet{Kratz1995}, AVHRR instrument channels
	\item[lowtran]
	Gas absorption parameterization from LOWTRAN; code adopted from SBDART \citep{Ricchiazzi1998b}; 
	please see the section on "Spectral resolution".
	\item[sbdart]
	Identical to LOWTRAN.
	\end{description}
	If \code{mol\_abs\_param} kato/kato2/kato2.96/fu/avhrr\_kratz is specified, the extraterrestrial 
	flux is taken from
	internally defined files specific for each parameterization, not 
	from \code{source solar file}. The output is the integrated irradiance for 
	each band. To get e.g. integrated shortwave irradiance, simply add all 
	bands of the \citet{Kato1999b} or the \citet{fu92,fu93}
	parameterization. The five AVHRR channels are weighted sums of the
	libRadtran output. Examples how to integrate the output in the
	\code{avhrr\_kratz} case are included in the {\sl uvspec} self check
	which is initiated with  
	\code{make check}.
		
\option{mol_file}
	Specify density profiles (or matrix, see below) of various trace gases to be included in the radiative 
	transfer calculation.
	\fcode{
	mol\_file gas\_species filename [unit]
	}
	At the moment following \code{gas\_species} are included: ozone (O3), nitrogen dioxide (NO2), water vapor (H2O),
	bromine oxide (BRO), chlorine dioxide (OCLO), formaldehyde (HCHO), carbon dioxide (CO2), 
	sulphur dioxide (SO2), oxygen (O2), the oxygen dimer (O4), methane (CH4), nitrous oxide (N20), 
        carbon monoxide (CO), and nitrogen (N2). 
	The gas species is identified 
	by their abbrevations given in the parenthesis above.

	The model expects a density file with two columns:
	\begin{description}
	\item[1] Altitude above sea level in km.
	\item[2] The density of trace gas [in the specified unit]
	\end{description}
	The altitude grid may be different from that in \code{atmosphere\_file}. All densities inside the range 
	of the \code{mol\_file} are replaced. For all other altitudes the values from the 
	\code{atmosphere\_file} are used. If the density is specified as -1 at a level, 
	the value from \code{atmosphere\_file} is used.
        Altitude ranges not covered by the \code{atmosphere\_file} are ignored.

	\code{unit} is an optional argument to define the unit of the density. The profiles can
	be given in particles per cm$^{3}$ (\code{cm\_3}), in particles per m$^{3}$ (\code{m\_3}), as volume mixing ratio (\code{vmr}), as mass mixing 
	ratio (\code{mmr}), or as relative humidity (\code{rh}) (only for water). The default for \code{unit} is cm$^{-3}$.
	
	To scale the profile to a total column value use \code{mol\_modify}.
	
	For airmass factor calculations it is for some species necessary to account for the 
	variation of the profile with sza. This may be accomplished by specifying a \code{mol\_file} 
	in the following format:
	\fcode{
	0.0       SZA1      SZA2 ...\\
	z(1)    dens(1,1)    ...\\
	z(2)      .           .\\
	 .        .           .
	}
	where z(i) are the altitude levels above sea level in km, SZA is the solar zenith 
	angle in degrees, and dens is the density [in the specified unit] of the trace gases as 
	function of solar zenith angle and altitude. 
	The matrix may only be specified for one species. It may however be combined with profiles
	of other species. 
	A density matrix can only be used in connection with \code{rte\_solver sdisort}!
		
\option{mol_modify}
	Set the total column of a density profile. The column is integrated between the 
	user-defined \code{altitude} and TOA (top of atmosphere). The syntax is
	\fcode{
	mol\_modify species column unit
	}
	where \code{species} is one of O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, O4, SO2, 
        CH4, N2O, CO, or N2, see also \code{mol\_file}.
	The second argument is the total column value, and the third argument is the unit, 
	in which the column is given. The unit can be DU (Dobson units), CM\_2 (molecules/cm$^2$) or MM.

	Please note that the unit MM is only valid for species H2O and specifies the precipitable water 
	in kg / m2 (which is approximately 1mm).The water vapor profile is scaled accordingly. The precipitable water 
        is integrated from the user-defined \code{altitude} to TOA (top of atmosphere).

	The default units are DU for O3, and CM\_2 for all other gases. It is possible to have
	several \code{mol\_modify} commands in the input file (maximum one per species). The following sets
	the NO$_2$ total column to 1.2 DU.
	\fcode{
	mol\_modify NO2 1.2 DU
	}
		
\option{mol_tau_file}
	Location of molecular scattering or absorption optical depth file.
	\fcode{
	mol\_tau\_file sca/abs filename
	}
	\begin{description}
	\parameter{sca}	Usually, the Rayleigh scattering
	cross section is calculated from the air pressure provided in \code{atmosphere\_file}
	(scaled with \code{pressure}). Use this parameter only if you really want to specify
	the optical depth directly (e.g. for a model intercomparison). The
	optical thickness profile may be either monochromatic or spectral. 
	\parameter{abs} Usually, molecular absorption
	is calculated from trace gas concentrations provided in \code{atmosphere\_file}
	(scaled with \code{mol\_modify O3}, etc.). Use this option only if you want to specify
	the optical depth directly (e.g. for a model intercomparison) or for a line-by-line
	calculation. If a spectral \code{mol\_tau\_file} is specified, the wavelength 
	grid defined there is used as the internal wavelength grid for the radiative transfer
	calculation, if not defined otherwise with \code{wavelength\_grid\_file}.
	\end{description}

	The file can be either of the following three formats:
	\begin{description}
	\parameter{Monochromatic} 
	Column 1 is the altitude in km 
	%AK: This is probably old stuff prior to redistribute
	%where the altitude grid must be exactly equal 
	%to the altitude grid specified in \code{atmosphere\_file}. 
	Column 2 is the absorption optical depth of each layer.
	\parameter{Spectral, ASCII}
	The first line contains the level altitudes in decreasing order; the following lines 
	contain the wavelength [nm] in the first column and then the absorption optical depths
	of each layer. 
	\parameter{Spectral,  netcdf}
	An example is available at the libRadtran homepage, 
	the file \file{UVSPEC.O2A.afglms.cdf} is a line-by-line spectrum of the oxygen A-Band 
	around 760nm, calculated for the mid-latitude summer 
	atmosphere. The advantage of
	netcdf compared to ASCII is that it is much faster to read, and that the file 
	is a self-contained, including data and a description of the variables and arrays.
	It is therefore particularly useful for line-by-line calculations where usually 
	many spectral data points are involved. 
	%netcdf is a common platform independent format; the description, a library to read and 
	%write netcdf including some tools to generate netcdf is available at 
	%http://www.unidata.ucar.edu/packages/netcdf/. A \code{mol\_tau\_file abs} must obey 
	%certain rules; 
	\end{description}
	Comments start with \code{\#}. Empty lines are ignored.
		
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

\option{no_absorption}
	Switch off absorption. Please note that this option simply sets
	the absorption optical thickness to 0.
	\fcode{
	no\_absorption [name]
	}
	If \code{name} is not set, all absorption (molecular, aerosol, cloud, ice cloud, and any profile) is switched off.

	If used together with \code{xxx\_modify set tau} this might be a bit confusing but 
	probably the most logical way. E.g. when using \code{aerosol\_default} and 
	\code{aerosol\_modify set tau 1}, the aerosol optical thickness is set to 1,
	with 0.940539 scattering and 0.059461 absorption. If \code{no\_absorption} 
	is added, the absorption optical thickness is set to 0 while the scattering
	optical thickness is preserved at 0.940539 (even though 1 was specified by 
	the user). We find this the most logical solution of the problem because 
	by switching \code{no\_absorption} off and on one tests the effect of the 
	absorber in an isolated way, rather than mixing absorption and scattering. 
	The same is true for water and ice clouds. Note, that thermal emission of 
	molecules is also switched off. 

	Possible choises for the optional argument \code{name} are:
	\begin{description}
	\item[mol] Switch off molecular absorption.
	\end{description}
		
\option{no_scattering}
	Switch scattering off.
	\fcode{
	no\_scattering [name]
	}
	If \code{name} is not set, all scattering (molecular, aerosol, cloud, ice cloud, and any profile) is switched off.

	Possible choises for the optional argument \code{name} are:
	\begin{description}
	\item[mol] Switch off molecular scattering.
	\item[aer] Switch off scattering by aerosols.
	\item[wc]  Switch off scattering by water clouds.
	\item[ic]  Switch off scattering by ice clouds.
	\item[profile] Switch off scattering by any profile defined in \code{profile typename}.
	\end{description}
		
\option{number_of_streams}
	Number of streams used to solve the radiative transfer equation.
	\fcode{
	number_of_streams value
	}
	Default is 6 for fluxes and 16 for radiances.
	(For \code{rte_solver fdisort1}, \code{fdisort2} and \code{disort} only
	even number of streams are possible.)
		
\option{output_file}
	\fcode {             
	    output\_file filename
	}
	uvspec output will be written to \code{filename}. The format can be set by the option 
	\code{output\_format} (default is ascii).  
		
\option{output_format}
	\fcode {
	    output\_format format
	}
	where \code{format} is either \code{ascii} (default) or \code{flexstor}. 
	Note that \code{flexstor} does not currently work when \code{umu} and/or \code{phi} is specified.
	There is also the
	possibility to write uvspec simulation results to an existing netCDF file. For
	that \code{format} must be \code{netCDF} and the option \code{output\_file} must be given
	and point to a file that contains a lat/lon/time grid.
	If \code{format} is set to \code{sat picture} then \code{output\_file} must be given and point
	to a NetCDF-File that contains a pixel x/pixel y/time grid.
		
\option{output_process}
	Decide how the output from \code{uvspec} is processed:
	\fcode{
	output\_process type
	}
	where type is one of
	\begin{description}
	\parameter{sum} 
	Sum output over wavelength. Useful in combination with the 
	\code{mol\_abs\_param} option (\code{kato}, \code{kato2}, \code{kato2.96}, 
	\code{Fu}, \code{avhrr\_kratz}, \code{reptran}).
        In case of \code{mol\_abs\_param reptran}, the units are automatically converted to \code{per\_band} before summation. 
	\parameter{integrate} 
	Integrate output over wavelength for solar and over wavenumber for thermal simulations. 
	Useful for spectral calculations and \code{mol\_abs\_param lowtran}.
        In case of \code{mol\_abs\_param reptran} the spectral range of the integration may slightly deviate 
        from the spectral range specified by the \code{wavelength} option.
	\parameter{per\_nm } 
	Output is given in W/(m$^2$ nm) or mW/(m$^2$ nm) (W or mW is determined by the extraterrestrial spectrum.)
	\parameter{per\_cm-1 } 
	Output is given in W/(m$^2$ cm$^{-1})$ or mW/(m$^2$ cm$^{-1}$).
	\parameter{per\_band } 
	Output is given in W/m$^2$ or mW/m$^2$ per correlated-k or reptran band. (This option can not be used for 
	spectral calculations and \code{mol\_abs\_param} LOWTRAN in the solar range.)
	\parameter{none} 
	No processing - output spectral information (default).
	\end{description}
		
\option{output_quantity}
	Convert radiances / irradiances to equivalent output quantity.
	\fcode{
	output\_quantity quantity
	}
	where \code{quantity} can be one of the following:
	\begin{description}
	\parameter{brightness}
	Convert radiances / irradiances to equivalent brightness temperatures.
	\parameter{reflectivity}
	Calculate transmission / reflectivity instead of absolute quantities. 
	For irradiances / actinic fluxes the transmission T is defined as 
	\begin{equation}
	   T = {{E}\over{E_0 \cos{\theta}}}
	\end{equation}
	  where $E$ is the irradiance / actinic flux, $E_0$ is the extraterrestrial flux,
	  and $\theta$ is the solar zenith angle.
	The reflectivity R is defined as 
	\begin{equation}
	   R = {{\pi \cdot L}\over{E_0 \cos{\theta}}}
	\end{equation}
	 where $L$ is the radiance, $E_0$ is the extraterrestrial flux,
	 and $\theta$ is the solar zenith angle.
	Obviously, reflectivities do not depend on Sun-Earth distance. Please 
	note the difference to \code{transmittance}.
	\parameter{transmittance}
	Calculate transmittance / reflectance instead of absolute quantities. 
	That is, set the extraterrestrial irradiance to 1 and do not correct for Sun-Earth distance:
	\begin{equation}
	   T = \frac{E}{E_0}
	\end{equation}
	  where $E$ is the irradiance / actinic flux / radiance and $E_0$ is the extraterrestrial flux.
	Please note the difference to \code{reflectivity}.
        \end{description}
		
\option{output_user}
	User defined output. 
	%This option is case sensitive.
	Here the user may specify the columns desired for output. 
	\fcode{
	output\_user format
	}
	where \code{format} is one or more of the following.
	\begin{description}
	\parameter{lambda} 
	Wavelength in nm.
	\parameter{wavenumber} 
	Wave number in cm$^{-1}$.
	\parameter{sza} 
	solar zenith angle
	\parameter{zout} 
	Output altitude in km.
	\parameter{edir, eglo, edn, eup, enet, esum}
	The direct, global, diffuse downward, and diffuse upward irradiance.
	Net is global - upward, sum is global + upward. 
	\parameter{uu} 
	Radiances:  uu(umu(0),phi(0)) ... uu(umu(0),phi(m)) ... uu(umu(n),phi(0)) ... uu(umu(n),phi(m))
	\parameter{fdir, fglo, fdn, fup, f}
	The direct, global, diffuse downward, diffuse upward, and total actinic flux.
	\parameter{uavgdir, uavgglo, uavgdn, uavgup, uavg} 
	The Direct, global, diffuse downward, diffuse upward, and total diffuse 
	mean intensity (= actinic flux / 4$\pi$).
	\parameter{spher\_alb}
	Spherical albedo of the complete atmosphere.
	\parameter{albedo} 
	Albedo.
	\parameter{heat} 
	Heating rate in K/day.
	\end{description}

	It is also possible to gain some information about the atmosphere and the clouds:
	\begin{description}
	\parameter{p} 
	pressure [hPa], ,
	\parameter{T, T\_d}
	temperature [K], dewpoint temperature [K]
	\parameter{T\_sur} 
	surface temperature [K]
	\parameter{theta } 
	potential temperature [K]
	\parameter{theta\_e} 
	equivalent potential temperature [K]
	\parameter{n\_xxx} 
	number density of the gas xxx [cm$^{-3}$]
	\parameter{rho\_xxx} 
	mass density of the gas xxx [kg/m$^3$]
	\parameter{mmr\_xxx} 
	mass mixing ratio of the gas xxx [kg/kg]
	\parameter{vmr\_xxx} 
	volume mixing ratio of the gas xxx [m$^3$/m$^3$]
	\parameter{rh} 
	relative humidity over water [percent]
	\parameter{rh\_ice} 
	relative humidity over ice   [percent]
	\parameter{c\_p} 
	specific heat capacity of the air (humidity and temperature dependent)
	\parameter{CLWC} 
	cloud liquid water content   [kg/kg]
	\parameter{CLWD} 
	cloud liquid water density   [g/m$^3$]
	\parameter{CIWC} 
	cloud ice water content      [kg/kg]
	\parameter{CIWD} 
	cloud ice water density      [g/m$^3$]
	\parameter{TCC} 
	total cloud cover            [0-1]
	\end{description}
	where \code{xxx} is one of 
	AIR, O3, O2, H2O, CO2, NO2, BRO, OCLO, HCHO, or O4.

	Default output is  
	\fcode{
	output\_user lambda lambda, edir, edn, eup, uavgdir, uavgdn, uavgup 
	}
	for \code{fdisort1}, \code{sdisort}, and \code{spsdisort}, whereas the default for 
	\code{twostr} is 
	\fcode{
	output\_user lambda, edir, edn, eup, uavg.
	}
	The lines containing radiances and the output of \code{rte\_solver polradtran} are not affected. 
		
\option{ozone_column}
Obsolete, use \code{mol_modify} instead.
%Set the ozone column to a given value. The column is integrated between the 
%user-defined \code{altitude} and TOA (top of atmosphere). The value must be in 
%Dobson units. If value is negative or \code{ozone_column}
%is not specified, the ozone column is not scaled.

%CE Commented for ESA version
\option{phi}
	Azimuth output angles (in degrees) in increasing order. 
	\fcode{
	phi values
	}
	The radiance is output at \code{phi} and \code{umu}.
	\begin{itemize}
	\item Sensor in the North (looking South):  0 deg 
	\item Sensor in the East  (looking West):   90 deg 
	\item Sensor in the South (looking North):  180 deg 
	\item Sensor in the West  (looking East):   270 deg 
	\end{itemize}
	For all one-dimensional solvers the absolute azimuth does not matter, but only the relative azimuth \code{phi}-\code{phi0}.
		
\option{phi0}
	Azimuth angle of the sun (0 to 360 degrees). 
	\fcode{phi0 value}
	\begin{itemize}
	\item Sun in the South:    0 degrees
	\item Sun in the West:    90 degrees
	\item Sun in the North: 180 degrees
	\item Sun in the East: 270 degrees
	\end{itemize}
	For all one-dimensional solvers the absolute azimuth does not matter, but only the relative azimuth \code{phi}-\code{phi0}.
		
\option{polradtran}
	Specify polradtran values. This option is only relevant for \code{rte_solver polradtran}.

	\fcode{
	polradtran aziorder value
	}
	Order of Fourier azimuth series.
	The value \code{0} (default for irradiance) is the azimuthally symmetric case.
	For radiance computation a higher order is required, thus the default for radiances is 4. 

	\fcode{
	polradtran nstokes value
	}
	Number of Stokes parameters
	where \code{value} is one of
	\begin{description}
	\item[1] for I (no polarization, default)
	\item[2] for I,Q,U (Since V is very small in the atmosphere, it makes sense
	           to compute only I,Q,U. This saves computation time and
		   memory).
	\item[3] for I,Q,U,V
	\end{description}
	Default is 1.

	\fcode{
	polrdatran src_code value
	}
	Radiation sources included
	which may be
	\begin{description}
	\item[0] none
	\item[1] solar
	\item[2] thermal
	\item[3] both
	\end{description}
	Default 1.
		
\option{polradtran_max_delta_tau}
	Initial layer thickness for doubling; governs accuracy, 10E-5 should be
	adequate. Do not go beyond half the real precision, i.e. 10e-8 for REAL*8.
	Default 1.e-05.
	\fcode{
	polradtran_max_delta_tau value
	}
	This option is only relevant for \code{rte_solver polradtran}.
		
\option{polradtran_quad_type}
	Type of quadrature used:
	\fcode{
	polradtran_quad_type type
	}
	where \code{type} is one of
	\begin{description}
	\item[G] (default) gaussian
	\item[D] double gaussian, 
	\item[L] Lobatto
	\item[E] extra-angle(s), this must be used of \code{polradtran} is used in 
	combination with \code{umu}. Will internally use Gaussian scheme (\code{G}). 
	See also radtran documentation (\file{libsrc_f/README.polRadtran}).
	\end{description}
	Default G.                                  
	This option is only relevant for \code{rte_solver polradtran}.
		
\option{pressure}
	The surface pressure (at the user-defined \code{altitude}) in hPa. 
	\fcode{
	pressure value
	}
	The pressure profile as well as air, O2 and CO2 density profiles 
	are scaled accordingly.
		
\option{pressure_out} 
	Specify the output levels in pressure coordinates. The syntax is   
	\fcode{ 
	pressure\_out p1 p2 ... 
	}
	where '\code{p1 p2 ...}' are the output levels in hPa.  
	The pressure output levels must be sorted in decreasing order.  
	Output pressure levels must be within the range defined in the 
	\code{atmosphere\_file}. You can also use \code{toa} for  
	\code{top of atmosphere} and \code{sur} for \code{surface altitude} and \code{cpt}  
	for \code{cold point tropopause}. 
		
\option{print_disort_info}
	Specify one or more integers between 1 and 7.
	\fcode{
	print\_disort\_info  value
	}
	Print various disort input and output in disorts own format. See 
	\file{libsrc\_f/DISORT2.doc} for more information.
	\strong{Warning:} Produces a lot of output. 
		
\option{profile_file}
	Define file containing properties of clouds, aerosols, hydrometeors, etc. 
	This option is a generalization of the options \code{wc\_file} and \code{ic\_file}.
	
	Usage:
	\fcode{
	profile\_file typename type file
	}
	
	\code{typename} describes the name of the profile; typically this 
	describes what kind of particles are dealt with here. Examples
	are \code{wc} (water clouds), \code{ic} (ice clouds), \code{aer} (aerosols), drizzle. 
	The \code{typename} is needed to refer to this profile when using other options, 
	such as \code{profile\_properties}. Note that \code{typename} "wc" and "ic"
	have special effects (i.e.~default properties, and "ic" properties are
	not allowed with "wc" files, and vice versa).
	
	\code{type} defines the file type, which is identical to \code{wc\_file type}.
	Please refer to \code{wc\_file} for choises of \code{type} and a detailed description 
	on the file structures.	
		
\option{profile_modify}
	Modify profile optical properties.
	\fcode{
	profile typename variable scale/set value
	}
	\code{typename} describes the name of the profile; it must be identical
	to the one defined in \code{profile\_file}.

	This option is identical to \code{wc\_modify}.
	Please refer to \code{wc\_modify} for a detailed description of \code{variable}.
		
\option{profile_properties}
	Define how liquid/ice water content/mass concentration and effective
	particle radius are translated to optical properties for profile
	\code{typename}. This option is a generalization of the options
	\code{wc\_properties} and \code{ic\_properties}.
	
	Usage:
	\fcode{
	profile\_properties typename property [interpolate]
	}
	
	\code{typename} describes the name of the profile; it must be identical
	to the one defined in \code{profile\_file}.
	
	Please refer to \code{wc\_properties} and \code{ic\_properties} for possible choices for \code{property}.
		
\option{pseudospherical}
	Invokes pseudo-spherical geometry in disort/twostr. Default is plane-parallel.
		
\option{quiet}
	If specified, informative messages are turned off. See also \code{verbose}.
		
\option{radiosonde}
	This option allows to change the temperature and pressure profile, and optionally to 
	specify one or more density profiles. The entry in the input file looks like this:
	\fcode{
	radiosonde filename [gas\_species] [unit] ...
	}
	Currently the following gas\_species are included: ozone (O3), nitrogen dioxide (NO2), 
	water vapor (H2O), bromine oxide (BRO), chlorine dioxide (OCLO), formaldehyde (HCHO), 
	carbon dioxide (CO2), sulphur dioxide (SO2), and the oxygen dimer (O4). 
	Each gas species is identified by its abbrevations given in parentheses above.
	Unit is an optional argument to defines the unit of the density. The profiles can
	be given in particles per cm$^3$ (CM-3), in particles per m$^3$ (M-3), as volume 
	mixing ratio (VMR), as mass mixing ratio in kg/kg (MMR), or as relative humidity 
	(RH) (only for water). The default unit is RH for water vapour, 
	MMR for ozone, and CM3 for all other gases.
	The radiosonde file must have (2 + number of gases) columns:
	\begin{description}
	\item[1] pressure in hPa
	\item[2] temperature in Kelvin
	\item[3, 4, ...] density of trace gas in the specified unit
	\end{description}
	A new z-grid will be calculated, starting at \code{altitude} and assuming a linear temperature variation
	between levels. The air density will be recalculated according to the ideal gas law, and the density of 
	the well mixed gases O2 and CO2 will be scaled accordingly.
	The atmospheric data above the radiosonde data is taken from the \code{atmosphere\_file} level by level, starting 
	at the first pressure level above the radiosonde data. The z-grid of the \code{atmosphere\_file} in 
	this height region is shifted accordingly.
	Also if the density in the radiosonde file is specified as -1 at a level, 
	the value from the \code{atmosphere\_file} is used.
	Possible calls are
	\fcode{
	radiosonde ../examples/radiosonde.dat
	}
	just in order to change the temperature and pressure profile, or
	\fcode{
	radiosonde ../examples/radiosonde2.dat H2O RH O3 MMR NO2
	}
	where water vapour density will be given as relative humidity, ozone as mass mixing ratio, 
	and NO2 in cm$^{-3}$ (default).
		
\option{radiosonde_levels_only}
	The atmosphere considered in the simulation has the same height range as the data in 
	the \code{radiosonde}-file. No further levels are added above those.
	This option has only an effect in combination with \code{radiosonde}.
		
\option{raman}
	\fcode{
	raman [original]
	}
	The \code{raman} option includes single order rotational Raman scattering in 
	the calculation. The solution treats Raman as a perturbation similar to the 
	approaches of \citet{Vountas1998} and \citet{Spurr2008}.
	
	The \code{raman} option may only be used for spectral calculation.

	The disort radiative transfer solver with a general source is needed to solve the radiative transfer
	equation including Raman scattering. This solver is automagically invoked when specifying 
	the \code{raman} option. It is thus not neccessary to set the \code{rte_solver}.

	The \code{raman} option is optimized with respect to speed. The optimized implementation
	should be just as accurate as the original version. To use the original version
	invode \code{raman original}. With the optional argument \code{original} each wavelength 
	is treated individually and is thus accurate, but computationally very expensive. 

	Please note that while the \code{raman} option has been extensively tested and verified, it 
	is nevertheless a new option, hence, use it with care.
		
\option{rayleigh_depol}
	Rayleigh depolarization factor.
	\fcode{
	rayleigh\_depol value
	}
	The Rayleigh scattering phase function is
	$p(\mu) = a + b  \mu^2$ where $a = 1.5{(1+\texttt{depol})/(2+\texttt{depol})}$ and  
	$b = 1.5{(1-\texttt{depol})/(2+\texttt{depol})}$. By default the depolarization is calculated
	using the expressions from \citet{Bodhaine1999}.
		
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

\option{reptran_file}
        Location of the representative wavelengths file. 
        \fcode{
          reptran\_file file
        }
        This option is useful together with 'mol\_abs\_param reptran' and 'mol\_abs\_param reptran\_channel' 
        options, if you want to use your own representative wavelengths parameterization.
        
\option{reverse_atmosphere}
	Option for the strong and bold. Reverses the atmospheric input to the radiative
	transfer solvers. That is, the atmosphere is turned on the head.
	Yes, that is actually useful for some purposes. If you think
	you need this contact the author. Otherwise, do not use.
		
\option{rte_solver}
	Set the radiative transfer equation solver to be used. 
	\fcode{
	rte_solver type
	}
	If not specified the default \code{rte_solver} is \code{disort}.
	Choices for \code{type} are
	\begin{description}

	\parameter{disort}
	C-version of the disort algorithm, translated from Fortran by Tim Dowling. 
	This is the recommended discrete ordinate code in {\sl libRadtran}.
	For documentation see \file{src\_f/DISORT2.doc} as well as the papers and
	the DISORT report at
	ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/.  The intensity
	correction can be performed according to \citet{nakajima1988} using
	\code{disort_intcor moments} (like in the original code), or with the 
	improvements described in (Buras, Dowling, Emde, in preparation; default). Can be run in  
	plane-parallel geometry (default) or in pseudo-spherical geometry 
	(using \code{pseudospherical}).

	\parameter{twostr}
	C-version of the two--stream radiative transfer solver described by
	\citet{Kylling1995}. Can be run in plane-parallel geometry (default)
	or in pseudo-spherical geometry (using \code{pseudospherical}
	).

	\parameter{fdisort1}
	The standard plane--parallel disort algorithm by \citet{Stamnes1988c},
	version 1.3 -- provided for compatibility reasons. Use only if you have troubles with the default			         
        \code{disort} or for historical reasons. For documentation see \file{src\_f/DISORT.doc} as well 
	as the papers and the DISORT report at
	ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/. To optimize
	for computational time and memory, please adjust the parameters in
	\file{src\_f/DISORT.MXD} for your application and re-compile. For your
	application please use \code{rte_solver fdisort2} which is the advanced
	version, unless you e.g. want to explore how a specific feature of
	\code{fdisort2} (e.g. the \citet{nakajima1988} intensity correction)
	improves the \code{fdisort1} result.

	\parameter{fdisort2}
	Version 2 of the Fortran algorithm disort -- provided for compatibility reasons.
	Use only if you have troubles with the default \code{disort} or for historical reasons. For documentation see
	\file{src\_f/DISORT2.doc} as well as the papers and the DISORT report
	at ftp://climate1.gsfc.nasa.gov/wiscombe/Multiple\_Scatt/
	\code{fdisort2} has several improvements compared to its 'ancestor'
	\code{fdisort1} (version 1.3). To optimize for computational time and
	memory, please adjust the parameters in \code{src_f/DISORT.MXD} for
	your application and re-compile.
	Note! \code{fdisort2} is a new version of the original disort code
	which was implemented in summer 2009. It uses phase functions to
	calculate the intensity corrections by \citet{nakajima1988} instead of
	Legendre moments. Hence it needs cloud properties files which contain
	the phase functions. It is still possible to use the old version of
	disort, you need to specify \code{disort_intcor moments}.

	\parameter{sdisort}
	Pseudospherical disort as described by \citet{Dahlback1991}. Double
	precision version. To optimize for computational time and memory,
	please adjust the parameters in \file{src\_f/DISORT.MXD} for your
	application and re-compile.  

	\parameter{spsdisort} Pseudospherical
	disort as described by \citet{Dahlback1991}, single precision
	version. \strong{Warning:} it is not recommended to use
	\code{spsdisort} for really large solar zenith angles nor for cloudy
	conditions. For large optical thickness it is numerically unstable and
	may produce wrong results. To optimize for computational time and
	memory, please adjust the parameters in \file{src_f/DISORT.MXD} for
	your application and re-compile.  

	\parameter{polradtran} The plane-parallel radiative transfer solver of \citet{Evans1991}.
	Includes polarization. The full implementation of the polRadtran
	solver in \code{uvspec} is quite new (version 1.4). If you find
	unusual behaviour, please contact the {\sl libRadtran} authors.

	\parameter{ftwostr}
	Original Fortran-version of the two--stream radiative transfer solver described by
	\citet{Kylling1995}, in pseudo-spherical geometry.

	\parameter{rodents}
	Delta-Eddington two--stream code (RObert's Delta-EddingtoN Two-Stream), plane-parallel.

	\undocumented{
	 \parameter{twostrebe}
	 Delta-Eddington two--stream solver by Bernhard Mayer. \code{zout_interpolate} is activated automatically.
	}

	\parameter{sslidar}
	A simple single scattering lidar simulator by Robert Buras.

	\parameter{sos}
	A scalar pseudospherical succesive orders of scattering code. Works
	for solar zenith angles smaller than 90 degrees. Can calculate
	azimuthally averaged radiances. Set \code{sos_nscat} to specify the order
	of scattering.

	\ifmystic{
	\parameter{montecarlo}
	The MYSTIC Monte Carlo code. Monte Carlo is
	the method of choice (1) for horizontally
	inhomogeneous problems; (2) whenever
	polarization is involved; (3) for applications
	where spherical geometry plays a role; and (4)
	whenever sharp features of the scattering
	phase function play a role, like for the 
	calculation of the backscatter glory or the
	aureole.  
	\parameter{mystic}
	Same as \code{montecarlo}.
	}

	\parameter{tzs} TZS stands for "thermal, zero scattering" and is a
	very fast analytical solution for the special case of thermal emission
	in a non-scattering atmosphere. Please note that TZS does only
	radiance calculations at top of the atmosphere. "Blackbody clouds" may be included 
        using the option \code{tzs_cloud_top_height}. 

	\parameter{sss}
	SSS stands for "solar, single scattering" and is an analytical single
	scattering approximation which might be reasonable for an optically
	thin atmosphere. Please note that SSS does only radiance calculations
	at top of the atmosphere. This is an experimental solver - be careful!

	\undocumented{
	\parameter{sssi} Missing documentation!!!
	}

	\parameter{null}
	The NULL solver does not solve the radiative transfer
	equation. However, it sets up the optical properties, and does the
	post-processing; useful if you are either interested in the overhead
	time required by a particular model input or if you are simply
	interested in the optical properties, as output by \code{verbose}.
	\end{description}
	Default: \code{disort}
		
\undocumented{
\option{satellite_geometry}
	With this option the satellite geometry is determinded. The argument for this option 
	\fcode{ 
	satellite_geometry  netCDF_file  
	} 
	is the location of a \code{netCDF_file}, which must contain latitude and longitude position 
	as well as zenith and azimuth viewing angle for each pixel. 
 
	%ak20100306 What is below does not belong in the general Users Guide yet. 
	%For format specification see the example file data/satellite/MSG\_seviri/MSG\_seviri\_geometry.nc. 
	%This option has to be used in combination with \code{satellite_pixel}. 
	%You can also use the abbreviation 
	%\code{satellite_geometry MSG} 
	%for MSG simulation, \emph{but in order to use this, you have to copy the netCDF file} 
	%\fcode{ 
	%cp -r /data/A3/satellite_geometry/satellites libRadtran/data/ 
	%} 
	%\emph{which is not in cvs version, as it is too large.} 
		}
\undocumented{
\option{satellite_pixel}
	This option specifies which pixel of the satellite image that should be simulated.  
	\fcode{ 
	satellite_pixel  pixel_x pixel_y 
	} 
	The arguments \code{pixel_nr_x} and \code{pixel_nr_y} specifies the pixel position in the native  
	system of the satellite, which is determinded by the option \code{satellite_geometry}. 
		}
\option{sdisort}
	Specify sdisort values. This option is only relevant for \code{rte_solver sdisort}.

	%CE: Commented for ESA version
	\undocumented{
	\fcode{
	sdisort ichapman value
	}
	Undocumented option for sdisort. Default is 1.
	}

	\fcode{
	sdisort nscat value
	}
	Set the order of scattering. 
	If \code{value} is set to 1 \code{sdisort} will run in single scattering mode 
	while if set to 2, \code{sdisort} runs in full multiple scattering mode.
	Default is 2 for \code{rte_solver sdisort}. 

	\fcode{
	sdisort nrefrac value
	}
	Include refraction where \code{value} has the meaning
	\begin{description}
	\item[0] No refraction, default.
	\item[1] Refraction included using fast, but harsh method.
	\item[2] Refraction included using slow, but accurate method.
	\end{description}
	If refraction is included also set parameter \code{refraction_file}.
		
\option{slit_function_file}
	If specified, the calculated spectrum is convolved with the function found
	in the \file{slit\_function\_file}. 
	\fcode{
	slit\_function\_file file 
	}
	The file must contain two columns. Column 1 is the wavelength, in nm, and relative to the
	center wavelength. Column 2 is the corresponding slit function value. It must
	be unity at the maximum. The wavelength steps in the slit function file must
	be equidistant. Comments start with \code{\#}. Empty lines are ignored. Please 
	note that prior to convolution the spectrum is interpolated to the wavelength
	steps of the slit function. For this reason, make sure that the resolution
	of the slit function is high enough even if the slit function is e.g. a 
	simple triangle which could in principle be described with 3 grid points. For an 
	example see \file{examples/TRI\_SLIT.DAT} and the \code{make\_slitfunction} tool.
		
\undocumented{
\option{sos_nscat}
	Set the order of scattering for \code{rte_solver sos}.
	\fcode{
	sos_nscat value
	}
	Default is 20.
		}
\option{source}
	Set the radiation source type
	\fcode{
	source  type 
	}
	where \code{type} is either \code{solar} or \code{thermal}.
	Solar radiation is per default output in W/(m2 nm) if no \code{mol\_abs\_param} is specified 
        or \code{mol\_abs\_param} options \code{crs}, \code{lowtran}, or \code{reptran} are specified.
        For all other \code{mol\_abs\_param} options
	the output is integrated over the wavelength band.
	Thermal radiation is per default output in W/(m2 cm-1), if REPTRAN is used or the bandwidth 
	is equal to 1 cm-1 (default for \code{mol\_abs\_param lowtran}). 
	Otherwise the output is the integrated flux over the
	wavenumber interval specified by \code{thermal\_bandwidth}, \code{thermal\_bands\_file},
	or by the \code{mol\_abs\_param} option (\code{kato}, \code{kato2}, \code{kato2.96}, 
	\code{fu}, or \code{avhrr\_kratz}).
	\fcode{
	source type [file] [unit]
	}
	The second argument \code{file} specifies the location of file holding the extraterrestrial spectrum.
        In general, \code{file} is required for solar calculations if \code{mol\_abs\_param} is not used.
	\code{file} is ignored if \code{mol\_abs\_param} other than \code{lowtran} oder \code{reptran} is specified.

	The file must contain two columns. Column 1 is the wavelength in nm, and column 2 
	the corresponding extraterrestrial flux. The user may freely use any units
	he/she wants for the extraterrestrial flux. The wavelength  specified grid
	defines the wavelength resolution at which results are returned. However, 
	the wavelength range is determined by \code{wavelength}. \code{file} may be 
	omitted for thermal radiation calculations (\code{source thermal}) as well as 
	\code{output_quantity transmittance} and \code{output_quantity reflectivity} calculations. If omitted, the 
	output resolution equals the internal wavelength grid which the model chooses 
	for the radiative transfer calculation.
	Comments start with \code{\#}. Empty lines are ignored.

	For some purposes it is useful to tell libRadtran the units of the spectrum.
	This can be done with the optional third argument. 
	Possible choises for \code{unit} are \code{per\_nm}, \code{per\_cm-1} or \code{per\_band}.
	If \code{unit} is set to \code{per\_nm}
	libRadtran assumes that the unit of the spectrum is W/(m2 nm), if set to \code{per\_cm-1}
	it assumes W/(m2 cm-1).

		
\option{spline}
	\fcode{
	spline lambda\_0 lambda\_1 lambda\_step
	}
	Spline interpolate the calculated spectrum between wavelengths \code{lambda\_0} 
	and \code{lambda\_1} in steps of \code{lambda\_step},
	in nm. Specified as e.g.
	\fcode{
	spline 290. 365. 0.5
	}
	Here, the calculated spectrum is interpolated to wavelengths 290.0, 290.5, 291.0, 
	..., 364.5, 365.0. For interpolation to arbitrary wavelengths use \code{spline\_file}.
	The specified wavelength interval must be within the one specified by \code{wavelength}.
		
\option{spline_file}
	Spline interpolate to arbitrary wavelengths, in nm, given as a single column in file 
	\file{spline\_file}. 
	\fcode{
	spline\_file file
	}
	The specified wavelengths must be within the range specified 
	by \code{wavelength}. Comments start with \code{\#}. Empty lines are ignored.
		
\option{sslidar}
	Set single scattering lidar parameters (\code{rte_solver sslidar}).
	\fcode{
	sslidar variable value
	}
	\code{variable} can be one of the following:
	\begin{description}
	\item[area] Set area of single scattering lidar in units of square meters (default: 1.0)
	\item[E0] Set Laser pulse energy for single scattering lidar in units of (default: 0.1)
	Joule (You can also use a \code{source solar file}
	instead... not yet implemented.)
	\item[eff] Set lidar efficiency for single scattering lidar (default: 0.5)
	\item[position] Set lidar position for single scattering lidar in units of km (default: 0.0)
	\item[range] Set lidar range bin width for single scattering lidar in units of km (default: 0.1)
	\end{description}
		
\option{sslidar_nranges}
	Set number of range bins for single scattering lidar (solver
	\code{sslidar}). Default is 100.
		
\option{sslidar_polarisation}
	Turn on polarisation measurement for lidar (solver \code{sslidar}). Default is without polarisation.
		
\option{sur_temperature}
	Surface temperature, used for thermal infrared calculations. 
	\fcode{
	sur\_temperature value
	}
	If not specified, the
	temperature of the lowest atmospheric level is used as surface temperature.
		
\undocumented{
\option{surface_temperature_map}
	Specify a surface\_temperature map with a \emph{netCDF} file which is used in combinition 
	with the options \code{latitude}, \code{longitude}, and \code{time}.
	\fcode{
	surface_temperature_map file [variable_name]
	}
	where \code{file} is the location of the \emph{netCDF} file. 
	\code{libRadtran} reads the value at the nearest pixel to the given \code{latitude} and 
	\code{longitude}. No spatial interpolation or averaging of the values is done.
		}
\undocumented{
\option{surface_type_map}
	Specify a surface type map, which is used in combinition with 
	\code{albedo_library}, \code{latitude}, and \code{longitude} in order to 
	select the surface type relevant for the simulation.
	No pixel interpolation is done.
	The format of the call is:
	\fcode{
	surface_type_map file [variable_name]
	}
	where \code{file} is the location of the surface type map file. The map is expected to 
	be in \emph{netCDF} format. The file must contain the variables \code{double lat(nlat)}, 
	\code{double lon(nlon)}, and \code{byte brdf_rpv_type (nlat, nlon)}. If the name of the 
	surface type variable is different, the optional argument can be used in order to specify 
	the variable name. For format specification see also 
	\file{data/albedo/IGBP_map/SURFACE_TYPE_IGBP_10min.cdf}.

	For using the IGBP map, the call is \code{surface_type_map IGBP}. This map has a resolution of 10 
	minutes and contains the surface types 1 to 18 defined in the \code{albedo_library IGBP}. 
	Fresh snow and sea ice are not included, as their extent is too variable. Attention:  
	That implies e.g. that the Arctic is considered ocean\_water and not sea\_ice!

	Locations on the pixel boundaries are interpreted as the pixel northward and
	eastward respectively. E.g. location 0 N, 0 E is interpreted like the pixel ranging from
	0 to 10min North and from 0 to 10min East.
		}
\option{sza}
	The solar zenith angle (degrees). 
	\fcode{
	sza value
	}
	The default solar zenith angle is 0.
		
\option{sza_file}
	Location of solar zenith angle file for wavelength-dependent solar
	zenith angle. 
	\fcode{
	sza_file file
	}
	This option is useful if you want to simulate an instrument 
	which scans so slowly that the solar zenith angle may change significantly during 
	the wavelength scan. 
	The file must have two or three columns. Column 1 is the wavelength, 
	in nm, and column 2 the corresponding  solar zenith angle. Optionally
	the third column may contain the corresponding solar azimuth angle. The solar 
	azimuth angle is only needed when calculating radiances. The wavelength 
	grid may be freely set. The solar zenith and azimuth angle will be interpolated to 
	the wavelength grid used for the radiation calculation.
	Comments start with \code{\#}. Empty lines are ignored.
		
\undocumented{
\option{test_optical_properties}
	Write optical properties to file. This is made for test suite. 
	rte_solver will be switched to null_solver.
		}
\option{thermal_bands_file}
	File with the center wavelengths and the wavelength band intervals to be used for 
	calculations in the thermal range. 
	\fcode{
	thermal\_bands\_file file
	}
	The following three columns are expected:
	center (or reference) wavelength, lower wavelength limit, upper wavelength limit [nm].
	\code{thermal\_bands\_file} defines the wavelength grid for the radiative transfer 
	calculation. The RTE solver is called for each of the wavelengths in the first column. 
	The atmospheric (scattering, absorption, etc) properties are also evaluated at these 
	wavelengths. For thermal radiation calculations, the Planck function is integrated 
	over the wavelength bands defined in the second and third columns. The result will
	therefore be a band-integrated irradiance which does only make sense when the 
	\code{source solar file} grid equals the \code{thermal\_bands\_file} grid.
				
\option{thermal_bandwidth}
	Specify a constant bandwidth in cm-1 for thermal calculations. 
	\fcode{
	thermal\_bandwidth value
	}
	The default is 1 cm-1.
	This option is ignored if used together with \code{mol\_abs\_param kato/kato2/kato2.96/fu/avhrr\_kratz}.
			
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

\undocumented{
\option{tipa}
	With the option \code{tipa dirdiff}, the tilted independent
	pixel (column) approximation (TIPA, the same as TICA) is carried out,
	where a 3D-cloud file serves as input and the columns are tilted
	towards the sun (see, e.g., Wapler and Mayer, 2008).  \code{tipa
	dirdiff} works with all 1D-rte-solvers.  When choosing
	tipa (\code{tipa dir}), from the "tilted clouds" only the optical
	depth for the calculation of the direct radiation is derived and used
	later to calculate the diffuse radiation. To calculate the diffuse
	radiation, the original cloud-matrix is still necessary /
	maintained. \code{tipa dir} works only with Roberts (Buras)
	delta--Eddington two--stream model \code{rodents}.  For both tipa
	options, a \code{mc_sample_grid} has to be specified as Nx Ny dx dy. A
	3d-cloud-file (with flag 3) serves as input (see
	\code{wc_file 3D}). The output is written into mc.flx(.spc),
	however, a basename via \code{mc_basename} can be given.
		}
\option{tzs_cloud_top_height}
	Define cloud top height of "blackbody clouds" for radiative transfer solver \code{tzs}.
		
\option{umu}
	Cosine of output polar angles in increasing order, starting with 
	negative (downwelling radiance, looking upward) values (if any) and on through 
	positive (upwelling radiance, looking downward) values. 
	Must not be zero.
	\fcode{
	umu values
	}
		
\option{verbose}
	If specified abundances of informative messages are output to stderr. To make
	use of this information, you may want to write the standard \code{uvspec} output to 
	one file and the diagnostic messages to another. To do so, try
	\code{(./uvspec < uvspec.inp > uvspec.out) >\& verbose.txt} (depending on your shell you 
	might need a slightly different syntax). The irradiances 
	and radiances will be written to \file{uvspec.out} while all 
	diagnostic messages go into \file{verbose.txt}. See also \code{quiet}.
		
\option{wavelength}
	Set the wavelength range by specifying first and last wavelength in nm. 
	\fcode{
	wavelength lambda\_0 lambda\_1
	}
	The default output wavelength grid is that defined in \code{source solar file}, 
	unless \code{spline} is specified. Note that the radiative transfer calculations 
	are done on an internal grid which can be influenced with \code{wavelength\_grid\_file}
	or \code{mol\_tau\_file abs file}
				
\option{wavelength_grid_file}
	Location of single column file that sets the wavelength grid used for the 
	internal transmittance calculations. 
	\fcode{
	wavelength\_grid\_file file
	}
	The wavelengths must be in nm.
	Do not use this option unless you know what you are doing.
	Comments start with \code{\#}. Empty lines are ignored.
				
\option{wavelength_index}
	Set the wavelengths to be selected. To be used together with predefined wavelength grids,
	such as \code{wavelength\_grid\_file}, \code{mol\_tau\_file abs file} and particularly 
	useful in combination with the \code{mol\_abs\_param} option where often only a 
	specified number of wavelength bands is required. E.g., in combination with
	\code{mol\_abs\_param AVHRR\_KRATZ}, \code{wavelength\_index 15 15} will select wavelength
	index 15 which corresponds to channel 4, or \code{wavelength\_index 10 14} will select 
	those bands required for channel 3. Indices start from 1.
				
\option{wc_file}
	Location of file defining water cloud properties. 
	\fcode{
	wc\_file type file
	}
	\code{type} defines the file type, which can be one of the
	following:
	\begin{description}

	\item[1D] Location of file defining one-dimensional profile.

	The file must contain three columns: Column 1 is the altitude in km,
	column 2 the liquid water content (LWC) in grams per cubic meter, and
	column 3 the effective droplet radius in micrometer. Empty lines are
	ignored. Comments start with \code{\#}.  Note that the definition of
	cloud altitudes in \code{wc\_file 1D} refers to sea level, not to altitude
	above ground. E.g., when \code{altitude} is set to 1.63km, and the
	first cloud level is defined at 3km, the cloud would start at 1.37km
	above ground.  An example of a water cloud is given in
	\file{examples/WC.DAT}.
	
	Per default the cloud properties are interpreted as layer
	properties. Before version 1.4 the default was level properties: The
	optical depth of a layer was calculated using information from the
	upper and lower levels defining the layer.  
	To switch to the old behaviour, use \code{interpret\_as\_level}. 
	See section~\ref{seq:wc} about water clouds for a
	realistic example how the contents of the \code{wc\_file 1D} are converted
	to optical properties.
	
	\ifthreedmystic{

	\item[3D] Define a MYSTIC 3D water cloud input file. 
	In case the \code{rte\_solver} is not \code{mystic},
	independent column calculations are performed for the 3D field.
	(See \file{examples/UVSPEC\_WC\_IPA.INP} for an example.)
	The expected format of the cloud file is: 
	\fcode{
	Nx  Ny  Nz flag \\
	dx  dy  z(1) z(2) ... z(n) \\
	ix  iy  iz  ext  g     ssa  (if flag == 1)\\
	ix  iy  iz  ext  reff       (if flag == 2)\\
	ix  iy  iz  LWC  reff       (if flag == 3)\\
	}
	where \code{Nx}, \code{Ny} and \code{Nz} are the number of grid boxes in 
	\code{x}, \code{y}, and \code{z}-direction.
	The parameter \code{flag} determines the format of the 3rd and 
	following lines. In the second line \code{dx} and \code{dy} are the sizes of the boxes in x- 
	and y-direction in km. 
	In the third and following lines the indices \code{ix}, \code{iy}, and \code{iz} specify cloudy pixels. 
	The optical properties of the cloud, are given by the other parameters in the line, 
	where \code{ext} is the extinction coefficient [1/km], \code{g} the asymmetry parameter,
	\code{reff} the effective radius [micrometer], and \code{ssa} the single scattering albedo.
	The conversion from microphysical to optical properties is defined by \code{wc\_properties} 
	(identical to the 1D case). Note that the first dimension (x)
	propagates east, and the second dimension (y) north. For more
	information see section~\ref{sec:mystic}.

%	\item[ipa] Define a 3D input file.
%	Independent column calculations are performed for the 3D field.
%	As argument a name of a 3D cloud file must be given. This file has to
%	be in the format as needed by MYSTIC, see \code{wc\_file 3D}.
%	(See \file{examples/UVSPEC\_WC\_IPA.INP} for an example.)
	}
	
	\item[ipa\_files] A two-column file, defining water cloud property files (see \code{wc\_file 1D}) in the first 
	column and the correspoding weights in the second column. 
	The radiative transfer calculation is performed independently for each
	cloud column and the result is the weighted average of all independent
	columns. If \code{ic\_file ipa\_files} and \code{wc\_file ipa\_files} are both
	defined, both must have the same columns in the same order, otherwise
	\code{uvspec} will complain. See
	\file{examples/UVSPEC\_WC\_IC\_IPA\_FILES.INP} for an example.
	
	\item[moments] A way to specify water cloud extinction coefficient, single
	scattering albedo, and scattering phase function for each layer.

	The file specified by 
	\code{wc\_file moments} has two columns where column 1 is the altitude in km. The second
	column is the name of a file which defines the optical properties of the layer 
	starting at the given altitude. The files specified in the second column must 
	have the following format:
	\begin{description}
	\item{Column 1: } 
	The wavelength in nm. These wavelengths may be different from those in \code{source solar filename}. 
	Optical properties are interpolated to the requested wavelengths.
	\item{Column 2:} 
	The extinction coefficient of the layer in units km-1. 
	\item{Column 3:} 
	The single scattering albedo of the layer.
	\item{Column 4-(nmom+4):} 
	The moments of the scattering phase function. 
	\end{description}
	Note that if using the \code{rte\_solver disort} or \code{rte\_solver
	fdisort2} it makes good sense to make the number of moments larger
	than \code{number\_of\_streams} because all moments are used in the calculation. For
	\code{rte\_solver fdisort1} and \code{rte\_solver polradtran} the number
	of moments included in the calculations will be \code{number\_of\_streams}+1. Higher
	order moments will be ignored for these solvers.  Please note that the
	uppermost line of the \code{wc\_file moments} denotes simply the top altitude
	of the uppermost layer. The optical properties of this line are
	consequently ignored. There are two options for this line: either an
	optical property file with zero optical thickness is specified
	or "NULL" instead.

	\end{description}
		
\option{wc_modify}
	Modify water cloud optical properties.
	\fcode{
	wc\_modify variable scale/set value
	}
	\code{variable} can be one of the following parameter:
	\begin{description}

	\parameter{gg}
	Modify the water cloud asymmetry factor for all wavelengths and altitudes.
	\begin{description}
	\item[set]
	\code{value} can be a float between -1.0 and 1.0. 
	Please note that this option is only applied if a Henyey-Greenstein 
	phase function is used but not if an explicit phase function is 
	defined e.g. with \code{wc\_file moments}. It doesn't make sense to modify only 
	the first moment of an explicit phase function.  
	This option is useful only for monochromatic 
	calculations or in wavelength regions where the optical 
	properties of water clouds can be considered constant, 
	e.g. the ultraviolet range. 
	\item[scale]
	Scale the water cloud asymmetry factor for all wavelengths and altitudes
	with \code{value} between 0.0 and 1.0.
	\end{description}
	
	\parameter{ssa}
	Modify the water cloud single scattering albedo for all wavelengths
	and altitudes.
	\begin{description}
	\item[set] \code{value} can be a float between 0.0 and 1.0. 
	This option is useful only for monochromatic 
	calculations or in wavelength regions where the optical properties of water clouds 
	can be considered constant, e.g. the ultraviolet range.
	\item[scale]
	Scale the water cloud single scattering albedo for all wavelengths and altitudes
	with \code{value} between 0.0 and 1.0.
	\end{description}

	\parameter{tau}
	Modify the total water cloud optical thickness.
	\begin{description}
	\item[set] Set optical thickness to a constant value for all wavelengths.
	The optical thickness defined here is the integral from the surface at the 
	user-defined \code{altitude} to TOA (top of atmosphere). This option is useful only 
	for monochromatic calculations or in wavelength regions where the optical properties 
	of water clouds can be considered constant, e.g. the ultraviolet range.
	\item[scale] Scale the water cloud optical thickness for all wavelengths and altitudes
	with \code{value} between 0.0 and 1000000.0. Also works for 3d clouds.
	\end{description}

	\parameter{tau550}
	Set the water cloud optical thickness at 550nm. 
	The optical thickness defined here 
	is the integral from the surface at the user-defined \code{altitude} 
	to TOA (top of atmosphere). Other wavelengths are scaled accordingly.
	Note that this option requires for technical reasons that the
	wavelength interval defined by \code{wavelength} does contain 550nm.
	\end{description}
	\fcode{
	wc\_modify tau550 set value
	}
			
\option{wc_properties}
	Define how liquid water content and effective radius are translated to optical properties.
	\fcode{
	wc\_properties property [interpolate]
	}
	Possible choices for \code{property} are:
	\begin{description}
	\parameter{hu} 
	Parameterization by \citet{Hu1993}; 
	this is the default setting. Note that 
	the parameterization is somewhat different for \code{mol\_abs\_param FU} than for all other 
	cases because in the latter case the parameterization from the newer (March 2000) 
	Fu and Liou code is used while otherwise the data are taken from the original 
	paper by \citet{Hu1993}. Note that this parameterization has been
	developed to calculate irradiances, hence it is less suitable for radiances.
	This is due to the use of the Henyey-Greenstein phase function as an approximation
	of the real Mie phase function.
	\parameter{echam4} 
	Use the very simple two-band parameterization of the ECHAM4 climate model, described 
	in \cite{Roeckner1996}; this is probably only meaningful if you want to compare
	your results with ECHAM4, the two bands are 0.2 - 0.68 micrometer and 0.68 - 4.0 micrometer;
	within these bands, the optical properties are assumed constant.
	\parameter{mie} 
	Use pre-calculated Mie tables; useful for \code{mol\_abs\_param}; 
	the tables are expected in \code{data\_files\_path}\file{/correlated\_k/.../}. \\
	For spectral or pseudo-spectral (\code{mol\_abs\_param sbdart}) calculations, 
	a set of pre-calculated tables is also available.
	%CE: commented some parts. Should be rewritten when new Mie tables are well
	% optimized.
	% the wavelength grid points of these 
	%data has been carefully selected such that the extinction cross section, 
	%single scattering albedo, and the asymmetry parameter are accurate to 1\% 
	%(compared to the fully-resolved Mie calculation) for all wavelengths
	%between 250nm and 100 micrometer. 
	For spectral or pseudo-spectral
	calculations optional argument \code{interpolate} has to be defined explicitely to
	initiate the interpolation of the optical properties to the internal wavelength grid.
	%Please note that this option may be extremely memory-consuming because for each 
	%internal wavelength a full set of Legendre moments of the phase function is 
	%stored (up to several thousands). 
	The Mie tables are not part of the standard distribution 
	(because of their large size) but they are freely available from http://www.libradtran.org. 
	This is the correct option to calculate radiances, to be preferred over the 
	Henyey-Greenstein approach of \citet{Hu1993}.
	\parameter{filename} 
	Read optical properties from specified filename; file format is as produced 
	by the \code{mie}-tool of the {\sl libRadtran} package (see \code{output\_user netcdf}).
	\end{description}

	With the optional argument \code{interpolate} the water cloud optical properties 
	are interpolated over wavelength; useful for precalculated optical property files. 
	Please note that this option may be extremely memory-consuming because for each internal wavelength 
	a full set of Legendre moments of the phase function is stored (up to several thousands). 
		
\option{z_interpolate}
\option{zout}
	This option is used to specify the output altitudes in km \emph{above surface altitude}. 
	One or more altitudes may be specified in increasing magnitude. 
	\fcode{
	zout 0 1 2 3 4 5 ...
	}
	Output altitudes must be within the range defined in the
	\code{atmosphere\_file}. Note that \code{zout} does not restructure the
	atmosphere model. Hence, if you specify \code{zout 0.730} and have your
	atmosphere model in \code{atmosphere\_file} go all the way down to 
	sea level, i.e. 0.0km., output is presented at 0.730km and calculations
	performed with an atmosphere between 0.0 and 0.730 km (and above of course). 
	If you want calculations done for e.g. an elevated site you have to restructure
	the atmosphere model and make sure it stops at the appropriate altitude.
	This you may either due by editing the atmosphere file or by using 
	\code{altitude}. Note that for \code{rte\_solver polradtran} the atmosphere file must 
	contain the altitudes specified by \code{zout}. You can also use \code{toa} for 
	\code{top of atmosphere} and \code{sur} for \code{surface altitude} and \code{cpt} 
	for \code{cold point tropopause}.

	Instead of specifying the altitudes in km, it is also possible to use keywords as 
	argument for this option. Possible keywords are \code{atm\_levels}, 
	\code{all\_levels}, \code{model\_levels}, \code{model\_layers}, and \code{model\_levels\_and\_layers}. 
	For \code{atm\_levels}, all levels from the \code{atmosphere\_file} are used as output levels.  
	For \code{all\_levels}, all levels (including levels from \code{atmosphere\_file}, 
	\code{mol\_file}, cloud files, \code{altitude} options) are used as output levels. 
	For \code{model\_levels}, \code{model\_layers}, \code{model\_levels\_and\_layers} 
	the levels, layers, or both from the \code{ECMWF\_atmosphere\_file} are used as output level.
	Usage e.g.:
	\fcode{
	zout model\_levels [nlev\_max]
	}
	With the optional argument \code{nlev\_max} the user may specify the number of zout layers 
	from the ground.
		
\option{zout_interpolate}
	The z-grid of optical properties is determined by the \code{atmosphere\_file},
	and, if specified, by other profile files like \code{mol\_file}, \code{rh\_file}, 
	or \code{refractive\_index\_file}. 
	Additional levels might be introduced by the \code{zout} 
	option and the second argument of the \code{altitude} option. By default 
	(if \code{zout\_interpolate} is not specified) levels introduced 
	by the \code{zout} option will not affect the optical property 
	profiles, that is, the optical properties are constant within the layers specified by 
	the \code{atmosphere\_file} and profile files.
	If \code{zout\_interpolate} is specified, the atmospheric profiles (tracegases, temperature ...) 
	are interpolated to the levels introduced by \code{zout}, and optical 
	properties are determined from the interpolated atmospheric properties.
	If \code{heating\_rate}, \code{rte\_solver polradtran}, \code{rte\_solver rodents},
	or \code{rte\_solver twostrebe}
	is specified, \code{zout\_interpolate} 
	will also be automatically activated. 
	\code{zout\_interpolate} generally causes smoother variation of the optical properties.
		
\option{zout_sea}
	like zout, but \emph{above sea surface}
		
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

  Input.UTC.tm_yday = NOT_DEFINED_INTEGER;
  Input.aer.alpha_0 = NOT_DEFINED_FLOAT;
  Input.aer.alpha_1 = NOT_DEFINED_FLOAT;
  Input.aer.alpha_2 = NOT_DEFINED_FLOAT;
  Input.aer.beta = NOT_DEFINED_FLOAT;
  Input.aer.n_species = NOT_DEFINED_INTEGER;
  Input.aer.nimag = NOT_DEFINED_FLOAT;
  Input.aer.spec = 0;
  Input.aer.tau_wvl_lambda = NOT_DEFINED_FLOAT;
  Input.aer.tau_wvl_tau = NOT_DEFINED_FLOAT;
  Input.aer.visibility = NOT_DEFINED_FLOAT;
  Input.alb.albedo = 0.0;
  Input.alb.source = NOT_DEFINED_INTEGER;
  Input.alb.surf_type_map = FALSE;
  Input.alb.surface = NOT_DEFINED_INTEGER;
  Input.alt.altitude = NOT_DEFINED_FLOAT;
  Input.alt.altitude_dz = 0.0;
  Input.alt.source = ALT_NOT_DEFINED;
  Input.atm.sza = NOT_DEFINED_FLOAT;
  Input.atm.sza_source = SZA_BY_TIME_AND_LOCATION;
  Input.atm.zout_interpolate = NO_ZOUT_INTERPOLATE;
  Input.cm.solar_wind = 0;
  Input.convolve = 0;
  Input.disort2_brdf = BRDF_NONE;
  Input.flu.fluorescence = 0.0;
  Input.flu.source = NOT_DEFINED_INTEGER;
  Input.hapke.source = HAPKE_CONSTANT;
  Input.ipa = 0;
  Input.pressure = NOT_DEFINED_FLOAT;
  Input.quiet = 0;
  Input.r_earth = 6370.0;
  Input.rayleigh_depol = NOT_DEFINED_FLOAT;
  Input.rossli.hotspot = BRDF_ROSSLI_HOTSPOT_OFF;
  Input.rossli.source = ROSSLI_CONSTANT;
  Input.rpv.source = RPV_CONSTANT;
  Input.rte.deltam = SWITCH_ON;
  Input.rte.mc.alpha = 5;
  Input.rte.mc.azimuth = MCAZIMUTH_NEW;
  Input.rte.mc.backward.islower = NOT_DEFINED_INTEGER;
  Input.rte.mc.backward.isupper = NOT_DEFINED_INTEGER;
  Input.rte.mc.backward.jslower = NOT_DEFINED_INTEGER;
  Input.rte.mc.backward.jsupper = NOT_DEFINED_INTEGER;
  Input.rte.mc.coherent_backscatter = 0;
  Input.rte.mc.concentration_is = FALSE;
  Input.rte.mc.dt = 1e-07;
  Input.rte.mc.ipa = 0;
  Input.rte.mc.ixmax = NOT_DEFINED_INTEGER;
  Input.rte.mc.ixmin = NOT_DEFINED_INTEGER;
  Input.rte.mc.iymax = NOT_DEFINED_INTEGER;
  Input.rte.mc.iymin = NOT_DEFINED_INTEGER;
  Input.rte.mc.maxscatters = NOT_DEFINED_INTEGER;
  Input.rte.mc.polarisation = 0;
  Input.rte.mc.polarisation_state = 0;
  Input.rte.mc.reflectalways = 0;
  Input.rte.mc.refraction = 0;
  Input.rte.mc.spectral = 0;
  Input.rte.mc.std = 0;
  Input.rte.nstr = 6;
  Input.rte.pol_max_delta_tau = 1e-05;
  Input.rte.reverse = 0;
  Input.spline = 0;
  Input.surface_temperature = NOT_DEFINED_FLOAT;
  Input.test_optical_properties = 0;
  Input.verbose = 0;
  Input.wl.end = NOT_DEFINED_FLOAT;
  Input.wl.end_index = NOT_DEFINED_INTEGER;
  Input.wl.start = NOT_DEFINED_FLOAT;
  Input.wl.start_index = NOT_DEFINED_INTEGER;





  
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


