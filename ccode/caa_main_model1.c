/*!
  \file caa_main_model1.c
  \brief Containing the outer main routine for fitting age and lga model
  \author Hanne Rognebakke

  Fitting the age and lga models are performed in the :main_model1 routine.
*/

/*Include Files:*/
#include "caa.h"
#include "caa_estimate.h"
#include "caa_read_binary.h"
#include "caa_read_write.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif
 





/*!
  \brief Fitting age and lga model through MCMC simulations.
  \author Hanne Rognebakke

  This routine gets all the input parameter from binary files. The input is first stored in
  input-structures. Then the routine ::main_model1 is called, where all the approperiate
  c-structures are made and other routines for performing the MCMC simulations
  are called.

  The simulations are printed to a binary file with all age parameters
  printed first, then all lga parameters. Simulated variables are only
  printed for every \em num_it_inner simulation after \em burn_in burnin iterations.
  For a further description, see the ::write_it routine in caa_routines.c
*/
int main(int argc, char *argv[])
{
  char     dirname[MAX_STR];
  char     filename[MAX_STR];
  char     buffer[MAX_STR];
  int      err = 0;
  FILE    *fp_log;

  /* Variables connected to input data*/
  Input_common *inCommon;
  Input_age    *inAge;
  Input_lga    *inLga;
  Input_wgl    *inHsz=NULL;
  Input_prior  *inPrior;
  Data_orig    *D_orig;         
  Data_CC      *D_CC=NULL;       /* Parameters for coastal cod */


  if(argc!= 2)
    {
      fprintf(stderr,"Usage: caa_main_model1 fit.directory\n");
      exit(1);
    }
  sscanf(argv[1],"%s",dirname);

  sprintf(filename,"%slog.txt",dirname);
  fp_log = fopen(filename,"w");
  fprintf(fp_log,"Start caa_main_model1.c\n");
  
  #ifdef LOG_FILE
  sprintf(buffer,"%s/caa_logfile_model1.txt",dirname);
  g_caa_log = fopen(buffer,"w");
  fprintf(g_caa_log,"Start C program: estimating age and length-given-age model\n");
  #endif
  fprintf(stderr,"Start C program: estimating age and length-given-age model\n\n");

  /* Allocate space for objects */
  err = alloc_objects_age_lga(&inCommon, &D_orig, &inAge, &inLga, &D_CC, &inHsz, &inPrior);
  if(err)
    {
      write_warning("caa_main_model1:Error calling alloc_objects\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
  
  /* Make struct for input data */
  sprintf(filename,"%scommon_par_fit_ascii", dirname);
  fprintf(stderr,"Read input data from file: %s\n", filename);
  err = readdata_common_ascii(filename, inCommon, D_orig, D_CC);
  if(err)
    {
      write_warning("caa_main_model1:Error calling readdata_common\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  /* Make struct for age data */
  sprintf(filename,"%sstoxdata_age", dirname);
  fprintf(stderr,"Read input data from file: %s\n", filename);
  err = readdata_age_lga(filename, D_orig, inAge, inLga);
  if(err)
    {
      write_warning("caa_main_model1:Error calling readdata_age\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  /* Make struct for lga data */
  sprintf(filename,"%sstoxdata_lga", dirname);
  fprintf(stderr,"Read input data from file: %s\n", filename);
  err = readdata_lga(filename, D_orig, inLga);
  if(err)
    {
      write_warning("caa_main_model1:Error calling readdata_lga\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
  #ifdef DEBUG_PROG
  int i;
  fprintf(stderr,"age:ncov=%d\n",inAge->cov->n_cov);
  for(i=0;i<inAge->cov->n_cov;i++)
    fprintf(stderr,"i=%d,nlev=%d\n",i,inAge->cov->n_lev[i]);
  fprintf(stderr,"lga int:ncov=%d\n",inLga->int_cov->n_cov);
  for(i=0;i<inLga->int_cov->n_cov;i++)
    fprintf(stderr,"i=%d,nlev=%d\n",i,inLga->int_cov->n_lev[i]);
  fprintf(stderr,"lga slp:ncov=%d\n",inLga->slp_cov->n_cov);
  for(i=0;i<inLga->slp_cov->n_cov;i++)
    fprintf(stderr,"i=%d,nlev=%d\n",i,inLga->slp_cov->n_lev[i]);
  #endif
  
  if(inCommon->inc_hsz==1)
    {
      /* Make struct for haulsize data */
      sprintf(filename,"%sstoxdata_hsz",dirname);
      fprintf(stderr,"Read input data from file: %s\n", filename);
      err = readdata_hsz(filename, D_orig, inHsz);
      if(err)
	{
	  write_warning("caa_main_model1:Error calling readdata_hsz\n");
          #ifdef LOG_FILE
	  fclose(g_caa_log);
          #endif
	  fclose(fp_log);
	  return(err);
	}
    }

  /* Make struct for prior data */
  sprintf(filename,"%spriorList",dirname);
  fprintf(stderr,"Read input data from file: %s\n", filename);
  err = readdata_prior(filename, inPrior);
  if(err)
    {
      write_warning("caa_main_model1:Error calling read_prior_age_lga\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  err = add_object_info_age_lga(inCommon, D_orig, inAge, inLga, inHsz);
  if(err)
    {
      write_warning("caa_main_model1:Error calling add_object_info\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Run main c-program\n");
  #endif
  /* Run the main c-program */
  err = main_model1(inCommon, inAge, inLga, inHsz, inPrior, D_orig, D_CC);
  if(err)
    {
      write_warning("caa_main_model1:Error calling main_model1\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
  

  /* Clean up */
  #ifdef LOG_FILE
  fprintf(g_caa_log,"Clean up\n");
  #endif

  #ifdef DEBUG_PROG
  fprintf(stderr,"caa_main_model1.c: CLEAN UP!!\n");
  #endif
if(0){
  err = re_alloc_objects_age_lga(&inCommon, &D_orig, &inAge, &inLga, &D_CC, &inHsz, &inPrior);  
  if(err)
    {
      write_warning("caa_main_model1:Error calling re_alloc_objects\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
 }  

  #ifdef LOG_FILE
  fprintf(g_caa_log,"caa_main_model1 finished\n");
  fclose(g_caa_log);
  #endif


  fprintf(fp_log,"End caa_main_model1.c\n");
  fprintf(fp_log,"OK\n");
  fclose(fp_log);

  return(0);
}		/* end of caa_main_model1 */
  
