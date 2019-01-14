/*!
  \file caa_main_model2.c
  \brief Containing the main routines for fitting wgl model
  \author Hanne Rognebakke

  Fitting the wgl models are performed in the ::caa_main_model2 routine. The
  main estimation is done in the :main_model2 routine.
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
  \brief Fitting wgl model through MCMC simulations.
  \author Hanne Rognebakke

  This routine gets all the input parameter from binary files. The input is first stored in
  input-structures. Then the routine ::main_model2 is called, where all the approperiate
  c-structures are made and other routines for performing the MCMC simulations
  are called.

  The simulations are printed to a binary file. Simulated variables are only
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
  
  /* Variables connected to input data */
  Input_common *inCommon;
  Input_wgl    *inWgl;
  Input_prior  *inPrior;
  Data_orig    *D_orig;  
  Data_CC      *D_CC=NULL;       /* Parameters for coastal cod */
  
  if(argc!= 2)
    {
      fprintf(stderr,"Usage: caa_main_model2 fit.directory\n");
      exit(1);
    }
  sscanf(argv[1],"%s",dirname);


  sprintf(filename,"%slog.txt",dirname);
  fp_log = fopen(filename,"w");
  fprintf(fp_log,"Start caa_main_model2.c\n");
  
  #ifdef LOG_FILE
  sprintf(buffer,"%s/caa_logfile_model2.txt",dirname);
  g_caa_log = fopen(buffer,"w");
  fprintf(g_caa_log,"Start C program: estimating weight-given-length model\n");
  #endif
  fprintf(stderr,"Start C program: estimating weight-given-length model\n");


  /* Allocate space for objects */
  err = alloc_objects_wgl(&inCommon, &D_orig, &inWgl, &D_CC, &inPrior);
  if(err)
    {
      write_warning("caa_main_model2:Error calling alloc_objects_wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  /* Make struct for input data */
  sprintf(filename,"%scommon_par_fit_ascii",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_common_ascii(filename, inCommon, D_orig, D_CC);
  if(err)
    {
      write_warning("caa_main_model2:Error calling readdata_common\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  /* Make struct for wgl data */
  sprintf(filename,"%sstoxdata_wgl",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_wgl(filename, D_orig, inWgl);
  if(err)
    {
      write_warning("caa_main_model2:Error calling readdata_wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  /* Make struct for prior data */
  sprintf(filename,"%spriorList",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_prior(filename, inPrior);
  if(err)
    {
      write_warning("caa_main_model2:Error calling readdata_prior\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
  
  err = add_object_info_wgl(inCommon, D_orig, inWgl);
  if(err)
    {
      write_warning("caa_main_model2:Error calling add_object_info_wgl\n");
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
  err = main_model2(inCommon, inWgl, inPrior, D_orig, 0);
  if(err)
    {
      write_warning("caa_main_model2:Error calling main_model2\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

   /* Clean up */
  #ifdef DEBUG_PROG
  fprintf(stderr,"caa_main_model1.c: CLEAN UP!!\n");
  #endif
if(0){
  err = re_alloc_objects_wgl(&inCommon, &D_orig, &inWgl, &D_CC, &inPrior);  
  if(err)
    {
      write_warning("caa_main_model2:Error calling re_alloc_objects_wgl\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
 }  

  #ifdef LOG_FILE
  fclose(g_caa_log);
  #endif


  fprintf(fp_log,"End caa_main_model2.c\n");
  fprintf(fp_log,"OK\n");
  fclose(fp_log);

  return(0);
}		/* end of caa_main_model2 */
  
