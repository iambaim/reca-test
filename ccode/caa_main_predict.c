/*!
  \file caa_main_predict.c
  \brief Containing the outer main routine for predicting catch at age
  \author Hanne Rognebakke

  The main estimation is done in the ::predict routine.
*/

/*Include Files:*/
#include "caa.h"
#include "caa_predict.h"
#include "caa_read_binary.h"
#include "caa_read_write.h"

#ifdef LOG_FILE
extern FILE     *g_caa_log; 
#endif
 

/*!
  \author Hanne Rognebakke
  \brief Estimates catch-at-age for simulated parameters and random effects.

  This routine gets all the input parameter from binary files. The input is first stored in
  input-structures. Then the routine ::predict is called, where all the approperiate
  c-structures are made and other routines for performing the MCMC simulations
  are called.

  The predictions are printed to a binary file with all
  predictions from one iteration in one sequential block. See the
  ::write_it_totcatch routine for the format of this block.
*/
int main(int argc, char *argv[])
{
  char     dirname[MAX_STR];
  char     filename[MAX_STR];
  char     buffer[MAX_STR];
  int      err = 0;
  FILE    *fp_log;

  /* Variables connected to input data*/
  Input_cell     *inCell;
  Input_totcatch *inCatch;
  Input_predict  *inPredict;

  if(argc!= 2)
    {
      fprintf(stderr,"Usage: caa_main_predict predict.directory\n");
      exit(1);
    }
  sscanf(argv[1],"%s",dirname);

  sprintf(filename,"%slog.txt",dirname);
  fp_log = fopen(filename,"w");
  fprintf(fp_log,"Start caa_main_predict.c\n");

  #ifdef LOG_FILE
  sprintf(buffer,"%s/caa_logfile_predict.txt",dirname);
  g_caa_log = fopen(buffer,"w");
  fprintf(g_caa_log,"Start C program: predicting catch-at-age\n");
  #endif
  fprintf(stderr,"Start C program: predicting catch-at-age\n");

  /* Make struct for input data */
  sprintf(filename,"%scommon_par_predict_ascii",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_predict_ascii(filename, &inPredict);
  if(err)
    {
      write_warning("caa_main_predict:Error calling readdata_predict\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }


  sprintf(filename,"%sdata_catch",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_catch(filename, &inCatch);
  if(err)
    {
      write_warning("caa_main_predict:Error calling readdata_catch\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }
  sprintf(filename,"%sdist_cell",dirname);
  fprintf(stderr,"Read input data from file: %s\n",filename);
  err = readdata_cell(filename, &inCell);
  if(err)
    {
      write_warning("caa_main_predict:Error calling readdata_cell\n");
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
  err = predict(inPredict, inCatch, inCell, NULL);
  if(err)
    {
      write_warning("caa_main_predict:Error calling predict\n");
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
  err = re_readdata_cell(&inCell);  
  if(err)
    {
      write_warning("caa_main_predict:Error calling re_readdata_cell\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  err = re_readdata_catch(&inCatch);  
  if(err)
    {
      write_warning("caa_main_predict:Error calling re_readdata_catch\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  err = re_readdata_predict(&inPredict);  
  if(err)
    {
      write_warning("caa_main_predict:Error calling re_readdata_predict\n");
      #ifdef LOG_FILE
      fclose(g_caa_log);
      #endif
      fclose(fp_log);
      return(err);
    }

  
  #ifdef LOG_FILE
  fprintf(g_caa_log,"caa_main_predict finished\n");
  fclose(g_caa_log);
  #endif


  fprintf(fp_log,"End caa_main_predict.c\n");
  fprintf(fp_log,"OK\n");
  fclose(fp_log);

  return(0);
}		/* end of caa_main_predict */
