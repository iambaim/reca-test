#INCLUDE = -I/usr/local/include -I/nr/prog/R/2.10.1_UBUNTU8.04/lib64/R/include/
INCLUDE = -I/usr/local/include -I/usr/lib 
#CFLAGS = -g -Wall -pg -O3 $(INCLUDE)
#DEBUG_CFLAGS = -DLOG_FILE -DDEBUG_INPUT -DDEBUG_PROG -DDEBUG_GAUSS_FILE 
DEBUG_CFLAGS =
#CFLAGS = $(DEBUG_CFLAGS) -g -pg -O0 $(INCLUDE) -L. 
CFLAGS = $(DEBUG_CFLAGS) -g -O0 $(INCLUDE) -L. 

LIB    = -ltaucs -lmetis -llapack -lblas -lgfortran -lquadmath -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm 
LFLAGS = -L/usr/lib -L../libs -L/usr/lib/atlas-base/atlas -L/usr/lib/atlas-base -static
BIN    = .
CC     = gcc -Wdiv-by-zero 
LD     = gcc 


TEST   = test.o 

OBJ1   = caa_main_model1.o 

OBJ2   = caa_main_model2.o 

OBJ3   = caa_main_predict.o 

OBJ  =  caa_predict.o \
	caa_estimate.o \
	caa_cell_constr.o \
        caa_chol.o \
        caa_evaluate.o \
	caa_init.o \
	caa_input.o \
	caa_lqp.o \
        caa_mcmc.o \
	caa_read_binary.o \
	caa_read_write.o \
        caa_routines.o \
        caa_sample_g_a.o \
        caa_sample_gauss.o \
        caa_sample_multi.o \
        caa_summaries.o \
        caa_util.o \
        caa_utl.o \
        blockupdate.o \
        com.o \
        d_brent.o \
        d_erfcc.o \
        d_erff.o \
        d_f1dim.o \
        d_gammln.o \
        d_gammp.o \
        d_gauher.o \
        d_gcf.o \
        d_gser.o \
        distributions.o \
	d_ksone.o \
        d_linmin.o \
        d_mnbrak.o \
        d_nrutil.o \
        d_powell.o \
	d_probks.o \
	d_sort.o \
        error-handler.o \
        gdens.o \
        globals.o \
        graph.o \
        hash.o \
        hidden-approx.o \
        lapack-interface.o \
        linpack.o \
        optimize.o \
        problem-setup.o \
        random.o \
        ranlib.o \
        rw.o \
        smtp-band.o \
        smtp-profile.o \
        smtp-taucs.o \
        sparse-interface.o \
        tabulate-Qfunc.o \
        timer.o \
        utl_malloc.o \
        wa.o \
	acm582.o \
        dcdflib.o \
        GMRFLib-fortran.o 

#.c.o : $(CC) -c $(CFLAGS) $<


caa: caa_main_model1 caa_main_model2 caa_main_predict

caa_main_model1:	$(OBJ1) $(OBJ)
	$(LD) $(CFLAGS)  $(OBJ1) $(OBJ) -o $(BIN)/caa_main_model1 $(LFLAGS) $(LIB)

caa_main_model2:	$(OBJ2) $(OBJ)
	$(LD) $(CFLAGS)  $(OBJ2) $(OBJ) -o $(BIN)/caa_main_model2 $(LFLAGS) $(LIB)

caa_main_predict:	$(OBJ3) $(OBJ)
	$(LD) $(CFLAGS)  $(OBJ3) $(OBJ) -o $(BIN)/caa_main_predict $(LFLAGS) $(LIB)
#	$(LD) $(CFLAGS)   $(OBJ) -o $(BIN)/predict $(LFLAGS) $(LIB)

caa_main_model1.o: caa.h caa_read_binary.h caa_estimate.h caa_main_model1.c
caa_main_model2.o: caa.h caa_read_binary.h caa_estimate.h caa_main_model2.c
caa_main_predict.o: caa.h caa_read_binary.h caa_predict.h caa_main_predict.c
caa_read_write.o: caa.h caa_read_write.h caa_read_write.c
caa_read_binary.o: caa.h caa_read_binary.h caa_read_binary.c
caa_routines.o: caa.h caa_routines.h caa_routines.c
caa_sample_gauss.o: caa.h caa_sample_gauss.h caa_sample_gauss.c
caa_sample_multi.o: caa.h caa_sample_multi.h caa_sample_multi.c
caa_sample_g_a.o: caa.h caa_sample_g_a.h caa_sample_g_a.c
caa_simulate.o: caa.h caa_simulate.h caa_simulate.c
caa_evaluate.o: caa.h caa_evaluate.h numrec.h caa_evaluate.c
caa_chol.o: caa.h caa_chol.h caa_chol.c
caa_lqp.o: caa.h caa_lqp.h caa_lqp.c
caa_util.o: caa.h caa_util.c
caa_utl.o: caa.h caa_utl.c
caa_init.o: caa.h caa_init.c
caa_input.o: caa.h caa_input.c
caa_predict.o: caa.h caa_predict.c
caa_summaries.o: caa.h caa_summaries.c
d_erff.o: d_erff.c
d_erfcc.o: d_erfcc.c
d_gammp.o: d_gammp.c
d_gser.o: d_gser.c
d_gcf.o: d_gcf.c
d_gammln.o: d_gammln.c
d_ksone.o: d_ksone.c
d_gauher.o: d_gauher.c
ranlib.o: ranlib.c
com.o: com.c
linpack.o: linpack.c
utl_malloc.o: utl_malloc.c
blockupdate.o: blockupdate.c
distributions.o: distributions.c
error-handler.o: error-handler.c
gdens.o: gdens.c
globals.o: globals.c
graph.o: graph.c
hidden-approx.o: hidden-approx.c
lapack-interface.o: lapack-interface.c
optimize.o: optimize.c
tabulate-Qfunc.o: tabulate-Qfunc.c
problem-setup.o: problem-setup.c
random.o: random.c
rw.o:rw.c
hash.o: hash.c
smtp-band.o:smtp-band.c
smtp-profile.o:smtp-profile.c
smtp-taucs.o:smtp-taucs.c
sparse-interface.o:sparse-interface.c
timer.o:timer.c
wa.o: wa.c
acm582.o: acm582.f
GMRFLib-fortran.o: GMRFLib-fortran.f
dcdflib.o: dcdflib.f
d_nrutil.o: d_nrutil.c
d_powell.o: d_powell.c
d_linmin.o: d_linmin.c
d_mnbrak.o: d_mnbrak.c
d_f1dim.o: d_f1dim.c
d_brent.o: d_brent.c

clean:
	rm -f *.o 
