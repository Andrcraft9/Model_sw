# Need before make:
# export PETSC_DIR=/path/to/petsc and export PETSC_ARCH=petsc-arch
#----------------------------------------------------------------------------

FFLAGS= -cpp -Wtabs ${PETSC_FC_INCLUDES} -I./Inc
CFLAGS= ${PETSC_CC_INCLUDES}
FCFLAGS= -I./Inc

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules




SRCCONTROL =	\
	Control/ctlwrite.f \
	Control/ocforcag.f \
	Control/ocforcog.f \
	Control/ocini.f    \
	Control/ocoutput.f \
	Control/ocstep.f   \
	Control/setlqbc.f 

SRCSERVICE =	    \
	Service/basinpar.f       \
	Service/gridcon.f        \
	Service/inoutd.f         \
	Service/ocalg.f          \
	Service/rw_ctl.f         \
	Service/rwparam.f        \
	Service/s_timu.f         \
	Service/s_treat.f        \
	Service/sigmaztr.f       \
	Service/vgrid.f          
         

SRCFUNCTION =	\
	Function/m_adapt.f  \
	Function/m_bottomfr.f \
	Function/m_divi4.f    \
	Function/m_seaice.f   \
	Function/m_sweeps.f      \
	Function/m_trants.f      \
	Function/m_tranuv.f      \
	Function/m_unidif.f      \
	Function/m_uvwfun.f      \
	Function/m_vertmo.f      \
	Function/m_np_filters.f      
 

compile: mpi_parallel_tools.o mod_shallowater.o m_sl_no_split.o octask.o 
	-${FLINKER} -o octask mpi_parallel_tools.o octask.o mod_shallowater.o $(SRCFUNCTION) m_sl_no_split.o \
	$(SRCCONTROL) $(SRCSERVICE) ${PETSC_KSP_LIB} 


clean_o: 
	rm *.o
	rm *.mod

run:
	-@${MPIEXEC} -n 16 ./octask -ksp_type gmres
