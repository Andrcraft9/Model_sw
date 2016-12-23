C(C) DIANSKY N.A. (dinar@inm.ras.ru), begin at 16.10.99 16:15
C INCLUDE MODULE FOR DESCRIPTION COMMON OGCM VARIABLES AND
C                                TASK CONTROL PARAMETERS
C
C STREAM FUNCTION ARRAYS
      REAL*8   EPS !PRECISION PARAMETER FOR FUNCTION OR SEA LEVEL ITERATION
      INTEGER NSMAX  !MAX NUMBER OF ITERATION
      PARAMETER(EPS=0.0005,NSMAX=999)      !for stream function
C     PARAMETER(EPS=0.00001, NSMAX=400)    !for sea level (GMRS)
C RDR - BOTTOM FRICTION FOR STREAM FUNCTION
C RMU - ITERATION PARAMETERS FOR STREAM FUNCTION OR SEA LEVEL ITERATION
C RIGID_LID = 1 - FOR RIGID LID CONDITION;
C                 RIGID_LID=0/1 IN TERMS "MEAN PRESSURE"/"SLH"
C PJEBARRDR = 1 - FOR JEBAR: =0, WITHOUT JEBAR: =1
      COMMON /SFSLPAR/RDR,RMU,RIGID_LID,PJEBAR
      REAL*8          RDR,RMU,RIGID_LID,PJEBAR
C BASE 2D-3D-ARRAYS:
      
C      INCLUDE '2SHALLOWATER.INC'

C---------------------------------------------------------------------72
C BASE 3D-ARRAYS:
      INCLUDE '2HYDRODYN.INC'           

C+++++++++++++++++++PARAMETER TO DEFINE TASK++++++++++++++++++++++++++72

      INCLUDE '2TASKCNTRL.INC'

C SEA SURFACE BOUNDARY CONDITION

      INCLUDE '2SSBCOND.INC'

C######################################################################
C DIFFUSION PARAMETERS

      INCLUDE '2DIFFCNTRL.INC'

C---------------------------------------------------------------------72
C COEFFICIENTS OF VISCOSITY  AND DIFFUSION

      INCLUDE '2VISDIFFCOEFF.INC'

C######################################################################
C SEA SURFACE BOUNDARY CONDITION
      
      INCLUDE '2SEASURF.INC'

C######################################################################
      COMMON /WORK/ WOR(NX,NY), ACCUM(NX,NY,9)
      REAL ACCUM,WOR

C Connection ATM TO ICE 
      INCLUDE '2ATMTOICE.INC'

      INCLUDE '2PASSTRAC.INC'

