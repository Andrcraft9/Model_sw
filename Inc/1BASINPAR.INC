C PARAMETERS OF THE OCEAN MODEL

C--------------------------------------------------------------------
C NX, NY - NUMBER OF TOTAL POINTS IN X-, Y-DIRECTION IN BASE ARRAYS
C NZ     - NUMBER OF T-LEVELS
C MMM,NNN - BEGINNING AND MM,NN - END OF CALCULATED T-POINTS.
C MDIM = MAX(NX,NY,NZ)
C MMD -  CICLIC PERIOD IN NUMBER OF GRID POINTS
C NBASINS - NUMBER OF REGIONAL BASINS(FIRST- WHOLE BASIN)
C FOR NON-PERIODIC CASE: SET PARAMETER AS MMD=0
      INTEGER NX, NY , NZ, MLR, MMM, NNN, MM, NN, MDIM
      INTEGER MMD, NBASINS
      PARAMETER (NX=208,NY=108,NZ=1,MLR=13,
     &           MMM=3, NNN=3, MM=NX-2,NN=NY-2,MDIM=MAX(NX,NY,NZ),
     &           MMD=0,NBASINS=1)

C--------------------------------------------------------------------
C RLAT - FIRST CALCUTABLE T-POINT (M=MMM) ON LATITUDE  (IN DEGREES)
C RLON - FIRST CALCUTABLE T-POINT (N=NNN) ON LONGITUDE (IN DEGREES)
C DYST & DXST - LATITUDE & LONGITUDE STEPS (IN DEGREES)
      REAL      RLAT, RLON, DXST, DYST
      PARAMETER(RLON=0.05,RLAT =42.55,
     &          DXST =0.05,DYST =0.05	)



C Parameters of rotation
      REAL ROTATION_ON_LON,ROTATION_ON_LAT
      PARAMETER (ROTATION_ON_LON = 0.0, ROTATION_ON_LAT = 0.0)

CCCC   ROTATION_ON LON, ROTATION_ON_LAT - geographic coordinates of
CCCC   zero of latitude and longitude on rotated grid


      REAL X_POLE,          !LON OF REAL NORTH POLE IN SPEC. SYSTEM
     &     Y_POLE,          !LAT OF REAL NORTH POLE IN SPEC. SYSTEM
     &     P_POLE,          !LON OF NEW NORT POLE IN GEO SYSTEM
     &     Q_POLE           !LAT OF NEW NORT POLE IN GEO SYSTEM


      PARAMETER(X_POLE=100.0,          !LON OF REAL NORTH POLE IN SPEC. SYSTEM
     &          Y_POLE=70.0,          !LAT OF REAL NORTH POLE IN SPEC. SYSTEM
     &          P_POLE=X_POLE,        !LON OF NEW NORT POLE IN GEO SYSTEM
     &          Q_POLE=Y_POLE)        !LAT OF NEW NORT POLE IN GEO SYSTEM
C----------------------------------------------------------------------
C NORTH POLAR FILTER PARAMETERS:

      REAL CRIT_LAT   !CRITICAL LATITUDE (DGR) FOR BEGIINING OF FILTRATION

      PARAMETER(CRIT_LAT = 180.00)

      INTEGER
     &  NUM_LAT_BEG_FILTER, !NUMBER OF LATITUDE ON T-GRID FOR BEGINNING FILTERING
     &  MAX_LONG_SFLT,      !MAX NUMBER OF COEFFICIENTS OF SIMMETRIC FILTER
     &  NUM_NPSF            !NUMBER OF GRIDED LATITUDE CIRCELE FOR NORTH POLE FILTRING

      PARAMETER(NUM_LAT_BEG_FILTER = MAX(NNN,
     &                 MIN(NINT((CRIT_LAT - RLAT)/DYST) +NNN,NN+1) ),
     &               MAX_LONG_SFLT = MAX((MM-MMM+1)/2-1,1),
     &                    NUM_NPSF = MAX(NN-NUM_LAT_BEG_FILTER+1,1)    )

C-----------------------------------------------------------------------

      REAL LAT_CRIT_4D       !CRITICAL LATITUDE (DGR) FOR BEGINING OF 4TH ORDER FILTRATION

      PARAMETER(LAT_CRIT_4D=180.0)

C-----------------------------------------------------------------------
