C################## ATMOSPHERIC PARAMETERS #############################
      INTEGER   NXA,NYA, NOFMONTH
      REAL      DXA,DYA, X0A, Y0A
      INTEGER   MMMA,MMA,NNNA,NNA
      PARAMETER(
     &      NOFMONTH=4,        !DIM. OF TIME LEVEL FOR TIME INTRPOALION
     &           NXA=194,      !ATM. ARRAY DIM. ON LON
     &           NYA= 95,      !ATM. ARRAY DIM. ON LAT
     &           DXA=1.875,    !STEP ON LON (USED IN CASE OF REGULAR GRID)
     &           DYA=1.875,    !STEP ON LAT (-#-)
     &           X0A= 0.0,     !INITIAL LON POINT FOR M=MMMA (-#-)
     &           Y0A=-88.5420, !INITIAL LAT POINT FOR N=NNNA (-#-)
     &          MMMA=2,        !LEFFT BOUNDARY OF USED POINTS
     &           MMA=NXA-1,    !RIGHT BOUNDARY OF USED POINTS
     &          NNNA=2,        !DOWN  BOUNDARY OF USED POINTS
     &           NNA=NYA)      !UPPER BOUNDARY OF USED POINTS

      REAL WPM2_GCMPS,MMPDAY_CMPS,UPRH0
      PARAMETER (
     &  WPM2_GCMPS=1.0/CPWRH0*1.0E-04,!CONVERSION FACTOR [W/M**2]=>[GR*CM/S]
     & MMPDAY_CMPS=1.0/8.64E+05,      !CONVERSION FACTOR [MM/DAY]=>[CM/S]
     &       UPRH0=1.0/RH0)           !DIVIDE BY RH0
C   HERE CPWRH0, RH0 - USED FROM 0COM.INC 

C   CONVERSION FACTORS FOR ATMOSPHERIC DATA      

      REAL  CONFAC_QBAL,       !CONVFACTOR FOR HEAT BALANCE 
     &      CONFAC_QSWR,       !CONVFACTOR FOR SW RADIATION 
     &      CONFAC_QSAL,       !CONVFACTOR FOR PRECIPITATION-EVAPORATION 
     &      CONFAC_T0B,        !CONVFACTOR FOR SEA SURFACE TEM
     &      CONFAC_SICE,       !CONVFACTOR FOR SEA ICE MASK 
     &      CONFAC_TAUX,       !CONVFACTOR FOR ZONAL WIND STRESS 
     &      CONFAC_TAUY        !CONVFACTOR FOR MERIDIONAL WIND STRESS 

C CONVERSION FACTORS FOR OCEAN MODEL ONLY
      PARAMETER(CONFAC_QBAL=WPM2_GCMPS,          ! [W/M**2]=>[GRAD*CM/S]   HEAT BALANCE
     &          CONFAC_QSWR=WPM2_GCMPS,          ! [W/M**2]=>[GRAD*CM/S]   SW RADIATION
     &          CONFAC_QSAL=100.0,           ! [M/S]=>[CM/S]!            PRECIPITATION-EVAPORATION
     &          CONFAC_T0B =1.0,             ! [GRAD C]=>[GRAD C]        SEA SURFACE TEM
     &          CONFAC_SICE=0.01,            ! [%]=>[RATIO]              SEA ICE MASK
     &          CONFAC_TAUX=10.0/RH0,            ! [Pa]=>[CM**2/S**2]         ZONAL WIND STRESS    
     &          CONFAC_TAUY=10.0/RH0)            ! [Pa]=>[CM**2/S**2]         MERIDIONAL WIND STRESS


      REAL  CONFAC_LW,       !CONVFACTOR FOR DW-LW-RADIATION
     &      CONFAC_SW,       !CONVFACTOR FOR DW-SW-RADIATION
     &      CONFAC_PR,       !CONVFACTOR FOR PRECIPITATION
     &      CONFAC_TA,       !CONVFACTOR FOR ATMOSPHERE TEMPERATURE
     &      CONFAC_QA,       !CONVFACTOR FOR HUMIDITY
     &      CONFAC_UWND,     !CONVFACTOR FOR ZONAL WIND SPEED
     &      CONFAC_VWND,     !CONVFACTOR FOR MERIDIONAL WIND SPEED
     &      CONFAC_SLPR      !CONVFACTOR FOR PRESSURE

C  CONVERSION FACTORS FOR OCEAN MODEL WITH ICE BLOCK
      PARAMETER(CONFAC_LW=1000.0,      ! [W/M**2]=>[ERG/CM**2/S]      DW-LW-RADIATION
     &          CONFAC_SW=1000.0,      ! [W/M**2]=>[ERG/CM**2/S]      DW-SW-RADIATION
c    &          CONFAC_PR=MMPDAY_CMPS, ! [MM/DAY]=>[CM/S]             PRECIPITATION(climate)
     &          CONFAC_PR=0.1,         ! [MM/S]=>[CM/S]               PRECIPITATION(real-year)
     &          CONFAC_TA=1.0,         ! [GRAD C]=>[GRAD C]      ATMOSPHERE TEMPERATURE
     &          CONFAC_QA=1.0,         ! [KG/KG]=>[G/G]      HUMIDITY
c !!!! Tempory correction factor of 1.5 is used for hourly to minutes wind variability correction !!! 
c     &          CONFAC_UWND=100.0*1.5, ! [M/S]=>[CM/S]      ZONAL WIND    
c     &          CONFAC_VWND=100.0*1.5, ! [M/S]=>[CM/S]      MERID WIND
c---------------------base values for wind convertion------------------
     &          CONFAC_UWND=100.0,     ! [M/S]=>[CM/S]      ZONAL WIND    
     &          CONFAC_VWND=100.0,     ! [M/S]=>[CM/S]      MERID WIND
     &          CONFAC_SLPR=1.0)       ! [Pa]=>[Pa]      PRESSURE

      REAL CONFAC_RUNOFF               !CONFACTOR FOR RIVER RUNOFF
!     PARAMETER(CONFAC_RUNOFF=100.0)   ![M/S]->[CM/S]        !omip data
      PARAMETER(CONFAC_RUNOFF=0.1)     ![MM/S]->[CM/S]        !core data

      INTEGER FTYPE_QBL,                ! TYPE OF ATMOSPHERIC FORCING,
     &        FTYPE_QSW,                ! 
     &        FTYPE_PME,                ! 0 - MONTHLY MEAN DATA
     &        FTYPE_TSA,                ! 
     &        FTYPE_SIC,                ! 1 - REAL-TIME DATA
     &        FTYPE_TXA,                ! 
     &        FTYPE_TYA,                ! 
     &        FTYPE_SLP,                ! 
     &        FTYPE_ROFF,               ! 
     &        FTYPE_LWAT,               ! 
     &        FTYPE_SWAT,               ! 
     &        FTYPE_PRAT,               !
     &        FTYPE_TAT,                ! 
     &        FTYPE_QAT,                ! 
     &        FTYPE_UAT,                !
     &        FTYPE_VAT                 !

      PARAMETER ( FTYPE_QBL = 0,   ! TYPE OF ATMOSPHERIC FORCING,
     &            FTYPE_QSW = 0,   ! 
     &            FTYPE_PME = 0,   ! 0 - MONTHLY MEAN DATA
     &            FTYPE_TSA = 0,   ! 
     &            FTYPE_SIC = 0,   ! 1 - REAL-TIME DATA
     &            FTYPE_TXA = 0,   ! 
     &            FTYPE_TYA = 0,   ! 
     &            FTYPE_SLP = 1,   ! 
!     &            FTYPE_ROFF= 1,   ! OMIP
     &            FTYPE_ROFF= 0,   !  CORE    
     &            FTYPE_LWAT= 1,   ! 
     &            FTYPE_SWAT= 1,   ! 
     &            FTYPE_PRAT= 0,   !
     &            FTYPE_TAT = 1,   ! 
     &            FTYPE_QAT = 1,   ! 
     &            FTYPE_UAT = 1,   !
     &            FTYPE_VAT = 1         
     &                           )        

      INTEGER TIME_RESOLUTION_QBL,          ! TIME RESOLUTION 
     &        TIME_RESOLUTION_QSW,          ! 
     &        TIME_RESOLUTION_PME,          ! OF ATMOSPHERIC FORCING
     &        TIME_RESOLUTION_TSA,          ! 
     &        TIME_RESOLUTION_SIC,          ! IN HOURS
     &        TIME_RESOLUTION_TXA,          ! 
     &        TIME_RESOLUTION_TYA,          ! (IMPORTANT FOR REAL-TIME DATA ONLY)
     &        TIME_RESOLUTION_SLP,          ! 
     &        TIME_RESOLUTION_ROFF,         ! 
     &        TIME_RESOLUTION_LWAT,         ! 
     &        TIME_RESOLUTION_SWAT,         ! 
     &        TIME_RESOLUTION_PRAT,         !
     &        TIME_RESOLUTION_TAT,          ! 
     &        TIME_RESOLUTION_QAT,          ! 
     &        TIME_RESOLUTION_UAT,          !
     &        TIME_RESOLUTION_VAT           !

      PARAMETER( TIME_RESOLUTION_QBL=720,    ! TIME RESOLUTION 
     &           TIME_RESOLUTION_QSW=720,    ! 
     &           TIME_RESOLUTION_PME=720,    ! OF ATMOSPHERIC FORCING
     &           TIME_RESOLUTION_TSA=720,    ! 
     &           TIME_RESOLUTION_SIC=720,    ! IN HOURS
     &           TIME_RESOLUTION_TXA=720,    ! 
     &           TIME_RESOLUTION_TYA=720,    ! (SIGNIFICANT FOR REAL-TIME DATA ONLY)
     &           TIME_RESOLUTION_SLP=6,      ! 
!     &           TIME_RESOLUTION_ROFF=24,    ! OMIP
     &           TIME_RESOLUTION_ROFF=720,   !  CORE    
     &           TIME_RESOLUTION_LWAT=24,    ! 
     &           TIME_RESOLUTION_SWAT=24,    ! 
     &           TIME_RESOLUTION_PRAT=720,   !
     &           TIME_RESOLUTION_TAT=6,      ! 
     &           TIME_RESOLUTION_QAT=6,      ! 
     &           TIME_RESOLUTION_UAT=6,      !
     &           TIME_RESOLUTION_VAT=6       !
     &                                   )


      INTEGER NTIME_QBL,                ! TIME DIMENSION OF ATMOSPHERIC DATA ARRAY,
     &        NTIME_QSW,                ! 
     &        NTIME_PME,                ! 4 - FOR MONTHLY MEAN DATA
     &        NTIME_TSA,                ! 
     &        NTIME_SIC,                ! 31*24/TIME_RESOLUTION+3 - FOR REAL-TIME DATA
     &        NTIME_TXA,                ! 
     &        NTIME_TYA,                ! (TO BE CALCULATED AUTOMATICALLY)
     &        NTIME_SLP,                ! 
     &        NTIME_ROFF,               ! 
     &        NTIME_LWAT,               ! 
     &        NTIME_SWAT,               ! 
     &        NTIME_PRAT,               !
     &        NTIME_TAT,                ! 
     &        NTIME_QAT,                ! 
     &        NTIME_UAT,                !
     &        NTIME_VAT                 !


      PARAMETER(NTIME_QBL= (31*24/TIME_RESOLUTION_QBL+3)*FTYPE_QBL
     &                                     +NOFMONTH*(1-FTYPE_QBL),   ! TIME DIMENSION OF ATMOSPHERIC DATA ARRAY,

     &          NTIME_QSW= (31*24/TIME_RESOLUTION_QSW+3)*FTYPE_QSW
     &                                     +NOFMONTH*(1-FTYPE_QSW),   ! 

     &          NTIME_PME= (31*24/TIME_RESOLUTION_PME+3)*FTYPE_PME
     &                                     +NOFMONTH*(1-FTYPE_PME),   ! 4 - FOR MONTHLY MEAN DATA

     &          NTIME_TSA= (31*24/TIME_RESOLUTION_TSA+3)*FTYPE_TSA
     &                                     +NOFMONTH*(1-FTYPE_TSA),   ! 

     &          NTIME_SIC= (31*24/TIME_RESOLUTION_SIC+3)*FTYPE_SIC
     &                                     +NOFMONTH*(1-FTYPE_SIC),   ! 31*24/TIME_RESOLUTION - FOR REAL-TIME DATA

     &          NTIME_TXA= (31*24/TIME_RESOLUTION_TXA+3)*FTYPE_TXA
     &                                     +NOFMONTH*(1-FTYPE_TXA),   ! 

     &          NTIME_TYA= (31*24/TIME_RESOLUTION_TYA+3)*FTYPE_TYA
     &                                     +NOFMONTH*(1-FTYPE_TYA),   !  (TO BE CALCULATED AUTOMATICALLY)

     &          NTIME_SLP= (31*24/TIME_RESOLUTION_SLP+3)*FTYPE_SLP
     &                                     +NOFMONTH*(1-FTYPE_SLP),   ! 

     &          NTIME_ROFF=(31*24/TIME_RESOLUTION_ROFF+3)*FTYPE_ROFF
     &                                     +NOFMONTH*(1-FTYPE_ROFF),  ! 

     &          NTIME_LWAT=(31*24/TIME_RESOLUTION_LWAT+3)*FTYPE_LWAT
     &                                     +NOFMONTH*(1-FTYPE_LWAT),  ! 

     &          NTIME_SWAT=(31*24/TIME_RESOLUTION_SWAT+3)*FTYPE_SWAT
     &                                     +NOFMONTH*(1-FTYPE_SWAT),  ! 

     &          NTIME_PRAT=(31*24/TIME_RESOLUTION_PRAT+3)*FTYPE_PRAT
     &                                     +NOFMONTH*(1-FTYPE_PRAT),  !

     &          NTIME_TAT =(31*24/TIME_RESOLUTION_TAT+3)*FTYPE_TAT
     &                                     +NOFMONTH*(1-FTYPE_TAT),   ! 

     &          NTIME_QAT =(31*24/TIME_RESOLUTION_QAT+3)*FTYPE_QAT
     &                                     +NOFMONTH*(1-FTYPE_QAT),   ! 

     &          NTIME_UAT =(31*24/TIME_RESOLUTION_UAT+3)*FTYPE_UAT
     &                                     +NOFMONTH*(1-FTYPE_UAT),   !

     &          NTIME_VAT =(31*24/TIME_RESOLUTION_VAT+3)*FTYPE_VAT
     &                                     +NOFMONTH*(1-FTYPE_VAT)    !

     &                               )

C PARAMETERS OF ATMOSPHERIC FORCING ON ATMOSPHERIC GRID:
      REAL QBL,  !HEAT BALANCE [W/M**2] 
     &     QSW,  !SW RADIATION BALANCE[W/M**2]
     &     PME,  !PRECIPITATION-EVAPORATION[M/S]
     &     TSA,  !SEA SURFACE TEM.[°C]
     &     SIC,  !SEA ICE MASK[%]  
     &     TXA,  !ZONAL WIND STRESS[Pa=N/M**2] 
     &     TYA,  !MERIDIONAL WIND STRESS[Pa=N/M**2]
     &     SLP,  !PRESSURE AT SEA SURFACE
     &     ROFF, !RIVER RUNOFF
     &     LWAT, !DW-LW-RAD[W/M**2]
     &     SWAT, !DW-SW-RAD[W/M**2]
     &     PRAT, !PRECIPIT[M/S]
     &      TAT, !TEMP OF ATMOSPHERE[°C]
     &      QAT, !HUMIDITY [G/KG]
     &      UAT, !U-WIND SPEED[M/S]
     &      VAT  !V-WIND SPEED[M/S]

      COMMON /ATMGRIDFORC/
     &      QBL(NXA,NYA,NTIME_QBL), !HEAT BALANCE [W/M**2] 
     &      QSW(NXA,NYA,NTIME_QSW), !SW RADIATION BALANCE[W/M**2] 
     &      PME(NXA,NYA,NTIME_PME), !PRECIPITATION-EVAPORATION[M/S] 
     &      TSA(NXA,NYA,NTIME_TSA), !SEA SURFACE TEM.[°C] 
     &      SIC(NXA,NYA,NTIME_SIC), !SEA ICE MASK[%] 
     &      TXA(NXA,NYA,NTIME_TXA), !ZONAL WIND STRESS[Pa=N/M**2]
     &      TYA(NXA,NYA,NTIME_TYA), !MERIDIONAL WIND STRESS[Pa=N/M**2] 
     &      SLP(NXA,NYA,NTIME_SLP), !PRESSURE AT SEA SURFACE
     &     ROFF(NXA,NYA,NTIME_ROFF),!RIVER RUNOFF
     &     LWAT(NXA,NYA,NTIME_LWAT),!DW-LW-RAD[W/M**2]
     &     SWAT(NXA,NYA,NTIME_SWAT),!DW-SW-RAD[W/M**2]
     &     PRAT(NXA,NYA,NTIME_PRAT),!PRECIPIT[M/S]
     &      TAT(NXA,NYA,NTIME_TAT), !TEMP OF ATMOSPHERE[°C]
     &      QAT(NXA,NYA,NTIME_QAT), !HUMIDITY [G/KG]
     &      UAT(NXA,NYA,NTIME_UAT), !U-WIND SPEED[M/S]
     &      VAT(NXA,NYA,NTIME_VAT)  !V-WIND SPEED[M/S]

      INTEGER IND_CLIM_QBL,                !  INDEX OF CLIMATISITY
     &        IND_CLIM_QSW,                ! 
     &        IND_CLIM_PME,                !  0 - DATA OF CLIMATIC YEAR
     &        IND_CLIM_TSA,                ! 
     &        IND_CLIM_SIC,                !  1 - DATA OF REAL YEAR
     &        IND_CLIM_TXA,                !  (PATH TO DATA DEPENDS ON CURRENT YEAR)
     &        IND_CLIM_TYA,                ! 
     &        IND_CLIM_SLP,                ! 
     &        IND_CLIM_ROFF,               ! 
     &        IND_CLIM_LWAT,               ! 
     &        IND_CLIM_SWAT,               ! 
     &        IND_CLIM_PRAT,               !
     &        IND_CLIM_TAT,                ! 
     &        IND_CLIM_QAT,                ! 
     &        IND_CLIM_UAT,                !
     &        IND_CLIM_VAT                 !

      PARAMETER( IND_CLIM_QBL=0,                !  INDEX OF DATA TIME PERIOD
     &           IND_CLIM_QSW=0,                ! 
     &           IND_CLIM_PME=0,                !  0 - DATA OF CLIMATIC YEAR
     &           IND_CLIM_TSA=0,                ! 
     &           IND_CLIM_SIC=0,                !  1 - DATA OF REAL YEAR
     &           IND_CLIM_TXA=0,                !  (PATH TO DATA DEPENDS ON CURRENT YEAR)
     &           IND_CLIM_TYA=0,                ! 
     &           IND_CLIM_SLP=1,                ! 
!     &           IND_CLIM_ROFF=0,               ! OMIP
     &           IND_CLIM_ROFF=1,               !  CORE    
     &           IND_CLIM_LWAT=1,               ! 
     &           IND_CLIM_SWAT=1,               ! 
     &           IND_CLIM_PRAT=1,               !
     &           IND_CLIM_TAT=1,                ! 
     &           IND_CLIM_QAT=1,                ! 
     &           IND_CLIM_UAT=1,                !
     &           IND_CLIM_VAT=1                 !
     &                        )                 !


