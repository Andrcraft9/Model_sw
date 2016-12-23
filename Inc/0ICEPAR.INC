*--------------------------------------------
*     PARAMETERS OF HEAT, MASS AND MOMENTUM
*         EXCHANGE OF WATER, ICE AND SNOW.
*
*            VERSION 04.01.2003.
C---------------------------------------------------------------
      REAL EW,ES,EI,                           !EMISSIVITIES
     &     AOPW,AIDRY,AIWET,ASDRY,ASWET,       !ALBEDOS
C    &     ALICED,ALICEW,ALOC,
     &     AI0,                                !PENETRATING RADIATION PARAMETERIZATION
     &     CKSNOW,CKICE,                       !HEAT TRANSFER COEFFICIENTS
     &     CPA,CPW8, QICE,QSNOW,QLI,QLW,       !HEAT CAPACITY AND HEAT OF FUSION
     &     SICE8,SIGMA,ROA,ROW,ROI,ROSDRY,
     &     ROSWET,TMELT,E0, CMELT,             !GENERAL CONSTANTS
     &     CDHW,CDHI,CDA,CDW,                  !TRANSFER COEFFICIENTS FOR BULK FORMULAS
     &     HREF,GAMMA,HSMAX, HIMIN,HSMIN,AIMIN,
     &     ASWP                                !PART OF SW-RAD,PENETRATED IN WATER

*     EMISSIVITY OF WATER, SNOW AND ICE.
      PARAMETER(EW= 0.97,
     &          ES= 0.99,
     &          EI= 0.97)

*     ALBEDO.

      PARAMETER(AOPW=0.10,  ! OPEN WATER
     &         AIDRY=0.60,  ! DRY ICE       0.60      !JAKOVLEV'S VARIANT
c    &         AIWET=0.58,  ! WET ICE       0.50
     &         AIWET=0.50,  ! WET ICE       0.50     
c    &         ASDRY=0.71,  ! DRY SNOW      0.80
     &         ASDRY=0.8,   ! DRY SNOW      0.80
     &         ASWET=0.7)   ! WET SNOW      0.70

c      PARAMETER(AOPW=0.10,  ! OPEN WATER
c     &         AIDRY=0.65,  ! DRY ICE       0.60      !JAKOVLEV'S VARIANT NEW
c     &         AIWET=0.55,  ! WET ICE       0.50
c     &         ASDRY=0.82,  ! DRY SNOW      0.80
c     &         ASWET=0.75)  ! WET SNOW      0.70


*     PENETRATING RADIATION PARAMETERIZATION.
      PARAMETER(AI0= 0.17)

*     HEAT TRANSFER COEFFICIENTS

      PARAMETER(CKSNOW=0.31E5,
     &          CKICE =2.04E5)

*     HEAT CAPACITY AND HEAT OF FUSION.

      PARAMETER(CPA= 1004.E4,     ! AIR SPEC. HEAT STORAGE (MASS) ERG/G
c    &         CPW8= 3990.E4,     ! WATER SPEC. HEAT STORAGE (MASS) ERG/G
     &         CPW8= 4190.E4,     ! WATER SPEC. HEAT STORAGE (MASS) ERG/G
     &         QICE= 2.679E9,     ! HEAT OF FUSION, ICE, ERG/CM3
     &        QSNOW= 1.097E9,     ! HEAT OF FUSION, SNOW, ERG/CM3
     &          QLI= 2.835E10,    ! LATENT HEAT OF SUBLIMATION, ERG/G
     &          QLW= 2.501E10)    ! LATENT HEAT OF EVAPORATION, ERG/G
C     &           E0= 0.   )!!!!!!!!!!!!!!!!!!!

*     GENERAL CONSTANTS.

      PARAMETER(SICE8 = 4. ,      ! ICE MEAN SALINITY, PROMILLE
     &          SIGMA = 5.67E-5,  ! STEFAN-BOLZMANN (CGS)
c    &            ROA = 1.27E-3,  ! AIR DENSITY, G/CM3
     &            ROA = 1.30E-3,  ! AIR DENSITY, G/CM3
c    &            ROW = 1.026  ,  ! WATER DENSITY, G/CM3
     &            ROW = 1.025  ,  ! WATER DENSITY, G/CM3     
     &            ROI = 0.900  ,  ! ICE DENSITY, G/CM3
     &         ROSDRY = 0.32   ,  ! DRY SNOW DENSITY, G/CM3
     &         ROSWET = 0.45   ,  ! WET SNOW DENSITY, G/CM3
     &          TMELT = 0.      ) ! ICE & SNOW MELTING POINT, GRAD C


*     TRANSFER COEFFICIENTS FOR BULK FORMULAS (HEAT AND MOMENTUM).

      PARAMETER( CDHW= 1.20E-3,    ! AIR/WATER HEAT
c    &           CDHI= 3.0E-3 ,    ! AIR/ICE   HEAT
     &           CDHI= 1.5E-3 ,    ! AIR/ICE   HEAT
     &           CDA = 2.75E-3,    ! AIR/ICE   MOMENTUM
     &           CDW = 5.5E-3 )    ! ICE/WATER MOMENTUM


*     HEAT EXCHANGE COEFFICIENT FOR WINDLESS CONDITIONS.
      PARAMETER(E0=100.)  !  ROCT = RO_*CDH_*(E0+WIND).

*     ICE COMPACTNESS CHANGE DUE TO SIDE MELTING:
      PARAMETER(CMELT=0.5)        ! HIBLER, 1979
*     PARAMETER(CMELT=0.7)        ! MELLOR AND KANTHA, 1989


*     RATE OF LEADS CLOSING, HIBLER'S 1979 H*, H=50CM,
*     OR OBERHUBER'S HI - INSTANTENEOUS ICE, H= 1.
      PARAMETER(HREF=50.)

*     SNOW OLDERING AND CONVERSION TO ICE, ONLY IF HS > HSMAX,
*     OBERHUBER: 1E-7SEC-1.
      PARAMETER(GAMMA =1.E-7,
     &          HSMAX = 20.)

*     CUTOFFS FOR ICE AND SNOW.

      PARAMETER(HIMIN= 1.E-6,    ! MINIMAL ICE MASS (H*A)
     &          HSMIN= 1.E-8,    ! MINIMAL SNOW MASS
     &          AIMIN= 1.E-3)    ! MINIMAL COMPACTNESS

      PARAMETER(ASWP=0.4)






