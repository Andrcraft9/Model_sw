C  PARAMETERS FOR MATRIX INTERPOLATION FROM MODEL TO GEOGRAPHICAL GRID
      REAL     WGHT_MTRX_M2G,  !NONZERO MATRIX ELEMENTS FOR M2G INTERPOLATION
     &         WGHT_MTRX_M2G_NSKIP,
     &      ROTVEC_COEFF_M2G,  !COS(1) AND SIN(2) OF ANGLE BETWEEN LAT-LINES
     &             RMASK_OUT, RMASK_NSKIP,
     &             XT_G,YT_G,DX_G,DY_G

      INTEGER    I_INPUT_M2G,  !I-NUMBERS OF MATRIX ELEMENTS
     &           J_INPUT_M2G,  !J-NUMBERS OF MATRIX ELEMENTS
     &           I_INPUT_M2G_NSKIP,  !I-NUMBERS OF MATRIX ELEMENTS
     &           J_INPUT_M2G_NSKIP,  !J-NUMBERS OF MATRIX ELEMENTS
     &           IMASK_OUT

      COMMON /M2GINTRP/
     &      WGHT_MTRX_M2G(NX_GEO,NY_GEO,4),
     &      WGHT_MTRX_M2G_NSKIP(NX_GEO,NY_GEO,4),     
     &   ROTVEC_COEFF_M2G(NX_GEO,NY_GEO,4),
     &          RMASK_OUT(NX_GEO,NY_GEO),
     &          RMASK_NSKIP(NX_GEO,NY_GEO),
     &               XT_G(NX_GEO),
     &               YT_G(NY_GEO),
     &               DX_G(NX_GEO,NY_GEO),
     &               DY_G(NX_GEO,NY_GEO),
     &        I_INPUT_M2G(NX_GEO,NY_GEO,4),
     &        J_INPUT_M2G(NX_GEO,NY_GEO,4),
     &        I_INPUT_M2G_NSKIP(NX_GEO,NY_GEO,4),
     &        J_INPUT_M2G_NSKIP(NX_GEO,NY_GEO,4),
     &          IMASK_OUT(NX_GEO,NY_GEO)
