C PARAMETERS OF EARTH GLOBE SPHERE.
      REAL      PI, RE, GRV, RH0, CPW, G, RN
      REAL      OMEGA, CPWRH0

      PARAMETER(PI = 3.1415926,
     &          RE = 6.37E+08,           !EARTH RADIUS[CM]
C    &          RN = 1.0,                !SCALE NUMBER (=1.0, IF GRID STEPS ARE IN CENTIMETERS)
     &          RN = RE,                 !SCALE NUMBER (= RE, IF GRID STEPS ARE IN RADIANS)
     &       OMEGA = 7.2921235E-05,      !EARTH ANGULAR VELOCITY[RAD/SEC]
     &         CPW = 4.1,       !HEAT CAPACITY OF WATER[J/GR/GRAD)]FOR 35%.
     &         RH0 = 1.025,     !REFERENCE DENSITY[GR/CM**3]
     &         GRV = 980.6,     !FREE FALL ACCELERATION[CM/SEC**2]
     &           G = GRV/RH0,   !FREE FALL ACCELERATION[CM/SEC**2]/RHO
     &      CPWRH0 = CPW*RH0)
