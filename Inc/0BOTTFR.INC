ccc
ccc   INClude file for bottom friction parametrization.
ccc   Here is defined BottomFriction array, some physics coeffs
ccc   and type_of_bott_friction parameter
ccc

c     Type_fric - flag:
c                    0 - no bottom fric
c                    1 - linear
c                    2 - nonlinear
      INTEGER, PARAMETER:: type_fric = 2

      REAL       BottomFriction ! on temperature mask
      COMMON  /BOTTOMFRICTION/BOTTOMFRICTION(NX,NY)
c     Cb_l  - coef of    linear bottom friction (typical value 5E-2 cm/s)
c     Cb_nl - coef of nonlinear bottom friction (typical value 1E-3 )
c     Ebottom - bottom turbulent kinetic Energy (typical value 0:5E+2cm2/s2)
      REAL, PARAMETER::    Cb_l    = 5e-2,
c    &                     Cb_nl   = 1e-3,     !OPA
     &                     Cb_nl   =2.5e-3,    !Popov S.K.
c    &                     Ebottom = 10.0      !FRAM
     &                     Ebottom = 25.0      !FRAM
c    &                     Ebottom = 0.0       !OPA
