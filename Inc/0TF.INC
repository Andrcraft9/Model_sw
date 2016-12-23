C MELTING TEMPERATURE AS FANCTION OF SALINITY:
       REAL      S, TF
       EXTERNAL  TF
C  S -SALINITY, TF -FREEZING POINT
C-----------------------------------------------------------------------
C  Makshtas A.P. The heat budget of Arctic ice in the winter.
C  Publ. by Int. Glaciological Soc. Cambridge CB2 1ER. UK. 1991. 77 p.:
C     TF(S) = -0.054*MAX(S,0.0)
C  Mellor&Kantha 1989:
C     TF(S) = -0.0543*S +7.59e-4(grad m-1)
C  Increasing of pressure causes increasing of melting temp.
C  Millero (UNESCO) 1978:
C     TF(S) = Tmelt-0.0575*S+1.71e-3 *S**1.5 - 2.155*S**2
C-----------------------------------------------------------------------
