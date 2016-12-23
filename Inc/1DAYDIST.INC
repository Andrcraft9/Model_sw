      REAL DAYS_IN_4YR(49) !DAY DISTRIBUTION IN 4-YEAR
      DATA DAYS_IN_4YR/
c     &    0.,  30.,  60.,  90., 120., 150.,      !1-ST HALF YEAR OF 1-ST YEAR
c     &  180., 210., 240., 270., 300., 330.,      !2-ND HALF YEAR OF 1-ST YEAR
c     &  360., 390., 420., 450., 480., 510.,      !1-ST HALF YEAR OF 2-ND YEAR
c     &  540., 570., 600., 630., 660., 690.,      !2-ND HALF YEAR OF 2-ND YEAR
c     &  720., 750., 780., 810., 840., 870.,      !1-ST HALF YEAR OF 3-RD YEAR
c     &  900., 930., 960., 990.,1020.,1050.,      !2-ND HALF YEAR OF 3-RD YEAR
c     & 1080.,1110.,1140.,1170.,1200.,1230.,      !1-ST HALF YEAR OF 4-TH YEAR
c     & 1260.,1290.,1320.,1350.,1380.,1410.,1440./!2-ND HALF YEAR OF 4-TH YEAR
C------------REAL DISTRIBUTION WITHOUT LEAP-YEAR------------------------------
     &    0.,  31.,  59.,  90., 120., 151.,      !1-ST HALF YEAR OF 1-ST YEAR
     &  181., 212., 243., 273., 304., 334.,      !2-ND HALF YEAR OF 1-ST YEAR
     &  365., 396., 424., 455., 485., 516.,      !1-ST HALF YEAR OF 2-ND YEAR
     &  546., 577., 608., 638., 669., 699.,      !2-ND HALF YEAR OF 2-ND YEAR
     &  730., 761., 789., 820., 850., 881.,      !1-ST HALF YEAR OF 3-RD YEAR
     &  911., 942,  973.,1003.,1034.,1064.,      !2-ND HALF YEAR OF 3-RD YEAR
     & 1095.,1126.,1154.,1185.,1215.,1246.,      !1-ST HALF YEAR OF 4-TH YEAR
     & 1276.,1307.,1338.,1368.,1399.,1429.,1460./!2-ND HALF YEAR OF 4-TH YEAR
c----------------------------------------------------------------------
!     &    0.,  31.,  59.,  90., 120., 151.,      !1-ST HALF YEAR OF 1-ST YEAR
!     &  181., 212., 243., 273., 304., 334.,      !2-ND HALF YEAR OF 1-ST YEAR
!     &  365., 396., 424., 455., 485., 516.,      !1-ST HALF YEAR OF 2-ND YEAR
!     &  546., 577., 608., 638., 669., 699.,      !2-ND HALF YEAR OF 2-ND YEAR
!     &  730., 761., 789., 820., 850., 881.,      !1-ST HALF YEAR OF 3-RD YEAR
!     &  911., 942,  973.,1003.,1034.,1064.,      !2-ND HALF YEAR OF 3-RD YEAR
!     & 1095.,1126.,1155.,1186.,1216.,1247.,      !1-ST HALF YEAR OF 4-TH YEAR
!     & 1277.,1308.,1339.,1369.,1400.,1430.,1461./!2-ND HALF YEAR OF 4-TH YEAR
