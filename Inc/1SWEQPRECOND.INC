        
        integer nzmax , precond, use_direct_solver
c
c  good value of nzmax for world ocean 2.5x2 is 550000
c  typical number for fill-in ration is about 5-7 (depends of size of matrice)
c 
c   remark: if time step tau and friction R are too small stability of ilut can be lost. 
c           you can treat that problem by decreasing of the time step, using direct solver,
c           increasing of nzmax parameter
        parameter (nzmax=10000000,  !size of data storage for precondition
     &             precond = 1,     !precond = 0 no preconditioner
     &             use_direct_solver = 0)      !if 0 then use iterative solver
C    &             use_direct_solver = 1)      !if 1 then use direct solver