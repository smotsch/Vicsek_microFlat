Program MicroVic_flat

  use toolkit
  use input_output_MicroVic_flat
  use boundary_MicroVic_flat
  use initial_condition_MicroVic_flat
  use grid_interaction
  use interaction_MicroVic_flat
  use stat_MicroVic_flat
  
  implicit none
  
  !--------------  Declaration of variables  --------------!
  !--------------------------------------------------------!
  TYPE(PARAM_MicroVic_flat)                     :: P
  TYPE(PARAM_init)                              :: Pinit
  Integer                                       :: N
  Double Precision, Dimension(:,:), Allocatable :: X, V, V_ave, V_old, B
  Double Precision, Dimension(:), Allocatable   :: theta,angleB,angleV_ave,Noise
  Double Precision, Dimension(:), Allocatable   :: cos_dTheta0,sin_dTheta0,cos_dTheta1,sin_dTheta1,alpha,beta
  Double Precision, Dimension(:), Allocatable   :: globalOrderList
  Integer, Dimension(:), Allocatable            :: posGrid, firstParticleGrid
  Integer, Dimension(:), Allocatable            :: verletListPrev, verletListNext
  Integer                                       :: iStep,nSteps
  Character(80)                                 :: nameFile
  Double Precision                              :: sq_2d_dt,stdNoise
  Integer                                       :: i0, i0PosGrid, i0PosGrid_old
  Real                                          :: start, finish
  
  !--  Lecture of the parameters   --!
  !----------------------------------!
  Call Lecture(P,Pinit)
  N = P%N
  sq_2d_dt = sqrt(2*P%d*P%dt)
  ! Std noise:  √( d/ν (1-exp(-2ν∆t)) )
  !    interpolation between √(2d∆t) and √(d/ν)
  stdNoise = sqrt( P%d/P%nu*(1-exp(-2*P%nu*P%dt)) )
  !- Number of steps
  nSteps = floor(P%Time/P%dt + .5)
  !- warning
  Call TestParameter(P)
  !-- Information
  Call PrintInfo(P)
  
  !--      Initialisation       --!
  !-------------------------------!
  !- Allocate
  Allocate(X(N,2),V(N,2),V_ave(N,2),theta(N))
  Allocate(V_old(N,2),B(N,2),angleB(N),angleV_ave(N),Noise(N))
  Allocate(globalOrderList(nSteps+1))
  If ( P%numericalScheme==2 ) Then
     Allocate(cos_dTheta0(N),sin_dTheta0(N),cos_dTheta1(N),sin_dTheta1(N))
     Allocate(alpha(N),beta(N))
  end If
  If ( P%isGrid ) Then
     Allocate(posGrid(N))
     Allocate(firstParticleGrid(P%nCaseX*P%nCaseY))
     Allocate(verletListPrev(N),verletListNext(N))
  End If
  !- Initial condition
  Call InitCond(X,theta, P,Pinit)
  V = reshape( (/ cos(theta), sin(theta) /) , (/N,2/) )
  !- Verlet list
  if ( P%isGrid ) Then
     Call InitVerletList(X,P,posGrid,firstParticleGrid,&
          verletListPrev,verletListNext)
  endif
  !- Let's count the time it takes
  Call Cpu_time(start)
  if (P%isTrajectorySave) then
     !- Write the intial condition
     if (P%isFormatVtk) then
        Call FilePrintParticle_vtk(X,V,"../data/particle_"//P%strSeed,0)
     else
        Call FilePrintVector(X(:,1),"../data/particleX_",.true.,0)
        Call FilePrintVector(X(:,2),"../data/particleY_",.true.,0)
        Call FilePrintVector(theta,"../data/particleTheta_",.true.,0)
     endif
  endif
  !- The moments
  if (P%dx/=0d0)     Then ; Call Moment1Dx(X(:,1),theta,P,0);    end if
  if (P%dxy/=0d0)    Then ; Call Moment2D(X,theta,P,0);          end if
  if (P%dtheta/=0d0) Then ; Call DensityAngle(theta,P%dtheta,0); end if


     
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!
  !-----------------------   The big loop   ----------------------!
  !---------------------------------------------------------------!
  !---------------------------------------------------------------!

  Do iStep=1,nSteps

     !- 1) Compute the average velocity around each particle -!
     !--------------------------------------------------------!
     If ( P%isGrid ) Then
        Call VelocityAverageFast(X,V,P,&
             posGrid,firstParticleGrid,verletListNext, V_ave)
     else
        Call VelocityAverageSlow(X,V,P, V_ave)
     End If
     
     !- 2) Compute the new velocity and position -!
     !--------------------------------------------!
     ! init
     Call RandNorm_N(Noise)
     V_old = V
     ! update V
     select case (P%numericalScheme)
     case (0)
        !-- Euler Forward method --!
        Call AngleVec_N(V_ave,angleV_ave)
        theta = theta + P%nu*P%dt*sin(angleV_ave - theta) + sq_2d_dt*Noise
        V     = reshape( (/ cos(theta), sin(theta) /) , (/N,2/) )
     case (1)
        !-- Technic Circle (Improved Euler method) --!
        B = V + P%nu*P%dt/2*( V_ave - V )
        Call AngleVec_N(B,angleB)
        theta = theta + 2*(angleB - theta) + sq_2d_dt*Noise
        V     = reshape( (/ cos(theta), sin(theta) /) , (/N,2/) )
     case(2)
        !-- Explicit solution
        ! explicit solution of ∆θ'= -μ sin ∆θ (here ∆θ_0 = θ_0 - θ_bar)
        !- init: ∆θ_0
        cos_dTheta0 = V(:,1)*V_ave(:,1) + V(:,2)*V_ave(:,2) ! cos ∆θ = cos(θ-θ_bar) = ...
        sin_dTheta0 = V(:,2)*V_ave(:,1) - V(:,1)*V_ave(:,2) ! sin ∆θ = sin(θ-θ_bar) = ...
        !- ∆θ_1
        alpha = (1+cos_dTheta0)*exp(2*P%Nu*P%Dt)
        beta  = (1-cos_dTheta0)
        cos_dTheta1 = (alpha-beta)/(alpha+beta)
        sin_dTheta1 = sqrt(abs(1-cos_dTheta1**2))*sign(1d0,sin_dTheta0)
        !- V = R_∆θ(t)⋅V_ave
        V(:,1) = cos_dTheta1*V_ave(:,1) - sin_dTheta1*V_ave(:,2)
        V(:,2) = sin_dTheta1*V_ave(:,1) + cos_dTheta1*V_ave(:,2)
        !- add noise
        theta  = atan2(V(:,2),V(:,1)) + sq_2d_dt*Noise
        V(:,1) = cos(theta)
        V(:,2) = sin(theta)
     end select
     !- update X
     X = X + P%c*P%dt/2*(V_old + V)
     !- The particles may have smashed the wall...
     if (P%boundCond==5) Then
        Call Wall_Circle(X,V,V_old,theta,P)
     else
        Call Wall(X,V,theta,P)
     endif
    
     !- 3) Update the Verlet_list -!
     !-----------------------------!
     If (P%isGrid) Then
        do i0=1,N
           !-- init
           i0PosGrid     = CellNumber(X(i0,:),P)
           i0PosGrid_old = posGrid(i0)
           !- Test if i0 has changed its position on the grid
           if (i0PosGrid /= i0PosGrid_old) then
              !- We modify the Verlet_list...
              Call ListRemove(i0,i0PosGrid_old,firstParticleGrid,&
                   verletListPrev,verletListNext)
              Call ListAdd(i0,i0PosGrid,firstParticleGrid,&
                   verletListPrev,verletListNext)
              !- ...and we update posGrid
              posGrid(i0) = i0PosGrid
           end if
        end do
     end If

     !- 4) Write the data and statistical analysis -!
     !----------------------------------------------!
     if (modulo(iStep,P%jumpPrint)==0 .or. iStep==nSteps) Then
        if (P%isTrajectorySave) then
           !- print particles
           if (P%isFormatVtk) then
              Call FilePrintParticle_vtk(X,V,"../data/particle_"//P%strSeed,iStep)
           else
              Call FilePrintVector(X(:,1),"../data/particleX_",.true.,iStep)
              Call FilePrintVector(X(:,2),"../data/particleY_",.true.,iStep)
              Call FilePrintVector(theta,"../data/particleTheta_",.true.,iStep)
           endif
        end if
        if (P%dx/=0d0)     Then ; Call Moment1Dx(X(:,1),theta,P,iStep);    end if
        if (P%dxy/=0d0)    Then ; Call Moment2D(X,theta,P,iStep);          end if
        if (P%dtheta/=0d0) Then ; Call DensityAngle(theta,P%dtheta,iStep); end if
     end if
     !- Statistical analysis
     globalOrderList(iStep+1) = GlobalOrder(theta)     

     !- progress...
     Call BarProgress(iStep,nSteps)
  End do
  
  !---------------------------------------------------------------!
  !---------------------   End loop in time   --------------------!
  !---------------------------------------------------------------!
  
  !- To conclude...
  Call FilePrintVector(globalOrderList,"../data/globalOrder",.false.)
  !- Deallocate  
  Deallocate(X,V,theta,V_ave,V_old,B)
  Deallocate(angleB,Noise,globalOrderList)
  If ( P%numericalScheme==2 ) Then
     Deallocate(cos_dTheta0,sin_dTheta0,cos_dTheta1,sin_dTheta1,alpha,beta)
  end If
  If ( P%isGrid ) Then
     Deallocate(posGrid,firstParticleGrid,verletListPrev,verletListNext)
  End If
  Call Cpu_time(finish)
  print *,""
  print "(A24,f11.3)"," Time to compute (s)   = ",finish-start
  print *,"******************************************************"
  print *,""

  
End Program MicroVic_flat
