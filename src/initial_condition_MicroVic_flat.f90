module initial_condition_MicroVic_flat

  !-- contains :
  !       InitCond

  
  use toolkit                      ! for random
  use input_output_MicroVic_flat   ! for PARAM_MicroVic_flat

contains

  Subroutine InitCond(X,theta,P,Pinit)
    !- initial condition
    implicit none
    Double Precision, Dimension(:,:), intent(out)   :: X
    Double Precision, Dimension(:), intent(out)     :: theta
    TYPE(PARAM_MicroVic_flat), intent(in)           :: P
    TYPE(PARAM_init), intent(in)                    :: Pinit
    Double Precision, Dimension(2)                  :: center
    Integer                                         :: i0
    Double Precision, Dimension(2)                  :: X_i0
    Double Precision                                :: r2_i0,a,b
    logical                                         :: isInsideDomain
    Double Precision, Dimension(1)                  :: thetaTemp,blockTemp
    Character(9)                                    :: Extension
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

    if (P%isInitRand) Then
       Call InitRandomSeed()
    endif
    ! init for ellipse IC
    a = P%Lx/2
    b = P%Ly/2

    !-- Compute the initial condition --!
    !-----------------------------------!

    if (Pinit%readPreviousSimu) then
       !-- Read the IC from a previous simulation --!
       !--------------------------------------------!
99     format(I9.9)
       write(Extension,99) Pinit%previousTimeStep
       ! The external file should be "data/particleX_*********.dat", etc. 
       Open(Unit=10,file=trim('data/particleX_')//trim(Extension)//trim('.udat'))
       Open(Unit=11,file=trim('data/particleY_')//trim(Extension)//trim('.udat'))
       Open(Unit=12,file=trim('data/particleTheta_')//trim(Extension)//trim('.udat'))
       read(10) X(:,1)
       read(11) X(:,2)
       read(12) theta(:)
       close(10)
       close(11)
       close(12)

    else
       !-- Generate the initial condition --!
       !------------------------------------!
       center = (/ P%Lx/2d0 , P%Ly/2d0 /)
       Do i0=1,P%N
          !----- Position -----!
          isInsideDomain = .false.
          Do while (.not. isInsideDomain)
             select case (Pinit%initCondX)
             case(1)
                !-  Uniform distribution
                Call Random_number(X_i0)
                X_i0 = (/ P%Lx*X_i0(1) , P%Ly*X_i0(2) /)
             case(2)
                !- Gaussian distribution
                X_i0(1) = Pinit%xMean + Pinit%xVar*RandNorm()
                X_i0(2) = Pinit%yMean + Pinit%yVar*RandNorm()
             case(3)
                !- Riemann problem
                !- multiply by f(x) with
                !-           x* (P%Lx/2 / proportion)                  if x<proportion
                !-   f(x) =
                !-           P%Lx/2 + (x-proportion)*(P%Lx/2 / 1-proportion)   if proportion<x
                Call Random_number(X_i0)
                X_i0(1) = f_dens(X_i0(1),P%Lx,Pinit%propLeft)
                X_i0(2) = P%Ly*X_i0(2)
             case(4)
                !-  4 blocks
                Call Random_number(X_i0)
                X_i0 = (/ P%Lx/2*X_i0(1) , P%Ly/2*X_i0(2) /)
                Call Random_number(blockTemp)
                if (blockTemp(1) > .5d0) then
                   X_i0(1) = X_i0(1)+P%Lx/2;
                   X_i0(2) = X_i0(2)+P%Ly/2;
                endif
             case default
                !- No good
                print *,"Wrong choice for initial position!!"
                print *,Pinit%initCondX," is not 1,2 or 3."
                stop
             end select
             !- Test if the position of the particle is in the box
             !----------------------------------------------------
             if (P%boundCond==4) Then
                r2_i0 = (X_i0(1)-a)**2/a**2 + (X_i0(2)-b)**2/b**2
                if (r2_i0<1) Then
                   isInsideDomain = .true.
                   X(i0,:) = X_i0
                endif
             else
                if (0d0<=X_i0(1) .and. X_i0(1)<P%Lx .and. 0d0<=X_i0(2) .and. X_i0(2)<P%Ly) Then
                   isInsideDomain = .true.
                   X(i0,:) = X_i0
                endif
             endif
          end Do
          !----- Angle θ -----!
          select case (Pinit%initCondTheta)
          case(1)
             !-  Uniform distribution
             Call Random_number(theta(i0))
             theta(i0) = 2*PI*theta(i0)
          case(2)
             !- Gaussian distribution
             theta(i0) = Pinit%thetaMean + Pinit%thetaVar*RandNorm()
          case(3)
             !- Riemann problem
             Call VonMisesRnd(thetaTemp,Pinit%temperature)
             if (X(i0,1)<P%Lx/2) Then
                theta(i0) = Pinit%thetaL + thetaTemp(1)
             else
                theta(i0) = Pinit%thetaR + thetaTemp(1)
             endif
          case default
             !- No good
             print *,"Wrong choice for initial theta !!"
             print *,Pinit%initCondTheta," is not 1,2 or 3."
             stop
          end select
       end Do                   !- end loop particle

    end if

  end subroutine InitCond

end module initial_condition_MicroVic_flat
