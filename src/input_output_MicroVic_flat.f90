module input_output_MicroVic_flat

  !-- contains :
  !--     TYPE PARAM_MicroVic_flat, PARAM_InitMicroVic_flat
  !--     Lecture, PrintInfo, TestParam
  !--     FilePrintVector, FilePrintArray2D, BarProgress

  use toolkit

  TYPE PARAM_MicroVic_flat
     Integer                            :: N
     Double Precision                   :: c, nu, d, rVision
     Double Precision                   :: Lx, Ly
     Integer                            :: boundCond
     Double Precision                   :: epBC
     Double Precision                   :: dt,Time
     Integer                            :: numericalScheme
     Logical                            :: isGrid
     Logical                            :: isInitRand
     Integer                            :: nbSeed
     Character(len=80)                  :: strSeed
     Logical                            :: isTrajectorySave
     Integer                            :: jumpPrint
     Double Precision                   :: dx, dxy, dtheta

     Integer                            :: nCaseX, nCaseY
  end TYPE PARAM_MicroVic_flat

  TYPE PARAM_init
     Logical                            :: readPreviousSimu
     Double Precision                   :: previousTimeStep
     Integer                            :: initCondX
     Double Precision                   :: xMean,yMean,xVar,yVar,propLeft
     Integer                            :: initCondTheta
     Double Precision                   :: thetaMean,thetaVar,thetaL,thetaR,temperature
  end TYPE PARAM_init


contains

  Subroutine Lecture(P,Pinit)
    !- Read the parameters of the simulation from
    !- the file PARAMETER_MicroVic_flat
    !
    implicit none
    TYPE(PARAM_MicroVic_flat), intent(out)  :: P
    TYPE(PARAM_init), intent(out)           :: Pinit
    !- temp
    character(15)                           :: temp,tpStr
    Double Precision, PARAMETER             :: PI = 3.14159265358979323846

    !----------------------------------------!
    !- 1) Digesting the parameters file...  -!
    !----------------------------------------!
    open(unit=15,file='PARAMETER_MicroVic_flat.txt',status='old')

    read(15,*)temp
    read(15,*)temp
    read(15,*)temp
    read(15,*)temp
    read(15,*)P%N
    read(15,*)temp
    read(15,*)P%c
    read(15,*)P%nu
    read(15,*)P%d
    read(15,*)P%rVision
    read(15,*)temp
    read(15,*)P%Lx
    read(15,*)P%Ly
    read(15,*)temp
    read(15,*)temp
    read(15,*)P%boundCond
    read(15,*)temp
    read(15,*)P%dt
    read(15,*)P%Time

    read(15,*)temp
    read(15,*)temp
    read(15,*)P%numericalScheme
    read(15,*)temp
    read(15,*)P%isGrid
    read(15,*)temp
    read(15,*)P%isInitRand
    read(15,*)P%nbSeed
    read(15,*)temp
    read(15,*)P%isTrajectorySave
    read(15,*)temp
    read(15,*)P%dx
    read(15,*)P%dxy
    read(15,*)P%dtheta
    read(15,*)temp
    read(15,*)P%jumpPrint

    close(unit=15)

    !- extra term
    P%nCaseX  = floor( P%Lx/P%rVision )
    P%nCaseY  = floor( P%Ly/P%rVision )
    if (P%isInitRand) then
       P%strSeed = ""
    else
       write(tpStr , *) P%nbSeed
       P%strSeed = "seed"//trim(adjustl(tpStr))//"_"
    end if

    !---------------------------!
    !- 2) Digesting the IC...  -!
    !---------------------------!
    open(unit=16,file='PARAMETER_init.txt',status='old')

    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%readPreviousSimu
    read(16,*)Pinit%previousTimeStep
    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%initCondX
    read(16,*)temp
    read(16,*)Pinit%xMean
    read(16,*)Pinit%yMean
    read(16,*)Pinit%xVar
    read(16,*)Pinit%yVar
    read(16,*)temp
    read(16,*)Pinit%propLeft
    read(16,*)temp
    read(16,*)temp
    read(16,*)temp
    read(16,*)Pinit%initCondTheta
    read(16,*)temp
    read(16,*)Pinit%thetaMean
    read(16,*)Pinit%thetaVar
    read(16,*)temp
    read(16,*)Pinit%thetaL
    read(16,*)Pinit%thetaR
    read(16,*)Pinit%temperature

    close(unit=16)

  end Subroutine Lecture


  Subroutine PrintInfo(P)
    !- Write the parameters in the terminal during the simulation
    !
    implicit none
    TYPE(PARAM_MicroVic_flat), intent(in)               :: P
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

    !-----  Information on the terminal  -----!
    !-----------------------------------------!
    print *,"******************************************************"
    print *,"***********  Parameters of the simulation  ***********"
    print *,"******************************************************"
    print *,"|--------  Parameters model  --------|"
    print "(A25,I8)",  " Number particles (N) = ",P%N
    print "(A25,f8.2)","      speed (c)        =  ",P%c
    print "(A25,f8.2)","      relax (nu)       =  ",P%nu
    print "(A25,f8.2)","      noise (d)        =  ",P%d
    print "(A25,f8.2)","      rVision          =  ",P%rVision
    print "(A25,f8.2)","     Nb neighbors      = ",PI*P%rVision**2/(P%Lx*P%Ly)*P%N
    !--  boundary shape
    print "(A25,f7.1)","         Lx            = ",P%Lx
    print "(A25,f7.1)","         Ly            = ",P%Ly
    !--  boundary condition
    select case (P%boundCond)
    case(1)
       print *," Boundary condtion    :  periodic"
    case(2)
       print *," Boundary condtion    :  reflexive"
    case(3)
       print *," Boundary condtion    :  periodic-reflexive"
    case(4)
       print *," Boundary condtion    :  reflexive-periodic"
    case(5)
       print *," Boundary condtion    :  ellipse"
    case default
       print *, "  Incorrect choice for the boundary condition."
       stop
    end select
    print "(A25,f9.3)","         dt            = ",P%dt    
    print "(A25,f7.1)","     Total time        = ",P%Time

    !-- Computational stuf
    print *,"|----  Parameters computation  ------|"

    select case (P%numericalScheme)
    case (0)
       print *," Numerical scheme     :    Euler forward"
    case (1)
       print *," Numerical scheme     :    Circle method"
    case (2)
       print *," Numerical scheme     :    Exact"
    end select
    print *," Use grid             :     ",P%isGrid
    print *," Init  Random         :     ",P%isInitRand
    print *," Trajectory save      :     ",P%isTrajectorySave
    if (P%dx>0) Then
       print *," Save macro 1D in x"
       print "(A25,f8.2)","   Meshgrid (dx) in x  =   ",P%dx
       print "(A25,f8.2)","   Nb particles in dx  =   ",P%N*P%dx/P%Lx
    end if
    if (P%dxy>0) Then
       print *," Save macro 2D"
       print "(A25,f8.2)","   Meshgrid (dx,dy)    =   ",P%dxy
       print "(A25,f8.2)","   Nb particles dxdy   =   ",P%N*P%dxy/P%Lx*P%dxy/P%Ly
    end if
    if (P%dtheta>0) Then
       print *," Save distribution angle"
       print "(A25,f8.2)","   Meshgrid (dTheta)   =   ",P%dTheta
    end if
    
  end subroutine PrintInfo


  Subroutine TestParameter(P)
    !- Check if the parameters are not crazy...
    TYPE(PARAM_MicroVic_flat), intent(in)      :: P
    !- Some tests
    If (P%N>1d6) Then
       print *,"Warning : the number of particles is too large for 'FilePrint' (maximum 1 million)"
    end If
    If ( P%isGrid .and. P%nCaseX<3 ) Then
       print *,"The number of cases in x is too small (nCaseX=",P%nCaseX,")"
    end If
    If ( P%isGrid .and. P%nCaseY<3 ) Then
       print *,"The number of cases in y is too small (nCaseY=",P%nCaseY,")"
    end If
  end subroutine TestParameter


  Subroutine FilePrintVector(X,fileName,isUnFormatted,nbIteration)
    !- Write the real vector X in a file called 'fileName'.
    !- The writing can be easier formatted (i.e. human readable)
    !- or unformatted (binary format).
    implicit none
    Double Precision, Dimension(:), intent(in) :: X
    Character (len=*), intent(in)              :: fileName
    Logical, intent(in)                        :: isUnFormatted
    Integer, intent(in), optional              :: nbIteration
    Character(9)                               :: Extension
    Character (len=80)                         :: fileName_ext
    !- Init: change the fileName
    if (present(nbIteration)) Then
       write(Extension,'(I9.9)') nbIteration
       fileName_ext = trim(fileName)//trim(Extension)
    else
       fileName_ext = trim(fileName)
    endif
    if (isUnFormatted) then
       fileName_ext = trim(fileName_ext)//'.udat'
    else
       fileName_ext = trim(fileName_ext)//'.dat'
    endif
    !- Open and write
    if (isUnFormatted) then
       ! unformatted format
       Open(Unit=10,file=fileName_ext,form='UNFORMATTED')
       write(10) X
    else
       ! formatted format
30     format(100000(ES 15.7))
       Open(Unit=10,file=fileName_ext)
       write(10,30) X
    end if
    close(10)
  end Subroutine FilePrintVector


  Subroutine FilePrintArray2D(density,fileName,isUnFormatted,nbIteration)
    !- Same as FilePrint for array
    implicit none
    Double Precision, Dimension(:,:), intent(in) :: density
    Character (len=*), intent(in)                :: fileName
    Logical, intent(in)                          :: isUnFormatted
    Integer, intent(in), optional                :: nbIteration
    Character(9)                                 :: Extension
    Character (len=80)                           :: fileName_ext
30  format(100000(ES 15.7))

    !- Init: change the fileName
    if (present(nbIteration)) Then
       write(Extension,'(I9.9)') nbIteration
       fileName_ext = trim(fileName)//trim(Extension)
    else
       fileName_ext = trim(fileName)
    endif
    if (isUnFormatted) then
       fileName_ext = trim(fileName_ext)//'.udat'
    else
       fileName_ext = trim(fileName_ext)//'.dat'
    endif
    !- Open and write
    if (isUnFormatted) then
       ! unformatted format
       Open(Unit=10,file=fileName_ext,form='UNFORMATTED')
       write(10) density
       close(10)
    else
       ! formatted format
       Open(Unit=10,file=fileName_ext)
       write(10,30) density
       close(10)
    end if
  end Subroutine FilePrintArray2D


  Subroutine BarProgress(i,imax)
    !- Inform the user of the progression of the simulation
    !
    implicit none
    integer          :: i,imax,k
    character(len=1) :: bar, back
    ! the subroutine
    !-init
    back = char(8)
    bar  = '='
    ! print the percentage and the bar
    write(*,'(256a1)', advance='no') (back, k =1,(30*i/imax)+9)
    write(*,'(2x,1i3,1a1,2x,1a1,256a1)', advance='no') &
         100*i/imax,'%','|', (bar, k =1,30*i/imax)
    if (i==imax) Then
       write(*,'(a)') '| done.'
    end if
  End subroutine BarProgress


end module input_output_MicroVic_flat
