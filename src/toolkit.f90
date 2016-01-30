module toolkit

  !-- contains :
  !--   . AngleVec, AngleVec_N, AngleVec_NN, AngleModulo, f_dens
  !--   . InitRandSeed, RandNorm, RandNorm_N, RandomExpCos, VonMisesRnd

contains

  !---------------------------------!
  !--   fonctions élémentaires    --!
  !---------------------------------!

  function AngleVec(Vec)
    !- Compute the angle of the vector Vec
    !-    Vec :      a 2D vector
    !-    AngleVec : a scalar (θ)
    !
    implicit none
    Double Precision, Dimension(2), intent(in) :: Vec
    Double Precision                           :: AngleVec
    Double Precision, PARAMETER                :: PI = 3.14159265358979323846
    !-- arctan
    if (Vec(1)==0) Then
       AngleVec = PI/2*sign(1d0,Vec(2))
    else
       AngleVec = atan( Vec(2)/Vec(1) );
    endif
    !-- on remet entre -pi et pi
    if (Vec(1)<0d0) Then
       AngleVec = AngleVec + PI*sign(1d0,Vec(2))
    endif
  endfunction AngleVec

  subroutine AngleVec_N(V,Angle_V)
    !- Compute the angles of N vectors
    !-    V :        a N×2 matrix
    !-    Angle_V :  a N vector (θ_1,...,θ_N)
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in) :: V
    Double Precision, Dimension(:), intent(out)  :: Angle_V
    Integer                                      :: i,N
    !- init
    N = size(Angle_V)
    !- loop
    Do i=1,N
       Angle_V(i) = AngleVec( (/ V(i,1) , V(i,2) /) )
    end Do
  end subroutine AngleVec_N

  subroutine AngleVec_NN(U,V,Angle_mat)
    !- Compute the angles of N^2 vectors
    !-    U,V :        N×N matrix
    !-    Angle_Mat :  N×N matrix (θ_ij)
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in)  :: U,V
    Double Precision, Dimension(:,:), intent(out) :: Angle_mat
    Integer                                       :: i,j,M,N
    !- init
    M = size(U(:,1))
    N = size(U(1,:))
    !- loop
    Do i=1,M
       Do j=1,N
          Angle_mat(i,j) = AngleVec( (/ u(i,j) , v(i,j) /) )
       end Do
    end Do
  end subroutine AngleVec_NN

  function AngleModulo(theta)
    !- Give the angle θ between -π and π
    !     theta       :  scalar (angle)
    !     AngleModulo :  scalar (θ in ]-π,π])
    !
    implicit none
    Double Precision, intent(in)   :: theta
    Double Precision, PARAMETER    :: PI = 3.14159265358979323846
    Double Precision               :: AngleModulo
    !-- the function
    AngleModulo = theta - floor( (theta+PI)/(2d0*PI) )*(2d0*PI)
  end function AngleModulo


  function f_dens(x,L,lambda)
    !- Modify the uniform distribution for a two values distribution
    !
    implicit none
    Double Precision, intent(in)     :: x,L,lambda
    Double Precision                 :: f_dens
    !-           x* (L/2 / lambda)                  if x<lambda
    !-   f(x) =
    !-           L/2 + (x-lambda)*(L/2 / 1-lambda)   if lambda<x
    if (x<lambda) Then
       f_dens = x*(L/2/lambda)
    else
       f_dens = L/2 + (x-lambda)*(L/2/(1-lambda))
    end if
  end function f_dens


  !--------------------------------!
  !--          Random            --!
  !--------------------------------!

  Subroutine InitRandomSeed(nbSeed)
    !- Initialise le random
    !
    implicit none
    Integer, intent(in), optional       :: nbSeed
    Integer                             :: i, n, clock
    Integer, Dimension(:), Allocatable  :: seed

    Call Random_seed(size = n)
    Allocate(seed(n))
    if (present(nbSeed)) then
       seed = nbSeed * (/ (i - 1, i = 1, n) /)
    else
       Call System_clock(COUNT=clock)
       seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    end if
    Call Random_seed(PUT = seed)
    Deallocate(seed)
  End Subroutine InitRandomSeed

  function RandNorm()
    !- Generate a Gaussian with zero mean and variance 1
    !     RandNorm : a scalar
    !
    implicit none
    Double Precision, dimension(2)     :: X_unif
    Double Precision                   :: RandNorm
    Double Precision, Parameter        :: PI = 3.14159265358979323846
    !- classic formula
    Call Random_number(X_unif)
    RandNorm = sqrt(-2*log(X_unif(1)))*cos(2*PI*X_unif(2));
  end function RandNorm

  Subroutine RandNorm_N(Noise)
    !- Generate a vector of Gaussian
    !     RandNorm_N : a N vector of Gaussian
    !
    implicit none
    Double Precision, dimension(:), intent(out)   :: Noise
    Integer                                       :: N
    Double Precision, dimension(:,:), Allocatable :: X_unif
    Double Precision, Parameter                   :: PI = 3.14159265358979323846
    !- Allocate
    N = size(Noise)
    Allocate(X_unif(N,2))
    !- the classic formula for the Gaussian
    Call Random_number(X_unif)
    Noise = sqrt(-2*log(X_unif(:,1)))*cos(2*PI*X_unif(:,2));
    !- Deallocate
    DeAllocate(X_unif)
  end Subroutine RandNorm_N

  Subroutine RandomExpCos(theta,d)
    !- Generate number according to the Von-Mises distribution
    !    d     : parameter of the distribution
    !    theta : the number generated
    !
    implicit none
    Double Precision, intent(out)   :: theta
    Double Precision, intent(in)    :: d
    Double Precision, PARAMETER     :: PI = 3.14159265358979323846
    Double Precision                :: M,U,Y
    Logical                         :: bool

    !- init
    M = exp(1d0/d)
    bool = .true.

    !- loop (Rejection sampling)
    Do while (bool)
       Call Random_number(Y)
       Call Random_number(U)
       Y = 2d0*PI*(Y-.5d0)
       !- test
       if (U*M<exp(cos(Y)/d)) Then
          bool = .false.
       endif
    end Do

    !- finally
    theta = Y
  end Subroutine RandomExpCos


  Subroutine VonMisesRnd(theta,d)
    ! Fisher algorithm to generate a Von-Mises distribution. The idea is to use an acceptance-rejection method
    ! with a function close to a Von Mises distribution and that can be generated easily. The best function
    ! is a wrapped Cauchy distribution (Cauchy law on the circle).

    ! adaptation from:
    !    http://www.google.com/codesearch#RqChUiMGv2Q/trunk/octave-forge/main/statistics/inst/vmrnd.m&q=%22von%20mises%22&type=cs
    !
    ! Remark : the implementation in Fortran is slightly more complicated than in Octave/Matlab. Applying a 'mask' to an array is
    !          more delicate. We have used the 'where' statement, but it does not allow to use the random_number routine, which
    !          is why we use the extra variables randUnif*.


    implicit none
    Double Precision, dimension(:), intent(out)   :: theta
    Double Precision, intent(in)                  :: d
    Integer                                       :: N
    Double Precision                              :: a,b,r
    Double Precision, dimension(:), allocatable   :: u1,u2,u3
    Double Precision, dimension(:), allocatable   :: randUnif1,randUnif2,randUnif3
    Double Precision, dimension(:), allocatable   :: z,f,c
    Logical, dimension(:), allocatable            :: notDone
    Double Precision, PARAMETER                   :: PI = 3.14159265358979323846

    if (d > 1e6) then
       ! D is large: sample uniformly on circle
       Call Random_number(theta)
       theta = 2*PI*theta - PI;
    else
       ! init
       N = size(theta)
       a = 1 + sqrt (1 + 4/d**2);
       b = (a - sqrt (2 * a))*d/2;
       r = (1 + b**2) / (2 * b);
       !- Allocate
       Allocate(u1(N),u2(N),u3(N))
       Allocate(randUnif1(N),randUnif2(N),randUnif3(N))
       Allocate(z(N),f(N),c(N),notDone(N))
       notDone = .true.

       ! algorithm    
       do while (any(notDone))
          Call Random_number(randUnif1)
          Call Random_number(randUnif2)
          Call Random_number(randUnif3)
          where (notDone)
             u1 = randUnif1
             u2 = randUnif2
             u3 = randUnif3
             z = cos (PI * u1)
             f = (1 + r * z) / (r + z)
             c = (r - f)/d
          endwhere
          ! check
          notDone = (u2  >= c*(2 - c)) .and. (log(c) - log (u2) + 1 - c < 0)
       end do
       ! finish
       theta = sign (1d0,u3-0.5) * acos (f)

       ! Deallocate
       Deallocate(u1,u2,u3,randUnif1,randUnif2,randUnif3,z,f,c,notDone)
    endif
  end Subroutine VonMisesRnd




end module toolkit
