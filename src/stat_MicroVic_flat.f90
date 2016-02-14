module stat_MicroVic_flat

  
  !-- contains 
  !-     GlobalOrder, DensityAngle
  !-     Moment1Dx (with Div0x)
  !-     Moment2D  (with Div0)
  !-     DensityNgp1D, DensityNgp2D
  !-     DensityPic1D, DensityPic2D

  use toolkit                      ! for angle_vec_N, AngleModulo
  use input_output_MicroVic_flat   ! for PARAM_MicroVic_flat, FilePrintVector, FilePrint2D
  
contains

  !####################################################!
  !################   Global measure   ################!
  !####################################################!

  function GlobalOrder(theta)
    !- Compute the global flux 1/N Σ v_i
    !     theta :       a N vectors of angle θ_i
    !     GlobalOrder : the sum 1/N Σ v_i
    !
    implicit none
    Double Precision, dimension(:), intent(in)  :: theta
    Integer                                     :: N
    Double Precision, dimension(2)              :: meanMicroVic_flat
    Double Precision                            :: GlobalOrder
    !- init
    N = size(theta)
    meanMicroVic_flat = 1d0/N*(/ sum(cos(theta)) , sum(sin(theta)) /)
    !- the value
    GlobalOrder = sqrt( meanMicroVic_flat(1)**2 + meanMicroVic_flat(2)**2 )

  endfunction GlobalOrder


  Subroutine DensityAngle(theta,dtheta,iTime)
    !- Compute the global density of angle θ at time t_i
    !
    implicit none
    Double Precision, dimension(:), intent(in)  :: theta
    Double Precision, intent(in)                :: dtheta
    Integer, intent(in)                         :: iTime
    Double Precision, PARAMETER                 :: PI = 3.14159265358979323846
    Double Precision, dimension(:), Allocatable :: thetaModulo
    Double Precision, dimension(:), Allocatable :: densTheta
    Integer                                     :: i,n_th
    !- Init
    n_th = size(theta)
    Allocate(thetaModulo(n_th))
    Allocate(densTheta( floor(2*PI/dtheta)+1 ))
    !- Angle between -pi and pi
    Do i=1,n_th
       thetaModulo(i) = AngleModulo(theta(i))
    end Do
    !- The function (very simple!)
    Call DensityPic1D(-PI,PI, thetaModulo, densTheta)
    Call FilePrintVector(densTheta,"../data/densTheta_",.true.,iTime)
    !- Deallocate
    Deallocate(thetaModulo,densTheta)
  end Subroutine DensityAngle



  !#####################################################!
  !################     Moments in x     ###############!
  !#####################################################!

  Subroutine Moment1Dx(X,theta,P,iTime)
    !- Compute the macroscopic density and velocity
    !- of the distribution of particle in x.
    !
    implicit none
    Double Precision, dimension(:), intent(in)    :: X     !!!!! X est unidimensionnel ici
    Double Precision, dimension(:), intent(in)    :: theta
    TYPE(PARAM_MicroVic_flat), intent(in)         :: P
    Integer, intent(in)                           :: iTime
    Double Precision, dimension(:), allocatable   :: dens_x,f_u_x,u_x,f_v_x,v_x,angle_V
    Integer                                       :: nX
    !- Init
    nX = floor(P%Lx/P%dx)+1
    Allocate(dens_x(nX),f_u_x(nX),u_x(nX),f_v_x(nX),v_x(nX),angle_V(nX))
    !- Compute densities with PIC method
    Call DensityPic1D(0d0,P%Lx,&
         X, dens_x)
    Call DensityPic1D(0d0,P%Lx,&
         X, f_u_x, cos(theta))
    Call DensityPic1D(0d0,P%Lx,&
         X, f_v_x, sin(theta))
    !- Change the values at the boundaries if necessary
    if (P%boundCond==2 .or. P%boundCond==3) Then
       dens_x(1)   = (dens_x(1)+dens_x(nX))/2d0
       dens_x(nX) = dens_x(1)
       f_u_x(1)   = (f_u_x(1)+f_u_x(nX))/2d0
       f_u_x(nX) = f_u_x(1)
       f_v_x(1)   = (f_v_x(1)+f_v_x(nX))/2d0
       f_v_x(nX) = f_v_x(1)
    endif
    Call Div0x(f_u_x,dens_x,u_x)
    Call Div0x(f_v_x,dens_x,v_x)
    !- Print the result
    Call FilePrintVector(dens_x,"../data/dens1Dx_"//P%strSeed,.true.,iTime)
    Call FilePrintVector(u_x,"../data/u1Dx_"//P%strSeed,.true.,iTime)
    Call FilePrintVector(v_x,"../data/v1Dx_"//P%strSeed,.true.,iTime)
    !- Angle
    Call AngleVec_N(reshape( (/ u_x , v_x /), (/nX,2/) ),angle_V)
    Call FilePrintVector(angle_V,"../data/theta1Dx_"//P%strSeed,.true.,iTime)
    !- Deallocate
    DeAllocate(dens_x,f_u_x,u_x,f_v_x,v_x,angle_V)

  end Subroutine Moment1Dx

  Subroutine Div0x(f_dens,n_dens,u_dens)
    ! Divide one vector by another one with possible zero values
    !
    implicit none
    Double Precision, Dimension(:), intent(in)  :: f_dens,n_dens
    Double Precision, Dimension(:), intent(out) :: u_dens
    Integer                                     :: i,m
    !- Tnit
    m = size(f_dens)
    u_dens = 0d0
    !- Loop
    Do i=1,m
       if (n_dens(i)>0d0) Then
          u_dens(i) = f_dens(i)/n_dens(i)
       endif
    end Do
  end Subroutine Div0x



  !######################################################!
  !################     Moments in 2D     ###############!
  !######################################################!

  Subroutine Moment2D(X,theta,P,iTime)
    !- Compute the macroscopic density and velocity
    !- of the distribution of particle in (x,y).
    !
    implicit none
    Double Precision, dimension(:,:), intent(in)  :: X     !!!!! X is 2D
    Double Precision, dimension(:), intent(in)    :: theta
    TYPE(PARAM_MicroVic_flat), intent(in)         :: P
    Integer                                       :: iTime
    Double Precision, dimension(:,:), allocatable :: f_u,f_v
    real(kind=8), dimension(:,:), allocatable     :: rho,u,v
    Integer                                       :: nX,nY

    !- Init
    nX = floor(P%Lx/P%dxy + .5) + 1   ! for PIC
    nY = floor(P%Ly/P%dxy + .5) + 1   ! for PIC
    allocate(rho(nX,nY),f_u(nX,nY),u(nX,nY),f_v(nX,nY),v(nX,nY))
    
    !- Compute densities with PIC method
    Call DensityPic2D(0d0,P%Lx,0d0,P%Ly, &
         X(:,1), X(:,2), rho)
    Call DensityPic2D(0d0,P%Lx,0d0,P%Ly, &
         X(:,1), X(:,2), f_u, cos(theta))
    Call DensityPic2D(0d0,P%Lx,0d0,P%Ly, &
         X(:,1), X(:,2), f_v, sin(theta))
    ! Call DensityNGP2D(0d0,P%Lx,0d0,P%Ly,&
    !      X(:,1), X(:,2), rho)
    ! Call DensityNGP2D(0d0,P%Lx,0d0,P%Ly,&
    !      X(:,1), X(:,2), f_u, cos(theta))
    ! Call DensityNGP2D(0d0,P%Lx,0d0,P%Ly,&
    !      X(:,1), X(:,2), f_v, sin(theta))
    
    !- Change the values at the boundaries if necessary
    !  1:periodic, 2:reflexive, 3:periodic-reflexive, 4:reflexive-periodic, 5:ellipse
    if (P%boundCond==1 .or. P%boundCond==3) Then
       ! ρ(x=0,.)=ρ(x=Lx,.)
       rho(1,:)  = (rho(1,:)+rho(nX,:))/2d0
       rho(nX,:) = rho(1,:)
       f_u(1,:)  = (f_u(1,:)+f_u(nX,:))/2d0
       f_u(nX,:) = f_u(1,:)
       f_v(1,:)  = (f_v(1,:)+f_v(nX,:))/2d0
       f_v(nX,:) = f_v(1,:)
    endif
    if (P%boundCond==1 .or. P%boundCond==4) Then
       ! ρ(.,y=0)=ρ(.,y=Ly)
       rho(:,1)  = (rho(:,1)+rho(:,nY))/2d0
       rho(:,nY) = rho(:,1)
       f_u(:,1)  = (f_u(:,1)+f_u(:,nY))/2d0
       f_u(:,nY) = f_u(:,1)
       f_v(:,1)  = (f_v(:,1)+f_v(:,nY))/2d0
       f_v(:,nY) = f_v(:,1)
    endif
    Call Div0(f_u,rho,u)
    Call Div0(f_v,rho,v)
    !- Print the moments
    Call FilePrintArray2D(rho,"../data/rho2D_",.true.,iTime)
    Call FilePrintArray2D(u  ,"../data/u2D_"  ,.true.,iTime)
    Call FilePrintArray2D(v  ,"../data/v2D_"  ,.true.,iTime)
    
    Deallocate(rho,f_u,u,f_v,v)

  end Subroutine Moment2D


  Subroutine Div0(f_dens,n_dens,u_dens)
    ! Divide one matrix by another one with possible zero values
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in)  :: f_dens,n_dens
    Double Precision, Dimension(:,:), intent(out) :: u_dens
    Integer                                       :: i,j,m,n
    !- Init
    m = size(f_dens(:,1))
    n = size(f_dens(1,:))
    u_dens = 0d0
    !- Loop
    Do i=1,m
       Do j=1,n
          if (n_dens(i,j)>0d0) Then
             u_dens(i,j) = f_dens(i,j)/n_dens(i,j)
          endif
       end Do
    end Do

  end Subroutine Div0




  !###################################################!
  !################    density_NGP   #################!
  !###################################################!

  Subroutine DensityNgp1D(x_left,x_right,X_pos, M_dens,Attribut)
    ! Compute the density of particle with a Nearest-Grid-Point method:
    !   n_ep = n*phi   phi with square support
    ! in the domain
    !  x_left --- x_right
    !
    implicit none
    Double Precision, intent(in)                  :: x_left,x_right
    Double Precision, Dimension(:), intent(in)    :: X_pos
    Double Precision, Dimension(:), intent(out)   :: M_dens
    Double Precision, Dimension(:), intent(in), optional :: Attribut
    Double Precision                              :: dx, x
    Integer                                       :: i0, i, nX, N
    !- Init
    N = size(X_pos)
    nX = size(M_dens)
    dx  = (x_right - x_left)/(nX-1)
    M_dens = 0d0

    !------!
    ! Grid !     (1)    -   (2)    -       ....     -    (nX-1)    -    (nX) 
    !------!

    !---------------------------!
    !------     Loop      ------!
    !---------------------------!
    do i0=1,N
       x = X_pos(i0)
       !- Position on the grid of the i0 particle (west to be exact...)
       i = floor((x-x_left + dx/2 )/dx) + 1
       if ( 1<=i .and. i<=nX ) Then
          !- Add mass 1 on the grid
          if (present(Attribut)) Then
             M_dens(i) = M_dens(i) + Attribut(i0)
          else
             M_dens(i) = M_dens(i) + 1d0
          end if
       endif
    end do
    !- We normalize...
    M_dens = M_dens/dx

  end Subroutine DensityNgp1D


  Subroutine DensityNgp2D(x_left,x_right,y_down,y_up,X_pos,Y_pos, M_dens,Attribut)
    ! Number of particles per unit surface (histogram) in the domain.
    !    |--   y_up   --|
    !    |              |
    !  x_left        x_right
    !    |              |
    !    |--  y_down  --|
    !
    ! This 1st order method estimates the density (and u,v) at the *center* of the cells.
    ! 
    implicit none
    Double Precision, intent(in)                  :: x_left,x_right,y_down,y_up
    Double Precision, Dimension(:), intent(in)    :: X_pos, Y_pos
    Double Precision, Dimension(:,:), intent(out) :: M_dens
    Double Precision, Dimension(:), intent(in), optional :: Attribut
    Double Precision                              :: dx, dy
    Double Precision                              :: x, y
    Integer                                       :: i0, i, j, nX, nY, N
    !- Init
    N   = size(X_pos)
    nX  = size(M_dens(:,1))
    nY  = size(M_dens(1,:))
    dx  = (x_right - x_left)/nX
    dy  = (y_up    - y_down)/nY
    M_dens = 0d0

    !------!
    ! Grid !
    !------!

    !   (1,nY)   -  (2,nY)   -       ....     -   (nX-1,nY)   -   (nX,nY)
    !                                                                      
    !      .                                                               .   
    !      .                                                               .   
    !      .                                                               .   
    !      .                                                               .   
    !                                                                        
    !    (1,1)    -   (2,1)    -       ....     -    (nX-1,1)    -    (nX,1) 
    !

    !- warning
    if (size(Y_pos)/=N) Then
       print *,"In density_NGP, must have the same length for X_pos and Y_pos."
       stop
    end if


    !---------------------------!
    !------     Loop      ------!
    !---------------------------!

    do i0=1,N
       x = X_pos(i0)
       y = Y_pos(i0)
       !- Position on the grid of the i0 particle (south-west to be exact...)
       i = floor((x-x_left)/dx) + 1
       j = floor((y-y_down)/dy) + 1
       if ( 1<=i .and. i<=nX .and. 1<=j .and. j<=nY ) Then
          !- Add mass 1 on the grid
          if (present(Attribut)) Then
             M_dens(i,j) = M_dens(i,j) + Attribut(i0)
          else
             M_dens(i,j) = M_dens(i,j) + 1d0
          end if
       endif
    end do

    !- We normalize...
    M_dens = M_dens/dx/dy

  end Subroutine DensityNgp2D



  !###################################################!
  !################    density_PIC   #################!
  !###################################################!

  Subroutine DensityPic1D(x_left,x_right,X_pos, M_dens,Attribut)
    ! Compute the density of particle with second method :
    !   n_ep = n*phi   phi with square support
    ! in the domain
    !  x_left        x_right
    !
    implicit none
    Double Precision, intent(in)                  :: x_left,x_right
    Double Precision, Dimension(:), intent(in)    :: X_pos
    Double Precision, Dimension(:), intent(out)   :: M_dens
    Double Precision, Dimension(:), intent(in), optional :: Attribut
    Double Precision                              :: dx
    Double Precision                              :: x, x_p
    Integer                                       :: i0, i, nX, N
    !- Init
    N  = size(X_pos)
    nX = size(M_dens)
    dx = (x_right - x_left)/(nX-1)
    M_dens = 0d0

    !------!
    ! Grid !           (1)        ...              ...           (nX-1)      
    !------!     (1)    -   (2)    -       ....     -    (nX-1)    -    (nX) 

    !---------------------------!
    !------     Loop      ------!
    !---------------------------!
    do i0=1,N
       x = X_pos(i0)
       !- Position on the grid of the i0 particle (west to be exact...)
       i = floor((x-x_left)/dx) + 1
       if ( 1<=i .and. i<=(nX-1) ) Then
          !---     (i-1)dx < x < idx
          x_p = x - (x_left + (i-1)*dx)
          !- Add mass 1 on the grid
          if (present(Attribut)) Then
             M_dens(i)   = M_dens(i) + (dx - x_p)*Attribut(i0)
             M_dens(i+1) = M_dens(i+1) + x_p*Attribut(i0)
          else
             M_dens(i)   = M_dens(i) + (dx - x_p)
             M_dens(i+1) = M_dens(i+1) + x_p
          end if
       endif
    end do
    !- we normalize...
    M_dens = 1d0/N/dx**2*M_dens
    !- ...and adjust the values at the boundary
    M_dens(1)   = 2*M_dens(1)
    M_dens(nX) = 2*M_dens(nX)

  end Subroutine DensityPic1D


  Subroutine DensityPic2D(x_left,x_right,y_down,y_up,X_pos,Y_pos, M_dens,Attribut)
    ! Compute the number of particles using a 2nd oder method:
    !   n_ep = n*phi   phi with square support
    ! in the domain
    !    |--   y_up   --|
    !    |              |
    !  x_left        x_right
    !    |              |
    !    |--  y_down  --|
    !
    ! This 2nd order method estimates the density (and u,v) at the *grid point* of the cells.
    
    implicit none
    Double Precision, intent(in)                  :: x_left,x_right,y_down,y_up
    Double Precision, Dimension(:), intent(in)    :: X_pos, Y_pos
    Double Precision, Dimension(:,:), intent(out) :: M_dens
    Double Precision, Dimension(:), intent(in), optional :: Attribut
    Double Precision                              :: dx, dy
    Double Precision                              :: x, y, x_p, y_p
    Integer                                       :: i0, i, j, nX, nY, N
    !- Init
    N   = size(X_pos)
    nX  = size(M_dens(:,1))
    nY  = size(M_dens(1,:))
    dx  = (x_right - x_left)/(nX-1)
    dy  = (y_up    - y_down)/(nY-1)
    M_dens = 0d0

    !------!
    ! Grid !
    !------!
    !   (1,nY)   -  (2,nY)   -       ....     -   (nX-1,nY)   -   (nX,nY)
    !          (1,nY-1)      ...              ...           (nX-1,nY-1)   
    !      |                                                               |   
    !      |      .                                               .        |   
    !      |      .                                               .        |   
    !      |                                                                   
    !           (1,1)         ...              ...            (nX-1,1)      
    !    (1,1)    -   (2,1)    -       ....     -    (nX-1,1)    -    (nX,1) 
    !
    !- warning
    if (size(Y_pos)/=N) Then
       print *,"In density_PIC, must have the same length for X_pos and Y_pos."
       print *,"We have size(X_pos)=",size(X_pos)," and size(Y_pos)=",size(Y_pos)
       stop
    end if

    !---------------------------!
    !------     Loop      ------!
    !---------------------------!
    do i0=1,N
       x = X_pos(i0)  ;   y = Y_pos(i0)
       !- Position on the grid of the i0 particle (south-west to be exact...)
       i = floor((x-x_left)/dx) + 1
       j = floor((y-y_down)/dy) + 1
       if ( 1<=i .and. i<=(nX-1) .and. 1<=j .and. j<=(nY-1) ) Then
          !---     (i-1)dx < x < idx  ,   (j-1)dy < y < jdy
          x_p = x - (x_left + (i-1)*dx)
          y_p = y - (y_down + (j-1)*dy)
          !- Add mass 1 on the grid
          if (present(Attribut)) Then
             M_dens(i,j)     = M_dens(i,j) + (dx - x_p)*(dy - y_p)*Attribut(i0)
             M_dens(i+1,j+1) = M_dens(i+1,j+1) + x_p*y_p*Attribut(i0)
             M_dens(i+1,j)   = M_dens(i+1,j) + x_p*(dy-y_p)*Attribut(i0)
             M_dens(i,j+1)   = M_dens(i,j+1) + (dx-x_p)*y_p*Attribut(i0)
          else
             M_dens(i,j)     = M_dens(i,j) + (dx - x_p)*(dy - y_p)
             M_dens(i+1,j+1) = M_dens(i+1,j+1) + x_p*y_p
             M_dens(i+1,j)   = M_dens(i+1,j) + x_p*(dy-y_p)
             M_dens(i,j+1)   = M_dens(i,j+1) + (dx-x_p)*y_p
          end if
       endif
    end do

    !---- We normalize...
    M_dens = M_dens/dx**2/dy**2
    !- ...and adjust the values at the boundary
    M_dens(:,1)  = 2*M_dens(:,1)
    M_dens(:,nY) = 2*M_dens(:,nY)
    M_dens(1,:)  = 2*M_dens(1,:)
    M_dens(nX,:) = 2*M_dens(nX,:)

  end Subroutine DensityPic2D



end module stat_MicroVic_flat
