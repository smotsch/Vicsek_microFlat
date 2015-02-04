module boundary_MicroVic_flat

  !-- contains :
  !       Wall

  use input_output_MicroVic_flat          ! for PARAM_MicroVic_flat

contains

  subroutine Wall(X,V,theta,P)
    !-   Check if the particles in X have reached the boundaries.
    !-    If this is the case, the particle makes a "rebund".
    !
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: X
    Double Precision, Dimension(:,:), intent(inout) :: V
    Double Precision, dimension(:), intent(inout)   :: theta
    TYPE(PARAM_MicroVic_flat), intent(in)           :: P
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

    select case(P%boundCond)
    case(1)
       ! periodic BC
       X(:,1) = X(:,1) + P%Lx * MERGE( 1, 0, X(:,1)<0 )
       X(:,1) = X(:,1) - P%Lx * MERGE( 1, 0, X(:,1)>P%Lx)
       X(:,2) = X(:,2) + P%Ly * MERGE( 1, 0, X(:,2)<0 )
       X(:,2) = X(:,2) - P%Ly * MERGE( 1, 0, X(:,2)>P%Ly)
    case(2)
       ! reflexive BC
       where (X(:,1)<0)
          X(:,1) = -X(:,1)
          V(:,1) = -V(:,1)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
       where (X(:,1)>P%Lx)
          X(:,1) = P%Lx - (X(:,1) - P%Lx)
          V(:,1) = -V(:,1)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
       where (X(:,2)<0)
          X(:,2) = -X(:,2)
          V(:,2) = -V(:,2)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
       where (X(:,2)>P%Ly)
          X(:,2) = P%Ly - (X(:,2) - P%Ly)
          V(:,2) = -V(:,2)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
    case(3)
       ! periodic in x, reflexive in y
       X(:,1) = X(:,1) + P%Lx * MERGE( 1, 0, X(:,1)<0 )
       X(:,1) = X(:,1) - P%Lx * MERGE( 1, 0, X(:,1)>P%Lx)
       where (X(:,2)<0)
          X(:,2) = -X(:,2)
          V(:,2) = -V(:,2)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
       where (X(:,2)>P%Ly)
          X(:,2) = P%Ly - (X(:,2) - P%Ly)
          V(:,2) = -V(:,2)
          theta  = atan2(V(:,2),V(:,1))
       endwhere
     
    end select
    
  end subroutine Wall

subroutine Wall_Circle(X,V,V_old,theta,P)
  !-   Check if the particles in X have reached the boundaries.
  !-    If this is the case, the particle makes a "rebund".
  !
  implicit none
  Double Precision, Dimension(:,:), intent(inout) :: X
  Double Precision, Dimension(:,:), intent(inout) :: V
  Double Precision, Dimension(:,:), intent(in)    :: V_old
  Double Precision, dimension(:), intent(inout)   :: theta
  TYPE(PARAM_MicroVic_flat), intent(in)               :: P
  Double Precision, dimension(2)                  :: X_i0,V_i0,V_m,X_star,eta
  Double Precision                                :: a,b,r2_i0,aa,bb,cc,t,scalar_prod
  Integer                                         :: i0
  Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

  ! init
  a = P%Lx/2
  b = P%Ly/2
  
  do i0=1,P%N
     !- Has the particule X_i0 collided with the wall?
     X_i0 = X(i0,:)
     V_i0 = V(i0,:)
     V_m  = (V(i0,:) + V_old(i0,:))/2 ! V middle (old+new)/2
     r2_i0 = (X_i0(1)-a)**2/a**2 + (X_i0(2)-b)**2/b**2
     !--  circle  --!
     if (r2_i0>1) Then
        ! time of impact
        aa = V_m(1)**2/a**2 + V_m(2)**2/b**2
        bb = V_m(1)*(X_i0(1)-a)/a**2 + V_m(2)*(X_i0(2)-b)/b**2
        cc = (X_i0(1)-a)**2/a**2 + (X_i0(2)-b)**2/b**2 - 1
        t = (-bb + sqrt(bb**2-aa*cc))/aa
        ! point of impact
        X_star = X_i0 + t*V_m
        ! normal η
        eta = (/ 2*(X_i0(1)-a)/a**2 , 2*(X_i0(2)-b)/b**2 /)
        eta = eta/sqrt(eta(1)**2+eta(2)**2)
        ! new position
        scalar_prod = (X_i0(1)-X_star(1))*eta(1) + (X_i0(2)-X_star(2))*eta(2)
        X(i0,:) = X_star - 2*scalar_prod*eta
        ! new velocity
        scalar_prod = V_i0(1)*eta(1) + V_i0(2)*eta(2)
        V(i0,:) = V_i0 - 2*scalar_prod*eta
        theta(i0)  = atan2(V(i0,2),V(i0,1))
     endif

  end do
end subroutine Wall_Circle

  
end module boundary_MicroVic_flat
