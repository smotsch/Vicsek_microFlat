module boundary_SppFlat2D

  !-- contains :
  !       Wall

  use input_output_SppFlat2D          ! for PARAM_SppFlat2D

contains

  subroutine Wall(X,V,theta,P)
    !-   Check if the particles in X have reached the boundaries.
    !-    If this is the case, the particle makes a "rebund".
    !
    implicit none
    Double Precision, Dimension(:,:), intent(inout) :: X
    Double Precision, Dimension(:,:), intent(inout) :: V
    Double Precision, dimension(:), intent(inout)   :: theta
    TYPE(PARAM_SppFlat2D), intent(in)               :: P
    Double Precision, dimension(2)                  :: X_i0,V_i0
    Integer                                         :: i0
    Double Precision, PARAMETER                     :: PI = 3.14159265358979323846

    do i0=1,P%N
       !- Has the particule X_i0 collided with the wall?
       X_i0 = X(i0,:)
       V_i0 = V(i0,:)
       !--  Left wall  --!
       if (X_i0(1)<0) Then
          if (P%boundCond==1 .or. P%boundCond==3) Then
             !- we shift the particle
             X(i0,1)    = X_i0(1) + P%Lx
          else
             !- rebound
             X(i0,1)    = -X_i0(1)
             V(i0,1)    = -V_i0(1)
             theta(i0)  = atan2(V(i0,2),V(i0,1))
          endif
       endif
       !--  Right wall  --!
       if (X_i0(1)>P%Lx) Then
          if (P%boundCond==1 .or. P%boundCond==3) Then
             !- we shift
             X(i0,1)    = X_i0(1) - P%Lx
          else
             !- rebound
             X(i0,1)    = P%Lx - (X_i0(1) - P%Lx)
             V(i0,1)    = -V_i0(1)
             theta(i0)  = atan2(V(i0,2),V(i0,1))
          endif
       endif
       !--  Down wall  --!
       if (X_i0(2)<0) Then
          if (P%boundCond==1) Then
             !- we shift the particle
             X(i0,2)    = X_i0(2) + P%Ly
          else
             !- rebound
             X(i0,2)    = -X_i0(2)
             V(i0,2)    = -V_i0(2)
             theta(i0)  = atan2(V(i0,2),V(i0,1))
          endif
       endif
       !--  Up wall  --!
       if (X_i0(2)>P%Ly) Then
          if (P%boundCond==1) Then
             !- we shift
             X(i0,2)    = X_i0(2) - P%Ly
          else
             !- rebound
             X(i0,2)    = P%Ly - (X_i0(2) - P%Ly)
             V(i0,2)    = -V_i0(2)
             theta(i0)  = atan2(V(i0,2),V(i0,1))
          endif
       endif

    end do
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
  TYPE(PARAM_SppFlat2D), intent(in)               :: P
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

  
end module boundary_SppFlat2D