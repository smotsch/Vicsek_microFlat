module interaction_MicroVic_flat

  !-- contains :
  !       TestNeighbor
  !       VelocityAverageSlow
  !       VelocityAverageFast

  use toolkit                     ! for AngleVec
  use input_output_MicroVic_flat  ! for PARAM_MicroVic_flat
  use grid_interaction

contains

  Function TestNeighbor(X1,X2,R)
    !- Check whether the two particles X1 and X2 interact
    !
    implicit none
    Double precision, dimension(2), intent(in)  :: X1,X2
    Double precision, intent(in)                :: R
    Logical                                     :: TestNeighbor
    !- A boolean
    TestNeighbor = ( ((X2(1)-X1(1))**2 + (X2(2)-X1(2))**2) <= R*R )
  End Function TestNeighbor


  Subroutine VelocityAverageSlow(X,V,P,V_ave)
    !- Compute the average velocity V_ave with the "slow" technic of order O(N^2)
    implicit none
    Double Precision, Dimension(:,:), intent(in)  :: X, V
    TYPE(PARAM_MicroVic_flat), intent(in)         :: P
    Double Precision, Dimension(:,:), intent(out) :: V_ave
    Double Precision, Dimension(2)                :: X_i0, X_j
    Integer                                       :: i0,j
    !- init
    V_ave = V                   ! a particle sees itself

    !----------  Look around particle i0  ----------!
    Do i0=1,P%N
       X_i0 = X(i0,:)
       !- Find the neighbor of the particle i0
       Do j=i0+1,P%N
          !- The jth particle
          X_j  = X(j,:)
          !- Take care of the periodic boundary condiiton
          select case (P%boundCond)
          case(2)
             !- Periodic only in x
             if (abs(X_j(1)-X_i0(1))>P%Lx/2) Then
                X_j(1) = X_j(1) - P%Lx*sign(1d0,X_j(1)-X_i0(1))
             endif
          case(3)
             !- Periodic in x
             if (abs(X_j(1)-X_i0(1))>P%Lx/2) Then
                X_j(1) = X_j(1) - P%Lx*sign(1d0,X_j(1)-X_i0(1))
             endif
             !- Periodic in y
             if (abs(X_j(2)-X_i0(2))>P%Ly/2) Then
                X_j(2) = X_j(2) - P%Ly*sign(1d0,X_j(2)-X_i0(2))
             endif
          end select

          !- X_j is now the position of j in the reference of i0
          if (TestNeighbor(X_i0,X_j,P%rVision)) Then
             !- we take into account the speed of j
             V_ave(i0,:) = V_ave(i0,:) + V(j,:)
             V_ave(j,:)  = V_ave(j,:)  + V(i0,:)
          Endif
       end Do
    end Do

    ! Normalize V_ave
    Do i0=1,P%N
       V_ave(i0,:) = V_ave(i0,:) / sqrt( sum(V_ave(i0,:)**2) )
    end do
       
  end subroutine VelocityAverageSlow


  Subroutine VelocityAverageFast(X,V,P,posGrid,firstParticleGrid,verletListNext,V_ave)
    !- Computation of the average velocity V_ave with the "fast" method
    !- which requires the Verlet List(s)
    !
    implicit none
    Double Precision, Dimension(:,:), intent(in)  :: X, V
    TYPE(PARAM_MicroVic_flat), intent(in)         :: P
    Integer, Dimension(:), intent(in)             :: posGrid, firstParticleGrid
    Integer, Dimension(:), intent(in)             :: verletListNext
    Double Precision, Dimension(:,:), intent(out) :: V_ave
    Double Precision, Dimension(2)                :: X_i0, X_j, V_j, V_i0
    Integer                                       :: i0,j
    Integer                                       :: k, k_i0,i_grid,j_grid,i_g,j_g
    !- boundary condition
    Integer                                       :: ii,jj
    Double Precision                              :: dec_x, dec_y

    !- init
    V_ave = V                   ! a particle sees itself

    !- For each particle
    Do i0=1,P%N
       !-- Init
       X_i0 = X(i0,:)
       V_i0 = V(i0,:)
       !- Position on the grid
       k_i0   = posGrid(i0)
       i_grid = modulo(k_i0-1,P%nCaseX) + 1
       j_grid =      (k_i0-1)/P%nCaseX  + 1
       !- Compute the sum of the velocity of its neighbors
       Do i_g=i_grid-1,i_grid+1
          Do j_g=j_grid-1,j_grid+1
             ! ******* boundary conditions
             ii = i_g     ;  jj = j_g
             dec_x = 0d0  ;  dec_y = 0d0
             if (ii==0 .and. (P%boundCond==2 .or. P%boundCond==3)) Then
                ii = P%nCaseX ; dec_x = -P%Lx
             endif
             if (ii==(P%nCaseX+1) .and. (P%boundCond==2 .or. P%boundCond==3)) Then
                ii = 1      ; dec_x =  P%Lx
             endif
             if (jj==0 .and. P%boundCond==3) Then
                jj = P%nCaseY ; dec_y = -P%Ly
             endif
             if (jj==(P%nCaseY+1) .and. P%boundCond==3) Then
                jj = 1      ; dec_y =  P%Ly
             endif
             !- Summing
             !---------
             if (ii>=1 .and. ii<=P%nCaseX .and. jj>=1 .and. jj<=P%nCaseY) Then
                ! If non-periodic, we skip cases out of the box
                !- init
                k = (jj-1)*P%nCaseX + ii           ! case (ii,jj) in the list
                j = firstParticleGrid(k)           ! first particle in the case k
                !- while there is somebody in this box...
                Do While (j/=0)
                   if (j>i0) Then
                      !- The neighbor
                      X_j  = X(j,:) + (/dec_x , dec_y/)
                      V_j  = V(j,:)
                      !- Is j a neighbor ?
                      if (TestNeighbor(X_i0,X_j,P%rVision)) Then
                         !- We take into account the velocity of j
                         V_ave(i0,:) = V_ave(i0,:) + V_j
                         V_ave(j,:)  = V_ave(j,:)  + V_i0
                      End if
                   end if
                   ! Next neighbour!
                   j = verletListNext(j)
                end Do
             endif
             
          end Do
       end Do
    end do

    ! Normalize V_ave
    Do i0=1,P%N
       V_ave(i0,:) = V_ave(i0,:) / sqrt( sum(V_ave(i0,:)**2) )
    end do

  end subroutine VelocityAverageFast
  

end module interaction_MicroVic_flat
