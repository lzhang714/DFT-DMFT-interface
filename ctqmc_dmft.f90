!
! MK: modified, at least to account for multiple atoms, potentially also in other ways.
! Original header:
!
!-------------------------------------------------------------------------
! project : azalea
! program : ctqmc_dmft_selfer
!           ctqmc_dmft_conver
!           ctqmc_dmft_mixer
! source  : ctqmc_dmft.f90
! type    : subroutine
! author  : li huang (email:huangli712@yahoo.com.cn)
! history : 09/16/2009 by li huang
!           09/18/2009 by li huang
!           09/21/2009 by li huang
!           09/22/2009 by li huang
!           10/28/2009 by li huang
!           12/01/2009 by li huang
!           12/06/2009 by li huang
!           12/24/2009 by li huang
!           01/13/2010 by li huang
! purpose : self-consistent engine for dynamical mean field theory (DMFT)
!           simulation. it is only suitable for hybridization expansion
!           version continuous time quantum Monte Carlo (CTQMC) quantum
!           impurity solver plus bethe lattice model.
! input   :
! output  :
! status  : unstable
! comment :
!-------------------------------------------------------------------------

! ========================================================================================= 
!>>> the self-consistent engine for continuous time quantum Monte Carlo
! quantum impurity solver plus dynamical mean field theory simulation

  subroutine ctqmc_dmft_selfer(iter)

     use constants
     use control
     use context
     use plo_quantities

     implicit none

! external arguments
! current iteration number
     integer, intent(in) :: iter

! local variables
! loop index over flavors
     integer :: i

! loop index over frequencies
     integer :: k

! atom. index.
     integer :: iat

! status flag
     integer :: istat

! dummy hybridization function, in matsubara frequency axis, matrix form
     complex(dp), allocatable :: htmp(:,:,:,:)

! allocate memory
     allocate(htmp(mfreq,norbs,norbs,natom), stat=istat)
     if ( istat /= 0 ) then
         call ctqmc_print_error('ctqmc_dmft_selfer','can not allocate enough memory')
     endif

! initialize htmp
     htmp = czero

! calculate new hybridization function using self-consistent condition
    call hybfunc_and_levels_sigma(htmp)

    do iat=1,natom
        call find_symmetries(iat)
    enddo

! mixing new and old hybridization function: htmp and hybf
     call ctqmc_dmft_mixer(hybf, htmp)

! update original hybridization function
     hybf = htmp

    select case (solvctrl)

    case('QMC  ')
!
! calculate new bath weiss's function in QMC case

    do iat=1,natom

     do i=1,norbs
         do k=1,mfreq
             wssf(k,i,i,iat) = cmesh(k) - eimp(i,iat) - hybf(k,i,i,iat)
         enddo ! over k={1,mfreq} loop
     enddo ! over i={1,norbs} loop

     do k=1,mfreq
         call ctqmc_zmat_inv(norbs, wssf(k,:,:,iat))
     enddo ! over k={1,mfreq} loop

! write out the new bath weiss's function
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_wssf(rmesh, wssf,iat)
     endif

! write out the new hybridization function
     if ( myid == master ) then ! only master node can do it
         call ctqmc_dump_hybf(rmesh, hybf,iat)
     endif

     enddo ! atoms

! print necessary self-consistent simulation information
     if ( myid == master ) then ! only master node can do it
         write(mystd,'(2X,a,i3)') 'DMFT hybridization function is updated at iteration ', iter
         write(mystd,*)
     endif

    end select

! deallocate memory
     deallocate(htmp)

     return
  end subroutine ctqmc_dmft_selfer

! ========================================================================================= 
!>>> check the convergence of self-energy function

  subroutine ctqmc_dmft_conver(convergence)

     use constants
     use control
     use context

     implicit none

! external arguments
! convergence flag
     logical, intent(inout) :: convergence

! local variables
! loop index over orbitals
     integer  :: i
     integer :: iat

! dummy variables
     real(dp) :: diff
     real(dp) :: norm
     real(dp) :: seps

! calculate diff and norm
! why not using the whole matrix? since the off-diagonal elementes may be NaN!

    do iat=1,natom

     diff = zero
     do i=1,norbs
         diff = diff + abs( sum( sig2(:,i,i,iat) - sig1(:,i,i,iat) ) )
     enddo ! over i={1,norbs} loop

     norm = zero
     do i=1,norbs
         norm = norm + abs( sum( sig2(:,i,i,iat) + sig1(:,i,i,iat) ) )
     enddo ! over i={1,norbs} loop
     norm = norm / two

! calculate seps
     seps = (diff / norm) / real(mfreq * norbs)

    end do ! natom

! write convergence information to screen
     if ( myid == master ) then ! only master node can do it 
         write(*,*)'----------------------------------------------------------------'
         write(mystd,'(2X,a,E13.6)') 'DMFT self-energy function convergence=', seps
         write(*,*)'----------------------------------------------------------------'
         write(mystd,*) 
     endif

! judge convergence status
     convergence = (seps <= eps8)

    do iat=1,natom
! update sig1
     sig1(:,:,:,iat) = sig1(:,:,:,iat) * (one - alpha) + sig2(:,:,:,iat) * alpha
    enddo

     return
  end subroutine ctqmc_dmft_conver

! ========================================================================================= 
!>>> complex(dp) version, mixing two vectors using linear mixing algorithm 

  subroutine ctqmc_dmft_mixer(vec1, vec2)

     use constants
     use control
     use plo_quantities, only : natom
     implicit none

! external arguments
! older green/weiss/sigma function in input
     complex(dp), intent(inout) :: vec1(mfreq,norbs,norbs,natom)

! newer green/weiss/sigma function in input, newest in output
     complex(dp), intent(inout) :: vec2(mfreq,norbs,norbs,natom)

! linear mixing scheme
     vec2 = vec1 * (one - alpha) + vec2 * alpha

     return

  end subroutine ctqmc_dmft_mixer




