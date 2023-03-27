!
! Calculates the initial hybridization function and impurity levels from the dft input.
! 
! MPI included.
!
! Long Zhang created, most recently modified 2017
!
!
  subroutine hybfunc_and_levels(atom)

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only : eimp, rmesh
   use ctqmc_wmat, only : hybf

   implicit none

    integer, intent(in) :: atom

    integer :: i
    integer :: is
    integer :: iom

    complex(dp) :: z

    call plo_calc_g0_mpi()

    do is = 1, nspin

    do iom = 1,mfreq

    z = dcmplx(0.00_dp,rmesh(iom))
!
!   Compute Hybridization function
!
      call ctqmc_zmat_inv( nband, G0loc(:,:,iom,is,atom) )

      do i = 1,nband
         G0loc(i,i,iom,is,atom)=G0loc(i,i,iom,is,atom)-z
      end do

         G0loc(:,:,iom,is,atom)=-(G0loc(:,:,iom,is,atom))

    enddo      ! mfreq

    enddo      ! spins
!
! Write tails into array eimp
!
  do i=1,nband
    eimp(i,atom) = Real(G0loc(i,i,mfreq,1,atom))
    eimp(i+nband,atom) = Real(G0loc(i,i,mfreq,2,atom))
  enddo
!
! Write Delta with tails -> 0 into array hybf
!
  do iom=1,mfreq
    do i = 1,nband
    hybf(iom,i,i,atom) = dcmplx( (Real(G0loc(i,i,iom,1,atom)-G0loc(i,i,mfreq,1,atom))),Imag(G0loc(i,i,iom,1,atom)) )
    hybf(iom,i+nband,i+nband,atom) = dcmplx( (Real(G0loc(i,i,iom,2,atom)-G0loc(i,i,mfreq,2,atom))),Imag(G0loc(i,i,iom,2,atom)) )
    enddo ! nband
  enddo ! mfreq

  end subroutine hybfunc_and_levels
!
! Calculate G0
!
  subroutine plo_calc_g0()

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only : eimp, rmesh
   use ctqmc_wmat, only : hybf

   implicit none

    integer :: ib
    integer :: is
    integer :: ikp
    integer :: iom
    integer :: iat

    real(dp) :: q0

    complex(dp) :: z
    complex(dp) :: g0(lda_nband,lda_nband)
    complex(dp) :: Gloc2(nband,nband,natom)

    G0loc(:,:,:,:,:) = czero

    q0 = 1.d0/nkp

    do is = 1, nspin

    do iat = 1, natom

    do iom = 1, mfreq

       z = dcmplx(efermi,rmesh(iom))
!
!    Calculation of local Green function.
!
       do ikp = 1, nkp

          g0(:,:) = czero

          forall( ib = 1:lda_nband ) g0(ib,ib) = z - dcmplx(hbloch(ib,ikp,is),0.00_dp)
!
!    Construct local Green's function from Bloch representation
!
      call ctqmc_zmat_inv( lda_nband, g0 )
!
!    non-interacting only

        call plo_downfold(g0, locpro(:,:,iat,ikp,is), Gloc2(:,:,iat))
    
          G0loc(:,:,iom,is,iat) = G0loc(:,:,iom,is,iat) + Gloc2(:,:,iat)

       enddo   ! nkp

      G0loc(:,:,iom,is,iat) = q0*G0loc(:,:,iom,is,iat)
      
      enddo ! mfreq
      
      enddo ! natom
      
      enddo ! nspin

  end subroutine plo_calc_g0

!
! Calculate G0
!
  subroutine plo_calc_g0_mpi()

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only : eimp, rmesh
   use ctqmc_wmat, only : hybf

   implicit none

    integer :: ib
    integer :: is
    integer :: ikp
    integer :: iom
    integer :: iat
    integer :: Iwpart
    integer :: Iwrest
    integer :: lim1
    integer :: lim2

    real(dp) :: q0

    complex(dp) :: z
    complex(dp) :: g0(lda_nband,lda_nband)
    complex(dp) :: Gloc2(nband,nband,natom)
    complex(dp) :: g0_loc(nband,nband,mfreq,nspin,natom)

    q0 = 1.d0/nkp

    g0_loc(:,:,:,:,:) = czero
!
!  Divide frequencies over processors
!
    Iwpart = (mfreq) / nprocs
    Iwrest = mod(mfreq,nprocs)
!
!  Setup frequency division for MPI
!
    if( myid == master )then
       lim1 = myid*Iwpart+1
       lim2 = (Iwpart*(myid+1))-1+Iwrest+1
    else
       lim1 = myid*Iwpart+Iwrest+1
       lim2 = (Iwpart*(myid+1))-1+Iwrest+1
    end if
!
!  The Loop
!
    do is = 1, nspin

    do iat = 1, natom

    do iom = lim1, lim2

       z = dcmplx(efermi,rmesh(iom))
!
!    Calculation of local Green function.
!
       do ikp = 1, nkp

          g0(:,:) = czero

          forall( ib = 1:lda_nband ) g0(ib,ib) = z - dcmplx(hbloch(ib,ikp,is),0.00_dp)
!
!    Construct local Green's function from Bloch representation
!
      call ctqmc_zmat_inv( lda_nband, g0 )
!
!    non-interacting only

        call plo_downfold(g0, locpro(:,:,iat,ikp,is), Gloc2(:,:,iat))

         g0_loc(:,:,iom,is,iat) = g0_loc(:,:,iom,is,iat) + Gloc2(:,:,iat)

       enddo   ! nkp

      g0_loc(:,:,iom,is,iat) = q0*g0_loc(:,:,iom,is,iat)

      enddo ! mfreq

      enddo ! natom

      enddo ! nspin

    call mp_barrier()

    call mp_allreduce(g0_loc,G0loc)

    call mp_barrier()

  end subroutine plo_calc_g0_mpi

!
! Calculates the initial hybridization function and impurity levels
! in general. With self energy.
! Since this tends to get lenghty, we have MPI.
! This calls plo_gk() and does some additional 
! Hybridization calculation in the end...
!
  subroutine hybfunc_and_levels_sigma(hybfunc)

   use constants
   use context, only : grnf
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only : eimp, rmesh
   use ctqmc_wmat
   use ctqmc_smat

   implicit none

   complex(dp), intent(out) :: hybfunc(mfreq,norbs,norbs,natom)

    integer :: i
    integer :: is
    integer :: iom
    integer :: iat

    complex(dp) :: Delta(nband,nband,mfreq,nspin,natom)

    Delta = czero

    call plo_gk()
!
! Compute Hybridization function
!
  do iat=1,natom

  do is=1,nspin
   do iom=1,mfreq

    Delta(:,:,iom,is,iat) = Gloc(:,:,iom,is,iat)

    call ctqmc_zmat_inv( nband, Delta(:,:,iom,is,iat) )

    Delta(:,:,iom,is,iat) = -Delta(:,:,iom,is,iat)-sigma(:,:,iom,is,iat)

    forall( i = 1:nband ) Delta(i,i,iom,is,iat) = Delta(i,i,iom,is,iat) + dcmplx(0.00_dp,rmesh(iom))

   enddo ! iom
  enddo ! spin

!
! Write tails into array eimp
!
  do i=1,nband
    eimp(i,iat) = Real(Delta(i,i,mfreq,1,iat))
    eimp(i+nband,iat) = Real(Delta(i,i,mfreq,2,iat))
  enddo
!
! Write Delta with tails -> 0 into array hybf
!
  do iom=1,mfreq
    do i = 1,nband
    hybfunc(iom,i,i,iat) = dcmplx( (Real(Delta(i,i,iom,1,iat)-Delta(i,i,mfreq,1,iat))),Imag(Delta(i,i,iom,1,iat)) )
    hybfunc(iom,i+nband,i+nband,iat) = dcmplx( (Real(Delta(i,i,iom,2,iat)-Delta(i,i,mfreq,2,iat))),Imag(Delta(i,i,iom,2,iat)) )
    enddo ! nband
  enddo ! mfreq
  
    enddo ! iatom

  end subroutine hybfunc_and_levels_sigma

!
! Find the symmetries of the system.
!
  subroutine find_symmetries(atom)

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only : eimp,symm

   implicit none

   integer :: i
   integer :: j
   
   integer, intent(in) :: atom

   do i=1,norbs
    do j=i,norbs

!   When we compare the same indices, do nothing
    if(i .eq. j)then
    continue

!   Compare and do nothing if they are 'degenerate'
    else if( (abs(eimp(i,atom)-eimp(j,atom))).lt.0.0001 )then
    continue

!   Compare and increase symm(j) by one if they are different
!   AND their indices are the same!
    else if( ( abs((eimp(i,atom)-eimp(j,atom))).gt.0.0001 ).and.(symm(i,atom) .eq. symm(j,atom)) )then
    symm(j,atom) = symm(j,atom) + 1

    end if

    enddo
   enddo

    if(myid==master) then
    write(*,*)'Impurity levels and Symmetries for Atom',atom,':'
    write(*,*)''
    write(*,"(30f13.8)")(eimp(i,atom), i=1,norbs)
    write(*,"(30i13.8)")(symm(i,atom), i=1,norbs)
    write(*,*)''
    end if
    
  end subroutine find_symmetries
