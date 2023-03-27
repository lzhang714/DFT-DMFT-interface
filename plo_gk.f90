!
! Calculates the local Green function
!   G(z) = Sum_k ( z - H(k) - Sigma(z) )^-1
!
! Long Zhang created, most recently modified 2017
!
!
  subroutine plo_gk()

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only :  rmesh,cmesh
   use ctqmc_smat

    implicit none

    integer :: i
    integer :: ib
    integer :: is
    integer :: ikp
    integer :: iom
    integer :: lim1
    integer :: lim2
    integer :: Iwpart
    integer :: Iwrest
    integer :: iat

    real(dp) :: q0

    complex(dp) :: z
    complex(dp) :: g0(lda_nband,lda_nband)
    complex(dp) :: g(lda_nband,lda_nband)
    complex(dp) :: g_loc2(nband,nband,natom)
    complex(dp) :: sigblochmat(lda_nband,lda_nband)
    complex(dp) :: g0_loc(nband,nband,mfreq,nspin,natom)
    complex(dp) :: g_loc(nband,nband,mfreq,nspin,natom)
    complex(dp), allocatable :: gb_loc(:,:,:)

    allocate(gb_loc(lda_nband,mfreq,nspin))

    gb_loc(:,:,:) = czero
    g0_loc(:,:,:,:,:) = czero
    g_loc(:,:,:,:,:) = czero

    q0=1.d0/nkp
!
!  Divide frequencies over processors
!
    Iwpart = (mfreq) / nprocs
    Iwrest = mod(mfreq,nprocs)
!
!   Get Sigma
!
    call plo_rewrite_array( sig2,sigma,'solver2plo' )
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
    do is=1,nspin

      do iom=lim1,lim2

!      z = dcmplx(efermi,rmesh(iom))
       z = cmesh(iom) + efermi

        do ikp=1,nkp

        g(:,:) = czero
        g0(:,:) = czero

        forall( ib = 1:lda_nband ) g0(ib,ib) = z - dcmplx(hbloch(ib,ikp,is),0.00_dp)
!
!    Transform self-energy to bloch basis (multiatom version).
!
        call plo_upfold( sigma(:,:,iom,is,:),locpro(:,:,:,ikp,is),sigblochmat(:,:) )
!
!    Construct local Greens function from Bloch representation
!
        g(:,:) =  g0(:,:) - sigblochmat(:,:)

        call ctqmc_zmat_inv( lda_nband, g0 )
        call ctqmc_zmat_inv( lda_nband, g )
!
!    non-interacting
!
    do iat=1,natom

!    Downfold
    call plo_downfold( g0,locpro(:,:,iat,ikp,is),g_loc2(:,:,iat) )
!    Local
    g0_loc(:,:,iom,is,iat) = g0_loc(:,:,iom,is,iat)+g_loc2(:,:,iat)

    enddo
!
!    interacting
!
!    Bloch (only store diagonal)
    forall(i=1:lda_nband) gb_loc(i,iom,is) = gb_loc(i,iom,is)+g(i,i)

    do iat=1,natom

!    Downfold
    call plo_downfold(g,locpro(:,:,iat,ikp,is),g_loc2(:,:,iat))
!    Local
    g_loc(:,:,iom,is,iat) = g_loc(:,:,iom,is,iat)+g_loc2(:,:,iat)

    enddo

      enddo !kpoints
!
!    Normalize
!
    g_loc(:,:,iom,is,:) = q0*g_loc(:,:,iom,is,:)

    gb_loc(:,iom,is) =  q0*gb_loc(:,iom,is)

    g0_loc(:,:,iom,is,:) =  q0*g0_loc(:,:,iom,is,:)

    enddo !iom

    enddo !spin
!
!  Wait for all to finish
!
    call mp_barrier()
!
!  Store the full functions 
!
    call mp_allreduce(g_loc,Gloc)

    call mp_allreduce(gb_loc,Gbloc)

    call mp_allreduce(g0_loc,G0loc)

    call mp_barrier()

    deallocate(gb_loc)

  end subroutine plo_gk

!
!  Density matrix for charge selfconsistency
!
  subroutine plo_dens()

   use constants
   use plo_quantities
   use control
   use mmpi
   use ctqmc_umat, only :  rmesh
   use ctqmc_smat

    implicit none

    integer :: i
    integer :: j
    integer :: ib
    integer :: is
    integer :: ikp
    integer :: iom
    integer :: lim1
    integer :: lim2
    integer :: Iwpart
    integer :: Iwrest

    complex(dp) :: z
    complex(dp) :: z_KS
    complex(dp) :: g0(lda_nband,lda_nband)
    complex(dp) :: g(lda_nband,lda_nband)
    complex(dp) :: g_KS(lda_nband,lda_nband)
    complex(dp) :: sigblochmat(lda_nband,lda_nband)
    complex(dp) :: densbuff(lda_nband,lda_nband)
    complex(dp) :: densbuff2(lda_nband,lda_nband)
    complex(dp) :: dmft_dens_part(lda_nband,lda_nband)
!
!  Divide frequencies over processors
!
    Iwpart = (mfreq) / nprocs
    Iwrest = mod(mfreq,nprocs)
!
!   Get Sigma
!
    call plo_rewrite_array( sig2,sigma,'solver2plo' )
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
!  The Loop, now over k first, since we want Delta N(k)
!
    is=1

!    do is=1,nspin

    do ikp=1,nkp

    dmft_dens_part(:,:) = czero

      do iom=lim1,lim2

      z = dcmplx(efermi,rmesh(iom))
      z_KS = dcmplx(efermi_dft,rmesh(iom))

        g(:,:) = czero
        g0(:,:) = czero
        g_KS(:,:) = czero

        forall( ib = 1:lda_nband ) g0(ib,ib) = z - dcmplx(hbloch(ib,ikp,is),0.00_dp)

        forall( ib = 1:lda_nband ) g_KS(ib,ib) = z_KS - dcmplx(hbloch(ib,ikp,is),0.00_dp)
!
!    Transform self-energy to bloch basis.
!
        call plo_upfold( sigma(:,:,iom,is,:),locpro(:,:,:,ikp,is),sigblochmat(:,:) )
!
!    Construct local Greens function from Bloch representation
!
        g(:,:) =  g0(:,:) - sigblochmat(:,:)

        call ctqmc_zmat_inv( lda_nband, g0 )
        call ctqmc_zmat_inv( lda_nband, g )
        call ctqmc_zmat_inv( lda_nband, g_KS )

    forall( ib = 1:lda_nband ) sigblochmat(ib,ib) = sigblochmat(ib,ib) - (efermi - efermi_dft)
!
! Delta N(k) (positive matsubara frequency part)
!
!
! G_KS * Sigma
!
    call zgemm('N','N',lda_nband,lda_nband,lda_nband,cone,g_KS,lda_nband,sigblochmat(1,1),lda_nband,czero,densbuff,lda_nband)

!
! (G_KS * Sigma) * G
!
    call zgemm('N','N',lda_nband,lda_nband,lda_nband,cone,densbuff,lda_nband,g,lda_nband,czero,densbuff2,lda_nband)

! Add it
    dmft_dens_part(:,:) = dmft_dens_part(:,:) + densbuff2(:,:)

!
! Delta N(k) (negative matsubara frequency part)
!
!
! G_KS * Sigma
!
    call zgemm('C','C',lda_nband,lda_nband,lda_nband,cone,g_KS,lda_nband,sigblochmat(1,1),lda_nband,czero,densbuff,lda_nband)

!
! (G_KS * Sigma) * G
!
    call zgemm('N','C',lda_nband,lda_nband,lda_nband,cone,densbuff,lda_nband,g,lda_nband,czero,densbuff2,lda_nband)

! Add it
    dmft_dens_part(:,:) = dmft_dens_part(:,:) + densbuff2(:,:)

    enddo !iom

    call mp_barrier()

! Collect all parts from all procesors

     call mp_reduce( dmft_dens_part(:,:),dmft_dens(:,:), master)

    if(myid == master)then ! only the master

! Normalize matsubara sum
     dmft_dens(:,:) = (1.00_dp/beta)*dmft_dens(:,:)

! Write it into the file
     write(24,*)''
     do j=1,lda_nband
         do i=1,lda_nband
         write(24,*)real(dmft_dens(i,j)),aimag(dmft_dens(i,j))
         end do
     end do

    end if ! Master

      enddo !kpoints

!    enddo !spin
!
!  Wait for all to finish
!
    call mp_barrier()

  end subroutine plo_dens
