!
! Calculates the Double Counting energy
!
! L.Zhang: created, then most recently modified 2017
!

subroutine plo_dcc(iter)

  use control
  use plo_quantities
  use ctqmc_umat
  use constants
  use mmpi

  implicit none

  integer, intent (in ) :: iter

  integer :: i
  integer :: j
  integer :: is
  integer :: iat

  real(dp) :: usum, corr_dos, step_dc

  ntot(:,:) = 0.00_dp

!
!    We have FLL, AMF, TrG=TrG0 (MET) and a number as choices
!

!
! Now atom dependent. Free of charge.
!

   do iat=1, natom

!     Compute average interaction parameters from umatrix

    Uavg(iat)=0.0

    Javg(iat)=0.0

      do i=(norbs/2)+1,norbs
        do j=1,norbs/2
          Uavg(iat)=Uavg(iat)+uumat(i,j,iat)
        enddo
      enddo

      Uavg(iat)=Uavg(iat)/((norbs/2)**2)

      do i=1,norbs/2
        do j=1,norbs/2
          Javg(iat)=Javg(iat)+uumat(i,j,iat)
        enddo
      enddo

      Javg(iat)=Javg(iat)/(((norbs/2)**2)-(norbs/2))

      Javg(iat)=(Uavg(iat)-Javg(iat))

    select case (dc_type(iat))
!
!     Fully localized limit
!
      case(0)

        do is = 1, nspin
          do i=1, nband
          ntot(is,iat) = ntot(is,iat) + occ(i,i,is,iat)
          enddo
        end do

        if( nspin == 1 ) ntot(2,iat) = ntot(1,iat)

        do is = 1, nspin
          do i = 1, nband
            dcc(is,iat) = - Uavg(iat)*( sum(ntot(:,iat)) - 0.5 ) +  Javg(iat)*( ntot(is,iat)  - 0.5 )
          end do
        end do
!
!	Around mean field
!
      case(1)
        do is = 1, nspin
          do i=1, nband
          ntot(is,iat) = ntot(is,iat) + occ(i,i,is,iat)
          enddo
        end do
        
!	Sum up one line or column of density-density U matrix
	usum=0.0
	j=1
	do i = 1,norbs
	    usum = usum + uumat(i,j,iat)
	end do

        do is = 1, nspin
          do i = 1, nband
            dcc(is,iat) = - (sum(ntot(:,iat))/(norbs))*usum
          end do
        end do

!
!	Tr G=Tr G0
!
    case(2)

    if(iter==1)then 
	do i = 1, nband
    		dcc(:,iat) = dc_start(iat)
    	end do
    else

!    measure dos at the fermi level
    corr_dos = 0.0
    step_dc = 0.0

    do i=1,nband
	corr_dos = corr_dos + (-aimag(Gloc(i,i,1,1,iat))/pi)
    end do

!    change dcc
    do i=1, nband

    if(corr_dos > 1)then
	step_dc = ( ncorr(1,1,iat) - ncorr(2,1,iat) )/corr_dos
    else 
	step_dc = ( ncorr(1,1,iat) - ncorr(2,1,iat) )
    end if
    
    end do
    
    do i = 1, nband
    	dcc(:,iat) = dcc(:,iat) - (alpha*step_dc)/nband
    end do
    
    end if
    
      case(3)
!
!    Double counting is a fixed constant
!
    do i = 1, nband
    	dcc(:,iat) = dc_start(iat)
    end do

    end select

  if(myid==master)then
  write(*,*)' The double counting is set to', dcc(1,iat),'eV for atom',iat
  end if

!
!  Apply the dcc to eimp
  eimp(:,iat)=eimp(:,iat)+dcc(1,iat)

  enddo   ! natom

  call mp_barrier()
  
end subroutine plo_dcc
