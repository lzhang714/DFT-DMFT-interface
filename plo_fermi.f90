!
! Compute (see below) and write particle numbers and fix chemical potential (\mu)
!
! L.Zhang: created, then most recently modified 2017
!

subroutine plo_fermi()

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none
!
!Call to Ridders method
!
  if( myid == master )then
  write(*,*)'==================== Fermi level and Filling  ===================='
  write(*,*)''
  end if

  if( muman==0 )then
    call find_fermi(efermi)
  else
    efermi = manfermi
  end if

  call mp_barrier()
!
! compute occupancies and write them
!
  call calculate_occ(efermi)
!
!    Print out
!
  if( myid == master )then
    write(*,*)'================================================================'
    write(*,*)''
  end if! master

      call mp_barrier()

end subroutine plo_fermi

!
! Compute (see below) and write particle numbers for DFT, Sigma=0
!
subroutine plo_fermi_dft()

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none
!
!call to Ridders method
!
  if( myid == master )then
  write(*,*)'==================== DFT Fermi level  ===================='
  write(*,*)''
  end if

  call find_fermi(efermi_dft)

  call mp_barrier()
!
! write it
  if( myid == master )then
    write(*,"(' DFT Fermi level  ',20x,f15.11)") efermi_dft
    write(*,*)''
  end if
!
!  call calculate_occ(efermi)
!
!    Print out
!
  if( myid == master )then
    write(*,*)'=========================================================='
    write(*,*)''
  end if! master

      call mp_barrier()

end subroutine plo_fermi_dft

!
!   This subroutine calculates the total Green function and then the occupancy
!   Using Ridders' Method:
!
!   Ridders, C.; , "A new algorithm for computing a single root of a real continuous function,"
!   Circuits and Systems, IEEE Transactions on , vol.26, no.11, pp. 979- 980, Nov 1979
!   doi: 10.1109/TCS.1979.1084580
!
!   From numerical recipes. thats why it is so 'beautiful'.
!
subroutine find_fermi(zriddr)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none
   
   real(dp), intent(out) :: zriddr

  integer :: j
  integer :: MAXIT

  real(dp) :: fh,fl,fm,fnew,s,xh,xl,xm,xnew

  real(dp) :: x1,x2,xacc,UNUSED

! Some crazy small value
  UNUSED=-1.11E30
! Maximum number of iterations
  MAXIT=60
! Accuracy
  xacc = 1.d-5
!
! The interval is -10, 10, should be big enough.
!
  x1=-10.00_dp
  x2=10.00_dp
!
! Compute at the borders
!
  call calculate_diffocc(x1,fl)

  call calculate_diffocc(x2,fh)

  if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
   xl=x1
   xh=x2
   zriddr=UNUSED

   do j=1,MAXIT
!
!  Evaluate at midpoint
!
    xm=0.5*(xl+xh)

    call calculate_diffocc(xm,fm)
!
! Update
!
    s=sqrt(fm**2-fl*fh)
    if(s.eq.0.) return
    xnew=xm+(xm-xl)*(sign(1.,fl-fh)*fm/s)
    if (abs(xnew-zriddr).le.xacc) return
    zriddr=xnew

    call calculate_diffocc(zriddr,fnew)

    if (fnew.eq.0.) return

    if(sign(fm,fnew).ne.fm) then

    xl=xm
    fl=fm
    xh=zriddr
    fh=fnew

    else if(sign(fl,fnew).ne.fl) then

    xh=zriddr
    fh=fnew

    else if(sign(fh,fnew).ne.fh) then

    xl=zriddr
    fl=fnew

    else

    stop 'never get here in zriddr'

    endif

   if(abs(xh-xl).le.xacc)return

   enddo ! iteration done

   stop 'zriddr exceed maximum iterations'

   else if (fl.eq.0.) then

   zriddr=x1

   else if (fh.eq.0.) then

   zriddr=x2

   else

   stop 'root must be bracketed in zriddr'

   endif
   
end subroutine find_fermi
!
! Calculates Occupancy for given fermi level and write it.
!
subroutine calculate_occ(fermilevel)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none

  real(dp), intent(in )  :: fermilevel

  integer :: i
  integer :: j
  integer :: is
  integer :: iat

  real(dp) :: q0,numparts

  real(dp) :: occtot, occref

  real(dp) :: occ_array(3)

!
!   Compute Full Green functions
!
    efermi = fermilevel

    call plo_gk()

    occ(:,:,:,:) = 0.00_dp
    occ_array(:) = 0.00_dp
    numparts = 0.00_dp
    occtot = 0.00_dp
    ncorr(:,:,:) = 0.00_dp
!
!   Calculate occupation matrix
!
    do is = 1, nspin
!
!	Total occupancy first
!

       do i = 1, lda_nband
          q0 = 1.00_dp
          call qlmp( mfreq,beta,q0,occtot,Gbloc(i,:,is),rmesh)
          numparts = numparts + occtot
          occ_array(is) = occ_array(is) + occtot
       enddo

!
!	Partial Occupancies atom by atom
!
    do iat = 1,natom

    do i = 1, nband
        q0 = 1.00_dp
          call qlmp( mfreq,beta,q0,occref,G0loc(i,i,:,is,iat),rmesh)
          ncorr(1,is,iat) = ncorr(1,is,iat) + occref
          do j = 1, nband
             q0 = 0.00_dp
             if( i == j ) q0 = 1.D0
             call qlmp( mfreq,beta,q0,occ(i,j,is,iat),Gloc(i,j,:,is,iat),rmesh )
          enddo
          ncorr(2,is,iat) = ncorr(2,is,iat) + occ(i,i,is,iat)
       enddo
       
       enddo ! natom

    enddo ! spin

    occ_array(3) = numparts

      if( myid == master )then
!
!  Print it.
!
    write(*,"(' Number of particles for spin 1 : ',f15.11,/,                        &
               ' Number of particles for spin 2 : ',f15.11,/,                        & 
               ' Total magnetization              ',f15.11,/,                        &
               ' Total number of particles        ',f15.11)") occ_array(1), occ_array(2),  &
                                                   occ_array(2)-occ_array(1), occ_array(3)
    write(*,"(' Fermi level  ',20x,f15.11)") efermi
    write(*,*)''
    write(*,*)'Impurity:'
    write(*,*)''

     do iat = 1, natom
     write(*,*)'================================================================'
     write(*,"(' Atom ',i2)") iat
     write(*,*)'================================================================'
      do is = 1, nspin
        write(*,"(' Spin ',i2)") is
        do i = 1, nband
          write(*,"(20(2x,f7.4))") ( occ(i,j,is,iat), j = 1,nband )
        end do
        write(*,"(' Number of electrons in spin channel (G0)',f12.8)") ncorr(1,is,iat)
        write(*,"(' Number of electrons in spin channel (G) ',f12.8)") ncorr(2,is,iat)
        write(*,*)''
      end do
      
        write(*,"(' Total Impurity Filling (G0)',f12.8)") sum(ncorr(1,:,iat))
        write(*,"(' Total Impurity Filling (G) ',f12.8)") sum(ncorr(2,:,iat))
        write(*,*)''

     end do ! natom

      end if

end subroutine calculate_occ

!
! Calculates Occupancy for given fermi level and does not write it.
!
subroutine calculate_occ_silent(fermilevel)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none

  real(dp), intent(in )  :: fermilevel

  integer :: i
  integer :: j
  integer :: is
  integer :: iat

  real(dp) :: q0

  real(dp) :: occref

!
!   Compute Full Green functions
!
    efermi = fermilevel

    call plo_gk()

    occ(:,:,:,:) = 0.00_dp
    ncorr(:,:,:) = 0.00_dp
!
!   Calculate occupation matrix
!
    do is = 1, nspin
!
!       Partial Occupancies atom by atom only
!
    do iat = 1,natom

    do i = 1, nband
        q0 = 1.00_dp
          call qlmp( mfreq,beta,q0,occref,G0loc(i,i,:,is,iat),rmesh)
          ncorr(1,is,iat) = ncorr(1,is,iat) + occref
          do j = 1, nband
             q0 = 0.00_dp
             if( i == j ) q0 = 1.D0
             call qlmp( mfreq,beta,q0,occ(i,j,is,iat),Gloc(i,j,:,is,iat),rmesh )
          enddo
          ncorr(2,is,iat) = ncorr(2,is,iat) + occ(i,i,is,iat)
       enddo

       enddo ! natom

    enddo ! spin

end subroutine calculate_occ_silent

!
! Calculates DIFFERENCE between input occupancy ntotal and
! the one we get using fermilevel.
!
subroutine calculate_diffocc(fermilevel,numparts)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat
   use mmpi

   implicit none
   
   real(dp), intent(in)  :: fermilevel
   real(dp), intent(out) :: numparts
  
  integer :: i
  integer :: is
  
  real(dp) :: q0

  real(dp) :: occtot
!
!   Compute Full Green functions
!
    efermi = fermilevel

    call mp_barrier()

    call plo_gk()

    call mp_barrier()

    numparts = 0.00_dp
!
!   Calculate occupation matrix
!
    do is = 1, nspin

       do i = 1, lda_nband
          q0 = 1.00_dp
          call qlmp( mfreq,beta,q0,occtot,Gbloc(i,:,is),rmesh)
          numparts = numparts + occtot
       enddo

    enddo ! spin

    numparts = numparts-ntotal

    call mp_barrier()

end subroutine calculate_diffocc
!
! Qlmp subroutine from v2 of A. Poteryaevs code.
!
subroutine qlmp( nepus,beta1,q0,qlm,Gf11,om )

    use mmpi
    use constants

    implicit none

    integer,  intent (in   ) :: nepus
    real(dp), intent (in   ) :: beta1
    real(dp), intent (inout) :: q0
    real(dp), intent (  out) :: qlm
    real(dp), intent (in   ) :: om(nepus)
    complex(dp),   intent (in   ) :: Gf11(nepus)

    integer :: ie
    real(dp) :: q1
    complex(dp) :: sum

    q1 = 0.00_dp

    if( q0 > 1.d-2 ) q1 = real( 1.00_dp / Gf11(nepus) )

    sum  = czero
    do ie = 1,nepus
      sum = sum + Gf11(ie) - q0 / cmplx( q1,om(ie) )
    end do

    qlm = 2.00_dp * real(sum) / beta1

    if( q0 < 1.d-2 ) return

    q0 = beta1*q1
    if( abs( q0 ) < 5.d2 )then
       qlm = qlm + 1.00_dp / ( 1.00_dp + exp( -q0 ) )
      elseif( q0 >= 5.d2 )then
       qlm = qlm + 1.00_dp
    end if 
    
  end subroutine qlmp
