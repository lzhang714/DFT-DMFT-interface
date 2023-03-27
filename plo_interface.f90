!
! Reads Bloch hamiltonian and projectors on correlated subspace from input file. 
!
! Long Zhang created, most recently modified 2017
!
!
  subroutine plo_read_projectors()

    use constants
    use plo_quantities
    use control

    implicit none

    integer :: iat
    integer :: is
    integer :: ib
    integer :: idum
    integer :: ikp
    integer :: ios1
    integer :: ios2
    integer :: ios3
    integer :: ios4
    integer :: m
    integer :: LDWORK
    integer :: INFOR 
    integer :: ind
    logical :: exists_1
    logical :: exists_2
    logical :: exists_3
    logical :: exists_4

    real(dp), allocatable :: eigenvals(:)
    real(dp), allocatable :: RDWORK(:)

    complex(dp), allocatable   :: cbuff(:,:)
    complex(dp), allocatable   :: DWORK(:)
    complex(dp), allocatable   :: dbuff(:,:)
    complex(dp), allocatable   :: ubuff(:,:)
    complex(dp), allocatable   :: obuff(:,:)
    complex(dp), allocatable   :: locprold(:,:,:,:,:)
    complex(dp), allocatable   :: locbuff(:,:,:,:)
    complex(dp), allocatable   :: locbuffnew(:,:,:,:)


!
! Check if we have the files.
!
    inquire (file = 'BAND_up', exist = exists_1)
    inquire (file = 'BAND_down', exist = exists_2)
    inquire (file = 'PROJECTOR_up', exist = exists_3)
    inquire (file = 'PROJECTOR_down', exist = exists_4)

        if((exists_1==.false. .or. exists_2==.false. .or. exists_3==.false. .or. exists_4==.false.))then
    	    stop 'One or more of BAND_up/BAND_down ; PROJECTOR_up/down is missing!'
        end if

!
!   Open files.
!
    open( 2,file='PROJECTOR_up',form='formatted',status='old',iostat=ios1,&
          action='read',position='rewind' )
    open( 3,file='PROJECTOR_down',form='formatted',status='old',iostat=ios2,&
          action='read',position='rewind' )
    open( 4,file='BAND_up',form='formatted',status='old',iostat=ios3,&
          action='read',position='rewind' )
    open( 5,file='BAND_down',form='formatted',status='old',iostat=ios4,&
          action='read',position='rewind' )

    read(2,*) natom, nband, nkp, lda_nband
    read(3,*) natom, nband, nkp, lda_nband

    allocate( hbloch(lda_nband,nkp,nspin) )
!
!   Reading Fermi level and eigenvalues.
!
    do is =1,nspin
     read(is+3,*) efermi_dft
      do ikp = 1, nkp
       do ib = 1, lda_nband
        read(is+3,*) idum, idum, hbloch(ib,ikp,is)
       enddo
      enddo
     close(is+3)
    enddo

        hbloch(:,:,:) = hbloch(:,:,:) - efermi_dft
!
!   Reading projectors.
!
    allocate( locprold(nband,lda_nband,natom,nkp,nspin) )

    do is =1,nspin

     do ikp = 1, nkp
      do iat = 1, natom
       do m = 1, nband
        do ib = 1, lda_nband
         read(is+1,*) idum, idum, idum
         read(is+1,*) locprold(m,ib,iat,ikp,is)
        enddo
       enddo
      enddo
     enddo

     close(is+1)

    enddo
!
!       Construct the normalization matrix
!
    allocate( eigenvals(nband*natom), ubuff((nband*natom),(nband*natom)), obuff((nband*natom),(nband*natom)),& 
              cbuff((nband*natom),(nband*natom)), &
              dbuff((nband*natom),(nband*natom)), DWORK((2*(nband*natom)-1)), RDWORK((3*(nband*natom)-2)), &
              locpro(nband,lda_nband,natom,nkp,nspin), locbuff((nband*natom),lda_nband,nkp,nspin), & 
              locbuffnew((nband*natom),lda_nband,nkp,nspin) )
!
!       copy locprold(:,:,:,:) into locbuff(:,:,:)
!
        forall( iat=1:natom, ikp=1:nkp, m=1:nband, ib=1:lda_nband, is=1:nspin )
            locbuff((m+(nband*(iat-1))),ib,ikp,is)=locprold(m,ib,iat,ikp,is)
        end forall
!
!       set length of workarray of zheev routine (unoptimized)
!
       LDWORK = (2*(nband*natom)-1)

    do is=1, nspin

       do ikp = 1, nkp
!
!    Construct Overlap-Matrix
!
       call zgemm('N','C',(nband*natom),(nband*natom),lda_nband,cone,locbuff(1,1,ikp,is),(nband*natom),&
                  locbuff(1,1,ikp,is),(nband*natom),czero,cbuff,(nband*natom))
!
!    Compute Eigenvalues and -vectors of Overlap-Matrix
!
         call zheev('V','U',(nband*natom), cbuff(1,1), (nband*natom), eigenvals,DWORK,LDWORK, RDWORK,INFOR )
!
!       Form diagonalized Overlap-Matrix^-1/2 
!
        dbuff(:,:) = czero
        do ind = 1,(nband*natom)
         dbuff(ind,ind) = (1.d0)/(sqrt(eigenvals(ind)))
        enddo
!
!    Perform unitary transnformation to obtain Overlap^-1/2
!
         call zgemm('N','N',(nband*natom),(nband*natom),(nband*natom),cone,cbuff(1,1),(nband*natom),&
                     dbuff(1,1),(nband*natom),czero,ubuff,(nband*natom))

         call zgemm('N','C',(nband*natom),(nband*natom),(nband*natom),cone,ubuff(1,1),(nband*natom),&
                     cbuff(1,1),(nband*natom),czero,obuff,(nband*natom))
!
!       Overwrite locbuff with properly normalized locbuffnew
!
        call zgemm('N','N',(nband*natom),lda_nband,(nband*natom),cone,obuff(1,1),(nband*natom),&
                   locbuff(1,1,ikp,is),(nband*natom),czero,locbuffnew(1,1,ikp,is),(nband*natom))

       enddo ! ikp

       enddo ! is
!
!       write new locpro
!
        forall( iat=1:natom, ikp=1:nkp, m=1:nband, ib=1:lda_nband, is=1:nspin )
            locpro(m,ib,iat,ikp,is) = locbuffnew((m+(nband*(iat-1))),ib,ikp,is)
        end forall

    deallocate( cbuff,eigenvals,ubuff,obuff,dbuff,DWORK,RDWORK,locbuff,locbuffnew )

 end subroutine plo_read_projectors

!
! Compute the Hopping matrix
!

subroutine plo_make_hopmat()
    use constants
    use plo_quantities
    use control

      implicit none

    character (LEN=30) :: cdum
    integer :: i,j,iat, ib, is, ikp, ios, m, LDWORK, idum, INFOR
    double precision, allocatable :: eigenvals(:),RDWORK(:),kpoints(:,:)
    double complex, allocatable   :: DWORK(:),hbuff(:,:),ebuff(:,:),Hk(:,:),locbuff(:,:,:,:)
    logical :: exists
!
!   Open files.
!

    inquire (file = 'solver.kpoints.in', exist = exists)
    
    if(exists==.false.)then
	stop 'solver.kpoints.in file not found!'
    end if

    open( 888,file='solver.kpoints.in',form='formatted',status='old',iostat=ios,    &
          action='read',position='rewind' )

    allocate( kpoints(nkp,3) )

!
!  Read Kpoints from kpoints file
!

!
!  header stuff
!
    read(888,*)cdum
    read(888,*)idum
    read(888,*)cdum

    do ikp = 1, nkp
       read(888,*)kpoints(ikp,1),kpoints(ikp,2),kpoints(ikp,3),idum
    enddo

    close(888)
!
!   Construct the normalization matrix..
!

    allocate( eigenvals(nband*natom), DWORK((2*(nband*natom)-1)), RDWORK((3*(nband*natom)-2)), &
	      locbuff((nband*natom),lda_nband,nkp,nspin), & 
	      hbuff((nband*natom),lda_nband),ebuff(lda_nband,lda_nband),Hk((nband*natom),(nband*natom)) )
!
!	set hopping matrix to zero
!
    hopmat(:,:) = czero
!
!	shove locpro(:,:,:,:,:) into locbuff(:,:,:,:)
!
	forall(iat=1:natom, ikp=1:nkp, m=1:nband, ib=1:lda_nband,is=1:nspin)
	    locbuff((m+(nband*(iat-1))),ib,ikp,is)=locpro(m,ib,iat,ikp,is)
	end forall
!
!   	set length of workarray of zheev (unoptimized)
!
    LDWORK = (2*(nband*natom)-1)

    do ikp = 1, nkp   
!
!	compute H(k)
!	
	ebuff(:,:) = czero
	    
	do ib = 1,lda_nband
	    ebuff(ib,ib) = hbloch(ib,ikp,1)
	enddo ! ib
	    
	DWORK(:) = czero
	RDWORK(:) = czero
	Hk(:,:) = czero
	hbuff(:,:) = czero
	eigenvals(:) = czero

!	Computing H(k)...
	call zgemm('N','N',(nband*natom),lda_nband,lda_nband,cone,locbuff(1,1,ikp,1),(nband*natom), &
		    ebuff(1,1),lda_nband,czero,hbuff(1,1),(nband*natom))

	call zgemm('N','C',(nband*natom),(nband*natom),lda_nband,cone,hbuff(1,1),(nband*natom), &
		    locbuff(1,1,ikp,1),(nband*natom),czero,Hk(1,1),(nband*natom))
!
!	compute on-site hopping matrix elements Hop=1/N sum H(k)
!
	forall(i=1:(nband*natom), j=1:(nband*natom))
	    hopmat(i,j) = hopmat(i,j)+Hk(i,j)
	end forall
!
!	Diagonalize H(k)
!	    
	call zheev('N','U',(nband*natom),Hk(1,1),(nband*natom),eigenvals,DWORK,LDWORK,RDWORK,INFOR)
	    
    enddo ! ikp

!
!	normalize and write the hopping matrix
!

	forall(i=1:(nband*natom), j=1:(nband*natom))
	    hopmat(i,j) = hopmat(i,j)/(nkp)
	end forall

    deallocate(eigenvals,DWORK,RDWORK,locbuff,hbuff,ebuff,Hk,kpoints)

50  format(x,10f8.5)
99  format(('(',F8.5$))
55  format((1X,F8.5,')'$))
77  format((1X,F8.5$))

end subroutine plo_make_hopmat

!
! Diagonalize the on site Hamiltonian to obtain the crystal fiel basis
!
 subroutine plo_make_cfb()
    use constants
    use plo_quantities
    use control
    use mmpi

    implicit none  
 
    integer :: i, j, iat, LDWORK, INFOR, m, ib, ikp
    real(dp) :: FindDet
    double precision, allocatable :: eigenvals(:),RDWORK(:)
    double complex, allocatable   :: DWORK(:), locpro_cfb(:,:,:,:),hopbuff(:,:),hopblock(:,:)

! characters for plotting
    character(len=12) :: lmchar(5)=(/'*dxy(u,v)   ','*dyz(u,v)   ','*d3z2m1(u,v)','*dzx(u,v)   ','*dx2my2(u,v)'/)

    if(myid==master)then
	open( 99,file='solver.hopmat.dat',form='formatted',action='write',position='append')
	open( 999,file='orbital_plot.dat',form='formatted',action='write',position='append')
    end if

	allocate( eigenvals(nband), DWORK((2*nband-1)), RDWORK((3*nband-2)), & 
		  hopbuff((nband*natom),(nband*natom)),locpro_cfb(nband,lda_nband,natom,nkp), hopblock(nband,nband))


	forall(i=1:(nband*natom), j=1:(nband*natom))
	    hopbuff(i,j) = hopmat(i,j)
	end forall

	LDWORK = (2*nband-1)

      do iat=1,natom

	hopblock(:,:) = czero

	forall(i=1:nband, j=1:nband)
		hopblock(i,j) = hopbuff(i+((iat-1)*nband),j+((iat-1)*nband))
	end forall

	eigenvals(:) = czero
	DWORK(:) = czero
	RDWORK(:) = czero
!
!	diagonalize onsite blocks in hopmat
!	    
	call zheev('V','U',nband,hopblock,nband,eigenvals,DWORK,LDWORK,RDWORK,INFOR )
	
	do i=1,nband
		do j=1,nband
		cfb_coeff(iat,i,j) = real(hopblock(j,i))
		cfb_det(iat,i,j) = real(hopblock(j,i))
		enddo
	enddo

    if(myid==master)then

	write(99,*)'eigenvectors from cfb_coeff for checking'
	do i=1,nband
		do j=1,nband
		write(99,*)iat,i,j, cfb_coeff(iat,i,j)
		enddo
		write(99,*)''
	enddo
	write(99,*)''
!
!	Eigenvectors for the plot. in line.
!
	write(999,*)'#Atom',iat
	write(999,*)'#============'
	do j=1,nband
	write(999,'(a)')'orb1(u,v)=\'
		do i=1,nband
			write(999,777)real(hopblock(i,j)),trim(lmchar(i))
		enddo
	write(999,*)
	write(999,'(a)')'pos1(u,v)= orb1(u,v)<0 ? 0 :  orb1(u,v)'
	write(999,'(a)')'neg1(u,v)= orb1(u,v)>0 ? 0 : -orb1(u,v)'
        write(999,'(a)')'set autoscale'
	write(999,'(a)')'splot pos1(u,v)*x(u,v), pos1(u,v)*y(u,v), pos1(u,v)*z(u,v),\'
	write(999,'(a)')'neg1(u,v)*x(u,v), neg1(u,v)*y(u,v), neg1(u,v)*z(u,v)'
	write(999,*)
	
        enddo
	write(999,*)''

    write(99,*)'eigenvectors (read in column order) '

	do i=1,nband
		do j=1,nband
			write(99,77)real(hopblock(i,j))
		enddo
		write(99,*)''
		write(99,*)
	enddo
	write(99,*)''

	    write(99,*)'crystal field splitting'
	    
		do i=1,nband
			write(99,77)eigenvals(i)
	    enddo
	    
	    write(99,*)''
	    write(99,*)''
	    
	end if ! master
	    
 enddo !iat

!
!	Transform into crystal field basis.
!
!	Be careful! This procedure rearranges the orbitals. The VASP order is destroyed
!	and the orbitals are ordered by their on-site eigenvalues.
!	This is OK, since the Coulomb interaction is rotated to the same basis later.
!

	locpro_cfb(:,:,:,:)=czero


	do i=1,nband
		forall(m=1:nband, ib=1:lda_nband, iat=1:natom, ikp=1:nkp)
	    	locpro_cfb(m,ib,iat,ikp) = locpro_cfb(m,ib,iat,ikp)+(cfb_coeff(iat,m,i)*locpro(i,ib,iat,ikp,1))
		end forall
	enddo !nband

    if(myid==master)then
    write(99,*)'old onsite energy matrix'
	do i=1,(nband*natom)
	    do j=1,(nband*natom)
		write(99,77)real(hopmat(i,j))
	    enddo
	    write(99,*)''
	enddo
    end if ! master

    locpro(:,:,:,:,1) = locpro_cfb(:,:,:,:)
    
    locpro(:,:,:,:,2) =  locpro(:,:,:,:,1)

    deallocate(eigenvals,DWORK,RDWORK,hopbuff,locpro_cfb,hopblock)

    if(myid==master) then
        write(99,*)''
	close(99)
	close(999)
!
! check if the transformation was unitary, warn if it was not
!
    write(*,*)'Making Crystal Field Basis'

    do iat=1,natom
	write(*,*)'Determinant of transformation matrix for atom',iat,abs(FindDet(cfb_det(iat,:,:), nband))
    enddo
    
    end if ! master

50  format(x,10f8.5)
99  format(('(',F10.7$))
55  format((1X,F10.7,')'$))
77  format((1X,F10.7$))
777  format((SP,F10.7,a$))

end subroutine plo_make_cfb

!
! Diagonalize  the on site density matrix to obtain something
!
 subroutine plo_diag_densmat()
    use constants
    use plo_quantities
    use control
    use mmpi

    implicit none  
 
    integer :: i, j, iat, LDWORK, INFOR, m, ib, ikp
    real(dp) :: FindDet
    double precision, allocatable :: eigenvals(:),RDWORK(:)
    double complex, allocatable   :: DWORK(:), locpro_cfb(:,:,:,:),hopbuff(:,:),hopblock(:,:)

! characters for plotting
    character(len=12) :: lmchar(5)=(/'*dxy(u,v)   ','*dyz(u,v)   ','*d3z2m1(u,v)','*dzx(u,v)   ','*dx2my2(u,v)'/)

    if(myid==master)then
        open( 99,file='solver.hopmat.dat',form='formatted',action='write',position='append')
        open( 999,file='orbital_plot.dat',form='formatted',action='write',position='append')
    end if

        allocate( eigenvals(nband), DWORK((2*nband-1)), RDWORK((3*nband-2)), & 
                  hopbuff((nband*natom),(nband*natom)),locpro_cfb(nband,lda_nband,natom,nkp), hopblock(nband,nband) )

        LDWORK = (2*nband-1)

        call calculate_occ_silent(0.0)

      do iat=1,natom

        hopblock(:,:) = czero

        forall(i=1:nband, j=1:nband)
                hopblock(i,j) = occ(i,j,1,iat)
        end forall

        eigenvals(:) = czero
        DWORK(:) = czero
        RDWORK(:) = czero
!
!       diagonalize onsite density matrix
!           
        call zheev('V','U',nband,hopblock,nband,eigenvals,DWORK,LDWORK,RDWORK,INFOR )

        do i=1,nband
                do j=1,nband
                cfb_coeff(iat,i,j) = real(hopblock(j,i))
                cfb_det(iat,i,j) = real(hopblock(j,i))
                enddo
        enddo

    if(myid==master)then

        write(99,*)'eigenvectors from cfb_coeff for checking'
        do i=1,nband
                do j=1,nband
                write(99,*)iat,i,j, cfb_coeff(iat,i,j)
                enddo
                write(99,*)''
        enddo
        write(99,*)''

!
!       Eigenvectors for the plot. in line.
!
        write(999,*)'#Atom',iat
        write(999,*)'#============'
        do j=1,nband
        write(999,'(a)')'orb1(u,v)=\'
                do i=1,nband
                        write(999,777)real(hopblock(i,j)),trim(lmchar(i))
                enddo
        write(999,*)
        write(999,'(a)')'pos1(u,v)= orb1(u,v)<0 ? 0 :  orb1(u,v)'
        write(999,'(a)')'neg1(u,v)= orb1(u,v)>0 ? 0 : -orb1(u,v)'
        write(999,'(a)')'set autoscale'
        write(999,'(a)')'splot pos1(u,v)*x(u,v), pos1(u,v)*y(u,v), pos1(u,v)*z(u,v),\'
        write(999,'(a)')'neg1(u,v)*x(u,v), neg1(u,v)*y(u,v), neg1(u,v)*z(u,v)'
        write(999,*)

        enddo
        write(999,*)''

    write(99,*)'eigenvectors (read in column order) '

        do i=1,nband
                do j=1,nband
                        write(99,77)real(hopblock(i,j))
                enddo
                write(99,*)''
                write(99,*)
        enddo
        write(99,*)''

            write(99,*)'crystal field splitting'
            
                do i=1,nband
                        write(99,77)eigenvals(i)
            enddo
            
            write(99,*)''
            write(99,*)''
            
        end if ! master
            
 enddo !iat

!
!       Transform into crystal field basis.
!
!
!       Be careful! This procedure rearranges the orbitals. The VASP order is destroyed
!       and the orbitals are ordered by their on-site eigenvalues.
!

        locpro_cfb(:,:,:,:)=czero


        do i=1,nband
                forall(m=1:nband, ib=1:lda_nband, iat=1:natom, ikp=1:nkp)
                locpro_cfb(m,ib,iat,ikp) = locpro_cfb(m,ib,iat,ikp)+(cfb_coeff(iat,m,i)*locpro(i,ib,iat,ikp,1))
                end forall
        enddo !nband

    if(myid==master)then
    write(99,*)'old onsite energy matrix'
        do i=1,(nband*natom)
            do j=1,(nband*natom)
                write(99,77)real(hopmat(i,j))
            enddo
            write(99,*)''
        enddo
    end if ! master

    locpro(:,:,:,:,1) = locpro_cfb(:,:,:,:)
    
    locpro(:,:,:,:,2) =  locpro(:,:,:,:,1)

    deallocate(eigenvals,DWORK,RDWORK,hopbuff,locpro_cfb,hopblock)

    if(myid==master) then
        write(99,*)''
        close(99)
        close(999)
!
! check if the transformation was unitary, warn if it was not
!
    write(*,*)'Diagonalizing local density matrix'

    do iat=1,natom
	write(*,*)'Determinant of transformation matrix for atom',iat,abs(FindDet(cfb_det(iat,:,:), nband))
    enddo

    end if ! master

50  format(x,10f8.5)
99  format(('(',F10.7$))
55  format((1X,F10.7,')'$))
77  format((1X,F10.7$))
777  format((SP,F10.7,a$))

end subroutine plo_diag_densmat

!
! Perform downfolding from Bloch to Local space via PLO projectors
!
  subroutine plo_downfold(bloch_array,plo_projectors,local_array)

   use constants
   use plo_quantities
   use control
   use mmpi

   implicit none

!External
   complex(dp), intent(in    ) :: bloch_array(lda_nband,lda_nband)
   complex(dp), intent(in    ) :: plo_projectors(nband,lda_nband)
   complex(dp), intent(inout ) :: local_array(nband,nband)

!Local
   integer :: bloch_dim
   integer :: local_dim
   complex(dp), allocatable :: buffer(:,:)

!Get dimensions. I know, redundant. Deal with it.
    bloch_dim=size(bloch_array,dim=1)
    local_dim=size(local_array,dim=1)

    local_array(:,:) = czero

    allocate( buffer(local_dim, bloch_dim) )

!Downfolding. Matrix multiplication...
    call zgemm('N','N',local_dim,bloch_dim,bloch_dim,cone,plo_projectors(1,1),local_dim,&
                bloch_array,bloch_dim,czero,buffer,local_dim)
    call zgemm('N','C',local_dim,local_dim,bloch_dim,cone,buffer,local_dim,&
                plo_projectors(:,:),local_dim,czero,local_array,local_dim)

  deallocate( buffer )

  end subroutine plo_downfold
!
! Perform upfolding from Local to Bloch space via PLO projectors
! (multiatom version)
!
  subroutine plo_upfold(local_array,plo_projectors,bloch_array)

   use constants
   use plo_quantities
   use control
   use mmpi

   implicit none

! External
   complex(dp), intent(in    ) :: local_array(nband,nband,natom)
   complex(dp), intent(in    ) :: plo_projectors(nband,lda_nband,natom)
   complex(dp), intent(inout ) :: bloch_array(lda_nband,lda_nband)

! Local
   integer :: bloch_dim
   integer :: local_dim
   integer :: iat
   integer :: m
   integer :: ib
   integer :: i
   integer :: j

   complex(dp), allocatable :: buffer(:,:)
   complex(dp), allocatable :: sigloctot(:,:)
   complex(dp), allocatable :: locpro_mult(:,:)

! Get dimensions. I know, redundant, but still useful. Deal with it.
    bloch_dim=size(bloch_array,dim=1)
    local_dim=size(local_array,dim=1)

    bloch_array(:,:) = czero

    allocate( buffer(bloch_dim, (local_dim*natom)),sigloctot( (local_dim*natom),(local_dim*natom) ),&
              locpro_mult((local_dim*natom),bloch_dim) )

    sigloctot(:,:) = czero
    buffer(:,:) = czero
    locpro_mult(:,:) = czero
!
! Build multiatomic locpro
!
     forall( iat=1:natom, m=1:local_dim, ib=1:bloch_dim)
         locpro_mult((m+(local_dim*(iat-1))),ib)=plo_projectors(m,ib,iat)
     end forall
!
! Build Sigma for all atoms
!
    forall( iat=1:natom, i=1:local_dim, j=1:local_dim )
        sigloctot(((iat-1)*local_dim)+i,((iat-1)*local_dim)+j) = local_array(i,j,iat)
    end forall
!
!    Transform local matrix to bloch basis.
!
    call zgemm('C','N',bloch_dim,(local_dim*natom),(local_dim*natom),cone,locpro_mult(1,1),(local_dim*natom),&
                 sigloctot(1,1),(local_dim*natom),czero,buffer,bloch_dim)
    call zgemm('N','N',bloch_dim,bloch_dim,(local_dim*natom),cone,buffer,bloch_dim,&
                locpro_mult(1,1),(local_dim*natom),czero,bloch_array(1,1),bloch_dim)

  deallocate( buffer,sigloctot,locpro_mult )

  end subroutine plo_upfold
!
! Rewrite Selfenergy output from the solver into plo friendly format
! or the other way around.
!
  subroutine plo_rewrite_array(solver_array,plo_array,controltag)

   use constants
   use plo_quantities
   use control
   use mmpi

   implicit none

! External
   complex(dp), intent(inout ) :: solver_array(mfreq,norbs,norbs,natom)
   complex(dp), intent(inout ) :: plo_array(nband,nband,mfreq,nspin,natom)

! Control string, input should be 'plo2solver' with quotes
   character(10), intent(in    ) :: controltag

   integer :: i
   integer :: iom
   integer :: iat

! Decide what to do
   select case(controltag)

! We want to go from the plo array to the solver conform array
   case('plo2solver')

   do iat=1,natom
    do iom=1,mfreq
     do i = 1,nband
      solver_array(iom,i,i,iat) = plo_array(i,i,iom,1,iat)
      solver_array(iom,i+nband,i+nband,iat) = plo_array(i,i,iom,2,iat)
     enddo ! nband
    enddo ! mfreq
   enddo ! natom

! We want to go from the solver array to the plo conform array
   case('solver2plo')

   do iat=1,natom
    do iom=1,mfreq
     do i = 1,nband
      plo_array(i,i,iom,1,iat) = solver_array(iom,i,i,iat)
      plo_array(i,i,iom,2,iat) = solver_array(iom,i+nband,i+nband,iat)
     enddo ! nband
    enddo ! mfreq
   enddo ! natom

   end select

  end subroutine plo_rewrite_array

  subroutine plo_read_selfenergy(sigf)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat, only : rmesh

   implicit none

! self-energy function
     complex(dp), intent(inout) :: sigf(nband,nband,mfreq,nspin)

! local variables
! loop index
     integer :: i
     integer :: j
     integer :: is
     integer :: iom
     integer :: iat

! open data file: solver.sgm.in
     open(mytmp, file='solver.sgm.in', form='formatted', status='old',action='read',&
     position='rewind')

! read it directly into sigma
    do is=1,nspin
     do i=1,nband
     do j=i,nband
      read(mytmp,*) ( sigf(i,j,iom,is), iom=1,mfreq )
     enddo
     enddo
    enddo

! close data file
     close(mytmp)

!
!    Make sigma symmetric and complete second spin
!
   do iat=1,natom
    do i = 1, nband-1
      do j = i+1, nband
       sigma(j,i,:,:,iat) = sigma(i,j,:,:,iat)
      end do
    end do
   end do

  end subroutine plo_read_selfenergy

!
! Fourier Transform from Gbloch(iw) to Gbloch(tau)
!
  subroutine plo_gbloc_gbtau(green_bloch_mat)

   use constants
   use plo_quantities
   use control
   use ctqmc_umat, only : cmesh,tmesh
   use mmpi
   !! use m_fft

   implicit none

!External
   complex(dp), intent(in    ) :: green_bloch_mat(lda_nband,mfreq,nspin)

!Local
   integer :: is
   integer :: i
   integer :: iom
   integer :: iflag

   complex(dp), allocatable :: gaux_w(:,:,:)
   real(dp), allocatable :: gaux_tau(:,:,:)
!
!  Printout Bloch Green's function in time domain. Only the diagonal.
!
    allocate( gaux_w(0:mfreq-1,lda_nband,nspin), gaux_tau(ntime,lda_nband,nspin) )

    do is = 1,nspin
       do iom = 0,mfreq-1
          forall( i = 1:lda_nband ) gaux_w(iom,i,is) = green_bloch_mat(i,iom+1,is)
       enddo
    enddo

    do is = 1,nspin
      do i = 1,lda_nband
        iflag = 0
         !call invfourier( gaux_w(:,i,is),mfreq-1,gaux_tau(:,i,is),ntime,beta,iflag,efermi )
      end do
    end do

    call out_gb( gaux_tau )

    deallocate(gaux_w,gaux_tau)

  end subroutine plo_gbloc_gbtau
!
!MK Stolen from A. Poteryaevs code.
!
  subroutine out_gb( g )

    use constants
    use control
    use plo_quantities
!
!   Write to file Gbtau_i,i.dat Bloch Green functions
!
    implicit none
    real(dp), intent (in   ) :: g(ntime,lda_nband,nspin)     !  Green function

    integer :: i,k,is,ityp
    character*6  :: pos
    character*20 :: fname,fname_i
    real(dp) :: abssum
    real(dp) :: f(ntime)

    pos = 'rewind'

    do i = 1,lda_nband
        ityp = i
        f = g(:,ityp,1)
        abssum = sum( abs(f) ) / float(ntime)
        if( abssum > 1.d-6 )then
          write(fname_i,*) i
          fname_i = trim(adjustl(fname_i))
          fname   = 'Gbtau_'//fname_i
          fname   = trim(adjustl(fname))//','
          fname   = trim(adjustl(fname))//fname_i
          fname   = trim(adjustl(fname))//'.dat'

          open( 10+ityp,file=fname,form='formatted',status='unknown',position=pos )

          do k = 1,ntime
             write(10+ityp,50) beta*( k - 1 ) / float(ntime),(g(k,ityp,is),is=1,nspin)
          enddo

!          if( i == j )then
             write(10+ityp,50) beta,(-1.D0-g(1,ityp,is),is=1,nspin)
!          else
!             write(10+ityp,50) beta,(-g(1,ityp,is),is=1,nspin)
!          end if
!          write(10+ityp,*)

          close(10+ityp)
        end if
    end do

50  format(20(1x,f14.8))

  end subroutine out_gb
!
! Compute Selfenergy the old fashioned way and write it into sig2 array
! agrees with ctqmc internat subroutine
  subroutine plo_calc_sigma(grnf,grnf0)

   use constants
   use plo_quantities
   use control
   use mmpi

   implicit none

!External
   complex(dp), intent(in    ) :: grnf(mfreq,norbs,norbs,natom)
   complex(dp), intent(in    ) :: grnf0(mfreq,norbs,norbs,natom)

!Local
   integer :: iom
   integer :: iat
   integer :: is
   integer :: i
 
   complex(dp) :: grnf_loc(nband,nband,mfreq,nspin,natom)
   complex(dp) :: grnf_loc0(nband,nband,mfreq,nspin,natom)

   grnf_loc = czero
   grnf_loc0 = czero

! Rewrite array in plo form
   call plo_rewrite_array(grnf,grnf_loc,'solver2plo')
   call plo_rewrite_array(grnf0,grnf_loc0,'solver2plo')

! Compute Sigma

   sigma = czero

  do iat=1,natom

   do is =1,nspin
    do iom=1,mfreq
    call ctqmc_zmat_inv( nband,grnf_loc(:,:,iom,is,iat) )
    call ctqmc_zmat_inv( nband,grnf_loc0(:,:,iom,is,iat) )
    sigma(:,:,iom,is,iat) = grnf_loc0(:,:,iom,is,iat) - grnf_loc(:,:,iom,is,iat)
    enddo
   enddo

  enddo ! atoms

  end subroutine plo_calc_sigma

!
! compute the determinant of a matrix
! stolen from http://web.hku.hk/~gdli/UsefulFiles/Example-Fortran-program.html
!
!Function to find the determinant of a square matrix
!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
 REAL(dp) FUNCTION FindDet(matrix, n)
    use constants

    IMPLICIT NONE
    REAL(dp), DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL(dp) :: m, temp
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
!Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO

!Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO

 END FUNCTION FindDet
