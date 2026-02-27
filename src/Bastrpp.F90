!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for Viscous Burgers equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2025-07-09 by Chensheng Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Bastrpp
  !
  use commvar
  use parallel
  use fftwlink
  use tool
  use solution
  use readwrite
  use utility
  use cmdefne
  include 'fftw3-mpi.f03'
  !
  character(len=16) :: cmd
  integer :: thefilenumb
  real(8) :: scale
  !
  call mpiinitial
  !
  call readinput
  !
  call fftw_mpi_init()
  if(mpirank==0)  print *, "fftw_mpi initialized"
  !
  call mpisizedis_fftw
  if(mpirank==0)  print*, '** mpisizedis & parapp done!'
  !
  call parallelini
  if(mpirank==0)  print*, '** parallelini done!'
  !
  call allocation
  if(mpirank==0)  print*, '** allocation done!'
  !
  ! Planning
  call prepareplan_fftw_plan
  if(mpirank==0)  print*, '** fftw planing done!'
  !
  if(mpirank==0)  print*, '** Entering pp process...'
  !
  if(mpirank == 0) then
    call readkeyboad(cmd)
    print*,' ** pp command: ',cmd
  endif
  !
  call bcast(cmd)
  !
  if(trim(cmd)=='addfield2d') then
    call add_field_2d
  elseif(trim(cmd)=='transfer2d')then
    if(mpirank == 0) then
      call readkeyboad(cmd)
      read(cmd,'(i4)') thefilenumb
    endif
    call bcast(thefilenumb) 
    call transfer2d(thefilenumb)
  elseif(trim(cmd)=='spectra2d')then
    if(mpirank == 0) then
      call readkeyboad(cmd)
      read(cmd,'(i4)') thefilenumb
    endif
    call bcast(thefilenumb) 
    call spectra2d(thefilenumb)
  elseif(trim(cmd)=='hitgen2d')then
    call create_initial_field_2d
  elseif(trim(cmd)=='hitgen2dbis')then
    call create_initial_field_2d_bis
  elseif(trim(cmd)=='scaleinit2d')then
    if(mpirank == 0) then
      call readkeyboad(cmd)
      read(cmd,*) scale
    endif
    call bcast(scale)
    call scale_initial_field_2d(scale)
  else
    if(mpirank==0) print *, ' Unknown pp command!'
  endif
  !
  call deallocation
  !
  if(mpirank==0)  print *, "** End!"
  !
end program Bastrpp

subroutine add_field_2d
  use stlaio,  only: get_unit
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  !
  character(len=256) :: infilename      ! Name of the input HDF5 file
  integer :: ierr                      ! Error code for MPI operations
  character(len=1) :: modeio           ! Mode for HDF5 read operations
  real(8), allocatable :: u1A(:,:,:), u2A(:,:,:), u1B(:,:,:),u2B(:,:,:)
  real(8) :: cA, cB
  character(len=64) :: inputfile
  integer :: n,fh,i
  !
  allocate(u1A(1:im,1:jm,0:km), u2A(1:im,1:jm,0:km), &
            u1B(1:im,1:jm,0:km), u2B(1:im,1:jm,0:km))
  ! 
  !
  if(mpirank==0) then
    !
    inputfile='datin/addfield.dat'
    !
    open(fh,file=trim(inputfile),action='read')
    read(fh,'(////)')
    read(fh,*)cA,cB
    print *,' << ',trim(inputfile),' ... done'
  endif
  !
  call bcast(cA)
  call bcast(cB)
  !
  if(mpirank==0) then
    !
    print *, 'cA=', cA, ' cB=', cB
    !
  endif
  !
  ! read field
  modeio ='h'
  infilename='datin/flowini2dA.h5'
  !
  ! Initial field read
  !
  call h5io_init(filename=infilename,mode='read')
  call h5read(varname='u1', var=u1A(1:im,1:jm,0:km),mode = modeio)
  call h5read(varname='u2', var=u2A(1:im,1:jm,0:km),mode = modeio)
  call h5io_end
  call mpi_barrier(mpi_comm_world,ierr)
  if(mpirank==0) print *, ' << ',trim(infilename),' ... done'
  !
  infilename='datin/flowini2dB.h5'
  !
  call h5io_init(filename=infilename,mode='read')
  call h5read(varname='u1', var=u1B(1:im,1:jm,0:km),mode = modeio)
  call h5read(varname='u2', var=u2B(1:im,1:jm,0:km),mode = modeio)
  call h5io_end
  call mpi_barrier(mpi_comm_world,ierr)
  if(mpirank==0) print *, ' << ',trim(infilename),' ... done'
  !
  do i = 1,im
    do j = 1,jm
        u1(i,j,0) = cA * u1A(i,j,0) + cB * u1B(i,j,0)
        u2(i,j,0) = cA * u2A(i,j,0) + cB * u2B(i,j,0)
    end do
  end do
  !
  infilename='datin/sumfield.'//modeio//'5'
  !
  call h5io_init(trim(infilename),mode='write')
  call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode=modeio)
  call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode=modeio)
  call h5io_end
  !
  if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
  !
end subroutine add_field_2d
!
subroutine transfer2d(thefilenumb)
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  use solution
  use utility, only: listinit, listwrite
  !
  implicit none
  integer, intent(in)::thefilenumb
  character(len=128) :: infilename
  character(len=4) :: stepname
  real(8), allocatable :: trans(:,:,:),trans1(:,:,:),trans2(:,:,:)
  real(8), allocatable :: T(:),T1(:),T2(:)
  integer :: kOrdinal
  real(8) :: kk
  integer :: hand_f, i,j
  !
  call quantity_prepare
  !
  allocate(trans(1:im,1:jm,0:km),trans1(1:im,1:jm,0:km),trans2(1:im,1:jm,0:km))
  allocate(T(0:allkmax),T1(0:allkmax),T2(0:allkmax))
  ! read field
  write(stepname,'(i4.4)')thefilenumb
  infilename='outdat/flowfield'//stepname//'.h5'
  !
  call h5io_init(filename=infilename,mode='read')
  call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = 'h')
  call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = 'h')
  call h5read(varname='time',var=time)
  call h5read(varname='nstep', var=nstep)
  call h5io_end
  if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
  !
  do j=1,jm
  do i=1,im
  u1spe(i,j,1)=CMPLX(u1(i,j,0),0.d0,C_INTPTR_T);
  u2spe(i,j,1)=CMPLX(u2(i,j,0),0.d0,C_INTPTR_T);
  end do
  end do
  !
  call fft2d(u1spe)
  call fft2d(u2spe)
  !
  ! In fourier space
  u1x1 = imag * u1spe * k1 
  u1x2 = imag * u1spe * k2
  u2x1 = imag * u2spe * k1
  u2x2 = imag * u2spe * k2
  !
  call ifft2d(u1x1)
  call ifft2d(u1x2)
  call ifft2d(u2x1)
  call ifft2d(u2x2)
  call ifft2d(u1spe)
  call ifft2d(u2spe)
  !
  ! In physical spac
  do j=1,jm
  do i=1,im
      !
      ! uitA = uj \partial ui/ \partial uj
      u1tA(i,j,1) = - u1spe(i,j,1) * dreal(u1x1(i,j,1)) - u2spe(i,j,1) * dreal(u1x2(i,j,1))
      u2tA(i,j,1) = - u1spe(i,j,1) * dreal(u2x1(i,j,1)) - u2spe(i,j,1) * dreal(u2x2(i,j,1)) 
      !
      ! uitC = ui \partial uj/ \partial uj
      u1tC(i,j,1) = u1spe(i,j,1) * dreal(u1x1(i,j,1) + u2x2(i,j,1))
      u2tC(i,j,1) = u2spe(i,j,1) * dreal(u1x1(i,j,1) + u2x2(i,j,1))
      !
      ! ATTENTION: Reallocation: uixj for uiuj
      u1x1(i,j,1) = u1spe(i,j,1)*u1spe(i,j,1)
      u1x2(i,j,1) = u1spe(i,j,1)*u2spe(i,j,1)
      u2x1(i,j,1) = u2spe(i,j,1)*u1spe(i,j,1)
      u2x2(i,j,1) = u2spe(i,j,1)*u2spe(i,j,1)
  end do
  end do
  !
  call fft2d(u1tA)
  call fft2d(u2tA)
  call fft2d(u1spe)
  call fft2d(u2spe)
  call fft2d(u1tC)
  call fft2d(u2tC)
  call fft2d(u1x1)
  call fft2d(u1x2)
  call fft2d(u2x1)
  call fft2d(u2x2)
  ! In Fourier space
  !
  T = 0.d0
  T1 = 0.d0
  T2 = 0.d0
  do j=1,jm
  do i=1,im
      !
      kk=dsqrt(k1(i,j,0)**2+k2(i,j,0)**2)
      !
      trans(i,j,0) = dreal(conjg(u1spe(i,j,1))*u1tA(i,j,1)) + dreal(conjg(u2spe(i,j,1))*u2tA(i,j,1))
      trans1(i,j,0) = dimag(k1(i,j,0)*u1x1(i,j,1)*conjg(u1spe(i,j,1))+&
                            k2(i,j,0)*u1x2(i,j,1)*conjg(u1spe(i,j,1))+&
                            k1(i,j,0)*u2x1(i,j,1)*conjg(u2spe(i,j,1))+&
                            k2(i,j,0)*u2x2(i,j,1)*conjg(u2spe(i,j,1)))
      trans2(i,j,0) = dreal(conjg(u1spe(i,j,1))*u1tC(i,j,1)) + dreal(conjg(u2spe(i,j,1))*u2tC(i,j,1))
      !
      kOrdinal = kint(kk,1.d0,2,1.d0)
      !
      if(kOrdinal <= allkmax)then
          T(kOrdinal) = T(kOrdinal) + trans(i,j,0)
          T1(kOrdinal) = T1(kOrdinal) + trans1(i,j,0)
          T2(kOrdinal) = T2(kOrdinal) + trans2(i,j,0)
      endif
      !
  end do
  end do
  !
  do i=1,allkmax
      T(i) = psum(T(i))
      T1(i) = psum(T1(i))
      T2(i) = psum(T2(i))
  enddo
  !
  if(lio) then
    print *, "Calculation finish!"
    infilename='pp/Transfer'//stepname//'.dat'
    call listinit(filename=infilename,handle=hand_f, &
                      firstline='ns ti k T T1 T2')
    do i=1,allkmax
      call listwrite(hand_f,real(i,8),T(i),T1(i),T2(i))
    enddo
    print *, ' << ',trim(infilename),' ... done'
  endif
  !
end subroutine transfer2d
!
subroutine spectra2d(thefilenumb)
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  use solution
  use utility, only: listinit, listwrite
  !
  implicit none
  integer, intent(in)::thefilenumb
  character(len=128) :: infilename
  character(len=4) :: stepname
  integer :: kOrdinal
  real(8) :: kk, dk
  real(8) :: k2Edspe, kM1Edspe
  integer :: hand_f, i,j
  complex(8) :: u1s,u2s,u1d,u2d
  complex(8) ::  usspe,udspe
  !
  call quantity_prepare
  !
  ! read field
  write(stepname,'(i4.4)')thefilenumb
  infilename='outdat/flowfield'//stepname//'.h5'
  !
  call h5io_init(filename=infilename,mode='read')
  call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = 'h')
  call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = 'h')
  call h5read(varname='time',var=time)
  call h5read(varname='nstep', var=nstep)
  call h5io_end
  if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
  !
  do j=1,jm
  do i=1,im
  u1spe(i,j,1)=CMPLX(u1(i,j,0),0.d0,C_INTPTR_T);
  u2spe(i,j,1)=CMPLX(u2(i,j,0),0.d0,C_INTPTR_T);
  end do
  end do
  !
  call fft2d(u1spe)
  call fft2d(u2spe)
  !
  dk = 1.d0
  Ed = 0.d0
  Es = 0.d0
  k2Edspe = 0.d0
  kM1Edspe = 0.d0
  !
  do j=1,jm
  do i=1,im
    kk=dsqrt(k1(i,j,0)**2+k2(i,j,0)**2)
    if(kk>dk/2)then
      usspe = u1spe(i,j,1)*k2(i,j,0)/kk - u2spe(i,j,1)*k1(i,j,0)/kk
      udspe = u1spe(i,j,1)*k1(i,j,0)/kk + u2spe(i,j,1)*k2(i,j,0)/kk
      u1d =  udspe*k1(i,j,0)/kk
      u2d =  udspe*k2(i,j,0)/kk
      u1s =  usspe*k2(i,j,0)/kk 
      u2s = -usspe*k1(i,j,0)/kk
      kM1Edspe = kM1Edspe + (udspe*dconjg(udspe)) / kk
    else
      usspe = 0
      udspe = 0
      u1d = 0
      u2d = 0
      u1s = 0
      u2s = 0
    endif
    !
    kOrdinal = kint(kk,dk,2,1.d0)
    !
    if(kOrdinal <= allkmax)then
      Es(kOrdinal) = Es(kOrdinal) + usspe*conjg(usspe)
      Ed(kOrdinal) = Ed(kOrdinal) + udspe*conjg(udspe)
    endif
    Esspe = Esspe + (usspe*dconjg(usspe))
    Edspe = Edspe + (udspe*dconjg(udspe))
    k2Edspe = k2Edspe + (udspe*dconjg(udspe)) * (kk ** 2)
    !
  end do
  end do
  !
  !
  do i=1,allkmax
    Es(i) = psum(Es(i))
    Ed(i) = psum(Ed(i))
  enddo
  !
  Edspe = psum(Edspe)
  Esspe = psum(Esspe)
  k2Edspe = psum(k2Edspe)
  kM1Edspe = psum(kM1Edspe)
  !
  if(lio) then
    print *, "Calculation finish!"
    infilename='pp/Spectra'//stepname//'.dat'
    call listinit(filename=infilename,handle=hand_f, &
                      firstline='ns ti k Es Ed')
    do i=1,allkmax
      call listwrite(hand_f,real(i,8),Es(i),Ed(i))
    enddo
    print *, ' << ',trim(infilename),' ... done'
    !
    infilename='pp/Spectra_aux'//stepname//'.dat'
    call listinit(filename=infilename,handle=hand_f, &
                      firstline='ns ti Esspe Edspe k2Edspe kM1Edspe')
    call listwrite(hand_f,Esspe,Edspe,k2Edspe,kM1Edspe)
    print *, ' << ',trim(infilename),' ... done'
  endif
  !
end subroutine spectra2d
!
subroutine create_initial_field_2d
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  use solution
  !
  implicit none
  character(len=128) :: infilename
  real(8), parameter :: PI = 3.14159265358979323846d0
  integer :: hand_f, i,j,nseed
  integer, allocatable :: seed(:)
  real(8) :: var1,var2
  real(8) :: rand1,rand2
  real(8) :: kx,ky,kk,energy,factor
  !
  call quantity_prepare
  !
  !
  do j=1,jm
  do i=1,im
  u1spe(i,j,1)=CMPLX(0.d0,0.d0,C_INTPTR_T);
  u2spe(i,j,1)=CMPLX(0.d0,0.d0,C_INTPTR_T);
  end do
  end do
  !
  call random_seed(size=nseed)
  allocate(seed(1:nseed))
  call system_clock(count=seed(1))
  do i = 2, nseed
      seed(i) = seed(i-1) * 1103515245 + 12345
  end do
  !
  call random_seed(put=seed)
  !
  deallocate(seed)
  !
  call random_number(random_angle)
  !
  random_complex(:,:,1) = CMPLX(random_angle(:,:,0),0.d0,C_INTPTR_T)
  !
  call fft2d(random_complex)
  !
  do j=1,jm
  do i=1,im
    kx = k1(i,j,0)
    ky = k2(i,j,0)
    kk=dsqrt(kx**2+ky**2)
    if(kx==0 .and. ky==0) then
      u1spe(i,j,1)=0.d0
      u2spe(i,j,1)=0.d0
    else
      ! ran1: random number distributied in (0,1)
      !
      var1=kk**4*exp(-2.d0*(kk/forcek)**2)
      var2=sqrt(var1/2.d0/PI/kk)
      !
      random_angle(i,j,0) = atan2(aimag(random_complex(i,j,1)),dreal(random_complex(i,j,1))) 
      u1spe(i,j,1) = var2*kx/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
      u2spe(i,j,1) = var2*ky/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
    end if
  enddo
  enddo
  !
  call ifft2d(u1spe)
  call ifft2d(u2spe)
  !
  u1(:,:,0) = dreal(u1spe(:,:,1))
  u2(:,:,0) = dreal(u2spe(:,:,1))
  !
  energy = 0.d0
  do j=1,jm
  do i=1,im
      energy = energy + (u1(i,j,0)**2 + u2(i,j,0)**2)
  end do
  end do
  energy = psum(energy)/(ia*ja)
  factor = dsqrt(target_energy/energy)
  ! scale the field to the target energy
  u1(:,:,0) = u1(:,:,0) * factor
  u2(:,:,0) = u2(:,:,0) * factor
  !
  infilename='datin/flowini2d.h5'
  !
  call h5io_init(trim(infilename),mode='write')
  call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode='h')
  call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode='h')
  call h5write(varname='random_angle',var=random_angle(1:im,1:jm,0:km),mode='h')
  call h5io_end
  !
  if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
  !
end subroutine create_initial_field_2d
!
subroutine create_initial_field_2d_bis
  ! This imposes the uiuitheta
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  use solution
  !
  implicit none
  character(len=128) :: infilename
  real(8), parameter :: PI = 3.14159265358979323846d0
  integer :: hand_f, i,j,nseed
  integer, allocatable :: seed(:)
  real(8) :: var1,var2
  real(8) :: rand1,rand2
  real(8) :: kx,ky,kk,energy,factor
  real(8), allocatable, dimension(:,:,:) :: thetaA, thetaB,u1A,u1B,u2A,u2B
  type(C_PTR) :: c_thetaspe
  complex(c_double_complex), pointer, dimension(:,:,:):: thetaspe
  real(8) :: QX,QY,QZ,QW,QP,QQ,QDELTA,temp1,temp2,alpha,Qvalue
  call quantity_prepare
  !
  call allocate_fftw_complex(thetaspe, c_thetaspe)
  allocate(u1A(1:im,1:jm,1:1),u1B(1:im,1:jm,1:1),u2A(1:im,1:jm,1:1),u2B(1:im,1:jm,1:1))
  allocate(thetaA(1:im,1:jm,1:1),thetaB(1:im,1:jm,1:1))
  !
  !
  call random_seed(size=nseed)
  allocate(seed(1:nseed))
  call system_clock(count=seed(1))
  do i = 2, nseed
      seed(i) = seed(i-1) * 1103515245 + 12345
  end do
  !
  call random_seed(put=seed)
  !
  deallocate(seed)
  !
  ! First field generation
  !
  call random_number(random_angle)
  !
  random_complex(:,:,1) = CMPLX(random_angle(:,:,0),0.d0,C_INTPTR_T)
  !
  call fft2d(random_complex)
  !
  do j=1,jm
  do i=1,im
    kx = k1(i,j,0)
    ky = k2(i,j,0)
    kk=dsqrt(kx**2+ky**2)
    if(kx==0 .and. ky==0) then
      u1spe(i,j,1)=0.d0
      u2spe(i,j,1)=0.d0
      thetaspe(i,j,1)=0.d0
    else
      ! ran1: random number distributied in (0,1)
      !
      var1=kk**4*exp(-2.d0*(kk/forcek)**2)
      var2=sqrt(var1/2.d0/PI/kk)
      !
      random_angle(i,j,0) = atan2(aimag(random_complex(i,j,1)),dreal(random_complex(i,j,1))) 
      u1spe(i,j,1) = var2*kx/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
      u2spe(i,j,1) = var2*ky/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
      thetaspe(i,j,1) = CMPLX(0.d0,1.d0,C_INTPTR_T) * (kx * u1spe(i,j,1) + ky * u2spe(i,j,1))
    end if
  enddo
  enddo
  !
  call ifft2d(u1spe)
  call ifft2d(u2spe)
  call ifft2d(thetaspe)
  !
  u1A(:,:,1) = dreal(u1spe(:,:,1))
  u2A(:,:,1) = dreal(u2spe(:,:,1))
  thetaA(:,:,1)= dreal(thetaspe(:,:,1))
  !
  !
  ! Second field generation
  !
  call random_number(random_angle)
  !
  random_complex(:,:,1) = CMPLX(random_angle(:,:,0),0.d0,C_INTPTR_T)
  !
  call fft2d(random_complex)
  !
  do j=1,jm
  do i=1,im
    kx = k1(i,j,0)
    ky = k2(i,j,0)
    kk=dsqrt(kx**2+ky**2)
    if(kx==0 .and. ky==0) then
      u1spe(i,j,1)=0.d0
      u2spe(i,j,1)=0.d0
      thetaspe(i,j,1)=0.d0
    else
      ! ran1: random number distributied in (0,1)
      !
      var1=kk**4*exp(-2.d0*(kk/forcek)**2)
      var2=sqrt(var1/2.d0/PI/kk)
      !
      random_angle(i,j,0) = atan2(aimag(random_complex(i,j,1)),dreal(random_complex(i,j,1))) 
      u1spe(i,j,1) = var2*kx/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
      u2spe(i,j,1) = var2*ky/kk * CMPLX(sin(random_angle(i,j,0)),cos(random_angle(i,j,0)),C_INTPTR_T)
      thetaspe(i,j,1) = CMPLX(0.d0,1.d0,C_INTPTR_T) * (kx * u1spe(i,j,1) + ky * u2spe(i,j,1))
    end if
  enddo
  enddo
  !
  call ifft2d(u1spe)
  call ifft2d(u2spe)
  call ifft2d(thetaspe)
  !
  u1B(:,:,1) = dreal(u1spe(:,:,1))
  u2B(:,:,1) = dreal(u2spe(:,:,1))
  thetaB(:,:,1)= dreal(thetaspe(:,:,1))
  !
  ! Factor solution
  QX = 0.d0
  QY = 0.d0
  QZ = 0.d0
  QW = 0.d0
  !
  do j=1,jm
  do i=1,im
    QX = QX + (u1A(i,j,1) * u1A(i,j,1) + u2A(i,j,1) * u2A(i,j,1)) * thetaA(i,j,1)
    QY = QY + (u1B(i,j,1) * u1A(i,j,1) + u2B(i,j,1) * u2A(i,j,1)) * thetaA(i,j,1) + &
              (u1A(i,j,1) * u1B(i,j,1) + u2A(i,j,1) * u2B(i,j,1)) * thetaA(i,j,1) + &
              (u1A(i,j,1) * u1A(i,j,1) + u2A(i,j,1) * u2A(i,j,1)) * thetaB(i,j,1)
    QZ = QZ + (u1B(i,j,1) * u1B(i,j,1) + u2B(i,j,1) * u2B(i,j,1)) * thetaA(i,j,1) + &
              (u1A(i,j,1) * u1B(i,j,1) + u2A(i,j,1) * u2B(i,j,1)) * thetaB(i,j,1) + &
              (u1B(i,j,1) * u1A(i,j,1) + u2B(i,j,1) * u2A(i,j,1)) * thetaB(i,j,1)
    QW = QW + (u1B(i,j,1) * u1B(i,j,1) + u2B(i,j,1) * u2B(i,j,1)) * thetaB(i,j,1)
  enddo
  enddo
  QX = psum(QX)/(1.d0*ia*ja)
  QY = psum(QY)/(1.d0*ia*ja)
  QZ = psum(QZ)/(1.d0*ia*ja)
  QW = psum(QW)/(1.d0*ia*ja)
  !
  if(mpirank==0) print *, 'QX', QX, 'QY', QY, 'QZ', QZ, 'QW', QW
  !
  ! Solution of 3rd order equation
  ! Depressed cubic coefficients
  QP = QY/QW - (1.d0/3.d0)*(QZ/QW)**2
  QQ = 2.d0/27.d0*(QZ/QW)**3 - (1.d0/3.d0)*(QZ*QY)/QW**2 + QX/QW
  !
  ! Discriminant
  QDELTA = (QQ/2.d0)**2 + (QP/3.d0)**3
  !
  ! Solve cubic
  if (QDELTA >= 0.d0) then
    ! Cardano formula (real cube roots)
    temp1 = -QQ/2.d0 + SQRT(QDELTA)
    temp2 = -QQ/2.d0 - SQRT(QDELTA)
    alpha = SIGN(ABS(temp1)**(1.d0/3.d0), temp1) + &
            SIGN(ABS(temp2)**(1.d0/3.d0), temp2) - QZ/(3.d0*QW)
  else
    ! Trigonometric formula
    alpha = 2.d0*SQRT(-QP/3.d0) * COS( (1.d0/3.d0) * ACOS(3.d0*QQ/(2.d0*QP) * SQRT(-3.d0/QP)) ) - &
            QZ/(3.d0*QW)
  end if

  ! Check the solution
  Qvalue = QX + alpha*QY + alpha**2*QZ + alpha**3*QW

  ! Print results
  if(mpirank==0) print *, "alpha = ", alpha, "Check equation (=0)", Qvalue
  !
  call mpi_barrier(mpi_comm_world,ierr)
  !
  ! Summation
  u1(:,:,0) = u1A(:,:,1) + alpha * u1B(:,:,1)
  u2(:,:,0) = u2A(:,:,1) + alpha * u2B(:,:,1)
  !
  ! Normalization
  energy = 0.d0
  do j=1,jm
  do i=1,im
      energy = energy + (u1(i,j,0)**2 + u2(i,j,0)**2)
  end do
  end do
  energy = psum(energy)/(ia*ja)
  factor = dsqrt(target_energy/energy)
  !
  ! scale the field to the target energy
  u1(:,:,0) = u1(:,:,0) * factor
  u2(:,:,0) = u2(:,:,0) * factor
  !
  infilename='datin/flowini2d.h5'
  !
  call h5io_init(trim(infilename),mode='write')
  call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode='h')
  call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode='h')
  call h5write(varname='random_angle',var=random_angle(1:im,1:jm,0:km),mode='h')
  call h5io_end
  !
  if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
end subroutine create_initial_field_2d_bis
!
subroutine scale_initial_field_2d(scale)
  use stlaio,  only: get_unit
  use hdf5io
  use commvar
  use parallel,  only : mpirank,bcast
  !
  real(8), intent(in) :: scale
  character(len=256) :: filename      ! Name of the input HDF5 file
  integer :: ierr                      ! Error code for MPI operations
  character(len=1) :: modeio           ! Mode for HDF5 read operations
  integer :: n,fh,i
  ! 
  !
  ! read field
  modeio ='h'
  filename='datin/flowini2d.h5'
  !
  ! Initial field read
  !
  call h5io_init(filename=filename,mode='read')
  call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = modeio)
  call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = modeio)
  call h5io_end
  call mpi_barrier(mpi_comm_world,ierr)
  if(mpirank==0) print *, ' << ',trim(filename),' ... done'
  !
  do i = 1,im
    do j = 1,jm
        u1(i,j,0) = scale * u1(i,j,0)
        u2(i,j,0) = scale * u2(i,j,0)
    end do
  end do
  !
  filename='datin/flowini2d.'//modeio//'5'
  !
  call h5io_init(trim(filename),mode='write')
  call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode=modeio)
  call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode=modeio)
  call h5io_end
  !
  if(mpirank==0) print *, ' >> ',trim(filename),' ... done'
end subroutine scale_initial_field_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the pp program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!