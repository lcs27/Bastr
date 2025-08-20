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
  include 'fftw3-mpi.f03'
  !
  integer :: hand_f, hand_g, hand_a, hand_fo, i,j
  !
  call mpiinitial
  !
  call readinput
  !
  ndims = 2
  !
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
  call add_field
  !
  if(mpirank==0)  print *, "** End!"
  !
  !
end program Bastrpp

subroutine add_field
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
end subroutine add_field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the pp program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!