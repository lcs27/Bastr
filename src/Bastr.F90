!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for Viscous Burgers equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2025-07-09 by Chensheng Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Bastr
  !
  ! TODO: To test adaptation for 3D
  !       - 3D to be tested
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
  integer :: hand_f, hand_g, hand_a, hand_fo, nseed
  integer, allocatable :: seed(:)
  character(len=16) :: cmd
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
  call quantity_prepare
  !
  if(mpirank == 0) then
    call readkeyboad(cmd)
  endif
  !
  call bcast(cmd)
  !
  if (len_trim(cmd) == 0) then
    lrestart = .false.
    if(mpirank == 0) print*,' start from zero.... '
    call read_initial_field
  else
    lrestart = .true.
    read(cmd,'(i4)') filenumb
    if(mpirank == 0) print*,'continue to run: start from ',filenumb
    call read_continue_field
  endif
  !
  ! Initial logs
  if(mpirank==0) then
    select case(ndims)
    case(2)
      call listinit(filename='log/stat2d_specall.dat',handle=hand_f, &
                      firstline='ns ti Es Ed k2Ed kMEd')
      call listinit(filename='log/stat2d_spect.dat',handle=hand_g, &
                      firstline='ns ti k Es Ed')
      call listinit(filename='log/stat2d_velgrad.dat',handle=hand_a, &
                      firstline='ns ti umumtheta2 umumijji u2theta dissp etamin Tay ReTay Kol')
    case(3)
      call listinit(filename='log/stat3d_specall.dat',handle=hand_f, &
                      firstline='ns ti Es Ed k2Ed kMEd')
      call listinit(filename='log/stat3d_spect.dat',handle=hand_g, &
                      firstline='ns ti k Es Ed')
      call listinit(filename='log/stat3d_velgrad.dat',handle=hand_a, &
                      firstline='ns ti umumtheta2 umumijji u2theta dissp etamin Tay ReTay Kol')
    end select
    if(forcemethod>0) then
      call listinit(filename='log/factor_output.dat',handle=hand_fo, &
                      firstline='ns ti factor')
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
    end if
  endif
  !
  ! Solution !
  !
  !
  call mainloop(hand_f, hand_g, hand_a, hand_fo)
  !
  call deallocation
  !
  if(mpirank==0)  print *, "** End!"
  !
  !
end program Bastr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
