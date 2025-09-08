!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for Viscous Burgers equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2025-07-09 by Chensheng Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Bastr
  !
  ! TODO: To test adaptation for 3D
  !       - 2D still march and gives the same result - OK
  !       - 3D to be tested
  use commvar
  use parallel
  use fftwlink
  use tool
  use solution
  use readwrite
  use utility
  include 'fftw3-mpi.f03'
  !
  integer :: hand_f, hand_g, hand_a, hand_fo
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
  call read_initial_field
  !
  call quantity_prepare
  !
  if(mpirank==0) then
    select case(ndims)
    case(2)
      call listinit(filename='log/stat2d_specall.dat',handle=hand_f, &
                      firstline='ns ti Es Ed k2Ed kdEd k2dEd')
      call listinit(filename='log/stat2d_spect.dat',handle=hand_g, &
                      firstline='ns ti k Es Ed')
      call listinit(filename='log/stat2d_velgrad.dat',handle=hand_a, &
                      firstline='ns ti umumtheta2 umumijji u2theta dissp etamin Tay ReTay Kol')
    case(3)
      call listinit(filename='log/stat3d_specall.dat',handle=hand_f, &
                      firstline='ns ti Es Ed k2Ed kdEd k2dEd')
      call listinit(filename='log/stat3d_spect.dat',handle=hand_g, &
                      firstline='ns ti k Es Ed')
      call listinit(filename='log/stat3d_velgrad.dat',handle=hand_a, &
                      firstline='ns ti umumtheta2 umumijji u2theta dissp etamin Tay ReTay Kol')
    end select
      if(forcemethod>0) then
        call listinit(filename='log/factor_output.dat',handle=hand_fo, &
                        firstline='ns ti factor')
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
