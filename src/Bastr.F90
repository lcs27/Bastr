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
  integer :: time_values(8),seed(8)
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
        call date_and_time(values=time_values)
        do i = 1, 8
            seed(i) = time_values(1) * 10000000 + &   ! 年
                    time_values(2) * 100000 + &     ! 月
                    time_values(3) * 1000 + &       ! 日
                    ABS(time_values(4)) * 100 + &   ! 时差
                    time_values(5) * 6000000 + &    ! 时
                    time_values(6) * 60000 + &      ! 分
                    time_values(7) * 600 + &        ! 秒
                    time_values(8) * 60  ! 毫秒
            seed(i) = seed(i) + i * 314159
            !
            if (seed(i) < 0) seed(i) = -seed(i)
            !
            seed(i) = mod(seed(i), 1000000000)
        end do
        !
        call random_seed(put=seed)
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
