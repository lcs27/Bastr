!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This Program is a DNS solver for Viscous Burgers equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2025-07-09 by Chensheng Luo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program Bastr
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
  allkmax = ceiling(real(sqrt(2.d0)/3*min(ia,ja)))
  allkmax = min(allkmax, kmax)
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
      call listinit(filename='log/stat2d_specall.dat',handle=hand_f, &
                      firstline='ns ti Es Ed k2Ed kdEd k2dEd')
      call listinit(filename='log/stat2d_spect.dat',handle=hand_g, &
                      firstline='ns ti k Es Ed')
      call listinit(filename='log/stat2d_velgrad.dat',handle=hand_a, &
                      firstline='ns ti umumtheta2 umumijji u2theta dissp etamin')
      if(lforce) then
        call listinit(filename='log/factor_output.dat',handle=hand_fo, &
                        firstline='ns ti factor')
      end if
  endif
  !
  ! Solution !
  !
  !
  do while(nstep<=maxstep)
    !
    do j=1,jm
    do i=1,im
      u1spe(i,j)=CMPLX(u1(i,j,0),0.d0,C_INTPTR_T);
      u2spe(i,j)=CMPLX(u2(i,j,0),0.d0,C_INTPTR_T);
    end do
    end do
    !
    call fft2d(u1spe)
    call fft2d(u2spe)
    !
    call spectra_compute(hand_f, hand_g)
    !
    call velgrad_calculate(hand_a)
    !
    if(mpirank==0)  print *, "nstep, Es, Ed, Eall= ", nstep, Esspe, Edspe, &
                        Esspe+Edspe
    !
    if(lwsequ .and. nstep==nxtwsequ) then
      !
      call writeflfed
      !
      filenumb = filenumb + 1
      !
      nxtwsequ = min(nxtwsequ + feqwsequ, maxstep)
      !
      if(mpirank==0)  print *, "Next print is ", nxtwsequ
      !
    endif
    !
    nstep = nstep + 1 
    time = time + deltat
    !
    call RK3(hand_fo)
    !
  end do
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
