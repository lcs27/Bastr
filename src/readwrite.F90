!+---------------------------------------------------------------------+
!| This module contains subroutines of reading and writing files.      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Oct-2018  | Created by J. Fang @ Warrington                      |
!+---------------------------------------------------------------------+
module readwrite
  !
  use commvar
  use parallel, only: mpirank,mpirankname,mpistop,lio,irk,jrkm,jrk,    &
                       ptime,bcast,mpirankmax,ig0,jg0,kg0
  use hdf5io
  !
  implicit none
  !
  contains
  !
  subroutine readinput
    use stlaio,  only: get_unit
    !
    logical :: lfex
    character(len=64) :: inputfile
    character(len=5) :: char
    integer :: n,fh,i
    !
    if(mpirank==0) then
      !
      inputfile='datin/input.dat'
      !
      !
      fh=get_unit()
      !
      open(fh,file=trim(inputfile),action='read')
      read(fh,'(////)')
      read(fh,*)ia,ja,ka
      read(fh,'(/)')
      read(fh,*)ref_tem,reynolds
      read(fh,'(/)')
      read(fh,*)maxstep,deltat,lwsequ,feqwsequ,lwspectra,feqwspe,kmax,timemethod
      read(fh,'(/)')
      read(fh,*)forcemethod, forcek, target_energy, lprojectd
      read(fh,'(/)')
      read(fh,*)initialmethod
      print *,' << ',trim(inputfile),' ... done'
      !
      print *, 'ia:', ia, 'ja:', ja, 'ka:', ka
    endif
    !
    call bcast(ia)
    call bcast(ja)
    call bcast(ka)
    call bcast(ref_tem)
    call bcast(reynolds)
    call bcast(lwsequ)
    call bcast(maxstep)
    call bcast(feqwsequ)
    call bcast(lwspectra)
    call bcast(feqwspe)
    call bcast(deltat)
    call bcast(forcemethod)
    call bcast(target_energy)
    call bcast(lprojectd)
    call bcast(kmax)
    call bcast(timemethod)
    call bcast(forcek)
    call bcast(initialmethod)
    !
    if(ka == 0) then
      if(ja == 0 ) then
        ndims = 1
        error stop '1D case is not implemented yet!'
      else
        ndims = 2
        allkmax = int(ia / 3)
      endif
    else
      ndims = 3
      allkmax = int(ia / 3)
    endif
    !
    !
    allkmax = min(allkmax, kmax)

    if(mpirank==0) then
      !
      print *, 'ia=', ia, ' ja=', ja, ' ka=', ka
      print *, 'ndims=', ndims
      print *, 'ref_tem=', ref_tem
      print *, 'reynolds=', reynolds
      print *, 'maxstep=', maxstep
      print *, 'deltat=', deltat
      print *, 'lwsequ=', lwsequ
      print *, 'feqwsequ=', feqwsequ
      print *, 'lwspectra=', lwspectra
      print *, 'feqwspe=', feqwspe
      print *, 'kmax=', allkmax
      if(timemethod==1)then
        print *, 'timemethod= RK3'
      elseif(timemethod==2)then
        print *, 'timemethod= CN'
      else
        stop 'Undefined time method'
      endif
      print *, 'forcemethod=', forcemethod
      print *, 'forcek=', forcek
      print *, 'target_energy=', target_energy
      print *, 'lprojectd=', lprojectd
      print *, 'initialmethod=', initialmethod
    endif
    !
  end subroutine readinput
  !
  subroutine read_initial_field
    !
    use commvar, only : ndims
    character(len=256) :: infilename      ! Name of the input HDF5 file
    integer :: ierr                      ! Error code for MPI operations
    character(len=1) :: modeio           ! Mode for HDF5 read operations
    !
    if(initialmethod==0)then
      u1=0.d0
      u2=0.d0
      if(ndims==3)then
        u3=0.d0
      endif
      if(mpirank==0) print *, ' give void initial velocity field ... done'
    elseif(initialmethod==1)then
      ! read initial field
      modeio ='h'
      select case(ndims)
      case(2)
        infilename='datin/flowini2d.h5'
        call h5io_init(filename=infilename,mode='read')
        call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = modeio)
        call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = modeio)
        call h5io_end
        call mpi_barrier(mpi_comm_world,ierr)
        if(mpirank==0) print *, ' << ',trim(infilename),' ... done'
      case(3)
        infilename='datin/flowini3d.h5'
        call h5io_init(filename=infilename,mode='read')
        call h5read(varname='u1', var=u1(1:im,1:jm,1:km),mode = modeio)
        call h5read(varname='u2', var=u2(1:im,1:jm,1:km),mode = modeio)
        call h5read(varname='u3', var=u3(1:im,1:jm,1:km),mode = modeio)
        call h5io_end
        call mpi_barrier(mpi_comm_world,ierr)
        if(mpirank==0) print *, ' << ',trim(infilename),' ... done'
      case default
        error stop 'ndims should be 2 or 3!'
      end select
    endif
  !
  end subroutine read_initial_field
  !
  subroutine read_continue_field
    use commvar, only : ndims
    !
    character(len=256) :: infilename      ! Name of the input HDF5 file
    integer :: ierr                      ! Error code for MPI operations
    character(len=1) :: modeio           ! Mode for HDF5 read operations
    character(len=4) :: stepname
    !
    ! read initial field
    modeio ='h'
    write(stepname,'(i4.4)')filenumb
    infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    !
    ! Initial field read
    !
    call h5io_init(filename=infilename,mode='read')
    select case(ndims)
    case(2)
      call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = modeio)
      call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = modeio)
    case(3)
      call h5read(varname='u1',var=u1(1:im,1:jm,1:km),mode=modeio)
      call h5read(varname='u2',var=u2(1:im,1:jm,1:km),mode=modeio)
      call h5read(varname='u3',var=u3(1:im,1:jm,1:km),mode=modeio)
    case default
      error stop 'ndims should be 2 or 3!'
    end select
    !
    call h5read(varname='nstep',var=nstep)
    call h5read(varname='time',var=time)
    !
    if(mpirank==0) print *, ' >> ',trim(infilename),' ... done'
    call h5io_end
    filenumb = filenumb + 1
    nxtwsequ = min(nstep + feqwsequ, maxstep)
    nxtwspe = min(nstep + feqwspe, maxstep)
  !
  end subroutine read_continue_field
  !
  subroutine writeflfed
    use commvar, only : ndims
    !
    character(len=4) :: stepname
    character(len=64) :: outfilename
    character(len=1) :: modeio       ! Mode for HDF5 read operations
    modeio ='h'
    !
    write(stepname,'(i4.4)')filenumb
    !
    outfilename='outdat/flowfield'//stepname//'.'//modeio//'5'
    !
    ! outfilename='outdat/flowfield.h5'
    call h5io_init(trim(outfilename),mode='write')
    !
    select case(ndims)
    case(2)
      call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode=modeio)
      call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode=modeio)
    case(3)
      call h5write(varname='u1',var=u1(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='u2',var=u2(1:im,1:jm,1:km),mode=modeio)
      call h5write(varname='u3',var=u3(1:im,1:jm,1:km),mode=modeio)
    case default
      error stop 'ndims should be 2 or 3!'
    end select
    !
    call h5io_end
    if(mpirank == 0) then
      call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
      call h5srite(varname='time',var=time,filename=trim(outfilename))
    endif
    !
    if(mpirank==0) print *, ' >> ',trim(outfilename),' ... done'
    !
  end subroutine writeflfed
  !
end module readwrite
!+---------------------------------------------------------------------+
!| The end of the module readwrite.                                    |
!+---------------------------------------------------------------------+
