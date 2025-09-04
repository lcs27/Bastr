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
      read(fh,*)maxstep,deltat,lwsequ,feqwsequ,lwspectra,feqwspe,kmax
      read(fh,'(/)')
      read(fh,*)forcemethod, target_energy, lprojectd
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
    !
    if(mpirank==0) then
      !
      print *, 'ia=', ia, ' ja=', ja, ' ka=', ka
      print *, 'ref_tem=', ref_tem
      print *, 'reynolds=', reynolds
      print *, 'maxstep=', maxstep
      print *, 'deltat=', deltat
      print *, 'lwsequ=', lwsequ
      print *, 'feqwsequ=', feqwsequ
      print *, 'lwspectra=', lwspectra
      print *, 'feqwspe=', feqwspe
      print *, 'forcemethod=', forcemethod
      print *, 'target_energy=', target_energy
      print *, 'lprojectd=', lprojectd
    endif
    !
  end subroutine
  !
  subroutine read_initial_field
    !
    character(len=256) :: infilename      ! Name of the input HDF5 file
    integer :: ierr                      ! Error code for MPI operations
    character(len=1) :: modeio           ! Mode for HDF5 read operations
    ! read initial field
    modeio ='h'
    infilename='datin/flowini2d.h5'
    !
    ! Initial field read
    !
    call h5io_init(filename=infilename,mode='read')
    call h5read(varname='u1', var=u1(1:im,1:jm,0:km),mode = modeio)
    call h5read(varname='u2', var=u2(1:im,1:jm,0:km),mode = modeio)
    call h5io_end
    call mpi_barrier(mpi_comm_world,ierr)
    if(mpirank==0) print *, ' << ',trim(infilename),' ... done'
  !
  end subroutine read_initial_field
  !
  subroutine writeflfed
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
    call h5write(varname='u1',var=u1(1:im,1:jm,0:km),mode=modeio)
    call h5write(varname='u2',var=u2(1:im,1:jm,0:km),mode=modeio)
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
