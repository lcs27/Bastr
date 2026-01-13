!+---------------------------------------------------------------------+
!| This module contains subroutines used to help the development of FFT with FFTW      |
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-Mar-2024  | Created by C.S. Luo @ Beihang University             |
!| 02-Dec-2024  | Add for forcing
!+---------------------------------------------------------------------+
module fftwlink
    !
    use, intrinsic :: iso_c_binding
    use commvar,   only : im,jm,km,ia,ja,ka, ndims
    use parallel,  only : isize,jsize,ksize,irkm,jrkm,krkm,mpirank,bcast,mpirankmax
    !
    implicit none
    !
    !
    integer(C_INTPTR_T) :: alloc_local,iafftw,jafftw,kafftw,imfftw,jmfftw,kmfftw
    integer :: nproc, myid, ierr
    integer :: i0f,j0f,k0f
    ! Reserved for forcing part
    integer :: rooti0, imf, jmf, kmf
    integer, allocatable, dimension(:) :: counts_grid,disp_grid,counts_fence,disp_fence
    type(C_PTR) :: forward_plan, backward_plan
    !
    contains
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This subroutine is used to assign the distributions size on the ranks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Writen by Chensheng Luo, 2024-03-06.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    subroutine mpisizedis_fftw
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        ! 
        iafftw = ia
        jafftw = ja
        kafftw = ka
        !
        !
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        isize = 1
        jsize = 1
        ksize = 1
        !
        if(ndims==3) then
            ! 
            alloc_local = fftw_mpi_local_size_3d(kafftw,jafftw,iafftw,MPI_COMM_WORLD, kmfftw,k0fftw)
            jmfftw = jafftw
            imfftw = iafftw
            krk = myid
            ksize = nproc
            !
        elseif(ndims==2) then 
            alloc_local = fftw_mpi_local_size_2d(jafftw,iafftw,MPI_COMM_WORLD, jmfftw, j0fftw)
            imfftw = iafftw
            kmfftw = 0
            jrk = myid
            jsize = nproc
            !
        else
            stop 'fftw_mpi: 1D not implemented!'
            !alloc_local = fftw_mpi_local_size(myid,iafftw,MPI_COMM_WORLD, imfftw, i0fftw)
            jmfftw = 0
            kmfftw = 0
            irk = myid
            isize = nproc
            !
        endif
        !
        if(mpirank == 0)then
            write(*,'(3(A,I0))')'  ** mpi size= ',isize,' x ',jsize,' x ',ksize
        endif
        !
        imf = imfftw
        jmf = jmfftw
        kmf = kmfftw
        im = imf
        jm = jmf
        km = kmf
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        allocate(ids(0:(nproc-1)),irks(0:(nproc-1)),jrks(0:(nproc-1)),krks(0:(nproc-1)),&
                ims(0:(nproc-1)),jms(0:(nproc-1)),kms(0:(nproc-1)),i0s(0:(nproc-1)),    &
                j0s(0:(nproc-1)),k0s(0:(nproc-1)))
        !
        do i=0,(nproc-1)
            if(myid == i) then
                ids(i) = myid
                irks(i) = irk
                jrks(i) = jrk
                krks(i) = krk
                ims(i) = im 
                jms(i) = jm
                kms(i) = km 
                i0s(i) = i0fftw
                j0s(i) = j0fftw
                k0s(i) = k0fftw
            endif
            call mpi_bcast(ids(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(irks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jrks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(krks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(ims(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(kms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(i0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(j0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(k0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
        enddo
        !
        if(mpirank == 0)then
            open(fh,file='datin/parallel.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")isize,jsize,ksize
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,(nproc-1)
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print*,' << parallel.info ... done !'
        endif
        !
        irkm=isize-1
        jrkm=jsize-1
        krkm=ksize-1
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(mpirank == 0) print*,' ** parallel processing ... done.'
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine mpisizedis_fftw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The end of the subroutine mpisizedis_fftw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    subroutine mpisizedis_half_fftw
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        ! 
        iafftw = ia
        jafftw = ja
        kafftw = ka
        !
        !
        call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        !
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        isize = 1
        jsize = 1
        ksize = 1
        !
        if(ndims==3) then
            ! 
            alloc_local = fftw_mpi_local_size_3d(kafftw,jafftw,iafftw/2+1,MPI_COMM_WORLD, kmfftw, k0fftw)
            k0fftw = k0fftw
            jmfftw = jafftw
            imfftw = iafftw/2+1
            krk = myid
            ksize = nproc
            !
        elseif(ndims==2) then 
            alloc_local = fftw_mpi_local_size_2d(jafftw,iafftw/2+1,MPI_COMM_WORLD, jmfftw, j0fftw)
            j0fftw = j0fftw
            imfftw = iafftw/2+1
            kmfftw = 0
            jrk = myid
            jsize = nproc
            !
        else
            stop 'fftw_mpi: 1D not implemented!'
            !alloc_local = fftw_mpi_local_size(myid,iafftw,MPI_COMM_WORLD, imfftw, i0fftw)
            jmfftw = 0
            kmfftw = 0
            irk = myid
            isize = nproc
            !
        endif
        !
        if(mpirank == 0)then
            write(*,'(3(A,I0))')'  ** mpi size= ',isize,' x ',jsize,' x ',ksize
        endif
        !
        imf = imfftw
        jmf = jmfftw
        kmf = kmfftw
        im = imf
        jm = jmf
        km = kmf
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        allocate(ids(0:(nproc-1)),irks(0:(nproc-1)),jrks(0:(nproc-1)),krks(0:(nproc-1)),&
                ims(0:(nproc-1)),jms(0:(nproc-1)),kms(0:(nproc-1)),i0s(0:(nproc-1)),    &
                j0s(0:(nproc-1)),k0s(0:(nproc-1)))
        !
        do i=0,(nproc-1)
            if(myid == i) then
                ids(i) = myid
                irks(i) = irk
                jrks(i) = jrk
                krks(i) = krk
                ims(i) = im
                jms(i) = jm
                kms(i) = km
                i0s(i) = i0fftw
                j0s(i) = j0fftw
                k0s(i) = k0fftw
            endif
            call mpi_bcast(ids(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(irks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jrks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(krks(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(ims(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(jms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(kms(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(i0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(j0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
            call mpi_bcast(k0s(i),1,mpi_integer,i,mpi_comm_world,ierr)
        enddo
        !
        if(mpirank == 0)then
            open(fh,file='datin/parallel.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")isize,jsize,ksize
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,(nproc-1)
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print*,' << parallel.info ... done !'
        endif
        !
        irkm=isize-1
        jrkm=jsize-1
        krkm=ksize-1
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(mpirank == 0) print*,' ** parallel processing ... done.'
        !
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine mpisizedis_half_fftw
    ! !
    !
    !
    subroutine prepareplan_fftw_plan
        !
        use commvar, only: u1spe
        use parallel
        !
        include 'fftw3-mpi.f03'
        !
        !
        select case(ndims)
        case(2)
            forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
            backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        case(3)
            forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
            backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        case default
            stop 'prepareplan_fftw: Not implemented for 1D!'
        end select
        !
    end subroutine prepareplan_fftw_plan
    !
    subroutine allocate_fftw_complex(array, c_pointer)
        use, intrinsic :: iso_c_binding
        include 'fftw3-mpi.f03'
        !
        type(C_PTR), intent(out) :: c_pointer
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: array(:,:,:)
        !
        c_pointer = fftw_alloc_complex(alloc_local)
        select case(ndims)
        case(2)
            call c_f_pointer(c_pointer, array, [imfftw, jmfftw, int(1,C_INTPTR_T)])
        case(3)
            call c_f_pointer(c_pointer, array, [imfftw, jmfftw, kmfftw])
        case default
            stop 'allocate_fftw_complex: Not implemented for 1D!'
        end select
        !
    end subroutine allocate_fftw_complex
  !
    subroutine fft2d(array)
        !
        use commvar, only: im,jm,ia,ja
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:,:)
        integer :: i,j
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(forward_plan,array,array)
        do j=1,jm
        do i=1,im
            array(i,j,1)=array(i,j,1)/(1.d0*ia*ja)
        end do
        end do
    end subroutine fft2d
    !
    subroutine fft3d(array)
        !
        use commvar, only: im,jm,km,ia,ja,ka
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:,:)
        integer :: i,j,k
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(forward_plan,array,array)
        do k=1,km
        do j=1,jm
        do i=1,im
            array(i,j,k)=array(i,j,k)/(1.d0*ia*ja*ka)
        end do
        end do
        end do
    end subroutine fft3d
    !
     subroutine ifft2d(array)
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:,:)
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(backward_plan,array,array)
        !
    end subroutine ifft2d
        !
     subroutine ifft3d(array)
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:,:)
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(backward_plan,array,array)
        !
    end subroutine ifft3d
    !
    subroutine allocation
        !
        use commvar
        include 'fftw3-mpi.f03'
        ! This subroutine allocates the common variables.
        !
        select case(ndims)
        case(2)
            call allocation2D
        case(3)
            call allocation3D
        case default
            stop 'allocation: Not implemented for 1D!'
        end select
        !
        allocate(Es(0:allkmax),Ed(0:allkmax),kn(0:allkmax),Ecount(0:allkmax))
        !
    end subroutine allocation
    !
    subroutine allocation2D
        !
        use commvar
        include 'fftw3-mpi.f03'
        allocate(u1(1:im,1:jm,0:km), u2(1:im,1:jm,0:km))
        allocate(k1(1:im,1:jm,0:km), k2(1:im,1:jm,0:km))
        !
        call allocate_fftw_complex(u1spe, c_u1spe)
        call allocate_fftw_complex(u2spe, c_u2spe)
        call allocate_fftw_complex(u1x1,   c_u1x1)
        call allocate_fftw_complex(u1x2,   c_u1x2)
        call allocate_fftw_complex(u2x1,   c_u2x1)
        call allocate_fftw_complex(u2x2,   c_u2x2)
        call allocate_fftw_complex(u1xixi, c_u1xixi)
        call allocate_fftw_complex(u2xixi, c_u2xixi)
        call allocate_fftw_complex(thetaxixi, c_thetaxixi)
        call allocate_fftw_complex(u1tA,     c_u1tA)
        call allocate_fftw_complex(u2tA,     c_u2tA)
        call allocate_fftw_complex(force1,   c_force1)
        call allocate_fftw_complex(force2,   c_force2)
        if(timemethod==1)then
            call allocate_fftw_complex(u1tB,     c_u1tB)
            call allocate_fftw_complex(u2tB,     c_u2tB)
            call allocate_fftw_complex(u1tC,     c_u1tC)
            call allocate_fftw_complex(u2tC,     c_u2tC)
        elseif(timemethod==2)then
            allocate(u1old(1:im,1:jm,0:km), u2old(1:im,1:jm,0:km))
        else
            stop 'allocation2D: timemethod not recognized!'
        endif
        !
    end subroutine allocation2D
    !
    subroutine allocation3D
        use commvar
        include 'fftw3-mpi.f03'
        allocate(u1(1:im,1:jm,1:km), u2(1:im,1:jm,1:km), u3(1:im,1:jm,1:km))
        allocate(k1(1:im,1:jm,1:km), k2(1:im,1:jm,1:km), k3(1:im,1:jm,1:km))
        !
        call allocate_fftw_complex(u1spe, c_u1spe)
        call allocate_fftw_complex(u2spe, c_u2spe)
        call allocate_fftw_complex(u3spe, c_u3spe)
        call allocate_fftw_complex(u1x1,   c_u1x1)
        call allocate_fftw_complex(u1x2,   c_u1x2)
        call allocate_fftw_complex(u1x3,   c_u1x3)
        call allocate_fftw_complex(u2x1,   c_u2x1)
        call allocate_fftw_complex(u2x2,   c_u2x2)
        call allocate_fftw_complex(u2x3,   c_u2x3)
        call allocate_fftw_complex(u3x1,   c_u3x1)
        call allocate_fftw_complex(u3x2,   c_u3x2)
        call allocate_fftw_complex(u3x3,   c_u3x3)
        call allocate_fftw_complex(u1xixi, c_u1xixi)
        call allocate_fftw_complex(u2xixi, c_u2xixi)
        call allocate_fftw_complex(u3xixi, c_u3xixi)
        call allocate_fftw_complex(thetaxixi, c_thetaxixi)
        call allocate_fftw_complex(u1tA,     c_u1tA)
        call allocate_fftw_complex(u2tA,     c_u2tA)
        call allocate_fftw_complex(u3tA,     c_u3tA)
        call allocate_fftw_complex(u1tB,     c_u1tB)
        call allocate_fftw_complex(u2tB,     c_u2tB)
        call allocate_fftw_complex(u3tB,     c_u3tB)
        call allocate_fftw_complex(u1tC,     c_u1tC)
        call allocate_fftw_complex(u2tC,     c_u2tC)
        call allocate_fftw_complex(u3tC,     c_u3tC)
        call allocate_fftw_complex(force1,   c_force1)
        call allocate_fftw_complex(force2,   c_force2)
        call allocate_fftw_complex(force3,   c_force3)
        !
    end subroutine allocation3D
    !
    subroutine deallocation
        !
        use commvar
        use parallel, only: mpistop
        include 'fftw3-mpi.f03'
        !
        deallocate(Es, Ed, kn,Ecount)
        !
        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        call fftw_mpi_cleanup()
        
        select case(ndims)
        case(3)
            call deallocation3D
        case(2)
            call deallocation2D
        end select
        !
        call mpistop
        !
    end subroutine deallocation
    !
    subroutine deallocation3D
        !
        use commvar
        include 'fftw3-mpi.f03'
        !
        deallocate(u1,u2,u3,k1,k2,k3)
        call fftw_free(c_u1spe)
        call fftw_free(c_u2spe)
        call fftw_free(c_u3spe)
        call fftw_free(c_u1x1)
        call fftw_free(c_u1x2)
        call fftw_free(c_u1x3)
        call fftw_free(c_u2x1)
        call fftw_free(c_u2x2)
        call fftw_free(c_u2x3)
        call fftw_free(c_u3x1)
        call fftw_free(c_u3x2)
        call fftw_free(c_u3x3)
        call fftw_free(c_u1xixi)
        call fftw_free(c_u2xixi)
        call fftw_free(c_u3xixi)
        call fftw_free(c_thetaxixi)
        call fftw_free(c_u1tA)
        call fftw_free(c_u2tA)
        call fftw_free(c_u3tA)
        call fftw_free(c_u1tB)
        call fftw_free(c_u2tB)
        call fftw_free(c_u3tB)
        call fftw_free(c_u1tC)
        call fftw_free(c_u2tC)
        call fftw_free(c_u3tC)
        call fftw_free(c_force1)
        call fftw_free(c_force2)
        call fftw_free(c_force3)
    end subroutine deallocation3D
    !
    subroutine deallocation2D
        !
        use commvar
        include 'fftw3-mpi.f03'
        !
        deallocate(u1,u2,k1,k2)
        call fftw_free(c_u1spe)
        call fftw_free(c_u2spe)
        call fftw_free(c_u1x1)
        call fftw_free(c_u1x2)
        call fftw_free(c_u2x1)
        call fftw_free(c_u2x2)
        call fftw_free(c_u1xixi)
        call fftw_free(c_u2xixi)
        call fftw_free(c_thetaxixi)
        call fftw_free(c_u1tA)
        call fftw_free(c_u2tA)
        call fftw_free(c_u1tB)
        call fftw_free(c_u2tB)
        call fftw_free(c_u1tC)
        call fftw_free(c_u2tC)
        call fftw_free(c_force1)
        call fftw_free(c_force2)
    end subroutine deallocation2D
end module fftwlink