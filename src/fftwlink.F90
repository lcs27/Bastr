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
    interface fftw_grid_fence
        module procedure fftw_grid_fence_2D
        module procedure fftw_grid_fence_3D
    end interface
    !
    interface fftw_fence_grid
        module procedure fftw_fence_grid_2D
        module procedure fftw_fence_grid_3D
    end interface
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
    subroutine fftwprepare_forcing
        !
        if(ndims == 2) then
            call fftwprepare_forcing2D
        else
            stop "Not implemented error! fftwprepare_forcing3D"
        endif
    end subroutine fftwprepare_forcing
    !
    subroutine fftwprepare_forcing2D
        !
        include 'fftw3-mpi.f03'
        include 'mpif.h' 
        !
        integer :: fh,irk,jrk,krk
        integer(C_INTPTR_T) :: i0fftw,j0fftw,k0fftw
        integer :: i,ierr
        integer, allocatable, dimension(:) :: ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        rooti0 = isize * (mpirank/isize)
        jmf = jm/isize
        !
        if((jmf*ia) .ne. (im*jm))then
            print *,'jmf = ', jmf, 'ia=', ia, 'im=', im, 'jm=', jm
            stop 'fftwprepare_forcing2D: Parallel ranking not match!'
        endif
        !
        allocate(counts_grid(0:mpirankmax),disp_grid(0:mpirankmax))
        !
        counts_grid = 0
        disp_grid = 0
        counts_grid(rooti0) = im*jm
        do i = (rooti0+1), (rooti0+isize-1)
            counts_grid(i) = im*jm
            disp_grid(i) =disp_grid(i-1) + counts_grid(i-1)
        enddo
        !
        allocate(counts_fence(0:mpirankmax),disp_fence(0:mpirankmax))
        counts_fence = 0
        disp_fence = 0
        counts_fence(rooti0) = ia*jmf
        do i = (rooti0+1), (rooti0+isize-1)
            counts_fence(i) = ia*jmf
            disp_fence(i) =disp_fence(i-1) + counts_fence(i-1)
        enddo
        !
        iafftw = ia
        jafftw = ja
        kafftw = ka
        i0fftw = 0
        j0fftw = 0
        k0fftw = 0
        irk = 0
        jrk = 0
        krk = 0
        !
        alloc_local = fftw_mpi_local_size_2d(jafftw, iafftw,MPI_COMM_WORLD, jmfftw, j0fftw)
        !
        imfftw = iafftw
        kmfftw = 0
        !
        i0f = i0fftw
        j0f = j0fftw
        k0f = k0fftw
        !
        imf = ia
        kmf = 0
        if(jmfftw .ne. jmf)then
            print *, 'mpirank = ', mpirank, 'jmfftw=',jmfftw, 'jmf=', jmf
            stop "fftwprepare_forcing2D: Rank not matching!"
        endif
        !
        allocate(ids(0:mpirankmax),irks(0:mpirankmax),jrks(0:mpirankmax),krks(0:mpirankmax),&
                ims(0:mpirankmax),jms(0:mpirankmax),kms(0:mpirankmax),i0s(0:mpirankmax),    &
                j0s(0:mpirankmax),k0s(0:mpirankmax))
        !
        do i=0,mpirankmax
            if(mpirank == i) then
                ids(i) = mpirank
                irks(i) = 0
                jrks(i) = mpirank
                krks(i) = 0
                ims(i) = imf
                jms(i) = jmf
                kms(i) = kmf 
                i0s(i) = i0f
                j0s(i) = j0f
                k0s(i) = k0f
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
            open(fh,file='datin/parallelfftw.info',form='formatted')
            write(fh,"(3(A9,1x))")'isize','jsize','ksize'
            write(fh,"(3(I9,1x))")1,mpirank,1
            write(fh,"(10(A9,1x))")'Rank','Irk','Jrk','Krk','IM','JM','KM', 'I0','J0','K0'
            do i=0,mpirankmax
                write(fh,"(10(I9,1x))")ids(i),irks(i),jrks(i),krks(i),ims(i),jms(i),kms(i),i0s(i),j0s(i),k0s(i)
            enddo
            close(fh)
            print *,' << parallelfftw.info ... done !'
        endif
        !
        deallocate(ids,irks,jrks,krks,ims,jms,kms,i0s,j0s,k0s)
        !
    end subroutine fftwprepare_forcing2D
    !
    subroutine fftw_grid_fence_2D(gridarray_t,fencearray)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:)  ::  gridarray_t 
        real(8), intent(out), dimension(:,:) ::  fencearray
        real(8), allocatable, dimension(:,:) :: global_t, global
        !
        allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        !
        !
        call mpi_gatherv(gridarray_t,im*jm, MPI_REAL8, &
        global_t,counts_grid,disp_grid,MPI_REAL8,rooti0, &
        MPI_COMM_WORLD, ierr)
        !
        global = transpose(global_t)
        !
        call mpi_scatterv(global,counts_fence,disp_fence,MPI_REAL8,&
                        fencearray,ia*jmf,MPI_REAL8,rooti0,&
                            MPI_COMM_WORLD,ierr)
        !
        deallocate(global_t,global)
        !
    end subroutine fftw_grid_fence_2D
    !
    subroutine fftw_fence_grid_2D(fencearray,gridarray_t)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:)  :: fencearray
        real(8), intent(out), dimension(:,:) :: gridarray_t
        real(8), allocatable, dimension(:,:) :: global_t, global
        !
        allocate(global_t(1:jm,1:ia),global(1:ia,1:jm))
        !
        call mpi_gatherv(fencearray,ia*jmf, MPI_REAL8, &
        global,counts_fence,disp_fence,MPI_REAL8,rooti0,  &
        MPI_COMM_WORLD, ierr)
        !
        global_t = transpose(global)
        !
        call mpi_scatterv(global_t,counts_grid,disp_grid,MPI_REAL8,&
                            gridarray_t,im*jm,MPI_REAL8,rooti0,  &
                            MPI_COMM_WORLD,ierr)
        !
        deallocate(global_t,global)
        !
    end subroutine fftw_fence_grid_2D
    !
    subroutine fftw_grid_fence_3D(gridarray_t13,fencearray)
        !
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:,:)  ::  gridarray_t13
        real(8), intent(out), dimension(:,:,:) ::  fencearray
        real(8), allocatable, dimension(:,:,:) :: global_t, global
        !
        ! Under implementation
        !
    end subroutine fftw_grid_fence_3D
    !
    subroutine fftw_fence_grid_3D(fencearray,gridarray_t13)
        !
        use parallel
        !
        ! gridarray_t is transposed!
        real(8), intent(in), dimension(:,:,:)  :: fencearray
        real(8), intent(out), dimension(:,:,:) :: gridarray_t13
        real(8), allocatable, dimension(:,:,:) :: global_t, global
        !
        ! Under implementation
        !
    end subroutine fftw_fence_grid_3D
    !
    subroutine prepareplan_fftw_plan
        !
        use commvar, only: u1spe
        use parallel
        !
        include 'fftw3-mpi.f03'
        !
        !
        if(ndims == 2) then
            forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
            backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        elseif(ndims == 3) then
            forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
            backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        else
            stop 'prepareplan_fftw: Not implemented for 1D!'
        endif
        !
    end subroutine prepareplan_fftw_plan
    !
    subroutine allocate_fftw_complex(array, c_pointer)
        use, intrinsic :: iso_c_binding
        include 'fftw3-mpi.f03'
        !
        type(C_PTR), intent(out) :: c_pointer
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: array(:,:)
        !
        c_pointer = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_pointer, array, [imfftw, jmfftw])
        !
    end subroutine allocate_fftw_complex
  !
    subroutine fft2d(array)
        !
        use commvar, only: im,jm,ia,ja
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:)
        integer :: i,j
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(forward_plan,array,array)
        do j=1,jm
        do i=1,im
            array(i,j)=array(i,j)/(1.d0*ia*ja)
        end do
        end do
    end subroutine fft2d
    !
     subroutine ifft2d(array)
        !
        use commvar, only: im,jm,ia,ja
        !
        complex(C_DOUBLE_COMPLEX), pointer, intent(inout) :: array(:,:)
        integer :: i,j
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(backward_plan,array,array)
        !
    end subroutine ifft2d
    !
    subroutine allocation
        !
        use commvar
        include 'fftw3-mpi.f03'
        ! This subroutine allocates the common variables.
        !
        allocate(u1(1:im,1:jm,0:km), u2(1:im,1:jm,0:km))
        allocate(k1(1:im,1:jm),k2(1:im,1:jm))
        !
        call allocate_fftw_complex(u1spe, c_u1spe)
        call allocate_fftw_complex(u2spe, c_u2spe)
        call allocate_fftw_complex(u1x1,   c_u1x1)
        call allocate_fftw_complex(u1x2,   c_u1x2)
        call allocate_fftw_complex(u2x1,   c_u2x1)
        call allocate_fftw_complex(u2x2,   c_u2x2)
        call allocate_fftw_complex(u1xixi, c_u1xixi)
        call allocate_fftw_complex(u2xixi, c_u2xixi)
        call allocate_fftw_complex(u1tA,     c_u1tA)
        call allocate_fftw_complex(u2tA,     c_u2tA)
        call allocate_fftw_complex(u1tB,     c_u1tB)
        call allocate_fftw_complex(u2tB,     c_u2tB)
        call allocate_fftw_complex(u1tC,     c_u1tC)
        call allocate_fftw_complex(u2tC,     c_u2tC)
        !
        allocate(Es(0:allkmax),Ed(0:allkmax),kn(0:allkmax),Ecount(0:allkmax))
        !
    end subroutine allocation
    !
    subroutine deallocation
        !
        use commvar
        include 'fftw3-mpi.f03'
        !
        deallocate(u1, u2, k1,k2,Es, Ed, kn,Ecount)
        !
        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        call fftw_mpi_cleanup()
        call fftw_free(c_u1spe)
        call fftw_free(c_u2spe)
        call fftw_free(c_u1x1)
        call fftw_free(c_u1x2)
        call fftw_free(c_u2x1)
        call fftw_free(c_u2x2)
        call fftw_free(c_u1xixi)
        call fftw_free(c_u2xixi)
        call fftw_free(c_u1tA)
        call fftw_free(c_u2tA)
        call fftw_free(c_u1tB)
        call fftw_free(c_u2tB)
        call fftw_free(c_u1tC)
        call fftw_free(c_u2tC)
    end subroutine deallocation
end module fftwlink