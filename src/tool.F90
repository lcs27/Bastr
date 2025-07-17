!+---------------------------------------------------------------------+
!+---------------------------------------------------------------------+
module tool
    use commvar, only: im,jm,km,ia,ja,ka
    interface GenerateWave
        module procedure GenerateWave_2D
        module procedure GenerateWave_3D
    end interface
    !
    contains
    !
    subroutine GenerateWave_2D(j0,k1,k2)
        implicit none
        integer, intent(in) :: j0
        real(8), dimension(:,:),intent(out) :: k1,k2
        integer :: i,j
        !
        do j=1,jm
        do i=1,im
        !
        if((im .ne. ia) .and. ((2*im-2) .ne. ia))then
            stop "GenerateWave Error! im /= ia  and (2*im-2) /= ia"
        endif
        !
        if(i <= (ia/2+1)) then
            k1(i,j) = real(i-1,8)
        else if(i<=(ia)) then
            k1(i,j) = real(i-ia-1,8)
        else
            stop "GenerateWave Error, no wave number possible, i must smaller than ia-1 !"
        end if
        !
        if((j+j0) <= (ja/2+1)) then
            k2(i,j) = real(j+j0-1,8)
        else if((j+j0)<=(ja)) then
            k2(i,j) = real(j+j0-ja-1,8)
        else
            stop "GenerateWave Error, no wave number possible, (j+j0) must smaller than ja-1 !"
        end if
        !
        end do
        end do
    end subroutine GenerateWave_2D
    !
    subroutine GenerateWave_3D(k0,k1,k2,k3)
        implicit none
        integer, intent(in) :: k0
        real(8), dimension(:,:,:),intent(out) :: k1,k2,k3
        integer :: i,j,k
        !
        do k=1,km
        do j=1,jm
        do i=1,im
        !
        if((im .ne. ia) .and. ((2*im-2) .ne. ia))then
            stop "GenerateWave Error! im /= ia and (2*im-2) /= ia"
        endif
        !
        if(jm .ne. ja)then
            stop "GenerateWave Error! jm /= ja"
        endif
        !
        if(i <= (ia/2+1)) then
            k1(i,j,k) = real(i-1,8)
        else if(i<=ia) then
            k1(i,j,k) = real(i-ia-1,8)
        else
            stop "GenerateWave Error, no wave number possible, i must smaller than ia-1 !"
        end if
        !
        if(j <= (ja/2+1)) then
            k2(i,j,k) = real(j-1,8)
        else if(j<=ja) then
            k2(i,j,k) = real(j-ja-1,8)
        else
            stop "GenerateWave Error, no wave number possible, j must smaller than ja-1 !"
        end if
        !
        if((k+k0) <= (ka/2+1)) then
            k3(i,j,k) = real(k+k0-1,8)
        else if((k+k0)<=ka) then
            k3(i,j,k) = real(k+k0-ka-1,8)
        else
            stop "GenerateWave Error, no wave number possible, (k+k0) must smaller than ka-1 !"
        end if
        !
        enddo
        enddo
        enddo
    end subroutine GenerateWave_3D
    !
    real(8) function wav(i)
    ! This function gives the wave number of index i
    ! with the maximum im
        implicit none
        integer, intent(in) :: i
        !
        if(i <= (im/2+1)) then
        wav = real(i-1,8)
        else if(i<=im) then
        wav = real(i-im-1,8)
        else
        print *,"Error, no wave number possible, i must smaller than im !"
        end if
    end function wav
    !
    integer function invwav(i)
    ! This function gives the wave number of index i
    ! with the maximum im
        implicit none
        integer, intent(in) :: i
        !
        if(i < 0) then
        invwav = i + im + 1;
        else
        invwav = i + 1;
        end if
    end function invwav
    !
    integer function kint(k,dk,method,lambda)
    ! This function gives the nearby k
    !!
        implicit none
        real(8), intent(in) :: k,dk
        integer, intent(in) :: method
        real(8), intent(in), optional :: lambda
        !
        if(method == 1)then
        if( (.not. present(lambda)) .or. (lambda <= 0))then
            stop "Error! kint method 1 with no lambda or lambda is negative!"
        else
            if(k<(dk/2))then
                kint = 0
            else
                kint = floor(log(k/dk)/log(lambda))+1
            endif
        endif
        else
        kint = nint(k/dk)
        endif
    end function kint
    !
    subroutine dealiasing(field)
        ! Applies 2/3-rule dealiasing to a 2D spectral field
        use commvar, only : k1, k2
        implicit none
        complex(8), dimension(:,:), intent(inout) :: field
        integer :: i, j, cutoff
        !
        cutoff = int(ia / 3)
        do j=1,jm
        do i=1,im
            if (abs(k1(i,j)) > cutoff .or. abs(k2(i,j)) > cutoff) then
                field(i,j) = 0.d0
            end if
        end do
        end do
    end subroutine dealiasing
end module tool