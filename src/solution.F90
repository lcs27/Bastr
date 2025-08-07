module solution
    !
    use commvar
    use parallel
    use fftwlink
    use tool
    use readwrite
    !
    implicit none
    !
    contains
    !
    subroutine RK3
        !
        implicit none
        !
        integer :: i,j
        real(8) :: energy, factor
        ! Warning : The entrance of RK3 must be a Fourier-Transformed u1spe and u2spe
        ! Step 1
        ! 
        !
        call compute_ut(u1tA, u2tA)
        !
        ! Step 2
        do j=1,jm
        do i=1,im
            u1spe(i,j)=CMPLX(u1(i,j,0) + deltat * dreal(u1tA(i,j)) / 2.d0 ,0.d0,C_INTPTR_T);
            u2spe(i,j)=CMPLX(u2(i,j,0) + deltat * dreal(u2tA(i,j)) / 2.d0 ,0.d0,C_INTPTR_T);
        end do
        end do
        call fft2d(u1spe)
        call fft2d(u2spe)
        !
        call compute_ut(u1tB, u2tB)
        !
        ! Step 3
        do j=1,jm
        do i=1,im
            u1spe(i,j) = CMPLX(u1(i,j,0) - deltat * dreal(u1tA(i,j)) + 2.d0*deltat * dreal(u1tB(i,j)), 0.d0, C_INTPTR_T)
            u2spe(i,j) = CMPLX(u2(i,j,0) - deltat * dreal(u2tA(i,j)) + 2.d0*deltat * dreal(u2tB(i,j)), 0.d0, C_INTPTR_T)
        end do
        end do
        call fft2d(u1spe)
        call fft2d(u2spe)
        !
        call compute_ut(u1tC, u2tC)
        !
        ! Final update
        do j=1,jm
        do i=1,im
            u1(i,j,0) = u1(i,j,0) + deltat * (dreal(u1tA(i,j)) + 4.d0 * dreal(u1tB(i,j)) + dreal(u1tC(i,j))) / 6.d0
            u2(i,j,0) = u2(i,j,0) + deltat * (dreal(u2tA(i,j)) + 4.d0 * dreal(u2tB(i,j)) + dreal(u2tC(i,j))) / 6.d0
        end do
        end do
        !
        if(lforce)then
            energy = 0.d0
            do j=1,jm
            do i=1,im
                energy = energy + (u1(i,j,0)**2 + u2(i,j,0)**2) / 2.d0
            end do
            end do
            energy = psum(energy)/(ia*ja)
            !
            factor = sqrt(target_energy/energy)
            !
            energy = 0.d0
            do j=1,jm
            do i=1,im
                u1(i,j,0) = factor * u1(i,j,0) 
                u2(i,j,0) = factor * u2(i,j,0) 
            end do
            end do
        endif
        !
    end subroutine RK3
    !
    subroutine compute_ut(u1t, u2t)
        !
        use commvar, only: u1x1, u1x2, u2x1, u2x2, u1xixi, u2xixi
        use fftwlink, only: fft2d, ifft2d
        use tool, only : dealiasing
        !
        implicit none
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: u1t(:,:), u2t(:,:)
        integer :: i,j
        !
        !
        u1x1 = imag * u1spe * k1 
        u1x2 = imag * u1spe * k2
        u2x1 = imag * u2spe * k1
        u2x2 = imag * u2spe * k2
        u1xixi = - u1spe * (k1*k1 + k2*k2)
        u2xixi = - u2spe * (k1*k1 + k2*k2)
        !
        call ifft2d(u1x1)
        call ifft2d(u1x2)
        call ifft2d(u2x1)
        call ifft2d(u2x2)
        call ifft2d(u1xixi)
        call ifft2d(u2xixi)
        !
        do j=1,jm
        do i=1,im
            !
            u1t(i,j) = - u1(i,j,0) * dreal(u1x1(i,j)) - u2(i,j,0) * dreal(u1x2(i,j)) + nu * dreal(u1xixi(i,j))
            u2t(i,j) = - u1(i,j,0) * dreal(u2x1(i,j)) + u2(i,j,0) * dreal(u2x2(i,j)) + nu * dreal(u2xixi(i,j))
            !
        end do
        end do
        !
        !!!! Do 2d FFT
        !
        call fft2d(u1t)
        call fft2d(u2t)
        !
        call dealiasing(u1t)
        call dealiasing(u2t)
        !
        call ifft2d(u1t)
        call ifft2d(u2t)
        !
    end subroutine compute_ut

    subroutine spectra_compute(hand_f, hand_g)
    !
        use commvar, only: Ed, Es, kn, Ecount, Esspe, Edspe, k1, k2
        use parallel, only: psum
        use utility, only: listwrite
        !
        implicit none
        integer, intent(in) :: hand_f, hand_g
        !
        integer :: i, j, kOrdinal
        real(8) :: kk, dk
        complex(8) :: u1s,u2s,u1d,u2d
        complex(8) ::  usspe,udspe
        !
        dk = 1.d0
        Ed = 0.d0
        Es = 0.d0
        kn = 0.d0
        Ecount = 0
        Edspe = 0.d0
        Esspe = 0.d0
        !
        do j=1,jm
        do i=1,im
            kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
            if(kk>dk/2)then
                usspe = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
                udspe = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
                u1d =  udspe*k1(i,j)/kk
                u2d =  udspe*k2(i,j)/kk
                u1s =  usspe*k2(i,j)/kk 
                u2s = -usspe*k1(i,j)/kk
            else
                usspe = 0
                udspe = 0
                u1d = 0
                u2d = 0
                u1s = 0
                u2s = 0
            endif
            !
            kOrdinal = kint(kk,dk,2,1.d0)
            !
            if(kOrdinal <= allkmax)then
                Ecount(kOrdinal) = Ecount(kOrdinal) + 1
                Es(kOrdinal) = Es(kOrdinal) + usspe*conjg(usspe)*kk/2
                Ed(kOrdinal) = Ed(kOrdinal) + udspe*conjg(udspe)*kk/2
                kn(kOrdinal) = kn(kOrdinal) + kk
            endif
            !
            Edspe = Edspe + (udspe*dconjg(udspe))/2
            Esspe = Esspe + (usspe*dconjg(usspe))/2
            !
        end do
        end do
        !
        !
        do i=1,allkmax
            Ecount(i) = psum(Ecount(i))
            Es(i) = psum(Es(i))/Ecount(i)
            Ed(i) = psum(Ed(i))/Ecount(i)
            kn(i) =  psum(kn(i))/Ecount(i)
        enddo
        Edspe = psum(Edspe)
        Esspe = psum(Esspe)
        !
        if(lio) then
            call listwrite(hand_f,Esspe,Edspe)
            do i=1,allkmax
            call listwrite(hand_g,kn(i),Es(i),Ed(i))
            enddo
        endif
        !
        !
    end subroutine spectra_compute
    !
    !
    subroutine velgrad_calculate(hand_a)
        !
        use commvar, only: u1spe, u2spe, k1, k2, u1x1, u1x2, u2x1, u2x2
        use fftwlink, only: ifft2d
        use utility, only: listwrite
        implicit none
        integer, intent(in) :: hand_a
        integer :: i,j
        real(8) :: div, umumtheta2, umumijji, u2theta
        !
        u1x1 = imag * u1spe * k1 
        u1x2 = imag * u1spe * k2
        u2x1 = imag * u2spe * k1
        u2x2 = imag * u2spe * k2
        !
        call ifft2d(u1x1)
        call ifft2d(u1x2)
        call ifft2d(u2x1)
        call ifft2d(u2x2)
        !
        umumtheta2 = 0.d0
        umumijji = 0.d0
        u2theta = 0.d0
        do j=1,jm
        do i=1,im
            div = dreal(u1x1(i,j)) + dreal(u2x2(i,j))
            umumtheta2 = umumtheta2 + div**2 * (u1(i,j,0)**2 + u2(i,j,0)**2)
            umumijji = umumijji + (dreal(u1x1(i,j))**2 + dreal(u2x2(i,j))**2 + &
                                  2.d0*dreal(u1x2(i,j))*dreal(u2x1(i,j))) * &
                                  (u1(i,j,0)**2 + u2(i,j,0)**2)
            u2theta = u2theta + div * (u1(i,j,0)**2 + u2(i,j,0)**2)
        end do
        end do
        !
        umumtheta2 = psum(umumtheta2)/(ia*ja)
        umumijji = psum(umumijji)/(ia*ja)
        u2theta = psum(u2theta)/(ia*ja)
        !
        if(lio) then
            call listwrite(hand_a,umumtheta2,umumijji,u2theta)
        endif
        !
    end subroutine velgrad_calculate
    !
    subroutine quantity_prepare
        use commvar, only: nstep, nxtwsequ, feqwsequ, filenumb, time, nu
        use tool, only: GenerateWave
        use fftwlink, only: j0f
        use parallel, only: mpirank
        nstep = 0
        nxtwsequ = min(feqwsequ, maxstep)
        filenumb = 1
        call GenerateWave(j0f,k1,k2)
        time = 0.d0
        nu = miucal(ref_tem)/Reynolds
        if(mpirank == 0) print *, 'nu=', nu
    end subroutine quantity_prepare
    !
    real(8) function miucal(temper)
        !
        use commvar, only: ref_tem
        real(8),intent(in) :: temper
        !
        !
        real(8) :: tempconst, tempconst1
        !
        tempconst=110.4d0/ref_tem
        tempconst1=1.d0+tempconst
        miucal=temper*sqrt(temper)*tempconst1/(temper+tempconst)
        !
        return 
        !
    end function miucal
end module