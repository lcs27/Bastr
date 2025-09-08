module solution
    !
    !
    use commvar
    use parallel
    use fftwlink
    use tool
    use readwrite
    !
    implicit none
    !
    interface compute_ut
        module procedure compute_ut2D
        module procedure compute_ut3D
    end interface
    !
    contains
    !
    !
    subroutine mainloop(hand_f, hand_g, hand_a, hand_fo)
        use parallel, only: mpirank
        !
        integer, intent(in) :: hand_f, hand_g, hand_a, hand_fo
        !
        select case(ndims)
        case(2)
            call mainloop2D(hand_f, hand_g, hand_a, hand_fo)
        case(3)
            call mainloop3D(hand_f, hand_g, hand_a, hand_fo) 
        case default
            error stop 'ndims should be 2 or 3!'
        end select
        !
    end subroutine mainloop
    !
    subroutine mainloop2D(hand_f, hand_g, hand_a, hand_fo)
        !
        integer :: i,j
        integer, intent(in) :: hand_f, hand_g, hand_a, hand_fo
        !
        do while(nstep<=maxstep)
            !
            do j=1,jm
            do i=1,im
            u1spe(i,j,1)=CMPLX(u1(i,j,0),0.d0,C_INTPTR_T);
            u2spe(i,j,1)=CMPLX(u2(i,j,0),0.d0,C_INTPTR_T);
            end do
            end do
            !
            call fft2d(u1spe)
            call fft2d(u2spe)
            !
            call spectra_compute2D(hand_f, hand_g)
            !
            call velgrad_calculate2D(hand_a)
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
            call RK32D
            !
            !
            call forcing2D(hand_fo)
            !
        end do
    end subroutine mainloop2D
    !
    subroutine mainloop3D(hand_f, hand_g, hand_a, hand_fo)
        !
        integer :: i,j,k
        integer, intent(in) :: hand_f, hand_g, hand_a, hand_fo
        !
        do while(nstep<=maxstep)
            !
            do k=1,km
            do j=1,jm
            do i=1,im
            u1spe(i,j,k)=CMPLX(u1(i,j,k),0.d0,C_INTPTR_T);
            u2spe(i,j,k)=CMPLX(u2(i,j,k),0.d0,C_INTPTR_T);
            u3spe(i,j,k)=CMPLX(u3(i,j,k),0.d0,C_INTPTR_T);
            end do
            end do
            end do
            !
            call fft3d(u1spe)
            call fft3d(u2spe)
            call fft3d(u3spe)
            !
            call spectra_compute3D(hand_f, hand_g)
            !
            call velgrad_calculate3D(hand_a)
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
            call RK33D
            !
            !
            call forcing3D(hand_fo)
            !
        end do
    end subroutine mainloop3D
    !
    subroutine RK32D
        !
        implicit none
        integer :: i,j
        ! Warning : The entrance of RK3 must be a Fourier-Transformed u1spe and u2spe
        ! Step 1
        ! 
        !
        call compute_ut(u1tA, u2tA)
        !
        ! Step 2
        do j=1,jm
        do i=1,im
            u1spe(i,j,1)=CMPLX(u1(i,j,0) + deltat * dreal(u1tA(i,j,1)) / 2.d0 ,0.d0,C_INTPTR_T);
            u2spe(i,j,1)=CMPLX(u2(i,j,0) + deltat * dreal(u2tA(i,j,1)) / 2.d0 ,0.d0,C_INTPTR_T);
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
            u1spe(i,j,1) = CMPLX(u1(i,j,0) - deltat * dreal(u1tA(i,j,1)) + 2.d0*deltat * dreal(u1tB(i,j,1)), 0.d0, C_INTPTR_T)
            u2spe(i,j,1) = CMPLX(u2(i,j,0) - deltat * dreal(u2tA(i,j,1)) + 2.d0*deltat * dreal(u2tB(i,j,1)), 0.d0, C_INTPTR_T)
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
            u1(i,j,0) = u1(i,j,0) + deltat * (dreal(u1tA(i,j,1)) + 4.d0 * dreal(u1tB(i,j,1)) + dreal(u1tC(i,j,1))) / 6.d0
            u2(i,j,0) = u2(i,j,0) + deltat * (dreal(u2tA(i,j,1)) + 4.d0 * dreal(u2tB(i,j,1)) + dreal(u2tC(i,j,1))) / 6.d0
        end do
        end do
        !
        !
    end subroutine RK32D
    !
    subroutine RK33D
        !
        implicit none
        integer :: i,j,k
        ! Warning : The entrance of RK3 must be a Fourier-Transformed u1spe and u2spe
        ! Step 1
        ! 
        !
        call compute_ut(u1tA, u2tA, u3tA)
        !
        ! Step 2
        do k=1,km
        do j=1,jm
        do i=1,im
            u1spe(i,j,k)=CMPLX(u1(i,j,k) + deltat * dreal(u1tA(i,j,k)) / 2.d0 ,0.d0,C_INTPTR_T);
            u2spe(i,j,k)=CMPLX(u2(i,j,k) + deltat * dreal(u2tA(i,j,k)) / 2.d0 ,0.d0,C_INTPTR_T);
            u3spe(i,j,k)=CMPLX(u3(i,j,k) + deltat * dreal(u3tA(i,j,k)) / 2.d0 ,0.d0,C_INTPTR_T);
        end do
        end do
        end do
        !
        call fft3d(u1spe)
        call fft3d(u2spe)
        call fft3d(u3spe)
        !
        call compute_ut(u1tB, u2tB, u3tB)
        !
        ! Step 3
        do k=1,km
        do j=1,jm
        do i=1,im
            u1spe(i,j,k) = CMPLX(u1(i,j,k) - deltat * dreal(u1tA(i,j,k)) + 2.d0*deltat * dreal(u1tB(i,j,k)), 0.d0, C_INTPTR_T)
            u2spe(i,j,k) = CMPLX(u2(i,j,k) - deltat * dreal(u2tA(i,j,k)) + 2.d0*deltat * dreal(u2tB(i,j,k)), 0.d0, C_INTPTR_T)
            u3spe(i,j,k) = CMPLX(u3(i,j,k) - deltat * dreal(u3tA(i,j,k)) + 2.d0*deltat * dreal(u3tB(i,j,k)), 0.d0, C_INTPTR_T)
        end do
        end do
        end do
        !
        call fft3d(u1spe)
        call fft3d(u2spe)
        call fft3d(u3spe)
        !
        call compute_ut(u1tC, u2tC, u3tC)
        !
        ! Final update
        do k=1,km
        do j=1,jm
        do i=1,im
            u1(i,j,k) = u1(i,j,k) + deltat * (dreal(u1tA(i,j,k)) + 4.d0 * dreal(u1tB(i,j,k)) + dreal(u1tC(i,j,k))) / 6.d0
            u2(i,j,k) = u2(i,j,k) + deltat * (dreal(u2tA(i,j,k)) + 4.d0 * dreal(u2tB(i,j,k)) + dreal(u2tC(i,j,k))) / 6.d0

        end do
        end do
        end do
        !
        !
    end subroutine RK33D
    !
    subroutine forcing2D(hand_fo)
        use utility, only: listwrite
        implicit none
        !
        integer, intent(in) :: hand_fo
        integer :: i,j
        real(8) :: energy, factor
        integer :: forcek
        real(8) :: kk,dk,E
        !
        if(forcemethod == 1)then
            ! Linear forcing with f_i = (factor-1) * u_i
            ! s.t. the kinetic energy is maintained at target_energy
            ! All in physical space
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
            do j=1,jm
            do i=1,im
                u1(i,j,0) = factor * u1(i,j,0) 
                u2(i,j,0) = factor * u2(i,j,0) 
            end do
            end do
            !
            if (lio) then
                call listwrite(hand_fo,factor)
            endif
            !
        elseif(forcemethod == 2) then
            ! Linear forcing in a band of wave numbers in spectral space
            forcek = 5
            dk = 0.5d0
            !
            do j=1,jm
            do i=1,im
                u1spe(i,j,1)=CMPLX(u1(i,j,0),0.d0,C_INTPTR_T);
                u2spe(i,j,1)=CMPLX(u2(i,j,0),0.d0,C_INTPTR_T);
            end do
            end do
            !
            call fft2d(u1spe)
            call fft2d(u2spe)
            !
            do j=1,jm
            do i=1,im
                kk=dsqrt(k1(i,j,0)**2+k2(i,j,0)**2)
                !
                if((kk - dk)<forcek .and. (kk + dk)>forcek) then
                    force1(i,j,1) = u1spe(i,j,1)
                    force2(i,j,1) = u2spe(i,j,1)
                else
                    force1(i,j,1) = 0.d0
                    force2(i,j,1) = 0.d0
                endif
            end do
            end do
            !
            call ifft2d(force1)
            call ifft2d(force2)
            !
            E = 0.d0
            do j=1,jm
            do i=1,im
                E = E + dreal(force1(i,j,1))**2 + dreal(force2(i,j,1))**2
            end do
            end do
            E = psum(E)/(ia*ja)
            factor = max(sqrt(target_energy/E) - 1.d0, 0.d0)
            !
            do j=1,jm
            do i=1,im
                u1(i,j,0) = u1(i,j,0) + factor * dreal(force1(i,j,1))
                u2(i,j,0) = u2(i,j,0) + factor * dreal(force2(i,j,1))
            end do
            end do
            !
            if (lio) then
                call listwrite(hand_fo,factor,E)
            endif
            !
        endif
    end subroutine forcing2D
    !
    subroutine forcing3D(hand_fo)
        use utility, only: listwrite
        implicit none
        !
        integer, intent(in) :: hand_fo
        integer :: i,j,k
        real(8) :: energy, factor
        integer :: forcek
        real(8) :: kk,dk,E
        !
        if(forcemethod == 1)then
            ! Linear forcing with f_i = (factor-1) * u_i
            ! s.t. the kinetic energy is maintained at target_energy
            ! All in physical space
            energy = 0.d0
            do k=1,km
            do j=1,jm
            do i=1,im
                energy = energy + (u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2) / 2.d0
            end do
            end do
            end do
            energy = psum(energy)/(ia*ja*ka)
            !
            factor = sqrt(target_energy/energy)
            !
            do k=1,km
            do j=1,jm
            do i=1,im
                u1(i,j,k) = factor * u1(i,j,k) 
                u2(i,j,k) = factor * u2(i,j,k)
                u3(i,j,k) = factor * u3(i,j,k) 
            end do
            end do
            end do
            !
            if (lio) then
                call listwrite(hand_fo,factor)
            endif
            !
        elseif(forcemethod == 2) then
            ! Linear forcing in a band of wave numbers in spectral space
            forcek = 5
            dk = 0.5d0
            !
            do k=1,km
            do j=1,jm
            do i=1,im
                u1spe(i,j,k)=CMPLX(u1(i,j,k),0.d0,C_INTPTR_T);
                u2spe(i,j,k)=CMPLX(u2(i,j,k),0.d0,C_INTPTR_T);
                u3spe(i,j,k)=CMPLX(u3(i,j,k),0.d0,C_INTPTR_T);
            end do
            end do
            end do
            !
            call fft3d(u1spe)
            call fft3d(u2spe)
            call fft3d(u3spe)
            !
            do k=1,km
            do j=1,jm
            do i=1,im
                kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)
                !
                if((kk - dk)<forcek .and. (kk + dk)>forcek) then
                    force1(i,j,k) = u1spe(i,j,k)
                    force2(i,j,k) = u2spe(i,j,k)
                    force3(i,j,k) = u3spe(i,j,k)
                else
                    force1(i,j,k) = 0.d0
                    force2(i,j,k) = 0.d0
                    force3(i,j,k) = 0.d0
                endif
            end do
            end do
            end do
            !
            call ifft3d(force1)
            call ifft3d(force2)
            call ifft3d(force3)
            !
            E = 0.d0
            do k=1,km
            do j=1,jm
            do i=1,im
                E = E + dreal(force1(i,j,k))**2 + dreal(force2(i,j,k))**2 + dreal(force3(i,j,k))**2
            end do
            end do
            end do
            E = psum(E)/(ia*ja*ka)
            factor = max(sqrt(target_energy/E) - 1.d0, 0.d0)
            !
            do k=1,km
            do j=1,jm
            do i=1,im
                u1(i,j,k) = u1(i,j,k) + factor * dreal(force1(i,j,k))
                u2(i,j,k) = u2(i,j,k) + factor * dreal(force2(i,j,k))
                u3(i,j,k) = u3(i,j,k) + factor * dreal(force3(i,j,k))
            end do
            end do
            end do
            !
            if (lio) then
                call listwrite(hand_fo,factor,E)
            endif
            !
        endif
    end subroutine forcing3D
    !
    subroutine compute_ut2D(u1t, u2t)
        !
        use fftwlink, only: fft2d, ifft2d
        use tool, only : dealiasing
        !
        implicit none
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: u1t(:,:,:), u2t(:,:,:)
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
        call ifft2d(u1spe)
        call ifft2d(u2spe)
        !
        do j=1,jm
        do i=1,im
            !
            u1t(i,j,1) = - u1spe(i,j,1) * dreal(u1x1(i,j,1)) - u2spe(i,j,1) * dreal(u1x2(i,j,1)) + nu * dreal(u1xixi(i,j,1))
            u2t(i,j,1) = - u1spe(i,j,1) * dreal(u2x1(i,j,1)) - u2spe(i,j,1) * dreal(u2x2(i,j,1)) + nu * dreal(u2xixi(i,j,1))
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
        if(lprojectd)then
            call projection(u1t, u2t)
        endif
        !
        call ifft2d(u1t)
        call ifft2d(u2t)
        !
    end subroutine compute_ut2D
    !
    subroutine compute_ut3D(u1t, u2t,u3t)
        !
        use fftwlink, only: fft3d, ifft3d
        use tool, only : dealiasing
        !
        implicit none
        complex(C_DOUBLE_COMPLEX), pointer, intent(out) :: u1t(:,:,:), u2t(:,:,:), u3t(:,:,:)
        integer :: i,j,k
        !
        !
        u1x1 = imag * u1spe * k1 
        u1x2 = imag * u1spe * k2
        u1x3 = imag * u1spe * k3
        u2x1 = imag * u2spe * k1
        u2x2 = imag * u2spe * k2
        u2x3 = imag * u2spe * k3
        u3x1 = imag * u3spe * k1
        u3x2 = imag * u3spe * k2
        u3x3 = imag * u3spe * k3
        u1xixi = - u1spe * (k1*k1 + k2*k2 + k3*k3)
        u2xixi = - u2spe * (k1*k1 + k2*k2 + k3*k3)
        u3xixi = - u3spe * (k1*k1 + k2*k2 + k3*k3)
        !
        call ifft3d(u1x1)
        call ifft3d(u1x2)
        call ifft3d(u1x3)
        call ifft3d(u2x1)
        call ifft3d(u2x2)
        call ifft3d(u2x3)
        call ifft3d(u3x1)
        call ifft3d(u3x2)
        call ifft3d(u3x3)
        call ifft3d(u1xixi)
        call ifft3d(u2xixi)
        call ifft3d(u3xixi)
        call ifft3d(u1spe)
        call ifft3d(u2spe)
        call ifft3d(u3spe)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
            !
            u1t(i,j,k) = - u1spe(i,j,k) * dreal(u1x1(i,j,k)) - u2spe(i,j,k) * dreal(u1x2(i,j,k)) &
                         - u3spe(i,j,k) * dreal(u1x3(i,j,k)) + nu * dreal(u1xixi(i,j,k))
            u2t(i,j,k) = - u1spe(i,j,k) * dreal(u2x1(i,j,k)) - u2spe(i,j,k) * dreal(u2x2(i,j,k)) &
                         - u3spe(i,j,k) * dreal(u2x3(i,j,k)) + nu * dreal(u2xixi(i,j,k))
            u3t(i,j,k) = - u1spe(i,j,k) * dreal(u3x1(i,j,k)) - u2spe(i,j,k) * dreal(u3x2(i,j,k)) &
                         - u3spe(i,j,k) * dreal(u3x3(i,j,k)) + nu * dreal(u3xixi(i,j,k))
            !
        end do
        end do
        end do
        !
        !!!! Do 2d FFT
        !
        call fft3d(u1t)
        call fft3d(u2t)
        call fft3d(u3t)
        !
        call dealiasing(u1t)
        call dealiasing(u2t)
        call dealiasing(u3t)
        !
        if(lprojectd)then
            call projection(u1t, u2t, u3t)
        endif
        !
        call ifft3d(u1t)
        call ifft3d(u2t)
        call ifft3d(u3t)
        !
    end subroutine compute_ut3D
    !
    subroutine spectra_compute2D(hand_f, hand_g)
    !
        use parallel, only: psum
        use utility, only: listwrite
        !
        implicit none
        integer, intent(in) :: hand_f, hand_g
        !
        integer :: i, j, kOrdinal
        real(8) :: kk, dk
        real(8) :: k2Edspe, kdEdspe, k2dEdspe
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
        k2Edspe = 0.d0
        kdEdspe = 0.d0
        k2dEdspe = 0.d0
        !
        do j=1,jm
        do i=1,im
            kk=dsqrt(k1(i,j,0)**2+k2(i,j,0)**2)
            if(kk>dk/2)then
                usspe = u1spe(i,j,1)*k2(i,j,0)/kk - u2spe(i,j,1)*k1(i,j,0)/kk
                udspe = u1spe(i,j,1)*k1(i,j,0)/kk + u2spe(i,j,1)*k2(i,j,0)/kk
                u1d =  udspe*k1(i,j,0)/kk
                u2d =  udspe*k2(i,j,0)/kk
                u1s =  usspe*k2(i,j,0)/kk 
                u2s = -usspe*k1(i,j,0)/kk
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
            k2Edspe = k2Edspe + (udspe*dconjg(udspe))/2 * (kk ** 2)
            kdEdspe = kdEdspe + (udspe*dconjg(udspe))/2 * (kk ** 0.5)
            k2dEdspe = k2dEdspe + (udspe*dconjg(udspe))/2 * (kk ** 2.5)
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
        k2Edspe = psum(k2Edspe)
        kdEdspe = psum(kdEdspe)
        k2dEdspe = psum(k2dEdspe)
        !
        if(lio) then
            call listwrite(hand_f,Esspe,Edspe,k2Edspe,kdEdspe,k2dEdspe)
            if(lwspectra .and. nstep==nxtwspe) then
                !
                do i=1,allkmax
                    call listwrite(hand_g,kn(i),Es(i),Ed(i))
                enddo
                !
                nxtwspe = min(nxtwspe + feqwspe, maxstep)
                !
            endif
        endif
        !
        !
    end subroutine spectra_compute2D
    !
    subroutine spectra_compute3D(hand_f, hand_g)
    !
        use parallel, only: psum
        use utility, only: listwrite
        !
        implicit none
        integer, intent(in) :: hand_f, hand_g
        !
        integer :: i, j, k, kOrdinal
        real(8) :: kk, dk
        real(8) :: k2Edspe, kdEdspe, k2dEdspe
        complex(8) :: u1s,u2s,u3s,u1d,u2d,u3d
        complex(8) ::  udspe
        !
        dk = 1.d0
        Ed = 0.d0
        Es = 0.d0
        kn = 0.d0
        Ecount = 0
        Edspe = 0.d0
        Esspe = 0.d0
        k2Edspe = 0.d0
        kdEdspe = 0.d0
        k2dEdspe = 0.d0
        !
        do k=1,km
        do j=1,jm
        do i=1,im
            kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)
            if(kk>dk/2)then
                udspe = u1spe(i,j,k)*k1(i,j,k)/kk + u2spe(i,j,k)*k2(i,j,k)/kk + u3spe(i,j,k)*k3(i,j,k)/kk
                u1d =  udspe*k1(i,j,k)/kk
                u2d =  udspe*k2(i,j,k)/kk
                u3d =  udspe*k3(i,j,k)/kk
                u1s =  u1spe(i,j,k) - u1d
                u2s =  u2spe(i,j,k) - u2d
                u3s =  u3spe(i,j,k) - u3d
            else
                udspe = 0
                u1d = 0
                u2d = 0
                u3d = 0
                u1s = 0
                u2s = 0
                u3s = 0
            endif
            !
            kOrdinal = kint(kk,dk,2,1.d0)
            !
            if(kOrdinal <= allkmax)then
                Ecount(kOrdinal) = Ecount(kOrdinal) + 1
                Es(kOrdinal) = Es(kOrdinal) + u1s*conjg(u1s)*kk/2 + u2s*conjg(u2s)*kk/2 + u3s*conjg(u3s)*kk/2
                Ed(kOrdinal) = Ed(kOrdinal) + udspe*conjg(udspe)*kk/2
                kn(kOrdinal) = kn(kOrdinal) + kk
            endif
            !
            Edspe = Edspe + (udspe*dconjg(udspe))/2
            Esspe = Esspe + u1s*conjg(u1s)/2 + u2s*conjg(u2s)/2 + u3s*conjg(u3s)/2
            k2Edspe = k2Edspe + (udspe*dconjg(udspe))/2 * (kk ** 2)
            kdEdspe = kdEdspe + (udspe*dconjg(udspe))/2 * (kk ** 0.5)
            k2dEdspe = k2dEdspe + (udspe*dconjg(udspe))/2 * (kk ** 2.5)
            !
        end do
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
        k2Edspe = psum(k2Edspe)
        kdEdspe = psum(kdEdspe)
        k2dEdspe = psum(k2dEdspe)
        !
        if(lio) then
            call listwrite(hand_f,Esspe,Edspe,k2Edspe,kdEdspe,k2dEdspe)
            if(lwspectra .and. nstep==nxtwspe) then
                !
                do i=1,allkmax
                    call listwrite(hand_g,kn(i),Es(i),Ed(i))
                enddo
                !
                nxtwspe = min(nxtwspe + feqwspe, maxstep)
                !
            endif
        endif
        !
        !
    end subroutine spectra_compute3D
    !
    subroutine velgrad_calculate2D(hand_a)
        !
        use commvar, only: u1spe, u2spe, k1, k2, u1x1, u1x2, u2x1, u2x2, &
                            u1xixi, u2xixi, thetaxixi, eta_min
        use fftwlink, only: ifft2d
        use utility, only: listwrite
        implicit none
        integer, intent(in) :: hand_a
        integer :: i,j
        real(8) :: div, umumtheta2, umumijji, u2theta, dissp,epsilon
        real(8) :: urms, ufluc, Taylength, ReTay, dudx2, Kollength
        real(8), allocatable, dimension(:,:,:) :: kk
        allocate(kk(1:im,1:jm,0:km))
        !
        u1x1 = imag * u1spe * k1 
        u1x2 = imag * u1spe * k2
        u2x1 = imag * u2spe * k1
        u2x2 = imag * u2spe * k2
        kk = k1*k1 + k2*k2
        u1xixi = - u1spe * kk
        u2xixi = - u2spe * kk
        thetaxixi = - imag * (u1spe * k1 * kk + u2spe * k2 * kk)
        !
        call ifft2d(u1x1)
        call ifft2d(u1x2)
        call ifft2d(u2x1)
        call ifft2d(u2x2)
        call ifft2d(u1xixi)
        call ifft2d(u2xixi)
        call ifft2d(thetaxixi)
        !
        umumtheta2 = 0.d0
        umumijji = 0.d0
        u2theta = 0.d0
        dissp = 0.d0
        eta_min = 2*3.14
        dudx2 = 0.d0
        urms = 0.d0
        epsilon = 0.d0
        do j=1,jm
        do i=1,im
            div = dreal(u1x1(i,j,1)) + dreal(u2x2(i,j,1))
            umumtheta2 = umumtheta2 + div**2 * (u1(i,j,0)**2 + u2(i,j,0)**2)
            umumijji = umumijji + (dreal(u1x1(i,j,1))**2 + dreal(u2x2(i,j,1))**2 + &
                                  2.d0*dreal(u1x2(i,j,1))*dreal(u2x1(i,j,1))) * &
                                  (u1(i,j,0)**2 + u2(i,j,0)**2)
            u2theta = u2theta + div * (u1(i,j,0)**2 + u2(i,j,0)**2)
            dissp = dissp + 2.d0 * u1(i,j,0) * div * dreal(u1xixi(i,j,1)) &
                          + 2.d0 * u2(i,j,0) * div * dreal(u2xixi(i,j,1)) &
                          + (u1(i,j,0)**2 + u2(i,j,0)**2) * dreal(thetaxixi(i,j,1))
            epsilon = epsilon + nu * (dreal(u1x2(i,j,1)) - dreal(u2x1(i,j,1)))**2 &
                              + 4.d0/3.d0 * nu * div**2
            eta_min = min(eta_min, dsqrt(nu/abs(div)))
            urms = urms+ u1(i,j,0)**2 + u2(i,j,0)**2
            dudx2 = dudx2 + dreal(u1x1(i,j,1))**2 + dreal(u2x2(i,j,1))**2
        end do
        end do
        !
        umumtheta2 = psum(umumtheta2)/(ia*ja)
        umumijji = psum(umumijji)/(ia*ja)
        u2theta = psum(u2theta)/(ia*ja)
        dissp = psum(dissp)/(ia*ja)
        eta_min = pmin(eta_min)
        urms = sqrt(psum(urms)/(ia*ja))
        dudx2 = psum(dudx2)/(ia*ja)
        epsilon = psum(epsilon)/(ia*ja)
        !
        ufluc     = urms/sqrt(2.d0)
        Taylength = ufluc/sqrt(dudx2/2.d0)
        ReTay     = ufluc * Taylength / nu
        Kollength = sqrt(sqrt(nu**3/epsilon))
        !
        if(lio) then
            call listwrite(hand_a,umumtheta2,umumijji,u2theta,dissp,eta_min, Taylength, ReTay, Kollength)
        endif
        !
    end subroutine velgrad_calculate2D
    !
    subroutine velgrad_calculate3D(hand_a)
        !
        use commvar, only: u1spe, u2spe, u3spe, k1, k2, k3, u1x1, u1x2, u1x3, &
                            u2x1, u2x2, u2x3, u3x1, u3x2, u3x3,&
                            u1xixi, u2xixi, u3xixi, thetaxixi, eta_min
        use fftwlink, only: ifft3d
        use utility, only: listwrite
        implicit none
        integer, intent(in) :: hand_a
        integer :: i,j,k
        real(8) :: div, umumtheta2, umumijji, u2theta, dissp,epsilon
        real(8) :: urms, ufluc, Taylength, ReTay, dudx2, Kollength
        real(8), allocatable, dimension(:,:,:) :: kk
        allocate(kk(1:im,1:jm,0:km))
        !
        u1x1 = imag * u1spe * k1 
        u1x2 = imag * u1spe * k2
        u1x3 = imag * u1spe * k3
        u2x1 = imag * u2spe * k1
        u2x2 = imag * u2spe * k2
        u2x3 = imag * u2spe * k3
        u3x1 = imag * u3spe * k1
        u3x2 = imag * u3spe * k2
        u3x3 = imag * u3spe * k3
        kk = k1*k1 + k2*k2 + k3*k3
        u1xixi = - u1spe * kk
        u2xixi = - u2spe * kk
        u3xixi = - u3spe * kk
        thetaxixi = - imag * (u1spe * k1 * kk + u2spe * k2 * kk + u3spe * k3 * kk)
        !
        call ifft3d(u1x1)
        call ifft3d(u1x2)
        call ifft3d(u1x3)
        call ifft3d(u2x1)
        call ifft3d(u2x2)
        call ifft3d(u2x3)
        call ifft3d(u3x1)
        call ifft3d(u3x2)
        call ifft3d(u3x3)
        call ifft3d(u1xixi)
        call ifft3d(u2xixi)
        call ifft3d(u3xixi)
        call ifft3d(thetaxixi)
        !
        umumtheta2 = 0.d0
        umumijji = 0.d0
        u2theta = 0.d0
        dissp = 0.d0
        eta_min = 2*3.14
        dudx2 = 0.d0
        urms = 0.d0
        epsilon = 0.d0
        do k=1,km
        do j=1,jm
        do i=1,im
            div = dreal(u1x1(i,j,k)) + dreal(u2x2(i,j,k)) + dreal(u3x3(i,j,k))
            umumtheta2 = umumtheta2 + div**2 * (u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2)
            umumijji = umumijji + (dreal(u1x1(i,j,k))**2 + &
                                   dreal(u2x2(i,j,k))**2 + &
                                   dreal(u3x3(i,j,k))**2 + &
                                   2.d0*dreal(u1x2(i,j,k))*dreal(u2x1(i,j,k)) + &
                                   2.d0*dreal(u1x3(i,j,k))*dreal(u3x1(i,j,k)) + &
                                   2.d0*dreal(u2x3(i,j,k))*dreal(u3x2(i,j,k))) *&                                  
                                   (u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2)
            u2theta = u2theta + div * (u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2)
            dissp = dissp + 2.d0 * u1(i,j,k) * div * dreal(u1xixi(i,j,k)) &
                          + 2.d0 * u2(i,j,k) * div * dreal(u2xixi(i,j,k)) &
                          + 2.d0 * u3(i,j,k) * div * dreal(u3xixi(i,j,k)) &
                          + (u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2) * dreal(thetaxixi(i,j,k))
            epsilon = epsilon + nu * (dreal(u1x2(i,j,k)) - dreal(u2x1(i,j,k)))**2 &
                              + nu * (dreal(u2x3(i,j,k)) - dreal(u3x2(i,j,k)))**2 &
                              + nu * (dreal(u1x3(i,j,k)) - dreal(u3x1(i,j,k)))**2 &
                              + 4.d0/3.d0 * nu * div**2
            eta_min = min(eta_min, dsqrt(nu/abs(div)))
            urms = urms+ u1(i,j,k)**2 + u2(i,j,k)**2 + u3(i,j,k)**2
            dudx2 = dudx2 + dreal(u1x1(i,j,k))**2 + dreal(u2x2(i,j,k))**2 + dreal(u3x3(i,j,k))**2
        end do
        end do
        end do
        !
        umumtheta2 = psum(umumtheta2)/(ia*ja*ka)
        umumijji = psum(umumijji)/(ia*ja*ka)
        u2theta = psum(u2theta)/(ia*ja*ka)
        dissp = psum(dissp)/(ia*ja*ka)
        eta_min = pmin(eta_min)
        urms = sqrt(psum(urms)/(ia*ja*ka))
        dudx2 = psum(dudx2)/(ia*ja*ka)
        epsilon = psum(epsilon)/(ia*ja*ka)
        !
        ufluc     = urms/sqrt(3.d0)
        Taylength = ufluc/sqrt(dudx2/3.d0)
        ReTay     = ufluc * Taylength / nu
        Kollength = sqrt(sqrt(nu**3/epsilon))
        !
        if(lio) then
            call listwrite(hand_a,umumtheta2,umumijji,u2theta,dissp,eta_min, Taylength, ReTay, Kollength)
        endif
        !
    end subroutine velgrad_calculate3D
    !
    subroutine quantity_prepare
        use commvar, only: nstep, nxtwsequ, feqwsequ, filenumb, time, nu, nxtwspe, feqwspe, allkmax, ia, ja
        use tool, only: GenerateWave
        use fftwlink, only: j0f,k0f
        use parallel, only: mpirank
        nstep = 0
        nxtwsequ = min(feqwsequ, maxstep)
        nxtwspe = 0
        filenumb = 1
        select case(ndims)
        case(2)
            call GenerateWave(j0f,k1,k2)
        case(3)
            call GenerateWave(k0f,k1,k2,k3)
        case default
            error stop 'ndims should be 2 or 3!'
        end select
        time = 0.d0
        nu = miucal(1.d0)/Reynolds
        if(mpirank == 0) print *, 'nu=', nu
        if(mpirank==0) print *, 'allkmax=', allkmax
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