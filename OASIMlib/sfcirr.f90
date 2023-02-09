submodule (oasim) oasim_sfcirr
    implicit none

contains
    module subroutine sfcirr(self, iday, sec_c, slp, wsm, oz, wv, rh, &
         taua, asymp, ssalb, ccov, rlwp, cdre, error)
        implicit none

        class(calc_unit) :: self
        integer, intent(in) :: iday
        real(kind=real_kind), intent(in) :: sec_c
        real(kind=real_kind), dimension(:), intent(in) :: slp, wsm, oz, wv, rh
        real(kind=real_kind), dimension(:,:), intent(in) :: taua, asymp, ssalb
        real(kind=real_kind), dimension(:), intent(in) :: ccov, rlwp, cdre
        logical, intent(out) :: error

        real(kind=real_kind), parameter :: daypersec = 1.0d0 / 86400.0d0

        real(kind=real_kind) :: rday, daycor, sunz, cosunz, pres, ws, ozone, wvapor, relhum
        real(kind=real_kind) :: cov, clwp, re !, sirr
        integer :: i !, j

        self%eda = 0.0d0
        self%esa = 0.0d0

        rday = real(iday, real_kind) + sec_c * daypersec
        daycor = 1.0 + 1.67d-2 * cos(pi2 * (rday - 3.0d0) / 365.0d0)
        daycor = daycor * daycor

        do i = 1, self%p_size
            cosunz = cos(self%solz(i) * rad_1)
            sunz = self%solz(i)

            if (sunz < 90.0d0) then
                pres = slp(i)
                ws = wsm(i)
                ozone = oz(i)
                wvapor = wv(i)
                relhum = rh(i)

                self%ta = taua(i,:)
                self%asym = asymp(i,:)
                self%wa = ssalb(i,:)

                cov = ccov(i)
                clwp = rlwp(i)
                re = cdre(i)

                call self%light(sunz, cosunz, daycor, pres, ws, ozone, wvapor, relhum, &
                                self%lib%init_parameters%am, self%lib%init_parameters%vi, &
                                cov, clwp, re, error)

                ! sirr = 0.0
                self%eda(i,:) = self%ed
                self%esa(i,:) = self%es

                ! sirr = sum(self%eda(i,:)) + sum(self%esa(i,:))
            else
                ! sirr = 0.0
                self%eda(i,:) = 0.0
                self%esa(i,:) = 0.0
            end if
        end do


    end subroutine sfcirr
end submodule oasim_sfcirr