submodule (oasim) oasim_sfcsolz
    implicit none

contains
    module subroutine sfcsolz(self, iyr, iday, sec_b, sec_e, points)
        implicit none

        class(calc_unit) :: self
        integer, intent(in) :: iyr, iday
        real(kind=real_kind), intent(in) :: sec_b, sec_e
        integer, dimension(:), intent(in) :: points

        real(kind=real_kind) :: delsec, deltat, rs, sec_c
        integer :: steps, tstep

        if (self%lib%init_parameters%zenith_avg) then
            delsec = sec_e - sec_b
            steps = int(delsec / self%lib%init_parameters%integration_step_secs)
            steps = min(steps, self%lib%init_parameters%max_integration_steps)
            deltat = delsec / real(steps, real_kind)        

            sec_c = sec_b

            call self%sunmod(iday, iyr, sec_c, points, rs)
            self%integ = 0.5 * cos(self%sunz * rad_1)
            sec_c = sec_c + deltat

            do tstep = 1, steps - 1
                call self%sunmod(iday, iyr, sec_c, points, rs)
                self%integ = self%integ + cos(self%sunz * rad_1)
                sec_c = sec_c + deltat
            end do

            call self%sunmod(iday, iyr, sec_c, points, rs)
            self%integ = self%integ + 0.5 * cos(self%sunz * rad_1)

            self%solz = rad * acos(self%integ * deltat / delsec)
        else
            call self%sunmod(iday, iyr, sec_b, points, rs)
            self%solz = self%sunz
        end if
        self%solz = min(self%solz, 90.0d0)
        self%solz = max(self%solz, 0.0d0)
    end subroutine sfcsolz
end submodule oasim_sfcsolz