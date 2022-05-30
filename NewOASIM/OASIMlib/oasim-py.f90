type(c_ptr) function py_init_lib(n, m, filename, lat, lon) bind(c) result(ptr)
    use :: iso_c_binding
    use :: oasim

    implicit none

    integer(c_int), value, intent(in) :: n, m
    character(c_char), dimension(n), intent(in) :: filename
    real(real_kind), dimension(m), intent(in) :: lat, lon

    logical :: error_i
    character(len=n) :: filename_i
    integer :: i

    type(oasim_lib), pointer :: lib

    do i = 1, n
        filename_i(i:i) = filename(i)
    end do

    allocate(lib)

    lib = oasim_lib(filename_i, lat, lon, error_i)

    if (error_i) then
        ptr = c_null_ptr
    else
        ptr = c_loc(lib)
    end if
end function py_init_lib

subroutine py_finalize_lib(ptr) bind(c)
    use :: iso_c_binding
    use :: oasim

    implicit none

    type(c_ptr), value, intent(in) :: ptr

    type(oasim_lib), pointer :: lib

    call c_f_pointer(ptr, lib)
    call lib%finalize()
    if (associated(lib)) deallocate(lib)
    nullify(lib)
end subroutine py_finalize_lib

type(c_ptr) function py_init_calc(p_size, lib_ptr) bind(c) result(ptr)
    use :: iso_c_binding
    use :: oasim

    implicit none

    integer(c_int), value, intent(in) :: p_size
    type(c_ptr), value, intent(in) :: lib_ptr

    type(calc_unit), pointer :: calc
    type(oasim_lib), pointer :: lib

    call c_f_pointer(lib_ptr, lib)
    allocate(calc)

    calc = calc_unit(p_size, lib)

    ptr = c_loc(calc)
end function py_init_calc

logical(c_bool) function py_monrad(ptr, p_size, rows, points, iyr, iday, sec_b, sec_e, &
                                  sp, msl, ws10, tco3, t2m, d2m, &
                                  tcc, tclw, cdrem, taua, asymp, ssalb, edout, esout) bind(c) result(error)
    use :: iso_c_binding
    use :: oasim

    type(c_ptr), value, intent(in) :: ptr
    integer(c_int), value, intent(in) :: rows, p_size
    integer(c_int), dimension(p_size), intent(in) :: points
    integer(c_int), value, intent(in) :: iyr, iday
    real(kind=real_kind), value, intent(in) :: sec_b, sec_e
    real(kind=real_kind), dimension(p_size), intent(in) :: sp, msl, ws10, tco3, t2m, d2m, tcc, tclw, cdrem
    real(kind=real_kind), dimension(p_size, rows), intent(in) :: taua, asymp, ssalb
    real(kind=real_kind), dimension(p_size, rows), intent(out) :: edout, esout

    type(calc_unit), pointer :: calc
    logical :: error_i

    call c_f_pointer(ptr, calc)
    call calc%monrad(points, iyr, iday, sec_b, sec_e, sp, msl, ws10, tco3, t2m, d2m, & 
                     tcc, tclw, cdrem, taua, asymp, ssalb, edout, esout, error_i)

    error = error_i
end function py_monrad

subroutine py_finalize_calc(ptr) bind(c)
    use :: iso_c_binding
    use :: oasim

    implicit none

    type(c_ptr), value, intent(in) :: ptr

    type(calc_unit), pointer :: calc

    call c_f_pointer(ptr, calc)
    call calc%finalize()
    if (associated(calc)) deallocate(calc)
    nullify(calc)
end subroutine py_finalize_calc