module matrix_read
    use, intrinsic :: iso_fortran_env
    use :: oasim_common

    implicit none

    integer, parameter :: line_length = 1024
    character, parameter :: comment_char = '#'

contains
   integer function matrix_rows(file_unit) result(rows)
        implicit none

        integer, intent(in) :: file_unit

        integer :: iostat
        character(len=line_length) :: line

        rows = 0
        iostat = 0
        do
            read(file_unit, '(A)', iostat=iostat) line
            if (iostat == iostat_end) exit
            line = trim(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) /= comment_char) rows = rows + 1
        end do
        
        rewind(file_unit)
    end function matrix_rows

    subroutine read_matrix(file_unit, error, lows, highs, cols, values)
        implicit none

        integer, intent(in) :: file_unit
        real(kind=real_kind), dimension(:), intent(out) :: lows, highs
        integer, intent(in), optional :: cols
        real(kind=real_kind), dimension(:,:), intent(out), optional :: values
        logical, intent(out) :: error

        integer :: i, j, iostat
        character(len=line_length) :: line
        
        error = .true.

        i = 1
        do
            read(file_unit, '(A)', iostat=iostat) line
            if (iostat == iostat_end) exit
            if (iostat /= 0) return
            line = trim(line)
            if (line(1:1) == comment_char) cycle
            if (len_trim(line) == 0) cycle
            if (present(values) .and. present(cols)) then 
                read(line, *, iostat=iostat) lows(i), highs(i), (values(i,j), j=1, cols)
                if (iostat /= 0) return
            else
                read(line, *, iostat=iostat) lows(i), highs(i)
                if (iostat /= 0) return
            end if
            i = i + 1
        end do

        error = .false.
    end subroutine read_matrix  

end module matrix_read