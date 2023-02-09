module tables
    use, intrinsic :: iso_fortran_env
    use :: oasim_common
    use :: matrix_read
    implicit none

    private

    real(kind=real_kind) :: epsilon = 1e-9

    type table
        private

        integer :: rows, cols
        logical :: filled
        real(kind=real_kind), allocatable :: row_limits(:)
        real(kind=real_kind), pointer, public :: tab(:,:)
    contains
        procedure :: from_table
        procedure :: get_rows
        procedure :: get_cols
        procedure :: get_high
        procedure :: get_low
        procedure :: print
        procedure :: finalize
        ! final :: destroy
    end type table

    interface table
        module procedure :: init_table
        module procedure :: init_from_file
    end interface table

    public :: real_kind, table

contains
        
    type(table) function init_table(lows, highs, error) result(t)
        implicit none
        
        real(kind=real_kind), intent(in) :: lows(:), highs(:)
        logical, intent(out) :: error

        integer :: input_rows, i, j

        t%filled = .false.
        error = .true.

        if (size(lows) /= size(highs)) return
        input_rows = size(lows)

        t%rows = input_rows
        do i = 1, input_rows - 1
            if (lows(i) > highs(i)) return
            if (lows(i + 1) < highs(i)) return
            if ((lows(i + 1) - highs(i)) >= epsilon) t%rows = t%rows + 1
        end do
        if (lows(input_rows) > highs(input_rows)) return

        allocate(t%row_limits(t%rows + 1))

        error = .false.
 
        t%row_limits(1) = lows(1)
        j = 2
        do i = 2, input_rows
            if ((lows(i) - highs(i - 1)) >= epsilon) then
                t%row_limits(j) = highs(i - 1)
                j = j + 1
            end if
            t%row_limits(j) = lows(i)
            j = j + 1
        end do
        t%row_limits(j) = highs(input_rows)
    end function init_table

    type(table) function init_from_file(filename, error, cols) result(t)
        implicit none

        character(len=*), intent(in) :: filename
        integer, intent(in), optional :: cols
        logical, intent(out) :: error
        
        integer :: file_unit, iostat, rows
        real(kind=real_kind), dimension(:), allocatable :: lows, highs
        real(kind=real_kind), dimension(:,:), pointer :: matrix

        error = .true.
        open(action="read", status="old", newunit=file_unit, file=filename, iostat=iostat)
        if (iostat /= 0) return

        rows = matrix_rows(file_unit)

        allocate(lows(rows))
        allocate(highs(rows))
        if (present(cols)) allocate(matrix(rows, cols))

        if (present(cols)) then
            call read_matrix(file_unit, error, lows, highs, cols, matrix)
        else
            call read_matrix(file_unit, error, lows, highs)
        end if

        t = table(lows, highs, error)

        if (present(cols)) then
            t%cols = cols
            t%filled = .true.
            t%tab => matrix
        end if

        if (allocated(lows)) deallocate(lows)
        if (allocated(highs)) deallocate(highs)

        close(file_unit)

        error = .false.
    end function init_from_file  

    real(kind=real_kind) function overlap_fraction(xmin, xmax, ymin, ymax)
        implicit none

        real(kind=real_kind), intent(in) :: xmin, xmax, ymin, ymax
        real(kind=real_kind) :: ydiam_1, first, last
        real(kind=real_kind), parameter :: real0 = 0.0

        ydiam_1 = 1 / (ymax - ymin)
        first = max(xmin, ymin)
        last = min(xmax, ymax)

        overlap_fraction = max(last - first, real0) * ydiam_1        
    end function overlap_fraction

    subroutine from_table(self, other)
        implicit none

        class(table) :: self
        type(table), intent(in) :: other

        integer :: other_rows, i, other_i
        real(kind=real_kind) :: first, other_first, last, other_last, conv_factor

        self%cols = other%get_cols()
        other_rows = other%get_rows()

        allocate(self%tab(self%rows, self%cols))
        self%filled = .true.

        self%tab = 0.0
        i = 1
        other_i = 1
        do while (i <= self%get_rows() .and. other_i <= other_rows)
            first = self%get_low(i)
            last = self%get_high(i)
            other_first = other%get_low(other_i)
            other_last = other%get_high(other_i)

            conv_factor = overlap_fraction(first, last, other_first, other_last)
            self%tab(i,:) = self%tab(i,:) + conv_factor * other%tab(other_i,:)

            if (last > other_last) then
                if (other_i <= other_rows) other_i = other_i + 1
            else
                if (i <= self%rows) i = i + 1
            end if
        end do
    end subroutine from_table

    subroutine print(self)
        implicit none

        class(table) :: self

        integer :: i

        do i = 1, self%rows
            write (*, '("[", F12.3, ", ", F12.3, "]: ", *(F12.3))') self%row_limits(i), self%row_limits(i + 1), self%tab(i,:)
        end do
    end subroutine print

    integer function get_rows(self)
        implicit none

        class(table) :: self

        get_rows = self%rows
    end function get_rows

    integer function get_cols(self)
        implicit none

        class(table) :: self

        get_cols = self%cols
    end function get_cols

    real(kind=real_kind) function get_high(self, bin)
        implicit none

        class(table) :: self
        integer, intent(in) :: bin

        get_high = self%row_limits(bin + 1)
    end function get_high

    real(kind=real_kind) function get_low(self, bin)
        implicit none

        class(table) :: self
        integer, intent(in) :: bin

        get_low = self%row_limits(bin)
    end function get_low

    subroutine finalize(self)
        implicit none

        class(table) :: self

        if(allocated(self%row_limits)) deallocate(self%row_limits)
        if(self%filled) deallocate(self%tab)
        self%filled = .false.
    end subroutine finalize

    ! subroutine destroy(self)
    !     implicit none

    !     type(table) :: self

    !     call self%finalize()
    ! end subroutine destroy
end module tables