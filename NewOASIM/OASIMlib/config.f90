module configuration
    use, intrinsic :: iso_fortran_env
    use oasim_common
    use yaml_types, only: type_node, type_dictionary, type_error
    use yaml, only: parse, error_length

    implicit none

    private

    type parameters
        character(len=string_length) :: atmo_file
        character(len=string_length) :: abso_file
        character(len=string_length) :: slingo_file
        real(kind=real_kind) :: integration_step_secs
        integer :: max_integration_steps
        logical :: zenith_avg
        logical :: local_time
        real(kind=real_kind) :: am
        real(kind=real_kind) :: vi
        character(len=string_length) :: bin_file
    end type parameters

    interface parameters
        module procedure :: init_from_file
    end interface parameters

    public :: parameters

contains
    type(parameters) function init_from_file(filename, error) result(params)
        implicit none
            
        character(len=*), intent(in) :: filename
        logical, intent(out) :: error

        character(len=error_length) :: error_message
        class(type_node), pointer :: root_node
        type(type_dictionary), pointer :: dict
        type(type_error), pointer :: parse_error

        error = .true.
        error_message = ''
        nullify(root_node)
        nullify(dict)
        nullify(parse_error)

        root_node => parse(filename, unit=100, error=error_message)
        if (.not. associated(root_node)) then
            write(error_unit, *) error_message
            return
        end if
        if (error_message /= '') then
            write(error_unit, *) "Error parsing configuration file."
            goto 999
        end if

        select type(root_node)
        class is(type_dictionary)
            dict => root_node%get_dictionary('input', required=.true., error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing input section of configuration file."
                goto 999
            end if

            params%atmo_file = dict%get_string('atmospheric_data_file', default='', error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing atmospheric_data_file parameter."
                goto 999
            end if

            params%abso_file = dict%get_string('absorption_data_file', default='', error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing absorption_data_file parameter."
                goto 999
            end if

            params%slingo_file = dict%get_string('cloud_slingo_file', default='', error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing cloud_slingo_file parameter."
                goto 999
            end if

            dict => root_node%get_dictionary('compute', required=.true., error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing compute section of configuration file."
                goto 999
            end if

            params%integration_step_secs = dict%get_real('integration_step_secs', default=real(0, real_kind), error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing integration_steps_secs parameter."
                goto 999
            end if

            params%max_integration_steps = dict%get_integer('max_integration_steps', default=0, error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing max_integration_steps parameter."
                goto 999
            end if

            params%zenith_avg = dict%get_logical('zenith_avg', default=.true., error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing zenith_avg parameter."
                goto 999
            end if

            params%local_time = dict%get_logical('local_time', default=.false., error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing local_time parameter."
                goto 999
            end if

            params%am = dict%get_real('am', default=real(0, real_kind), error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing am parameter."
                goto 999
            end if            

            params%vi = dict%get_real('vi', default=real(0, real_kind), error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing am parameter."
                goto 999
            end if 

            dict => root_node%get_dictionary('output', required=.true., error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing output section of configuration file."
                goto 999
            end if            

            params%bin_file = dict%get_string('bin_file', default='', error=parse_error)
            if (associated(parse_error)) then
                write(error_unit, *) "Error parsing bin_atmospheric_file parameter."
                goto 999
            end if        

        class default
            write(error_unit, *) "Wrong format of configuration file."
            goto 999
        end select

        error = .false.

        
999     if (associated(parse_error)) deallocate(parse_error)
        call root_node%finalize()
        deallocate(root_node)
    end function init_from_file
end module configuration