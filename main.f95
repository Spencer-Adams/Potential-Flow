program main
    use Potential_flows
    implicit none
    type(cylinder) :: my_cyl
    real, dimension(:,:), allocatable :: up_stream_vals, low_stream_vals
    real, dimension(2, 600) :: geometry_values
    integer :: start_count, end_count 
    real :: count_rate, runtime

    call system_clock(start_count, count_rate)

 ! Initialize variables (including json)
! |//////////////////////////////////////////////////////|
! |              FIGURE A WAY TO PARSE JSONS             |
! |   FOR NOW, JSON FILE VALUES CAN BE VIEWED BELOW      |      
! |//////////////////////////////////////////////////////|

    my_cyl%cylinder_radius = 2.0
    my_cyl%x_start         = -5.0
    my_cyl%x_lower_limit   = -5.0
    my_cyl%x_upper_limit   = 5.0
    my_cyl%delta_s         = 0.01    
    my_cyl%delta_y         = 0.4
   
    my_cyl%x_leading_edge  = -2.0
    my_cyl%x_trailing_edge = 2.0
    my_cyl%x_cent          = 0.0

    my_cyl%vel_inf   = 10.0 
    my_cyl%number_pi = 3.1415926535897932384626433
    my_cyl%alpha     = 15.0*my_cyl%number_pi/180
    my_cyl%vortex_strength = 80.0
    my_cyl%n_lines  = 20

! |//////////////////////////////////////////////////////|
! |              FIGURE A WAY TO PARSE JSONS             |
! |   FOR NOW, JSON FILE VALUES CAN BE VIEWED ABOVE      |      
! |//////////////////////////////////////////////////////|

    geometry_values = my_cyl%calc_geom_values()

    up_stream_vals = my_cyl%calc_upper_stream_values()

    low_stream_vals = my_cyl%calc_lower_stream_values()

    call system_clock(end_count)

    runtime = real(end_count - start_count)/count_rate

    write(*,*) "Fortran__cylinder_calculations_time:", runtime, "s"

end program main


