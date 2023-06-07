
program main
    use Potential_flows
    implicit none
    type(cylinder) :: my_cyl
    real :: x_value
    real, dimension(:),allocatable :: ARB_POINT, other_point
    integer :: direction ! direction will be used later when figuring out how to integrate from the leading edge stagnation line.

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

    x_value         = 1.5
    ! allocate(ARB_POINT(2))
    ARB_POINT       = [0.0 , 2.5]
    direction = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS IS A PLACEHOLDER VALUE FOR NOW. USE TO DEBUG rk4  
    other_point     = [-4.8, 2.5]
    ! |//////////////////////////////////////////////////////|
    ! |              FIGURE A WAY TO PARSE JSONS             |
    ! |   FOR NOW, JSON FILE VALUES CAN BE VIEWED ABOVE      |      
    ! |//////////////////////////////////////////////////////|

    write(*,*) "Cylinder radius:", my_cyl%cylinder_radius
    write(*,*) "Starting x:", my_cyl%x_start
    write(*,*) "x_min on plot:", my_cyl%x_lower_limit
    write(*,*) "x_max on plot:", my_cyl%x_upper_limit
    write(*,*) "Space step:", my_cyl%delta_s
    write(*,*) "Number of streamlines:", my_cyl%n_lines                  
    write(*,*) "Y distance between lines at starting x:", my_cyl%delta_y
    write(*,*) "Leading edge x location:", my_cyl%x_leading_edge
    write(*,*) "Trailing edge x location:", my_cyl%x_trailing_edge
    write(*,*) "x center of the geometry:", my_cyl%x_cent
    write(*,*) "V infinity:", my_cyl%vel_inf
    write(*,*) "Angle of attack:", my_cyl%alpha
    write(*,*) "Vortex strength:", my_cyl%vortex_strength
    write(*,*) "x_value:", x_value
    write(*,*) "ARB_POINT:", ARB_POINT

    write(*,*) "[Upper, Lower, Camber]: ", my_cyl%calc_geometry(x_value) ! call is only used for subroutines.
   
    call my_cyl%normal_vec(x_value)
    write(*,*) "Upper unit normal vector at x_value:", my_cyl%unit_normal_upper
    write(*,*) "Lower unit normal vector at x_value:", my_cyl%unit_normal_lower

    call my_cyl%tangent_vec(x_value)
    write(*,*) "Upper unit tangent vector at x_value:", my_cyl%unit_tangent_upper
    write(*,*) "Lower unit tangent vector at x_value:", my_cyl%unit_tangent_lower

    write(*,*) "Cartesian velocity at ARB_POINT:", my_cyl%calc_vel(ARB_POINT)

    write(*,*) "[Upper_surf_tan_vel, Lower_surf_tan_vel]: ", my_cyl%calc_surf_tan(ARB_POINT)

    write(*,*) "Vel norm by mag: ", my_cyl%calc_vel_norm(ARB_POINT)
   
    write(*,*) "Cylinder Stagnation Points: ", my_cyl%calc_stag()

    write(*,*) "Rk4 next step: ", my_cyl%calc_rk4(ARB_POINT, 1)

    ! write(*,*) "Arb_point", ARB_POINT

    write(*, '(1(2F10.2))')  my_cyl%calc_streamline(ARB_POINT)

end program main


