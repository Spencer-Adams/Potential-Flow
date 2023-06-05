
program main
    use Potential_flows
    implicit none
    type(cylinder) :: my_cyl
    real :: upper_surface, lower_surface, camber 
    real :: cylinder_radius, x_start, x_lower_limit, x_upper_limit, delta_s, delta_y  
    real :: x_value, x_leading_edge, x_trailing_edge, x_cent, vel_inf, alpha, vortex_strength
    real :: number_pi
    real, dimension(:),allocatable :: ARB_POINT 
    integer :: direction, n_lines  ! direction will be used later when figuring out how to integrate from the leading edge stagnation line.

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
    ARB_POINT       = [0.0 , 2.5]
    ! direction = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THIS IS A PLACEHOLDER VALUE FOR NOW. 
    
    ! |//////////////////////////////////////////////////////|
    ! |              FIGURE A WAY TO PARSE JSONS             |
    ! |   FOR NOW, JSON FILE VALUES CAN BE VIEWED ABOVE      |       
    ! |//////////////////////////////////////////////////////|

    write(*,*) "Cylinder radius:", cylinder_radius
    write(*,*) "Starting x:", x_start
    write(*,*) "x_min on plot:", x_lower_limit
    write(*,*) "x_max on plot:", x_upper_limit
    write(*,*) "Space step:", delta_s
    write(*,*) "Number of streamlines:", n_lines                  
    write(*,*) "Y distance between lines at starting x:", delta_y
    write(*,*) "Leading edge x location:", x_leading_edge 
    write(*,*) "Trailing edge x location:", x_trailing_edge
    write(*,*) "x center of the geometry:", x_cent
    write(*,*) "V infinity:", vel_inf
    write(*,*) "Angle of attack:", alpha
    write(*,*) "Vortex strength:", vortex_strength
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
    
    ! call my_cyl%stag()
    ! write(*,*) "Stagnation points: ", my_cyl%stagnation 

end program main 


