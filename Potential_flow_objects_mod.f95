module Potential_flows
    implicit none
    type cylinder
    ! This is a cylinder object or type.
        real :: cylinder_radius, x_start, x_lower_limit, x_upper_limit, delta_s, delta_y
        real :: x_leading_edge, x_trailing_edge, x_cent, vel_inf, alpha, vortex_strength
        real :: number_pi
        real,dimension(:),allocatable :: unit_normal_upper, unit_normal_lower, unit_tangent_upper, unit_tangent_lower
        real,dimension(:),allocatable ::  vel_cyl, forward_stag, aft_stag
        integer :: n_lines
        contains
        procedure :: calc_geometry => geometry
        procedure :: normal_vec => normal
        procedure :: tangent_vec => tangent
        procedure :: calc_vel => velocity
        procedure :: calc_surf_tan => surf_tan_vel
        procedure :: calc_vel_norm => vel_norm_by_mag  
        procedure :: calc_stag => stagnation  
        procedure :: calc_rk4 => rk4 
        procedure :: calc_streamline => streamline 

    end type cylinder 
    ! Here are the functions.

    contains
        subroutine normal(this, x)
            ! This subroutine takes in an x value, and spits out the normal vectors at the surface of the cylinder at that x value.
            implicit none 
            real :: angle_up, angle_low
            real, intent(in) :: x
            real, dimension(:), allocatable :: geometry_vec
            class(cylinder),intent(inout) :: this
            geometry_vec = this%calc_geometry(x)
            angle_up = atan2(geometry_vec(1), x)
            angle_low = atan2(geometry_vec(2), x)
            this%unit_normal_upper = [cos(angle_up),sin(angle_up)]
            this%unit_normal_lower = [cos(angle_low), sin(angle_low)]
        end subroutine normal

        subroutine tangent(this, x)
            ! This subroutine takes in an x value, and spits out the tangent vectors at the surface of the cylinder at that x value.
            implicit none 
            real :: angle_up, angle_low
            real, intent(in) :: x
            real, dimension(:), allocatable :: geometry_vec
            class(cylinder),intent(inout) :: this
            geometry_vec = this%calc_geometry(x)
            angle_up = atan2(geometry_vec(1), x)
            angle_low = atan2(geometry_vec(2), x)
            this%unit_tangent_upper = [sin(angle_up),-cos(angle_up)]
            this%unit_tangent_lower = [-sin(angle_low),cos(angle_low)]
        end subroutine tangent 

        function geometry(this, x) result(geometry_vec)
            implicit none
            real, intent(in) :: x
            real :: camber, upper_surface, lower_surface
            real, dimension(:), allocatable :: geometry_vec
            class(cylinder), intent(inout) :: this
            camber = 0.0
            upper_surface = sqrt((this%cylinder_radius**2 - x**2))
            lower_surface = -sqrt((this%cylinder_radius**2 - x**2))
            geometry_vec = [upper_surface, lower_surface, camber]
        end function geometry 

        function velocity(this, point) result(vel_cart)
        ! This function turns a cylindrical velocity at a point into a cartesian velocity at that point.
            implicit none 
            real :: angle, r
            real,dimension(:),allocatable :: vel_cart 
            real,dimension(:),allocatable, intent(in) :: point 
            class(cylinder),intent(inout) :: this

            angle = atan2(point(2), point(1))
            r = sqrt((point(1))**2 + (point(2))**2)

            this%vel_cyl = [this%vel_inf*(1-(this%cylinder_radius**2/(r**2)))*cos(angle-this%alpha), &
             -(this%vel_inf*(1+(this%cylinder_radius**2/(r**2)))*sin(angle-this%alpha)+ &
             (this%vortex_strength)/(2*this%number_pi*r))]
               
            vel_cart = [cos(angle)*this%vel_cyl(1)-sin(angle)*this%vel_cyl(2), &
            sin(angle)*this%vel_cyl(1) + cos(angle)*this%vel_cyl(2)]
        end function velocity 

        function surf_tan_vel(this, point) result(tan_vel_list)
        ! This function finds the tangent velocity at a point. This will help with finding the stagnation points in the stagnation subroutine
            implicit none
            real, dimension(:), allocatable :: tan_vel_list, vel_cart
            real, dimension(:), allocatable, intent(in) :: point
            class(cylinder), intent(inout) :: this
            real :: tan_vel_up, tan_vel_low
            ! Call calc_vel and tangent_vec subroutines to populate the required arrays
            vel_cart = this%calc_vel(point)
            call this%tangent_vec(point(1))
            tan_vel_up = dot_product(vel_cart, this%unit_tangent_upper)
            tan_vel_low = dot_product(vel_cart, this%unit_tangent_lower)
            tan_vel_list = [tan_vel_up, tan_vel_low]
        end function surf_tan_vel

        function magnitude(vector) result(mag)
            implicit none 
            real :: mag
            real, dimension(:), allocatable, intent(in) :: vector
            mag = sqrt(vector(1)**2 + vector(2)**2)
        end function magnitude 

        function vel_norm_by_mag(this, point) result(vel_norm)
        ! This subroutine simply takes a point, and returns the unit velocity vector at that point
            implicit none
            real ::  mag
            real, dimension(:), allocatable :: velocity_at_point, vel_norm
            real, dimension(:), allocatable, intent(in) :: point
            class(cylinder), intent(inout) :: this  
            velocity_at_point = this%calc_vel(point)
            mag = magnitude(velocity_at_point)
            vel_norm = velocity_at_point / mag
        end function vel_norm_by_mag

        function linspace(x_start, x_end, x_len) result(x)
            real, dimension(:), allocatable :: x
            real :: dx
            real, intent(in) :: x_start, x_end
            integer :: i
            integer, intent(in) :: x_len
            dx = (x_end-x_start)/(x_len-1)
            allocate(x(1:x_len))
            x = [(x_start + ((i-1)*dx), i=1, x_len)]    
        end function linspace 

        function stagnation(this) result(stags)
        !This subroutine returns the stagnation point of a cylinder given the velocity field properties.
            implicit none  
            integer :: x_num_steps, i, j
            real, dimension(:), allocatable :: x_values, x_values_2, point_up, point_low
            real, dimension(:), allocatable :: point_up_2, point_low_2, geom_vec, geom_vec_2
            real, dimension(:), allocatable :: surface_tan_up, surface_tan_low
            real, dimension(:), allocatable :: surface_tan_up_2, surface_tan_low_2
            real, dimension(:), allocatable :: forward_stag, aft_stag, stags
            class(cylinder), intent(inout) :: this
            x_num_steps = int((this%x_trailing_edge-this%x_leading_edge)/this%delta_s)
            ! write(*,*) x_num_steps
            x_values = linspace(this%x_leading_edge, this%x_trailing_edge/2, x_num_steps)
            x_values_2 = linspace(this%x_leading_edge+(this%x_trailing_edge/2), this%x_trailing_edge, int(x_num_steps/2))
            do i = 1, size(x_values)
                geom_vec = this%calc_geometry(x_values(i))
                point_up = [x_values(i), geom_vec(1)]
                point_low = [x_values(i), geom_vec(2)]
                surface_tan_up = this%calc_surf_tan(point_up)
                surface_tan_low = this%calc_surf_tan(point_low)
                if (surface_tan_up(1) <= 0.001 .and. surface_tan_low(2) >= 0.001) then 
                    forward_stag = [x_values(i), geom_vec(1)]
                else if (surface_tan_low(2) <= 0.001 .and. surface_tan_up(1) >= 0.001) then
                    forward_stag = [x_values(i), geom_vec(2)]
                end if
            end do 
            do j = 1, size(x_values_2)
                geom_vec_2 = this%calc_geometry(x_values_2(j))
                point_up_2 = [x_values_2(j), geom_vec_2(1)]
                point_low_2 = [x_values_2(j), geom_vec_2(2)]
                surface_tan_up_2 = this%calc_surf_tan(point_up_2)
                surface_tan_low_2 = this%calc_surf_tan(point_low_2)
                if (surface_tan_up_2(1) <= 0.001 .and. surface_tan_low_2(2) >= 0.001) then 
                    aft_stag = [x_values_2(j), geom_vec_2(1)]
                else if (surface_tan_low_2(2) <= 0.001 .and. surface_tan_up_2(1) >= 0.001) then
                    aft_stag = [x_values_2(j), geom_vec_2(2)]
                end if 
            end do 
            stags = [forward_stag, aft_stag]
        end function stagnation 
       
        function rk4(this, point, direction) result(point_new)
            ! This function takes in a starting position, and integrates it one time to get a new position using rk4 integration.
            implicit none
            real, dimension(:), allocatable :: k2_val, k3_val, k4_val
            real, dimension(:), allocatable :: point_new, k1, k2
            real, dimension(:), allocatable :: k3, k4 
            real, dimension(:), allocatable, intent(in) :: point
            integer, intent(in) :: direction
            class(cylinder), intent(inout) :: this
            k1 = this%calc_vel_norm(point)
            k2_val = (point+((direction*this%delta_s)/2)*k1)
            k2 = this%calc_vel_norm(k2_val)
            k3_val = point+((direction*this%delta_s)/2)*k2
            k3 = this%calc_vel_norm(k3_val)
            k4_val = point+((direction*this%delta_s))*k3
            k4 = this%calc_vel_norm(k4_val)
            point_new = point + ((direction*this%delta_s)/6)*(k1+2*k2+2*k3+k4)
        end function rk4 
       
        function streamline(this, point) result(streamline_array)
            implicit none 
            real :: stopper_1, stopper_2
            logical :: stopper
            integer :: direction, i  
            real, dimension(:), allocatable :: point_new
            real, dimension(:), allocatable, intent(inout) :: point
            real, dimension(:,:), allocatable :: streamline_array
            class(cylinder), intent(inout) :: this
            ! write(*,*)  "Point"
            if (point(1) >= this%x_leading_edge - 0.001 .and. point(1) <= this%x_cent - 0.001) then 
                direction = -1
            else 
                direction = 1
            end if
            allocate(streamline_array(500,2))
            allocate(point_new(2))
            ! allocate(point(2))
            stopper_1 = this%x_upper_limit - 0.001
            stopper_2 = this%x_lower_limit + 0.001    !!!!!!!!!!! why will this not work? please just why?
            i = 1
            streamline_array(i, :) = point
            stopper = .false.
            do while (.not. stopper)
                i = i + 1 
                point_new = this%calc_rk4(point, direction)
                streamline_array(i, :) = point_new
                point = point_new
                if (direction == 1 .and. point(1) >= stopper_1 .or. direction == -1 .and. point(1) <= stopper_2) then
                    stopper = .true. 
                end if    
            end do 
        end function streamline 

        ! function start(vecA, vecB) result(vecC)
        !     implicit none
        !     integer, dimension(:), allocatable :: vecA, vecB
        !     integer, dimension(:, :), allocatable :: vecC
       
        !     ! Allocate vecC with appropriate shape
        !     allocate(vecC(size(vecA) + size(vecB), 2))
       
        !     ! Copy the elements from vecA and vecB to vecC
        !     vecC(:, 1) = vecA
        !     vecC(:, 2) = vecB
       
        ! end function start

end module Potential_flows
