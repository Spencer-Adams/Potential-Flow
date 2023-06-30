module Potential_flows
    implicit none
    type cylinder
    ! This is a cylinder object or type.
        real :: cylinder_radius, x_start, x_lower_limit, x_upper_limit, delta_s, delta_y
        real :: x_leading_edge, x_trailing_edge, x_cent, vel_inf, alpha, vortex_strength
        real :: number_pi
        real,dimension(:),allocatable :: unit_normal_upper, unit_normal_lower, unit_tangent_upper, unit_tangent_lower
        real,dimension(:),allocatable ::  vel_cyl, stags
        real,dimension(:,:), allocatable :: stag_streams
        integer :: n_lines
        contains
        procedure :: calc_geometry => geometry
        procedure :: calc_geom_values => geom_values
        procedure :: normal_vec => normal
        procedure :: tangent_vec => tangent
        procedure :: calc_vel => velocity
        procedure :: calc_surf_tan => surf_tan_vel
        procedure :: calc_vel_norm => vel_norm_by_mag  
        procedure :: calc_stag => stagnation 
        procedure :: calc_stag_values => stagnation_values 
        procedure :: calc_rk4 => rk4 
        procedure :: calc_streamline => streamline
        procedure :: calc_upper_stream_values => upper_streamline_values
        procedure :: calc_lower_stream_values => lower_streamline_values

        end type cylinder 
   
    ! Here are the functions and subroutines 

    contains

        function geometry(this, x) result(geometry_vec)
            ! This function returns the y values of the upper surf, lower surf, and camber at given x inputs
            implicit none
            real, intent(in) :: x
            real :: camber, upper_surface, lower_surface
            real, dimension(:), allocatable :: geometry_vec
            class(cylinder), intent(inout) :: this
            camber = 0.0
            allocate(geometry_vec(3))
            upper_surface = sqrt((this%cylinder_radius**2 - x**2))
            lower_surface = -sqrt((this%cylinder_radius**2 - x**2))
            geometry_vec = [upper_surface, lower_surface, camber]
        end function geometry 

        subroutine normal(this, x)
        ! This subroutine takes in an x value, and returns the normal vectors at the surface of the cylinder at that x value.
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
        ! This subroutine takes in an x value, and returns the tangent vectors at the surface of the cylinder at that x value.
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

            ! Call velocity function and tangent_vec subroutine to populate the required arrays
            vel_cart = this%calc_vel(point)
            call this%tangent_vec(point(1))
            tan_vel_up = dot_product(vel_cart, this%unit_tangent_upper)
            tan_vel_low = dot_product(vel_cart, this%unit_tangent_lower)
            tan_vel_list = [tan_vel_up, tan_vel_low]
        end function surf_tan_vel

        function magnitude(vector) result(mag) !!!!!!!!!!!!!!!!! work on making a common module and stick stuff like this there.
        ! This function calculates and returns the magnitude of a given vector
            implicit none 
            real :: mag
            real, dimension(:), allocatable, intent(in) :: vector
            mag = sqrt(vector(1)**2 + vector(2)**2)
        end function magnitude 

        function vel_norm_by_mag(this, point) result(vel_norm)
        ! This subroutine takes a point, and returns the unit velocity vector at that point
            implicit none
            real ::  mag
            real, dimension(:), allocatable :: velocity_at_point, vel_norm
            real, dimension(:), allocatable, intent(in) :: point
            class(cylinder), intent(inout) :: this  
            velocity_at_point = this%calc_vel(point)
            mag = magnitude(velocity_at_point)
            vel_norm = velocity_at_point / mag
        end function vel_norm_by_mag

        function linspace(x_start, x_end, x_len) result(x) !!!!!!!!!!!!!!!!! work on making a common module and stick stuff like this there. 
        ! This works exactly like np.linspace in python. Notice it isn't part of the cylinder structure. It can be used anywhere. 
            real, dimension(:), allocatable :: x
            real :: dx
            real, intent(in) :: x_start, x_end
            integer :: i
            integer, intent(in) :: x_len
            ! class(cylinder), intent(inout) :: this  
            dx = (x_end-x_start)/(x_len-1)
            allocate(x(1:x_len))
            x = [(x_start + ((i-1)*dx), i=1, x_len)]    
        end function linspace 

        subroutine stagnation(this) 
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
            x_values = linspace(this%x_leading_edge, this%x_trailing_edge/2, x_num_steps)
            ! x_values_2 is used for the trailing stagnation line
            x_values_2 = linspace(this%x_leading_edge+(this%x_trailing_edge/2), this%x_trailing_edge, int(x_num_steps/2))
             
            allocate(stags(2))
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
            this%stags = [forward_stag, aft_stag]
        end subroutine stagnation 
       
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
        ! This function traces a streamline by looping the rk4 function until we get enough information to make a pretty graph.
            implicit none 
            real :: stopper_1, stopper_2
            logical :: stopper
            integer :: direction, i  
            real, dimension(:), allocatable :: point_new
            real, dimension(:), allocatable, intent(inout) :: point
            real, dimension(:,:), allocatable :: streamline_array
            class(cylinder), intent(inout) :: this
            
            if (point(1) >= this%x_leading_edge - 0.001 .and. point(1) <= this%x_cent - 0.001) then 
                direction = -1
            else 
                direction = 1
            end if
            
            allocate(streamline_array(1300,2), source=0.)
            allocate(point_new(2))
            stopper_1 = this%x_upper_limit - 0.001
            stopper_2 = this%x_lower_limit + 0.001 
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
        
        function geom_values(this) result(geoms)
        ! This function loops the geometry function to output all of the x and y values of the cylinder into a csv. 
            implicit none
            integer :: x_num_steps, i, j, k 
            real, dimension(:,:), allocatable :: geoms
            real, dimension(:), allocatable ::  geom_vec_1, geom_vec_2, geom_vec_3
            real, dimension(:), allocatable :: x_values, point
            class(cylinder), intent(inout) :: this
            
            ! file-related variables
            character(len = 256) :: file_path
            integer :: file_unit
            
            ! Get the current directory
            call getcwd(file_path)

            ! Add the file name to the path
            file_path = trim(file_path) // "/output_geom_values.csv"

            x_num_steps = int((this%x_trailing_edge-this%x_leading_edge)/this%delta_s)
            
            x_values = linspace(this%x_leading_edge, this%x_trailing_edge, x_num_steps) 

            ! Allocate geoms matrix 
            allocate(geoms(x_num_steps*3, 2), source = 0.)

            ! Open the file for writing
            open(newunit=file_unit, file=file_path, status="replace")

            do i = 1, size(x_values)
                geom_vec_1 = this%calc_geometry(x_values(i))
                point = [x_values(i), geom_vec_1(1)]
                ! fill the geoms matrix with the upper cylinder values
                geoms(i, :) = point
                !write the i'th row of the geoms matrix into the csv file 
                write(file_unit, *) geoms(i, 1), geoms(i,2)
            end do 

            do j = 1, size(x_values)
                geom_vec_2 = this%calc_geometry(x_values(j))
                point = [x_values(j), geom_vec_2(2)]
                ! fill matrix with the lower cylinder values 
                geoms(j+size(x_values), :) = point ! the j + size part is so the csv will be 2 columns and enough rows to fit upper, lower, and camber 
                !write the j'th row of the geoms matrix into the csv file 
                write(file_unit, *) geoms(j+size(x_values), 1), geoms(j+size(x_values),2)
            end do
            
            do k = 1, size(x_values)
                geom_vec_3 = this%calc_geometry(x_values(k))
                point = [x_values(k), geom_vec_2(3)]
                ! fill matrix with camber values 
                geoms(j+2*size(x_values), :) = point ! the j + 2*size part is so the csv will be 2 columns and enough rows to fit upper, lower, and camber
                !write the k'th row of the geoms matrix into the csv file 
                write(file_unit, *) geoms(j+2*size(x_values), 1), geoms(j+2*size(x_values),2)
            end do 
            ! Close the file
            close(file_unit)
            
        end function geom_values

        subroutine stagnation_values(this) 
        ! This function uses the streamline function and the stagnation function to calculate streamlines starting at the stagnation points.
            implicit none 
            integer :: i, j
            real, dimension(4) :: stag_points
            real, dimension(1) :: stag_point_forward_x, stag_point_forward_y
            real, dimension(:), allocatable :: point_forward, point_aft
            real, dimension(1) :: stag_point_aft_x, stag_point_aft_y
            real, dimension(:,:), allocatable :: stag_streamlines
            real, dimension(:,:), allocatable :: stag_forward, stag_aft
            class(cylinder), intent(inout) :: this

            ! file-related variables
            character(len = 256) :: file_path
            integer :: file_unit

            ! Get the current directory
            call getcwd(file_path)

            ! Add the file name to the path
            file_path = trim(file_path) // "/output_stag_streamlines.csv"
            
            ! Allocate arrays 
            allocate(stag_forward(400, 2), source = 0.)
            allocate(stag_aft(400, 2), source = 0.)
            allocate(stag_streamlines(1600,2), source = 0.)

            ! Open the file for writing
            open(newunit=file_unit, file=file_path, status="replace")

            ! calculate the stagnation points
            call this%calc_stag
            stag_points = this%stags

            stag_point_forward_x = stag_points(1)
            stag_point_forward_y = stag_points(2)
            point_forward = (/stag_point_forward_x, stag_point_forward_y/)
            
            stag_point_aft_x = stag_points(3)
            stag_point_aft_y = stag_points(4)
            point_aft = (/stag_point_aft_x, stag_point_aft_y/)
            
            ! calculate the stagnation streamlines starting at the forward and aft stagnation points
            stag_forward = this%calc_streamline(point_forward) 
            stag_aft = this%calc_streamline(point_aft)
            do i = 1, 384 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fix the hard coding later
                stag_streamlines(i, 1:2) = stag_forward(i,:)
                write(file_unit, *) stag_streamlines(i,1), stag_streamlines(i,2)  ! Write each data row to the CSV file
            end do

            do j = 385, 701 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fix the hard coding later
                stag_streamlines(j, 1:2) = stag_aft(j-384,:)
                write(file_unit, *) stag_streamlines(j,1), stag_streamlines(j,2)  ! Write each data row to the CSV file
                this%stag_streams = stag_streamlines
            end do  

            ! Close the csv file
            close(file_unit)   
        end subroutine stagnation_values
        
        function upper_streamline_values(this) result(streamline_vals)
        ! This function traces streamlines above the stagnation streamlines using the streamline function. 
            implicit none 
            integer :: i, j
            real :: y, q
            real, dimension(:,:), allocatable :: stag_points
            real, dimension(:,:), allocatable :: streamline_vals
            real, dimension(:,:), allocatable :: stag_forward, stag_aft
            real, dimension(:), allocatable :: point 
            class(cylinder), intent(inout) :: this

            ! csv file-related variables
            character(len = 256) :: file_path
            integer :: file_unit
        
            ! Get the current directory for csv
            call getcwd(file_path)
        
            ! Add the file name to the csv path
            file_path = trim(file_path) // "/output_streamlines_upper.csv"
        
            ! Allocate arrays
            allocate(stag_forward(400, 2), source = 0.)
            allocate(stag_aft(400, 2), source = 0.)
            allocate(streamline_vals(1300,20), source=0.)
            allocate(point(2))
            
            ! call stagnation_values subroutine so we can define which y values above which to trace streamlines
            call this%calc_stag_values()
            stag_points = this%stag_streams
            
            ! Pull last value of the forward stagnation streamline as a start loop value
            y = stag_points(383, 2)
            q = stag_points(383, 2)
        
            ! Open the csv file for writing
            open(newunit=file_unit, file=file_path, status="replace")

            i = 1
            j = 1
            
            do while (y <= 4.5)
                y = y + this%delta_y
                point = [this%x_start, y]
                ! Fill streamline vals with the i'th row 
                streamline_vals(:, i:i+1) = this%calc_streamline(point)
                do j = 1, size(streamline_vals(:,i))
                    ! write x, y to the csv
                    write(file_unit, *) streamline_vals(j, i), streamline_vals(j,i+1) ! write to a csv
                end do 
                i = i + 1
            end do

            ! Close the csv file
            close(file_unit)
        end function upper_streamline_values

        function lower_streamline_values(this) result(streamline_vals)
            implicit none 
            integer :: i, j
            real :: y, q, neg_delta_y 
            real, dimension(:,:), allocatable :: stag_points
            real, dimension(:,:), allocatable :: streamline_vals
            real, dimension(:,:), allocatable :: stag_forward, stag_aft
            real, dimension(:), allocatable :: point 
            class(cylinder), intent(inout) :: this

            ! csv file-related variables
            character(len = 256) :: file_path
            integer :: file_unit
        
            ! Get the current directory for csv
            call getcwd(file_path)
        
            ! Add the file name to the path
            file_path = trim(file_path) // "/output_streamlines_lower.csv"
        
            ! Allocate arrays
            allocate(stag_forward(400, 2), source = 0.)
            allocate(stag_aft(400, 2), source = 0.)
            allocate(streamline_vals(1300,20), source=0.)
            allocate(point(2))
            
            ! call stagnation_values subroutine so we can define which y values below which to trace streamlines
            call this%calc_stag_values()
            stag_points = this%stag_streams
            
            ! Pull last value of the forward stagnation streamline as a start loop value
            y = stag_points(383, 2)
            q = stag_points(383, 2)
        
            ! Open the csv file for writing
            open(newunit=file_unit, file=file_path, status="replace")

            i = 1
            j = 1
            neg_delta_y = -this%delta_y
            
            do while (y >= -7.5)
                y = y + neg_delta_y
                point = [this%x_start, y]
                ! Fill streamline vals with the i'th row 
                streamline_vals(:, i:i+1) = this%calc_streamline(point)
                do j = 1, size(streamline_vals(:,i))
                    ! write x, y to the csv
                    write(file_unit, *) streamline_vals(j, i), streamline_vals(j,i+1) ! write to a csv
                end do 
                i = i + 1
            end do

            ! Close the file
            close(file_unit)
        end function lower_streamline_values
        
end module Potential_flows
