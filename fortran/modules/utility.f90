module utility

    implicit none

    integer , parameter :: fp = selected_real_kind( 15 , 307 )
    integer , parameter :: maxFileLen = 50
    integer , parameter :: maxStrLen = 100

    contains

    function CROSS(a,b)
        implicit none
        real (fp), dimension(3) :: CROSS, a, b 

        CROSS(1) = a(2)*b(3) - a(3)*b(2)
        CROSS(2) = a(3)*b(1) - a(1)*b(3)
        CROSS(3) = a(1)*b(2) - a(2)*b(1)

    end function CROSS

    function POINT( diff_x , diff_y , diff_z )
        implicit none

        real (fp), dimension(3) :: POINT

        real (fp) :: diff_x , diff_y , diff_z

        POINT(1) = cos(atan2( sqrt( diff_y**2 + diff_z**2 ) , diff_x ))
        POINT(2) = cos(atan2( sqrt( diff_x**2 + diff_z**2 ) , diff_y ))
        POINT(3) = cos(atan2( sqrt( diff_x**2 + diff_y**2 ) , diff_z ))

    end function POINT

end module utility
    
