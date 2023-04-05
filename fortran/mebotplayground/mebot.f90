PROGRAM MAIN

    ! DECLARATIONS
    
        IMPLICIT NONE

        INTEGER , PARAMETER :: fp = SELECTED_REAL_KIND( 15 , 307 )
        INTEGER , PARAMETER :: ip = SELECTED_INT_KIND(9)

        INTEGER(ip) :: n , count0 , count1 , count2

        INTEGER(ip) , PARAMETER :: ne = 1 , np = 1

        REAL(fp) , DIMENSION(3,ne) :: e_r , e_v , e_a , e_a2 , e_a3 , e_a4 , e_a5 , e_a6 , e_a7 , e_a8

        REAL(fp) , DIMENSION(3,np) :: p_r , p_v , p_a , p_a2 , p_a3 , p_a4 , p_a5 , p_a6 , p_a7 , p_a8

        REAL(fp) :: T , dT

        REAL(fp) , DIMENSION(3) :: Fe , Fg , Fm
        REAL(fp) , DIMENSION(3) :: r

        REAL(fp) :: KE_e , KE_p , U_g , U_e , U_m

    !

    ! DEFINITIONS

        REAL(fp) , PARAMETER :: pi = 4*atan(1.0_fp) , c = 299792458 , k = c**2*1.0E-7_fp , G = 6.67408E-11_fp , M = 1.0E-7_fp

        REAL(fp) , PARAMETER :: m_e = 9.1093837015E-31_fp , m_p = 1.67262192369E-27_fp
        REAL(fp) , PARAMETER :: q_e = -1.6021766208E-19_fp , q_p = 1.6021766208E-19_fp

        T = 1.0E+0_fp
        dT = 1.0E-4_fp

        n = NINT(T/dT)

    !

    OPEN( unit = 1 , file = "nbody_e.txt" )
    OPEN( unit = 2 , file = "nbody_p.txt" )
    OPEN( unit = 3 , file = "energy.txt" )

    CALL QUANTILE( e_r , 0.0_fp , (/1.0_fp,0.0_fp,0.0_fp/) )
    CALL QUANTILE( e_v , 0.0_fp , (/0.0_fp,14.0_fp,0.0_fp/) )
    
    CALL QUANTILE( p_r , 0.0_fp , (/-1.0_fp,0.0_fp,0.0_fp/)*m_e/m_p )
    CALL QUANTILE( p_v , 0.0_fp , (/0.0_fp,-14.0_fp,0.0_fp/)*m_e/m_p )

    DO count0 = 1 , n

        WRITE( unit = 1 , fmt = * ) ( e_r(:,count1) , count1 = 1 , ne )
        WRITE( unit = 2 , fmt = * ) ( p_r(:,count1) , count1 = 1 , np )

        KE_e = SUM( m_e/2 * NORM2(e_v,1)**2 )
        KE_p = SUM( m_p/2 * NORM2(p_v,1)**2 )

        U_g = 0
        U_e = 0
        U_m = 0

        DO count1 = 1 , ne

            DO count2 = 1 , count1 - 1

                r = e_r(:,count1) - e_r(:,count2)

                U_g = U_g - G * m_e * m_e / NORM2(r)**2
                U_e = U_e + k * q_e * q_e / NORM2(r)**2
                U_m = U_m + DOT_PRODUCT( -9.2847647043E-24_fp*UNIT(e_v(:,count1)) , M*q_e/NORM2(r)**2*CROSS( e_v(:,count2) , r ) )

            END DO

            DO count2 = 1 , np

                r = e_r(:,count1) - p_r(:,count2)

                U_g = U_g - G * m_e * m_p / NORM2(r)**2
                U_e = U_e + k * q_e * q_p / NORM2(r)**2
                U_m = U_m + DOT_PRODUCT( -9.2847647043E-24_fp*UNIT(e_v(:,count1)) , M*q_p/NORM2(r)**2*CROSS( p_v(:,count2) , r ) )

            END DO

        END DO

        DO count1 = 1 , np

            DO count2 = 1 , count1 - 1

                r = p_r(:,count1) - p_r(:,count2)

                U_g = U_g - G * m_p * m_p / NORM2(r)**2
                U_e = U_e + k * q_p * q_p / NORM2(r)**2
                U_m = U_m + DOT_PRODUCT( 1.41060679736E-26_fp*UNIT(p_v(:,count1)) , M*q_p/NORM2(r)**2*CROSS( p_v(:,count2) , r ) )

            END DO

        END DO

        WRITE( unit = 3 , fmt = * ) KE_e , KE_p , U_g , U_e , U_m

        e_a = 0.0_fp
        e_a2 = 0.0_fp
        e_a3 = 0.0_fp
        e_a4 = 0.0_fp
        e_a5 = 0.0_fp
        e_a6 = 0.0_fp
        e_a7 = 0.0_fp
        e_a8 = 0.0_fp

        p_a = 0.0_fp
        p_a2 = 0.0_fp
        p_a3 = 0.0_fp
        p_a4 = 0.0_fp
        p_a5 = 0.0_fp
        p_a6 = 0.0_fp
        p_a7 = 0.0_fp
        p_a8 = 0.0_fp

        CALL CALC( e_r , e_v , e_a , p_r , p_v , p_a )
        CALL CALC( e_r + 0.25 * e_a * dT**2 , e_v , e_a2 , p_r + 0.25 * p_a * dT**2 , p_v , p_a2 )
        CALL CALC( e_r + 0.25 * e_a2 * dT**2 , e_v , e_a3 , p_r + 0.25 * p_a2 * dT**2 , p_v , p_a3 )
        CALL CALC( e_r + 0.25 * e_a3 * dT**2 , e_v , e_a4 , p_r + 0.25 * p_a3 * dT**2 , p_v , p_a4 )

        e_r = e_r + e_v * dT + ( e_a + 2*e_a2 + 2*e_a3 + e_a4 ) * dT**2 / 12

        p_r = p_r + p_v * dT + ( p_a + 2*p_a2 + 2*p_a3 + p_a4 ) * dT**2 / 12

        CALL CALC( e_r , e_v , e_a5 , p_r , p_v , p_a5 )
        CALL CALC( e_r , e_v + 0.5 * e_a5 * dT , e_a6 , p_r , p_v + 0.5 * p_a5 * dT , p_a6 )
        CALL CALC( e_r , e_v + 0.5 * e_a6 * dT , e_a7 , p_r , p_v + 0.5 * p_a6 * dT , p_a7 )
        CALL CALC( e_r , e_v + 0.5 * e_a7 * dT , e_a8 , p_r , p_v + 0.5 * p_a7 * dT , p_a8 )

        e_v = e_v + ( e_a + 2*e_a2 + 2*e_a3 + e_a4 + e_a5 + 2*e_a6 + 2*e_a7 + e_a8 ) * dT / 12

        p_v = p_v + ( p_a + 2*p_a2 + 2*p_a3 + p_a4 + p_a5 + 2*p_a6 + 2*p_a7 + p_a8 ) * dT / 12
    
    END DO

    CONTAINS

        FUNCTION UNIT( vector )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3) :: UNIT , vector

            UNIT(1) = cos(atan2( sqrt( vector(2)**2 + vector(3)**2 ) , vector(1) ))
            UNIT(2) = cos(atan2( sqrt( vector(1)**2 + vector(3)**2 ) , vector(2) ))
            UNIT(3) = cos(atan2( sqrt( vector(1)**2 + vector(2)**2 ) , vector(3) ))

        END FUNCTION UNIT

        FUNCTION CROSS( a , b )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3) :: CROSS , a , b

            CROSS(1) = a(2)*b(3) - a(3)*b(2)
            CROSS(2) = a(3)*b(1) - a(1)*b(3)
            CROSS(3) = a(1)*b(2) - a(2)*b(1)
            
        END FUNCTION CROSS

        SUBROUTINE QUANTILE( array , width , offset )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(:,:) :: array

            REAL(fp) :: width

            REAL(fp) , DIMENSION(3) :: offset

            CALL RANDOM_NUMBER(array)

            array = width * ( 5*LOG(array/(1-array))/12 + 0.9*( array - 0.5 ) )

            array(1,:) = array(1,:) + offset(1)
            array(2,:) = array(2,:) + offset(2)
            array(3,:) = array(3,:) + offset(3)

        END SUBROUTINE QUANTILE

        SUBROUTINE CALC( e_pos , e_vel , e_acc , p_pos , p_vel , p_acc )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3,ne) :: e_pos , e_vel , e_acc
            REAL(fp) , DIMENSION(3,np) :: p_pos , p_vel , p_acc

            DO count1 = 1 , ne

                DO count2 = 1 , count1 - 1

                    r = e_pos(:,count1) - e_pos(:,count2)

                    Fe = k * q_e * q_e / NORM2(r)**2 * UNIT(r)
                    Fg = -G * m_e * m_e / NORM2(r)**2 * UNIT(r)
                    Fm = M * q_e * q_e * CROSS( e_vel(:,count1) , CROSS( e_vel(:,count2) , r ) )

                    e_acc(:,count1) = e_acc(:,count1) + ( Fe + Fg + Fm ) / m_e

                    e_acc(:,count2) = e_acc(:,count2) - ( Fe + Fg + Fm ) / m_e

                END DO

                DO count2 = 1 , np

                    r = e_pos(:,count1) - p_pos(:,count2)

                    Fe = k * q_e * q_p / NORM2(r)**2 * UNIT(r)
                    Fg = -G * m_e * m_p / NORM2(r)**2 * UNIT(r)
                    Fm = M * q_e * q_p * CROSS( e_vel(:,count1) , CROSS( p_vel(:,count2) , r ) )

                    e_acc(:,count1) = e_acc(:,count1) + ( Fe + Fg + Fm ) / m_e

                    p_acc(:,count2) = p_acc(:,count2) - ( Fe + Fg + Fm ) / m_p

                END DO

            END DO

            DO count1 = 1 , np

                DO count2 = 1 , count1 - 1

                    r = p_pos(:,count1) - p_pos(:,count2)

                    Fe = k * q_p * q_p / NORM2(r)**2 * UNIT(r)
                    Fg = -G * m_p * m_p / NORM2(r)**2 * UNIT(r)
                    Fm = M * q_p * q_p * CROSS( p_vel(:,count1) , CROSS( p_vel(:,count2) , r ) )

                    p_acc(:,count1) = p_acc(:,count1) + ( Fe + Fg + Fm ) / m_p

                    p_acc(:,count2) = p_acc(:,count2) - ( Fe + Fg + Fm ) / m_p

                END DO

            END DO

        END SUBROUTINE CALC

    !

END PROGRAM MAIN
