PROGRAM MAIN

    ! DECLARATIONS

        IMPLICIT NONE

        INTEGER , PARAMETER :: fp = SELECTED_REAL_KIND( 15 , 307 )

        REAL(fp) :: m_e , m_p , q_e , q_p
        REAL(fp) :: k , G , M
        REAL(fp) , DIMENSION(3) :: Fe , Fg , Fm
        
        REAL(fp) , DIMENSION(3,1) :: Electron_pos , Electron_vel , Electron_acc , Electron_acc2
        REAL(fp) , DIMENSION(3,1) :: Proton_pos , Proton_vel , Proton_acc , Proton_acc2


        REAL(fp) :: T , dT
        REAL(fp) :: dx , dy , dz , r

        INTEGER :: n , count0 , count1 , count2

    !

    ! DEFINITIONS

        m_e = 9.1093837015E-31_fp
        m_p = 1.67262192369E-27_fp

        q_e = -1.602176634E-19_fp
        q_p = 1.602176634E-19_fp

        k = 8.987551782E9_fp
        G = 6.6743E-11_fp
        M = 1.0E-7_fp

        T = 1E+0_fp
        dT = 1E-6_fp

        n = NINT(T/dT)

    !

    OPEN( unit = 1 , file = "nbody_e.txt" )
    OPEN( unit = 2 , file = "nbody_p.txt" )

    Electron_pos = QUANTILE( Electron_pos , 0.0_fp )
    Electron_vel = QUANTILE( Electron_vel , 0.0_fp )

    Electron_pos(2,1) = 1.0_fp
    Electron_vel(1,1) = 14.0_fp

    Proton_pos = QUANTILE( Proton_pos , 0.0_fp )
    Proton_vel = QUANTILE( Proton_vel , 0.0_fp )

    Proton_pos(2,1) = -1.0_fp * m_e / m_p
    Proton_vel(1,1) = 14.0_fp * m_e / m_p

    DO count0 = 0 , n

        IF ( MOD(count0,n/1000) == 0 ) THEN

            PRINT * , count0

        END IF

        WRITE( unit = 1 , fmt = * ) ( electron_pos(1,count1) , electron_pos(2,count1) , electron_pos(3,count1) , &
        & count1 = 1 , SIZE( Electron_pos , 2 ) )
        WRITE( unit = 2 , fmt = * ) ( proton_pos(1,count1) , proton_pos(2,count1) , proton_pos(3,count1) , &
        & count1 = 1 , SIZE( Proton_pos , 2 ) )

        Electron_acc(:,:) = 0._fp
        Electron_acc2(:,:) = 0._fp
        Proton_acc(:,:) = 0._fp
        Proton_acc2(:,:) = 0._fp

        DO count1 = 1 , SIZE( Electron_pos , 2 )

            DO count2 = 1 , count1 - 1

                dx = Electron_pos(1,count1) - Electron_pos(1,count2)
                dy = Electron_pos(2,count1) - Electron_pos(2,count2)
                dz = Electron_pos(3,count1) - Electron_pos(3,count2)

                r = sqrt( dx**2 + dy**2 + dz**2 )

                Fe = k * q_e * q_e / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_e * m_e / r**2 * POINT( dx , dy , dz )
                Fm = M * q_e * q_e / r**2 * CROSS( Electron_vel(:,count2) , &
                    & CROSS( Electron_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Electron_acc(:,count1) = Electron_acc(:,count1) + ( Fe + Fg + Fm ) / m_e

                Electron_acc(:,count2) = Electron_acc(:,count2) - ( Fe + Fg + Fm ) / m_e

            END DO

            DO count2 = 1 , SIZE( Proton_pos , 2 )

                dx = Electron_pos(1,count1) - Proton_pos(1,count2)
                dy = Electron_pos(2,count1) - Proton_pos(2,count2)
                dz = Electron_pos(3,count1) - Proton_pos(3,count2)

                r = sqrt( dx**2 + dy**2 + dz**2 )

                Fe = k * q_e * q_p / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_e * m_p / r**2 * POINT( dx , dy , dz )
                Fm = M * q_e * q_p / r**2 * CROSS( Proton_vel(:,count2) , &
                    & CROSS( Electron_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Electron_acc(:,count1) = Electron_acc(:,count1) + ( Fe + Fg + Fm ) / m_e

                Proton_acc(:,count2) = Proton_acc(:,count2) - ( Fe + Fg + Fm ) / m_p

            END DO

        END DO

        DO count1 = 1 , SIZE( Proton_pos , 2 )

            DO count2 = 1 , count1 - 1

                dx = Proton_pos(1,count1) - Proton_pos(1,count2)
                dy = Proton_pos(2,count1) - Proton_pos(2,count2)
                dz = Proton_pos(3,count1) - Proton_pos(3,count2)
    
                r = sqrt( dx**2 + dy**2 + dz**2 )
    
                Fe = k * q_p * q_p / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_p * m_p / r**2 * POINT( dx , dy , dz )
                Fm = M * q_p * q_p / r**2 * CROSS( Proton_vel(:,count2) , &
                    & CROSS( Proton_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Proton_acc(:,count1) = Proton_acc(:,count1) + ( Fe + Fg + Fm ) / m_p

                Proton_acc(:,count2) = Proton_acc(:,count2) - ( Fe + Fg + Fm ) / m_p

            END DO

        END DO

        Electron_pos = Electron_pos + Electron_vel * dT + 0.5 * Electron_acc * dT**2

        Proton_pos = Proton_pos + Proton_vel * dT + 0.5 * Proton_acc * dT**2

        DO count1 = 1 , SIZE( Electron_pos , 2 )

            DO count2 = 1 , count1 - 1

                dx = Electron_pos(1,count1) - Electron_pos(1,count2)
                dy = Electron_pos(2,count1) - Electron_pos(2,count2)
                dz = Electron_pos(3,count1) - Electron_pos(3,count2)

                r = sqrt( dx**2 + dy**2 + dz**2 )

                Fe = k * q_e * q_e / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_e * m_e / r**2 * POINT( dx , dy , dz )
                Fm = M * q_e * q_e / r**2 * CROSS( Electron_vel(:,count2) , &
                    & CROSS( Electron_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Electron_acc2(:,count1) = Electron_acc2(:,count1) + ( Fe + Fg + Fm ) / m_e

                Electron_acc2(:,count2) = Electron_acc2(:,count2) - ( Fe + Fg + Fm ) / m_e

            END DO

            DO count2 = 1 , SIZE( Proton_pos , 2 )

                dx = Electron_pos(1,count1) - Proton_pos(1,count2)
                dy = Electron_pos(2,count1) - Proton_pos(2,count2)
                dz = Electron_pos(3,count1) - Proton_pos(3,count2)

                r = sqrt( dx**2 + dy**2 + dz**2 )

                Fe = k * q_e * q_p / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_e * m_p / r**2 * POINT( dx , dy , dz )
                Fm = M * q_e * q_p / r**2 * CROSS( Proton_vel(:,count2) , &
                    & CROSS( Electron_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Electron_acc2(:,count1) = Electron_acc2(:,count1) + ( Fe + Fg + Fm ) / m_e

                Proton_acc2(:,count2) = Proton_acc2(:,count2) - ( Fe + Fg + Fm ) / m_p

            END DO

        END DO

        DO count1 = 1 , SIZE( Proton_pos , 2 )

            DO count2 = 1 , count1 - 1

                dx = Proton_pos(1,count1) - Proton_pos(1,count2)
                dy = Proton_pos(2,count1) - Proton_pos(2,count2)
                dz = Proton_pos(3,count1) - Proton_pos(3,count2)
    
                r = sqrt( dx**2 + dy**2 + dz**2 )
    
                Fe = k * q_p * q_p / r**2 * POINT( dx , dy , dz )
                Fg = -G * m_p * m_p / r**2 * POINT( dx , dy , dz )
                Fm = M * q_p * q_p / r**2 * CROSS( Proton_vel(:,count2) , &
                    & CROSS( Proton_vel(:,count1) , (/dx,dy,dz/)/r ) &
                & )

                Proton_acc2(:,count1) = Proton_acc2(:,count1) + ( Fe + Fg + Fm ) / m_p

                Proton_acc2(:,count2) = Proton_acc2(:,count2) - ( Fe + Fg + Fm ) / m_p

            END DO

        END DO

        Electron_vel = Electron_vel + 0.5 * ( Electron_acc + Electron_acc2 ) * dT

        Proton_vel = Proton_vel + 0.5 * ( Proton_acc + Proton_acc2 ) * dT

    END DO

    CONTAINS
    
        FUNCTION QUANTILE( arg , std_dev )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(SIZE(Electron_pos,1),SIZE(Electron_pos,2)) :: QUANTILE , arg
        
            REAL(fp) :: std_dev

            CALL RANDOM_NUMBER(arg)
        
            QUANTILE = std_dev * ( 5*LOG(arg/(1-arg))/12 + 0.9 * ( arg - 0.5 ) )
        
        END FUNCTION QUANTILE

        FUNCTION POINT( diff_x , diff_y , diff_z )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3) :: POINT
            
            REAL(fp) :: diff_x , diff_y , diff_z

            POINT(1) = cos(atan2( sqrt( diff_y**2 + diff_z**2 ) , diff_x ))
            POINT(2) = cos(atan2( sqrt( diff_x**2 + diff_z**2 ) , diff_y ))
            POINT(3) = cos(atan2( sqrt( diff_x**2 + diff_y**2 ) , diff_z ))

        END FUNCTION POINT

        FUNCTION CROSS( a , b )

            REAL(fp) , DIMENSION(3) :: CROSS , a , b

            CROSS(1) = a(2)*b(3) - a(3)*b(2)
            CROSS(2) = a(3)*b(1) - a(1)*b(3)
            CROSS(3) = a(1)*b(2) - a(2)*b(1)
            
        END FUNCTION CROSS

    !

END PROGRAM MAIN
