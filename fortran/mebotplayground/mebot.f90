PROGRAM MAIN

    ! DECLARATIONS

        IMPLICIT NONE

        INTEGER , PARAMETER :: fp = SELECTED_REAL_KIND( 15 , 307 )

        INTEGER , PARAMETER :: ne = 40 , np = 40

        REAL(fp) :: m_e , m_p , q_e , q_p
        REAL(fp) :: k , G , M
        REAL(fp) , DIMENSION(3) :: Fe , Fg , Fm
        
        REAL(fp) , DIMENSION(3,ne) :: Electron_pos , Electron_vel
        REAL(fp) , DIMENSION(3,ne) :: Electron_acc , Electron_acc2 , Electron_acc3 , Electron_acc4
        REAL(fp) , DIMENSION(3,ne) :: Electron_acc5 , Electron_acc6 , Electron_acc7 , Electron_acc8

        REAL(fp) , DIMENSION(3,np) :: Proton_pos , Proton_vel
        REAL(fp) , DIMENSION(3,np) :: Proton_acc , Proton_acc2 , Proton_acc3 , Proton_acc4
        REAL(fp) , DIMENSION(3,np) :: Proton_acc5 , Proton_acc6 , Proton_acc7 , Proton_acc8


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
        dT = 1E-4_fp

        n = NINT(T/dT)

    !

    OPEN( unit = 1 , file = "nbody_e.txt" )
    OPEN( unit = 2 , file = "nbody_p.txt" )

    CALL QUANTILE( Electron_pos , 0.4_fp , (/1.0_fp,0.0_fp,0.0_fp/) )
    CALL QUANTILE( Electron_vel , 0.4_fp , (/0.0_fp,20.0_fp,0.0_fp/) )

    CALL QUANTILE( Proton_pos , 0.4_fp , (/0.0_fp,0.0_fp,0.0_fp/) )
    CALL QUANTILE( Proton_vel , 0.4_fp , (/0.0_fp,0.0_fp,0.0_fp/) )

    DO count0 = 0 , n

        IF ( MOD(count0,n/1000) == 0 ) THEN

            PRINT * , count0

        END IF

        WRITE( unit = 1 , fmt = * ) ( electron_pos(1,count1) , electron_pos(2,count1) , electron_pos(3,count1) , &
        & count1 = 1 , ne )
        WRITE( unit = 2 , fmt = * ) ( proton_pos(1,count1) , proton_pos(2,count1) , proton_pos(3,count1) , &
        & count1 = 1 , np )

        Electron_acc = 0.0_fp
        Electron_acc2 = 0.0_fp
        Proton_acc = 0.0_fp
        Proton_acc2 = 0.0_fp

        CALL CALC( Electron_pos , Electron_vel , Electron_acc , &
        & Proton_pos , Proton_vel , Proton_acc )

        CALL CALC( Electron_pos + Electron_acc * dT**2 /4 , Electron_vel , Electron_acc2 , &
        & Proton_pos + Proton_acc * dT**2 /4, Proton_vel , Proton_acc2 )

        CALL CALC( Electron_pos + Electron_acc2 * dT**2 /4 , Electron_vel , Electron_acc3 , &
        & Proton_pos + Proton_acc2 * dT**2 /4, Proton_vel , Proton_acc3 )

        CALL CALC( Electron_pos + Electron_acc3 * dT**2 /4 , Electron_vel , Electron_acc4 , &
        & Proton_pos + Proton_acc3 * dT**2 /4, Proton_vel , Proton_acc4 )

        CALL CALC( Electron_pos , Electron_vel , Electron_acc5 , &
        & Proton_pos , Proton_vel , Proton_acc5 )

        CALL CALC( Electron_pos , Electron_vel + Electron_acc5 * dT/2 , Electron_acc6 , &
        & Proton_pos , Proton_vel + Proton_acc5 * dT/2 , Proton_acc6 )

        CALL CALC( Electron_pos , Electron_vel + Electron_acc6 * dT/2 , Electron_acc7 , &
        & Proton_pos , Proton_vel + Proton_acc6 * dT/2 , Proton_acc7 )

        CALL CALC( Electron_pos , Electron_vel + Electron_acc7 * dT/2 , Electron_acc8 , &
        & Proton_pos , Proton_vel + Proton_acc7 * dT/2 , Proton_acc8 )

        Electron_pos = Electron_pos + Electron_vel * dT + &
        &( Electron_acc + Electron_acc2 + Electron_acc3 + Electron_acc4 ) * dT**2 /6

        Proton_pos = Proton_pos + Proton_vel * dT + &
        &( Proton_acc + Proton_acc2 + Proton_acc3 + Proton_acc4 ) * dT**2 /6

        Electron_vel = Electron_vel + &
        &( Electron_acc5 + Electron_acc6 + Electron_acc7 + Electron_acc8 ) * dT**2 /6

        Proton_vel = Proton_vel + &
        &( Proton_acc5 + Proton_acc6 + Proton_acc7 + Proton_acc8 ) * dT**2 /6

    END DO

    CONTAINS

        FUNCTION POINT( diff_x , diff_y , diff_z )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3) :: POINT
            
            REAL(fp) :: diff_x , diff_y , diff_z

            POINT(1) = cos(atan2( sqrt( diff_y**2 + diff_z**2 ) , diff_x ))
            POINT(2) = cos(atan2( sqrt( diff_x**2 + diff_z**2 ) , diff_y ))
            POINT(3) = cos(atan2( sqrt( diff_x**2 + diff_y**2 ) , diff_z ))

        END FUNCTION POINT

        FUNCTION CROSS( a , b )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(3) :: CROSS , a , b

            CROSS(1) = a(2)*b(3) - a(3)*b(2)
            CROSS(2) = a(3)*b(1) - a(1)*b(3)
            CROSS(3) = a(1)*b(2) - a(2)*b(1)
            
        END FUNCTION CROSS

        SUBROUTINE QUANTILE( arg , std_dev , offset )

            IMPLICIT NONE

            REAL(fp) , DIMENSION(:,:) :: arg
        
            REAL(fp) , DIMENSION(3) :: offset

            REAL(fp) :: std_dev

            CALL RANDOM_NUMBER(arg)

            arg = std_dev * ( 5*LOG(arg/(1-arg))/12 + 0.9 * ( arg - 0.5 ) )

            arg(1,:) = arg(1,:) + offset(1)
            arg(2,:) = arg(2,:) + offset(2)
            arg(3,:) = arg(3,:) + offset(3)

        END SUBROUTINE QUANTILE

        SUBROUTINE CALC( arg1 , arg2 , arg3 , arg4 , arg5 , arg6 )

            REAL(fp) , DIMENSION(3,ne) :: arg1 , arg2 , arg3
            REAL(fp) , DIMENSION(3,np) :: arg4 , arg5 , arg6

            DO count1 = 1 , ne

                DO count2 = 1 , count1 - 1
    
                    dx = arg1(1,count1) - arg1(1,count2)
                    dy = arg1(2,count1) - arg1(2,count2)
                    dz = arg1(3,count1) - arg1(3,count2)
    
                    r = sqrt( dx**2 + dy**2 + dz**2 )
    
                    Fe = k * q_e * q_e / r**2 * POINT( dx , dy , dz )
                    Fg = -G * m_e * m_e / r**2 * POINT( dx , dy , dz )
                    Fm = M * q_e * q_e / r**2 * CROSS( arg2(:,count2) , &
                        & CROSS( arg2(:,count1) , (/dx,dy,dz/)/r ) &
                    & )
    
                    arg3(:,count1) = arg3(:,count1) + ( Fe + Fg + Fm ) / m_e
    
                    arg3(:,count2) = arg3(:,count2) - ( Fe + Fg + Fm ) / m_e
    
                END DO
    
                DO count2 = 1 , np
    
                    dx = arg1(1,count1) - arg4(1,count2)
                    dy = arg1(2,count1) - arg4(2,count2)
                    dz = arg1(3,count1) - arg4(3,count2)
    
                    r = sqrt( dx**2 + dy**2 + dz**2 )
    
                    Fe = k * q_e * q_p / r**2 * POINT( dx , dy , dz )
                    Fg = -G * m_e * m_p / r**2 * POINT( dx , dy , dz )
                    Fm = M * q_e * q_p / r**2 * CROSS( arg5(:,count2) , &
                        & CROSS( arg2(:,count1) , (/dx,dy,dz/)/r ) &
                    & )
    
                    arg3(:,count1) = arg3(:,count1) + ( Fe + Fg + Fm ) / m_e
    
                    arg6(:,count2) = arg6(:,count2) - ( Fe + Fg + Fm ) / m_p
    
                END DO
    
            END DO
    
            DO count1 = 1 , np
    
                DO count2 = 1 , count1 - 1
    
                    dx = arg4(1,count1) - arg4(1,count2)
                    dy = arg4(2,count1) - arg4(2,count2)
                    dz = arg4(3,count1) - arg4(3,count2)
        
                    r = sqrt( dx**2 + dy**2 + dz**2 )
        
                    Fe = k * q_p * q_p / r**2 * POINT( dx , dy , dz )
                    Fg = -G * m_p * m_p / r**2 * POINT( dx , dy , dz )
                    Fm = M * q_p * q_p / r**2 * CROSS( arg5(:,count2) , &
                        & CROSS( arg5(:,count1) , (/dx,dy,dz/)/r ) &
                    & )
    
                    arg6(:,count1) = arg6(:,count1) + ( Fe + Fg + Fm ) / m_p
    
                    arg6(:,count2) = arg6(:,count2) - ( Fe + Fg + Fm ) / m_p
    
                END DO
    
            END DO

        END SUBROUTINE CALC

    !

END PROGRAM MAIN
