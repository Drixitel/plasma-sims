

program testing

    implicit none
    integer , parameter :: fp = selected_real_kind( 15 , 307 )
    integer :: i,k,m,l
    real :: a, b
    INTEGER , PARAMETER :: ne = 1 , np = 2
    real(fp), dimension(3,ne) :: Electron_pos
    real(fp), dimension(3,np) :: Proton_pos
    real (fp), dimension(3) :: shift, shift2

    shift = (/0.0_fp,5.0_fp,0.0_fp/) 
    shift2 = (/6.0_fp,5.0_fp,0.0_fp/) 

    call quantile( Electron_pos, 0.2_fp, shift)
    call quantile( Proton_pos, 0.2_fp, shift2)


    print *, Electron_pos
    print *, "space"
    print *, Proton_pos

    contains

    subroutine quantile(arg, std_dev, offset)
        implicit none
        real (fp), dimension(:,:) :: arg
        real (fp) :: std_dev
        real (fp), dimension(3) :: offset


        call random_number(arg)

        arg = std_dev * ( 5*LOG(arg/(1-arg))/12 + 0.9 * ( arg - 0.5 ) )

        do i = 1, 3
            arg(i,:) = arg(i,:) + offset(i)
        end do

        ! arg(1,:) = arg(1,:) + offset(1)
        ! arg(2,:) = arg(2,:) + offset(2)
        ! arg(3,:) = arg(3,:) + offset(3)
        

    end subroutine quantile

  
  end program testing
  