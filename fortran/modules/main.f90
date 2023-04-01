program main

    use utility, only: fp, CROSS, POINT
    use problemsetup
    ! USE TIMESTEP, ONLY: TAKE_STEP

    implicit none

    real (fp) :: eposoffset(nDim)      ! offset position
    real (fp) :: eveloffset(nDim)        ! offset velocity
    real (fp) :: pposoffset(nDim)        ! offset position
    real (fp) :: pveloffset(nDim)        ! offset velocity

    ! Set loop variable for time step 
    integer :: n

    !SET OFFSETS
    eposoffset = (/0.0_fp, 0.0_fp, 0.0_fp/)
    eveloffset = (/0.0_fp, 0.0_fp, 0.0_fp/)

    pposoffset = (/0.0_fp, 0.0_fp, 0.0_fp/)
    pveloffset = (/0.0_fp, 0.0_fp, 0.0_fp/)


    call problemsetup_Init()

    print "(' Selected number of electrons: ',i2)", nElectrons
    print "(' Selected number of protons: ',i2)", nProtons
    print "(' Selected dt: ',i2)", dt

    call set_ics()

    time(1) = 0
    ! Fill rest of the array by integrating in time
    do n=1,nSteps-1
        call take_step(dt,mass,pos(:,:,n),mom(:,:,n),pos(:,:,n + 1),mom(:,:,n + 1))
        time(n+1) = n*dt
    end do

    call write_data()

    call deallocate_data()



    contains

    !INITIAL CONDITIONS
    subroutine set_ics()
        implicit none

        if (nElectrons == 1 .AND. nProtons == 1) then 

            epos(:,1,1) = (/0.0_fp,1.0_fp,0.0_fp/)
            evel(:,1,1) = (/14.0_fp,0.0_fp,0.0_fp/)

            ppos(:,1,1) = (/0.0_fp,-1.0_fp,0.0_fp/) * m_e / m_p 
            pvel(:,1,1) = (/-14.0_fp,0.0_fp,0.0_fp/) * m_e / m_p 
            
        else
            call quantile(epos(:,:,1), stdEpos, eposoffset)
            call quantile(evel(:,:,1), stdEpos, eveloffset)

            call quantile(ppos(:,:,1), stdEpos, pposoffset)
            call quantile(pvel(:,:,1), stdEpos, pveloffset)
        end if

    end subroutine set_ics

    
    subroutine allocate_data()
        implicit none
        allocate(time(nSteps))
        ! allocate(mass(nParticles))
        ! allocate(pos(nDim,nParticles,nSteps))
        ! allocate(mom(nDim,nParticles,nSteps))
    end subroutine allocate_data

    subroutine deallocate_data()
        implicit none
        deallocate(time)
        ! deallocate(mass)
        ! deallocate(pos)
        ! deallocate(mom)
    end subroutine deallocate_data
    subroutine write_data()
        implicit none
        print *, 'Writing file: ', outFile
        ! open(20,file = "sol.dat", status = "replace")
        open(20,file = outFile, status = "replace")
        do n =1, nSteps
            ! write(20,*) time(n), pos(:,:,n), mom(:,:,n)
        end do
        close(20)
    end subroutine write_data


end program main 