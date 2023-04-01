module problemsetup

    use utility, only: fp, maxFileLen, maxStrLen
    implicit none
    private

    ! DECLARATIONS 

    ! INIT FILE
    character(len=maxStrLen) :: inFile
    character(len=maxStrLen), public :: runName, outFile  ! Name of the run and the output file the code should write

    ! CONSTANTS
    integer, parameter, public :: nDim = 3                ! R3 vectors 
    integer, public :: nElectrons, nProtons, nSteps    

    real (fp), public :: tFinal, dt
    real (fp), public :: stdEpos, stdEvel
    real (fp), public :: stdPpos, stdPvel



    real (fp), parameter, public :: m_e = 9.1093837015E-31_fp       ! mass: electron                
    real (fp), parameter, public :: m_p = 1.67262192369E-27_fp      ! mass: proton           
    real (fp), parameter, public :: q_e = -1.602176634E-19_fp       ! charge: electron          
    real (fp), parameter, public :: q_p = 1.602176634E-19_fp        ! charge: proton  
    
    real (fp), parameter, public :: k = 8.987551782E9_fp            ! electrical constant         
    real (fp), parameter, public :: G = 6.6743E-11_fp               ! gravitational constant
    real (fp), parameter, public :: M = 1.0E-7_fp                   ! magnetic constant

    
    ! Arrays 
    real (fp), allocatable, public :: mass(:)              ! Mass array
    real (fp), allocatable, public :: epos(:,:,:)          ! Position array
    real (fp), allocatable, public :: ppos(:,:,:)          ! Position array
    real (fp), allocatable, public :: emom(:,:,:)          ! Momentum array 
    real (fp), allocatable, public :: pmom(:,:,:)          ! Momentum array 
    real (fp), allocatable, public :: evel(:,:,:)          ! velocity array 
    real (fp), allocatable, public :: pvel(:,:,:)          ! velocity array 
    real (fp), allocatable, public :: time(:)              ! Time array

    

    ! Subroutines
    public :: problemsetup_Init                           ! Subroutine for problem, reads files
    public :: quantile
    


    
contains

    subroutine quantile(arg, std_dev, offset)
        implicit none

        real (fp), dimension(:,:) :: arg
        real (fp), dimension(nDim) :: offset
        real (fp) :: std_dev
        integer :: i 

        call random_number(arg)

        arg = std_dev * ( 5*LOG(arg/(1-arg))/12 + 0.9 * ( arg - 0.5 ) )

        do i = 1, 3
            arg(i,:) = arg(i,:) + offset(i)
        end do

        ! arg(1,:) = arg(1,:) + offset(1)
        ! arg(2,:) = arg(2,:) + offset(2)
        ! arg(3,:) = arg(3,:) + offset(3)
    end subroutine quantile


!INIT FILE SETUP 
    ! subroutine: problemsetup_Init---------------------------------------------------------------------
    ! purpose: Set module variables to values defined in an input file
    ! inputs: <none>
    ! outputs: <none>, only sets internal variables
    subroutine problemsetup_Init()
        ! Calls all 3 functions to read init files 
        implicit none

        ! Get name of input file
        call get_command_argument(1,inFile)
        print *, "Reading from ",inFile
        
        ! Fill in default values, this is overwritten by init files if different 
        nElectrons = 1
        nProtons = 1
        stdEpos = 0.0_fp
        stdEvel = 0.0_fp
        stdPpos = 0.0_fp
        stdPvel = 0.0_fp
        nSteps = 10000
        tFinal = 50.0_fp
        
        ! Read problem settings from the input file
        nElectrons = read_initFileInt('num_electrons')
        nProtons = read_initFileInt('num_protons')
        stdEpos = read_initFileReal('std_width_e_pos')
        stdEvel = read_initFileReal('std_width_e_vel')
        stdPpos = read_initFileReal('std_width_p_pos')
        stdPvel = read_initFileReal('std_width_p_vel')
        nSteps = read_initFileInt('num_steps')
        tFinal = read_initFileReal('time_final')
        dt = tFinal/nSteps
        
        ! Set the name of the run and echo it out
        call read_initFileChar('run_name',runName)
        print *, 'Running problem: ', runName
        
        ! Set the output file, note that // does string concatenation
        outFile = trim(runName) // '.dat'
    end subroutine problemsetup_Init

    ! function: read_initFileInt -------------------------------------------------------------------
    ! purpose: Pull one integer value from an input file
    ! inputs: varName -- String that names the variable, this must be first entry on a line
    ! outputs: varValue -- Integer value that will hold the result from the input file
    function read_initFileInt(varName) result(varValue)

        implicit none
        character(len=*),intent(IN) :: varName
        integer :: varValue

        integer :: i,openStatus,inputStatus
        integer :: simInitVars
        character(len=maxStrLen) :: simCharVars
        integer :: pos1,pos2

        open(unit = 11, file=inFile, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')

        do i=1,maxFileLen
        read(11, FMT = 101, IOSTAT=inputStatus) simCharVars
        pos1 = index(simCharVars,varName)
        pos2 = pos1+len_trim(varName)
        if (pos2 > len_trim(varName)) then
            read(simCharVars(pos2+1:),*)simInitVars
            varValue = simInitVars
        endif
        end do

        close(11)

    101 FORMAT(A, 1X, I5)

    end function read_initFileInt

    ! function: read_initFileReal
    ! purpose: Pull one real value from an input file
    ! inputs: varName -- String that names the variable, this must be first entry on a line
    ! outputs: varValue -- Real value that will hold the result from the input file
    function read_initFileReal(varName) result(varValue)
        ! Requires no edit 
        implicit none
        character(len=*),intent(IN) :: varName
        real (fp) :: varValue

        integer :: i,openStatus,inputStatus
        real :: simInitVars
        character(len=maxStrLen) :: simCharVars
        integer :: pos1,pos2

        open(unit = 10, file=inFile, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')

        do i=1,maxFileLen
        read(10, FMT = 100, IOSTAT=inputStatus) simCharVars
        pos1 = index(simCharVars,varName)
        pos2 = pos1+len_trim(varName)
        if (pos2 > len_trim(varName)) then
            read(simCharVars(pos2+1:),*)simInitVars
            !print*,varName,len_trim(varName)
            !print*,simCharVars
            !print*,pos1,pos2,simCharVars(pos2+1:),simInitVars;stop
            varValue = simInitVars
        endif
        end do

        close(10)

    100 FORMAT(A, 1X, F3.1)

    end function read_initFileReal

    ! subroutine: read_initFileChar
    ! purpose: Pull one string with no spaces from an input file
    ! inputs: varName -- String that names the variable, this must be first entry on a line
    ! outputs: varValue -- String that will hold the result from the input file
    subroutine read_initFileChar(varName,varValue)
        ! Requires no edit 
        implicit none
        character(len=*),intent(IN)  :: varName
        character(len=*),intent(OUT) :: varValue

        integer :: i,openStatus,inputStatus
        character(len=maxStrLen) :: simInitVars
        character(len=maxStrLen) :: simCharVars
        integer :: pos1,pos2

        open(unit = 13, file=inFile, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')

        do i=1,maxFileLen
        read(13, FMT = 103, IOSTAT=inputStatus) simCharVars
        pos1 = index(simCharVars,varName)
        pos2 = pos1+len_trim(varName)

        if (pos2 > len_trim(varName)) then
            read(simCharVars(pos2+1:),*)simInitVars
            varValue = simInitVars
        endif
        end do

        close(13)

    103 FORMAT(A, 1X, A)

    end subroutine read_initFileChar
    
end module problemsetup



