!---------------------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------DOCUMENTATION-------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------------------
! Module Name: error_handling
! Contents: support module for the program Test_performances. Includes a very general subroutine for debugging and a function for 
!           checking if two matrices are the same
!---------------------------------------------------------------------------------------------------------------------------------------------------------
! Subroutine Name: Check
! Description: by passing the argument DEBUG = .TRUE. you can check whether a condition, important for the correctedness of 
!              your program is true ( and output an optional message) or false. In this case, you can output an optional error message 
!              and stop the program if STOPPAGE = .TRUE. 
!              Optionally, you can pass in a variable to check its type
! Input: 
!        debug               =  bool, if false the check subroutine is turned off
!        condition_to check  =  bool, the condition to check in the debugging procedure, i.e. matrix_rows > 0
!        verbose             =  optional bool, not yet implemented. In the future it might suppress messages or give shorter ones
!        stoppage            =  optional bool, if .TRUE. the program gets stopped if condition_to_check is false. Default is .TRUE.
!        message             =  optional character, message in the case the program works fine (condition_to_check true)
!        error_message       =  optional character, error message (condition_to_check false)
!        variable            =  optional class(*), can be whatever fortran native type. Pass it to check the type of something in your program
!---------------------------------------------------------------------------------------------------------------------------------------------------------



MODULE ERROR_HANDLING

    implicit none


    contains

    SUBROUTINE CHECK (debug, condition_to_check, verbose, stoppage, message, error_message, variable)

    implicit none

        logical, intent(IN) :: debug, condition_to_check
        logical, intent(IN), optional :: verbose, stoppage
        class(*), optional :: variable
        character(*), optional :: message, error_message

        !ALT is dummy used within the subroutine for clarity, has the role of stoppage if specified
        
        logical :: ALT

        if (.NOT.PRESENT(stoppage)) then
            ALT = .TRUE.
        else
            ALT = stoppage
        end if


        !*******************************************************************
        ! if the program behaves as it should because the condition is true
        !*******************************************************************
        
        IF (condition_to_check) then
            IF (debug) then

                if (present(message)) then
                    print *, message
                end if

                if (present(variable)) then
                    select type(variable)

                        type is (integer(2))
                            print *, variable, 'integer(2)'
                        type is (integer(4))
                            print *, variable, 'integer(4)'
                        type is (real(4))
                            print *, variable, 'real(4): single precision'
                        type is (real(8))
                            print *, variable, 'real(8): double precision'
                        type is (complex(4))
                            print *, variable, 'complex(4)'
                        type is (complex(8))
                            print *, variable, 'complex(8)'
                        type is (logical)
                            print *, variable, 'logical'
                        type is (character(*))
                            print *, variable, 'character(*)'

                    end select
                end if
            END IF
            


        !***************************************************************
        ! if the condition is false and the program has something wrong
        !***************************************************************

        else

            IF (debug) then

                if (present(error_message)) then
                    print *, error_message
                end if

                if (present(variable)) then
                    select type(variable)

                        type is (integer(2))
                            print *, variable, 'integer(2)'
                        type is (integer(4))
                            print *, variable, 'integer(4)'
                        type is (real(4))
                            print *, variable, 'real(4): single precision'
                        type is (real(8))
                            print *, variable, 'real(8): double precision'
                        type is (complex(4))
                            print *, variable, 'complex(4)'
                        type is (complex(8))
                            print *, variable, 'complex(8)'
                        type is (logical)
                            print *, variable, 'logical'
                        type is (character(*))
                            print *, variable, 'character(*)'

                    end select
                end if
            END IF


            ! if ALT is enabled (default or specified with stoppage), stop the program
            if (ALT) then 
                STOP
            end if

        END IF
            
    END SUBROUTINE CHECK



    !-----------------------------------------------------------------------------------------------------------------------
    ! Function that returns .TRUE. if TrA == sum eigenvalues. It's specific to this program with hermitian square matrix A 
    !-----------------------------------------------------------------------------------------------------------------------
    FUNCTION trace_check(A, eig) result(is_trace_invariant)

        complex*16, dimension(:, :) :: A
        real*8, dimension(:)        :: eig
        logical                     :: is_trace_invariant
        !indexes and sum
        integer                     :: ii
        real*8                      :: Tr_A, Tr_eig

        Tr_A = 0
        Tr_eig = 0

        do ii = 1, size(A,1)
            Tr_A = Tr_A + A(ii, ii)
        end do

        do ii = 1, size(eig,1)
            Tr_eig = Tr_eig + eig(ii)
        end do

        is_trace_invariant = (ABS(Tr_A-Tr_eig) < 0.00001)

    END function


END MODULE ERROR_HANDLING








MODULE HERMITIAN_RANDOM
    
    !this module doesn't know the other one, explicit interface required error otherwise
    use error_handling

    implicit none

    CONTAINS


    !-----------------------------------------------------------------------------------------------------------------------------
    ! Function should create a real diagonal square NxN matrix, given N as input. Random on [-1,1]
    ! Since it already has its eigenvalues on the diagonal, we can directly put them in a 1D array
    ! and order it (abscending), gaining computational efficiency
    !Uses Dlasrt (needs initialization of integer INFO:
    ! http://www.netlib.org/lapack/explore-html/d2/d24/group__aux_o_t_h_e_rcomputational_ga88a9fb57b6459176b68f21720d3d2ad2.html
    !-----------------------------------------------------------------------------------------------------------------------------
    FUNCTION Rand_Diagonal(N) result(eig_diagonal)

        integer*8              :: N, ii, jj, INFO
        real*8, dimension(N)   :: eig_diagonal

        !check if input dimension is valid

        call CHECK(debug = .TRUE.,&
                   condition_to_check = (N >= 1), &
                   stoppage =.TRUE., &
                   error_message = "Error: invalid diagonal random matrix dimension, it should be positive")


        do ii = 1, N
            eig_diagonal(ii) = 2*RAND(0) - 1
        end do

        !'I': orders third argument, the array, in incresing order
        call DLASRT('I', N, eig_diagonal, INFO)

    END function


    !------------------------------------------------------------------------------------------
    ! Function that creates an Hermitian square NxN matrix, given N as input. Random on [-1,1]
    !-------------------------------------------------------------------------------------------
    FUNCTION Rand_Hermitian(N) result(A)

        integer*8                   :: N, ii, jj
        complex*16, dimension(N, N) :: A


        !check if input dimension is valid

        call CHECK(debug = .TRUE.,&
                   condition_to_check = (N >= 1), &
                   stoppage =.TRUE., &
                   error_message = "Error: invalid hermitian random matrix dimension, it should be positive")



        !fill the matrix

        do ii = 1,N
            do jj = 1, ii

                if (ii == jj) then
                    A(ii, ii) = complex(RAND(0)*2-1, 0)
                else
                    A(ii, jj) = complex(RAND(0)*2 - 1, RAND(0)*2 - 1)
                    A(jj, ii) = conjg(A(ii, jj))
                end if

            end do
        end do

    end function Rand_Hermitian


    !------------------------------------------------------------------------------------------------------------------------------
    ! Function that creates the vector of spacings given real eigenvalues. Eigenvalues are in crescent order. s has size(eig,1)-1
    !------------------------------------------------------------------------------------------------------------------------------
    FUNCTION SPACINGS(eig) result(s)

    real(8), dimension(:) :: eig
    !index
    integer(8) :: ii, N

    !dummy variable for normalized spacings
    real(8) :: mean_lambda

    ! the vector of spacings has one less value than eig
    real(8), dimension(lbound(eig,dim=1):(ubound(eig,dim=1)-1)) :: s

    N = size(eig, 1)

    do ii = 1, N -1
        s(ii) = eig(ii+1) - eig(ii)
    end do

    !telescopic sum of delta lambda leaves last - first, average means divide by size(s) = N-1
    mean_lambda = (eig(N) - eig(1))/ (N-1)

    do ii = 1, N -1
        s(ii) = s(ii)/mean_lambda
    end do

    END function spacings



    !-----------------------------------------------------------------------------------
    ! Subroutine that prints matrix in a file given the complex matrix and filename
    !------------------------------------------------------------------------------------

    SUBROUTINE SAVE_cmplx_MATRIX(in_mat, filename)

        implicit none

        integer(8)    :: ii, jj
        character(*)  :: filename
        complex*16, dimension(:,:) :: in_mat

        !using unit < 20 seems to be discouraged online
        open(unit = 29, file = filename)

        DO ii = 1, size(in_mat, 1)
            DO jj = 1, size(in_mat, 2)

                write(29, '(1X,F20.10,F20.10)', advance = 'no') in_mat(ii,jj)

            END DO
            ! newline
            write(29, *) ''
        END DO

        close(29)

    END SUBROUTINE save_cmplx_matrix


    !-------------------------------------------------------------------------------------------------------------
    ! Subroutine that prints vector in a file given the real vector and filename. Appends if file already exists
    !-------------------------------------------------------------------------------------------------------------
    SUBROUTINE SAVE_VECTOR(in_vec, filename)

        implicit none

        integer*8 :: ii
        character(*) :: filename
        real*8, dimension(:) :: in_vec

        open(unit = 30, file = filename, position = 'append', action = 'write')

        do ii = 1, size(in_vec,1)
            write(30, '(1X,F20.10)') in_vec(ii)
        end do

        close(30)

    END subroutine save_vector


END MODULE hermitian_random




PROGRAM Eigenvalues
    
    use error_handling
    use hermitian_random

    implicit none


    integer*8     :: N, ii, repeat_program
    character(50) :: filename
    character(50) :: filename_diagonal
 
    complex*16, dimension(:,:), allocatable :: A
    real*8, dimension(:), allocatable       :: eig, s

    !the same vectors for the diagonal matrix
    real*8, dimension(:), allocatable       :: eig_diagonal, s_diagonal


    !declare zheev auxiliary variables, in caps lock
    integer                                     :: INFO, LWORK
    real*8, dimension(:), allocatable :: RWORK
    complex*16, dimension(:), allocatable     :: WORK

    filename = 'hermitian_spacings.txt'
    filename_diagonal = 'diagonal_spacings.txt'

    !append spacings to file for 15 different matrices
    DO repeat_program = 1, 30

        !doing 2001x2001 matrices so I have 2000 spacings
        N = 2001

        LWORK = max(1, 2*N-1)
        
        allocate(A(N,N))
        allocate(eig(N))
        allocate(RWORK(3*N-2))
        allocate(WORK(LWORK))

        !create random hermitian matrix
        A = Rand_Hermitian(N)

        !eig contains the eigenvalues
        call ZHEEV('N', 'U', N, A, N, eig, WORK, LWORK, RWORK, INFO)


        !check that TrA = sum eig = Tr diagonal A, trace should be invariant
        call CHECK( debug = .TRUE., &
                    condition_to_check = TRACE_CHECK(A, eig), &
                    stoppage = .TRUE., &
                    error_message = "Error: the sum of eigenvalues isn't equal to matrix trace")



        !calculate spacings
        s = spacings(eig)

        !save spacings to file
        call SAVE_VECTOR(s, TRIM(filename))

        deallocate(A)
        deallocate(eig)
        deallocate(RWORK)
        deallocate(WORK)
        deallocate(s)



        !part for generating diagonal random matrix spacings
        

        allocate(eig_diagonal(N))

        eig_diagonal = Rand_Diagonal(N)
        s_diagonal   = spacings(eig_diagonal)

        !check that normalized spacings sum to N-1. Switched it off because it was a trivial condition 
        call CHECK(debug = .FALSE., &
                   condition_to_check = ( ABS( (N-1) - SUM(s_diagonal) ) < 0.000001   ), &
                   stoppage = .TRUE., &
                   error_message = "Error: normalized spacings should sum to N-1, instead they sum to", &
                   variable = SUM(s_diagonal))

        call SAVE_VECTOR(s_diagonal, TRIM(filename_diagonal))
        
        deallocate(eig_diagonal)
        deallocate(s_diagonal)

    END DO 


end program eigenvalues