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


END MODULE ERROR_HANDLING












module TI_Harmonic
	implicit none

contains
	
	!------------------------------------------------------------------------------------------
    ! Function that creates grid [-a,a] as a vector, given N and spacing dx, with a = N*dx
    ! Watch out: L = 2*N*dx is total length
    !-------------------------------------------------------------------------------------------
	FUNCTION GRID_1D(N, dx) result(grid)

		integer*4 :: N, ii    									!*4 up to 2147483647
		real*8    :: dx, a  								! double precision, 15 digits accuracy, order 10**[-308, 308]
		real*8, dimension(2*N+1) :: grid

		!interval from -a to a
		a = N*dx

		do ii = 0, 2*N
			grid(ii+1) = -a + ii*dx
		end do

	END function



	!----------------------------------------------------------------------------------------------------------------------------
    ! Subroutine (functions can't return two results) that should crete a symmetric tridiagonal matrix. Since DSTEV needs only
    ! an array with the diagonal terms and one for the subdiagonal, d and d_sub arrays get modified to contain those.
    !--------------------------------------------------------------------------------------------------------------------------
	SUBROUTINE TRI_H( d, sub_d, omega, grid)

 		
 		real*8, dimension(:), intent(OUT)  :: d, sub_d !OUT used in subroutines to not take them initialized
 		real*8 :: dx
 		real*8, dimension(:), intent(IN) :: grid        !do not modify it
 		integer*4 :: ii, omega

 		dx = grid(2) - grid(1)

 		do ii = 1, size(grid,1)
 			d(ii) = 1d0/(dx**2d0) +  (omega**2d0 * grid(ii)**2d0)/2d0 
 		end do

 		!subdiagonal has 1 less value than diagonal
 		do ii = 1, size(grid,1)-1
 			sub_d(ii) = -1d0/(2d0*dx**2d0) 
 		end do


	END subroutine


END module TI_Harmonic





MODULE PRINTING

	implicit none

contains
	

	!-------------------------------------------------------------------------------------------------------------
    ! Subroutine that prints vector in a file given the real vector and filename. CLEARS it if file exists
    !-------------------------------------------------------------------------------------------------------------
    SUBROUTINE SAVE_VECTOR(in_vec, filename)

        implicit none

        integer*8 :: ii
        character(*) :: filename
        real*8, dimension(:) :: in_vec

        open(unit = 30, file = filename, action = 'write')

        do ii = 1, size(in_vec,1)
            write(30, '(1X,F20.10)') in_vec(ii)
        end do

        close(30)

    END subroutine save_vector


    !-----------------------------------------------------------------------------------
    ! Subroutine that prints matrix in a file given the real matrix and filename
    !------------------------------------------------------------------------------------

    SUBROUTINE SAVE_MATRIX(in_mat, filename)

        implicit none

        integer(8)    :: ii, jj
        character(*)  :: filename
        real*8, dimension(:,:) :: in_mat

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

    END SUBROUTINE save_matrix


END module printing





PROGRAM continuous_ti_se

	use TI_Harmonic
	use error_handling
	use printing

	implicit none

	!omega integer to handle filenames better
	integer*4 :: N, ii, dimension, omega, a
	real*8    :: dx
	
	!declare 1D grid, d = H diagonal, sub_d = H subdiagonal
	real*8, dimension(:), allocatable :: lattice, d, sub_d

	character(50) :: filename_values
	character(50) :: filename_functions
	!grid is superfluos, included if somebody else wants to use the data without python script
	character(50) :: filename_grid        
	character(10)  :: STRING_N, STRING_OMEGA, STRING_A

	!declare dstev auxiliary variables, in caps lock
    integer                                     :: INFO, LWORK    !lwork is superfluous
    real*8, dimension(:), allocatable 			:: WORK
    !will contain orthonormal eigenvalues, dim(LDZ, size(d,1)). LDZ = size(d,1) for us
    real*8, dimension(:, :), allocatable 			:: Z
    
    !variables for debugging: trace of tridiagonal matrix shouldn't change after dstev
    real*8 :: tr_before, tr_after

	
	!ask user for input
	print *, '*************************************************************************'
	print *, 'The interval will be [-a, a], with dx = a/N'
    print *, 'Insert edge a'
    read(*, *) a
	print *, 'Insert integer N'
    read(*, *) N
    print *, 'Insert pulsation omega' 
    read(*, *) omega
    print *, "*************************************************************************"


    !convert integers to double prec. floats before /
    dx = DFLOAT(a)/DFLOAT(N)
    print *, dx, a, N


    !write filename with N, omega that varies but no trailing zeroes
    if (N < 10) then
    	WRITE(STRING_N, '(I1)') N 
    else if (N < 100) then
    	WRITE(STRING_N, '(I2)') N 
    else if (N < 1000) then
    	WRITE(STRING_N, '(I3)') N 
    else if (N < 10000) then
    	WRITE(STRING_N, '(I4)') N 
    else if (N < 100000) then
    	WRITE(STRING_N, '(I5)') N 
    end if

    if (A < 10) then
        WRITE(STRING_A, '(I1)') A 
    else if (A < 100) then
        WRITE(STRING_A, '(I2)') A 
    else if (A < 1000) then
        WRITE(STRING_A, '(I3)') A 
    else if (A < 10000) then
        WRITE(STRING_A, '(I4)') A 
    else if (A < 100000) then
        WRITE(STRING_A, '(I5)') A 
    end if

    if (OMEGA < 10) then
    	WRITE(STRING_OMEGA, '(I1)') OMEGA 
    else if (OMEGA < 100) then
    	WRITE(STRING_OMEGA, '(I2)') OMEGA 
    else if (OMEGA < 1000) then
    	WRITE(STRING_OMEGA, '(I3)') OMEGA 
    else if (OMEGA < 10000) then
    	WRITE(STRING_OMEGA, '(I4)') OMEGA 
    else if (OMEGA < 100000) then
    	WRITE(STRING_OMEGA, '(I5)') OMEGA 
    else if (OMEGA < 1000000) then
    	WRITE(STRING_OMEGA, '(I6)') OMEGA 
    else if (OMEGA < 10000000) then
    	WRITE(STRING_OMEGA, '(I7)') OMEGA 
    else if (OMEGA < 100000000) then
    	WRITE(STRING_OMEGA, '(I8)') OMEGA 
    end if


    ! // is string concatenation
    filename_values = "eigenvalues"//trim(STRING_N)//"_"//trim(STRING_A)//"_"//trim(STRING_OMEGA)//".txt"
    filename_functions = "eigenfunctions"//trim(STRING_N)//"_"//trim(STRING_A)//"_"//trim(STRING_OMEGA)//".txt"

    


	!not necessary: function GRID_1D still works w/out
	allocate(lattice(2*N+1))     
	
	!necessary: w/out invalid memory ref in TRI_H call
	allocate(d(2*N+1))
	allocate(sub_d(2*N))


	lattice = GRID_1D(N,dx)

	!diagonal and subdiagonal for stdev ready
	call TRI_H (d, sub_d, omega, lattice)

	!initialize dstev auxiliary variables, 2N+1 is the size of d
	dimension = 2*N + 1
	LWORK = max(1, 2*dimension-1)
	allocate(WORK(LWORK))
	allocate(Z(dimension, dimension))


	tr_before = sum(d)

	! exit: d = eigenvalues ascending, Z = eigenvectors
	call DSTEV('V', dimension, d, sub_d, Z, dimension, WORK, INFO)

	tr_after = sum(D)


	call CHECK( debug = .TRUE., &
                condition_to_check = ( ABS( tr_before - tr_after ) < 0.000001   ), &
                stoppage = .TRUE., &
                error_message = "Error: Tr(H) before dstev differs from Sum(Eigenvalues), which equals", &
                variable = tr_after )

    !renormalize output if needed (dstev already does it)
    !do ii = 1, dimension
    !   Z(:, ii) = Z(:, ii) / (dot_product(Z(:, ii), Z(:, ii)))
    !end do

    !normalize eigenvectors: sum dx*psi*psi = dx sum psi*psi = dx
    Z = Z * sqrt(1/dx)

    
    

    call SAVE_VECTOR(D, filename_values)
    call SAVE_MATRIX(Z, filename_functions)


    deallocate(lattice)
    deallocate(d)
	deallocate(sub_d)
	deallocate(WORK)
	deallocate(Z)

end program