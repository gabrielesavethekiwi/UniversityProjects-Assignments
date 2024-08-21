
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
! Function Name: matrix_distance
! Description:  evaluates if two matrices have the same entries up to a certain threshold, with debug options derived from the subroutine Check
! Input:
!        mat_1, mat_2 = real matrices to confront
! Output:
!        error = real number, sum of absolute values of the elementwise difference of the input matrices
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



	FUNCTION matrix_distance(mat_1, mat_2, DEBUG, STOPPAGE) result(error)

		real*8, dimension(:,:) :: mat_1, mat_2
		logical :: DEBUG, STOPPAGE
		integer*4 :: ii,jj
		real*8 :: error

		error = 0

		call CHECK ( debug = DEBUG, &
    			 condition_to_check = ( (size(mat_1,1) == size(mat_2,1)) .AND. (size(mat_1,2) == size(mat_2,2)) ), &   !matrix sizes should match
    			 stoppage = STOPPAGE, &
    			 message = 'The resulting matrices can be compared', &
    			 error_message = 'Error: the different methods have produced different matrix dimensions')


    	do ii = 1, size(mat_1,1)
    		do jj = 1, size(mat_1,2)

    			error = error + ABS(mat_1(ii,jj) - mat_2(ii,jj))

    		end do
    	end do

    END FUNCTION matrix_distance


END MODULE ERROR_HANDLING





!---------------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------DOCUMENTATION-------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------------
! Module Name: naive_mult
! Contents: two inefficient functions for matrix matrix multiplication
!---------------------------------------------------------------------------------------------------------------------------------------------------
! Function Name: mult_usual
! Description: usual row by column multiplication with the common first matrix column and second matrix row index as innermost loop
! Input:
!	mat_1   = first real matrix to multiply (leftmost one)
!	mat_2   = second real matrix (rightmost one)
! Output:
!         mat_3   = resulting real matrix
!---------------------------------------------------------------------------------------------------------------------------------------------------
! Function Name: mult_cache
! Description: matrix matrix multiplication, but if first matrix has elements a_(i,j) and the second one b_(j,k) the looping is, from outermost
!              to innermost (i,j,k). When changing k, I am changing b_(j,k) by the columns: since the cache stores matrix elements in row major
!              order, they are accessed more rapidly than in the function mult_usual. The fastest order would be (j, k, i)
! Input:
!	      mat_1   = first real matrix to multiply (leftmost one)
!	      mat_2   = second real matrix (rightmost one)
! Output:
!         mat_3   = resulting real matrix
!---------------------------------------------------------------------------------------------------------------------------------------------------

MODULE NAIVE_MULT

	implicit none


contains

	!**********************************************************************************************************************
	! Usual row by column matrix multiplication. The first matrix has elements indexed by (i,j) and the second one (j,k).
	!**********************************************************************************************************************

	FUNCTION MULT_USUAL(mat_1, mat_2) result(mat_3)

	implicit none

	INTEGER*4 :: ii, jj, kk 

	REAL*8, dimension(:,:) :: mat_1, mat_2
	REAL*8, dimension(size(mat_1,1), size(mat_2,2)) :: mat_3

	do ii = 1, size(mat_1,1)
		do kk = 1, size(mat_2,2)
			do jj = 1, size(mat_1,2)

				mat_3(ii,kk) = mat_1(ii,jj)*mat_2(jj,kk) + mat_3(ii, kk)
			end do
		end do
	end do
	END FUNCTION MULT_USUAL




	!**************************************************************************************************************************************************************************************************
	! change the order of looping. Probably faster because before mat_2(jj,kk) in the innermost loop jumps rows and matrix elements are stored in the cache in row major order, so it has to be loaded slightly slower. Here mat_2(jj,kk) varies first ok k, column, probably won't miss cache. The fastest way would be to loop in the order (jj, kk, ii).
	!***************************************************************************************************************************************************************************************************

	FUNCTION MULT_CACHE(mat_1, mat_2) result(mat_3)

	implicit none

	INTEGER*4 :: ii, jj, kk

	REAL*8, dimension(:,:) :: mat_1, mat_2
	REAL*8, dimension(size(mat_1, 1), size(mat_2, 2)) :: mat_3

	do ii = 1, size(mat_1,1)
		do jj = 1, size(mat_1,2)
			do kk = 1, size(mat_2,2)

				mat_3(ii,kk) = mat_1(ii,jj)*mat_2(jj,kk) + mat_3(ii, kk)
			end do
		end do
	end do

	END FUNCTION mult_cache




END MODULE NAIVE_MULT





!---------------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------DOCUMENTATION-------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------------
! Program Name: test_performance
! Contents: given the size of two matrices the program checks whether they can be multiplied and carries the multiplication out with three methods:
!           naive_mult, cache_mult from module NAIVE_MULT and the fortran built-in matmul. Their performance is compared through the difference in
!           required CPU_TIME which is written after the matrix row and column dimension in a different text file for each method.
!
! Input: integer (row, column) dimensions for the first matrix to multiply, then another prompt asks for (row, column) of the second matrix. 
!        If you want to enable debugging run "./this_program_executable -d" with the debug flag "-d"
!
! Output: usual_matrix_multiplication.txt, cache-friendly_matrix_multiplication.txt, fortran_matrix_multiplication.txt
!         prints error messages if the program runs in debug mode and encounters some problems
!---------------------------------------------------------------------------------------------------------------------------------------------------



program Test_performance

	use NAIVE_MULT
	use ERROR_HANDLING

	implicit none

	INTEGER*4 N_1, M_1, N_2, M_2 !matrix dimensions, will receive them from user
	REAL*8 :: START, STOP !to store CPU times.

	REAL*8, dimension(:,:), allocatable :: A, B, out_usual, out_cache, out_fortran 
	!A,B are mat_1 and mat_2. Out-s are AB but changing method 


	logical :: DEBUG, STOPPAGE
	real*8 :: mat_distance_a, mat_distance_b


	! debug_flag will determine if DEBUG=YES in checkpoints via user input in command line

	CHARACTER(len = :), allocatable :: debug_flag    !deferred length objects from Fortran 2003. Character(*) assumes its length from something else so it won't work here
	integer*4 :: command_length
	call GET_COMMAND_ARGUMENT(1, length = command_length)

    allocate(character(command_length) :: debug_flag)
    
    call GET_COMMAND_ARGUMENT(1, value = debug_flag)


    !if -d for debug follows file name then the debug bool in checkpoint is true
   

    if (debug_flag == '-d') then
    	DEBUG = .TRUE.
    	STOPPAGE = .TRUE.
    else
    	DEBUG = .FALSE.
    	STOPPAGE = .FALSE.
    end if


    ! user input
    print *, '***********************************************************************************'
    print *, 'To enable debugging run the program with -d flag, as in:'
    print *, './program -d'
    print *, '***********************************************************************************'
    print *, 'Insert (row, column) dimensions of the first matrix'
    read(*, *) N_1, M_1
    print *, 'Insert (row, column) dimensions of the second matrix'
    read(*, *) N_2, M_2




    ! check if dimensions are positive

    call CHECK ( debug = DEBUG, &
    			 condition_to_check = (N_1 > 0) .AND. (M_1 > 0) .AND. (N_2 > 0) .AND. (M_2 > 0), &
    			 stoppage = STOPPAGE, &
    			 message = 'The dimensions are valid (positive integers)', &
    			 error_message = 'Error: matrix dimensions should be positive integers')



    ! check if the matrices can be multiplied: #columns_1 == rows_2

	call CHECK ( debug = DEBUG, &
    			 condition_to_check = (N_2 == M_1), &
    			 stoppage = STOPPAGE, &
    			 message = 'The matrices can be multiplied', &
    			 error_message = 'Error: the column and row of first and second matrix do not match')


	allocate(A(N_1,M_1))
	allocate(B(N_2,M_2))

	A = 0.d0
	B = 0.d0


	CALL RANDOM_NUMBER(A)  !uniform between 0,1 for each entry
	CALL RANDOM_NUMBER(B)


	allocate(out_usual(N_1, M_2))
	allocate(out_cache(N_1,M_2))
	allocate(out_fortran(N_1,M_2))


	!initialize entries with zero because the naive function use mat_3 = mat_3 +... in innermost loop

	out_usual = 0.d0
	out_cache = 0.d0
	out_fortran = 0.d0 

		 
	!data for usual matrix moltiplication
	CALL CPU_TIME(start)
	out_usual = MULT_USUAL(A,B)
	CALL CPU_TIME(stop)

	!in fortran unit = N_1, M_2, file = ... N is used to refer to that same file
	open(unit = 0, file = 'usual_matrix_multiplication.txt', position = 'append', action = 'write')
	write(0,*) N_1, M_2, stop-start
	close(0)



	!data for matrix multiplication changing innermost loop, now theoretically a bit faster
	CALL CPU_TIME(start)
	out_cache = MULT_CACHE(A,B)
	CALL CPU_TIME(stop)

	OPEN(unit = 1, file = 'cache-friendly_matrix_multiplication.txt', position = 'append', action = 'write')
	WRITE(1,*) N_1, M_2, stop-start
	CLOSE(1)
	!using 1 instead of repeating 0 was not necessary



	!fortran built in function
	CALL CPU_TIME(start)
	out_fortran = MATMUL(A,B) !fortran built in
	CALL CPU_TIME(stop)
	open(unit = 2, file = 'fortran_matrix_multiplication.txt', position = 'append', action = 'write')
	write(2,*) N_1, M_2, stop-start
	close(2)


	!post condition: check that the different multiplication methods don't obtain too different matrices (or even different shapes)


	mat_distance_a = MATRIX_DISTANCE(out_usual, out_cache, DEBUG, STOPPAGE)    !usual and cache
	mat_distance_b = MATRIX_DISTANCE(out_usual, out_fortran, DEBUG, STOPPAGE)  !usual and fortran multiplication


	call CHECK ( debug = DEBUG, &
    			 condition_to_check = (mat_distance_a < 0.000000001), &
    			 stoppage = STOPPAGE, &
    			 message = 'Same matrix with usual multiplication and cache optimized version', &
    			 error_message = 'Error: results with usual and cache optimized multiplication do not match')

    call CHECK ( debug = DEBUG, &
    			 condition_to_check = (mat_distance_b < 0.000000001), &
    			 stoppage = STOPPAGE, &
    			 message = 'Same matrix with usual multiplication and fortran one', &
    			 error_message = 'Error: results with usual and fortran multiplication do not match')


	deallocate(A)
	deallocate(B)
	deallocate(out_usual)
	deallocate(out_fortran)
	deallocate(out_cache)


END PROGRAM Test_performance