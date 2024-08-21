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













!explain assumptions on ordering




MODULE Quantum_States


	use error_handling
	implicit none

	!N-particles state, d is local physical space dimension. separable = .TRUE. if the state is separable
	type STATE

		!single particle Hilbert space dimension
		integer :: d

		!# particles
		integer :: N


		logical :: separable

		complex*16, dimension(:), allocatable :: psi

	end type

contains

	!system state = sys
	SUBROUTINE init_random_state(sys) 

		type(STATE) :: sys
		integer*4 :: total_dim
		real*8, dimension(:), allocatable    :: A, B
		real*8 :: psi_norm

		if (sys%separable .eqv. .TRUE.) then

			total_dim = sys%d * sys%N 

		else

			total_dim = sys%d ** sys%N 

		end if


		allocate(sys%PSI( total_dim ))
		allocate(A(total_dim))
		allocate(B(total_dim))

		call Random_Number(A)
		call Random_Number(B)

		sys%psi = cmplx(A,B)

		!normalize. dot_prod conjugates one of the complex vectors
		psi_norm =  SQRT(dot_product(sys%psi, sys%psi))

		sys%psi = sys%psi/psi_norm

		deallocate(A,B)


	end subroutine	



	!comment
	FUNCTION state_to_dmatrix(sys) result(density_matrix)

		type(STATE)  :: sys
		real*8 :: trace
		complex*16, dimension(:,:), allocatable :: density_matrix, transposed, diff

		!total number of elements in vector
		integer*4 :: NN, ii


		NN = size(sys%psi)
		allocate (density_matrix(NN,NN))


		!https://emeryyi.wordpress.com/2010/05/22/outer-product-of-two-vectors-in-fortran/
		!copy a vector ncopies times making it the rows (dim = 1) or the columns (dim = 2) of a matrix
		! A * B multiplies elementwise the matrices, NO matmul


		density_matrix = spread(sys%psi(1:NN), dim = 2, ncopies = NN)  *  spread( CONJG(sys%psi(1:NN)), dim = 1, ncopies = NN) 



		!check Tr rho = 1. Think about Tr rho^2 = 1 for pure states
		trace = 0

		do ii = 1, NN

			trace = trace + REAL(density_matrix(ii, ii))

		end do


		call CHECK ( debug = .FALSE., &
    			 condition_to_check = ( ABS(trace - 1) < 0.0000001  ), &
    			 stoppage = .TRUE., &
    			 error_message = 'Error: Tr(rho) \= 1', &
    			 variable = trace)


    	!check Hermitianity

    	transposed = conjg(transpose(density_matrix))

		diff = density_matrix - transposed

		do ii = 1, size(density_matrix,1)
			call CHECK ( debug = .FALSE., &
	    		condition_to_check = ( SUM(ABS(diff(ii, :))) < 0.0000001  ), &
	    		stoppage = .TRUE., &
	    		variable = density_matrix(ii, ii))
		end do

		deallocate(transposed, diff)


	end function

	


	function PARTIAL_TRACE(R, S, in_mat) result(out_mat)

		integer*4 :: R,S, ii, jj, uu
		complex*16    in_mat(R*S, R*S), reshaped_mat(R,S,R,S), out_mat(R,R)
		complex*16 :: trace

		![,] = \(  \) since Fortran 2003


		!(0.d0,0.d0) is padding, double complex. Needed to invert 1,2 and 3,4 because of column major order
		! [1,2,3,4,5,6,7,8,9] -> 	1 4 7     because columns are filled first. Just change order in [2,1]
		!							2 5 8 
		!							3 6 9

		!https://stackoverflow.com/questions/3708307/how-to-initialize-two-dimensional-arrays-in-fortran

		reshaped_mat = RESHAPE(in_mat, [R,S,R,S], [(0.d0,0.d0)], [2,1,4,3])

		!sum over indexes 2,4 to produce rho_a

		do ii = 1, R
			do jj = 1, R

				trace = 0

				do uu = 1, S

					trace = trace + reshaped_mat(ii, uu, jj, uu )

				end do

				out_mat(ii, jj) = trace

			end do
		end do


	end function



	function PARTIAL_TRACE_B(R, S, in_mat) result(out_mat)

		integer*4 :: R,S, ii, jj, uu
		complex*16    in_mat(R*S, R*S), reshaped_mat(R,S,R,S), out_mat(R,R)
		complex*16 :: trace

		reshaped_mat = RESHAPE(in_mat, [R,S,R,S], [(0.d0,0.d0)], [2,1,4,3])

		!sum over indexes 2,4 to produce rho_a

		do ii = 1, R
			do jj = 1, R

				trace = 0

				do uu = 1, S

					trace = trace + reshaped_mat(uu, ii, uu, jj )

				end do

				out_mat(ii, jj) = trace

			end do
		end do


	end function




END MODULE quantum_states





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


    !https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
    !third answer; this works because "i0" format is assumed == print with least number of digits require -> much simpler than my previous implementation

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str


end module printing















program density_test

	use error_handling
	use quantum_states
	use printing

	implicit none


	type(STATE) :: test
	real*8 :: start, stop
	character(1) :: is_separable
	integer*4    :: dim, number
	


	test%separable = .TRUE.
	open(20, file = "times_separable.txt")

	do dim = 2, 10, 2

		test%d = dim

		do number = 1, 5 * 10**4
		
			test%N = number

			call cpu_time(start)
			call INIT_RANDOM_STATE(test)
			call cpu_time(stop)

			write(20, *) dim, number, stop-start

			deallocate(test%psi)

		end do
	end do 

	close(20)


	!test%separable = .FALSE.
	!open(21, file = "times_NOT_separable.txt")

	!do dim = 2, 9, 2

	!	test%d = dim

	!	do number = 1, 8

	!		test%N = number

	!		call cpu_time(start)
	!		call INIT_RANDOM_STATE(test)
	!		call cpu_time(stop)

	!		write(21, *) dim, number, stop-start

	!		deallocate(test%psi)

	!	end do
	!end do 

	!close(21)

	

end program