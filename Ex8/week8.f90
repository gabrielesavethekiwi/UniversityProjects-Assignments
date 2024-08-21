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


    !useful for filenames
    !https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
    !third answer; this works because "i0" format is assumed == print with least number of digits require -> much simpler than my previous implementation

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str


END MODULE ERROR_HANDLING



!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************



module TENSOR_PRODUCTS

	use ERROR_HANDLING

	implicit none

contains



	!-------------------------------------------------------------------------------------------------------------------------
    ! Function that takes a (complex*16) matrix A (m x n) and B (p x q) and returns their tensor product C (mp x nq) 
    ! The tensor product is A_11 * B, A_12*B, ...
    !-------------------------------------------------------------------------------------------------------------------------

	function TPROD(A,B) result(C)

		complex*16, dimension(:,:), intent(IN) :: A,B
		complex*16, dimension(:,:), allocatable:: C

		!matrix dimensions for readability
		integer*4 :: m,n,p,q
		!indexes for loops
		integer*4 :: ii, jj
		!where to start and finish assigning values to C for Aij*B varying ij (start row, finish row, ... column)
		integer*4 :: srow, frow, scol, fcol

		m = size(A,1)
		n = size(A,2)
		p = size(B,1)
		q = size(B,2)

		allocate( C(m*p, n*q))

		!loop over A to get Aij that multiplies a whole B
		do ii = 1, m
			do jj = 1, n

			!start, finish: C[start:finish, ...] 
			srow = (jj - 1) * q + 1
			frow =  jj * q
			scol = (ii - 1) * p + 1
			fcol =  ii * p 

			C(srow:frow, scol:fcol) = A(ii,jj) * B

			end do
		end do


	end function




	!don't set 0 and 1 separately https://rosettacode.org/wiki/Identity_matrix#Notorious_trick
	!-------------------------------------------------------------------------------------------------------------------------
    ! Function that returns an identity matrix (complex*16) given its dimension N
    ! Modular arithmetic: i/j = 0 if i < j. Then, i/j * j/i is 1 iff i == j. I don't know if it's actually faster
    !-------------------------------------------------------------------------------------------------------------------------
	
	function cEYE(N) result(C)

		integer*4, intent(IN) :: N
		integer*4			  :: ii, jj

		!imaginary and complex parts
		real*8, dimension(:,:), allocatable :: A, B

		complex*16, dimension(N,N) :: C

		!made to deallocate them later. Could I have used dimension(N,N)?
		allocate(A(N,N))
		allocate(B(N,N))

		B = 0

		ForAll( ii = 1:N, jj = 1:N) A(ii, jj) = (ii/jj) * (jj/ii)

		C = CMPLX(A, B)

		deallocate(A,B)

	end function





	!---------------------------------------------------------------------------------------------------------------------------
    ! Wrapper for ZHEEV, subroutine that takes a hermitian matrix and eigenvalues array and fills the latter with eigs. 
    ! Warning: eigenvalues vector needs to be allocated, either before being passed or inside the wrapper (in this case inside)
    !---------------------------------------------------------------------------------------------------------------------------
	subroutine C_DIAGONALIZE(ham, eigs)

		complex*16, dimension(:,:), allocatable :: ham
		real*8, dimension(:), allocatable :: eigs
		integer*4 :: N

		!declare zheev auxiliary variables, in caps lock
	    integer                                     :: INFO, LWORK
	    real*8, dimension(:), allocatable :: RWORK
	    complex*16, dimension(:), allocatable     :: WORK

	    N = size(ham,1)
	    LWORK = max(1, 2*N-1)

	    !important: eigs needs to be allocated before zheev
	    allocate(eigs(N))
	    allocate(RWORK(3*N-2))
	    allocate(WORK(LWORK))

	    !eig contains the eigenvalues
	    call ZHEEV('N', 'U', N, ham, N, eigs, WORK, LWORK, RWORK, INFO)


	    deallocate(RWORK, WORK)


	end subroutine



	!EIGENVECTORS

	!---------------------------------------------------------------------------------------------------------------------------
    ! Wrapper for ZHEEV, subroutine that takes K, a hermitian matrix, eigenvalues array and fills the latter with eigs. 
    ! Introduced eig_vec to avoid modifying the hamiltonian
    ! They can be used as columns of a projector P to truncate operators as P^T Operator P after removing some of them
    ! Warning: eigenvalues vector needs to be allocated, either before being passed or inside the wrapper (in this case inside)
    !---------------------------------------------------------------------------------------------------------------------------


	subroutine C_EIGENVECTORS(ham, eig_vec, eigs)

		complex*16, dimension(:,:), allocatable :: ham, eig_vec
		real*8, dimension(:), allocatable :: eigs
		integer*4 :: N, k

		!declare zheev auxiliary variables, in caps lock
	    integer                                     :: INFO, LWORK
	    real*8, dimension(:), allocatable :: RWORK
	    complex*16, dimension(:), allocatable     :: WORK

	    N = size(ham,1)
	    LWORK = max(1, 2*N-1)

	    !important: eigs needs to be allocated before zheev
	    allocate(eigs(N))
	    allocate(RWORK(3*N-2))
	    allocate(WORK(LWORK))

	    eig_vec = ham

	    !eig contains the eigenvalues, eig_vec the eigenvectors
	    call ZHEEV('V', 'U', N, eig_vec, N, eigs, WORK, LWORK, RWORK, INFO)


	    deallocate(RWORK, WORK)


	end subroutine



	!-------------------------------------------------------------------------------------------------------------------------
    ! Subroutine that projects the matrix M using P, as out_matrix = P^T M P. Using dynamic arrays, allocation is automatic
    ! Check: M is square matrix, MP is valid (#M columns == #P rows)
    !-------------------------------------------------------------------------------------------------------------------------

    function PROJECT(M,P) result(M_new)

    	complex*16, dimension(:,:) :: M,P
    	complex*16, dimension(:,:), allocatable :: M_new


    	call CHECK( debug = .FALSE., &
    				condition_to_check = ( size(M,2) == size(P,1) ), &
    				stoppage = .TRUE., &
    				message = 'P rows are', &
    				error_message = 'Error: matrix #columns differs from projector #rows ', &
    				variable = size(P,1) )



    	call CHECK( debug = .FALSE., &
    				condition_to_check = ( size(M,1) == size(M,2) ), &
    				stoppage = .TRUE., &
    				message = 'Input square matrix size is', &
    				error_message = 'Error: matrix to project is not a square matrix ', &
    				variable = size(M,1) )
    	

    	allocate(M_new(size(P,2), size(P,2)))

    	!the projection
    	M_new = MATMUL( MATMUL(  transpose(P), M)  ,P)


    end function



end module




!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************



module Quantum_Ising

	use Tensor_Products
	use ERROR_HANDLING

	implicit none


contains
	
	
	
	!-----------------------------------------------------------------------------------------------------------------------------------------------
    ! Given system size N and uniform magnetic strength h  (easily generalizable), returns sum_i (I ø  I ø ... \sigma_z^i ø ... )
    ! Uses fortran's dynamical array magic that allows H_single = TPROD(H_single, sigma_z), where H_single appears with two different dimensions
    !-----------------------------------------------------------------------------------------------------------------------------------------------

	!h is magnetic strength, N system dimension,
	function MAGNETIC_TERM (N, h) result(H_magnetic)

		integer*4, intent(IN) :: N
		integer*4             :: ii, jj
		real*8   , intent(IN) :: h
		complex*16, dimension(2,2) :: sigma_z, id
		complex*16, dimension(:, :), allocatable :: H_single, H_magnetic



		!So I won't initialize it everytime inside a cycle
		id = cEYE(2)

		sigma_z = reshape([1, 0, 0, -1], shape(sigma_z), order = [2,1] )

		allocate(H_magnetic(2**N, 2**N))
		H_magnetic = 0


		do ii = 1, N

			!starting term
			if (ii == 1) then
				H_single = sigma_z
			else
				H_single = id
			end if


			do jj = 1, N - 1

				IF ( jj == ii - 1) then

					H_single = TPROD(H_single, sigma_z)

				else

					H_single = TPROD(H_single, id)

				END IF

			end do

			!can generalize to h(i)
			H_magnetic = H_magnetic + h * H_single

		end do



	end function
	



	!---------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Given system size N and pairwise interaction strength J  (easily generalizable to Jij), returns sum_i (I ø  I ø ... \sigma_x^{i-1} ø \sigma_x^{i}... )
    ! H_tmp is a single I ... sigma x sigma x ... I to be added multiplied by J to H_pair, which is the returned matrix
    !---------------------------------------------------------------------------------------------------------------------------------------------------------
	function PAIR_SPIN_TERM(N, J) result(H_pair)

		integer*4, intent(IN) :: N
		real*8, intent(IN)    :: J
		integer*4             :: ii, jj


		complex*16, dimension(:,:), allocatable :: H_pair, H_tmp
		complex*16, dimension(2,2) :: sigma_x, id


		id = cEYE(2)
		sigma_x = reshape([0, 1, 1, 0], shape(sigma_x), order = [2,1] )

		allocate(H_pair(2**N, 2**N))
		H_pair = 0

		!position of sigma_x^(i-1) sigma_X^I
		do ii = 1, N-1

		!sigma_x sigma_x I I I ...
			IF (ii == 1) then

				H_tmp = TPROD(sigma_x, sigma_x)

				do jj = 3, N

					H_tmp = TPROD(H_tmp, id)

				end do
				H_pair = H_pair + J * H_tmp

			ELSE 

			! I ... and then sigma_x-s are at ii
				H_tmp = id

				do jj = 2, N - 1

				!place sigma_x sigma_x	
					if (jj == ii) then

						H_tmp = TPROD(TPROD(H_tmp, sigma_x), sigma_x)

					else
				!add another I
						H_tmp = TPROD(H_tmp, id)
					end if
				
				end do
				H_pair = H_pair + J * H_tmp

			END IF


		end do


	end function


	!*****************************************************************************
	!----------- DOUBLING THE SYSTEM SIZE ----------------------------------------
	!*****************************************************************************


	!---------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Given system size N and single site Hilbert space physical dimension d (a bit silly: I use sigma_x 2x2), returns A = Left System Interaction Term = (I ø  I ø ...  ø \sigma_x^{N})
    !                              B = Left System Interaction Term = (\sigma_x^{1} ø I ø  I ø ...  ø I)
    ! The interaction term H_AB in the 2N system Hamiltonian is A ø B
    ! Idea: TPROD(I,I) is still the identity, with doubled dimension
    !---------------------------------------------------------------------------------------------------------------------------------------------------------
    subroutine INIT_INTERACTION(N, d, A, B)

    	integer*4, intent(IN) :: N, d
    	complex*16, dimension(:,:), allocatable :: A, B
    	complex*16, dimension(2,2) :: sigma_x

    	if (allocated(A)) then
    		deallocate(A)
    	end if

    	if (allocated(B)) then
    		deallocate(B)
    	end if

    	sigma_x = reshape([0, 1, 1, 0], shape(sigma_x), order = [2,1] )


    	A = TPROD( cEYE(d**(N-1)), sigma_x)
    	B = TPROD( sigma_x, cEYE(d**(N-1)))

    end subroutine


    !---------------------------------------------------------------------------------------------------------------------------------------------------------
    ! Given H of N sites and H_AB interaction builds H_2N = TPROD(H_A, I^N) + TPROD(I^N, H_B) + TPROD(A,B)
    ! It's always the initial step, since you should pass H_tilde, A_tilde and B_tilde after P^T M P 
    ! Instead of using d**N with N = #spins for cEYE I used the matrix H size directly, to check that it doubles.
    !---------------------------------------------------------------------------------------------------------------------------------------------------------

    function Ham_Double(H, A, B) result(H_2N)

    	complex*16, dimension(:,:), intent(IN)  :: H, A, B
    	complex*16, dimension(:,:), allocatable :: H_2N
    	integer*4 :: N, t

    	!total number 2**N
    	t = size(H,1)


    	H_2N = TPROD(H, cEYE(t)) + TPROD(cEYE(t), H) - TPROD(A,B)
    	!- sign to have -sigma_x x sigma_x interaction

    	call CHECK( debug = .FALSE., &
    				condition_to_check = ( size(H_2N,1) == t**2 ), &
    				stoppage = .TRUE., &
    				message = 'Doubled Hamiltonian size is', &
    				error_message = 'Error in Ham_Double: Hamiltonian of 2N spins does not have the size of H_N squared ; instead it has', &
    				variable = size(H_2N,1) )    	


    end function


end module



!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************
!**********************************************************************************************************************************************************************************




program Real_Space_RG


use Quantum_Ising
use Tensor_Products
use Error_Handling

	implicit none


	complex*16, dimension(:,:), allocatable :: H_single, H_double, eig_vec, A, B, A_double, B_double, P
	!H_single = Hamiltonian of the system I have to double, which will have H_double
	!eig_vec to contain eigenvectors, P will contain the first 2**N of them and be used for projections
	!A and B form the interaction term in H_double via tensor product, both then have to dìbe doubled and projected

	!eigenvalues might risk overflow: E_gs = -JN (lambda = 0) and N is increasing. In that case we could divide H by something at each step 

	integer*4 :: N, ii, iteration_count
	integer*8 :: system_size
	real*8 :: h, J, e_gs, old_egs, threshold, step

	!N is the number of Ising spins in the initial system, ii an iterator. 
	!iteration_count is used to count how many iterations are required for convergence
	!system size starts from N and doubles, in order to get e = E_gs / (2^k*N) where 2^k*N are the total spins now
	!h = magnetic field, also denoted with lambda
	!J = pair spin interaction
	!e_gs = ground state energy
	!old_egs = used for the gs of the previous iteration, in order to check convergence
	!threshold = convergence threshold
	!step = size of step for going from 0 to 3 magnetic field and execute tha algorithm

	real*8, dimension(:), allocatable :: eigs
	!eigenvalues of Hamiltonian, initialized by diagonalization subroutine C_eigenvalues

	character(50) :: fmtString 
	fmtString = '(9999("{",F16.5,",",1X,F16.5,"}",:,1X))'
	!useful to print complex matrices



	N = 1
	J = 1.
	h = 1.

	threshold = 1e-13
	!how many digits of the result should be the same as the previous one


	!lambda = h from 0 to 3
	step = 0.005

	DO ii = 1, int(3./step) + 1

		h = step * ii

		!---------------------------------------------------------------------------
		!starting Hamiltonian for N spins
		!--------------------------------------------------------------------------

		H_single = magnetic_term(N, h) - pair_spin_term(N, J)
		call INIT_INTERACTION(N, 2, A, B)

		!--------------------------------------------------------------------------------------------
		!d = 2 physical space dim is second argument. Now assign values to make the while loop start
		!-----------------------------------------------------------------------------------------

		old_egs =  1.
		e_gs    = -1.
		system_size = N


		!-------------------------------------------------------------------------
		!RSRG algorithm
		!--------------------------------------------------------------------------
		iteration_count = 0
		do while( ABS(1 - old_egs/e_gs) > threshold  )

			old_egs = e_gs

			H_double = HAM_DOUBLE(H_single, A, B)
			A_double = TPROD( A, cEYE(2**N))
			B_double = TPROD(cEYE(2**N), B)
			
			call C_EIGENVECTORS(H_double, eig_vec, eigs )	
			P = dcmplx(eig_vec(:, :2**N))

			!---------------------------------------------------------------
			!the 2N system becomes the one to double in the next iteration
			!---------------------------------------------------------------

			H_single = PROJECT(H_double, P)
			A = PROJECT(A_double, P)
			B = PROJECT(B_double, P)
			system_size = 2 * system_size
			

			e_gs = eigs(1)/ system_size

			deallocate(eigs, eig_vec)
			iteration_count = iteration_count + 1

		end do

	!now the program converged to e_gs in it_count iterations for the current lambda; we write that down
	!write to file magnetic field, ground state energy density
	
	open(20, file = "e_gs_"//trim(str(int(N)))//".txt", position = 'append', action = 'write')
		write (20,*) h, e_gs, iteration_count
	close(20)

	END DO 

end program
