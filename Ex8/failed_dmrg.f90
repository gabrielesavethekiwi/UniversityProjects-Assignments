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
    	logical :: debug_flag

    	debug_flag = .FALSE.


    	call CHECK( debug = debug_flag, &
    				condition_to_check = ( size(M,2) == size(P,1) ), &
    				stoppage = .TRUE., &
    				message = 'P rows are', &
    				error_message = 'Error: matrix #columns differs from projector #rows ', &
    				variable = size(P,1) )



    	call CHECK( debug = debug_flag, &
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


		call CHECK ( debug = .TRUE., &
    			 condition_to_check = ( ABS(trace - 1) < 0.0000001  ), &
    			 stoppage = .TRUE., &
    			 message = 'This is the trace', &
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




!added later for convenience when you know what kind of state you are dealing with

	FUNCTION wavefunction_to_dmatrix(sys_psi) result(density_matrix)

		complex*16, dimension(:), allocatable :: sys_psi
		real*8 :: trace
		complex*16, dimension(:,:), allocatable :: density_matrix, transposed, diff

		!total number of elements in vector
		integer*4 :: NN, ii


		NN = size(sys_psi)
		allocate (density_matrix(NN,NN))


		!https://emeryyi.wordpress.com/2010/05/22/outer-product-of-two-vectors-in-fortran/
		!copy a vector ncopies times making it the rows (dim = 1) or the columns (dim = 2) of a matrix
		! A * B multiplies elementwise the matrices, NO matmul


		density_matrix = spread(sys_psi(1:NN), dim = 2, ncopies = NN)  *  spread( CONJG(sys_psi(1:NN)), dim = 1, ncopies = NN) 



		!check Tr rho = 1. Think about Tr rho^2 = 1 for pure states
		trace = 0

		do ii = 1, NN

			trace = trace + REAL(density_matrix(ii, ii))

		end do


		call CHECK ( debug = .FALSE., &
    			 condition_to_check = ( ABS(trace - 1) < 0.0000001  ), &
    			 stoppage = .TRUE., &
    			 message = 'This is the trace', &
    			 error_message = 'Error: Tr(rho) \= 1', &
    			 variable = trace)


	end function






END MODULE quantum_states



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


    	H_2N = TPROD(H, cEYE(t)) + TPROD(cEYE(t), H) + TPROD(A,B)

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
! m = 3 spins  1 spin   b= block sb = super block
!   -------              -------   
!   |o o o|     o |  o   |o o o| 
!   -------              -------
!   H_b      sigma_z      H_b
!------------------------------------------------------------
!
! 1. |psi>_GS of H_sb = H_left + H_right + H_int
!    H_left = H_b ø I +  I ø ... ø I ø sigma_z + I^1 ø ...I^{m-1} ø sigma_x ø sigma_x
!    H_right = I ø H_b +  sigma_z ø I ø ... ø I + sigma_x ø sigma_x ø I^1 ø ...I^{m-1}
!    H_int = I_(2^mx2^m) ø sigma_x ø sigma_x ø I_(2^mx2^m) make the two   o | o   interact  
!
! 2. rho = |psi>_GS <psi|_GS, trace out a subsystem. 
! 3. Trace out WLOG right subsystem, diagonalize reduced density matrix rho_l and keep the first m eigenvectors in a matrix P
! 4. H_b_tilde = P^T H_left P, so I aumented by 1 spin but kept a 2^mx2^m matrix. Calculate physical quantities, then restart with this new H_b
!------------------------------------------------------------------------------------

program Real_Space_RG


use Quantum_Ising
use Tensor_Products
use Error_Handling
use Quantum_States

	implicit none


	complex*16, dimension(:,:), allocatable :: H_b, H_left, H_right, H_SB , eig_vec, rho, rho_L, P, A, B, H_b_left, H_b_right
	complex*16, dimension(:,:), allocatable :: H_edge_L, H_edge_R

	!A, B represent rightmost and leftmost interaction if you tensor product them. So they remain at the edge
	
	!eig_vec to contain eigenvectors, P will contain the first 2**N of them and be used for projections

	!eigenvalues might risk overflow: E_gs = -JN (lambda = 0) and N is increasing. In that case we could divide H by something at each step 

	complex*16, dimension(:), allocatable :: psi_gs
	complex*16, dimension(2,2) :: sigma_x, sigma_z


	integer*4 :: m, ii, kk, iteration_count, k
	integer*8 :: system_size
	real*8 :: h, J, e_gs, old_egs, threshold, step
 
	!iteration_count is used to count how many iterations are required for convergence
	!system size starts from 
	!h = magnetic field, also denoted with lambda
	!J = pair spin interaction
	!e_gs = ground state energy
	!old_egs = used for the gs of the previous iteration, in order to check convergence
	!threshold = convergence threshold
	!step = size of step for going from 0 to 3 magnetic field and execute tha algorithm

	real*8, dimension(:), allocatable :: eigs
	!eigenvalues holder, initialized by diagonalization subroutine C_eigenvalues

	character(50) :: fmtString 
	fmtString = '(9999("{",F16.5,",",1X,F16.5,"}",:,1X))'
	!useful to print complex matrices



	m = 1
	J = 1.
	h = 1.
	
	k = 4

	threshold = 1e-6
	!how many digits of the result should be the same as the previous one

	sigma_x = reshape([0, 1, 1, 0], shape(sigma_x), order = [2,1] )
	sigma_z = reshape([1, 0, 0, -1],shape(sigma_z), order = [2,1] )


	!lambda = h from 0.005 to 3; this avoids degeneracy 
	step = 0.005

	DO ii = 1, int(3./step) + 1

		h = step * ii


		!---------------------------------------------------------------------------
		!starting block Hamiltonian for m spins
		!--------------------------------------------------------------------------

		H_b = magnetic_term(m, h) - pair_spin_term(m, J)

		sigma_z = reshape([1, 0, 0, -1],shape(sigma_z), order = [2,1] )
		sigma_z = h *reshape([1, 0, 0, -1],shape(sigma_z), order = [2,1] )
		!don't forget the h and J! I won't worry about J since J = 1, otherwise multiply sigmax by sqrt(J) maybe

		!--------------------------------------------------------------------------------------------
		!Now assign values to make the while loop start
		!-----------------------------------------------------------------------------------------

		old_egs =  1.
		e_gs    = -1.
		system_size = 2*m + 2

!****************************************************************************************
!-------------------- FIRST STEP -------------------------------------------------------
!****************************************************************************************
! A single block has a 2**m x 2**m Hamiltonian

		A = TPROD(cEYE(2**m), sigma_x)
		B = TPROD(sigma_x, cEYE(2**m))
		!to build interaction term in the superblock hamiltonian 

		H_b_left = H_b
		H_b_right = H_b

		H_edge_L = TPROD(cEYE(2**(m-1)),  sigma_x)
		H_edge_R = TPROD(sigma_x, cEYE(2**(m-1)))
		!The block can be connected to the new spin with this operator TPROD sigma_x of the rogue spin


		H_left = TPROD(H_b_left, cEYE(2)) + TPROD(cEYE(2**m ), sigma_z) - TPROD( H_edge_L, sigma_x)   
		H_right= TPROD(cEYE(2), H_b_right) + TPROD(sigma_z, cEYE(2**m )) -TPROD( sigma_x, H_edge_R)
		!You could use the same H_b due to symmetry
		!Imagine a block with Hamiltonian H_b_left interacting with a single spin. Their spaces have dimensions m' = m_truncation and 2, so this is 2m'x2m' 

		H_SB = TPROD(H_left, cEYE(2**(m+1))) + TPROD(cEYE(2**(m+1)), H_right) - TPROD(A,B)		


		call C_EIGENVECTORS(H_SB, eig_vec, eigs )	
		psi_gs = eig_vec(:, 1)
		!diagonalize the super-block Hamiltonian
			
		e_gs = eigs(1)/ system_size
		system_size = system_size + 2 
		!every step adds two total spins to the superblock

		rho = wavefunction_to_dmatrix(psi_gs)
		rho_L = PARTIAL_TRACE(2**(m+1), 2**(m+1), rho)
		!density matrix from GS (avoid degeneracy h = 0, or you will get a linear combination of degenerate states)

		deallocate(eigs, eig_vec)
		!make space for reduced density matrix eigs

		call C_EIGENVECTORS(rho_L, eig_vec, eigs )
		!P = eig_vec(:, size(eig_vec)-min(2**(m+1), k) +1: size(eig_vec))

 	P = eig_vec(:, size(eig_vec)-k +1: size(eig_vec))


		!keep last (bigger) k eigenvectors of P; if k is too big it will keep 2^m


		!-----------------------------------------------------------------------------------
		!m+1 block Hamiltonian projected in 2**m x 2**m subspace; project also interactions
		!-----------------------------------------------------------------------------------

		if (allocated(H_b_left)) then
			deallocate(H_b_left)
		end if

		if (allocated(H_b_right)) then
			deallocate(H_b_right)
		end if


		H_b_left = PROJECT(H_left, P)
		H_b_right = PROJECT(H_right, P)

!useful to check if matrix is dying due to lack of eigenvalues			
print*, 'This is H_b_left 1', size(H_b_left,1), size(H_b_left,2), 'with m', m
		do kk = 1, size(H_b_left,1)
 			print*, H_b_left(kk, :)
		end do 
	print*, '-----------------------------------------'


		!you need to create an operator representing the new m+1 block, that allows its
		!interaction with new spins or whatever. This is I x... x sigma_x on the rightmost (for left)
		!This will however be projected, as H_b is now a projection of a m+1 operator too, and used with H_b for the next left block Hamiltonian


		deallocate(H_edge_L, H_edge_R)

		H_edge_L = PROJECT(A, P)
		H_edge_R = PROJECT(B, P)

!useful to check if matrix is dying due to lack of eigenvalues			
print*, 'This is H_edge left 1'
		do kk = 1, size(H_edge_L,1)
 			print*, H_edge_L(kk, :)
		end do 

		deallocate(eigs, eig_vec, rho, H_left, H_right, A, B, H_SB)
		!Watch out when leaving allocated dynamic arrays which will be set equal to matrices with a lower dimension
		!H_left would mantain 2**m, 2**m with lots of 0 instead of becoming k x k !

		iteration_count = iteration_count + 1

!****************************************************************************************
!-------------------- SECOND STEP -------------------------------------------------------
!****************************************************************************************
! Now you have a block, with effective Hamiltonian H_left k x k, interacting with a new spin
!	     _________          _________
!		|____L____| o  | o |____R____|
!            1      2    3      4  
!
!    H_left = H_1 x I_2 + I_1 x H_2 + H_12
!    And I_1 means a kxk identity 
!    
!    H_12 is tricky: L block interacts with a rightmost sigma_x, but it had to be made k x k with P
		
		A = TPROD(cEYE(k), sigma_x)
		B = TPROD(sigma_x, cEYE(k))
		! Superblock Interaction happens between 2 and 3, so kxk identities on 1,4

		iteration_count = 0
		do while( ABS(1 - old_egs/e_gs) > threshold  )
			
			old_egs = e_gs


			H_left = TPROD(H_b_left, cEYE(2)) + TPROD(cEYE(k), sigma_z) -  TPROD( H_edge_L, sigma_x)  
			H_right= TPROD(cEYE(2), H_b_right) + TPROD(sigma_z, cEYE(k)) - TPROD( sigma_x, H_edge_R)


		do kk = 1, size(rho,1)
 			print*, H_left(kk, :), 'with size of H left that should be kxk', size(H_left,1), size(H_left,2)
		end do 




			H_SB = TPROD(H_left, cEYE(2*k)) + TPROD(cEYE(2*k), H_right) - TPROD(A,B)

			call C_EIGENVECTORS(H_SB, eig_vec, eigs )	
			psi_gs = eig_vec(:, 1)
			!diagonalize the super-block Hamiltonian
			
			e_gs = eigs(1)/ system_size
			system_size = system_size + 2 
			!every step adds two total spins to the superblock


			rho = wavefunction_to_dmatrix(psi_gs)
		do kk = 1, size(rho,1)
 			print*, rho_L(kk, :), 'with size of rho full density mat', size(rho,1), size(rho, 2)
		end do 

			rho_L = PARTIAL_TRACE(size(H_left,1), size(H_left,1), rho)

		do kk = 1, size(rho_L,1)
 			print*, rho_L(kk, :), 'with size', size(rho_L,1), size(rho_L, 2 )
		end do 


			deallocate(eigs, eig_vec)
			!make space for reduced density matrix eigs

			call C_EIGENVECTORS(rho_L, eig_vec, eigs )
			!P = eig_vec(:, size(eig_vec)-min(2**(m+1),k) +1: size(eig_vec))

		P = eig_vec(:, size(eig_vec)-k +1: size(eig_vec))

		do kk = 1, size(P,1)
 			print*, P(kk, :)
		end do 

			!-----------------------------------------------------------------------------------
			!m+1 block Hamiltonian projected in k x k subspace; project also interactions
			!-----------------------------------------------------------------------------------

			if (allocated(H_b_left)) then
				deallocate(H_b_left)
			end if

			if (allocated(H_b_right)) then
				deallocate(H_b_right)
			end if


			H_b_left = PROJECT(H_left, P)
			H_b_right = PROJECT(H_right, P)

!useful to check if matrix is dying due to lack of eigenvalues			
print*, 'This is H_b_left 2'
		do kk = 1, size(H_b_left,1)
 			print*, H_b_left(kk, :)
		end do 



			!create operator to make block interact with other blocks (for us a new rogue spin)
			!You need an identity on all the left block, 1, and a rightmost sigma_x on the new spin it acquired

			H_edge_L = PROJECT(A, P)
			H_edge_R = PROJECT(B, P)


			deallocate(eigs, eig_vec, rho)
			iteration_count = iteration_count + 1

			print*, e_gs

		end do

	!now the program converged to e_gs in it_count iterations for the current lambda; we write that down
	!write to file magnetic field, ground state energy density
	
	open(20, file = "infiniterg_"//trim(str(int(m)))//trim(str(int(k)))//".txt", position = 'append', action = 'write')
		write (20,*) h, e_gs, iteration_count
	close(20)

	END DO 

end program
