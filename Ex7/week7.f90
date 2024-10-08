module TENSOR_PRODUCTS

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



end module






module Quantum_Ising

	use Tensor_Products
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



	!https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
    !third answer; this works because "i0" format is assumed == print with least number of digits require -> much simpler than my previous implementation

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str


end module






program Q_Ising

	use TENSOR_PRODUCTS
	use Quantum_Ising

	implicit none


	integer*4 :: N, ii, k
	real*8 :: h, J
	complex*16, dimension(:,:), allocatable :: Ham_Ising

	!watch out: now the subroutine for diagonalization will initialize it
	real*8, dimension(:), allocatable :: eigs

	

	!https://stackoverflow.com/questions/40547599/format-specifier-for-a-complex-matrix
	character(50) :: fmtString  
	fmtString = '(9999("{",F6.3,",",1X,F6.3,"}",:,1X))'


    !ask user for input
	print *, '********************************************************************************************'
	print *, 'This program returns the first k eigenvalues of a 1-D Ising spin chain, composed of N spins.'
	print *, 'The spins are under a magnetic field h and have pairwise interaction strength J.'
    print *, 'Insert number of spins N'
    read(*, *) N
	print *, 'Insert field strength h'
    read(*, *) h
    print *, 'Insert pairwise interaction strength J' 
    read(*, *) J
    print *, 'Insert eigenvalue cutoff k'
    read(*,*) k
    print *, "*********************************************************************************************"


	Ham_Ising = magnetic_term(N, h) - pair_spin_term(N, J)

	call C_DIAGONALIZE(Ham_Ising, eigs)

   
	open(20, file = "eigs_"//trim(str(int(N)))//".txt", position = 'append', action = 'write')

		do ii = 1, k
			write(20, *) h, J, eigs(ii)
		end do

	close(20)

	!not elegant; gets allocated by c_diagonalize
	deallocate(eigs)

end program




