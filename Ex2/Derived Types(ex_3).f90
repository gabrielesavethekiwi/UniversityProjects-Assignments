MODULE complex_matrices

	implicit none

	type CMATRIX
		!matrix dimensions
		integer(8), dimension(2) :: dim

		!matrix elements
		complex(16), dimension(:,:), allocatable :: elements

		!store determinant and trace
		complex(16) :: Det, Tr
	end type 

	|-----------------------------------------------------------------------------------------------
	!interface: you can call these operators as an alias for the corresponding function/subroutine
	|-----------------------------------------------------------------------------------------------

	interface operator(.TR.)
		module procedure CMATRIX_TRACE
	end interface

	interface operator(.ADJ.)
		module procedure ADJOINT
	end interface


CONTAINS

	!-----------------------------------------------------------------
	!initialize cmatrix given rows N and columns M
	!-----------------------------------------------------------------
	function INIT_CMATRIX(N,M) result(mat)

		integer, intent(IN) :: N, M

		!two real matrices to fill with random numbers and then use cmplx
		real(16), dimension(N,M) :: A, B
		

		!define the complex matrix, access the dim vector and give it input rows, columns 
		type(cmatrix) :: mat 

		mat%dim(1) = N
		mat%dim(2) = M

		!error handling
		if ( (N <= 0) .OR. (M <= 0) ) then
			print*, 'Error: matrix initialization requires positive integer dimensions'
			STOP
		end if


		allocate(mat%elements(N, M))

		CALL RANDOM_NUMBER(A)
		CALL RANDOM_NUMBER(B)
		mat%elements = cmplx(A,B)

		!CALL RANDOM_NUMBER(mat%elements)

		!only square matrices have determinant and trace
		if ( size(mat%elements, 1) == size(mat%elements, 2) )  then

			!initialize complex determinant
			mat%det = (0d0, 0d0) 

			!complex trace via function acting on this type
			mat%tr = CMATRIX_TRACE(mat)
		end if


	END function init_cmatrix


	!----------------------------------------------------------------------
	!adjoint of a type CMATRIX
	!----------------------------------------------------------------------
	FUNCTION ADJOINT(in_mat) result (mat_dagger)

	type(CMATRIX), intent(IN) :: in_mat
	type(CMATRIX) :: mat_dagger

	mat_dagger%dim(1) = in_mat%dim(2)
	mat_dagger%dim(2) = in_mat%dim(1)

	!fortran has built in transposed function for 2-d arrays and conjugation
	!given an adjoint matrix, det and tr are obtained by complex conjugation

	mat_dagger%elements = conjg( transpose( in_mat%elements))
	mat_dagger%tr = conjg( in_mat%tr)
	mat_dagger%det = conjg( in_mat%det)

	END FUNCTION


	!-------------------------------------------
	! Trace of a type CMATRIX
	!-------------------------------------------

	FUNCTION CMATRIX_TRACE(in_mat) result(trace)

	complex(16) :: trace
	integer(8)  :: ii
	type(CMATRIX), intent(IN) :: in_mat


	trace = (0d0, 0d0)

	if ( in_mat%dim(1) == in_mat%dim(2) ) then

		do ii = 1, size(in_mat%elements, 1)
			trace = trace + in_mat%elements(ii, ii)
		end do
	else 
		print *, 'You asked for the trace of a non squared matrix'

	end if

	END FUNCTION cmatrix_trace


	!------------------------------------------------------------------------------
	! Subroutine that prints the matrix in a file given the CMATRIX and filename
	!------------------------------------------------------------------------------

	SUBROUTINE SAVE_CMATRIX(in_mat, filename)

		implicit none

		integer(8)    :: ii, jj
		character(*)  :: filename
		type(CMATRIX) :: in_mat

		!using unit < 20 seems to be discouraged online
		open(unit = 29, file = filename)

		DO ii = 1, size(in_mat%elements, 1)
			DO jj = 1, size(in_mat%elements, 2)

				write(29, '(1X,F20.10,F20.10)', advance = 'no') in_mat%elements(ii,jj)

			END DO
			! newline
			write(29, *) ''
		END DO

		!only write trace and determinant for square matrices
		if (size(in_mat%elements, 1) == size(in_mat%elements, 2)) then

			write(29, *) 'Trace'
			write(29, '(1X,F20.10,F20.10)', advance = 'no') in_mat%tr
			write(29, *) ''
			write(29, *) 'Determinant'
			write(29, '(1X,F20.10,F20.10)', advance = 'no') in_mat%det

		end if


		close(29)

	END SUBROUTINE save_cmatrix


END module complex_matrices









PROGRAM derived_types

	use complex_matrices
	implicit none


	!let's test on a matrix A and its transposed
	type(CMATRIX) :: A, A_dagger
	integer(4)    :: N, M                                      !dim is a vector with two components, int*8, dim(2) means you have two 4 bytes numbers
	complex(16)   :: trace, trace_dagger
	!character(len = :),allocatable :: filename_matrix
	!character(len = :),allocatable :: filename_transposed

	character(10) :: filename_matrix = 'matrix.txt'
	character(21) :: filename_transposed = 'transposed_matrix.txt'


	!ask for user input
	print*, 'Input the number of matrix rows, columns'
	read*, N, M


	!initialize complex matrix
	A = INIT_CMATRIX(N, M)
	trace = .TR.A

	A_dagger = .ADJ.A
	trace_dagger = .TR.A_dagger

	!print the traces, then check they are the same as in .txt
	print*, 'Trace:', trace
	print*, 'Adjoint trace', trace_dagger


	CALL save_cmatrix(A, filename_matrix)
	CALL save_cmatrix(A_dagger, filename_transposed)


END PROGRAM derived_types