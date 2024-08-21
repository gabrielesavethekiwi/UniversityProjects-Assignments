MODULE NAIVE_MULT
	implicit none


contains

	FUNCTION MULT_USUAL(mat_1, mat_2) result(mat_3)
	! usual row by column matrix multiplication. The first matrix has elements indexed by (i,j) and the second one (j,k).

	implicit none

	INTEGER*4 :: ii, jj, kk 

	REAL*8, dimension(:,:) :: mat_1, mat_2
	REAL*8, dimension(size(mat_1,1), size(mat_2,2)) :: mat_3

	do ii = 1, size(mat_1,1)
		do kk = 1, size(mat_2,2)
			do jj = 1, size(mat_1,2)

				mat_3(ii,kk) = mat_1(ii,jj)*mat_2(jj,kk) + mat_3(ii, jj)
			end do
		end do
	end do
	END FUNCTION MULT_USUAL



	FUNCTION MULT_CACHE(mat_1, mat_2) result(mat_3)
	! change the order of looping. Probably faster because before mat_2(jj,kk) in the innermost loop jumps rows and matrix elements are stored in the cache in row major order, so it has to be loaded slightly slower. Here mat_2(jj,kk) varies first ok k, column, probably won't miss cache. The fastest way would be to loop in the order (jj, kk, ii).

	implicit none

	INTEGER*4 :: ii, jj, kk

	REAL*8, dimension(:,:) :: mat_1, mat_2
	REAL*8, dimension(size(mat_1, 1), size(mat_2, 2)) :: mat_3

	do ii = 1, size(mat_1,1)
		do jj = 1, size(mat_1,2)
			do kk = 1, size(mat_2,2)

				mat_3(ii,kk) = mat_1(ii,jj)*mat_2(jj,kk) + mat_3(ii, jj)
			end do
		end do
	end do

	END FUNCTION mult_cache


END MODULE NAIVE_MULT


program Test_performance

use NAIVE_MULT
implicit none

INTEGER*2 N !matrix dimension, assuming only square matrices
REAL*8 :: START, STOP !to store CPU times.

REAL*8, dimension(:,:), allocatable :: A, B, out_usual, out_cache, out_fortran 
!A,B are mat_1 and mat_2. Out-s are AB changing method 






DO N = 10, 1000, 20



	allocate(A(N,N))
	allocate(B(N,N))

	CALL RANDOM_NUMBER(A)  !uniform between 0,1 for each entry
	CALL RANDOM_NUMBER(B)


	allocate(out_usual(N,N))
	allocate(out_cache(N,N))
	allocate(out_fortran(N,N))


	!initialize entries with zero because the naive function use mat_3 = mat_3 +... in innermost loop

	out_usual = 0.
	out_cache = 0.
	out_fortran = 0. 

	 
	!data for usual matrix moltiplication
	CALL CPU_TIME(start)
	out_usual = MULT_USUAL(A,B)
	CALL CPU_TIME(stop)

	!in fortran unit = N, file = ... N is used to refer to that same file

	open(unit = 0, file = 'usual_matrix_multiplication.txt', position = 'append', action = 'write')
	write(0,*) N, stop-start
	close(0)



	!data for matrix multiplication changing innermost loop, now theoretically a bit faster
	CALL CPU_TIME(start)
	out_cache = MULT_CACHE(A,B)
	CALL CPU_TIME(stop)

	OPEN(unit = 1, file = 'cache-friendly_matrix_multiplication.txt', position = 'append', action = 'write')
	WRITE(1,*) N, stop-start
	CLOSE(1)
	!using 1 instead of repeating 0 was not necessary



	!fortran built in function
	CALL CPU_TIME(start)
	out_fortran = MATMUL(A,B) !fortran built in
	CALL CPU_TIME(stop)

	open(unit = 2, file = 'fortran_matrix_multiplication.txt', position = 'append', action = 'write')
	write(2,*) N, stop-start
	close(2)


	deallocate(A)
	deallocate(B)
	deallocate(out_usual)
	deallocate(out_fortran)
	deallocate(out_cache)


END DO

END PROGRAM Test_performance