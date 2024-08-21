!---------------------------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------DOCUMENTATION-------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------------------------
! Module Name: error_handling
! CoNteNts: support module for the program Test_performances. Includes a very general subroutine for debugging and a function for 
!           checking if two matrices are the same
!---------------------------------------------------------------------------------------------------------------------------------------------------------
! subroutine name: Check
! DescriptioN: by passing the argument DEBUG = .TRUE. you caNncheck whether a condition, important for the correctedness of 
!              your program is true ( aNd output an optional message) or false. in this case, you can output an optional error message 
!              and stop the program if STOPPAGE = .TRUE. 
!              optionally, you can pass in a variable to check its type
! input: 
!        debug               =  bool, if false the check subroutine is turned off
!        condition_to check  =  bool, the condition to check in the debugging procedure, i.e. matrix_rows > 0
!        verbose             =  optional bool, Not yet implemeNted. in the future it might suppress messages or give shorter oNes
!        stoppage            =  optional bool, if .TRUE. the program gets stopped if condition_to_check is false. Default is .TRUE.
!        message             =  optional character, message in the case the program works fine (condition_to_check true)
!        error_message       =  optional character, error message (condition_to_check false)
!        variable            =  optional class(*), caN be whatever fortraN Native type. Pass it to check the type of something in your program
!---------------------------------------------------------------------------------------------------------------------------------------------------------



MODULE ERROR_handling

    implicit none


    contains

    subroutine CHECK (debug, condition_to_check, verbose, stoppage, message, error_message, variable)

    implicit none

        logical, intent(in) :: debug, condition_to_check
        logical, intent(in), optional :: verbose, stoppage
        class(*), optional :: variable
        character(*), optional :: message, error_message

        !ALT is dummy used within the subroutine for clarity, has the role of stoppage if specified
        
        logical :: ALT

        if (.NOT.present(stoppage)) then
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
                            print *, variable, 'real(8): real*8'
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
            end IF
            


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
                            print *, variable, 'real(8): real*8'
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
            end IF


            ! if ALT is enabled (default or specified with stoppage), stop the program
            if (ALT) then 
                STOP
            end if

        end IF
            
    end subroutine CHECK


end MODULE ERROR_handling












module schroedinger

	implicit none


    TYPE GRID_OBJECT

        integer*4  :: dim
        real*8     :: lower
        real*8     :: upper 
        real*8     :: dx

        real*8, dimension(:), allocatable :: lattice

    end TYPE


contains
	

     


	!-----------------------------------------------------------------------------------------------------------
    ! function that creates grid given the number of points N_tot of the lattice [minimum, maximum] (arguments)
    !-----------------------------------------------------------------------------------------------------------
	function GRID_1D(N_tot, minimum, maximum) result(grid)

		integer*4          :: N_tot, ii    						!*4 up to 2147483647
		real*8, intent(in) :: minimum, maximum				! real*8, 15 digits accuracy, order 10**[-308, 308]
        type(GRID_OBJECT)  :: grid

        grid%dim   = N_tot 
		grid%lower = minimum
        grid%upper = maximum
        grid%dx    = (maximum-minimum)/N_tot
        
        allocate(grid%lattice(grid%dim))

		do ii = 0, N_tot - 1
			grid%lattice(ii+1) = minimum + (ii * grid%dx)
		end do

	end function




    !-----------------------------------------------------------------------------------
    !Potential 1/2 (q-t/T)^2, q,t,TT are real numbers
    !--------------------------------------------------------------------------------------
    function Calculate_V (q, t, TT) result(V_t)

        real*8, intent(in) :: q, t, TT
        real*8             :: V_t

        V_t = ((q - t/TT)**2d0 )/2d0


    end function


    !https://www.mcs.anl.gov/~itf/dbpp/text/node84.html
    !-------------------------------------------------------------------------------------------
    !Calculates the expected value of an operator O, using array operations like numpy
    !Performs int (dl |psi(l)|^2 * O(l) ) where l = x,p for example.
    !psi = wavefunction, c_vector; O = real vector containing function (ex: x**2), dx = spacing
    !Warning: operator must be diagonal in l-representation
    !--------------------------------------------------------------------------------------------
    function MEAN_VALUE(psi, O , dx) result(mv)

        complex*16, dimension(:), allocatable   :: psi
        real*8                                  :: dx, mv
        real*8, dimension(:), allocatable       :: O
        integer*4                               :: ii

        mv = 0

        do ii = 1, size(psi)
            mv = mv + ABS(psi(ii))**2.d0 * O(ii) * dx
        end do

    end function 



    !------------------------------------------------------------------------------------------------------------
    !Subroutine that takes a complex vector, the lattice (ex x=[dx, 2dx, ...]) and spacing dx, also the total 
    !number of points N just to be sure, and puts in the complex vector the harmonic oscillator groundstate
    !------------------------------------------------------------------------------------------------------------
    subroutine INIT_PSI(psi, xlattice, dxlattice, N)

        integer*4 :: N, ii
        complex*16, dimension(:), allocatable   :: psi
        real*8, dimension(:), allocatable :: xlattice
        real*8 ::  dxlattice

        do ii = 1, N
            psi(ii) = EXP(-(xlattice(ii)**2d0)/2d0)
        end do

        psi = psi / Norm(psi, dxlattice)

    end subroutine


    !-------------------------------------------------------------------
    !Calculates norm of wavefunction given complex vector and spacing
    !-------------------------------------------------------------------
    function NORM(psi, dx) result(psi_Norm)

        complex*16, dimension(:), allocatable :: psi
        real*8                                :: dx, psi_Norm

        !use SUM on array |psi_i|^2 * dx, elementwise operatioNs
        psi_Norm = SQRT( SUM(  ABS(psi)**2d0 * dx  ) )

    end function




    !---------------------------------------------------------------------------------------------------
    !Following two functions calculate the fast fourier transform (x -> p) and the inverse one given
    !complex vector wavefunction. Not normalized
    !-----------------------------------------------------------------------------------------------------
        
        function fft(psi) result(psi_p)
            
            complex*16, dimension(:) :: psi
            complex*16, dimension(size(psi)) :: psi_p
            
            integer*8 :: plan


            call dfftw_plan_dft_1d(plan, size(psi), psi, psi_p, -1)
            call dfftw_execute_dft(plan, psi, psi_p)
            call dfftw_destroy_plan(plan)
        
        end function 


        
        function anti_fft(psi) result(psi_x)
            
            complex*16, dimension(:) :: psi

            complex*16, dimension(size(psi)) :: psi_x
            integer*8 :: plan


            call dfftw_plan_dft_1d(plan, size(psi), psi, psi_x, +1)
            call dfftw_execute_dft(plan, psi, psi_x)
            call dfftw_destroy_plan(plan)
        
        end function 










    !MODIFY
	!----------------------------------------------------------------------------------------------------------------------------
    ! subroutine (fuNctioNs caN't returN two results) that should create a symmetric tridiagoNal matrix. Since DSTEV Needs oNly
    ! aN array with the diagoNal terms aNd oNe for the subdiagoNal, d aNd d_sub arrays get modified to coNtain those.
    !--------------------------------------------------------------------------------------------------------------------------
	subroutine TRI_H( d, sub_d, omega, grid)

 		
 		real*8, dimension(:), intent(OUT)  :: d, sub_d !OUT used in subroutines to Not take them initialized
 		real*8 :: dx
 		real*8, dimension(:), intent(in) :: grid        !do Not modify it
 		integer*4 :: ii, omega

 		dx = grid(2) - grid(1)

 		do ii = 1, size(grid,1)
 			d(ii) = 1d0/(dx**2d0) +  (omega**2d0 * grid(ii)**2d0)/2d0 
 		end do

 		!subdiagonal has 1 less value than diagonal
 		do ii = 1, size(grid,1)-1
 			sub_d(ii) = -1d0/(2d0*dx**2d0) 
 		end do


	end subroutine


end module schroedinger





MODULE printing

	implicit none

contains
	

	!-------------------------------------------------------------------------------------------------------------
    ! subroutine that prints vector in a file giveN the real vector aNd filename. CLEARS it if file exists
    !-------------------------------------------------------------------------------------------------------------
    subroutine SAVE_VECTOR(in_vec, filename)

        implicit none

        integer*8 :: ii
        character(*) :: filename
        real*8, dimension(:) :: in_vec

        open(unit = 30, file = filename, action = 'write')

        do ii = 1, size(in_vec,1)
            write(30, '(1X,F20.10)') in_vec(ii)
        end do

        close(30)

    end subroutine save_vector


    !-----------------------------------------------------------------------------------
    ! subroutine that prints matrix in a file giveN the real matrix aNd filename
    !------------------------------------------------------------------------------------

    subroutine SAVE_MATRIX(in_mat, filename)

        implicit none

        integer(8)    :: ii, jj
        character(*)  :: filename
        real*8, dimension(:,:) :: in_mat

        !using unit < 20 seems to be discouraged online
        open(unit = 29, file = filename)

        DO ii = 1, size(in_mat, 1)
            DO jj = 1, size(in_mat, 2)

                write(29, '(1X,F20.10,F20.10)', advance = 'no') in_mat(ii,jj)

            end DO
            ! newline
            write(29, *) ''
        end DO

        close(29)

    end subroutine save_matrix


    !https://stackoverflow.com/questions/1262695/convert-integers-to-strings-to-create-output-filenames-at-run-time
    !third answer; this works because "i0" format is assumed == print with least number of digits require -> much simpler than my previous implementation

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str


end module printing





PROGRAM t_ev_harmonic

	use schroedinger
	use error_handling
	use printing


    implicit none

	


    complex*16, dimension(:), allocatable   :: psi_x, psi_momentum

    type(GRID_OBJECT)                       :: x_grid, p_grid, t_grid

    !kinetic and potential; P momenta vector, Q2 = (q1^2,...,qN^2)
    complex*16, dimension(:), allocatable   :: K, V 
    real*8, dimension(:), allocatable :: P, Q2, P2

    real*8                            :: L, pi, dp
    real*8                            :: mean_x, mean_x2, mean_p, mean_p2
    integer                           :: N_x, N_t, ii, jj, ll, kk

    !will contain the times to pass to the program
    real*8, dimension(6)              :: TT


    !most precise way, if you can't trust architecture put the what you need by hand
    pi    = 4.d0 * datan(1.d0)

    
    ! L should be > 20, you want 2L/N_x < 0.01 but I wouln't go below 0.001
    L   = 40.d0
    N_x = 12000

    N_t = 1000

    ! different values of TT
    TT  = (/ 1.d0, 2.d0, 5.d0, 8.d0, 25.d0, 50.d0 /)

    do ll = 1, size(TT)
        

        x_grid = Grid_1D(N_x,-L, L)
        !momentum grid [-pi/L, pi/L]
        p_grid = Grid_1D(N_x, -pi/x_grid%dx, pi/x_grid%dx)
        t_grid = Grid_1D(N_t, 0.d0, TT(ll))
        
        !prepare vector with all p, in momentum basis, for kinetic term
        allocate(P(N_x))

        dp = p_grid%dx

        !fftw3 wants [0,...,pi/L, -pi/L, ..., - dp ]: increasing order of positive eigenvalues, then negative ones
        do ii = 0, int(N_x/2) - 1
            P(ii + 1) = dp * dble(ii)
        end do

        do ii = int(N_x/2), N_x -1 
            P(ii + 1) = - dp * dble(N_x - ii) 
        end do

        
        !will be the Observable in Mean_Value function
        allocate( Q2(N_x), P2(N_x))
        Q2 = x_grid%lattice**2
        P2 = P**2




        allocate( psi_x(N_x), psi_momentum(N_x),K(N_x), V(N_x)  )
        

        !with this implementation you don't need to declare filename variables

        open(20, file = "T_"//trim(str(int(TT(ll))))//"_times.txt")
        open(21, file = "T_"//trim(str(int(TT(ll))))//"_x.txt")
        open(22, file = "T_"//trim(str(int(TT(ll))))//"_sigma_x.txt")
        open(23, file = "T_"//trim(str(int(TT(ll))))//"_p.txt")
        open(24, file = "T_"//trim(str(int(TT(ll))))//"_sigma_p.txt")
        open(25, file = "T_"//trim(str(int(TT(ll))))//"_wavex.txt")


        ! |psi(0)> = |0> of Harmonic Oscillator, gaussian
        call INIT_PSI(psi_x, x_grid%lattice, x_grid%dx, N_x)

        !perform step-by-step time simulation, with split operator method
        DO ii = 1, N_t
            do jj = 1, N_x

                ! exp(- i/2 V dt)
                V(jj) = EXP(   -0.5 * dcmplx(0d0,1d0) * t_grid%dx * CALCULATE_V(x_grid%lattice(jj), t_grid%lattice(ii), TT(ll))   )

                ! exp(- i T dt) with T = p^2/2
                K(jj) = EXP(-0.5 * dcmplx(0d0,1d0) * t_grid%dx * P(jj)**2d0 )
                
            end do

            !apply potential term exp(-i/2 V dt)
            do jj = 1, N_x
                psi_x(jj) = V(jj) * psi_x(jj)
            end do

            !go to momentum space, I understand psi_x might be misleading
            psi_x = fft(psi_x)/SQRT( dble(N_x) )
            !psi_x = psi_x / NORM(psi_x) done with /sqrtNx


            !apply kinetic term in p space
            do jj = 1, N_x
                psi_x(jj) = K(jj) * psi_x(jj)
            end do

            !go back to x representation
            psi_x = anti_fft(psi_x)/SQRT( dble(N_x) )

            do jj = 1, N_x
                psi_x(jj) = V(jj) * psi_x(jj)
            end do

            !normalize
            psi_x = psi_x/Norm(psi_x, x_grid%dx)


            !p representation of wavefunction for momenta related observables
            psi_momentum = fft(psi_x)

            psi_momentum = psi_momentum/NORM(psi_momentum, dp)

            

            mean_x  = Mean_Value(psi_x, x_grid%lattice, x_grid%dx)
            mean_x2 = Mean_Value(psi_x, q2, x_grid%dx)
            mean_p  = Mean_Value(psi_momentum, P, dp)
            mean_p2 = Mean_Value(psi_momentum, P2, dp)

            write(20, *) t_grid%lattice(ii)
            write(21, *) mean_x
            write(22, *) SQRT(mean_x2 - mean_x**2)
            write(23, *) mean_p
            write(24, *) SQRT(mean_p2 - mean_p**2)
            


            if( MOD(ii, 100) == 0) then

                do kk = 1, size(psi_x)
                    write(25, *) x_grid%lattice(kk), dble(psi_x(kk) * CONJG(psi_x(kk)))
                end do

           end if


        END do

        deallocate( psi_x, psi_momentum, P)
        deallocate( V, K, q2, p2)
        deallocate( x_grid%lattice, p_grid%lattice, t_grid%lattice )

        close(20)
        close(21)
        close(22)
        close(23)
        close(24)
        close(25)

    end do


end program