SUBROUTINE CHECK (debug, condition_to_check, verbose, stoppage, message, error_message, variable)

implicit none

	logical, intent(IN) :: debug, condition_to_check
	logical, intent(IN), optional :: verbose, stoppage
	class(*), optional :: variable
	character(*), optional :: message, error_message

	! used within subroutine for clarity if stoppage is not specified
	logical :: ALT

	if (.NOT.PRESENT(stoppage)) then
		ALT = .TRUE.
	else
		ALT = stoppage
	end if



	if (condition_to_check) then
		if (debug) then

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

				END SELECT
			END IF
			! ended variable type block
		END IF
		! ended debug block

		! if ALT is enabled (default or specified with stoppage), stop the program
		if (ALT) then 
			STOP
		end if

	! if the condition is false
	else
		if (present(error_message)) then
			print *, error_message
		end if
		! if ALT is enabled (default or specified with stoppage), stop the program
		if (ALT) then 
			STOP
		end if

	END IF
		

END SUBROUTINE CHECK