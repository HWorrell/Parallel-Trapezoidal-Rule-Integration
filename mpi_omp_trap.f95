program trapMpiOmp

        use mpi

        use omp_lib

!declare variables

        integer(kind=4) :: status(MPI_Status_Size)

	integer(kind=8) :: intervals, local_intervals, counts = 0
	
	integer(kind=4) :: worldsize, myrank = 0

	real(kind=8) :: a, b, start, finish, elapsed, interval_range = 0

	real(kind=8) :: local_a, local_range, local_sum, x, global_sum = 0
	
    character(len=640) :: arg

	integer :: mpi_error_code
	
!mpi init

	CALL MPI_Init(mpi_error_code)
	
	start = MPI_Wtime()

	call MPI_Comm_size ( MPI_COMM_WORLD, worldsize, mpi_error_code )

	call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, mpi_error_code )
	

!if process 0

	if(myrank == 0) then

!get number of intervals

	call getarg(1, arg)
	
	read(arg, *) a
	
	call getarg(2, arg)
	
	read(arg, *) b

	call getarg(3, arg)	
	
	read(arg, '(i100)') intervals

        print *, "A is ", a

        print *, "B is ", b

        print *, "N is ", intervals

	call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, status, mpi_error_code)
	
	call MPI_Bcast(b, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, status, mpi_error_code)
		
	call MPI_Bcast(intervals, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status, mpi_error_code)

	!calculate interval size

	interval_range = b - a
	
	local_intervals = intervals / worldsize
	
	local_range = (b - a) / intervals
	
	local_a = a + myrank * local_intervals * local_range
	
	call omp_set_num_threads(10)

        local_sum = (f(local_a) + f(local_a + local_range * local_intervals)) / 2.0

		counts = 0
		
	!$OMP PARALLEL DO PRIVATE(x) REDUCTION(+:local_sum)
	
	do counts = 0, (local_intervals - 1)
	
		x = local_a + counts * local_range
		
		local_sum = local_sum + f(x)
	
	end do
	
	!$OMP END PARALLEL DO

        local_sum = local_sum * local_range

        print *, myrank, local_sum

        call MPI_Reduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_error_code)
	
	print *, "The area is ", global_sum
	
!if not process 0
	else
	
	call MPI_Bcast(a, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, status, mpi_error_code)
	
	call MPI_Bcast(b, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, status, mpi_error_code)
		
	call MPI_Bcast(intervals, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, status, mpi_error_code)

	!calculate interval size

	interval_range = b - a
	
	local_intervals = intervals / worldsize
	
	local_range = (b - a) / intervals
	
	local_a = a + myrank * local_intervals * local_range

	call omp_set_num_threads(10)

        local_sum = (f(local_a) + f(local_a + local_range * local_intervals)) / 2.0

		counts = 0
		
	!$OMP PARALLEL DO PRIVATE(x) REDUCTION(+:local_sum)
	
	do counts = 0, local_intervals
	
		x = local_a + counts * local_range
		
		local_sum = local_sum + f(x)
	
	end do
	
	!$OMP END PARALLEL DO

        local_sum = local_sum * local_range

        print *, myrank, local_sum

        call MPI_Reduce(local_sum, global_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, mpi_error_code)

	
	end if
!end

	finish = MPI_Wtime()

	CALL MPI_Finalize(mpi_error_code)


end program trapMpiOmp

real function f(x)

	real(kind=8) :: x
	
	f = x
	
	return
	
end function f
