    subroutine resetseed()
    integer :: size,seed_size
    integer, allocatable :: seed_array(:)
    write (6,*) 'Resetting seed'
    call random_seed( size = seed_size )
    allocate( seed_array( seed_size ) )
    seed_array = (/ (i, i=1,seed_size) /)
    call random_seed( put = seed_array )
    end subroutine resetseed
