!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! merge_fields
! J. Cipriani
! August / 2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Accepts 3 arguments
!     1. Source NetCDF file
!     2. Destination NetCDF file
!     3. Variable name
! The requested variable must be in the source NetCDF file, but not necessarily
! in the destination NetCDF file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program merge_fields

    use netcdf

    implicit none

    !Arguments
    character(len = 500) :: source_file  ! DA output (source)
    character(len = 500) :: dest_file    ! init.nc output (destination)
    character(len = 500) :: var          ! Reqeusted variable/field

    ! Variables
    integer  :: ncid, ncid_out
    integer  :: nPoints, nRecords, nz
    integer  :: numdims, numatts
    integer  :: varid, varid_out
    integer, dimension(nf90_max_var_dims) :: dimids
    logical  :: dest_has_var
    real, dimension(:,:,:), allocatable :: tmp_3d
    real, dimension(:,:), allocatable   :: tmp_2d
    real, dimension(:), allocatable     :: tmp_1d

    ! only used if creating a new variable
    integer  :: xtype
    integer  :: i
    character(len = 500) :: dimname
    integer, dimension(nf90_max_var_dims) :: dimids_dest


    ! Process arguments
    call getarg(1, source_file)    ! DA output (source)
    call getarg(2, dest_file)      ! init.nc output (destination)
    call getarg(3, var)            ! Requested variable/field
    if (trim(var) .eq. "xtime") then
        write(0, *) 'xtime not supported'
        stop 3
    end if

    ! Open NetCDF Files
    call check( nf90_open(trim(source_file), NF90_NOWRITE, ncid) )
    call check( nf90_open(trim(dest_file), NF90_WRITE, ncid_out) )

    ! Get IDs for requested variable
    call check( nf90_inq_varid(ncid, trim(var), varid) )

    ! Determine if we're writing a new variable to the destination
    dest_has_var = .true.
    if (nf90_inq_varid(ncid_out, trim(var), varid_out) .ne. NF90_NOERR) then
        dest_has_var = .false.
    end if

    ! Get dimensions/attributes for source variable
    call check( nf90_inquire_variable(ncid, varid, ndims = numdims, natts = numatts) )
    call check( nf90_inquire_variable(ncid, varid, dimids = dimids) )
    if (.not.dest_has_var) then

        ! Get some more information necessary for variable creation
        !   xtype of the source data
        !   dimension information (need to figure out what the source dimensions
        !       are and get the proper IDs for those dimensions in the dest.)
        !   dimension IDs of the dest. data (to create var using proper IDs)
        call check( nf90_inquire_variable(ncid, varid, xtype = xtype) )

        ! Get dimension names of the source data
        ! Get dimension ids of destination (Using dimension names of source)
        do i = 1, numdims
            call check( nf90_inquire_dimension(ncid, dimids(i), dimname) )
            call check( nf90_inq_dimid(ncid_out, trim(dimname), dimids_dest(i) ) )
        end do

        call check( nf90_redef(ncid_out) )  ! Go into redefine mode to create variable
        call check( nf90_def_var(ncid_out, trim(var), xtype, dimids_dest(1:numdims), varid_out) )
        call check( nf90_enddef(ncid_out) )

    end if

    ! Perform merging (copying from source to destination)
    ! 3 dimensions
    if (numdims .eq. 3) then
        call check( nf90_inquire_dimension(ncid, dimids(1), len = nz) )
        call check( nf90_inquire_dimension(ncid, dimids(2), len = nPoints) )
        call check( nf90_inquire_dimension(ncid, dimids(3), len = nRecords) )
        allocate(tmp_3d(nz,nPoints,nRecords))
        call check( nf90_get_var(ncid, varid, tmp_3d) )
        call check( nf90_put_var(ncid_out, varid_out, tmp_3d) )
        deallocate(tmp_3d)
    ! 2 dimensions
    else if (numdims .eq. 2) then
        call check( nf90_inquire_dimension(ncid, dimids(1), len = nPoints) )
        call check( nf90_inquire_dimension(ncid, dimids(2), len = nRecords) )
        allocate(tmp_2d(nPoints,nRecords))
        call check( nf90_get_var(ncid, varid, tmp_2d) )
        call check( nf90_put_var(ncid_out, varid_out, tmp_2d) )
        deallocate(tmp_2d)
    ! 1 dimension
    else
        call check( nf90_inquire_dimension(ncid, dimids(1), len = nPoints) )
        allocate(tmp_1d(nPoints))
        call check( nf90_get_var(ncid, varid, tmp_1d) )
        call check( nf90_put_var(ncid_out, varid_out, tmp_1d) )
        deallocate(tmp_1d)
    end if

    ! Close netCDF files
    call check( nf90_close(ncid) )
    call check( nf90_close(ncid_out) )


contains

    subroutine check(status)
        integer, intent (in) :: status
        if(status .ne. nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop 2
        end if
    end subroutine check

end program merge_fields
