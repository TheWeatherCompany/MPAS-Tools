PROGRAM merge_fields

  use netcdf

  IMPLICIT NONE

  !Arguments
  character(len = 500)                                 :: INCOMING_FCST
  character(len = 500)                                 :: INTERP_INPUT
  character(len = 500)                                 :: var

  integer                                              :: ncid, ncid_out
  integer                                              :: nPoints, records, nz
  integer                                              :: numdims, numatts
  integer                                              :: varid, varid_out
  integer, dimension(nf90_max_var_dims)                :: dimids
  
  real, dimension(:,:,:), allocatable                  :: tmp_3d
  real, dimension(:,:), allocatable                    :: tmp_2d
  real, dimension(:), allocatable                      :: tmp_1d


! Get file name as argument
 call getarg(1,INTERP_INPUT)
 call getarg(2,INCOMING_FCST)
 call getarg(3,var)

! Open netCDF files
 call check( nf90_open(trim(INTERP_INPUT), NF90_NOWRITE, ncid) )
 call check( nf90_open(trim(INCOMING_FCST), NF90_WRITE, ncid_out) )

! Get varid's and VARIABLES from netCDF
 call check( nf90_inq_varid(ncid, trim(var), varid) )
 call check( nf90_inq_varid(ncid_out, trim(var), varid_out) )
 call check( nf90_inquire_variable(ncid, varid, ndims = numdims, natts = numatts) )
 call check( nf90_inquire_variable(ncid, varid, dimids = dimids) )

 if (trim(var) .ne. "xtime") then
   ! 3 dimensions
   if (numdims .eq. 3) then
       call check( nf90_inquire_dimension(ncid, dimids(1), len = nz) )
       call check( nf90_inquire_dimension(ncid, dimids(2), len = nPoints) )
       call check( nf90_inquire_dimension(ncid, dimids(3), len = records) )
       allocate(tmp_3d(nz,nPoints,records))
       call check( nf90_get_var(ncid, varid, tmp_3d) )
       call check( nf90_put_var(ncid_out,varid_out,tmp_3d) )
       deallocate(tmp_3d)
   ! 2 dimensions
   else if (numdims .eq. 2) then
       call check( nf90_inquire_dimension(ncid, dimids(1), len = nPoints) )
       call check( nf90_inquire_dimension(ncid, dimids(2), len = records) )
       allocate(tmp_2d(nPoints,records))
       call check( nf90_get_var(ncid, varid, tmp_2d) )
       call check( nf90_put_var(ncid_out,varid_out,tmp_2d) )
       deallocate(tmp_2d)
   ! 1 dimension
   else
       call check( nf90_inquire_dimension(ncid, dimids(1), len = nPoints) )
       allocate(tmp_1d(nPoints))
       call check( nf90_get_var(ncid, varid, tmp_1d) )
       call check( nf90_put_var(ncid_out,varid_out,tmp_1d) )
       deallocate(tmp_1d)
   end if
end if


! Close netCDF files
 call check( nf90_close(ncid) )
 call check( nf90_close(ncid_out) )

! "Check" subroutine
   CONTAINS

       SUBROUTINE check(status)
         integer, intent (in) :: status
         if(status .ne. nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       END SUBROUTINE check


END PROGRAM merge_fields

