! this header contains the interface blocks to call
! functions of the cuda runtime from Fortran 90

 interface
   subroutine cudamallochost(hostptr, size) bind (c)
     use iso_c_binding
     type (c_ptr), intent(out) :: hostptr
     integer (c_int), intent(in) :: size
   end subroutine cudamallochost
 end interface




