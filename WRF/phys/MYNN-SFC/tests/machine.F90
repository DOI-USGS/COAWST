!> \file machine.f90
!!
module machine
!! \section arg_table_machine
!! \htmlinclude machine.html
!!
  implicit none
  
  integer, parameter :: kind_sngl_prec = 4
  integer, parameter :: kind_dbl_prec = 8
#ifdef __PGI
  integer, parameter :: kind_qdt_prec = 8
#else
  integer, parameter :: kind_qdt_prec = 16
#endif
  integer, parameter :: kind_integer = 4
  integer, parameter :: kind_logical = 4
  integer, parameter :: kind_io4 = kind_sngl_prec
  integer, parameter :: kind_ior = kind_dbl_prec
  integer, parameter :: kind_grid = kind_dbl_prec
  
  ! Physics single precision flag
#ifndef SINGLE_PREC
  integer, parameter :: kind_phys = kind_dbl_prec
#else
  integer, parameter :: kind_phys = kind_sngl_prec
#endif
  
  ! Note kind_io8 is not always 8 bytes
  integer, parameter :: kind_io8 = kind_phys
  
  ! Dynamics single precision flag
#ifdef OVERLOAD_R4
  integer, parameter :: kind_dyn = kind_sngl_prec
#else
  integer, parameter :: kind_dyn = kind_dbl_prec
#endif
  
  ! Machine precision to restrict dep
  real(kind=kind_phys), parameter :: mprec = 1.0e-12_kind_phys
  
  ! GRIB undefined value
  real(kind=kind_phys), parameter :: grib_undef = 9.99e20_kind_phys

end module machine