module machine
  !! sets the precision if not set by a host model

  implicit none
  private

  public :: kind_phys, kind_sngl_prec, kind_dbl_prec, kind_io8

  integer, parameter :: kind_sngl_prec = 4, &
    kind_dbl_prec = 8

! physics single precision flag
#ifndef SINGLE_PREC 
  integer, parameter :: kind_phys = kind_dbl_prec
#else
  integer, parameter :: kind_phys = kind_sngl_prec
#endif
  integer, parameter :: kind_io8 = kind_phys

end module machine