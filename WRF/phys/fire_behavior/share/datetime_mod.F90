  module datetime_mod

    use, intrinsic :: iso_fortran_env, only : OUTPUT_UNIT

    implicit none

    private

    public :: datetime_t

    integer, parameter :: DATETIME_LEN = 19, MONTHS_IN_YEAR = 12, MINS_IN_AN_HOUR = 60, SECONDS_IN_A_MIN = 60, &
        HOURS_IN_A_DAY = 24
    integer, dimension(MONTHS_IN_YEAR), parameter :: DAYS_IN_MONTHS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

    type :: datetime_t 
      character (len = DATETIME_LEN) :: datetime
      real :: fractional_seconds
    contains
      procedure :: Add_seconds_int, Add_seconds_real
      generic, public :: Add_seconds => Add_seconds_int
      generic, public :: Add_seconds => Add_seconds_real
      procedure, public :: Calc_days_in_month => Calc_days_in_month
      procedure, public :: Get_datetime_as_ints => Get_datetime_as_ints
      procedure, public :: Get_dd => Get_dd
      procedure, public :: Get_frac_ss => Get_frac_ss
      procedure, public :: Get_hh => Get_hh
      procedure, public :: Get_min => Get_min
      procedure, public :: Get_mm => Get_mm
      procedure, public :: Get_ss => Get_ss
      procedure, public :: Get_yyyy => Get_yyyy
      procedure, public :: Is_leap => Is_leap
      procedure, public :: Print_datetime => Print_datetime
      procedure, public :: Update_datetime => Update_datetime
        ! Operators
      procedure, public :: Different
      generic :: operator(/=) => Different
      procedure, public :: Lower_or_equal
      generic :: operator(<=) => Lower_or_equal
      procedure, public :: Equal
      generic :: operator(==) => Equal
      procedure, public :: Lower
      generic :: operator(<) => Lower
      procedure, pass(this) :: Copy
      generic, public       :: assignment(=) => Copy
    end type datetime_t


    interface datetime_t
      module procedure Datetime_t_const
    end interface datetime_t

  contains

    pure subroutine Add_seconds_int (this, delta_seconds)

      implicit none

      class (datetime_t), intent (in out) :: this
      integer, intent (in) :: delta_seconds

      type (datetime_t) :: datetime
      integer :: yyyy, mm, dd, hh, mins, secs, days_in_this_month


      days_in_this_month = this%Calc_days_in_month ()
      call this%Get_datetime_as_ints (yyyy, mm, dd, hh, mins, secs)

      secs = secs + delta_seconds
      Loop_secs: do while (secs >= SECONDS_IN_A_MIN)
        secs = secs - 60
        mins = mins + 1
        if_mins: if (mins == MINS_IN_AN_HOUR) then
          mins = 0
          hh = hh + 1
          if_hh: if (hh == HOURS_IN_A_DAY) then
            hh = 0
            dd = dd + 1
            if_dd: if (dd == days_in_this_month + 1) then
              dd = 1
              mm = mm + 1
              if_mm: if (mm == MONTHS_IN_YEAR + 1) then
                mm = 1
                yyyy = yyyy + 1
                datetime = datetime_t (yyyy, mm, dd, hh, mins, secs)
              end if if_mm
              days_in_this_month = this%Calc_days_in_month ()
            end if if_dd
          end if if_hh
        end if if_mins
      end do Loop_secs

      call this%Update_datetime (yyyy, mm, dd, hh, mins, secs)

    end subroutine Add_seconds_int

    pure subroutine Add_seconds_real (this, delta_seconds)

      implicit none

      class (datetime_t), intent (in out) :: this
      real, intent (in) :: delta_seconds

      integer :: delta_seconds_int
      real :: delta_seconds_frac, frac_ss

      delta_seconds_int = int (delta_seconds)
      delta_seconds_frac = delta_seconds - delta_seconds_int

      frac_ss = this%Get_frac_ss ()
      frac_ss = frac_ss + delta_seconds_frac
      if (frac_ss >= 1.0) then
        delta_seconds_int = delta_seconds_int + 1
        frac_ss = frac_ss - 1.0
      end if

      this%fractional_seconds = frac_ss

      call this%Add_seconds_int (delta_seconds_int)

    end subroutine Add_seconds_real

    pure function Calc_days_in_month (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      integer :: mm

      mm = this%Get_mm ()
      return_value = DAYS_IN_MONTHS(mm)
      if (mm == 2 .and. this%Is_leap ()) return_value = 29

    end function Calc_days_in_month

    pure subroutine Copy (this, rhs)

      implicit none

      class (datetime_t), intent(in out) :: this
      class(datetime_t), intent(in) :: rhs

      this%datetime = rhs%datetime
      this%fractional_seconds = rhs%fractional_seconds

    end subroutine Copy

    pure function Datetime_t_const (yyyy, mm, dd, hh, minutes, ss, frac_ss) result (return_value)

      implicit none

      integer, intent (in) :: yyyy, mm, dd, hh, minutes, ss
      integer, intent (in), optional :: frac_ss
      type (datetime_t) :: return_value


      if (present (frac_ss)) then
        return_value%fractional_seconds = frac_ss
      else
        return_value%fractional_seconds = 0.0
      end if

      write (return_value%datetime, "(i4, 5(a,i2.2))") yyyy, '-', mm, '-', dd, '_', hh, ':', minutes, ':', ss

    end function Datetime_t_const

    pure function Different (lhs, rhs) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: lhs, rhs
      logical :: return_value


      if (lhs%datetime /= rhs%datetime .or. lhs%fractional_seconds /= rhs%fractional_seconds) then
        return_value = .true.
      else
        return_value = .false.
      end if

    end function Different

    pure function Equal (lhs, rhs) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: lhs, rhs
      logical :: return_value


      if (lhs%datetime == rhs%datetime .and. lhs%fractional_seconds == rhs%fractional_seconds) then
        return_value = .true.
      else
        return_value = .false.
      end if

    end function Equal

    pure subroutine Get_datetime_as_ints (this, yyyy, mm, dd, hh, minutes, seconds)

      implicit none

      class (datetime_t), intent (in) :: this
      integer, intent(out) :: yyyy, mm, dd, hh, minutes, seconds


      yyyy = this%Get_yyyy ()
      mm = this%Get_mm ()
      dd = this%Get_dd ()
      hh = this%Get_hh ()
      minutes = this%Get_min ()
      seconds = this%Get_ss ()

    end subroutine Get_datetime_as_ints

    pure function Get_dd (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read (this%datetime(9:10),'(i2.2)') return_value

    end function Get_dd

    pure function Get_frac_ss (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      real :: return_value

      return_value = this%fractional_seconds

    end function Get_frac_ss

    pure function Get_hh (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read (this%datetime(12:13),'(i2.2)') return_value

    end function Get_hh

    pure function Get_min (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read (this%datetime(15:16),'(i2.2)') return_value

    end function Get_min

    pure function Get_mm (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read (this%datetime(6:7),'(i2.2)') return_value

    end function Get_mm

    pure function Get_ss (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read (this%datetime(18:19),'(i2.2)') return_value

    end function Get_ss

    pure function Get_yyyy (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      integer :: return_value

      read(this%datetime(1:4),'(i4.4)') return_value

    end function Get_yyyy

    pure function Is_leap (this) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: this
      logical :: return_value

      integer :: yyyy

      yyyy = this%Get_yyyy ()
      if (mod (yyyy, 100) /= 0 .and. mod (yyyy, 4) == 0) then
        return_value = .true.
      elseif (mod(yyyy, 400) == 0) then
        return_value = .true.
      else
        return_value = .false.
      end if

    end function Is_leap

    pure function Lower (lhs, rhs) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: lhs, rhs
      logical :: return_value


      if (lhs%datetime < rhs%datetime) then
        return_value = .true.
      else
        if (lhs%datetime == rhs%datetime .and. lhs%fractional_seconds < rhs%fractional_seconds) then
          return_value = .true.
        else
          return_value = .false.
        end if
      end if

    end function Lower

    pure function Lower_or_equal (lhs, rhs) result (return_value)

      implicit none

      class (datetime_t), intent (in) :: lhs, rhs
      logical :: return_value


      if (lhs%datetime < rhs%datetime) then
        return_value = .true.
      else
        if (lhs%datetime == rhs%datetime .and. lhs%fractional_seconds <= rhs%fractional_seconds) then
          return_value = .true.
        else
          return_value = .false.
        end if
      end if

    end function Lower_or_equal

    subroutine Print_datetime (this)

      implicit none

      class (datetime_t), intent (in) :: this

      write (OUTPUT_UNIT, *) this%datetime, this%fractional_seconds

    end subroutine Print_datetime

    pure subroutine Update_datetime (this, yyyy, mm, dd, hh, minutes, ss, frac_ss)

      implicit none

      class (datetime_t), intent (in out) :: this
      integer, intent (in) :: yyyy, mm, dd, hh, minutes, ss
      real, intent (in), optional :: frac_ss

      if (present (frac_ss)) this%fractional_seconds = frac_ss
      write (this%datetime, "(i4, 5(a,i2.2))") yyyy, '-', mm, '-', dd, '_', hh, ':', minutes, ':', ss

    end subroutine Update_datetime

  end module datetime_mod

