module units

  implicit none
  private

  public :: toev, tokcal, toang

  ! --- Mathematical constants ---
  double precision, parameter, public :: PI      = 3.141592653589793d0
  double precision, parameter, public :: PI2     = PI**2
  double precision, parameter, public :: PI12    = sqrt(PI)

  ! --- Physical constants / unit conversions (CODATA) ---

  !> Hartree energy in eV.
  !> CODATA 2018: https://physics.nist.gov/cgi-bin/cuu/Value?hrev
  double precision, parameter, public :: hartree2ev = 27.211386245981d0

  !> Hartree energy in Joule.
  !> https://physics.nist.gov/cgi-bin/cuu/Value?hr
  double precision, parameter, public :: hartree2joule =   4.3597447222060d-18

  !> Avogadro constant
  !> https://physics.nist.gov/cgi-bin/cuu/Value?na
  double precision, parameter ::  Na =  6.02214076d23

  !> cal definition
  !> https://www.nist.gov/pml/special-publication-811
  double precision, parameter ::  cal =  4.184d0

  double precision, parameter, public :: hartree2kcal = hartree2joule * Na / cal


  !> Bohr radius in Angstrom.
  !> CODATA 2022: https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
  double precision, parameter, public :: bohr2ang = 0.529177210544d0


contains

  pure function toev(e_hartree) result(e_ev)
    double precision, intent(in) :: e_hartree
    double precision :: e_ev
    e_ev = e_hartree * hartree2ev
  end function toev


  pure function tokcal(e_hartree) result(e_kcal)
        double precision, intent(in) :: e_hartree
        double precision :: e_kcal

        e_kcal = e_hartree * hartree2kcal
        
  end function tokcal


  pure function toang(r_bohr) result(r_ang)
    double precision, intent(in) :: r_bohr
    double precision :: r_ang
    r_ang = r_bohr * bohr2ang
  end function toang

end module units
