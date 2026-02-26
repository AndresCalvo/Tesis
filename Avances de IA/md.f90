module nrtype
  ! Defines numerical types
  implicit none
  integer, parameter :: i4b = selected_int_kind(9)
  integer, parameter :: i2b = selected_int_kind(4)
  integer, parameter :: i1b = selected_int_kind(1)
  integer, parameter :: dp  = kind(0.d0)
  integer, parameter :: cdp = kind((0.0_dp, 1.0_dp))
end module nrtype
!_________________________________________________________________   
module declare
  use nrtype
  implicit none

  ! mathematical and physical constants
  complex(cdp), parameter ::  im = (0.0_dp, 1.0_dp)        ! imaginary unit
  real(dp), parameter     ::  halfpi = 1.57079632679489_dp ! pi/2
  real(dp), parameter     ::  pi = 3.14159265358979_dp     ! pi
  real(dp), parameter     ::  pi2 = 6.28318530717959_dp    ! 2*pi
  real(dp), parameter     ::  hbar = 1.054571817E-34_dp    ! ħ = h/(2*pi)  [J·s]
  real(dp), parameter     ::  eV_to_J = 1.602176634E-19_dp ! eV to Joules
  real(dp), parameter     ::  c_SI = 299792458.0_dp        ! speed of light [m/s] (for conversions only)
  real(dp), parameter     ::  lambda0 = 245.0E-9_dp        ! carrier wavelength [m]
  real(dp), parameter     ::  omega0_SI = pi2*c_SI/lambda0 ! carrier frequency [rad/s], ~7.682e15
  ! --- Normalized units: omega0=1, k0=1, c=1 ---
  real(dp), parameter     :: c    = 1.0_dp                 ! normalized speed of light
  real(dp), parameter     :: mu0  = 1.0_dp                 ! normalized permeability
  real(dp), parameter     :: eps0 = 1.0_dp                 ! normalized permittivity

  ! waveguide and problem parameteres
  real(dp), parameter     :: Tmax = 40.0_dp            ! half-width of time window [1/omega0]; covers ~5 pulse widths
  integer(i4b), parameter :: Nt = 256                 !size of T window
  real(dp), parameter     :: dz = 1E-1_dp              ! size of z step
  real(dp), parameter     :: dt = Tmax / Nt            ! time step size
  real(dp), parameter     :: dr = 1E-1_dp              ! size of r step
  real(dp), parameter     :: Rmax = 15.0_dp            ! radius of the waveguide
  real(dp), parameter     :: Lmax = 5.0_dp             ! total fiber length [m]
  integer(i4b), parameter :: Nz = floor(Lmax / dz) + 1 ! steps in z
  integer(i4b), parameter :: Nr = floor(Rmax / dr) + 1 ! steps in r
  real(dp), parameter     :: k0 = 1.0_dp               ! wavenumber k_0 (normalized: k0=omega0/c=1)
  real(dp), parameter     :: decay = 0.74_dp           ! decay rate gamma in the Liouville equations.
  real(dp), parameter     :: omega0 = 1.0_dp           ! center frequency (normalized)
  real(dp), parameter     :: beam_waist = 1.0_dp       ! initial gaussian beam waist [normalized]
  real(dp), parameter     :: v_g = 1.0_dp              ! group velocity = c (normalized)
  real(dp), parameter     :: dw = pi / Tmax            ! delta omega, frequency step size

  ! atomic properties
  real(dp), parameter               :: NDensity = 2.7E10_dp ! atomic density
  real(dp), parameter, dimension(7) :: EnergyLevels = [ 6.19921, 8.43651, 9.68565, 9.81955, 9.93349, 10.401, 11.0151] ! [eV]
  real(dp), parameter, dimension(7) :: DecayRates = [ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ] ! decay gammas
  ! FrequencyLevels in SI [rad/s]: eV -> J -> rad/s
  real(dp), parameter, dimension(7) :: FrequencyLevels = EnergyLevels * eV_to_J / hbar
  ! Delta: normalized detuning from omega0. (omega_i - omega0_SI) / omega0_SI
  real(dp), parameter, dimension(7) :: Delta = FrequencyLevels / omega0_SI - 1.0_dp

  ! coordinate grids
  real(dp), dimension(Nt) :: t(Nt), freq(Nt) ! time, space and frequency domains
  real(dp), dimension(Nz) :: z(Nz)           ! z (propagation direction) domain
  real(dp), dimension(Nr) :: r(Nr)           ! radial domain
  real(dp), dimension(7,7), parameter :: u = reshape( &  !u dipole moment [cm] only u_ij where i<j. symmetry is solved in equations
    [ 0.0_dp, 3.4E-9_dp, 0.0_dp    , 0.0_dp   , 0.0_dp   , -2.6E-9_dp , 0.0_dp     , &   ! row 1
      0.0_dp, 0.0_dp   , 1.35E-8_dp, 8.0E-9_dp, 9.9E-9_dp, 0.0_dp     , -4.1E-8_dp , &   ! row 2
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -9.2E-9_dp , 0.0_dp     , &   ! row 3
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -4.7E-9_dp , 0.0_dp     , &   ! row 4
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -2.8E-8_dp , 0._dp      , &   ! row 5
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , 0.0_dp     , -1.85E-8_dp, &   ! row 6
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , 0.0_dp     , 0.0_dp     ], &  ! row 7
       shape=[7, 7])

! physical variables (field, envelope, density matrix elements, etc)
  complex(cdp), dimension(Nt, Nr) :: Afield, Afieldp1, Afield_f, Afield_f_p1, Pfield_f
  complex(cdp), dimension(Nt, Nr, 7, 7) :: rho, rabi
  complex(cdp), dimension(Nt) :: initial_a_t
  complex(cdp), dimension(Nt) :: tdata, diff, sol
  complex(cdp), dimension(Nt, Nr) :: Pfield
end module declare
!_____________________________________________________________________
module globals
  use nrtype
  implicit none
  contains

    function complex_to_real(arr, m)
      implicit none
      integer(i4b), intent(in) :: m
      complex(dp), intent(in) :: arr(m)
      real(dp), dimension(2 * m) :: complex_to_real(2 * m)
      integer i

      do i=1, m
        complex_to_real(2 * i - 1)=real(arr(i))
        complex_to_real(2 * i)=aimag(arr(i))
      end do

      return
    end function complex_to_real
    !_______________________________________________________________
    function real_to_complex(arr, m)
      implicit none
      integer(i4b), intent(in) :: m
      real(dp), intent(in) :: arr(2 * m)
      complex(cdp) real_to_complex(m)
      integer ::  i

      do i=1, m
        real_to_complex(i) = cmplx(arr(2 * i - 1), arr(2 * i), kind=cdp)
      end do

      return
    end function real_to_complex
    !___________________________________________________________________
    function sech(x)
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: sech
      sech=1. / cosh(x)

      return
    end function sech
    !___________________________________________________________________
    function soliton(x, y) !x for space, y for time
      use declare, only : im
      implicit none
      complex(cdp) :: soliton
      real(dp), intent(in) :: x, y

      soliton = 4 * ( cosh(3 * y) + 3 * exp(im * 4 * x) * cosh(y) ) * exp(im * x) &
                              / ( cosh(4 * y) + 4 * cosh(2 * y) + 3 * cos(4 * x) )

      return
    end function soliton
    !___________________________________________________________________
    function lin_sus(w)
      use declare
      implicit none
      
      real(dp), intent(in) :: w
      real(dp) :: lin_sus
      lin_sus = ( NDensity / ( eps0 * hbar ) )
      return
    end function lin_sus
  end module globals
!_____________________________________________________________________
subroutine setup !sets range arrays for time and space coordinates
  use nrtype
  use globals
  use declare, only : Nt, Nr, Nz, t, r, z, freq, Tmax, Lmax, Rmax, pi2, Afield, initial_a_t, beam_waist
  implicit none
  integer(i4b) :: i

  do i=1, Nt/2
    freq(i) = pi2 * dble(i-1) / (2 * Tmax) !positive and zero frequencies
    freq(i + Nt / 2)=pi2 * dble(i-1-Nt/2)/(2*Tmax) !neg and max/min feq
  end do

  do i=1, Nt
    t(i) = (i-1) * 2.0_dp * Tmax / Nt - Tmax

    initial_a_t(i) = cmplx(sech(t(i)), 0.0_dp, kind=cdp)

  end do

  do i=1, Nz
    z(i) = (i - 1) * Lmax / Nz
  end do

  !DO i=1, Nr
  !  r(i) = (i - 1) * Rmax / Nr
  !END DO

  do i=1, Nr
    r(i) = (i - 1) * Rmax / Nr
    Afield(:, i) = initial_a_t(:) * exp( -( r(i) ** 2 / beam_waist ) )
  end do

end subroutine setup
!___________________________________________
program main
  use nrtype
  use declare
  use globals
  implicit none
  integer :: i, j, k, ii, count
  real(dp), dimension(2*Nt) :: temp, temp1, temp2

  call setup
  write(*, *) 'Finished setup'
  write(*, *) 'Starting density matrix evolution'

  open(9, file='initial_pulse.dat')

  ! Write initial pulse
  do k=1, Nr
    do j=1, Nt
      write(9, *) t(j), r(k), real(Afield(j, k)) ** 2
    end do
  end do
  write(9, *)
  write(9, *)
  close(9)

  ! main loop
  open(10, file='results.dat')
  do i=1, Nz !main loop over z
    count=mod(i, 50)
    !if (count == 0) then !loop for writing into file
        do j=1, Nt
          write(10, *) t(j), z(i), real(Afield(j, 1)) ** 2
        end do
      write(10, *)
      write(10, *)
    !end if
    
    ! main program
    call matrixDensityEvol(Afield, rho)
    call TotalPolarization(rho, Pfield)
  
    do ii=1, Nr
      temp1 = complex_to_real(Pfield(:, ii), Nt)
      call four1(temp1, Nt, 1)
      Pfield_f(:, ii) = real_to_complex(temp1, Nt)
    
      temp2 = complex_to_real(Afield(:, ii), Nt)
      call four1(temp2, Nt, 1)
      Afield_f(:, ii) = real_to_complex(temp2, Nt)
    end do

    call PropagationStep(Afield_f, Pfield_f, Afield_f_p1)

    do ii=1, Nr
      temp2 = complex_to_real(Afield_f_p1(:, ii), Nt)
      call four1(temp2, Nt, -1)
      Afield(:, ii) = real_to_complex(temp2, Nt) / Nt
    end do
  end do

  close(10)
  write(*, *) 'Finished writing data'
end program main
!________________________________________________________________________
subroutine RabiFreq(a, rabi) !Takes electric field to yield rabbi frequencies array per spatial point
  use nrtype
  use declare, only : Nt, Nr, hbar, u
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in) :: a ! slowly varying complex amplitude
  complex(cdp), dimension(Nt, Nr, 7, 7), intent(out) :: rabi ! rabbi frequencies array

  rabi(:, :, 1, 2) = u(1, 2) * a(:, :) / hbar
  rabi(:, :, 1, 6) = u(1, 6) * a(:, :) / hbar
  rabi(:, :, 2, 3) = u(2, 3) * a(:, :) / hbar
  rabi(:, :, 2, 4) = u(2, 4) * a(:, :) / hbar
  rabi(:, :, 2, 5) = u(2, 5) * a(:, :) / hbar
  rabi(:, :, 2, 7) = u(2, 7) * a(:, :) / hbar
  rabi(:, :, 3, 6) = u(3, 6) * a(:, :) / hbar
  rabi(:, :, 4, 6) = u(4, 6) * a(:, :) / hbar
  rabi(:, :, 5, 6) = u(5, 6) * a(:, :) / hbar
  rabi(:, :, 6, 7) = u(6, 7) * a(:, :) / hbar

end subroutine RabiFreq
!________________________________________________________________________
subroutine matrixDensityEvol(a, r) !solves for rho in time with A field as input
  use nrtype
  use declare, only : dt, decay, omega0, im, Nr, Nt, Delta
  implicit none
  complex(cdp), dimension(Nt, Nr),intent(in) :: a
  complex(cdp), dimension(Nt, Nr, 7, 7), intent(out) :: r ! r for rho, density matrix
  ! COMPLEX(CDP), DIMENSION(Nr, Nt, 7, 7) :: rabi ! rabi values
  complex(cdp), dimension(Nt, Nr, 7, 7) :: rabi
  integer :: i
  real(dp) :: Tr12, Tr16, Tr23, Tr24, Tr25, Tr27, Tr36, Tr46, Tr56, Tr67

  Tr12 = decay
  Tr16 = decay
  Tr23 = decay
  Tr24 = decay
  Tr25 = decay
  Tr27 = decay
  Tr36 = decay
  Tr46 = decay
  Tr56 = decay
  Tr67 = decay

  call RabiFreq(a, rabi)

  r(1, :, 1, 1) = 1.0_dp
  r(1, :, 1, 2) = 0.0_dp
  r(1, :, 1, 6) = 0.0_dp
  r(1, :, 2, 3) = 0.0_dp
  r(1, :, 2, 4) = 0.0_dp
  r(1, :, 2, 5) = 0.0_dp
  r(1, :, 2, 7) = 0.0_dp
  r(1, :, 3, 6) = 0.0_dp
  r(1, :, 4, 6) = 0.0_dp
  r(1, :, 5, 6) = 0.0_dp
  r(1, :, 6, 7) = 0.0_dp

  do i=1, Nt - 1
    ! diagonal elements
    r(i+1, :, 1, 1) = r(i, :, 1, 1) + ( decay * ( r(i, :, 2, 2) + r(i, :, 6, 6) ) &
                                        + 2.0_dp * ( conjg(rabi(i, :, 1, 2)) * aimag(r(i, :, 1, 2)) &
                                                   + conjg(rabi(i, :, 1, 6)) * aimag(r(i, :, 1, 6)) ) ) * dt

    r(i+1, :, 2, 2) = r(i, :, 2, 2) + ( decay * ( -r(i, :, 2, 2) + r(i, :, 3, 3) + r(i, :, 4, 4) &
                                                  + r(i, :, 5, 5) + r(i, :, 7, 7) ) &
                                        - 2.0_dp * ( conjg(rabi(i, :, 1, 2)) * aimag(r(i, :, 1, 2)) &
                                                   - conjg(rabi(i, :, 2, 3)) * aimag(r(i, :, 2, 3)) &
                                                   - conjg(rabi(i, :, 2, 4)) * aimag(r(i, :, 2, 4)) &
                                                   - conjg(rabi(i, :, 2, 5)) * aimag(r(i, :, 2, 5)) &
                                                   - conjg(rabi(i, :, 2, 7)) * aimag(r(i, :, 2, 7)) ) ) * dt

    r(i+1, :, 3, 3) = r(i, :, 3, 3) + ( decay * ( -r(i, :, 3, 3) + r(i, :, 6, 6) ) &
                                        - 2.0_dp * ( conjg(rabi(i, :, 2, 3)) * aimag(r(i, :, 2, 3)) &
                                                   - conjg(rabi(i, :, 3, 6)) * aimag(r(i, :, 3, 6)) ) ) * dt

    r(i+1, :, 4, 4) = r(i, :, 4, 4) + ( decay * ( -r(i, :, 4, 4) + r(i, :, 6, 6) ) &
                                        - 2.0_dp * ( conjg(rabi(i, :, 2, 4)) * aimag(r(i, :, 2, 4)) &
                                                   - conjg(rabi(i, :, 4, 6)) * aimag(r(i, :, 4, 6)) ) ) * dt

    r(i+1, :, 5, 5) = r(i, :, 5, 5) + ( decay * ( -r(i, :, 5, 5) + r(i, :, 6, 6) ) &
                                        - 2.0_dp * ( conjg(rabi(i, :, 2, 5)) * aimag(r(i, :, 2, 5)) &
                                                   - conjg(rabi(i, :, 5, 6)) * aimag(r(i, :, 5, 6))) ) * dt

    r(i+1, :, 6, 6) = r(i, :, 6, 6) + ( decay * (-4.0 * r(i, :, 6, 6) + r(i, :, 7, 7) ) &
                                        - 2.0_dp * ( conjg(rabi(i, :, 1, 6)) * aimag(r(i, :, 1, 6)) &
                                                   + conjg(rabi(i, :, 3, 6)) * aimag(r(i, :, 3, 6)) &
                                                   + conjg(rabi(i, :, 4, 6)) * aimag(r(i, :, 4, 6)) &
                                                   + conjg(rabi(i, :, 5, 6)) * aimag(r(i, :, 5, 6)) &
                                                   - conjg(rabi(i, :, 6, 7)) * aimag(r(i, :, 6, 7)) ) ) * dt

    r(i+1, :, 7, 7) = r(i, :, 7, 7) + ( decay * ( -2.0 * r(i, :, 7, 7) ) &
                                       - 2.0_dp * ( conjg(rabi(i, :, 2, 7)) * aimag(r(i, :, 2, 7)) &
                                                   + conjg(rabi(i, :, 6, 7)) * aimag(r(i, :, 6, 7)) ) ) * dt
    ! non diagonal elements
    r(i+1, :, 1, 2) = r(i, :, 1, 2) + im * ( rabi(i, :, 1, 2) * ( -r(i, :, 1, 1) + r(i, :, 2, 2) ) &
                                            + r(i, :, 1, 2) * ( Delta(2) + omega0 - 0.5_dp * (Tr12) ) ) * dt

    r(i+1, :, 1, 6) = r(i, :, 1, 6) + im * ( rabi(i, :, 1, 6) * ( -r(i, :, 1, 1) + r(i, :, 6, 6) ) &
                                            + r(i, :, 1, 6) * ( Delta(6) + omega0 - 0.5_dp * (Tr16 + Tr36 + Tr46 + Tr56)) ) * dt

    r(i+1, :, 2, 3) = r(i, :, 2, 3) + im * ( rabi(i, :, 2, 3) * ( -r(i, :, 2, 2) + r(i, :, 3, 3) ) &
                                            + r(i, :, 2, 3) * ( Delta(2) + Delta(3) + omega0 - 0.5_dp * (Tr12 + Tr23)) ) * dt

    r(i+1, :, 2, 4) = r(i, :, 2, 4) + im * ( rabi(i, :, 2, 4) * ( -r(i, :, 2, 2) + r(i, :, 4, 4) ) &
                                            + r(i, :, 2, 4) * ( -Delta(2) + Delta(4) + omega0 - 0.5_dp * (Tr12 + Tr24)) ) * dt

    r(i+1, :, 2, 5) = r(i, :, 2, 5) + im * ( rabi(i, :, 2, 5) * ( -r(i, :, 2, 2) + r(i, :, 5, 5) ) &
                                            + r(i, :, 2, 5) * ( -Delta(2) + Delta(5) + omega0 - 0.5_dp * (Tr12 + Tr25)) ) * dt

    r(i+1, :, 2, 7) = r(i, :, 2, 7) + im * ( rabi(i, :, 2, 7) * ( -r(i, :, 2, 2) + r(i, :, 7, 7) ) &
                                            + r(i, :, 2, 7) * ( -Delta(2) + Delta(7) + omega0 &
                                                               - 0.5_dp * (Tr12 + Tr27 + Tr67)) ) * dt

    r(i+1, :, 3, 6) = r(i, :, 3, 6) + im * ( rabi(i, :, 3, 6) * ( -r(i, :, 3, 3) + r(i, :, 6, 6) ) &
                                            + r(i, :, 3, 6) * ( -Delta(3) + Delta(6) + omega0 &
                                                               - 0.5_dp * (Tr16 + Tr23 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(i+1, :, 4, 6) = r(i, :, 4, 6) + im * ( rabi(i, :, 4, 6) * ( -r(i, :, 4, 4) + r(i, :, 6, 6) ) &
                                            + r(i, :, 4, 6) * ( -Delta(4) + Delta(6) + omega0 &
                                                               - 0.5_dp * (Tr16 + Tr24 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(i+1, :, 5, 6) = r(i, :, 5, 6) + im * ( rabi(i, :, 5, 6) * ( -r(i, :, 5, 5) + r(i, :, 6, 6) ) &
                                            + r(i, :, 5, 6) * ( -Delta(5) + Delta(6) + omega0 &
                                                               - 0.5_dp * (Tr16 + Tr25 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(i+1, :, 6, 7) = r(i, :, 6, 7) + im * ( rabi(i, :, 6, 7) * ( -r(i, :, 6, 6) + r(i, :, 7, 7)) &
                                            + r(i, :, 6, 7) * ( -Delta(6) + Delta(7) + omega0 &
                                                               - 0.5_dp * (Tr16 + Tr27 + Tr36 + Tr46 + Tr56 + Tr67)) ) * dt
  end do

end subroutine matrixDensityEvol
!________________________________________________________________________
subroutine four1(data,nn,isign)
  implicit none
  integer isign,nn
  double precision data(2*nn)
  integer i,istep,j,m,mmax,n
  double precision tempi,tempr
  double precision theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
  j=1
  do i=1,n,2
    !This is the bit-reversal section of the routine.
    if(j.gt.i)then  
      tempr=data(j)
      !Exchange the two complex numbers.
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
    endif
    m=n/2
    do while ((m.ge.2).and.(j.gt.m))
      j=j-m
      m=m/2
    end do
    j=j+m
  enddo
  mmax=2
  !Here begins the Danielson-Lanczos section of the routine.
  do while (n.gt.mmax)
    !Outer loop executed log 2 nn times.
    istep=2*mmax
    theta=6.28318530717959d0/(isign*mmax)
    !Initialize for the trigonometric recurrence
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do m=1,mmax,2
      do i=m,n,istep
        j=i+mmax
        !This is the Danielson-Lanczos formula:
        tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
        tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
        data(j)=data(i)-tempr
        data(j+1)=data(i+1)-tempi
        data(i)=data(i)+tempr
        data(i+1)=data(i+1)+tempi
      enddo
      !Trigonometric recurrence
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    enddo
    mmax=istep
  end do
  return
end subroutine four1
!_________________________________________
subroutine TotalPolarization(rho, pol) ! takes rho (density matrix) to calculate total polarization
  use nrtype
  use declare, only : NDensity, Nt, Nr, im, omega0, u
  implicit none
  complex(cdp), dimension(Nr, Nt, 7, 7), intent(in) :: rho
  complex(cdp), dimension(Nr, Nt), intent(out) :: pol

  pol(:, :) = 2 * NDensity * ( u(1, 2) * real(rho(:, :, 1, 2), kind=dp) + u(1, 6) * real(rho(:, :, 1, 6), kind=dp) &
                             + u(2, 3) * real(rho(:, :, 2, 3), kind=dp) + u(2, 4) * real(rho(:, :, 2, 4), kind=dp) &
                             + u(2, 5) * real(rho(:, :, 2, 5), kind=dp) + u(2, 7) * real(rho(:, :, 2, 7), kind=dp) &
                             + u(3, 6) * real(rho(:, :, 3, 6), kind=dp) + u(4, 6) * real(rho(:, :, 4, 6), kind=dp) &
                             + u(5, 6) * real(rho(:, :, 5, 6), kind=dp) + u(6, 7) * real(rho(:, :, 6, 7), kind=dp) )

end subroutine TotalPolarization
!__________________________________
subroutine polarization_term(p, pterm)
  ! DEPRECATED. DO NOT USE.
  use nrtype
  use globals, only : complex_to_real, real_to_complex
  use declare, only : Nt, Nr, im, omega0, k0, mu0, freq
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in) :: p
  complex(cdp), dimension(Nt, Nr), intent(out) :: pterm
  complex(cdp), dimension(Nt, Nr) :: p_fourier_c
  real(cdp), dimension(2*Nt,Nr) :: p_fourier_r
  integer :: i
  
  do i=1, Nr 
    p_fourier_r(:, i) = complex_to_real(p(:, i), Nt)
  end do

  do i=1, Nr
    call four1(p_fourier_r(:, i), Nt, 1)
  end do

  do i=1, Nr
  p_fourier_c(:, i) = real_to_complex(p_fourier_r(:, i), Nt)
  end do

  do i=1, Nr
    pterm(:, i) = ( im * mu0 / (2 * k0) ) * ( omega0 * ( omega0 + 2 * freq(:) ) - freq(:) ** 2 ) * p_fourier_c(:, i)
  end do

  do i=1, Nr
    p_fourier_r(:, i) = complex_to_real(pterm(:, i), Nt)
  end do

  do i=1, Nr
    call four1(p_fourier_r(:, i), Nt, -1)
  end do

  do i=1, Nr
    pterm(:, i) = real_to_complex(p_fourier_r(:, i), Nt)
  end do

end subroutine polarization_term
!________________________________________________________________________
subroutine ddaperp(arr, darr) !radial part of laplacian for field
  use nrtype
  use declare, only : Nt, r, Nr, dr
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in) :: arr
  complex(cdp), dimension(Nt, Nr), intent(out) :: darr
  integer :: i

  darr(:, 1) = 4.0 * ( arr(:, 2) - arr(:, 1) ) / (dr ** 2)

  do i = 2, Nr - 1
    darr(:, i) = ( arr(:, i + 1) - 2.0 * arr(:, i) + arr(:, i - 1) ) / (dr ** 2) &
               + ( arr(:, i + 1) - arr(:, i - 1) ) / (2.0_dp * dr *r(i))
  end do

  darr(:, Nr) = ( - 2.0 * arr(:, Nr) + arr(:, Nr - 1) ) / (dr ** 2) &
               + (  - arr(:, Nr - 1) ) / (2.0_dp * dr *r(Nr))

end subroutine ddaperp
!_________________________________________________________________
subroutine radial_laplacian(arr, darr) !radial part of laplacian for field
  use nrtype
  use declare, only : r, Nr, dr
  implicit none
  complex(cdp), dimension(Nr), intent(in) :: arr
  complex(cdp), dimension(Nr), intent(out) :: darr
  integer :: i

  darr(1) = 4.0 * ( arr(2) - arr(1) ) / (dr ** 2)

  do i = 2, Nr - 1
    darr(i) = ( arr(i + 1) - 2.0 * arr(i) + arr(i - 1) ) / (dr ** 2) &
               + ( arr(i + 1) - arr(i - 1) ) / (2.0_dp * dr *r(i))
  end do

  darr(Nr) = ( - 2.0 * arr(Nr) + arr(Nr - 1) ) / (dr ** 2) &
               + (  - arr(Nr - 1) ) / (2.0_dp * dr *r(Nr))

end subroutine radial_laplacian
!_________________________________________________________________
subroutine RightHandFieldTerms(FAfield, omega, RHFT)
  use nrtype
  use declare, only : im, Nt, Nr, k0, omega0, dz, v_g, c
  implicit none
  integer(i4b) :: i
  complex(cdp), dimension(Nt, Nr) :: perpFA ! radial laplacian of FAfield
  complex(cdp), intent(in) :: FAfield(Nt, Nr) ! Fourier transform of Afield
  real(dp), dimension(Nt), intent(in) :: omega
  complex(cdp), intent(out) :: RHFT(Nt, Nr) ! Right Hand Field Terms of propagation equation
  
  call ddaperp(FAfield, perpFA)
  do i=1, Nt
    RHFT(i, :) = k0 * FAfield(i, :) ! (i)
    RHFT(i, :) = RHFT(i, :) - (2 * k0 / v_g) * omega(i) * FAfield(i, :) ! (iii)
    RHFT(i, :) = RHFT(i, :) + (1 / v_g ** 2) * (omega(i) ** 2) * FAfield(i, :) ! (v)
    RHFT(i, :) = RHFT(i, :) - perpFA(i, :) ! (vi)
    RHFT(i, :) = RHFT(i, :)  + ((omega0 / c) ** 2) * FAfield(i, :) ! (vii)
    RHFT(i, :) = RHFT(i, :)  - (2 * omega0 / (c ** 2)) * omega(i) * FAfield(i, :) ! (viii)
    RHFT(i, :) = RHFT(i, :) - omega(i) ** 2 * FAfield(i, :) ! (ix)
  end do
end subroutine RightHandFieldTerms
!_________________________________________________________________
subroutine PolarizationTerms(FPField, omega, RHPT)
  use nrtype
  use declare, only : im, Nt, Nr, hbar, u, k0, mu0, omega0, dz, v_g
  implicit none
  complex(cdp), intent(in) :: FPfield(Nt, Nr) ! Fourier transform of Polarization field (PField)
  real(dp), dimension(Nt),intent(in) :: omega
  complex(cdp), intent(out) :: RHPT(Nt, Nr) ! Right Hand Polarization Terms of propagation equation
  integer(i4b) :: i
  do i=1, Nt
    RHPT(i, :) = im * ((omega0 ** 2 * mu0) / (2*k0)) * FPfield(i, :) ! (i)
    RHPT(i, :) = RHPT(i, :) + im * (omega0 * mu0 / k0) * omega(i) * FPfield(i, :) ! (ii)
    RHPT(i, :) = RHPT(i, :) - im * (mu0 / (2 * k0)) * (omega(i) ** 2) * FPfield(i, :) ! (iii)
  end do
end subroutine PolarizationTerms
!_________________________________________________________________
subroutine PropagationStep(A_n, P_n, A_np1)
  use nrtype
  use declare, only : dt, Nt, Nr, im, dz, v_g, k0, freq
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in) :: A_n ! Fourier transform of Afield at step z_n
  complex(cdp), dimension(Nt, Nr), intent(in) :: P_n ! Fourier transform of Polarization field (PField)
  complex(cdp), dimension(Nt, Nr), intent(out) :: A_np1 ! Fourier transform of Afield at step z_n+dz
  complex(cdp), dimension(Nt, Nr) :: RHFT, RHPT
  real(dp), dimension(Nt) :: omega
  integer(i4b) :: i

  omega = freq  ! use frequency grid from declare module
  call RightHandFieldTerms(A_n, omega, RHFT)
  call PolarizationTerms(P_n, omega, RHPT)

  do i=1, Nt
    A_np1(i, :) = A_n(i, :) - (im * dz / 2) * (v_g / (k0 * v_g - omega(i))) * (RHFT(i, :) + RHPT(i, :))
  end do
end subroutine PropagationStep