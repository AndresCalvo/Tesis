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

  ! mathematical constants
  complex(cdp), parameter ::  im      = (0.0_dp, 1.0_dp)        ! imaginary unit
  real(dp), parameter     ::  halfpi  = 1.57079632679489_dp      ! pi/2
  real(dp), parameter     ::  pi      = 3.14159265358979_dp      ! pi
  real(dp), parameter     ::  pi2     = 6.28318530717959_dp      ! 2*pi

  ! SI constants (used only for parameter conversion, not in propagation)
  real(dp), parameter     ::  hbar_SI   = 1.054571817E-34_dp     ! ħ = h/(2*pi)  [J·s]
  real(dp), parameter     ::  eV_to_J   = 1.602176634E-19_dp     ! eV -> Joules
  real(dp), parameter     ::  c_SI      = 299792458.0_dp         ! speed of light [m/s]
  real(dp), parameter     ::  lambda0   = 245.0E-9_dp            ! carrier wavelength [m]
  real(dp), parameter     ::  omega0_SI = pi2 * c_SI / lambda0   ! ~7.682e15 rad/s

  ! ---------------------------------------------------------------
  ! Normalized units: omega0 = 1, k0 = 1, c = 1
  !   Time unit  : 1/omega0_SI  ~  1.3e-16 s  (sub-optical-cycle)
  !   Length unit: 1/k0_SI      =  lambda0/(2*pi) ~ 39 nm
  ! ---------------------------------------------------------------
  real(dp), parameter     :: c     = 1.0_dp   ! normalized speed of light
  real(dp), parameter     :: mu0   = 1.0_dp   ! normalized permeability
  real(dp), parameter     :: eps0  = 1.0_dp   ! normalized permittivity
  real(dp), parameter     :: omega0 = 1.0_dp  ! normalized carrier frequency
  real(dp), parameter     :: k0    = 1.0_dp   ! normalized wavenumber
  real(dp), parameter     :: v_g   = 1.0_dp   ! group velocity (= c in vacuum)

  ! Grid parameters
  real(dp), parameter     :: Tmax  = 40.0_dp            ! half-window in normalized time (~5 pulse widths)
  integer(i4b), parameter :: Nt    = 512                 ! time-grid points (larger -> smaller dt -> more stable)
  real(dp), parameter     :: dt    = 2.0_dp * Tmax / Nt ! time step
  real(dp), parameter     :: dz    = 1.0E-2_dp           ! z step (normalized)
  real(dp), parameter     :: dr    = 5.0E-1_dp           ! r step (normalized)
  real(dp), parameter     :: Rmax  = 50.0_dp             ! waveguide radius (normalized)
  real(dp), parameter     :: Lmax  = 5.0_dp              ! propagation length (normalized)
  integer(i4b), parameter :: Nz    = floor(Lmax / dz) + 1
  integer(i4b), parameter :: Nr    = floor(Rmax / dr) + 1
  real(dp), parameter     :: beam_waist = 100.0_dp       ! w^2: exp(-r^2/beam_waist), w=10, z_R=50
  real(dp), parameter     :: dw    = pi / Tmax           ! frequency step

  ! Atomic decay rate (normalized to omega0)
  real(dp), parameter     :: decay = 0.74_dp / omega0_SI * omega0_SI  ! kept as given; units: normalized

  ! Xenon energy levels [eV] -> normalized frequencies omega_i/omega0_SI
  real(dp), parameter, dimension(7) :: EnergyLevels = &
      [ 6.19921_dp, 8.43651_dp, 9.68565_dp, 9.81955_dp, 9.93349_dp, 10.401_dp, 11.0151_dp ]  ! [eV]
  real(dp), parameter, dimension(7) :: FrequencyLevels = &
      EnergyLevels * eV_to_J / hbar_SI   ! absolute frequencies [rad/s]
  ! Normalized detuning: Delta_i = (omega_i - omega0_SI) / omega0_SI
  real(dp), parameter, dimension(7) :: Delta = FrequencyLevels / omega0_SI - 1.0_dp

  real(dp), parameter, dimension(7) :: DecayRates = &
      [ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ]

  ! Dipole-moment coupling matrix u_ij [cm, e implicit]
  ! Physical dipole = e * u_cm * 1e-2  [C·m]
  ! Values from literature for Xe transitions.
  real(dp), dimension(7,7), parameter :: u = reshape( &
    [ 0.0_dp, 3.4E-9_dp,  0.0_dp,     0.0_dp,    0.0_dp,    -2.6E-9_dp,  0.0_dp,     &  ! col 1
      0.0_dp, 0.0_dp,     1.35E-8_dp, 8.0E-9_dp, 9.9E-9_dp, 0.0_dp,     -4.1E-8_dp, &  ! col 2
      0.0_dp, 0.0_dp,     0.0_dp,     0.0_dp,    0.0_dp,    -9.2E-9_dp,  0.0_dp,     &  ! col 3
      0.0_dp, 0.0_dp,     0.0_dp,     0.0_dp,    0.0_dp,    -4.7E-9_dp,  0.0_dp,     &  ! col 4
      0.0_dp, 0.0_dp,     0.0_dp,     0.0_dp,    0.0_dp,    -2.8E-8_dp,  0.0_dp,     &  ! col 5
      0.0_dp, 0.0_dp,     0.0_dp,     0.0_dp,    0.0_dp,     0.0_dp,    -1.85E-8_dp, &  ! col 6
      0.0_dp, 0.0_dp,     0.0_dp,     0.0_dp,    0.0_dp,     0.0_dp,     0.0_dp    ], &  ! col 7
       shape=[7, 7])

  ! Reference field amplitude and unit conversions
  real(dp), parameter :: E0       = 1.0E9_dp              ! peak field amplitude [V/m]
  real(dp), parameter :: e_charge = 1.602176634E-19_dp     ! electron charge [C]
  real(dp), parameter :: cm_to_m  = 1.0E-2_dp             ! cm -> m

  ! Normalised coupling matrix (dimensionless, ~O(1e-3)):
  !   u_norm_ij = (e * u_cm_ij[cm] * cm_to_m * E0) / (hbar_SI * omega0_SI)
  ! Rabi frequency in normalised units = u_norm_ij * A_normalised
  ! With A~1 and u_norm~6e-3, the time-step driving term ~ 6e-3 * dt ~ 1e-3  (safe)
  real(dp), dimension(7,7), parameter :: u_norm = &
      (e_charge * cm_to_m * E0 / (hbar_SI * omega0_SI)) * u

  ! Atomic number density [m^-3]
  real(dp), parameter :: NDensity = 2.7E10_dp

  ! Normalized polarization prefactor.
  ! In SI: P = 2*N * d * Re(rho), with d = e * u_cm * cm_to_m.
  ! The propagation eq. has: source ~ mu0 * omega^2 * P / (2*k0*(k0*vg - omega)).
  ! In normalized units (mu0=1, k0=1, omega0=1, c=1), the polarization
  ! must be expressed in the SAME normalized system as the field.
  ! The conversion factor is:
  !   P_norm = P_SI / (eps0_SI * E0)     [polarization in E0 units]
  ! where eps0_SI = 8.854187817e-12 F/m.
  ! Combined: P_norm = 2*N*d*Re(rho) / (eps0_SI*E0)
  !                  = 2*N*(e*u_cm*cm_to_m)/(eps0_SI*E0) * Re(rho)
  ! But u_norm_ij = (e*cm_to_m*E0)/(hbar_SI*omega0_SI) * u_ij,
  ! so e*u_cm*cm_to_m = u_norm_ij * hbar_SI * omega0_SI / E0.
  ! Thus: P_norm = 2*N* [u_norm * hbar_SI * omega0_SI / E0] / (eps0_SI*E0) * Re(rho)
  !              = 2*N*hbar_SI*omega0_SI / (eps0_SI*E0**2) * u_norm * Re(rho)
  real(dp), parameter :: eps0_SI  = 8.854187817E-12_dp   ! vacuum permittivity [F/m]
  real(dp), parameter :: pol_prefactor = &
      2.0_dp * NDensity * hbar_SI * omega0_SI / (eps0_SI * E0**2)

  ! Coordinate grids (allocated at runtime via setup)
  real(dp), dimension(Nt) :: t(Nt), freq(Nt)
  real(dp), dimension(Nz) :: z(Nz)
  real(dp), dimension(Nr) :: r(Nr)

  ! Field arrays
  complex(cdp), dimension(Nt, Nr) :: Afield, Afieldp1, Afield_f, Afield_f_p1, Pfield_f
  complex(cdp), dimension(Nt, Nr, 7, 7) :: rho, rabi
  complex(cdp), dimension(Nt)     :: initial_a_t
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
      complex(dp), intent(in)  :: arr(m)
      real(dp), dimension(2*m) :: complex_to_real(2*m)
      integer :: i
      do i = 1, m
        complex_to_real(2*i-1) = real(arr(i))
        complex_to_real(2*i)   = aimag(arr(i))
      end do
      return
    end function complex_to_real

    function real_to_complex(arr, m)
      implicit none
      integer(i4b), intent(in) :: m
      real(dp), intent(in)     :: arr(2*m)
      complex(cdp)             :: real_to_complex(m)
      integer :: i
      do i = 1, m
        real_to_complex(i) = cmplx(arr(2*i-1), arr(2*i), kind=cdp)
      end do
      return
    end function real_to_complex

    function sech(x)
      implicit none
      real(dp), intent(in) :: x
      real(dp)             :: sech
      sech = 1.0_dp / cosh(x)
      return
    end function sech

    function soliton(x, y)
      use declare, only : im
      implicit none
      complex(cdp) :: soliton
      real(dp), intent(in) :: x, y
      soliton = 4.0_dp * ( cosh(3.0_dp*y) + 3.0_dp*exp(im*4.0_dp*x)*cosh(y) ) * exp(im*x) &
                        / ( cosh(4.0_dp*y) + 4.0_dp*cosh(2.0_dp*y) + 3.0_dp*cos(4.0_dp*x) )
      return
    end function soliton

end module globals
!_____________________________________________________________________
subroutine setup
  use nrtype
  use globals
  use declare, only : Nt, Nr, Nz, t, r, z, freq, Tmax, Lmax, Rmax, pi2, Afield, initial_a_t, beam_waist, dt
  implicit none
  integer(i4b) :: i

  ! Frequency grid (FFT order: positive then negative)
  do i = 1, Nt/2
    freq(i)        =  pi2 * dble(i-1)      / (2.0_dp * Tmax)
    freq(i + Nt/2) =  pi2 * dble(i-1-Nt/2) / (2.0_dp * Tmax)
  end do

  ! Time grid: symmetric around 0
  do i = 1, Nt
    t(i) = (i-1) * 2.0_dp * Tmax / Nt - Tmax
    initial_a_t(i) = cmplx(sech(t(i)), 0.0_dp, kind=cdp)
  end do

  ! z grid
  do i = 1, Nz
    z(i) = (i-1) * Lmax / real(Nz-1, dp)
  end do

  ! r grid and initial Gaussian beam profile
  do i = 1, Nr
    r(i) = (i-1) * Rmax / real(Nr-1, dp)
    Afield(:, i) = initial_a_t(:) * exp( -(r(i)**2 / beam_waist) )
  end do

end subroutine setup
!_____________________________________________________________________
program main
  use nrtype
  use declare
  use globals
  implicit none
  integer :: i, j, ii, count
  integer, parameter :: Nwrite = 10           ! number of z snapshots to save
  real(dp), dimension(2*Nt) :: temp1, temp2

  call setup
  write(*, *) 'Setup complete.'
  write(*, '(A,I6,A,I6,A,I6)') '  Grid: Nt=', Nt, '  Nr=', Nr, '  Nz=', Nz
  write(*, '(A,F10.4,A,F10.6)') '  dt=', dt, '  dz=', dz

  ! Write initial pulse to file
  open(9, file='initial_pulse.dat')
  do ii = 1, Nr
    do j = 1, Nt
      write(9, *) t(j), r(ii), real(Afield(j, ii))**2
    end do
  end do
  write(9, *); write(9, *)
  close(9)

  write(*, *) 'Starting propagation...'

  open(10, file='results.dat')
  open(11, file='results2.dat')
  open(12, file='results3.dat')
  do i = 1, Nz

    ! --- Write snapshot ---
    count = mod(i, max(1, Nz/Nwrite))
    if (count == 0 .or. i == 1 .or. i == Nz) then
      do j = 1, Nt
        write(10, *) t(j), z(i), real(Afield(j, 1))**2
        write(11, *) t(j), z(i), real(Afield(j, Nr/2))**2
        write(12, *) t(j), z(i), real(Afield(j, Nr))**2
      end do
      write(10, *); write(10, *)
      write(11, *); write(11, *)
      write(12, *); write(12, *)
    end if

    ! ============================================================
    ! Symmetric split-step propagation:
    !   1. Half-step free dispersion in frequency domain
    !   2. Full-step diffraction + nonlinear (density matrix, polarization)
    !      in TIME domain (no 1/omega_env singularity)
    !   3. Half-step free dispersion in frequency domain
    ! ============================================================

    ! -- Step 1: FFT → half-step dispersion → IFFT --
    do ii = 1, Nr
      temp2 = complex_to_real(Afield(:, ii), Nt)
      call four1(temp2, Nt, 1)
      Afield_f(:, ii) = real_to_complex(temp2, Nt)
    end do
    ! Apply half-step phase: exp(i * freq(j) * dz/2)
    do j = 1, Nt
      Afield_f(j, :) = Afield_f(j, :) * exp(im * freq(j) * dz * 0.5_dp)
    end do
    do ii = 1, Nr
      temp2 = complex_to_real(Afield_f(:, ii), Nt)
      call four1(temp2, Nt, -1)
      Afield(:, ii) = real_to_complex(temp2, Nt) / real(Nt, dp)
    end do

    ! -- Step 2: Diffraction + nonlinear in TIME domain --
    ! 2a. Crank-Nicolson implicit diffraction step (unconditionally stable)
    !     Solves: dA/dz = (i/(2k0)) * nabla_perp^2 A
    call diffraction_step_CN(Afield, dz)

    ! 2b. Density matrix evolution
    call matrixDensityEvol(Afield, rho)
    ! 2c. Polarization
    call TotalPolarization(rho, Pfield)
    ! 2d. Add polarization source (small perturbation, Euler is fine here)
    Afield(:,:) = Afield(:,:) &
      - (im * dz * mu0 * omega0**2 / (2.0_dp * k0)) * Pfield(:,:)

    ! -- Step 3: FFT → half-step dispersion → IFFT --
    do ii = 1, Nr
      temp2 = complex_to_real(Afield(:, ii), Nt)
      call four1(temp2, Nt, 1)
      Afield_f(:, ii) = real_to_complex(temp2, Nt)
    end do
    do j = 1, Nt
      Afield_f(j, :) = Afield_f(j, :) * exp(im * freq(j) * dz * 0.5_dp)
    end do
    do ii = 1, Nr
      temp2 = complex_to_real(Afield_f(:, ii), Nt)
      call four1(temp2, Nt, -1)
      Afield(:, ii) = real_to_complex(temp2, Nt) / real(Nt, dp)
    end do

  end do

  close(10)
  close(11)
  close(12)
  write(*, *) 'Done. Output written to results.dat'

end program main
!_____________________________________________________________________
subroutine RabiFreq(a, rabi_out)
  use nrtype
  use declare, only : Nt, Nr, u_norm
  implicit none
  complex(cdp), dimension(Nt, Nr),       intent(in)  :: a
  complex(cdp), dimension(Nt, Nr, 7, 7), intent(out) :: rabi_out
  ! Normalised Rabi frequency: Omega_norm_ij = u_norm_ij * A_norm
  ! u_norm_ij = e * u_cm_ij * cm_to_m * E0 / (hbar_SI * omega0_SI)  (~6e-3 for strongest lines)
  ! A_norm ~ 1 at pulse peak  ->  Omega_norm ~ 6e-3  (well within stability region)
  rabi_out = 0.0_dp
  rabi_out(:, :, 1, 2) = u_norm(1,2) * a(:,:)
  rabi_out(:, :, 1, 6) = u_norm(1,6) * a(:,:)
  rabi_out(:, :, 2, 3) = u_norm(2,3) * a(:,:)
  rabi_out(:, :, 2, 4) = u_norm(2,4) * a(:,:)
  rabi_out(:, :, 2, 5) = u_norm(2,5) * a(:,:)
  rabi_out(:, :, 2, 7) = u_norm(2,7) * a(:,:)
  rabi_out(:, :, 3, 6) = u_norm(3,6) * a(:,:)
  rabi_out(:, :, 4, 6) = u_norm(4,6) * a(:,:)
  rabi_out(:, :, 5, 6) = u_norm(5,6) * a(:,:)
  rabi_out(:, :, 6, 7) = u_norm(6,7) * a(:,:)
end subroutine RabiFreq
!_____________________________________________________________________
subroutine matrixDensityEvol(a, rho_out)
  !
  ! Solves the Liouville-von Neumann equations for the density matrix
  ! using forward Euler for DIAGONAL (population) elements, and
  ! EXACT PHASE ROTATION for OFF-DIAGONAL (coherence) elements.
  !
  ! For the off-diagonals, the free-oscillation part is solved exactly:
  !   rho_ij(t+dt) = rho_ij(t) * exp((i*Omega_ij - gamma_ij/2)*dt)
  !                + [driving term] * dt
  ! This is unconditionally stable regardless of dt or detuning magnitude.
  !
  use nrtype
  use declare, only : dt, decay, omega0, im, Nr, Nt, Delta
  implicit none
  complex(cdp), dimension(Nt, Nr),       intent(in)  :: a
  complex(cdp), dimension(Nt, Nr, 7, 7), intent(out) :: rho_out

  complex(cdp), dimension(Nt, Nr, 7, 7) :: rabi
  integer :: i
  real(dp) :: g12, g16, g23, g24, g25, g27, g36, g46, g56, g67
  ! Phase-rotation factors for each coherence: exp((i*Delta_ij + i*omega0 - gamma/2)*dt)
  complex(cdp) :: ph12, ph16, ph23, ph24, ph25, ph27, ph36, ph46, ph56, ph67

  ! Decay rates per transition (all set to 'decay' for now)
  g12 = decay;  g16 = decay;  g23 = decay;  g24 = decay;  g25 = decay
  g27 = decay;  g36 = decay;  g46 = decay;  g56 = decay;  g67 = decay

  ! Precompute (scalar) phase-rotation factors — computed once, applied each step
  ! For rho_ij the free-oscillation exponent is (i*(Delta_i - Delta_j + omega0) - gamma_ij/2)
  ! where Delta here is the normalized detuning of each level from omega0.
  ! NOTE: Delta(1)=detuning of level 1, etc.  For ground state (level 1), Delta(1)≈0.224.
  ph12  = exp( (im*(Delta(2) - Delta(1) + omega0) - 0.5_dp*g12 ) * dt )
  ph16  = exp( (im*(Delta(6) - Delta(1) + omega0) - 0.5_dp*(g16+g36+g46+g56)) * dt )
  ph23  = exp( (im*(Delta(3) - Delta(2) + omega0) - 0.5_dp*(g12+g23)) * dt )
  ph24  = exp( (im*(Delta(4) - Delta(2) + omega0) - 0.5_dp*(g12+g24)) * dt )
  ph25  = exp( (im*(Delta(5) - Delta(2) + omega0) - 0.5_dp*(g12+g25)) * dt )
  ph27  = exp( (im*(Delta(7) - Delta(2) + omega0) - 0.5_dp*(g12+g27+g67)) * dt )
  ph36  = exp( (im*(Delta(6) - Delta(3) + omega0) - 0.5_dp*(g16+g23+g36+g46+g56)) * dt )
  ph46  = exp( (im*(Delta(6) - Delta(4) + omega0) - 0.5_dp*(g16+g24+g36+g46+g56)) * dt )
  ph56  = exp( (im*(Delta(6) - Delta(5) + omega0) - 0.5_dp*(g16+g25+g36+g46+g56)) * dt )
  ph67  = exp( (im*(Delta(7) - Delta(6) + omega0) - 0.5_dp*(g16+g27+g36+g46+g56+g67)) * dt )

  call RabiFreq(a, rabi)

  ! --- Initial conditions (all atoms in ground state, no coherences) ---
  rho_out(1, :, :, :) = (0.0_dp, 0.0_dp)
  rho_out(1, :, 1, 1) = 1.0_dp

  ! --- Time integration ---
  do i = 1, Nt - 1

    ! -- Diagonal (populations): forward Euler --
    rho_out(i+1, :, 1, 1) = rho_out(i, :, 1, 1) &
      + ( decay * ( rho_out(i,:,2,2) + rho_out(i,:,6,6) ) &
        + 2.0_dp * ( conjg(rabi(i,:,1,2)) * aimag(rho_out(i,:,1,2)) &
                   + conjg(rabi(i,:,1,6)) * aimag(rho_out(i,:,1,6)) ) ) * dt

    rho_out(i+1, :, 2, 2) = rho_out(i, :, 2, 2) &
      + ( decay * ( -rho_out(i,:,2,2) + rho_out(i,:,3,3) + rho_out(i,:,4,4) &
                  + rho_out(i,:,5,5) + rho_out(i,:,7,7) ) &
        - 2.0_dp * ( conjg(rabi(i,:,1,2)) * aimag(rho_out(i,:,1,2)) &
                   - conjg(rabi(i,:,2,3)) * aimag(rho_out(i,:,2,3)) &
                   - conjg(rabi(i,:,2,4)) * aimag(rho_out(i,:,2,4)) &
                   - conjg(rabi(i,:,2,5)) * aimag(rho_out(i,:,2,5)) &
                   - conjg(rabi(i,:,2,7)) * aimag(rho_out(i,:,2,7)) ) ) * dt

    rho_out(i+1, :, 3, 3) = rho_out(i, :, 3, 3) &
      + ( decay * ( -rho_out(i,:,3,3) + rho_out(i,:,6,6) ) &
        - 2.0_dp * ( conjg(rabi(i,:,2,3)) * aimag(rho_out(i,:,2,3)) &
                   - conjg(rabi(i,:,3,6)) * aimag(rho_out(i,:,3,6)) ) ) * dt

    rho_out(i+1, :, 4, 4) = rho_out(i, :, 4, 4) &
      + ( decay * ( -rho_out(i,:,4,4) + rho_out(i,:,6,6) ) &
        - 2.0_dp * ( conjg(rabi(i,:,2,4)) * aimag(rho_out(i,:,2,4)) &
                   - conjg(rabi(i,:,4,6)) * aimag(rho_out(i,:,4,6)) ) ) * dt

    rho_out(i+1, :, 5, 5) = rho_out(i, :, 5, 5) &
      + ( decay * ( -rho_out(i,:,5,5) + rho_out(i,:,6,6) ) &
        - 2.0_dp * ( conjg(rabi(i,:,2,5)) * aimag(rho_out(i,:,2,5)) &
                   - conjg(rabi(i,:,5,6)) * aimag(rho_out(i,:,5,6)) ) ) * dt

    rho_out(i+1, :, 6, 6) = rho_out(i, :, 6, 6) &
      + ( decay * ( -4.0_dp*rho_out(i,:,6,6) + rho_out(i,:,7,7) ) &
        - 2.0_dp * ( conjg(rabi(i,:,1,6)) * aimag(rho_out(i,:,1,6)) &
                   + conjg(rabi(i,:,3,6)) * aimag(rho_out(i,:,3,6)) &
                   + conjg(rabi(i,:,4,6)) * aimag(rho_out(i,:,4,6)) &
                   + conjg(rabi(i,:,5,6)) * aimag(rho_out(i,:,5,6)) &
                   - conjg(rabi(i,:,6,7)) * aimag(rho_out(i,:,6,7)) ) ) * dt

    rho_out(i+1, :, 7, 7) = rho_out(i, :, 7, 7) &
      + ( decay * ( -2.0_dp*rho_out(i,:,7,7) ) &
        - 2.0_dp * ( conjg(rabi(i,:,2,7)) * aimag(rho_out(i,:,2,7)) &
                   + conjg(rabi(i,:,6,7)) * aimag(rho_out(i,:,6,7)) ) ) * dt

    ! -- Off-diagonal (coherences): exact phase rotation + Euler driving term --
    ! Formula: rho_ij(t+dt) = rho_ij(t) * exp_factor
    !                        + i * Omega_ij * (rho_jj - rho_ii) * dt
    ! The exp_factor = exp((i*freq_ij - gamma_ij/2)*dt) is precomputed above.

    rho_out(i+1, :, 1, 2) = rho_out(i, :, 1, 2) * ph12 &
      + im * rabi(i, :, 1, 2) * (-rho_out(i,:,1,1) + rho_out(i,:,2,2)) * dt

    rho_out(i+1, :, 1, 6) = rho_out(i, :, 1, 6) * ph16 &
      + im * rabi(i, :, 1, 6) * (-rho_out(i,:,1,1) + rho_out(i,:,6,6)) * dt

    rho_out(i+1, :, 2, 3) = rho_out(i, :, 2, 3) * ph23 &
      + im * rabi(i, :, 2, 3) * (-rho_out(i,:,2,2) + rho_out(i,:,3,3)) * dt

    rho_out(i+1, :, 2, 4) = rho_out(i, :, 2, 4) * ph24 &
      + im * rabi(i, :, 2, 4) * (-rho_out(i,:,2,2) + rho_out(i,:,4,4)) * dt

    rho_out(i+1, :, 2, 5) = rho_out(i, :, 2, 5) * ph25 &
      + im * rabi(i, :, 2, 5) * (-rho_out(i,:,2,2) + rho_out(i,:,5,5)) * dt

    rho_out(i+1, :, 2, 7) = rho_out(i, :, 2, 7) * ph27 &
      + im * rabi(i, :, 2, 7) * (-rho_out(i,:,2,2) + rho_out(i,:,7,7)) * dt

    rho_out(i+1, :, 3, 6) = rho_out(i, :, 3, 6) * ph36 &
      + im * rabi(i, :, 3, 6) * (-rho_out(i,:,3,3) + rho_out(i,:,6,6)) * dt

    rho_out(i+1, :, 4, 6) = rho_out(i, :, 4, 6) * ph46 &
      + im * rabi(i, :, 4, 6) * (-rho_out(i,:,4,4) + rho_out(i,:,6,6)) * dt

    rho_out(i+1, :, 5, 6) = rho_out(i, :, 5, 6) * ph56 &
      + im * rabi(i, :, 5, 6) * (-rho_out(i,:,5,5) + rho_out(i,:,6,6)) * dt

    rho_out(i+1, :, 6, 7) = rho_out(i, :, 6, 7) * ph67 &
      + im * rabi(i, :, 6, 7) * (-rho_out(i,:,6,6) + rho_out(i,:,7,7)) * dt

  end do

end subroutine matrixDensityEvol
!_____________________________________________________________________
subroutine four1(data, nn, isign)
  ! Numerical Recipes FFT. isign=+1 forward, isign=-1 inverse (NOT normalized).
  implicit none
  integer  :: isign, nn
  double precision :: data(2*nn)
  integer  :: i, istep, j, m, mmax, n
  double precision :: tempi, tempr, theta, wi, wpi, wpr, wr, wtemp
  n = 2*nn
  j = 1
  do i = 1, n, 2
    if (j > i) then
      tempr    = data(j);   tempi    = data(j+1)
      data(j)  = data(i);   data(j+1)= data(i+1)
      data(i)  = tempr;     data(i+1)= tempi
    end if
    m = n/2
    do while ((m >= 2) .and. (j > m))
      j = j - m;  m = m/2
    end do
    j = j + m
  end do
  mmax = 2
  do while (n > mmax)
    istep = 2*mmax
    theta = 6.28318530717959d0 / (isign*mmax)
    wpr   = -2.d0*sin(0.5d0*theta)**2
    wpi   = sin(theta)
    wr    = 1.d0;  wi = 0.d0
    do m = 1, mmax, 2
      do i = m, n, istep
        j     = i + mmax
        tempr = sngl(wr)*data(j)   - sngl(wi)*data(j+1)
        tempi = sngl(wr)*data(j+1) + sngl(wi)*data(j)
        data(j)   = data(i)   - tempr
        data(j+1) = data(i+1) - tempi
        data(i)   = data(i)   + tempr
        data(i+1) = data(i+1) + tempi
      end do
      wtemp = wr
      wr = wr*wpr - wi*wpi + wr
      wi = wi*wpr + wtemp*wpi + wi
    end do
    mmax = istep
  end do
  return
end subroutine four1
!_____________________________________________________________________
subroutine TotalPolarization(rho_in, pol)
  use nrtype
  use declare, only : pol_prefactor, Nt, Nr, u_norm
  implicit none
  complex(cdp), dimension(Nt, Nr, 7, 7), intent(in)  :: rho_in
  complex(cdp), dimension(Nt, Nr),       intent(out) :: pol
  ! P_norm = pol_prefactor * sum_ij u_norm_ij * Re(rho_ij)
  ! pol_prefactor = 2*N*hbar_SI*omega0_SI / (eps0_SI * E0^2)   ~ O(1e-16)
  ! u_norm ~ O(1e-3), so P_norm ~ O(1e-19 * rho), vanishingly small feedback.
  pol(:,:) = pol_prefactor * &
    ( u_norm(1,2) * real(rho_in(:,:,1,2), dp) + u_norm(1,6) * real(rho_in(:,:,1,6), dp) &
    + u_norm(2,3) * real(rho_in(:,:,2,3), dp) + u_norm(2,4) * real(rho_in(:,:,2,4), dp) &
    + u_norm(2,5) * real(rho_in(:,:,2,5), dp) + u_norm(2,7) * real(rho_in(:,:,2,7), dp) &
    + u_norm(3,6) * real(rho_in(:,:,3,6), dp) + u_norm(4,6) * real(rho_in(:,:,4,6), dp) &
    + u_norm(5,6) * real(rho_in(:,:,5,6), dp) + u_norm(6,7) * real(rho_in(:,:,6,7), dp) )
end subroutine TotalPolarization
!_____________________________________________________________________
subroutine ddaperp(arr, darr)
  ! Radial part of transverse Laplacian: (1/r) d/dr (r d/dr) arr
  use nrtype
  use declare, only : Nt, r, Nr, dr
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in)  :: arr
  complex(cdp), dimension(Nt, Nr), intent(out) :: darr
  integer :: i
  ! On-axis (r=0): use L'Hopital limit -> 4*(arr(:,2)-arr(:,1))/dr^2
  darr(:, 1) = 4.0_dp * ( arr(:,2) - arr(:,1) ) / (dr**2)
  do i = 2, Nr - 1
    darr(:, i) = ( arr(:,i+1) - 2.0_dp*arr(:,i) + arr(:,i-1) ) / (dr**2) &
               + ( arr(:,i+1) - arr(:,i-1) ) / (2.0_dp * dr * r(i))
  end do
  ! Boundary at r=Rmax (Dirichlet: arr beyond Rmax = 0)
  darr(:, Nr) = ( -2.0_dp*arr(:,Nr) + arr(:,Nr-1) ) / (dr**2) &
              + ( -arr(:,Nr-1) ) / (2.0_dp * dr * r(Nr))
end subroutine ddaperp
!_____________________________________________________________________
subroutine PolarizationTerms(FPfield, omega_abs, RHPT)
  ! 𝒫 = −μ₀(ω₀+ω_env)²P  where omega_abs = omega0 + omega_env
  ! So 𝒫 = −μ₀ * omega_abs² * P
  use nrtype
  use declare, only : Nt, Nr, mu0
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in)  :: FPfield
  real(dp),     dimension(Nt),     intent(in)  :: omega_abs
  complex(cdp), dimension(Nt, Nr), intent(out) :: RHPT
  integer(i4b) :: i
  do i = 1, Nt
    RHPT(i, :) = -mu0 * omega_abs(i)**2 * FPfield(i, :)
  end do
end subroutine PolarizationTerms
!_____________________________________________________________________
subroutine PropagationStep(A_n, P_n, A_np1)
  !-----------------------------------------------------------------
  ! Split-step propagator for the Brabec-Krausz-type equation:
  !   2i(k₀ − Ω/vg) ∂A/∂z = 𝒜 + 𝒫
  ! where Ω = ω₀ + ω_env is the ABSOLUTE frequency and ω_env = freq(i)
  ! is the baseband (envelope) Fourier frequency.
  !
  ! With k₀=1, vg=1, ω₀=1:
  !   k₀vg − Ω = 1 − (1+ω_env) = −ω_env
  !
  ! Step 1 — Free dispersion (exact):
  !   ∂A/∂z = −i(−ω_env) A = i·ω_env·A
  !   Exact: A → A × exp(i·ω_env·dz)
  !   This is the standard free-propagation phase for the envelope.
  !
  ! Step 2 — Diffraction + polarization (perturbative Euler):
  !   ΔA = −(i·dz/2)(vg/(−ω_env)) × (−Δ⊥A + 𝒫)
  !       = (i·dz/2)(1/ω_env) × (−Δ⊥A + 𝒫)
  !   At ω_env≈0 (DC of envelope): skip (Δ⊥A=0 for DC, 𝒫 negligible)
  !-----------------------------------------------------------------
  use nrtype
  use declare, only : Nt, Nr, im, dz, v_g, k0, freq, omega0
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(in)  :: A_n   ! Fourier field at z
  complex(cdp), dimension(Nt, Nr), intent(in)  :: P_n   ! Fourier polarization at z
  complex(cdp), dimension(Nt, Nr), intent(out) :: A_np1 ! Fourier field at z+dz
  complex(cdp), dimension(Nt, Nr) :: RHPT, perpFA
  real(dp),     dimension(Nt)     :: omega_abs   ! absolute frequency
  real(dp) :: denom, w_env
  integer(i4b) :: i

  ! Absolute frequency = carrier + envelope Fourier freq
  omega_abs = omega0 + freq

  ! Diffraction of the field in frequency domain (acts on r, at each ω)
  call ddaperp(A_n, perpFA)

  ! Polarization source term: 𝒫 = −μ₀ Ω² P̃
  call PolarizationTerms(P_n, omega_abs, RHPT)

  do i = 1, Nt
    w_env = freq(i)                         ! envelope frequency
    denom = k0*v_g - omega_abs(i)            ! = −w_env

    ! -- Step 1: exact phase propagator for free dispersion --
    ! exp(−i·denom·dz) = exp(i·w_env·dz)
    A_np1(i, :) = A_n(i, :) * exp(-im * denom * dz)

    ! -- Step 2: diffraction + polarization correction --
    ! Skip DC of envelope (ω_env≈0): correction diverges but is physically zero
    if (abs(w_env) < 1.0E-6_dp) cycle

    A_np1(i, :) = A_np1(i, :) &
      + (im * dz / 2.0_dp) * (v_g / denom) * perpFA(i, :) &   ! +Δ⊥A
      - (im * dz / 2.0_dp) * (v_g / denom) * RHPT(i, :)        ! −𝒫

  end do

end subroutine PropagationStep
!_____________________________________________________________________
subroutine diffraction_step_CN(A, hz)
  !-----------------------------------------------------------------
  ! Crank-Nicolson step for the paraxial diffraction equation:
  !   dA/dz = (i / (2*k0)) * nabla_perp^2 A
  ! where nabla_perp^2 = d^2/dr^2 + (1/r) d/dr  (cylindrical coord.)
  !
  ! CN scheme: (I - alpha/2 * L) A^{n+1} = (I + alpha/2 * L) A^n
  !   alpha = i*hz/(2*k0),  L = tridiagonal radial Laplacian
  !
  ! Solved per time-point using the Thomas algorithm for complex
  ! tridiagonal systems. Unconditionally stable.
  !-----------------------------------------------------------------
  use nrtype
  use declare, only : Nt, Nr, im, k0, dr, r
  implicit none
  complex(cdp), dimension(Nt, Nr), intent(inout) :: A
  real(dp), intent(in) :: hz

  complex(cdp) :: alpha
  ! Tridiagonal coefficients
  complex(cdp), dimension(Nr) :: diag_lhs, rhs_vec
  complex(cdp), dimension(Nr) :: sub, sup   ! sub-diagonal and super-diagonal
  ! Work arrays for Thomas algorithm
  complex(cdp), dimension(Nr) :: cp, dp_th
  complex(cdp) :: m_fac
  integer :: it, ir
  real(dp) :: dr2, a_sub, a_diag, a_sup

  alpha = im * hz / (2.0_dp * k0)
  dr2 = dr**2

  do it = 1, Nt

    ! -- Build RHS = (I + alpha/2 * L) * A --
    ! and LHS tridiagonal = (I - alpha/2 * L)

    ! On-axis (ir=1, r=0): L'Hôpital -> L = 4*(A(2)-A(1))/dr^2
    ! Stencil: L(1) = 4/dr^2 * [A(2) - A(1)]
    a_diag = -4.0_dp / dr2
    a_sup  =  4.0_dp / dr2

    diag_lhs(1) = (1.0_dp, 0.0_dp) - (alpha/2.0_dp) * a_diag
    sup(1)      =                   - (alpha/2.0_dp) * a_sup
    rhs_vec(1)  = ((1.0_dp, 0.0_dp) + (alpha/2.0_dp) * a_diag) * A(it, 1) &
                + (alpha/2.0_dp) * a_sup * A(it, 2)

    ! Interior points (ir=2..Nr-1)
    do ir = 2, Nr - 1
      a_sub  = 1.0_dp/dr2 - 1.0_dp/(2.0_dp*dr*r(ir))
      a_diag = -2.0_dp/dr2
      a_sup  = 1.0_dp/dr2 + 1.0_dp/(2.0_dp*dr*r(ir))

      sub(ir)      =                   - (alpha/2.0_dp) * a_sub
      diag_lhs(ir) = (1.0_dp, 0.0_dp) - (alpha/2.0_dp) * a_diag
      sup(ir)      =                   - (alpha/2.0_dp) * a_sup

      rhs_vec(ir)  = (alpha/2.0_dp) * a_sub  * A(it, ir-1) &
                   + ((1.0_dp, 0.0_dp) + (alpha/2.0_dp) * a_diag) * A(it, ir) &
                   + (alpha/2.0_dp) * a_sup  * A(it, ir+1)
    end do

    ! Boundary (ir=Nr): Dirichlet A(Nr+1) = 0
    a_sub  = 1.0_dp/dr2 - 1.0_dp/(2.0_dp*dr*r(Nr))
    a_diag = -2.0_dp/dr2
    ! a_sup = 0 (ghost point = 0)

    sub(Nr)      =                   - (alpha/2.0_dp) * a_sub
    diag_lhs(Nr) = (1.0_dp, 0.0_dp) - (alpha/2.0_dp) * a_diag
    rhs_vec(Nr)  = (alpha/2.0_dp) * a_sub  * A(it, Nr-1) &
                 + ((1.0_dp, 0.0_dp) + (alpha/2.0_dp) * a_diag) * A(it, Nr)

    ! -- Thomas algorithm (complex tridiagonal solver) --
    ! Forward sweep
    cp(1) = sup(1) / diag_lhs(1)
    dp_th(1) = rhs_vec(1) / diag_lhs(1)
    do ir = 2, Nr
      m_fac = diag_lhs(ir) - sub(ir) * cp(ir-1)
      cp(ir) = sup(ir) / m_fac        ! sup(Nr)=0 so cp(Nr)=0
      dp_th(ir) = (rhs_vec(ir) - sub(ir) * dp_th(ir-1)) / m_fac
    end do

    ! Back substitution
    A(it, Nr) = dp_th(Nr)
    do ir = Nr-1, 1, -1
      A(it, ir) = dp_th(ir) - cp(ir) * A(it, ir+1)
    end do

  end do

end subroutine diffraction_step_CN
