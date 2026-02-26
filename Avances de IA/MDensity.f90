MODULE nrtype
  ! Defines numerical types
  IMPLICIT NONE
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(1)
  INTEGER, PARAMETER :: DP  = KIND(0.D0)
  INTEGER, PARAMETER :: CDP = KIND((0.0_dp, 1.0_dp))
END MODULE nrtype
!_________________________________________________________________   
MODULE declare
  IMPLICIT NONE
  USE nrtype

  ! mathematical and physical constants
  COMPLEX(CDP), PARAMETER ::  im = (0.0_dp, 1.0_dp)        ! imaginary unit
  REAL(DP), PARAMETER     ::  halfpi = 1.57079632679489_dp ! pi/2
  REAL(DP), PARAMETER     ::  pi = 3.14159265358979_dp     ! pi
  REAL(DP), PARAMETER     ::  pi2 = 6.28318530717959_dp    ! 2*pi
  REAL(DP), PARAMETER     ::  hbar = 6.62607015E-34        ! reduced planks constant
  REAL(DP), PARAMETER     :: c = 299792458_dp              ! speed of light [m/s]
  REAL(DP), PARAMETER     :: mu0 = 1.2566370612720E-6_dp   ! permeability of free space
  REAL(DP), PARAMETER     :: eps0 = 8.854187817E-12_dp     ! permittivity of free space

  ! waveguide and problem parameteres
  REAL(DP), PARAMETER     :: Tmax = 10.0_dp            ! max time in frame [s]
  INTEGER(I4B), PARAMETER :: Nt = 1024                 !size of T window
  REAL(DP), PARAMETER     :: dz = 1E-3_dp              ! size of z step
  REAL(DP), PARAMETER     :: dt = Tmax / Nt            ! time step size
  REAL(DP), PARAMETER     :: dr = 1E-3_dp              ! size of r step
  REAL(DP), PARAMETER     :: Rmax = 15.0_dp            ! radius of the waveguide
  REAL(DP), PARAMETER     :: Lmax = 5.0_dp             ! total fiber length [m]
  INTEGER(I4B), PARAMETER :: Nz = floor(Lmax / dz) + 1 ! steps in z
  INTEGER(I4B), PARAMETER :: Nr = floor(Rmax / dr) + 1 ! steps in r
  REAL(DP), PARAMETER     :: k0 = 1.0_dp               ! wavenumber k_0 in envelope
  REAL(DP), PARAMETER     :: decay = 0.74_dp           ! decay rate gamma in the Liouville equations.
  REAL(DP), PARAMETER     :: omega0 = 1.0_dp           ! center frequency
  REAL(DP), PARAMETER     :: beam_waist = 1.0_dp       ! initial gaussian beam waist
  REAL(DP), PARAMETER     :: v_g = c                   ! group velocity
  REAL(DP), PARAMETER     :: dw = pi / Tmax            ! delta omega, frequency step size
  ! atomic properties
  REAL(DP), PARAMETER               :: NDensity = 2.7E10_dp ! atomic density
  REAL(DP), PARAMETER, DIMENSION(7) :: EnergyLevels = [ 6.19921, 8.43651, 9.68565, 9.81955, 9.93349, 10.401, 11.0151] ! [eV]
  REAL(DP), PARAMETER, DIMENSION(7) :: DecayRates = [ 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp ] ! decay gammas
  REAL(DP), PARAMETER, DIMENSION(7) :: FrequencyLevels = EnergyLevels / hbar   ! energy frequency levels
  REAL(DP), PARAMETER, DIMENSION(7) :: Delta = FrequencyLevels  ! frequency shifts

  ! coordinate grids
  REAL(DP), DIMENSION(Nt) :: t(Nt), freq(Nt) ! time, space and frequency domains
  REAL(DP), DIMENSION(Nz) :: z(Nz)           ! z (propagation direction) domain
  REAL(DP), DIMENSION(Nr) :: r(Nr)           ! radial domain
  REAL(DP), DIMENSION(7,7), PARAMETER :: u = RESHAPE( &  !u dipole moment [cm] only u_ij where i<j. symmetry is solved in equations
    [ 0.0_dp, 3.4E-9_dp, 0.0_dp    , 0.0_dp   , 0.0_dp   , -2.6E-9_dp , 0.0_dp     , &   ! row 1
      0.0_dp, 0.0_dp   , 1.35E-8_dp, 8.0E-9_dp, 9.9E-9_dp, 0.0_dp     , -4.1E-8_dp , &   ! row 2
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -9.2E-9_dp , 0.0_dp     , &   ! row 3
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -4.7E-9_dp , 0.0_dp     , &   ! row 4
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , -2.8E-8_dp , 0._dp      , &   ! row 5
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , 0.0_dp     , -1.85E-8_dp, &   ! row 6
      0.0_dp, 0.0_dp   , 0.0_dp    , 0.0_dp   , 0.0_dp   , 0.0_dp     , 0.0_dp     ], &  ! row 7
       shape=[7, 7])

! physical variables (field, envelope, density matrix elements, etc)
  COMPLEX(CDP), DIMENSION(Nt, Nr) :: Afield, Efield, RFarr, NLpol, dAfield
  !COMPLEX(CDP), DIMENSION(Nt, Nr, 7, 7) :: rho, RabiFreq
  COMPLEX(CDP), DIMENSION(Nt) :: initial_a_t
  COMPLEX(CDP), DIMENSION(Nt) :: tdata, diff, sol
end module declare
!_____________________________________________________________________
MODULE globals
  USE nrtype
  IMPLICIT NONE

  CONTAINS

    FUNCTION complex_to_real(arr, m)
      IMPLICIT NONE
      INTEGER(I4B), INTENT(IN) :: m
      COMPLEX(DP), INTENT(IN) :: arr(m)
      REAL(DP), DIMENSION(2 * m) :: complex_to_real(2 * m)
      INTEGER i

      DO i=1, m
        complex_to_real(2 * i - 1)=real(arr(i))
        complex_to_real(2 * i)=aimag(arr(i))
      END DO

      RETURN
    END FUNCTION complex_to_real
    !_______________________________________________________________
    FUNCTION real_to_complex(arr, m)
      IMPLICIT NONE
      INTEGER(I4B), INTENT(IN) :: m
      REAL(DP), iNTENT(IN) :: arr(2 * m)
      COMPLEX(CDP) real_to_complex(m)
      INTEGER ::  i

      DO i=1, m
        real_to_complex(i) = cmplx(arr(2 * i - 1), arr(2 * i), kind=CDP)
      END DO

      RETURN
    end function real_to_complex
    !___________________________________________________________________
    FUNCTION sech(x)
      IMPLICIT NONE
      REAL(DP), intent(in) :: x
      REAL(DP) :: sech

      sech=1. / cosh(x)

      RETURN
    END FUNCTION sech
    !___________________________________________________________________
    FUNCTION soliton(x, y) !x for space, y for time
      IMPLICIT NONE
      USE declare, ONLY : im
      COMPLEX(CDP) :: soliton
      REAL(DP), intent(in) :: x, y

      soliton = 4 * ( cosh(3 * y) + 3 * exp(im * 4 * x) * cosh(y) ) * exp(im * x) &
                              / ( cosh(4 * y) + 4 * cosh(2 * y) + 3 * cos(4 * x) )

      RETURN
    END FUNCTION soliton
    !___________________________________________________________________
    FUNCTION lin_sus(w)
      IMPLICIT NONE
      USE declare
      REAL(DP), INTENT(IN) :: w
      REAL(DP) :: lin_sus
      lin_sus = ( NDensity / ( eps0 * hbar ) )
      RETURN
    END FUNCTION lin_sus
  END MODULE globals
!_____________________________________________________________________
SUBROUTINE setup !sets range arrays for time and space coordinates
  IMPLICIT NONE
  USE nrtype
  USE globals
  USE declare, ONLY : Nt, Nr, Nz, t, r, z, freq, Tmax, Lmax, Rmax, pi2, Afield, initial_a_t, beam_waist
  IMPLICIT NONE
  INTEGER(I4B) :: i

  DO i=1, Nt/2
    freq(i) = pi2 * dble(i-1) / (2 * Tmax) !positive and zero frequencies
    freq(i + Nt / 2)=pi2 * dble(i-1-Nt/2)/(2*Tmax) !neg and max/min feq
  END DO

  DO i=1, Nt
    t(i) = (i-1) * 2.0_dp * Tmax / Nt - Tmax

    initial_a_t(i) = cmplx(sech(t(i)), 0.0_dp, kind=CDP)

  END DO

  DO i=1, Nz
    z(i) = (i - 1) * Lmax / Nz
  END DO

  !DO i=1, Nr
  !  r(i) = (i - 1) * Rmax / Nr
  !END DO

  DO i=1, Nr
    r(i) = (i - 1) * Rmax / Nr
    Afield(:, i) = initial_a_t(:) * exp( -( r(i) ** 2 / beam_waist ) )
  END DO



END SUBROUTINE setup
!___________________________________________
PROGRAM main
  IMPLICIT NONE
  USE nrtype
  USE declare
  IMPLICIT NONE
  INTEGER :: i, j, count

  CALL setup
  WRITE(*, *) 'Finished setup'
  WRITE(*, *) 'Starting density matrix evolution'
  open(10, file='initial_pulse.dat')

 DO i=1, Nr !main loop over z
    count=MOD(i, 50)
    IF (count == 0) then !loop for writing into file
      DO j=1, Nt !cycle per time value
        WRITE(10, *) r(i), t(j), REAL(Afield(j, i)) ** 2
      END DO

      WRITE(10, *)
      WRITE(10, *)
    END IF
  END DO
  
  CLOSE(10)
  WRITE(*, *) 'Finished writing data'
  call sleep(1)
END PROGRAM main
!________________________________________________________________________
SUBROUTINE RabiFreq(a, rabi) !Takes electric field to yield rabbi frequencies array per spatial point
  IMPLICIT NONE
  use nrtype
  use declare, only : Nt, Nr, hbar, u
  IMPLICIT NONE
  COMPLEX(CDP), DIMENSION(Nr, Nt), INTENT(IN) :: a ! slowly varying complex amplitude
  COMPLEX(CDP), ALLOCATABLE, INTENT(OUT) :: rabi(:, :, :, :) ! rabbi frequencies array

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

END SUBROUTINE RabiFreq
!________________________________________________________________________
SUBROUTINE matrixDensityEvol(a, r) !solves for rho in time
  IMPLICIT NONE
  USE nrtype
  USE declare, ONLY : dt, decay, omega0, im, Nr, Nt, Delta
  IMPLICIT NONE
  COMPLEX(CDP), DIMENSION(Nr, Nt),INTENT(IN) :: a
  COMPLEX(CDP), DIMENSION(Nr, Nt, 7, 7), INTENT(OUT) :: r ! r for rho, density matrix
  ! COMPLEX(CDP), DIMENSION(Nr, Nt, 7, 7) :: rabi ! rabi values
  COMPLEX(CDP), ALLOCATABLE :: rabi(Nr, Nt, 7, 7)
  INTEGER :: i
  REAL(CDP) :: Tr12, Tr16, Tr23, Tr24, Tr25, Tr27, Tr36, Tr46, Tr56, Tr67

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
  !CALL RabiFreq(a, rabi)


  r(:, 1, 1, 1) = 1.0_dp
  r(:, 1, 1, 2) = 0.0_dp
  r(:, 1, 1, 6) = 0.0_dp
  r(:, 1, 2, 3) = 0.0_dp
  r(:, 1, 2, 4) = 0.0_dp
  r(:, 1, 2, 5) = 0.0_dp
  r(:, 1, 2, 7) = 0.0_dp
  r(:, 1, 3, 6) = 0.0_dp
  r(:, 1, 4, 6) = 0.0_dp
  r(:, 1, 5, 6) = 0.0_dp
  r(:, 1, 6, 7) = 0.0_dp
 
  DO i=1, Nt - 1
    ! diagonal elements
    r(:, i + 1, 1, 1) = r(:, i, 1, 1) + ( decay * ( r(:, i, 2, 2) + r(:, i, 6, 6) ) &
                                          + 2.0_dp * ( CONJG(rabi(:, i, 1, 2)) * AIMAG(r(:, i, 1, 2)) &
                                                     + CONJG(rabi(:, i, 1, 6)) * AIMAG(r(:, i, 1, 6)) ) ) * dt

    r(:, i + 1, 2, 2) = r(:, i, 2, 2) + ( decay * ( -r(:, i, 2, 2) + r(:, i, 3, 3) + r(:, i, 4, 4) &
                                                   + r(:, i, 5, 5) + r(:, i, 7, 7) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 1, 2)) * AIMAG(r(:, i, 1, 2)) &
                                                     - CONJG(rabi(:, i, 2, 3)) * AIMAG(r(:, i, 2, 3)) &
                                                     - CONJG(rabi(:, i, 2, 4)) * AIMAG(r(:, i, 2, 4)) &
                                                     - CONJG(rabi(:, i, 2, 5)) * AIMAG(r(:, i, 2, 5)) &
                                                     - CONJG(rabi(:, i, 2, 7)) * AIMAG(r(:, i, 2, 7)) ) ) * dt

    r(:, i + 1, 3, 3) = r(:, i, 3, 3) + ( decay * ( -r(:, i, 3, 3) + r(:, i, 6, 6) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 2, 3)) * AIMAG(r(:, i, 2, 3)) &
                                                     - CONJG(rabi(:, i, 3, 6)) * AIMAG(r(:, i, 3, 6)) ) ) * dt

    r(:, i + 1, 4, 4) = r(:, i, 4, 4) + ( decay * ( -r(:, i, 4, 4) + r(:, i, 6, 6) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 2, 4)) * AIMAG(r(:, i, 2, 4)) &
                                                     - CONJG(rabi(:, i, 4, 6)) * AIMAG(r(:, i, 4, 6)) ) )  * dt

    r(:, i + 1, 5, 5) = r(:, i, 5, 5) + ( decay * ( -r(:, i, 5, 5) + r(:, i, 6, 6) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 2, 5)) * AIMAG(r(:, i, 2, 5)) &
                                                     - CONJG(rabi(:, i, 5, 6)) * AIMAG(r(:, i, 5, 6))) ) * dt

    r(:, i + 1, 6, 6) = r(:, i, 6, 6) + ( decay * (-4.0 * r(:, i, 6, 6) + r(:, i, 7, 7) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 1, 6)) * AIMAG(r(:, i, 1, 6)) &
                                                     + CONJG(rabi(:, i, 3, 6)) * AIMAG(r(:, i, 3, 6)) &
                                                     + CONJG(rabi(:, i, 4, 6)) * AIMAG(r(:, i, 4, 6)) &
                                                     + CONJG(rabi(:, i, 5, 6)) * AIMAG(r(:, i, 5, 6)) &
                                                     - CONJG(rabi(:, i, 6, 7)) * AIMAG(r(:, i, 6, 7)) ) ) * dt

    r(:, i + 1, 7, 7) = r(:, i, 7, 7) + ( decay * ( -2.0 * r(:, i, 7, 7) ) &
                                          - 2.0_dp * ( CONJG(rabi(:, i, 2, 7)) * AIMAG(r(:, i, 2, 7)) &
                                                     + CONJG(rabi(:, i, 6, 7)) * AIMAG(r(:, i, 6, 7)) ) ) * dt
    ! non diagonal elements
    r(:, i + 1, 1, 2) = r(:, i, 1, 2) + im * ( rabi(:, i, 1, 2) * ( -r(:, i, 1, 1) + r(:, i, 2, 2) ) &
                                                + r(:, i, 1, 2) * ( Delta(2) + omega0 - 0.5_dp * (Tr12) ) ) * dt

    r(:, i + 1, 1, 6) = r(:, i, 1, 6) + im * ( rabi(:, i, 1, 6) * ( -r(:, i, 1, 1) + r(:, i, 6, 6) ) &
                                                + r(:, i, 1, 6) * ( Delta(6) + omega0 - 0.5_dp * (Tr16 + Tr36 + Tr46 + Tr56)) ) * dt

    r(:, i + 1, 2, 3) = r(:, i, 2, 3) + im * ( rabi(:, i, 2, 3) * ( -r(:, i, 2, 2) + r(:, i, 3, 3) ) &
                                                + r(:, i, 2, 3) * ( Delta(2) + Delta(3) + omega0 - 0.5_dp * (Tr12 + Tr23)) ) * dt

    r(:, i + 1, 2, 4) = r(:, i, 2, 4) + im * ( rabi(:, i, 2, 4) * ( -r(:, i, 2, 2) + r(:, i, 4, 4) ) &
                                                + r(:, i, 2, 4) * ( -Delta(2) + Delta(4) + omega0 - 0.5_dp * (Tr12 + Tr24)) ) * dt

    r(:, i + 1, 2, 5) = r(:, i, 2, 5) + im * ( rabi(:, i, 2, 5) * ( -r(:, i, 2, 2) + r(:, i, 5, 5) ) &
                                                + r(:, i, 2, 5) * ( -Delta(2) + Delta(5) + omega0 - 0.5_dp * (Tr12 + Tr25)) ) * dt

    r(:, i + 1, 2, 7) = r(:, i, 2, 7) + im * ( rabi(:, i, 2, 7) * ( -r(:, i, 2, 2) + r(:, i, 7, 7) ) &
                                                + r(:, i, 2, 7) * ( -Delta(2) + Delta(7) + omega0 &
                                                                    - 0.5_dp * (Tr12 + Tr27 + Tr67)) ) * dt

    r(:, i + 1, 3, 6) = r(:, i, 3, 6) + im * ( rabi(:, i, 3, 6) * ( -r(:, i, 3, 3) + r(:, i, 6, 6) ) &
                                                + r(:, i, 3, 6) * ( -Delta(3) + Delta(6) + omega0 &
                                                                    - 0.5_dp * (Tr16 + Tr23 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(:, i + 1, 4, 6) = r(:, i, 4, 6) + im * ( rabi(:, i, 4, 6) * ( -r(:, i, 4, 4) + r(:, i, 6, 6) ) &
                                                + r(:, i, 4, 6) * ( -Delta(4) + Delta(6) + omega0 &
                                                                    - 0.5_dp * (Tr16 + Tr24 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(:, i + 1, 5, 6) = r(:, i, 5, 6) + im * ( rabi(:, i, 5, 6) * ( -r(:, i, 5, 5) + r(:, i, 6, 6) ) &
                                                + r(:, i, 5, 6) * ( -Delta(5) + Delta(6) + omega0 &
                                                                    - 0.5_dp * (Tr16 + Tr25 + Tr36 + Tr46 + Tr56) ) ) * dt

    r(:, i + 1, 6, 7) = r(:, i, 6, 7) + im * ( rabi(:, i, 6, 7) * ( -r(:, i, 6, 6) + r(:, i, 7, 7)) &
                                                + r(:, i, 6, 7) * ( -Delta(6) + Delta(7) + omega0 &
                                                                    - 0.5_dp * (Tr16 + Tr27 + Tr36 + Tr46 + Tr56 + Tr67)) ) * dt
  END DO

END SUBROUTINE matrixDensityEvol
!________________________________________________________________________
SUBROUTINE four1(data,nn,isign)
  IMPLICIT NONE
  INTEGER isign,nn
  double precision data(2*nn)
  INTEGER i,istep,j,m,mmax,n
  double precision tempi,tempr
  DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
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
END SUBROUTINE four1
!_________________________________________

SUBROUTINE TotalPolarization(rho, pol)
  IMPLICIT NONE
  USE nrtype
  USE declare, ONLY : NDensity, Nt, Nr, im, omega0, u
  IMPLICIT NONE
  COMPLEX(CDP), DIMENSION(Nr, Nt, 7, 7), INTENT(IN) :: rho
  REAL(CDP), DIMENSION(Nr, Nt), INTENT(OUT) :: pol

  pol(:, :) = 2 * NDensity * ( u(1, 2) * REAL(rho(:, :, 1, 2), kind=DP) + u(1, 6) * REAL(rho(:, :, 1, 6), kind=DP) &
                             + u(2, 3) * REAL(rho(:, :, 2, 3), kind=DP) + u(2, 4) * REAL(rho(:, :, 2, 4), kind=DP) &
                             + u(2, 5) * REAL(rho(:, :, 2, 5), kind=DP) + u(2, 7) * REAL(rho(:, :, 2, 7), kind=DP) &
                             + u(3, 6) * REAL(rho(:, :, 3, 6), kind=DP) + u(4, 6) * REAL(rho(:, :, 4, 6), kind=DP) &
                             + u(5, 6) * REAL(rho(:, :, 5, 6), kind=DP) + u(6, 7) * REAL(rho(:, :, 6, 7), kind=DP) )

END SUBROUTINE TotalPolarization
!__________________________________
SUBROUTINE polarization_term(p, pterm)
  IMPLICIT NONE
  USE nrtype
  use globals, only : complex_to_real, real_to_complex
  USE declare, ONLY : Nt, Nr, im, omega0, k0, mu0, freq
  COMPLEX(CDP), DIMENSION(Nr, Nt), INTENT(IN) :: p
  COMPLEX(CDP), DIMENSION(Nr, Nt), INTENT(OUT) :: pterm
  COMPLEX(CDP), DIMENSION(Nr, Nt) :: p_fourier_c
  REAL(CDP), DIMENSION(Nr, 2 * Nt) :: p_fourier_r
  INTEGER :: i
  
  DO i=1, Nr 
    p_fourier_r(i, :) = complex_to_real(p(i, :), Nt)
  END DO

  DO i=1, Nr
    CALL four1(p_fourier_r(i, :), Nt, 1)
  END DO

  DO i=1, Nr
  p_fourier_c(i, :) = real_to_complex(p_fourier_r(i, :), Nt)
  END DO

  DO i=1, Nr
    pterm(i, :) = ( im * mu0 / (2 * k0) ) * ( omega0 * ( omega0 + 2 * freq(:) ) - freq(:) ** 2 ) * p_fourier_c(i, :)
  END DO

  DO i=1, Nr
    p_fourier_r(i, :) = complex_to_real(pterm(i, :), Nt)
  END DO

  DO i=1, Nr
    CALL four1(p_fourier_r(i, :), Nt, -1)
  END DO

  DO i=1, Nr
    pterm(i, :) = real_to_complex(p_fourier_r(i, :), Nt)
  END DO

END SUBROUTINE polarization_term
!________________________________________________________________________
SUBROUTINE ddaperp(arr, darr) !radial part of laplacian for field
IMPLICIT NONE
USE nrtype
USE declare, ONLY : Nt, r, Nr, dr
COMPLEX(CDP), DIMENSION(Nt, Nr), INTENT(IN) :: arr
COMPLEX(CDP), DIMENSION(Nt, Nr), INTENT(OUT) :: darr
INTEGER :: i

darr(:, 1) = 4.0 * ( arr(:, 2) - arr(:, 1) ) / (dr ** 2)

DO i = 2, Nr - 1
  darr(:, i) = ( arr(:, i + 1) - 2.0 * arr(:, i) + arr(:, i - 1) ) / (dr ** 2) &
             + ( arr(:, i + 1) - arr(:, i - 1) ) / (2.0_dp * dr *r(i))
END DO

darr(:, Nr) = ( - 2.0 * arr(:, Nr) + arr(:, Nr - 1) ) / (dr ** 2) &
             + (  - arr(:, Nr - 1) ) / (2.0_dp * dr *r(Nr))

END SUBROUTINE ddaperp
!_________________________________________________________________
SUBROUTINE radial_laplacian(arr, darr) !radial part of laplacian for field
IMPLICIT NONE
USE nrtype
USE declare, ONLY : r, Nr, dr
COMPLEX(CDP), DIMENSION(Nr), INTENT(IN) :: arr
COMPLEX(CDP), DIMENSION(Nr), INTENT(OUT) :: darr
INTEGER :: i

darr(1) = 4.0 * ( arr(2) - arr(1) ) / (dr ** 2)

DO i = 2, Nr - 1
  darr(i) = ( arr(i + 1) - 2.0 * arr(i) + arr(i - 1) ) / (dr ** 2) &
             + ( arr(i + 1) - arr(i - 1) ) / (2.0_dp * dr *r(i))
END DO

darr(Nr) = ( - 2.0 * arr(Nr) + arr(Nr - 1) ) / (dr ** 2) &
             + (  - arr(Nr - 1) ) / (2.0_dp * dr *r(Nr))

END SUBROUTINE radial_laplacian
!_________________________________________________________________
SUBROUTINE RightHandFieldTerms(FAfield, omega, RHFT)
  IMPLICIT NONE
  USE nrtype
  USE declare, ONLY : im, Nt, Nr, k0, omega0, dz, v_g
  COMPLEX(CDP), DIMENSION(Nt, Nr) :: perpFA ! radial laplacian of FAfield
  COMPLEX(CDP), INTENT(IN) :: FAfield(Nt, Nr) ! Fourier transform of Afield
  REAL(CDP), INTENT(IN) :: omega
  COMPLEX(CDP), INTENT(OUT) :: RHFT(Nt, Nr) ! Right Hand Field Terms of propagation equation
  
  call ddaperp(FAfield, perpFA)
  
  RHFT(:, :) = k0 * FAfield(:, :) ! (i)
  RHFT(:, :) = RHFT(:, :) - (2 * k0 / v_g) * omega * FAfield(:, :) ! (iii)
  RHFT(:, :) = RHFT(:, :) + (1 / v_g ** 2) * (omega ** 2) * FAfield(:, :) ! (v)
  RHFT(:, :) = RHFT(:, :) - perpFA(:, :) ! (vi)
  RHFT(:, :) = RHFT(:, :)  + ((omega0 / c) ** 2) * FAfield(:, :) ! (vii)
  RHFT(:, :) = RHFT(:, :)  - (2 * omega0 / (c ** 2)) * omega * FAfield(:, :) ! (viii)
  RHFT(:, :) = RHFT(:, :) - omega ** 2 * FAfield(:, :) ! (ix)
END SUBROUTINE RightHandFieldTerms
!_________________________________________________________________
SUBROUTINE PolarizationTerms(PField, omega, RHPT)
  IMPLICIT NONE
  USE nrtype
  USE declare, ONLY : im, Nt, Nr, hbar, u, k0, mu0, omega0, dz, v_g
  COMPLEX(CDP), INTENT(IN) :: FPfield(Nt, Nr) ! Fourier transform of Polarization field (PField)
  REAL(CDP), INTENT(IN) :: omega
  COMPLEX(CDP), INTENT(OUT) :: RHPT(Nt, Nr) ! Right Hand Polarization Terms of propagation equation

  RHPT(:, :) = im * ((omega0 ** 2 * mu0) / (2*k0)) * FPfield(:, :) ! (i)
  RHPT(:, :) = RHPT(:, :) + im * (omega0 * mu0 / k0) * omega * FPfield(:, :) ! (ii)
  RHPT(:, :) = RHPT(:, :) - im * (mu0 / (2 * k0)) * (omega ** 2) * FPfield(:, :) ! (iii)
END SUBROUTINE PolarizationTerms