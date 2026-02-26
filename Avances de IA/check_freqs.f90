program check_freqs
  use nrtype
  use declare
  use globals
  implicit none
  integer :: i, ii
  real(dp), dimension(2*Nt) :: temp1, temp2
  real(dp) :: denom, ratio, maxratio, maxperp, maxrhpt

  call setup
  write(*,'(A,F12.6)') 'dt        = ', dt
  write(*,'(A,F12.6)') 'dz        = ', dz
  write(*,'(A,F12.6)') 'k0*v_g    = ', k0*v_g
  write(*,'(A,F12.6)') 'freq(1)   = ', freq(1)
  write(*,'(A,F12.6)') 'freq(2)   = ', freq(2)
  write(*,'(A,F12.6)') 'freq(Nt/2)= ', freq(Nt/2)
  write(*,'(A,F12.6)') 'freq(Nt/2+1)= ', freq(Nt/2+1)
  write(*,*)

  ! Check max |vg/(k0*vg - omega)| ratio
  maxratio = 0.0_dp
  do i = 1, Nt
    denom = k0*v_g - freq(i)
    if (abs(denom) < 1.0E-6_dp) cycle
    ratio = abs(v_g / denom)
    if (ratio > maxratio) maxratio = ratio
  end do
  write(*,'(A,E14.6)') 'Max |vg/(k0vg-w)|  = ', maxratio

  ! Print frequencies near omega=1
  write(*,*) 'Frequencies near omega=k0*v_g=1:'
  do i = 1, Nt
    denom = k0*v_g - freq(i)
    if (abs(denom) < 0.2_dp) then
      write(*,'(A,I4,A,F12.6,A,E14.6)') '  i=', i, '  freq=', freq(i), '  1/(1-w)=', 1.0_dp/denom
    end if
  end do

  ! Now check actual field magnitudes after one propagation step
  write(*,*)
  write(*,'(A,E14.6)') 'Initial max|A|^2 = ', maxval(abs(Afield)**2)

  ! density matrix
  call matrixDensityEvol(Afield, rho)
  write(*,'(A,E14.6)') 'max|rho(1,2)|    = ', maxval(abs(rho(:,:,1,2)))

  ! total polarization
  call TotalPolarization(rho, Pfield)
  write(*,'(A,E14.6)') 'max|Pfield|      = ', maxval(abs(Pfield))

  ! FFT
  do ii = 1, Nr
    temp1 = complex_to_real(Pfield(:, ii), Nt)
    call four1(temp1, Nt, 1)
    Pfield_f(:, ii) = real_to_complex(temp1, Nt)

    temp2 = complex_to_real(Afield(:, ii), Nt)
    call four1(temp2, Nt, 1)
    Afield_f(:, ii) = real_to_complex(temp2, Nt)
  end do
  write(*,'(A,E14.6)') 'max|Afield_f|    = ', maxval(abs(Afield_f))
  write(*,'(A,E14.6)') 'max|Pfield_f|    = ', maxval(abs(Pfield_f))

  ! Check perpFA magnitude
  call ddaperp(Afield_f, Afield_f_p1)
  maxperp = maxval(abs(Afield_f_p1))
  write(*,'(A,E14.6)') 'max|perpFA|      = ', maxperp

  ! Check that correction is bounded
  write(*,'(A,E14.6)') 'max correction   = ', maxperp * maxratio * dz / 2.0_dp

  ! Propagation step
  call PropagationStep(Afield_f, Pfield_f, Afield_f_p1)
  write(*,'(A,E14.6)') 'max|A_f_p1|      = ', maxval(abs(Afield_f_p1))

  ! IFFT
  do ii = 1, Nr
    temp2 = complex_to_real(Afield_f_p1(:, ii), Nt)
    call four1(temp2, Nt, -1)
    Afield(:, ii) = real_to_complex(temp2, Nt) / real(Nt, dp)
  end do
  write(*,'(A,E14.6)') 'max|A|^2 after z1= ', maxval(abs(Afield)**2)

end program check_freqs
