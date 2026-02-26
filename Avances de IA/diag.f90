! Quick diagnostic: print max field value after each major step
! to find exactly where the blowup originates
program diag
  use nrtype
  use declare
  use globals
  implicit none
  integer :: ii
  real(dp), dimension(2*Nt) :: temp1, temp2
  real(dp) :: maxval_field

  call setup
  write(*,'(A,F12.6)') 'Initial max|A|^2: ', maxval(abs(Afield)**2)

  ! Step 1: density matrix
  call matrixDensityEvol(Afield, rho)
  write(*,*) 'After matrixDensityEvol: OK'
  write(*,'(A,E14.6)') '  max|rho(1,1)|   = ', maxval(abs(rho(:,:,1,1)))
  write(*,'(A,E14.6)') '  max|rho(1,2)|   = ', maxval(abs(rho(:,:,1,2)))

  ! Step 2: polarization
  call TotalPolarization(rho, Pfield)
  write(*,*) 'After TotalPolarization: OK'
  write(*,'(A,E14.6)') '  max|Pfield|     = ', maxval(abs(Pfield))

  ! Step 3: FFT
  do ii = 1, Nr
    temp1 = complex_to_real(Pfield(:, ii), Nt)
    call four1(temp1, Nt, 1)
    Pfield_f(:, ii) = real_to_complex(temp1, Nt)

    temp2 = complex_to_real(Afield(:, ii), Nt)
    call four1(temp2, Nt, 1)
    Afield_f(:, ii) = real_to_complex(temp2, Nt)
  end do
  write(*,*) 'After FFT: OK'
  write(*,'(A,E14.6)') '  max|Afield_f|   = ', maxval(abs(Afield_f))
  write(*,'(A,E14.6)') '  max|Pfield_f|   = ', maxval(abs(Pfield_f))

  ! Step 4: propagation
  call PropagationStep(Afield_f, Pfield_f, Afield_f_p1)
  write(*,*) 'After PropagationStep: OK'
  write(*,'(A,E14.6)') '  max|Afield_f_p1|= ', maxval(abs(Afield_f_p1))

  ! Step 5: IFFT
  do ii = 1, Nr
    temp2 = complex_to_real(Afield_f_p1(:, ii), Nt)
    call four1(temp2, Nt, -1)
    Afield(:, ii) = real_to_complex(temp2, Nt) / real(Nt, dp)
  end do
  write(*,*) 'After IFFT: OK'
  write(*,'(A,E14.6)') '  max|A| after z1 = ', maxval(abs(Afield))
  write(*,'(A,E14.6)') '  max|A|^2 af z1  = ', maxval(abs(Afield)**2)

end program diag
