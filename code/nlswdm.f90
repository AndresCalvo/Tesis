!_________________________________________________________________
MODULE nrtype
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(1)
    INTEGER, PARAMETER :: DP = KIND(1.0D0)
  END MODULE nrtype
!_________________________________________________________________ 
MODULE global
  USE nrtype
  
  CONTAINS
  FUNCTION sech(x)
    use nrtype
    implicit none
    REAL(DP) :: sech
    REAL(DP), intent(in) :: x
    sech=1./cosh(x)
    return
  end Function sech
  FUNCTION soliton(x, y) !x for space, y for time
    use nrtype
    use declare, only : im
    implicit none
    COMPLEX(DP) :: soliton
    REAL(DP), intent(in) :: x,y
    soliton = 4 * (cosh(3*y)+3*exp(4*im*x)*cosh(y)) * exp(im*x) / (cosh(4*y)+4*cosh(2*y)+3*cos(4*x))
    return
  end Function soliton
END MODULE global
!_________________________________________________________________
module declare

  use nrtype
  implicit none

  character(7), parameter :: filename='sol.dat' !save calculated solution
  character(8), parameter :: filename1='tsol.dat' !saves theoretical sol. for 2nd order soliton
  character(8), parameter :: filename2='diff.dat' !saves difference between the above solutions

  REAL(DP), parameter ::  pi2 = 6.28318530717959
  REAL(DP), parameter :: beta2 = -1.0d0 !second coeff in beta exp
  REAL(DP), parameter :: beta3 = 0.0 !secsplot ond coeff in beta exp
  REAL(DP), parameter :: alpha = 0.0 !absorption coeff
  REAL(DP), parameter :: gamma = 1.0 !nonlinear coeff

  REAL(DP), parameter :: c = 299792458.0d0
  REAL(DP), parameter :: ng = 1.465
  REAL(DP), parameter :: vg = c / ng

  REAL(DP), parameter :: L = 5.0d0 !total fiber length
  REAL(DP), parameter :: Tmax = 10.0d0 !max time in frame
  INTEGER(I4B), parameter :: N = 4096 !size of T window
  REAL(DP), parameter :: dz = 0.001 !size of z step
  INTEGER(I4B), parameter :: zsteps = floor(L / dz) !steps in z

  REAL(DP), parameter :: dt = Tmax / N !time step size
  INTEGER(I4B), parameter :: Ord = 3
  INTEGER(I4B)  j, k
  integer count
  REAL(DP) t(N), z(zsteps), freq(N) !time, space and frequency domains

  double complex, parameter :: im = dcmplx(0.D0,1.D0)
  double complex a(N), tdata(N), diff(N), sol
end module declare
!_____________________________________________________________________
!module oldmain
!  !former main progran=m
!  use declare
!  implicit none
!
!  INTERFACE
!
!    FUNCTION sech(x)
!      double precision :: sech
!      double precision, intent(in) :: x
!    end Function sech
!    
!    FUNCTION soliton(x, y)
!      double precision, INTENT(IN) :: x, y
!      double complex soliton
!  END INTERFACE
!
! 
!  
!  open(10, file=filename)
!  open(11, file=filename1)
!  open(12, file=filename2)
!  do j=1, zsteps !main loop over z
!    count=MOD(j, j)
!    if (count == 0) then !loop for writing into file
!      do k=1, N !cycle per time value
!        sol=soliton(z(j), t(k)) !temp var to avoid constantly fetching value in array and perform calculations ??
!        write(10, *) t(k), z(j), abs(a(k))
!        write(11, *) t(k), z(j), abs(sol)
!        write(12, *) t(k), z(j), abs(sol-a(k))
!      end do
!      write(10, *)
!      write(10, *)
!      write(11, *)
!      write(11, *)
!      write(12, *)
!      write(12, *)
!    end if
!
!    call SSFstep(a, dz, 3)
!
!  enddo
!
!  close(10)
!  close(11)
!  close(12)
!end module oldmain
!________________________________________________________________________
program main
  use declare
  use nrtype
  IMPLICIT NONE
  INTEGER :: i
  double precision h
  COMPLEX(DP) :: r11, r22, r33, r44, r55, r66, r77 !r for rho, density matrix elements
  COMPLEX(DP) :: r12, r16, r23, r24, r25, r27, r36, r46, r56, r67 !non diagonal rho elements
  COMPLEX(DP) :: delta1, delta2, delta3, delta4, delta5, delta6, delta7
  COMPLEX(DP) :: RF12, RF16, RF23, RF24, RF25, RF27, RF36, RF46, RF56, RF67 !RF for Rabi frequency
  REAL(DP), DIMENSION(N) :: u12, u16, u23, u24, u25, u27, u36, u46, u56, u67 !u for mu, dipole moment elements. Ordered elements

  call setup !sets time, space and frequency ranges
  u12(1) = 3.4D-9
  u16(1) = -2.6D-9
  u23(1) = 9.9D-9
  u24(1) = 8D-9
  u25(1) = 1.35D-8
  u27(1) = -4.1D-8
  u36(1) = -9.2D-9
  u46(1) = -4.7D-9
  u56(1) = -2.8D-8
  u67(1) = -1.85D-8


END

!________________________________________________________________________
SUBROUTINE rhodevs(rho, drhodt, Afield)
  USE nrtype
  USE global
  IMPLICIT NONE
  COMPLEX(DP) :: rho, drhodt, Afield
END SUBROUTINE rhodevs
!________________________________________________________________________
SUBROUTINE four1(data,nn,isign)
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
SUBROUTINE setup !sets range arrays for time and space coordinates
  use declare
  REAL(DP):: sech
  REAL(DP), DIMENSION(N) :: t, a
  do j=1, n/2
    freq(j)=pi2*dble(j-1)/(2*Tmax) !positive and zero frequencies
    freq(j+N/2)=pi2*dble(j-1-N/2)/(2*Tmax) !neg and max/min feq
  end do
  do j=1, N
    t(j) = (j-1) * 2*Tmax / N - Tmax
    a(j) = Ord * dcmplx(sech(t(j)))
  end do
  do j=1, zsteps
    z(j) = (j-1) * L / zsteps
  end do
END SUBROUTINE setup
!___________________________________________
SUBROUTINE dispersionExp(arr, h) !dispersion operator in fourier space, h is step length
  use declare
  IMPLICIT NONE
  double complex arr(N)
  double precision h, tarr(2*N)

  INTERFACE

    FUNCTION dr_to_c(ar, m)
      integer, intent(in) :: m
      double precision, iNTENT(IN) :: ar(2*m)
      double complex  dr_to_c(m)
    END FUNCTION dr_to_c

    FUNCTION c_to_dr(ar, m)
      integer, intent(in) :: m
      double complex, intent(in) :: ar(m)
      double precision c_to_dr(2*m)
    END FUNCTION c_to_dr

  END INTERFACE

  tarr = c_to_dr(arr, N)
  call four2(tarr, N, 1)
  arr = dr_to_c(tarr, N)
  arr = exp(h*(im*beta2/2*(freq**2)+im*beta3/6*(freq**3)-alpha/2))*arr
  tarr = c_to_dr(arr, N)
  call four2(tarr, N, -1)
  tarr = tarr / N
  arr = dr_to_c(tarr, N)

END SUBROUTINE dispersionExp
!__________________________________
SUBROUTINE SSFstep(arr, h, cl) !cl is for total aproximation cycles for A(z+dz), must be int >= 0. 0 means simply using A(z)=A(z+dz) and no iterations.
  use declare
  double complex arr(N), temp(N), temp2(N)
  double precision h
  integer cl, p
  temp = arr !for storing aproximations of A(z+dz)
  temp2 = arr !dummy starting point for each iteration A(z)

  call dispersionExp(temp, h/2)
  temp = exp(h*i*gamma*(abs(arr**2)))*temp
  call dispersionExp(temp, h/2) !here temp is first aprox of A(z+dz)

  if (cl >= 1) then
   do p=1, cl
      temp2 = arr !need aditional copy for operating dispersion and maintain values for A(z) and A(z+dz)
      call dispersionExp(temp2, h/2) !temp2 = arr is A(z)
      temp = exp(h/2*i*gamma*(abs(arr**2)+abs(temp)**2))*temp2 !now we dont need previous values of A(z+dz), so we rewrite the variable
      call dispersionExp(temp, h/2) !here temp is next aprox of A(z+dz)
    end do
  end if
   arr = temp
END SUBROUTINE SSFstep
!_________________________________________________________________
FUNCTION c_to_dr(arr, m)
  use nrtype
  IMPLICIT NONE
  integer, intent(in) :: m
  double complex , intent(IN) :: arr(m)
  double precision c_to_dr(2*m)
  integer l
  do l=1, m
    c_to_dr(2*l-1)=real(arr(l))
    c_to_dr(2*l)=aimag(arr(l))
  end do
  return
end function c_to_dr
!_______________________________________________________________
FUNCTION dr_to_c(arr, m)
  use nrtype
  IMPLICIT NONE
  integer, intent(in) :: m
  double precision, iNTENT(IN) :: arr(2*m)
  double complex dr_to_c(m)
  integer l
  do l=1, m
    dr_to_c(l)=dcmplx(arr(2*l-1), arr(2*l))
  end do
  return
end function dr_to_c
!___________________________________________________________________

!___________________________________________________________________

!___________________________________________________________________
