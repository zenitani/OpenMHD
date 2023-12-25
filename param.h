!-*- mode: f90 -*-
  integer, parameter :: version = 20240101
  real(8), parameter :: gamma = 5.d0 / 3.d0
!  real(8), parameter :: gamma = 2.d0
! do not edit
  integer, parameter :: var1 = 9, var2 = 4  ! number of variables
  integer, parameter :: mx = 1, vx = 1  ! momentum (Mx) / velocity (vx)
  integer, parameter :: my = 2, vy = 2  ! momentum (My) / velocity (vy)
  integer, parameter :: mz = 3, vz = 3  ! momentum (Mz) / velocity (vz)
  integer, parameter :: en = 4, pr = 4  ! total energy (E) / pressure (p)
  integer, parameter :: ro = 5  ! plasma density (rho)
  integer, parameter :: bx = 6  ! magnetic field (Bx)
  integer, parameter :: by = 7  ! magnetic field (By)
  integer, parameter :: bz = 8  ! magnetic field (Bz)
  integer, parameter :: ps = 9  ! virtual potential (psi) for div B
! do not edit
! end
