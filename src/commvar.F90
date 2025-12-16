!+---------------------------------------------------------------------+
!| This module is to define common variables.                          |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 06-02-2021  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module commvar
  !
  use, intrinsic :: iso_c_binding
  use commtype
  !
  !
  implicit none
  !
  integer :: ia,ja,ka,im,jm,km,is,ie,js,je,ks,ke,ndims
  logical :: lihomo,ljhomo,lkhomo
  parameter (lihomo=.true., ljhomo=.true., lkhomo=.true.)
  integer :: npdci,npdcj,npdck
  logical :: lfftk,lreport,ltimrpt
  real(8) :: preptime=0.d0
  integer :: islice,jslice,kslice
  integer :: numq
  character(len=10) :: turbmode
  complex(8) :: imag = CMPLX(0.d0,1.d0,8)
  !
  real(8), allocatable :: u1(:,:,:), u2(:,:,:), u3(:,:,:)
  real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
  !
  !
  integer :: nstep,maxstep,feqwsequ,filenumb, nxtwsequ,feqwspe,nxtwspe, kmax
  real(8) :: nu, ref_tem, reynolds, deltat, time
  logical :: lwsequ,lwspectra
  !
  ! Middle process for ut calculate
  type(C_PTR) :: c_u1spe, c_u2spe, c_u3spe, c_u1x1, c_u1x2, c_u1x3, &
                 c_u2x1, c_u2x2, c_u2x3, c_u3x1, c_u3x2, c_u3x3,    &
                 c_u1xixi, c_u2xixi, c_u3xixi, c_thetaxixi
  complex(c_double_complex), pointer, dimension(:,:,:):: &
                u1spe, u2spe, u3spe, u1x1, u1x2, u1x3, u2x1, u2x2, u2x3,&
                u3x1, u3x2, u3x3, u1xixi, u2xixi, u3xixi, thetaxixi
  real(8) :: eta_min
  !
  ! Middle process for RK3
  type(C_PTR) :: c_u1tA, c_u2tA, c_u3tA, c_u1tB, c_u2tB, c_u3tB, c_u1tC, c_u2tC, &
                 c_u3tC, c_force1, c_force2, c_force3
  complex(c_double_complex), pointer, dimension(:,:,:):: &
                 u1tA, u2tA, u3tA, u1tB, u2tB, u3tB, u1tC, u2tC, u3tC, &
                 force1, force2, force3
  !
  real(8), allocatable, dimension(:) :: Es,Ed,kn
  integer, allocatable, dimension(:) :: Ecount
  real(8) :: Esspe, Edspe
  integer :: allkmax
  !
  logical :: lprojectd
  real(8) :: target_energy
  integer :: forcemethod, forcek
  !
  integer :: initialmethod
  !
  contains
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+