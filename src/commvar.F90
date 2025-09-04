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
  real(8), allocatable :: u1(:,:,:), u2(:,:,:)
  real(8), allocatable, dimension(:,:) :: k1,k2
  !
  !
  integer :: nstep,maxstep,feqwsequ,filenumb, nxtwsequ,feqwspe,nxtwspe, kmax
  real(8) :: nu, ref_tem, reynolds, deltat, time
  logical :: lwsequ
  !
  ! Middle process for ut calculate
  type(C_PTR) :: c_u1spe, c_u2spe, c_u1x1, c_u1x2, c_u2x1, c_u2x2, c_u1xixi, c_u2xixi, c_thetaxixi
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe, u1x1, u1x2, u2x1, u2x2, u1xixi, u2xixi
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: thetaxixi
  real(8) :: eta_min
  !
  ! Middle process for RK3
  type(C_PTR) :: c_u1tA, c_u2tA,c_u1tB, c_u2tB,c_u1tC, c_u2tC
  complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1tA, u2tA, u1tB, u2tB, u1tC, u2tC
  !
  real(8), allocatable, dimension(:) :: Es,Ed,kn
  integer, allocatable, dimension(:) :: Ecount
  real(8) :: Esspe, Edspe
  integer :: allkmax
  !
  logical :: lforce, lprojectd
  real(8) :: target_energy
  !
  contains
  !
end module commvar
!+---------------------------------------------------------------------+
!| The end of the module commvar.                                      |
!+---------------------------------------------------------------------+