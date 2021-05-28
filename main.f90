module global_variables
! math constants
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)

! grid system
  integer :: nx
  real(8) :: length_x, dx
  real(8),allocatable :: xx(:)

! quantum system
  real(8) :: R_dist
  real(8),allocatable :: phi_gs(:), rho_gs(:)
  real(8),allocatable :: v_ee(:), v_en(:)
  real(8),allocatable :: ham_mat(:,:), kinetic_op(:,:)
  real(8),allocatable :: hf_orbitals(:,:)
  real(8),allocatable :: ham_sp_mat(:,:)
  real(8),allocatable :: v_hx(:)


end module global_variables
!-------------------------------------------------------------------------------
program main
  implicit none

  call initialize
  call hartree_fock


end program main
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix

  R_dist = 7.0d0
  length_x = 40d0
  dx = 0.15d0
  nx = nint(length_x/dx)

  write(*,*)'length_x =',length_x
  write(*,*)'dx       =',dx
  write(*,*)'nx       =',nx


  allocate(xx(nx), phi_gs(nx), rho_gs(nx), v_hx(nx))
  allocate(v_ee(0:nx), v_en(nx))
  allocate(ham_mat(nx,nx), kinetic_op(nx,nx))
  allocate(hf_orbitals(nx,nx),ham_sp_mat(nx,nx))

  do ix = 1, nx
    xx(ix) = dx*ix-0.5d0*length_x
  end do

! construct interaction
  do ix = 0, nx
    v_ee(ix) = 1d0/sqrt(abs(dx*ix)**2+2d0)
  end do

  do ix = 1, nx
    v_en(ix) = -1d0/sqrt(abs(xx(ix)-R_dist/2d0)**2 + 1d0) &
               -1d0/sqrt(abs(xx(ix)+R_dist/2d0)**2 + 1d0)
  end do

  call construct_kinetic_energy_operator

! set the single-particle hamiltonian
  ham_sp_mat = kinetic_op
  do ix = 1,nx
    ham_sp_mat(ix,ix) = ham_sp_mat(ix,ix) + v_en(ix)
  end do

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine construct_kinetic_energy_operator
  use global_variables
  implicit none
  real(8) :: c0,c1,c2
  real(8) :: mass_e, mass_p
  integer :: i,j

  mass_p = 1836d0
  mass_e = 2d0*mass_p/(2d0*mass_p+1d0)

  c0 = -0.5d0/(mass_e*dx**2)*(-5d0/2d0)
  c1 = -0.5d0/(mass_e*dx**2)*(4d0/3d0)
  c2 = -0.5d0/(mass_e*dx**2)*(-1d0/12d0)

  kinetic_op = 0d0

  do i = 1, nx

    j = i-2
    if(j>= 1 .and. j<= nx)kinetic_op(i,j)=c2

    j = i-1
    if(j>= 1 .and. j<= nx)kinetic_op(i,j)=c1

    j = i
    if(j>= 1 .and. j<= nx)kinetic_op(i,j)=c0

    j = i+1
    if(j>= 1 .and. j<= nx)kinetic_op(i,j)=c1

    j = i+2
    if(j>= 1 .and. j<= nx)kinetic_op(i,j)=c2

  end do

end subroutine construct_kinetic_energy_operator
!-------------------------------------------------------------------------------
subroutine hartree_fock
  use global_variables
  implicit none
  real(8),parameter :: mixing_rate = 0.01d0
  real(8) :: Egs
  integer :: iscf, nscf, ix

  nscf = 1000

! initial guess
  rho_gs = 0d0

  do iscf = 1, nscf

    call calc_hartree_fock_potential
    call calc_hartree_fock_orbitals

! calc one-body density
    rho_gs = mixing_rate*phi_gs**2 + (1d0-mixing_rate)*rho_gs

    call calc_gs_energy(Egs)
    write(*,"(A,2x,I7,e26.16e3)")'iscf, Egs=',iscf, Egs
  end do

  open(20,file='phi_rho_x.out')
  do ix = 1,nx
    write(20,"(999e26.16e3)")xx(ix),phi_gs(ix),rho_gs(ix)
  end do
  close(20)

end subroutine hartree_fock
!-------------------------------------------------------------------------------
! Here, we compute the half of Hartree potential
subroutine calc_hartree_fock_potential
  use global_variables
  implicit none
  real(8) :: ss
  integer :: ix, jx

  v_hx = 0d0

  do ix = 1,nx
    ss = 0d0
    do jx = 1,nx
      ss = ss + rho_gs(jx)*v_ee(abs(ix-jx))
    end do
    ss = ss*dx
    v_hx(ix) = ss
  end do


end subroutine calc_hartree_fock_potential
!-------------------------------------------------------------------------------
subroutine calc_hartree_fock_orbitals
  use global_variables
  implicit none
  integer :: ix
  real(8) :: ss
!==LAPACK
    integer :: nmax
    integer :: lwork
    real(8),allocatable :: work_lp(:)
    real(8),allocatable :: rwork(:),w(:)
    integer :: info

    nmax = nx
    lwork = 6*nmax**2
    allocate(work_lp(lwork),rwork(3*nmax-2),w(nmax))
!==LAPACK



  ham_mat = kinetic_op
  do ix = 1,nx
    ham_mat(ix,ix) = ham_mat(ix,ix) + v_en(ix) + v_hx(ix)
  end do

  hf_orbitals=ham_mat
  call dsyev('V', 'U', nmax, hf_orbitals, nmax, w, work_lp, lwork, info)
  hf_orbitals = hf_orbitals/sqrt(dx)

  phi_gs(:) = hf_orbitals(:,1)
  
  ss = sum(phi_gs(:)**2)*dx
!  write(*,*)'norm-check = ',ss
  phi_gs = phi_gs/sqrt(ss)

end subroutine calc_hartree_fock_orbitals
!-------------------------------------------------------------------------------
subroutine calc_gs_energy(Egs)
  use global_variables
  implicit none
  real(8),intent(out) :: Egs
  real(8) :: hphi(nx), rho_tmp(nx), E_coulomb, ss
  integer :: ix, jx

  Egs = 0d0

  hphi = matmul(ham_sp_mat, phi_gs)
  Egs = 2d0*sum(phi_gs*hphi)*dx


  rho_tmp = phi_gs**2
  E_coulomb = 0d0
  do ix = 1, nx
    ss = 0d0
    do jx = 1,nx
      ss = ss + rho_tmp(jx)*v_ee(abs(ix-jx))
    end do
    ss = ss*dx
    E_coulomb = E_coulomb + ss*rho_tmp(ix)
  end do
  E_coulomb = E_coulomb*dx

  Egs= Egs + E_coulomb + 1d0/R_dist


end subroutine calc_gs_energy
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
