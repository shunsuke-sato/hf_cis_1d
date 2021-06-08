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

! symmetrization
  logical,parameter :: if_v_hx_symmetrize = .true.


end module global_variables
!-------------------------------------------------------------------------------
program main
  implicit none

  call initialize
  call hartree_fock
  call tamm_dancoff

end program main
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix
  real(8) :: ss

  R_dist = 2.0d0
  length_x = 150d0
  dx = 0.075d0
  nx = nint(length_x/dx)

  write(*,*)'length_x =',length_x
  write(*,*)'dx       =',dx
  write(*,*)'nx       =',nx


  allocate(xx(nx), phi_gs(nx), rho_gs(nx), v_hx(nx))
  allocate(v_ee(0:nx), v_en(nx))
  allocate(ham_mat(nx,nx), kinetic_op(nx,nx))
  allocate(hf_orbitals(nx,nx),ham_sp_mat(nx,nx))

  do ix = 1, nx
    xx(ix) = dx*ix
  end do
! symmetrize
  ss = sum(xx)/nx
  xx = xx -ss

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
  real(8),parameter :: mixing_rate = 0.75d0
  real(8) :: Egs
  integer :: iscf, nscf, ix

  nscf = 20

! initial guess
  rho_gs = 0d0

  do iscf = 1, nscf

    call calc_hartree_fock_potential
    call calc_hartree_fock_orbitals

! calc one-body density
    
    rho_gs = mixing_rate*phi_gs**2 + (1d0-mixing_rate)*rho_gs
    if(iscf == 1)rho_gs = phi_gs**2

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
  real(8) :: ss, v_hx_tmp(nx)
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

! symemtrize
  if(if_v_hx_symmetrize)then
    v_hx_tmp = v_hx
    do ix = 1,nx
      v_hx(ix) = 0.5d0*(v_hx_tmp(ix) + v_hx_tmp(nx+1-ix))
    end do
  end if

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
subroutine tamm_dancoff
  use global_variables
  implicit none
  real(8) :: ham_cis(nx,nx), hphi(nx),Egs
  real(8) :: rho00(nx), rho_0j(nx), rho_jk(nx), rho_0k(nx)
  real(8) :: vh00(nx), vh_0j(nx)
  real(8) :: h00, hjk
  integer :: j,k, nst_write
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

    ham_cis = 0d0
    call calc_gs_energy(Egs)    
    ham_cis(1,1) = Egs

    hphi(:) = matmul(ham_sp_mat,hf_orbitals(:,1))
    h00 = sum(hf_orbitals(:,1)*hphi(:))*dx

    rho00 = hf_orbitals(:,1)**2
    call calc_pot(rho00, vh00)

    do j = 2, nx
      hphi(:) = matmul(ham_sp_mat,hf_orbitals(:,j))
      rho_0j(:) = hf_orbitals(:,1)*hf_orbitals(:,j)
      call calc_pot(rho00, vh_0j)
      do k = j, nx
        ham_cis(j,k)= sum(hphi(:)*hf_orbitals(:,k))*dx
        if(j == k)ham_cis(j,k)=ham_cis(j,k)+h00+1d0/R_dist

        rho_jk = hf_orbitals(:,j)*hf_orbitals(:,k)
        rho_0k = hf_orbitals(:,1)*hf_orbitals(:,k)
        ham_cis(j,k)= ham_cis(j,k) + sum(vh00*rho_jk)*dx &
                                   + sum(vh_0j*rho_0k)*dx

        ham_cis(k,j)= ham_cis(j,k)
      end do
    end do

    call dsyev('V', 'U', nmax, ham_cis, nmax, w, work_lp, lwork, info)


!    nst_write = min(20,nx)

    open(20,file='data_CIS.out')
    write(20,"(999e26.16e3)")R_dist, w(1:min(20,nx))
    close(20)


    contains
      subroutine calc_pot(rho_in, vpot_out)
        implicit none
        real(8),intent(in) :: rho_in(nx)
        real(8),intent(out) :: vpot_out(nx)
        integer :: ix, jx
        real(8) :: ss
        
        do ix = 1,nx
          ss = 0d0
          do jx = 1,nx
            ss = ss + rho_in(jx)*v_ee(abs(ix-jx))
          end do
          ss = ss*dx
          vpot_out(ix) = ss
        end do

      end subroutine calc_pot
end subroutine tamm_dancoff
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
