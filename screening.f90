program screening
  use spglib_f08

  implicit none

  include "fftw3.f"

  integer, parameter :: dp=selected_real_kind(15,300)
  integer, parameter :: long=selected_int_kind(18)

  real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: Ha2eV=27.21138342902473_dp
  real(dp), parameter :: bohr2A=0.52917720859_dp
  real(dp), parameter :: alpha_xc=0.0246867_dp
  real(dp), parameter :: eps_sigma=2.4_dp
  real(dp), parameter :: d_sigma=2.8_dp/bohr2A

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  character(10) :: par_n_char, par_k_char
  character(11) :: symbol_int
  character(50) :: file_name

  logical :: par_NL
  logical :: par_inter
  logical :: par_sigma
  logical :: par_xc
  logical :: calc_k_symm
  logical :: calc_at_symm
  logical :: conv
  logical :: found
  logical, allocatable :: write_rot(:)

  integer :: i
  integer :: j
  integer :: k
  integer :: a
  integer :: b
  integer :: c
  integer :: N
  integer :: num_imp
  integer :: iter
  integer :: max_iter
  integer :: par_n
  integer :: par_f
  integer :: size_H
  integer :: num_weights
  integer :: num_kpoints
  integer :: mesh(3)
  integer :: is_shift(3)
  integer :: par_k
  integer :: max_size
  integer :: num_symops
  integer :: weight
  integer, allocatable :: iweights(:)
  integer, allocatable :: grid_point(:,:)
  integer, allocatable :: map(:)
  integer, allocatable :: atom_types(:)
  integer, allocatable :: rotation(:,:,:)
  integer, allocatable :: at_t(:)
  integer(long) :: plan_r2c
  integer(long) :: plan_c2r

  real(dp) :: lat_par
  real(dp) :: par_l
  real(dp) :: cell(2,2)
  real(dp) :: scell(2,2)
  real(dp) :: scellr(2,2)
  real(dp) :: shift(2)
  real(dp) :: eps
  real(dp) :: mu
  real(dp) :: vF
  real(dp) :: area
  real(dp) :: arear
  real(dp) :: amix
  real(dp) :: d
  real(dp) :: rms
  real(dp) :: charge2D
  real(dp) :: av_pot
  real(dp) :: kF
  real(dp) :: x
  real(dp) :: rms_tol
  real(dp) :: at_pos(2,2)
  real(dp) :: d_min
  real(dp) :: lattice(3,3)
  real(dp) :: at_posf(2,2)
  real(dp) :: shiftf(2)
  real(dp) :: new_pos(2)
  real(dp) :: k_point(2)
  real(dp) :: k_point2(2)
  real(dp), allocatable :: at_poss(:,:)
  real(dp), allocatable :: at_possf(:,:)
  real(dp), allocatable :: imp_pos(:,:)
  real(dp), allocatable :: par_h(:)
  real(dp), allocatable :: Z(:)
  real(dp), allocatable :: positions(:,:)
  real(dp), allocatable :: rweights(:,:)
  real(dp), allocatable :: translation(:,:)
  real(dp), allocatable :: qpoints(:,:)
  real(dp), allocatable :: r(:,:,:)
  real(dp), allocatable :: rr(:,:,:)
  real(dp), allocatable :: rrp(:,:,:)
  real(dp), allocatable :: rd(:,:)
  real(dp), allocatable :: rdr(:,:)
  real(dp), allocatable :: V(:,:)
  real(dp), allocatable :: Vmix(:,:)
  real(dp), allocatable :: V0(:,:)
  real(dp), allocatable :: rho(:,:)
  real(dp), allocatable :: Vind(:,:)
  real(dp), allocatable :: Vind_old(:,:)
  real(dp), allocatable :: Vxc(:,:)

  complex(dp), allocatable :: V0r(:,:)
  complex(dp), allocatable :: V0r_comp(:,:)
  complex(dp), allocatable :: Vindr(:,:)

  type(SpglibDataset) :: dset

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Select model                !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  open(10,file='screening.in')
  read(10,*) par_NL    ! Include non-linearity?
  read(10,*) par_inter ! Include interband transitions?
  read(10,*) par_sigma ! Include effective sigma-band screening? (if par_NL)
  read(10,*) par_xc    ! Include exchange-correlation potential? (if par_NL)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Input parameters (physical) !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  read(10,*) par_n   ! Supercell size (n x n)
  read(10,*) mu      ! Chemical potential [energy] (eV)
  read(10,*) eps     ! Background dielectric constant (if .not. par_sigma)
  read(10,*) num_imp ! Number of impurities
  if (num_imp>0) then
    allocate(imp_pos(2,num_imp))
    allocate(par_h(num_imp))
    allocate(Z(num_imp))
    read(10,*) imp_pos ! In-plane position of impurities [length] (Ang)
    read(10,*) par_h   ! Height of impurities above plane [length] (Ang)
    read(10,*) Z       ! Impurity charges [charge] (e)
    !Z=0.1_dp*0.5_dp*vF
  end if
  read(10,*) lat_par ! Carbon-carbon distance [length] (Ang)
  read(10,*) vF      ! Fermi velocity [length] [time^-1] (a.u.)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Input parameters (numerical) !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  read(10,*) par_f    ! Fineness of real-space grid (m x m, where m=3*n*f)
  read(10,*) max_iter ! Maximum number of self-consistency iterations (if par_NL)
  read(10,*) amix     ! Potential mixing factor for self-consistency loop (if par_NL)
  read(10,*) rms_tol  ! Convergence tolerance for self-consistency loop (if par_NL) [energy] (Ha)

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Input parameters (symmetry) !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  read(10,*) par_k        ! Monkhorst-Pack k-point mesh (k x k including Gamma)
  read(10,*) calc_k_symm  ! Calculate k-point symmetry?
  read(10,*) calc_at_symm ! Calculate site symmetry?
  close(10)

  !~~~~~~~~~~~~~~~~~~!
  ! Setup parameters !
  !~~~~~~~~~~~~~~~~~~!

  if (.not. par_NL) par_sigma=.false.
  N=par_n*3*par_f
  mu=mu/Ha2eV
  kF=abs(mu)/vF
  if (par_sigma) eps=1.0_dp
  par_h=par_h/bohr2A
  imp_pos=imp_pos/bohr2A
  lat_par=lat_par/bohr2A
  Z=-Z

  !~~~~~~~~~~~~~~~~!
  ! Setup geometry !
  !~~~~~~~~~~~~~~~~!

  par_l=sqrt(3.0_dp)*lat_par
  cell(1:2,1)=(/sqrt(3.0_dp)*0.5_dp*par_l, 0.5_dp*par_l/)
  cell(1:2,2)=(/sqrt(3.0_dp)*0.5_dp*par_l,-0.5_dp*par_l/)
  scell=par_n*cell
  area=0.5_dp*sqrt(3.0_dp)*(par_l*par_n)**2
  scellr(1:2,1)=(/Pi/scell(1,1),Pi/scell(2,1)/)
  scellr(1:2,2)=(/Pi/scell(1,2),Pi/scell(2,2)/)
  arear=0.5_dp*sqrt(3.0_dp)*(scellr(1,1)**2+scellr(2,1)**2)

  size_H=2*(par_n**2)
  allocate(at_poss(2,size_H))
  allocate(at_possf(2,size_H))
  at_pos(1:2,1)=(/       sqrt(3.0_dp)*par_l/3.0_dp,0.0_dp/)
  at_pos(1:2,2)=(/2.0_dp*sqrt(3.0_dp)*par_l/3.0_dp,0.0_dp/)
  at_posf(1:2,1)=(/1.0_dp/3.0_dp,1.0_dp/3.0_dp/)
  at_posf(1:2,2)=(/2.0_dp/3.0_dp,2.0_dp/3.0_dp/)
  c=0
  do a=0,par_n-1
    do b=0,par_n-1
      shift(1:2)=a*cell(1:2,1)+&
                 b*cell(1:2,2)
      shiftf(1:2)=(/real(a,dp)/real(par_n,dp),&
                    real(b,dp)/real(par_n,dp)/)
      c=c+1
      at_poss(1:2,c)=at_pos(1:2,1)+shift(1:2)
      at_possf(1:2,c)=at_posf(1:2,1)/real(par_n,dp)+shiftf(1:2)
      c=c+1
      at_poss(1:2,c)=at_pos(1:2,2)+shift(1:2)
      at_possf(1:2,c)=at_posf(1:2,2)/real(par_n,dp)+shiftf(1:2)
    end do
  end do

  !~~~~~~~~~~~~~~~~~~~!
  ! Symmetry analysis !
  !~~~~~~~~~~~~~~~~~~~!

  if (calc_k_symm .or. calc_at_symm) then
    lattice=0.0_dp
    do i=1,2
      do j=1,2
        lattice(j,i)=scell(i,j)
      end do
    end do
    lattice(3,3)=10.0_dp
    allocate(positions(3,size_H+num_imp))
    do i=1,size_H
      positions(1:3,i)=(/modulo(at_possf(1,i),1.0_dp),modulo(at_possf(2,i),1.0_dp),0.0_dp/)
    end do
    allocate(atom_types(size_H+num_imp))
    atom_types=1
    do a=1,num_imp
      positions(1:3,size_H+a)=(/imp_pos(1,a),imp_pos(2,a),par_h(a)/)
      atom_types(size_H+a)=1+a
      do b=1,a-1
        if (abs(Z(a)-Z(b))<1.0d-6) then
          atom_types(size_H+a)=atom_types(size_H+b)
          exit
        end if
      end do
    end do
  end if

  if (calc_at_symm) then
  if (0==0) then
    write(par_n_char,'(i10)') par_n
    file_name='weights'//trim(adjustl(par_n_char))
    open(20,file=trim(adjustl(file_name)))
    dset=spg_get_dataset(lattice,positions,atom_types,size_H+num_imp,1.0d-6)
    c=0
    do i=1,size_H
      if (i-1==dset%equivalent_atoms(i)) c=c+1
    end do
    write(20,'(i7)') c
    do i=1,size_H
      if (i-1==dset%equivalent_atoms(i)) then
        weight=count(dset%equivalent_atoms==i-1)
        write(20,'(i7,1x,i7)') i, weight-1
        do j=1,size_H
          if ((i/=j) .and. (dset%equivalent_atoms(i)==dset%equivalent_atoms(j))) write(20,'(i7)') j
        end do
      end if
    end do
    close(20)
  else
    at_poss=at_poss*bohr2A
    scell=scell*bohr2A
    write(par_n_char,'(i10)') par_n
    file_name='weights'//trim(adjustl(par_n_char))
    open(20,file=trim(adjustl(file_name)))
    num_weights=0
    allocate(iweights(size_H))
    allocate(rweights(2,size_H))
    iweights=0
    rweights=0.0_dp
    do i=1,size_H
      d_min=9999999.9_dp
      do a=-1,1
        do b=-1,1
          shift(1:2)=a*scell(1:2,1)+&
                     b*scell(1:2,2)
          d=(at_poss(1,i)-shift(1))**2+&
            (at_poss(2,i)-shift(2))**2
          if (abs(d-d_min)<1.0d-5) then
            if (abs(a)+abs(b)==0) then
              c=a
              j=b
            end if
          else if (d<d_min) then
            d_min=d
            c=a
            j=b
          end if
        end do
      end do
      if (abs(c)+abs(j)==0) then
        if (at_poss(2,i)>=0.0_dp) then
          d=sqrt(at_poss(1,i)**2+&
                 at_poss(2,i)**2)
          if (at_poss(2,i)==0.0_dp) then
            num_weights=num_weights+1
            iweights(num_weights)=i
            rweights(1,num_weights)=d
            rweights(2,num_weights)=6.0_dp
          else
            num_weights=num_weights+1
            iweights(num_weights)=i
            rweights(1,num_weights)=d
            rweights(2,num_weights)=12.0_dp
          end if
        end if
      end if
    end do
    write(20,'(i7)') num_weights
    do i=1,num_weights
      write(20,'(i7,2(1x,es22.15e2))') iweights(i), rweights(1:2,i)
    end do
    close(20)
    deallocate(rweights)
    deallocate(iweights)
    at_poss=at_poss/bohr2A
    scell=scell/bohr2A
  end if
  else
    write(par_n_char,'(i10)') par_n
    file_name='weights'//trim(adjustl(par_n_char))
    open(20,file=trim(adjustl(file_name)))
    write(20,'(i7)') size_H
    do i=1,size_H
      write(20,'(i7,1x,i7)') i, 0
    end do
    close(20)
  end if

  if (calc_k_symm) then

    mesh(1:3)=(/par_k,par_k,1/)
    allocate(grid_point(3,product(mesh)))
    grid_point=0
    allocate(map(product(mesh)))
    map=0
    is_shift=(/0,0,0/)
    num_kpoints=spg_get_ir_reciprocal_mesh(grid_point,map,mesh,is_shift,0,lattice,&
                                           positions,atom_types,size_H+num_imp,1.0d-6)

    max_size=24*par_n*par_n
    allocate(rotation(3,3,max_size))
    allocate(translation(3,max_size))

    num_symops=spg_get_symmetry(rotation,translation,max_size,lattice,positions,&
                                atom_types,size_H+num_imp,1.0d-6)

    allocate(write_rot(num_symops))
    write_rot=.false.
    write(par_k_char,'(i10)') par_k
    file_name='kweights'//trim(adjustl(par_n_char))//'_'//trim(adjustl(par_k_char))
    open(20,file=trim(adjustl(file_name)))
    write(20,'(i7)') num_kpoints
    do i=1,product(mesh)
      if (i-1==map(i)) then
        weight=count(map==i-1)
        k_point(1:2)=real(grid_point(1:2,i),dp)/real(mesh(1:2),dp)
        write(20,'(3(es22.15e2,1x),i7)') k_point(1:2), real(weight,dp)/real(product(mesh),dp), weight-1
        c=0
        do j=1,product(mesh)
          if (map(i)==map(j)) then
            found=.false.
            k_point2(1:2)=real(grid_point(1:2,j),dp)/real(mesh(1:2),dp)
            do k=1,num_symops
              new_pos(1:2)=rotation(1:2,1,k)*k_point(1)+&
                           rotation(1:2,2,k)*k_point(2)
              new_pos=modulo(new_pos,1.0_dp)
              if ((abs(new_pos(1)-modulo(k_point2(1),1.0_dp))<1.0d-5) .and. &
                  (abs(new_pos(2)-modulo(k_point2(2),1.0_dp))<1.0d-5)) then
                c=c+1
                if ((j/=i) .and. &
                    (.not. found)) then
                  found=.true.
                  write_rot(k)=.true.
                  write(20,'(i7)') k
                end if
              end if
            end do
          end if
        end do
        if (c/=num_symops) then
          print*, "kSYMOPS ERROR!", i, c
          stop
        end if
      end if
    end do
    close(20)

    if (any(write_rot)) then
      file_name='sym_ops'//trim(adjustl(par_n_char))//'_'//trim(adjustl(par_k_char))
      open(20,file=trim(adjustl(file_name)))
      write(20,'(i7,2(1x,i7))') num_symops, count(write_rot), size_H
    end if
    do i=1,num_symops
      if (write_rot(i)) then
        allocate(at_t(size_H))
        do j=1,size_H
          new_pos(1:2)=rotation(1,1:2,i)*at_possf(1,j)+&
                       rotation(2,1:2,i)*at_possf(2,j)
          new_pos(1:2)=new_pos(1:2)+translation(1:2,i)
          new_pos=modulo(new_pos,1.0_dp)
          c=0
          do k=1,size_H
            if ((abs(new_pos(1)-modulo(at_possf(1,k),1.0_dp))<1.0d-5) .and. &
                (abs(new_pos(2)-modulo(at_possf(2,k),1.0_dp))<1.0d-5)) then
              at_t(j)=k
              c=c+1
            end if
          end do
          if (c/=1) then
            print*, "aSYMOPS ERROR!", i, j, c
            stop
          end if
        end do
        write(20,'(i7)') i
        do j=1,size_H
          write(20,'(i7)') at_t(j)
        end do
        deallocate(at_t)
      end if
    end do
    if (any(write_rot)) close(20)

  else

    num_kpoints=par_k*par_k
    write(par_k_char,'(i10)') par_k
    file_name='kweights'//trim(adjustl(par_n_char))//'_'//trim(adjustl(par_k_char))
    open(20,file=trim(adjustl(file_name)))
    write(20,'(i7)') num_kpoints
    do i=0,par_k-1
      do j=0,par_k-1
        write(20,'(3(es22.15e2,1x),i7)') real(i,dp)/par_k, &
                                         real(j,dp)/par_k, &
                                         1.0_dp/num_kpoints, &
                                         0
      end do
    end do
    close(20)

  end if

  if (allocated(write_rot)) deallocate(write_rot)
  if (allocated(translation)) deallocate(translation)
  if (allocated(rotation)) deallocate(rotation)
  if (allocated(atom_types)) deallocate(atom_types)
  if (allocated(positions)) deallocate(positions)
  if (allocated(map)) deallocate(map)
  if (allocated(grid_point)) deallocate(grid_point)

  !~~~~~~~~~~~~~!
  ! Setup grids !
  !~~~~~~~~~~~~~!

  allocate(r(2,N,N))
  allocate(rd(N,N))
  allocate(rr(2,N/2+1,N))
  allocate(rrp(2,N/2+1,N))
  allocate(rdr(N/2+1,N))
  allocate(V0(N,N))
  allocate(V0r(N/2+1,N))
  allocate(V0r_comp(N/2+1,N))
  if (par_NL) then
    allocate(V(N,N))
    allocate(Vmix(N,N))
    allocate(rho(N,N))
    allocate(Vind(N,N))
    allocate(Vind_old(N,N))
    allocate(Vindr(N/2+1,N))
    if (par_xc) allocate(Vxc(N,N))
  end if
  call dfftw_plan_dft_r2c_2d(plan_r2c,N,N,V0,V0r,FFTW_ESTIMATE)
  call dfftw_plan_dft_c2r_2d(plan_c2r,N,N,V0r,V0,FFTW_ESTIMATE)

  do i=1,N
    do j=1,N
      r(1:2,i,j)=(real(i-1,dp)/N)*scell(1:2,1)+&
                 (real(j-1,dp)/N)*scell(1:2,2)
      rd(i,j)=9999999.9_dp
      do a=-1,1
        do b=-1,1
          shift(1:2)=a*scell(1:2,1)+&
                     b*scell(1:2,2)
          d=sqrt((r(1,i,j)+shift(1))**2+&
                 (r(2,i,j)+shift(2))**2)
          if (d<rd(i,j)) rd(i,j)=d
        end do
      end do
    end do
  end do

  do i=1,N/2+1
    do j=1,N
      rr(1:2,i,j)=real(i-1,dp)*scellr(1:2,1)+&
                  real(j-1,dp)*scellr(1:2,2)
      rdr(i,j)=9999999.9_dp
      do a=-1,1
        do b=-1,1
          shift(1:2)=a*N*scellr(1:2,1)+&
                     b*N*scellr(1:2,2)
          d=sqrt((rr(1,i,j)+shift(1))**2+&
                 (rr(2,i,j)+shift(2))**2)
          if (d<rdr(i,j)) then
            rrp(1:2,i,j)=rr(1:2,i,j)+shift(1:2)
            rdr(i,j)=d
          end if
        end do
      end do
    end do
  end do

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Calculate impurity potential with linear-response screening !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  V0r=cmplx_0
  do a=1,num_imp
    do i=1,N/2+1
      do j=1,N
        if (par_sigma) then
          d=(eps_sigma*(eps_sigma+1.0_dp-(eps_sigma-1.0_dp)*exp(-rdr(i,j)*d_sigma)))/&
            (eps_sigma+1.0_dp+(eps_sigma-1.0_dp)*exp(-rdr(i,j)*d_sigma))
        else
          d=1.0_dp
        end if
        if (par_NL) then
          V0r_comp(i,j)=cmplx_1
        else
          V0r_comp(i,j)=cmplx(1.0_dp+4.0_dp*kF/(eps*vF*rdr(i,j)),0.0_dp)
        end if
        if (par_inter .and. (rdr(i,j)>2.0*kF)) then
          x=2.0_dp*kF/rdr(i,j)
          V0r_comp(i,j)=exp(-rdr(i,j)*par_h(a))/((V0r_comp(i,j)*d-(x*sqrt(1.0_dp-x**2)-acos(x))/(eps*vF))*rdr(i,j))
        else
          V0r_comp(i,j)=exp(-rdr(i,j)*par_h(a))/(V0r_comp(i,j)*d*rdr(i,j))
        end if
        V0r_comp(i,j)=V0r_comp(i,j)*Z(a)*exp(-cmplx_i*(imp_pos(1,a)*rrp(1,i,j)+&
                                                       imp_pos(2,a)*rrp(2,i,j)))
      end do
    end do
    V0r=V0r+V0r_comp
  end do
  V0r(1,1)=cmplx_0
  call dfftw_execute_dft_c2r(plan_c2r,V0r,V0)
  V0=V0*2.0_dp*Pi/(eps*area)

  if (.not. par_NL) then
    open(20,file='potential.dat')
    do k=1,size_H
      do i=1,N/par_f
        a=par_f*(i-1)+1
        do j=1,N/par_f
          b=par_f*(j-1)+1
          d=(r(1,a,b)-at_poss(1,k))**2+&
            (r(2,a,b)-at_poss(2,k))**2
          if (d<1.0d-5) write(20,'(es22.15e2,2(1x,es22.15e2))') r(1:2,a,b)*bohr2A, V0(a,b)*Ha2eV
        end do
      end do
    end do
    close(20)
  end if

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Self-consistency loop for non-linear screening !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  if (par_NL) then

    Vmix=V0/(eps**2)
    Vind_old=0.0_dp
    conv=.false.
    do iter=1,max_iter

      rho=(mu-Vmix)*abs(mu-Vmix)
      if (par_xc) then
        do i=1,N
          do j=1,N
            Vxc(i,j)=-sign(1.0_dp,rho(i,j))*alpha_xc*sqrt(abs(rho(i,j)/(Pi*(vF**2))))*log(abs(rho(i,j)/(Pi*(vF**2))))
          end do
        end do
      end if
      call dfftw_execute_dft_r2c(plan_r2c,rho,Vindr)
      do i=1,N/2+1
        do j=1,N
          if (par_sigma) then
            d=(eps_sigma*(eps_sigma+1.0_dp-(eps_sigma-1.0_dp)*exp(-rdr(i,j)*d_sigma)))/&
              (eps_sigma+1.0_dp+(eps_sigma-1.0_dp)*exp(-rdr(i,j)*d_sigma))
          else
            d=1.0_dp
          end if
          if (par_inter .and. (rdr(i,j)>2.0*kF)) then
            x=2.0_dp*kF/rdr(i,j)
            Vindr(i,j)=Vindr(i,j)/((d-(x*sqrt(1.0_dp-x**2)-acos(x))/(eps*vF))*rdr(i,j))
          else
            Vindr(i,j)=Vindr(i,j)/(d*rdr(i,j))
          end if
        end do
      end do
      Vindr(1,1)=cmplx_0
      call dfftw_execute_dft_c2r(plan_c2r,Vindr,Vind)
      Vind=Vind*2.0_dp/(eps*(vF**2)*(N**2))
      rms=0.0_dp
      if (par_xc) then
        do i=1,N
          do j=1,N
            rms=rms+(Vind(i,j)+Vxc(i,j)-Vind_old(i,j))**2
          end do
        end do
      else
        do i=1,N
          do j=1,N
            rms=rms+(Vind(i,j)-Vind_old(i,j))**2
          end do
        end do
      end if
      rms=sqrt(rms/(N**2))
      print('(i5,1x,es12.5e2)'), iter, rms
      if (rms<rms_tol) conv=.true.
      if (par_xc) then
        Vind_old=Vind+Vxc
        V=V0+Vind+Vxc
      else
        Vind_old=Vind
        V=V0+Vind
      end if
      Vmix=amix*V+(1.0_dp-amix)*Vmix

      if (conv .or. (iter==max_iter)) then
        charge2D=0.0_dp
        do i=1,N
          do j=1,N
            charge2D=charge2D+rho(i,j)
          end do
        end do
        charge2D=charge2D*area/((N**2)*Pi*(vF**2))
        print('(a1,5x,es12.5e2)'), '#', -charge2D
        open(20,file='potential.dat')
        do k=1,size_H
          do i=1,N/par_f
            a=par_f*(i-1)+1
            do j=1,N/par_f
              b=par_f*(j-1)+1
              d=(r(1,a,b)-at_poss(1,k))**2+&
                (r(2,a,b)-at_poss(2,k))**2
              if (d<1.0d-5) write(20,'(es22.15e2,2(1x,es22.15e2))') r(1:2,a,b)*bohr2A, V(a,b)*Ha2eV
            end do
          end do
        end do
        close(20)
      end if

      if (conv) exit

    end do

  end if

  !~~~~~~~~~~!
  ! Clean up !
  !~~~~~~~~~~!

  call dfftw_destroy_plan(plan_c2r)
  call dfftw_destroy_plan(plan_r2c)
  if (par_NL) then
    if (par_xc) deallocate(Vxc)
    deallocate(Vindr)
    deallocate(Vind_old)
    deallocate(Vind)
    deallocate(rho)
    deallocate(Vmix)
    deallocate(V)
  end if
  deallocate(V0r_comp)
  deallocate(V0r)
  deallocate(V0)
  deallocate(rdr)
  deallocate(rrp)
  deallocate(rr)
  deallocate(rd)
  deallocate(r)
  deallocate(at_possf)
  deallocate(at_poss)
  if (num_imp>0) then
    deallocate(Z)
    deallocate(par_h)
    deallocate(imp_pos)
  end if

end program screening
