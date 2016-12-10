program TB
  use MatrixSwitch
  use fermi_dirac_stats

  implicit none
  include 'mpif.h'

  integer, parameter :: dp=selected_real_kind(15,300)

  type k_info
    integer :: num
    integer, allocatable :: list(:)
  end type k_info

  real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: Ha2eV=27.21138342902473_dp
  real(dp), parameter :: bohr2A=0.52917720859_dp
  real(dp), parameter :: k_B=8.6173324d-5

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp), cmplx_i=(0.0_dp,1.0_dp), cmplx_0=(0.0_dp,0.0_dp)

  character(1) :: c1
  character(5) :: m_storage, m_dstorage, m_zstorage
  character(3) :: m_operation
  character(10) :: at_name
  character(50) :: weights_file, kweights_file, sym_ops_file, ldos_file

  logical :: par_pot, selected_eigs, is_gamma, have_symm

  integer :: mpi_err, mpi_size, mpi_rank
  integer :: i, j, k, k1, k2, a, b, c, s
  integer :: par_n(2), num_kpoints, num_symops, size_H, ldwork, lzwork, lawork, info, size_dos, par_2D, par_bs, num_ldos
  integer :: pot_n, rmax
  integer :: num_eigs, nz, liwork
  integer, allocatable :: nn(:,:,:), ldos_at(:)
  integer, allocatable :: iwork(:), ifail(:), iclustr(:)

  real(dp) :: cell(2,2), scell(2,2), scellr(2,2), at_pos(2,2), shift(2)
  real(dp) :: par_T, par_l, par_b, par_c, par_lr, kp(2), dos_E1, dos_E2, dos_dE, dos_w, E, grad_fd, grad2_fd, E_F
  real(dp) :: par_mu, pot_dr, area, area_per_atom, r1, r2
  real(dp) :: occ_tot, occ_err, dist, kdotr, norm_dist1, norm_dist2
  real(dp) :: abstol, orfac, del
  real(dp) :: dr, r, pot_point, pot_av, sum
  real(dp), allocatable :: eigvals(:), occs(:), dos(:), at_poss(:,:), at_dist(:), at_dos(:,:), pot(:), density(:)
  real(dp), allocatable :: gap(:), dwork(:)
  real(dp), allocatable :: av(:,:), kpoints(:,:)

  complex(dp) :: zel
  complex(dp), allocatable :: zwork(:), awork(:)

  type(matrix) :: H, eigvecs

  type(k_info), allocatable :: sym_kpoints(:), at_t(:)

  call mpi_init(mpi_err)
  call mpi_comm_size(mpi_comm_world,mpi_size,mpi_err)
  call mpi_comm_rank(mpi_comm_world,mpi_rank,mpi_err)

  m_dstorage='pddbc'
  m_zstorage='pzdbc'
  m_operation='lap'

  if (mpi_rank==0) then
    open(10,file='TB.in')
    read(10,*) par_n(1:2)
    read(10,*) par_l
    read(10,*) par_b
    read(10,*) par_pot
    read(10,*) par_T
    read(10,*) par_mu
    read(10,*) selected_eigs
    read(10,*) dos_E1
    read(10,*) dos_E2
    read(10,*) dos_dE
    read(10,*) dos_w
    read(10,*) weights_file
    open(20,file=trim(adjustl(weights_file)))
    read(20,'(i7)') num_ldos
    allocate(ldos_at(num_ldos))
    do i=1,num_ldos
      read(20,'(i7,1x,i7)') ldos_at(i), a
      do j=1,a
        read(20,'()')
      end do
    end do
    close(20)
    read(10,*) kweights_file
    open(20,file=trim(adjustl(kweights_file)))
    read(20,'(i7)') num_kpoints
    allocate(kpoints(1:4,num_kpoints))
    allocate(sym_kpoints(num_kpoints))
    have_symm=.false.
    do i=1,num_kpoints
      read(20,'(3(es22.15e2,1x),i7)') kpoints(1:3,i), sym_kpoints(i)%num
      kpoints(4,i)=kpoints(3,i)/(sym_kpoints(i)%num+1)
      if (sym_kpoints(i)%num/=0) then
        have_symm=.true.
        allocate(sym_kpoints(i)%list(sym_kpoints(i)%num))
        do j=1,sym_kpoints(i)%num
          read(20,'(i7)') sym_kpoints(i)%list(j)
        end do
      end if
    end do
    close(20)
    if (have_symm) then
      read(10,*) sym_ops_file
      open(20,file=trim(adjustl(sym_ops_file)))
      read(20,'(i7,2(1x,i7))') num_symops, b, c
      allocate(at_t(num_symops))
      do i=1,num_symops
        at_t(i)%num=0
      end do
      do i=1,b
        read(20,'(i7)') k
        at_t(k)%num=c
        allocate(at_t(k)%list(c))
        do j=1,c
          read(20,'(i7)') at_t(k)%list(j)
        end do
      end do
      close(20)
    end if
    read(10,*) par_2D
    read(10,*) par_bs
    close(10)
  end if
  call mpi_bcast(par_n,2,mpi_int,0,mpi_comm_world,info)
  call mpi_bcast(par_l,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(par_b,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(par_pot,1,mpi_logical,0,mpi_comm_world,info)
  call mpi_bcast(par_T,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(par_mu,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(selected_eigs,1,mpi_logical,0,mpi_comm_world,info)
  call mpi_bcast(dos_E1,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(dos_E2,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(dos_dE,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(dos_w,1,mpi_double_precision,0,mpi_comm_world,info)
  call mpi_bcast(num_ldos,1,mpi_int,0,mpi_comm_world,info)
  if (.not. allocated(ldos_at)) allocate(ldos_at(num_ldos))
  call mpi_bcast(ldos_at,num_ldos,mpi_int,0,mpi_comm_world,info)
  call mpi_bcast(num_kpoints,1,mpi_int,0,mpi_comm_world,info)
  if (.not. allocated(kpoints)) allocate(kpoints(1:4,num_kpoints))
  call mpi_bcast(kpoints,4*num_kpoints,mpi_double_precision,0,mpi_comm_world,info)
  if (mpi_rank/=0) allocate(sym_kpoints(num_kpoints))
  do i=1,num_kpoints
    call mpi_bcast(sym_kpoints(i)%num,1,mpi_int,0,mpi_comm_world,info)
    if (sym_kpoints(i)%num/=0) then
      if (mpi_rank/=0) allocate(sym_kpoints(i)%list(sym_kpoints(i)%num))
      call mpi_bcast(sym_kpoints(i)%list,sym_kpoints(i)%num,mpi_int,0,mpi_comm_world,info)
    end if
  end do
  call mpi_bcast(have_symm,1,mpi_logical,0,mpi_comm_world,info)
  if (have_symm) then
    call mpi_bcast(num_symops,1,mpi_int,0,mpi_comm_world,info)
    if (mpi_rank/=0) allocate(at_t(num_symops))
    do i=1,num_symops
      call mpi_bcast(at_t(i)%num,1,mpi_int,0,mpi_comm_world,info)
      if (at_t(i)%num/=0) then
        if (mpi_rank/=0) allocate(at_t(i)%list(at_t(i)%num))
        call mpi_bcast(at_t(i)%list,at_t(i)%num,mpi_int,0,mpi_comm_world,info)
      end if
    end do
  end if
  call mpi_bcast(par_2D,1,mpi_int,0,mpi_comm_world,info)
  call mpi_bcast(par_bs,1,mpi_int,0,mpi_comm_world,info)

  call ms_scalapack_setup(mpi_size,par_2D,'c',par_bs)

  norm_dist1=2.0_dp*(dos_w**2)
  norm_dist2=dos_w*sqrt(2.0_dp*Pi)

  par_T=k_B*par_T
  par_l=sqrt(3.0_dp)*par_l
  cell(1:2,1)=(/sqrt(3.0_dp)*0.5_dp*par_l, 0.5_dp*par_l/)
  cell(1:2,2)=(/sqrt(3.0_dp)*0.5_dp*par_l,-0.5_dp*par_l/)
  scell(1:2,1)=par_n(1)*cell(1:2,1)
  scell(1:2,2)=par_n(2)*cell(1:2,2)
  scellr(1:2,1)=(/Pi/scell(1,1),Pi/scell(2,1)/)
  scellr(1:2,2)=(/Pi/scell(1,2),Pi/scell(2,2)/)
  size_H=2*(par_n(1)*par_n(2))
  allocate(at_poss(2,size_H))
  allocate(at_dist(size_H))
  allocate(nn(3,3,size_H))
  allocate(eigvals(size_H))
  allocate(occs(size_H))
  allocate(density(num_ldos))

  at_pos(1:2,1)=(/       sqrt(3.0_dp)*par_l/3.0_dp,0.0_dp/)
  at_pos(1:2,2)=(/2.0_dp*sqrt(3.0_dp)*par_l/3.0_dp,0.0_dp/)
  c=0
  do a=0,par_n(1)-1
    do b=0,par_n(2)-1
      shift(1:2)=a*cell(1:2,1)+&
                 b*cell(1:2,2)
      c=c+1
      at_poss(1:2,c)=at_pos(1:2,1)+shift(1:2)
      c=c+1
      at_poss(1:2,c)=at_pos(1:2,2)+shift(1:2)
    end do
  end do
  area=0.5_dp*sqrt(3.0_dp)*(par_l**2)*par_n(1)*par_n(2)
  area_per_atom=area/size_H

  do i=1,size_H
    c=0
    do j=1,size_H
      do a=-1,1
        do b=-1,1
          shift(1:2)=a*scell(1:2,1)+&
                     b*scell(1:2,2)
          dist=(at_poss(1,i)-at_poss(1,j)-shift(1))**2+&
               (at_poss(2,i)-at_poss(2,j)-shift(2))**2
          if (abs(dist-((par_l**2)/3.0_dp))<1.0d-5) then
            c=c+1
            if (c>3) then
              print*, 'ERROR1!', i, c
              stop
            else
              nn(1:3,c,i)=(/j,a,b/)
            end if
          end if
        end do
      end do
    end do
  end do

  at_dist=0.0_dp
  if (par_pot) then
    if (mpi_rank==0) then
      open(10,file='potential.dat')
      do i=1,size_H
        read(10,'(es22.15e2,2(1x,es22.15e2))') r1, r2, at_dist(i)
      end do
      close(10)
    end if
    call mpi_bcast(at_dist,size_H,mpi_double_precision,0,mpi_comm_world,info)
    del=0.0_dp
    do i=1,size_H
      del=del+at_dist(i)
    end do
    del=del/real(size_H,dp)
    at_dist=at_dist-del
  end if

  if (selected_eigs) then
    abstol=0.0_dp
    orfac=-1.0_dp
    allocate(ifail(size_H))
    allocate(iclustr(2*mpi_size))
    allocate(gap(mpi_size))
  end if

  size_dos=nint((dos_E2-dos_E1)/dos_dE)
  allocate(dos(size_dos))
  allocate(at_dos(size_dos,num_ldos))
  dos=0.0_dp
  at_dos=0.0_dp
  density=0.0_dp

  do k1=1,num_kpoints
      if ((kpoints(1,k1)<1.0d-12) .and. &
          (kpoints(2,k1)<1.0d-12)) then
        is_gamma=.true.
        m_storage=m_dstorage
      else
        is_gamma=.false.
        m_storage=m_zstorage
        kp(1:2)=kpoints(1,k1)*scellr(1:2,1)+&
                kpoints(2,k1)*scellr(1:2,2)
      end if
      call m_allocate(H,size_H,size_H,m_storage)
      call m_allocate(eigvecs,size_H,size_H,m_storage)
      call m_set(H,'a',0.0_dp,0.0_dp,m_operation)
      do i=1,size_H
        call m_set_element(H,i,i,at_dist(i),m_operation)
        do j=1,3
          if (is_gamma) then
            call m_get_element(H,i,nn(1,j,i),del,m_operation)
            call m_set_element(H,i,nn(1,j,i),del+par_b,m_operation)
          else
            shift(1:2)=nn(2,j,i)*scell(1:2,1)+&
                       nn(3,j,i)*scell(1:2,2)
            kdotr=kp(1)*shift(1)+&
                  kp(2)*shift(2)
            call m_get_element(H,i,nn(1,j,i),zel,m_operation)
            call m_set_element(H,i,nn(1,j,i),zel+cmplx(par_b,0.0_dp,dp)*exp(cmplx_i*kdotr),m_operation)
          end if
        end do
      end do
      if (selected_eigs) then
        if (is_gamma) then
          allocate(dwork(1))
          allocate(iwork(1))
          ldwork=-1
          liwork=-1
          call pdsyevx('v','v','u',size_H,H%dval,1,1,H%iaux1,dos_E1,dos_E2,0,0,abstol,&
                       num_eigs,nz,eigvals,orfac,eigvecs%dval,1,1,eigvecs%iaux1,&
                       dwork,ldwork,iwork,liwork,ifail,iclustr,gap,info)
          if (info/=0) print*, 'pdsyevx ERROR!'
          ldwork=nint(real(dwork(1),dp))
          liwork=iwork(1)
          deallocate(iwork)
          deallocate(dwork)
          allocate(dwork(max(1,ldwork)))
          allocate(iwork(max(1,liwork)))
          call pdsyevx('v','v','u',size_H,H%dval,1,1,H%iaux1,dos_E1,dos_E2,0,0,abstol,&
                       num_eigs,nz,eigvals,orfac,eigvecs%dval,1,1,eigvecs%iaux1,&
                       dwork,ldwork,iwork,liwork,ifail,iclustr,gap,info)
          if (info/=0) print*, 'pdsyevx ERROR!'
          deallocate(iwork)
          deallocate(dwork)
        else
          allocate(zwork(1))
          allocate(dwork(1))
          allocate(iwork(1))
          lzwork=-1
          ldwork=-1
          liwork=-1
          call pzheevx('v','v','u',size_H,H%zval,1,1,H%iaux1,dos_E1,dos_E2,0,0,abstol,&
                       num_eigs,nz,eigvals,orfac,eigvecs%zval,1,1,eigvecs%iaux1,&
                       zwork,lzwork,dwork,ldwork,iwork,liwork,ifail,iclustr,gap,info)
          if (info/=0) print*, 'pzheevx ERROR!'
          lzwork=nint(real(zwork(1),dp))
          ldwork=nint(real(dwork(1),dp))
          liwork=iwork(1)
          deallocate(iwork)
          deallocate(dwork)
          deallocate(zwork)
          allocate(zwork(max(1,lzwork)))
          allocate(dwork(max(1,ldwork)))
          allocate(iwork(max(1,liwork)))
          call pzheevx('v','v','u',size_H,H%zval,1,1,H%iaux1,dos_E1,dos_E2,0,0,abstol,&
                       num_eigs,nz,eigvals,orfac,eigvecs%zval,1,1,eigvecs%iaux1,&
                       zwork,lzwork,dwork,ldwork,iwork,liwork,ifail,iclustr,gap,info)
          if (info/=0) print*, 'pzheevx ERROR!'
          deallocate(iwork)
          deallocate(dwork)
          deallocate(zwork)
        end if
      else
        if (is_gamma) then
          allocate(dwork(1))
          ldwork=-1
          call pdsyev('v','u',size_H,H%dval,1,1,H%iaux1,eigvals,eigvecs%dval,1,1,eigvecs%iaux1,&
                      dwork,ldwork,info)
          if (info/=0) print*, 'pdsyev ERROR!'
          ldwork=nint(real(dwork(1),dp))
          deallocate(dwork)
          allocate(dwork(max(1,ldwork)))
          call pdsyev('v','u',size_H,H%dval,1,1,H%iaux1,eigvals,eigvecs%dval,1,1,eigvecs%iaux1,&
                      dwork,ldwork,info)
          if (info/=0) print*, 'pdsyev ERROR!'
          deallocate(dwork)
        else
          allocate(zwork(1))
          allocate(awork(1))
          lzwork=-1
          lawork=-1
          call pzheev('v','u',size_H,H%zval,1,1,H%iaux1,eigvals,eigvecs%zval,1,1,eigvecs%iaux1,&
                    zwork,lzwork,awork,lawork,info)
          if (info/=0) print*, 'pzheev ERROR!'
          lzwork=nint(real(zwork(1),dp))
          lawork=nint(real(awork(1),dp))
          deallocate(awork)
          deallocate(zwork)
          allocate(zwork(max(1,lzwork)))
          allocate(awork(max(1,lawork)))
          call pzheev('v','u',size_H,H%zval,1,1,H%iaux1,eigvals,eigvecs%zval,1,1,eigvecs%iaux1,&
                      zwork,lzwork,awork,lawork,info)
          if (info/=0) print*, 'pzheev ERROR!'
          deallocate(awork)
          deallocate(zwork)
        end if
        num_eigs=size_H
      end if
      occ_tot=FD_occ_tot_eigs(par_T,par_mu,size_H,eigvals,occs)
      do i=1,num_eigs
        do j=1,size_dos
          E=dos_E1+(j-1)*dos_dE
          dos(j)=dos(j)+kpoints(3,k1)*exp(-(E-eigvals(i))**2/norm_dist1)
        end do
        do k=1,num_ldos
          call m_get_element(eigvecs,ldos_at(k),i,zel,m_operation)
          density(k)=density(k)+kpoints(4,k1)*occs(i)*(abs(zel)**2)
          do j=1,size_dos
            E=dos_E1+(j-1)*dos_dE
            at_dos(j,k)=at_dos(j,k)+kpoints(4,k1)*(abs(zel)**2)*exp(-(E-eigvals(i))**2/norm_dist1)
          end do
          do s=1,sym_kpoints(k1)%num
            call m_get_element(eigvecs,at_t(sym_kpoints(k1)%list(s))%list(ldos_at(k)),i,zel,m_operation)
            density(k)=density(k)+kpoints(4,k1)*occs(i)*(abs(zel)**2)
            do j=1,size_dos
              E=dos_E1+(j-1)*dos_dE
              at_dos(j,k)=at_dos(j,k)+kpoints(4,k1)*(abs(zel)**2)*exp(-(E-eigvals(i))**2/norm_dist1)
            end do
          end do
        end do
      end do
      call m_deallocate(eigvecs)
      call m_deallocate(H)
      if (mpi_rank==0) call progress_bar(1,num_kpoints,k1)
  end do
  dos=dos/norm_dist2
  at_dos=at_dos/norm_dist2

  if (mpi_rank==0) then
    open(11,file='density.dat')
    do k=1,num_ldos
      write(11,'(es22.15e2,3(1x,es22.15e2))') at_poss(1:2,ldos_at(k)), density(k), 2.0_dp*density(k)/area_per_atom
    end do
    close(11)
  end if

  if (mpi_rank==0) then
    open(11,file='dos.dat')
    do k=1,num_ldos
      write(at_name,'(i10)') ldos_at(k)
      ldos_file='ldos.'//trim(adjustl(at_name))//'.dat'
      open(10000+k,file=trim(adjustl(ldos_file)))
      write(10000+k,'(a1,3(1x,es22.15e2))') '#', at_poss(1:2,ldos_at(k)), at_dist(ldos_at(k))
    end do
    do j=1,size_dos
      E=dos_E1+(j-1)*dos_dE
      write(11,'(es22.15e2,1x,es22.15e2)') E, dos(j)
      do k=1,num_ldos
        write(10000+k,'(es22.15e2,1x,es22.15e2)') E, at_dos(j,k)
      end do
    end do
    do k=1,num_ldos
      close(10000+k)
    end do
    close(11)
  end if

  deallocate(at_dos)
  deallocate(dos)
  if (selected_eigs) then
    deallocate(gap)
    deallocate(iclustr)
    deallocate(ifail)
  end if
  deallocate(density)
  deallocate(occs)
  deallocate(eigvals)
  deallocate(nn)
  deallocate(at_dist)
  deallocate(at_poss)
  if (have_symm) then
    do i=1,num_symops
      if (at_t(i)%num/=0) deallocate(at_t(i)%list)
    end do
    deallocate(at_t)
  end if
  do i=1,num_kpoints
    if (sym_kpoints(i)%num/=0) deallocate(sym_kpoints(i)%list)
  end do
  deallocate(sym_kpoints)
  deallocate(kpoints)
  deallocate(ldos_at)

  call mpi_finalize(mpi_err)

end program TB
