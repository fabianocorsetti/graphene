program adatom_post
  implicit none

  integer, parameter :: dp=selected_real_kind(15,300)

  type k_info
    integer :: num
    integer, allocatable :: list(:)
  end type k_info

  real(dp), parameter :: Pi=3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: Ha2eV=27.21138342902473_dp
  real(dp), parameter :: bohr2A=0.52917720859_dp
  real(dp), parameter :: k_B=8.6173324d-5

  character(10) :: at_name
  character(50) :: weights_file, ldos_file, rdos_file

  logical :: par_pot, selected_eigs

  integer :: i, j, k, a, b, c
  integer :: par_n(2), num_ldos, size_H, num_bins, size_dos
  integer, allocatable :: ldos_at(:), r_bins(:), r_count(:)

  real(dp) :: par_l, par_b, par_T, par_mu, dos_E1, dos_E2, dos_dE, dos_w, E
  real(dp) :: r1, r2, cell(2,2), scell(2,2), at_pos(2,2), shift(2), density(2), potential
  real(dp), allocatable :: rdensity(:,:), at_poss(:,:), at_dist(:), ldos(:), rdos(:,:), rpotential(:)

  type(k_info), allocatable :: at_sym(:)

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
  allocate(at_sym(num_ldos))
  do i=1,num_ldos
    read(20,'(i7,1x,i7)') ldos_at(i), at_sym(i)%num
    allocate(at_sym(i)%list(0:at_sym(i)%num))
    at_sym(i)%list(0)=ldos_at(i)
    do j=1,at_sym(i)%num
      read(20,'(i7)') at_sym(i)%list(j)
    end do
  end do
  close(20)

  if (par_n(1)/=par_n(2)) then
    print*, "ERROR! Supercell not n x n"
    stop
  end if

  par_l=sqrt(3.0_dp)*par_l
  cell(1:2,1)=(/sqrt(3.0_dp)*0.5_dp*par_l, 0.5_dp*par_l/)
  cell(1:2,2)=(/sqrt(3.0_dp)*0.5_dp*par_l,-0.5_dp*par_l/)
  scell(1:2,1)=par_n(1)*cell(1:2,1)
  scell(1:2,2)=par_n(2)*cell(1:2,2)
  size_H=2*(par_n(1)*par_n(2))
  allocate(at_poss(2,size_H))
  allocate(at_dist(size_H))

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

  allocate(r_bins(size_H))
  r_bins=0

  do i=1,size_H
    at_dist=999999.9_dp
    do a=-1,1
      do b=-1,1
        shift(1:2)=a*scell(1:2,1)+&
                   b*scell(1:2,2)
        r1=sqrt((at_poss(1,i)+shift(1))**2+&
                (at_poss(2,i)+shift(2))**2)
        if (r1<at_dist(i)) then
          at_dist(i)=r1
        end if
      end do
    end do
    r_bins(i)=nint(at_dist(i))
  end do

  num_bins=maxval(r_bins)
  allocate(r_count(num_bins))
  r_count=0
  do i=1,num_bins
    r_count(i)=count(r_bins==i)
  end do

  allocate(rdensity(2,num_bins))
  rdensity=0.0_dp
  open(11,file='density.dat')
  do i=1,num_ldos
    read(11,'(es22.15e2,3(1x,es22.15e2))') r1, r2, density(1:2)
    do j=0,at_sym(i)%num
      a=r_bins(at_sym(i)%list(j))
      rdensity(1:2,a)=rdensity(1:2,a)+density(1:2)
    end do
  end do
  close(11)

  open(12,file='rdensity.dat')
  do k=1,num_bins
    if (r_count(k)/=0) write(12,'(es22.15e2,2(1x,es22.15e2))') real(k,dp), rdensity(1:2,k)/r_count(k)
  end do
  close(12)
  deallocate(rdensity)

  size_dos=nint((dos_E2-dos_E1)/dos_dE)
  allocate(ldos(size_dos))
  allocate(rdos(size_dos,num_bins))
  rdos=0.0_dp
  allocate(rpotential(num_bins))
  rpotential=0.0_dp
  do i=1,num_ldos
    write(at_name,'(i10)') ldos_at(i)
    ldos_file='ldos.'//trim(adjustl(at_name))//'.dat'
    open(11,file=trim(adjustl(ldos_file)))
    read(11,'(1x,3(1x,es22.15e2))') r1, r2, potential
    do j=1,size_dos
      read(11,'(es22.15e2,1x,es22.15e2)') r1, ldos(j)
    end do
    close(11)
    do j=0,at_sym(i)%num
      a=r_bins(at_sym(i)%list(j))
      rdos(1:size_dos,a)=rdos(1:size_dos,a)+ldos(1:size_dos)
      rpotential(a)=rpotential(a)+potential
    end do
  end do

  do i=1,num_bins
    if (r_count(i)/=0) then
      write(at_name,'(i10)') i
      rdos_file='rdos.'//trim(adjustl(at_name))//'.dat'
      open(12,file=trim(adjustl(rdos_file)))
      write(12,'(a1,2(1x,es22.15e2))') '#', real(i,dp), rpotential(i)/r_count(i)
      do j=1,size_dos
        E=dos_E1+(j-1)*dos_dE
        write(12,'(es22.15e2,1x,es22.15e2)') E, rdos(j,i)/r_count(i)
      end do
      close(12)
    end if
  end do
  deallocate(rdos)
  deallocate(ldos)

  deallocate(r_count)
  deallocate(r_bins)
  deallocate(at_dist)
  deallocate(at_poss)
  do i=1,num_ldos
    deallocate(at_sym(i)%list)
  end do
  deallocate(at_sym)
  deallocate(ldos_at)

end program adatom_post
