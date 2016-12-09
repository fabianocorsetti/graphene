module fermi_dirac_stats
implicit none

private

integer, parameter :: dp=selected_real_kind(15,300)
integer, parameter :: max_iter=10000

real(dp), parameter :: err_tol=1.0d-10
real(dp), parameter :: gr=0.5_dp*(sqrt(5.0_dp)-1.0_dp)

interface FD_occ_tot_eigs
  module procedure occ_tot_eigs
  module procedure occ_tot_eigs_spin
  module procedure occ_tot_eigs_kpoints
  module procedure occ_tot_eigs_spin_kpoints
end interface FD_occ_tot_eigs

interface FD_occ_tot_DOS
  module procedure occ_tot_DOS
  module procedure occ_tot_DOS_spin
  module procedure occ_tot_DOS_kpoints
  module procedure occ_tot_DOS_spin_kpoints
end interface FD_occ_tot_DOS

interface FD_calc_E_F_eigs
  module procedure calc_E_F_eigs
  module procedure calc_E_F_eigs_spin
  module procedure calc_E_F_eigs_kpoints
  module procedure calc_E_F_eigs_spin_kpoints
end interface FD_calc_E_F_eigs

interface FD_calc_E_F_DOS
  module procedure calc_E_F_DOS
  module procedure calc_E_F_DOS_spin
  module procedure calc_E_F_DOS_kpoints
  module procedure calc_E_F_DOS_spin_kpoints
end interface FD_calc_E_F_DOS

public :: FD_occ_tot_eigs
public :: FD_occ_tot_DOS
public :: FD_calc_E_F_eigs
public :: FD_calc_E_F_DOS

contains

real(dp) function occ_tot_eigs(T,E_F,n,eigs,occs)
  implicit none

  integer, intent(in) :: n

  real(dp), intent(in) :: T, E_F, eigs(n)
  real(dp), intent(out), optional :: occs(n)

  integer :: i

  real(dp) :: beta, occ_temp

  beta=1.0_dp/T

  occ_tot_eigs=0.0_dp
  if (present(occs)) then
    do i=1,n
      occs(i)=1.0_dp/(exp(beta*(eigs(i)-E_F))+1.0_dp)
      occ_tot_eigs=occ_tot_eigs+occs(i)
    end do
  else
    do i=1,n
      occ_temp=1.0_dp/(exp(beta*(eigs(i)-E_F))+1.0_dp)
      occ_tot_eigs=occ_tot_eigs+occ_temp
    end do
  end if

end function occ_tot_eigs

real(dp) function occ_tot_eigs_spin(T,E_F,n,eigs,occs)
  implicit none

  integer, intent(in) :: n

  real(dp), intent(in) :: T, E_F, eigs(n,2)
  real(dp), intent(out), optional :: occs(n,2)

  integer :: i, j

  real(dp) :: beta, occ_temp

  beta=1.0_dp/T

  occ_tot_eigs_spin=0.0_dp
  if (present(occs)) then
    do j=1,2
      do i=1,n
        occs(i,j)=1.0_dp/(exp(beta*(eigs(i,j)-E_F))+1.0_dp)
        occ_tot_eigs_spin=occ_tot_eigs_spin+occs(i,j)
      end do
    end do
  else
    do j=1,2
      do i=1,n
        occ_temp=1.0_dp/(exp(beta*(eigs(i,j)-E_F))+1.0_dp)
        occ_tot_eigs_spin=occ_tot_eigs_spin+occ_temp
      end do
    end do
  end if

end function occ_tot_eigs_spin

real(dp) function occ_tot_eigs_kpoints(T,E_F,n,nk,eigs,weights,occs)
  implicit none

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, E_F, eigs(n,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out), optional :: occs(n,nk)

  integer :: i, k

  real(dp) :: beta, weight_temp, occ_temp

  beta=1.0_dp/T
  if (.not. present(weights)) weight_temp=1.0_dp/real(nk,dp)

  occ_tot_eigs_kpoints=0.0_dp
  if (present(occs)) then
    if (present(weights)) then
      do k=1,nk
        do i=1,n
          occs(i,k)=1.0_dp/(exp(beta*(eigs(i,k)-E_F))+1.0_dp)
          occ_tot_eigs_kpoints=occ_tot_eigs_kpoints+weights(k)*occs(i,k)
        end do
      end do
    else
      do k=1,nk
        do i=1,n
          occs(i,k)=1.0_dp/(exp(beta*(eigs(i,k)-E_F))+1.0_dp)
          occ_tot_eigs_kpoints=occ_tot_eigs_kpoints+weight_temp*occs(i,k)
        end do
      end do
    end if
  else
    if (present(weights)) then
      do k=1,nk
        do i=1,n
          occ_temp=1.0_dp/(exp(beta*(eigs(i,k)-E_F))+1.0_dp)
          occ_tot_eigs_kpoints=occ_tot_eigs_kpoints+weights(k)*occ_temp
        end do
      end do
    else
      do k=1,nk
        do i=1,n
          occ_temp=1.0_dp/(exp(beta*(eigs(i,k)-E_F))+1.0_dp)
          occ_tot_eigs_kpoints=occ_tot_eigs_kpoints+weight_temp*occ_temp
        end do
      end do
    end if
  end if

end function occ_tot_eigs_kpoints

real(dp) function occ_tot_eigs_spin_kpoints(T,E_F,n,nk,eigs,weights,occs)
  implicit none

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, E_F, eigs(n,2,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out), optional :: occs(n,2,nk)

  integer :: i, j, k

  real(dp) :: beta, weight_temp, occ_temp

  beta=1.0_dp/T
  if (.not. present(weights)) weight_temp=1.0_dp/real(nk,dp)

  occ_tot_eigs_spin_kpoints=0.0_dp
  if (present(occs)) then
    if (present(weights)) then
      do k=1,nk
        do j=1,2
          do i=1,n
            occs(i,j,k)=1.0_dp/(exp(beta*(eigs(i,j,k)-E_F))+1.0_dp)
            occ_tot_eigs_spin_kpoints=occ_tot_eigs_spin_kpoints+weights(k)*occs(i,j,k)
          end do
        end do
      end do
    else
      do k=1,nk
        do j=1,2
          do i=1,n
            occs(i,j,k)=1.0_dp/(exp(beta*(eigs(i,j,k)-E_F))+1.0_dp)
            occ_tot_eigs_spin_kpoints=occ_tot_eigs_spin_kpoints+weight_temp*occs(i,j,k)
          end do
        end do
      end do
    end if
  else
    if (present(weights)) then
      do k=1,nk
        do j=1,2
          do i=1,n
            occ_temp=1.0_dp/(exp(beta*(eigs(i,j,k)-E_F))+1.0_dp)
            occ_tot_eigs_spin_kpoints=occ_tot_eigs_spin_kpoints+weights(k)*occ_temp
          end do
        end do
      end do
    else
      do k=1,nk
        do j=1,2
          do i=1,n
            occ_temp=1.0_dp/(exp(beta*(eigs(i,j,k)-E_F))+1.0_dp)
            occ_tot_eigs_spin_kpoints=occ_tot_eigs_spin_kpoints+weight_temp*occ_temp
          end do
        end do
      end do
    end if
  end if

end function occ_tot_eigs_spin_kpoints

real(dp) function occ_tot_DOS(T,E_F,n,DOS,occs)
  implicit none

  integer, intent(in) :: n

  real(dp), intent(in) :: T, E_F, DOS(2,n)
  real(dp), intent(out), optional :: occs(n)

  integer :: i

  real(dp) :: beta, dE, occ_temp

  beta=1.0_dp/T
  dE=DOS(1,2)-DOS(1,1)

  occ_tot_DOS=0.0_dp
  if (present(occs)) then
    do i=1,n
      occs(i)=DOS(2,i)/(exp(beta*(DOS(1,i)-E_F))+1.0_dp)
      occ_tot_DOS=occ_tot_DOS+occs(i)
    end do
  else
    do i=1,n
      occ_temp=DOS(2,i)/(exp(beta*(DOS(1,i)-E_F))+1.0_dp)
      occ_tot_DOS=occ_tot_DOS+occ_temp
    end do
  end if

  occ_tot_DOS=dE*occ_tot_DOS

end function occ_tot_DOS

real(dp) function occ_tot_DOS_spin(T,E_F,n,DOS,occs)
  implicit none

  integer, intent(in) :: n

  real(dp), intent(in) :: T, E_F, DOS(2,n,2)
  real(dp), intent(out), optional :: occs(n,2)

  integer :: i, j

  real(dp) :: beta, dE, occ_temp

  beta=1.0_dp/T
  dE=DOS(1,2,1)-DOS(1,1,1)

  occ_tot_DOS_spin=0.0_dp
  if (present(occs)) then
    do j=1,2
      do i=1,n
        occs(i,j)=DOS(2,i,j)/(exp(beta*(DOS(1,i,j)-E_F))+1.0_dp)
        occ_tot_DOS_spin=occ_tot_DOS_spin+occs(i,j)
      end do
    end do
  else
    do j=1,2
      do i=1,n
        occ_temp=DOS(2,i,j)/(exp(beta*(DOS(1,i,j)-E_F))+1.0_dp)
        occ_tot_DOS_spin=occ_tot_DOS_spin+occ_temp
      end do
    end do
  end if

  occ_tot_DOS_spin=dE*occ_tot_DOS_spin

end function occ_tot_DOS_spin

real(dp) function occ_tot_DOS_kpoints(T,E_F,n,nk,DOS,weights,occs)
  implicit none

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, E_F, DOS(2,n,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out), optional :: occs(n,nk)

  integer :: i, k

  real(dp) :: beta, dE, weight_temp, occ_temp

  beta=1.0_dp/T
  if (.not. present(weights)) weight_temp=1.0_dp/real(nk,dp)
  dE=DOS(1,2,1)-DOS(1,1,1)

  occ_tot_DOS_kpoints=0.0_dp
  if (present(occs)) then
    if (present(weights)) then
      do k=1,nk
        do i=1,n
          occs(i,k)=DOS(2,i,k)/(exp(beta*(DOS(1,i,k)-E_F))+1.0_dp)
          occ_tot_DOS_kpoints=occ_tot_DOS_kpoints+weights(k)*occs(i,k)
        end do
      end do
    else
      do k=1,nk
        do i=1,n
          occs(i,k)=DOS(2,i,k)/(exp(beta*(DOS(1,i,k)-E_F))+1.0_dp)
          occ_tot_DOS_kpoints=occ_tot_DOS_kpoints+weight_temp*occs(i,k)
        end do
      end do
    end if
  else
    if (present(weights)) then
      do k=1,nk
        do i=1,n
          occ_temp=DOS(2,i,k)/(exp(beta*(DOS(1,i,k)-E_F))+1.0_dp)
          occ_tot_DOS_kpoints=occ_tot_DOS_kpoints+weights(k)*occ_temp
        end do
      end do
    else
      do k=1,nk
        do i=1,n
          occ_temp=DOS(2,i,k)/(exp(beta*(DOS(1,i,k)-E_F))+1.0_dp)
          occ_tot_DOS_kpoints=occ_tot_DOS_kpoints+weight_temp*occ_temp
        end do
      end do
    end if
  end if

  occ_tot_DOS_kpoints=dE*occ_tot_DOS_kpoints

end function occ_tot_DOS_kpoints

real(dp) function occ_tot_DOS_spin_kpoints(T,E_F,n,nk,DOS,weights,occs)
  implicit none

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, E_F, DOS(2,n,2,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out), optional :: occs(n,2,nk)

  integer :: i, j, k

  real(dp) :: beta, dE, weight_temp, occ_temp

  beta=1.0_dp/T
  if (.not. present(weights)) weight_temp=1.0_dp/real(nk,dp)
  dE=DOS(1,2,1,1)-DOS(1,1,1,1)

  occ_tot_DOS_spin_kpoints=0.0_dp
  if (present(occs)) then
    if (present(weights)) then
      do k=1,nk
        do j=1,2
          do i=1,n
            occs(i,j,k)=DOS(2,i,j,k)/(exp(beta*(DOS(1,i,j,k)-E_F))+1.0_dp)
            occ_tot_DOS_spin_kpoints=occ_tot_DOS_spin_kpoints+weights(k)*occs(i,j,k)
          end do
        end do
      end do
    else
      do k=1,nk
        do j=1,2
          do i=1,n
            occs(i,j,k)=DOS(2,i,j,k)/(exp(beta*(DOS(1,i,j,k)-E_F))+1.0_dp)
            occ_tot_DOS_spin_kpoints=occ_tot_DOS_spin_kpoints+weight_temp*occs(i,j,k)
          end do
        end do
      end do
    end if
  else
    if (present(weights)) then
      do k=1,nk
        do j=1,2
          do i=1,n
            occ_temp=DOS(2,i,j,k)/(exp(beta*(DOS(1,i,j,k)-E_F))+1.0_dp)
            occ_tot_DOS_spin_kpoints=occ_tot_DOS_spin_kpoints+weights(k)*occ_temp
          end do
        end do
      end do
    else
      do k=1,nk
        do j=1,2
          do i=1,n
            occ_temp=DOS(2,i,j,k)/(exp(beta*(DOS(1,i,j,k)-E_F))+1.0_dp)
            occ_tot_DOS_spin_kpoints=occ_tot_DOS_spin_kpoints+weight_temp*occ_temp
          end do
        end do
      end do
    end if
  end if

  occ_tot_DOS_spin_kpoints=dE*occ_tot_DOS_spin_kpoints

end function occ_tot_DOS_spin_kpoints

subroutine calc_E_F_eigs(T,E_F,guess_E_F,occ_tot,occ_err,n,eigs,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n

  real(dp), intent(in) :: T, occ_tot, eigs(n)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  if (guess_E_F) then
    xc=E_F
  else
    i=nint(occ_tot)
    if (i==0) then
      xc=eigs(1)
    else if (i==n) then
      xc=eigs(n)
    else
      xc=0.5_dp*(eigs(i)+eigs(i+1))
    end if
  end if
  yc=(FD_occ_tot_eigs(T,xc,n,eigs)-occ_tot)**2

  xb=xc+1.0_dp
  yb=(FD_occ_tot_eigs(T,xb,n,eigs)-occ_tot)**2

  xa=(xc+xb*(gr-1.0_dp))/gr
  ya=(FD_occ_tot_eigs(T,xa,n,eigs)-occ_tot)**2

  xd=xa+gr*(xb-xa)
  yd=(FD_occ_tot_eigs(T,xd,n,eigs)-occ_tot)**2

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      yc=(FD_occ_tot_eigs(T,xc,n,eigs)-occ_tot)**2
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      yd=(FD_occ_tot_eigs(T,xd,n,eigs)-occ_tot)**2
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,eigs,occs)
  else
    occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,eigs)
  end if

end subroutine calc_E_F_eigs

subroutine calc_E_F_eigs_spin(T,E_F,guess_E_F,occ_tot,occ_err,n,eigs,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n

  real(dp), intent(in) :: T, occ_tot, eigs(n,2)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,2)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  if (guess_E_F) then
    xc=E_F
  else
    i=nint(occ_tot)
    if (i==0) then
      xc=eigs(1,1)
    else if (i==n) then
      xc=eigs(n,1)
    else
      xc=0.5_dp*(eigs(i,1)+eigs(i+1,1))
    end if
  end if
  yc=(FD_occ_tot_eigs(T,xc,n,eigs)-occ_tot)**2

  xb=xc+1.0_dp
  yb=(FD_occ_tot_eigs(T,xb,n,eigs)-occ_tot)**2

  xa=(xc+xb*(gr-1.0_dp))/gr
  ya=(FD_occ_tot_eigs(T,xa,n,eigs)-occ_tot)**2

  xd=xa+gr*(xb-xa)
  yd=(FD_occ_tot_eigs(T,xd,n,eigs)-occ_tot)**2

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      yc=(FD_occ_tot_eigs(T,xc,n,eigs)-occ_tot)**2
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      yd=(FD_occ_tot_eigs(T,xd,n,eigs)-occ_tot)**2
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,eigs,occs)
  else
    occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,eigs)
  end if

end subroutine calc_E_F_eigs_spin

subroutine calc_E_F_eigs_kpoints(T,E_F,guess_E_F,occ_tot,occ_err,n,nk,eigs,weights,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, occ_tot, eigs(n,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,nk)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  if (guess_E_F) then
    xc=E_F
  else
    i=nint(occ_tot)
    if (i==0) then
      xc=eigs(1,1)
    else if (i==n) then
      xc=eigs(n,1)
    else
      xc=0.5_dp*(eigs(i,1)+eigs(i+1,1))
    end if
  end if
  if (present(weights)) then
    yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs,weights)-occ_tot)**2
  else
    yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs)-occ_tot)**2
  end if

  xb=xc+1.0_dp
  if (present(weights)) then
    yb=(FD_occ_tot_eigs(T,xb,n,nk,eigs,weights)-occ_tot)**2
  else
    yb=(FD_occ_tot_eigs(T,xb,n,nk,eigs)-occ_tot)**2
  end if

  xa=(xc+xb*(gr-1.0_dp))/gr
  if (present(weights)) then
    ya=(FD_occ_tot_eigs(T,xa,n,nk,eigs,weights)-occ_tot)**2
  else
    ya=(FD_occ_tot_eigs(T,xa,n,nk,eigs)-occ_tot)**2
  end if

  xd=xa+gr*(xb-xa)
  if (present(weights)) then
    yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs,weights)-occ_tot)**2
  else
    yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs)-occ_tot)**2
  end if

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      if (present(weights)) then
        yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs,weights)-occ_tot)**2
      else
        yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs)-occ_tot)**2
      end if
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      if (present(weights)) then
        yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs,weights)-occ_tot)**2
      else
        yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs)-occ_tot)**2
      end if
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,weights,occs)
    else
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,occs=occs)
    end if
  else
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,weights)
    else
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs)
    end if
  end if

end subroutine calc_E_F_eigs_kpoints

subroutine calc_E_F_eigs_spin_kpoints(T,E_F,guess_E_F,occ_tot,occ_err,n,nk,eigs,weights,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, occ_tot, eigs(n,2,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,2,nk)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  if (guess_E_F) then
    xc=E_F
  else
    i=nint(occ_tot)
    if (i==0) then
      xc=eigs(1,1,1)
    else if (i==n) then
      xc=eigs(n,1,1)
    else
      xc=0.5_dp*(eigs(i,1,1)+eigs(i+1,1,1))
    end if
  end if
  if (present(weights)) then
    yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs,weights)-occ_tot)**2
  else
    yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs)-occ_tot)**2
  end if

  xb=xc+1.0_dp
  if (present(weights)) then
    yb=(FD_occ_tot_eigs(T,xb,n,nk,eigs,weights)-occ_tot)**2
  else
    yb=(FD_occ_tot_eigs(T,xb,n,nk,eigs)-occ_tot)**2
  end if

  xa=(xc+xb*(gr-1.0_dp))/gr
  if (present(weights)) then
    ya=(FD_occ_tot_eigs(T,xa,n,nk,eigs,weights)-occ_tot)**2
  else
    ya=(FD_occ_tot_eigs(T,xa,n,nk,eigs)-occ_tot)**2
  end if

  xd=xa+gr*(xb-xa)
  if (present(weights)) then
    yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs,weights)-occ_tot)**2
  else
    yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs)-occ_tot)**2
  end if

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      if (present(weights)) then
        yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs,weights)-occ_tot)**2
      else
        yc=(FD_occ_tot_eigs(T,xc,n,nk,eigs)-occ_tot)**2
      end if
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      if (present(weights)) then
        yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs,weights)-occ_tot)**2
      else
        yd=(FD_occ_tot_eigs(T,xd,n,nk,eigs)-occ_tot)**2
      end if
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,weights,occs)
    else
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,occs=occs)
    end if
  else
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs,weights)
    else
      occ_err=occ_tot-FD_occ_tot_eigs(T,E_F,n,nk,eigs)
    end if
  end if

end subroutine calc_E_F_eigs_spin_kpoints

subroutine calc_E_F_DOS(T,E_F,guess_E_F,occ_tot,occ_err,n,eigs,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n

  real(dp), intent(in) :: T, occ_tot, eigs(2,n)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  !if ((occ_tot<0.0_dp) .or. &
  !    (occ_tot>real(n,dp))) then
  !  print*, "ERROR in calc_E_F!"
  !  stop
  !end if

  !if (guess_E_F) then
    xc=E_F
  !else
  !  i=nint(occ_tot)
  !  if (i==0) then
  !    xc=eigs(1)
  !  else if (i==n) then
  !    xc=eigs(n)
  !  else
  !    xc=0.5_dp*(eigs(i)+eigs(i+1))
  !  end if
  !end if
  yc=(FD_occ_tot_DOS(T,xc,n,eigs)-occ_tot)**2

  xb=xc+1.0_dp
  yb=(FD_occ_tot_DOS(T,xb,n,eigs)-occ_tot)**2

  xa=(xc+xb*(gr-1.0_dp))/gr
  ya=(FD_occ_tot_DOS(T,xa,n,eigs)-occ_tot)**2

  xd=xa+gr*(xb-xa)
  yd=(FD_occ_tot_DOS(T,xd,n,eigs)-occ_tot)**2

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      yc=(FD_occ_tot_DOS(T,xc,n,eigs)-occ_tot)**2
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      yd=(FD_occ_tot_DOS(T,xd,n,eigs)-occ_tot)**2
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,eigs,occs)
  else
    occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,eigs)
  end if

end subroutine calc_E_F_DOS

subroutine calc_E_F_DOS_spin(T,E_F,guess_E_F,occ_tot,occ_err,n,eigs,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n

  real(dp), intent(in) :: T, occ_tot, eigs(2,n,2)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,2)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  !if ((occ_tot<0.0_dp) .or. &
  !    (occ_tot>real(n,dp))) then
  !  print*, "ERROR in calc_E_F!"
  !  stop
  !end if

  !if (guess_E_F) then
    xc=E_F
  !else
  !  i=nint(occ_tot)
  !  if (i==0) then
  !    xc=eigs(1)
  !  else if (i==n) then
  !    xc=eigs(n)
  !  else
  !    xc=0.5_dp*(eigs(i)+eigs(i+1))
  !  end if
  !end if
  yc=(FD_occ_tot_DOS(T,xc,n,eigs)-occ_tot)**2

  xb=xc+1.0_dp
  yb=(FD_occ_tot_DOS(T,xb,n,eigs)-occ_tot)**2

  xa=(xc+xb*(gr-1.0_dp))/gr
  ya=(FD_occ_tot_DOS(T,xa,n,eigs)-occ_tot)**2

  xd=xa+gr*(xb-xa)
  yd=(FD_occ_tot_DOS(T,xd,n,eigs)-occ_tot)**2

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      yc=(FD_occ_tot_DOS(T,xc,n,eigs)-occ_tot)**2
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      yd=(FD_occ_tot_DOS(T,xd,n,eigs)-occ_tot)**2
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,eigs,occs)
  else
    occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,eigs)
  end if

end subroutine calc_E_F_DOS_spin

subroutine calc_E_F_DOS_kpoints(T,E_F,guess_E_F,occ_tot,occ_err,n,nk,dos,weights,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, occ_tot, dos(2,n,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,nk)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  !if (guess_E_F) then
    xc=E_F
  !else
  !  i=nint(occ_tot)
  !  if (i==0) then
  !    xc=eigs(1,1)
  !  else if (i==n) then
  !    xc=eigs(n,1)
  !  else
  !    xc=0.5_dp*(eigs(i,1)+eigs(i+1,1))
  !  end if
  !end if
  if (present(weights)) then
    yc=(FD_occ_tot_DOS(T,xc,n,nk,dos,weights)-occ_tot)**2
  else
    yc=(FD_occ_tot_DOS(T,xc,n,nk,dos)-occ_tot)**2
  end if

  xb=xc+1.0_dp
  if (present(weights)) then
    yb=(FD_occ_tot_DOS(T,xb,n,nk,dos,weights)-occ_tot)**2
  else
    yb=(FD_occ_tot_DOS(T,xb,n,nk,dos)-occ_tot)**2
  end if

  xa=(xc+xb*(gr-1.0_dp))/gr
  if (present(weights)) then
    ya=(FD_occ_tot_DOS(T,xa,n,nk,dos,weights)-occ_tot)**2
  else
    ya=(FD_occ_tot_DOS(T,xa,n,nk,dos)-occ_tot)**2
  end if

  xd=xa+gr*(xb-xa)
  if (present(weights)) then
    yd=(FD_occ_tot_DOS(T,xd,n,nk,dos,weights)-occ_tot)**2
  else
    yd=(FD_occ_tot_DOS(T,xd,n,nk,dos)-occ_tot)**2
  end if

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      if (present(weights)) then
        yc=(FD_occ_tot_DOS(T,xc,n,nk,dos,weights)-occ_tot)**2
      else
        yc=(FD_occ_tot_DOS(T,xc,n,nk,dos)-occ_tot)**2
      end if
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      if (present(weights)) then
        yd=(FD_occ_tot_DOS(T,xd,n,nk,dos,weights)-occ_tot)**2
      else
        yd=(FD_occ_tot_DOS(T,xd,n,nk,dos)-occ_tot)**2
      end if
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,weights,occs)
    else
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,occs=occs)
    end if
  else
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,weights)
    else
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos)
    end if
  end if

end subroutine calc_E_F_DOS_kpoints

subroutine calc_E_F_DOS_spin_kpoints(T,E_F,guess_E_F,occ_tot,occ_err,n,nk,dos,weights,occs)
  implicit none

  logical, intent(in) :: guess_E_F

  integer, intent(in) :: n, nk

  real(dp), intent(in) :: T, occ_tot, dos(2,n,2,nk)
  real(dp), intent(in), optional :: weights(nk)
  real(dp), intent(out) :: occ_err
  real(dp), intent(out), optional :: occs(n,2,nk)
  real(dp), intent(inout) :: E_F

  logical :: conv

  integer :: i

  real(dp) :: xa, xb, xc, xd
  real(dp) :: ya, yb, yc, yd

  if ((occ_tot<0.0_dp) .or. &
      (occ_tot>real(n,dp))) then
    print*, "ERROR in calc_E_F!"
    stop
  end if

  !if (guess_E_F) then
    xc=E_F
  !else
  !  i=nint(occ_tot)
  !  if (i==0) then
  !    xc=eigs(1,1,1)
  !  else if (i==n) then
  !    xc=eigs(n,1,1)
  !  else
  !    xc=0.5_dp*(eigs(i,1,1)+eigs(i+1,1,1))
  !  end if
  !end if
  if (present(weights)) then
    yc=(FD_occ_tot_DOS(T,xc,n,nk,dos,weights)-occ_tot)**2
  else
    yc=(FD_occ_tot_DOS(T,xc,n,nk,dos)-occ_tot)**2
  end if

  xb=xc+1.0_dp
  if (present(weights)) then
    yb=(FD_occ_tot_DOS(T,xb,n,nk,dos,weights)-occ_tot)**2
  else
    yb=(FD_occ_tot_DOS(T,xb,n,nk,dos)-occ_tot)**2
  end if

  xa=(xc+xb*(gr-1.0_dp))/gr
  if (present(weights)) then
    ya=(FD_occ_tot_DOS(T,xa,n,nk,dos,weights)-occ_tot)**2
  else
    ya=(FD_occ_tot_DOS(T,xa,n,nk,dos)-occ_tot)**2
  end if

  xd=xa+gr*(xb-xa)
  if (present(weights)) then
    yd=(FD_occ_tot_DOS(T,xd,n,nk,dos,weights)-occ_tot)**2
  else
    yd=(FD_occ_tot_DOS(T,xd,n,nk,dos)-occ_tot)**2
  end if

  conv=.false.
  do i=1,max_iter
    if (yc<yd) then
      E_F=xc
      if (yc<err_tol) then
        conv=.true.
        exit
      end if
      xb=xd
      yb=yd
      xd=xc
      yd=yc
      xc=xb-gr*(xb-xa)
      if (present(weights)) then
        yc=(FD_occ_tot_DOS(T,xc,n,nk,dos,weights)-occ_tot)**2
      else
        yc=(FD_occ_tot_DOS(T,xc,n,nk,dos)-occ_tot)**2
      end if
    else
      E_F=xd
      if (yd<err_tol) then
        conv=.true.
        exit
      end if
      xa=xc
      ya=yc
      xc=xd
      yc=yd
      xd=xa+gr*(xb-xa)
      if (present(weights)) then
        yd=(FD_occ_tot_DOS(T,xd,n,nk,dos,weights)-occ_tot)**2
      else
        yd=(FD_occ_tot_DOS(T,xd,n,nk,dos)-occ_tot)**2
      end if
    end if
  end do

  if (.not. conv) then
    print*, "WARNING: calc_E_F not converged!"
  end if

  if (present(occs)) then
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,weights,occs)
    else
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,occs=occs)
    end if
  else
    if (present(weights)) then
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos,weights)
    else
      occ_err=occ_tot-FD_occ_tot_DOS(T,E_F,n,nk,dos)
    end if
  end if

end subroutine calc_E_F_DOS_spin_kpoints

end module fermi_dirac_stats
