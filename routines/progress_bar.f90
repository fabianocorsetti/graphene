subroutine progress_bar(n1,n2,n)
  implicit none

  integer, parameter :: l=50
  integer, parameter :: unit_num=6
  integer, parameter :: dp=selected_real_kind(15,300)

  integer, intent(in) :: n1, n2, n

  integer :: i, p
  integer, save :: marker=0

  if (n==n1) then
    marker=0
    write(unit_num,'(a3)',advance='no') '#0%'
    do i=1+2,l/2-1
      write(unit_num,'(a1)',advance='no') '='
    end do
    write(unit_num,'(a3)',advance='no') '50%'
    do i=l/2+3,l-4
      write(unit_num,'(a1)',advance='no') '='
    end do
    write(unit_num,'(a5)',advance='yes') '100%#'
    write(unit_num,'(a1)',advance='no') '#'
    flush(6)
  end if
  p=nint(real(l,dp)*real(n-n1+1,dp)/(real(n2-n1+1,dp)))
  if (marker<p) then
    do i=1,p-marker
      write(unit_num,'(a1)',advance='no') '.'
      flush(6)
    end do
    marker=p
  end if
  if (n==n2) then
    write(unit_num,'(a1)',advance='yes') '#'
    marker=0
  end if

end subroutine progress_bar
