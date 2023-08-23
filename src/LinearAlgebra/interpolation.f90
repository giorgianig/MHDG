MODULE interpolation
  USE printutils
contains

  function binarysearch(length, array, value, delta)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "delta" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif

    implicit none
    integer, intent(in) :: length
    real, dimension(length), intent(in) :: array
    !f2py depend(length) array
    real, intent(in) :: value
    real, intent(in), optional :: delta

    integer :: binarysearch

    integer :: left, middle, right
    real :: d

    if (present(delta) .eqv. .true.) then
      d = delta
    else
      d = 1e-9
    endif

    left = 1
    right = length
    do
      if (left > right) then
        exit
      endif
      middle = nint((left + right)/2.0)
      if (abs(array(middle) - value) <= d) then
        binarySearch = middle
        return
      else if (array(middle) > value) then
        right = middle - 1
      else
        left = middle + 1
      end if
    end do
    binarysearch = right

  end function binarysearch

  real function interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
    implicit none
    integer, intent(in) :: x_len, y_len
    real, dimension(x_len), intent(in) :: x_array
    real, dimension(y_len), intent(in) :: y_array
    real, dimension(x_len, y_len), intent(in) :: f
    real, intent(in) :: x, y
    real, intent(in), optional :: delta
    !f2py depend(x_len) x_array, f
    !f2py depend(y_len) y_array, f

    real :: denom, x1, x2, y1, y2
    integer :: i, j

    i = binarysearch(x_len, x_array, x, delta)
    j = binarysearch(y_len, y_array, y, delta)

    IF (i == x_len) THEN
      WRITE (6, *) "Problem in the binary search"
      WRITE (6, *) "x", x
      WRITE (6, *) "max(x_array)", maxval(x_array)
      WRITE (6, *) "min(x_array)", minval(x_array)
    END IF

    IF (j == y_len) THEN
      WRITE (6, *) "Problem in the binary search"
      WRITE (6, *) "y", y
      WRITE (6, *) "max(y_array)", maxval(y_array)
      WRITE (6, *) "min(y_array)", minval(y_array)
    END IF

    x1 = x_array(i)
    x2 = x_array(i + 1)

    y1 = y_array(j)
    y2 = y_array(j + 1)

    denom = (x2 - x1)*(y2 - y1)

    interpolate = (f(i, j)*(x2 - x)*(y2 - y) + f(i + 1, j)*(x - x1)*(y2 - y) + &
      f(i, j + 1)*(x2 - x)*(y - y1) + f(i + 1, j + 1)*(x - x1)*(y - y1))/denom

  end function interpolate
    
   function nodesearch(x, y, xy_len, x_array, y_array)
    ! Given a  point (x,y), returns the index of the closest node in 2D
    implicit none
    integer,intent(in) :: xy_len
    real,intent(in) :: x, y
    real,dimension(xy_len),intent(in):: x_array, y_array
    real*8 :: d(xy_len)
    
    integer :: nodesearch
   
    d = sqrt( (x_array - x)**2 + (y_array - y)**2 )
    nodesearch = minloc(d, DIM = 1)
 
  end function nodesearch
  
  SUBROUTINE lineintegration(qp_len, x_vec, y_vec, f, Xc, T, nodes2D, nli) 
    ! Given a set of points (x, y) along a line of sight, retruns the line integration of
    ! f. The function f is evaluated in the closest nodes to points (x, y).
    implicit none
    integer,intent(in) :: qp_len, nodes2D, T(:,:)
    real*8,intent(in) :: x_vec(:), y_vec(:), f(:), Xc(:,:)
    real*8,intent(out) :: nli
    integer :: i, iel, inp, n_ind, f_ind(2), np_len
    real*8 :: x, y, dl(qp_len-1), x_qp(qp_len), f_qp(qp_len)
   
    ! Search closest node and evaluate f
    DO i = 1, qp_len
       x = x_vec(i)
       y = y_vec(i)
       np_len = size(Xc,1)
       n_ind = nodesearch(x, y, np_len, Xc(:,1), Xc(:,2))
       f_ind = findloc(T,n_ind)
       iel = f_ind(1)
       inp = f_ind(2)
       x_qp(i) = Xc(n_ind,1)
       f_qp(i) = f((iel-1)*nodes2D + inp)
    END DO

    dl = (x_qp(2:qp_len) - x_qp(1:qp_len - 1))*1.901e-3
    nli = sum(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len - 1))*dl)
     
 END SUBROUTINE lineintegration

  !real function lineintegration(qp_len, x_vec, y_vec, np_len, x_array, y_array ,f) 
    ! Given a set of points (x, y) along a line of sight, retruns the line integration of
    ! f. The function f is evaluated in the closest nodes to points (x, y).
  !  implicit none
  !  real,intent(in) :: qp_len, np_len
  !  real,dimension(qp_len),intent(in) :: x_vec,  y_vec
  !  real,dimension(np_len),intent(in) :: x_array,  y_array
  !  real, intent(in) :: f
  !  integer :: i, iel, inp, n_ind, f_ind(2)
  !  real*8 :: x, y, dl(qp_len-1), x_qp(qp_len), f_qp(qp_len)
   
    ! Search closest node and evaluate f
  !  DO i = 1, qp_len
  !     x = x_vec(i)
  !     y = y_vec(i)
  !     n_ind = nodesearch(x, y, np_len, x_array, y_array)
  !     f_ind = findloc(Mesh%T(n_ind))
  !     iel = f_ind(1)
  !     ip = f_ind(2)
  !     x_qp(i) = x_array(n_ind)
  !     f_qp(i) = f((iel-1)*refElPol%N2D + ip)
  !  END DO

    ! Integrate
  !  dl = x_qp(2:qp_len) - x_qp(1:qp_len-1)
  !  lineintegration = sum(0.5*(f_qp(2:qp_len) + f_qp(1:qp_len-1))*dl)
     
 !end function lineintegration
 
END MODULE interpolation
