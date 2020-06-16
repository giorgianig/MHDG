MODULE matrices_types

   TYPE MAT_CSR_TYP ! A type to store matrices in CSR format
      logical                        :: start  ! keeps track if it is the first solve
      integer :: n ! Number of (locally owned) columns in the matrix
      integer :: nnz ! Number of (locally owned) non-zeros
      integer, dimension(:), pointer :: rowptr => null() ! Index of first element of each row in cols and vals
      integer, dimension(:), pointer :: cols => null() ! Column of each element
      real(8), dimension(:), pointer :: vals => null()  ! Value of each element
      integer, dimension(:), pointer :: loc2glob => null() ! Global column number of the local columns
   END TYPE

   TYPE RHS_TYP
      integer :: n ! Number of (locally owned) columns in the matrix
      real(8), dimension(:), pointer :: vals => null()  ! Value of each element
      integer, dimension(:), pointer :: loc2glob => null() ! Global column number of the local columns
   END TYPE

   TYPE MAT_IJV_TYP ! A type to store matrices in IJV format
      integer :: n ! Number of (locally owned) columns in the matrix
      integer :: nnz ! Number of (locally owned) non-zeros
      integer, dimension(:), pointer :: rows => null() ! Row of each element
      integer, dimension(:), pointer :: cols => null() ! Column of each element
      real(8), dimension(:), pointer :: vals => null()  ! Value of each element
   END TYPE

END MODULE matrices_types
