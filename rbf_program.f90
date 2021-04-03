!*> =================================================================
!*> Radial Basis Functions (RBF) are frequently used when it is needed to construct an
!*> approximation of a certain multivariable function by interpolation between the
!*> existing data. This section will give a very brief description of the method to
!*> make the reader familiar with main principles and to show how it can be effectively
!*> combined with previously described POD theory in the structural context here of
!*> interest.
!*>
!*> For more detailed: Buhmann, M.D.: Radial Basis Functions.
!*>                    Cambridge University Press, Cambridge (2003).
!*>
!*> This source code uses the ludcmp and lubksb routines from the Numerical Recipes in Fortran 90, 
!*> Second Edition (1996) free, in our bookreader format.
!*> http://numerical.recipes/oldverswitcher.html
!*>
!*> Author: Bruno R. Jaccoud
!*> Version: 2.0
!*> Date: 2021-03-22
!*>
!*> =================================================================
!*
!*> =================================================================
!*>                     MODULE FOR GLOBAL VARIABLES
!*> =================================================================
module var_globais
  !* ... VARIABLES FOR COMPUTATIONAL TIME ...
  real(8) time_total,time_begin,time_end
  !* ... PARAMETERS AND VARIABLES TO RBF ...
  integer :: nsamples, npoints, ntraining
  integer, parameter :: ndim = 3
  real(8) :: shape_parameter
  real(8),allocatable,dimension(:,:) :: trainingMeas, trainingPoints, trainingPointsNormalize
  real(8),allocatable,dimension(:,:) :: evaluationPoints, evaluationFunction
  real(8),allocatable,dimension(:,:) :: basisFunctionMatrix, basisFunctionMatrixINV, interpolationCoefficientsMatrix
  real(8),allocatable,dimension(:,:) :: newInterpCoeffFun, appFunctionRBF
  real(8),allocatable,dimension(:) :: minP, maxP
end module var_globais
!*
!*> =================================================================
!*>                          MAIN PRINCIPAL
!*> =================================================================
program rbf_program
  use var_globais
  implicit none
  integer :: i, j
  real(8) :: func_y1, func_y2, func_y3

  ! Reading variables for training points matrix
  open(999, file = 'points.dat', status = 'old')
  read(999,*) npoints, nsamples

  ! Allocate variables
  call allocation
  
  do i = 1, npoints
    ! Training Points Matrix
    read(999,*) (trainingPoints(i,j), j = 1, nsamples)
  enddo
  ! Close training file
  close(999)
  
  ! Building the solutions matrix with a training points
  do j = 1, nsamples
    trainingMeas(1,j) = func_y1(trainingPoints(1,j))
    trainingMeas(2,j) = func_y2(trainingPoints(2,j))
    trainingMeas(3,j) = func_y3(trainingPoints(1,j),trainingPoints(2,j))
  enddo

  ! Constructing the Basis Function Matrix
  shape_parameter = 1.d2/dble(nsamples)
  call construct_rbf(npoints,nsamples,shape_parameter,trainingPoints,trainingPointsNormalize,minP,maxP,basisFunctionMatrix)

  ! Determinate a inverve of Basis Function Matrix
  basisFunctionMatrixINV = 0.d0  
  call pseudo_inverse(nsamples,nsamples,basisFunctionMatrix,basisFunctionMatrixINV)

  ! Determine a Interpolation Coefficients Matrix
  interpolationCoefficientsMatrix = matmul(trainingMeas,basisFunctionMatrixINV)

  ! Reading new points for evaluation
  open(998, file = 'exacts.dat', status = 'old')
  read(998,*) npoints !,ndim
  do i = 1, npoints
    ! Evaluation Points Matrix
    read(998,*) (evaluationPoints(i,j), j = 1, ndim)
  enddo
  close(998)

  do j = 1, ndim
    call forward_problem(npoints,nsamples,shape_parameter,evaluationPoints(:,j), &
        trainingPointsNormalize,minP,maxP,newInterpCoeffFun(:,j))
  enddo

  appFunctionRBF = matmul(interpolationCoefficientsMatrix,newInterpCoeffFun)

  ! Printing RBF Answer
  call printingRBF()

  ! ! Evaluation function Matrix
  ! real(8) :: evaluationFunction(ndim,ndim)
  ! do j = 1, ndim
  !   evaluationFunction(1,j) = func_y1(evaluationPoints(1,j))
  !   evaluationFunction(2,j) = func_y2(evaluationPoints(2,j))
  !   evaluationFunction(3,j) = func_y3(evaluationPoints(1,j),evaluationPoints(2,j))
  ! enddo

  ! real(8)
  ! open(unit = 1, file = "comparisonAnalyticalRBF.dat")
  ! evaluationPoints(npoints,ndim), evaluationFunction(ndim,ndim), appFunctionRBF(ndim,ndim), error

end program rbf_program
!*
!*> =================================================================
!*>                   ALLOCATE DEPENDENT VARIABLES
!*> =================================================================
subroutine allocation
  use var_globais
  allocate(trainingPoints(npoints,nsamples), trainingPointsNormalize(npoints,nsamples), trainingMeas(ndim,nsamples))
  allocate(basisFunctionMatrix(nsamples,nsamples),basisFunctionMatrixINV(nsamples,nsamples))
  allocate(minP(npoints), maxP(npoints))
  allocate(interpolationCoefficientsMatrix(ndim,nsamples))
  allocate(evaluationPoints(npoints,ndim))
  allocate(newInterpCoeffFun(nsamples,ndim), appFunctionRBF(ndim,ndim))
end subroutine allocation
!*
!*> =================================================================
!*>                   DEALLOCATE DEPENDENT VARIABLES
!*> =================================================================
subroutine deallocation
  use var_globais
  deallocate(trainingPoints, trainingPointsNormalize, trainingMeas)
  deallocate(basisFunctionMatrix,basisFunctionMatrixINV)
  deallocate(minP, maxP)
  deallocate(interpolationCoefficientsMatrix)
  deallocate(evaluationPoints, evaluationFunction)
  deallocate(newInterpCoeffFun)
end subroutine deallocation
!*
!*> =================================================================
!*>                        PRINT OUTPUT RESULTS
!*> =================================================================
subroutine output(m,n,k,name_variable,x)
  implicit none
  integer :: i, j, k
  character(len=*):: name_variable 
  integer :: m, n
  real(8) :: x(m,n)
  !*  --------------
  open(unit = k, file = trim(name_variable)//".dat")
  !*  --------------
  do i = 1, m
    write(k,*) (x(i,j), j = 1, n)
  enddo
  close(k)

end subroutine output
!*
!*> =================================================================
!*>              CONSTRUCT A SET OF RADIAL BASIS FUNCTION
!*> =================================================================
subroutine construct_rbf(m,n,shapeP,P,x,minParameters,maxParameters,phi)
  implicit none
  !*------------------------------------------------
  integer :: i,j
  integer :: m,n
  real(8) :: shapeP
  real(8) :: x(m,n),phi(n,n),xi(m),xj(m)
  real(8) :: P(m,n),minParameters(m),maxParameters(m)
  !*------------------------------------------------
  !*
  do i=1,m
    do j=1,n
      x(i,j)=(P(i,j)-minParameters(i))/(maxParameters(i)-minParameters(i))
    enddo
  enddo
  !*
  x = P
  do i=1,n
    xi(:)=x(:,i)
    do j=1,n
      xj(:)=x(:,j)-xi(:)
      phi(i,j)=dsqrt(dot_product(xj,xj)) !Linearsplines
      !phi(i,j)=dsqrt(dot_product(xj,xj)+c*c) !Multiquadrics
      ! phi(i,j)=1.d0/(dsqrt(dot_product(xj,xj)+c*c)) !Inverse Multiquadrics
    enddo
  enddo
end subroutine construct_rbf
!*
!*> =================================================================
!*>                FORWARD RBF PROBLEM - COMPUTE VALUES
!*> =================================================================
subroutine forward_problem(rows,cols,shapeP,evalPoints,setPoints,minPoints,maxPoints,funcPoints)
  implicit none
  integer :: i
  integer :: rows, cols
  real(8) :: shapeP
  real(8) :: evalPoints(rows), setPoints(rows,cols), normPoints(rows)
  real(8) :: minPoints(rows), maxPoints(rows), funcPoints(cols)
  
  do i = 1, rows
    normPoints(i) = (evalPoints(i) - minPoints(i)) / (maxPoints(i) - minPoints(i))
  enddo
  normPoints = evalPoints

  call solve_rbf(rows,cols,shapeP,setPoints,normPoints,funcPoints)

end subroutine forward_problem
!*
!*> =================================================================
!*>     INTERPOLATE A FUNCTION USING A PRE-CONSTRUCTED RADIAL BASIS
!*>                        FUNCTION COEFFICIENTS
!*> =================================================================
subroutine solve_rbf(rows,cols,shapeP,setPoints,normPoints,funcPoints)
  implicit none
  integer :: j
  integer :: rows, cols
  real(8) :: f1, shapeP
  real(8), dimension(:) :: normPoints(rows)
  real(8), dimension(:) :: setPoints(rows,cols)
  real(8), dimension(:) :: funcPoints(cols)

    f1 = 0.d0
    do j = 1, cols
      f1 = dsqrt(dot_product(normPoints(:) - setPoints(:,j),normPoints(:) - setPoints(:,j))) !Linearsplines
      ! f1 = dsqrt(dot_product(xj,xj)+c*c) !Multiquadrics
      ! f1 = 1.d0 / (dsqrt(dot_product(xj,xj)+c*c)) !Inverse Multiquadrics
      funcPoints(j) = f1
    enddo
end subroutine solve_rbf
!*
!*> =================================================================
!*>                       PRINT RBF INFORMATIONS
!*>   Training Points, Training Measurements, Basis Function Matrix
!*>   and Interpolation Coefficients Matrix
!*> =================================================================
subroutine printingRBF()
  use var_globais
  implicit none
  integer :: i, j

  open(unit = 1000, file = "printingRBF.dat")

  write(1000,'(/,A)') " =============================================="
  write(1000,'(A,I4,A,I4)') " Training Points : i = ", npoints, ", j = ", nsamples
  write(1000,'(A)') " ==============================================="
  do i = 1, npoints
    write(1000,*) (trainingPoints(i,j), j = 1, nsamples)
  enddo

  write(1000,'(/,A)') " =============================================="
  write(1000,'(A,I4,A,I4)') " Solution with the Training Points : i = ", ndim, ", j = ", nsamples
  write(1000,'(A)') " ==============================================="
  do i = 1, ndim
    write(1000,*) (trainingMeas(i,j), j = 1, nsamples)
  enddo

  write(1000,'(/,A)') " =============================================="
  write(1000,'(A,I4,A,I4)') " Basis Function Matrix : i = ", nsamples, ", j = ", nsamples
  write(1000,'(A)') " RADIAL BASIS FUNCTION : ..."
  write(1000,'(A)') " ==============================================="
  do i = 1, nsamples
    write(1000,*) (basisFunctionMatrix(i,j), j = 1, nsamples)
  enddo

  write(1000,'(/,A)') " =============================================="
  write(1000,'(A,I4,A,I4)') " Interpolation Coefficients Matrix : i = ", ndim, ", j = ", nsamples
  write(1000,'(A)') " ==============================================="
  do i = 1, ndim
    write(1000,*) (interpolationCoefficientsMatrix(i,j), j = 1, nsamples)
  enddo

  write(1000,'(/,A)') " =============================================="
  write(1000,'(A,I4,A,I4)') " New Evaluating - ", ndim, " new points"
  write(1000,'(A)') " ==============================================="
  write(1000,'(A)') " Evaluating Points "
  do i = 1, npoints
    write(1000,*) (evaluationPoints(i,j), j = 1, ndim)
  enddo
  write(1000,'(/,A)') " * Function Evaluating Points "
  do i = 1, ndim
    write(1000,*) (evaluationFunction(i,j), j = 1, ndim)
  enddo
  write(1000,'(/,A)') " * Interpolation Coefficients Matrix "
  do i = 1, nsamples
    write(1000,*) (newInterpCoeffFun(i,j), j = 1, ndim)
  enddo
  write(1000,'(/,A)') " * Interpolated results Matrix "
  do i = 1, ndim
    write(1000,*) (appFunctionRBF(i,j), j = 1, ndim)
  enddo

  close(1000)
end subroutine printingRBF
!*
!*> =================================================================
!*  FUNCTIONS
!*> =================================================================
function func_y1(x)
  real(8) :: func_y1, x
  func_y1 = sqrt(x)
end function func_y1
function func_y2(x)
  real(8) :: func_y2, x
  func_y2 = x**2.d0
end function func_y2
function func_y3(x1,x2)
  real(8) :: func_y3, x1, x2
  func_y3 = x1 + x2
end function func_y3
!*
!*> =================================================================
!*>           SOLVED SVD PROBLEM WITH DGESDD LAPACK SUBROUTINE
!*> =================================================================
subroutine svd(m,n,a,snapshotMatrixU,eigenMatrix,vMatrix)
  implicit none
  external dgesdd
  !*------------------------------------------------
  integer ::  m, n
  real(8) :: a(m,n), aold(m,n)
  real(8) :: snapshotMatrixU(m,m),vMatrix(n,n),eigenMatrix(m)
  !*
  real(8), dimension(:,:), allocatable :: vt
  real(8), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: iwork
  integer :: ldu,lda,lwork
  integer :: ldvt,info,mx,mn,lwmax
  character(len=1) :: jobz
  !*------------------------------------------------
  !*
  mx=max(m,n)
  mn=min(m,n)
  !*
  aold = a
  !*
  jobz='a'
  lda=m
  ldu=m
  ldvt=max(1,n)
  !*
  lwork=-1
  lwmax=3*mn*mn+max(mx,4*mn*(mn+1))+500
  !*
  allocate(vt(n,n))
  allocate(work(lwmax))
  allocate(iwork(8*mn))
  !  print*, 'a', shape(a)
  !  print*, 'vt', shape(vt)
  !  print*, 'snapshotMatrixU', shape(snapshotMatrixU)
  !  print*, 'work', shape(work)
  !  print*, 'iwork', shape(iwork)

  call dgesdd( jobz, m, n, aold, lda, eigenMatrix, snapshotMatrixU, ldu, vt, ldvt, work, lwork, iwork, info )
  lwork = min( lwmax, int( work( 1 ) ) )
  !write(*,*) lwork,lwmax,work(1)
  call dgesdd(jobz, m, n, aold, lda, eigenMatrix, snapshotMatrixU, ldu, vt, ldvt, work, lwork, iwork, info )
  !*
  vMatrix=transpose(vt)

end subroutine svd
!*> =================================================================
!*>                           PSEUDO INVERSE
!*> =================================================================
subroutine pseudo_inverse(m,n,x,a_pseudo)
  implicit none
  integer :: i,j
  integer :: m,n
  integer :: Ener
  real(8) :: a_pseudo(n,m),x(m,n)
  real(8) :: eigenMatrix(m,n),snapshotMatrixU(m,m),vMatrix(n,n),sp(n,m)
  real(8) :: sing_val(m)
  real(8) :: tol=1.d-6
  
  call svd(m,n,x,snapshotMatrixU,sing_val,vMatrix)
  
  eigenMatrix=0.d0
  Ener=0
  do j=1,m
    if(sing_val(j).gt.tol) Ener=Ener+1
  enddo
  !*
  do j=1,Ener
    eigenMatrix(j,j)=sing_val(j)
  enddo
  
  sp(1:n,1:m)=0.d0
  do i=1,min(m,n)
    if (eigenMatrix(i,i).ne.0.d0)then
      sp(i,i)=1.d0/eigenMatrix(i,i)
    endif
  enddo
  
  a_pseudo(1:n,1:m) = matmul(vMatrix(1:n,1:n), &
    matmul(sp(1:n,1:m),transpose(snapshotMatrixU(1:m,1:m))))

end subroutine pseudo_inverse