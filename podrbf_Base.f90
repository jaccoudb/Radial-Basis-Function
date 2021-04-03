program pod_rbf
implicit none
  integer :: i,j,cond
  integer :: info
  integer :: nsamples,nsnaps,nparam,Ener
  real(8) :: tol,c
  real(8),allocatable,dimension(:,:) :: xfield,Pfield
  real(8),allocatable,dimension(:,:) :: xfield_res_t,xfield_res_t_2
  real(8),allocatable,dimension(:,:) :: u,v,s
  real(8),allocatable,dimension(:,:) :: Phi,G,A,B
  real(8),allocatable,dimension(:,:) :: gP,xfield_res_t_3
  real(8),allocatable,dimension(:) :: e,sing_val,pX,minP,maxP
  !*===============================================
  parameter(tol=1.d-4)
  !*===============================================
  !*  READ AND ALLOCATE VARIABLES FOR SNAPSHOTS
  !*===============================================
  open(999,file='snapshots.dat',status='old')
  read(999,*)  nsamples,nsnaps
  !*
  allocate(xfield(nsamples,nsnaps),e(nsamples))
  allocate(xfield_res_t(nsamples,nsnaps),xfield_res_t_2(nsamples,nsnaps))
  allocate(u(nsamples,nsamples),sing_val(nsamples),v(nsnaps,nsnaps),s(nsamples,nsnaps))
  allocate(Phi(nsnaps,nsnaps))
  allocate(G(nsnaps,nsnaps),A(nsamples,nsnaps),B(nsamples,nsnaps))
  allocate(gP(nsnaps,1),xfield_res_t_3(nsamples,1))
  !*
  do i=1,nsamples
    read(999,*) (xfield(i,j),j=1,nsnaps)
  enddo
  !*
  !*===============================================
  !* DETERMINATE SVD RESULTS
  !*===============================================
  !*
  call dsvdc(xfield,nsamples,nsnaps,sing_val,e,u,v,11,info)
  !*
write(997,*) (sing_val(j),j=1,nsamples)
  !*===============================================
  !*  READ AND ALLOCATE VARIABLES FOR PARAMETERS
  !*===============================================
  open(998,file='parameters.dat',status='old')
  read(998,*) nparam
  !*
  allocate(Pfield(nparam,nsnaps),pX(nparam))
  allocate(minP(nparam),maxP(nparam))
  do i=1,nparam
    read(998,*) (Pfield(i,j),j=1,nsnaps)
  enddo
  c=0.d2/dble(nsnaps)
  !*====================================================
  !* CALCULATE ENERGY
  !*====================================================
  Ener=0
  do j=1,nsamples
    if(sing_val(j).gt.tol) Ener=Ener+1
  enddo
  !*
  s=0.d0
  do i=1,Ener
    s(i,i)=sing_val(i)
  enddo
  !*
  !*===============================================
  !* FIELD TEST SVD
  !*===============================================
!  xfield_res_t(1:nsamples,1:nsnaps)=matmul(u(1:nsamples,1:nsamples),&
!    matmul(s(1:nsamples,1:nsnaps),transpose(v(1:nsnaps,1:nsnaps))))
  xfield_res_t(1:nsamples,1:nsnaps)=matmul(u(1:nsamples,1:Ener),&
    matmul(s(1:Ener,1:Ener),transpose(v(1:nsnaps,1:Ener))))
  !*===============================================
  !* PRINT SVD RESULTS
  !*===============================================
  call output(nsamples,nsamples,1,'matrix_U',u)
  call output(nsnaps,nsnaps,1,'matrix_V',v)
  call output(nsamples,nsnaps,1,'matrix_S',s)
  call output(nsamples,1,1,'vetor_S',sing_val)
  call output(nsamples,nsnaps,1,'xfield_res_t',xfield_res_t) !ok!
  !*===============================================
  !* DETERMINATE RBF FOR PARAMETERS
  !*===============================================
  call construct_rbf(nparam,nsnaps,c,Pfield,minP,maxP,Phi)
  call output(nsnaps,nsnaps,1,'Phi',Phi)
  !*===============================================
  !* DETERMINATE RBF FOR TRAINING PARAMETERS
  !*===============================================
  A=matmul(s,transpose(v))
  call output(nsamples,nsnaps,1,'A',A)
  !*
  G=0.d0  
  call pseudo_inverse(nsnaps,nsnaps,Phi,G)   !ok!
  call output(nsnaps,nsnaps,1,'G',G)
  !*
  B=matmul(A,G)
  call output(nsamples,nsnaps,1,'B',B)
  !*
  xfield_res_t_2(1:nsamples,1:nsnaps)=matmul(u(1:nsamples,1:nsamples),&
    matmul(B(1:nsamples,1:nsnaps),G(1:nsnaps,1:nsnaps)))
  !*
  call output(nsamples,nsnaps,1,'xfield_res_t_2',xfield_res_t_2)
  !*
  !*===============================================
  !* RUN RESPONSE SURFACE
  !*===============================================
  pX(1)=(135.d3-minP(1))/(maxP(1)-minP(1))
  pX(2)=(11.d0-minP(2))/(maxP(2)-minP(2))
  !*
  call forward_problem(nparam,nsnaps,c,pX,Pfield,gP)
  call output(nsnaps,1,1,'gP',gP)
  !*
  xfield_res_t_3(1:nsamples,1)=matmul(u(1:nsamples,1:nsamples),&
    matmul(B(1:nsamples,1:nsnaps),gP(1:nsnaps,1)))

  call output(nsamples,1,1,'xfield_res_t_3',xfield_res_t_3)
!*
end program pod_rbf
!*
!*
!####################################################################
!####################################################################
!#  SUBROTINAS SVD
!####################################################################
!*===================================================================
!*
!*===================================================================
subroutine drotg(da,db,dc,ds)
implicit none
!*===================================================================
!* DESIGNED BY C.L.LAWSON, JPL, 1977 SEPT 08
!* CONSTRUCT THE GIVENS TRANSFORMATION
!*
!*     ( DC  DS )
!* G = (        ) ,    DC**2 + DS**2 = 1 ,
!*     (-DS  DC )
!*
!* WHICH ZEROS THE SECOND ENTRY OF THE 2-VECTOR  (DA,DB)**T .
!*
!* THE QUANTITY R = (+/-)SQRT(DA**2 + DB**2) OVERWRITES DA IN
!* STORAGE.  THE VALUE OF DB IS OVERWRITTEN BY A VALUE Z WHICH
!* ALLOWS DC AND DS TO BE RECOVERED BY THE FOLLOWING ALGORITHM:
!*
!*   IF Z=1  SET  DC=0.D0  AND  DS=1.D0
!*   IF DABS(Z) < 1  SET  DC=SQRT(1-Z**2)  AND  DS=Z
!*   IF DABS(Z) > 1  SET  DC=1/Z  AND  DS=SQRT(1-DC**2)
!*
!* NORMALLY, THE SUBPROGRAM DROT(N,DX,INCX,DY,INCY,DC,DS) WILL
!* NEXT BE CALLED TO APPLY THE TRANSFORMATION TO A 2 BY N MATRIX.
!*===================================================================
real(8), intent(in out) :: da
real(8), intent(in out) :: db
real(8), intent(out) :: dc
real(8), intent(out) :: ds
real(8)  :: u, v, r

if (dabs(da).le.dabs(db)) go to 10

!* HERE ABS(DA) > ABS(DB)
u = da + da
v = db / u

!* NOTE THAT U AND R HAVE THE SIGN OF DA
r = dsqrt(.25d0 + v**2) * u

!* NOTE THAT DC IS POSITIVE
dc = da / r
ds = v * (dc + dc)
db = ds
da = r
return

!* HERE ABS(DA) <= ABS(DB)
10 if (db.eq.0.d0) go to 20
u = db + db
v = da / u

!* NOTE THAT U AND R HAVE THE SIGN OF DB
!* (R IS IMMEDIATELY STORED IN DA)

da = dsqrt(0.25D0 + v**2) * u

!* NOTE THAT DS IS POSITIVE

ds = db / da
dc = v * (ds + ds)
if (dc.eq.0.d0) go to 15
db = 1.d0 / dc
return
15 db = 1.d0
return

!* HERE DA = DB = 0.D0

20 dc = 1.d0
ds = 0.d0

end subroutine drotg
!*
!*===================================================================
!*
!*===================================================================
subroutine dswap1(n, dx, dy)
implicit none
!* INTERCHANGES TWO VECTORS.
!* USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
!* JACK DONGARRA, LINPACK, 3/11/78.
!* This version is for increments = 1.
integer, intent(in) :: n
real(8), intent(in out) :: dx(*)
real(8), intent(in out) :: dy(*)
real(8) :: dtemp
integer :: i, m, mp1

if(n.le.0) return

!* CODE FOR BOTH INCREMENTS EQUAL TO 1
!* CLEAN-UP LOOP

m = mod(n,3)
if(m.eq.0) go to 40
do  i = 1,m
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
end do

if( n.lt.3 ) return
40 mp1 = m + 1
do  i = mp1,n,3
  dtemp = dx(i)
  dx(i) = dy(i)
  dy(i) = dtemp
  dtemp = dx(i + 1)
  dx(i + 1) = dy(i + 1)
  dy(i + 1) = dtemp
  dtemp = dx(i + 2)
  dx(i + 2) = dy(i + 2)
  dy(i + 2) = dtemp
end do

end subroutine dswap1
!*
!*===================================================================
!*
!*===================================================================
subroutine drot1(n, dx, dy, c, s)
implicit none
!* APPLIES A PLANE ROTATION.
!* JACK DONGARRA, LINPACK, 3/11/78.
!* This version is for increments = 1.
integer, intent(in) :: n
real(8), intent(in out) :: dx(*)
real(8), intent(in out) :: dy(*)
real(8), intent(in) :: c
real(8), intent(in) :: s
real(8) :: dtemp
integer :: i

if(n.le.0) return

!* CODE FOR BOTH INCREMENTS EQUAL TO 1
do  i = 1,n
  dtemp = c*dx(i) + s*dy(i)
  dy(i) = c*dy(i) - s*dx(i)
  dx(i) = dtemp
end do

end subroutine drot1
!*
!*===================================================================
!*
!*===================================================================
!module sub_globais
!implicit none
!contains
!* ...


!end module sub_globais

subroutine dsvdc(x,n,p,s,e,u,v,job,info)
implicit none
!integer, intent(in) :: n,p
!real(8), intent(in out) :: x(:,:)
!real(8), intent(out) :: s(:),e(:),u(:,:),v(:,:)
!integer, intent(in) :: job
!integer, intent(out) :: info

integer :: n,p
real(8) :: x(n,p)
real(8) :: s(n),e(n),u(n,n),v(p,p)
integer :: job
integer :: info

!* DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X
!* BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
!* DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
!* COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
!* AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
!*
!* ON ENTRY
!*
!* X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N.
!*           X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
!*           DECOMPOSITION IS TO BE COMPUTED.
!*           X IS DESTROYED BY DSVDC.
!*
!* LDX       INTEGER.
!*           LDX IS THE LEADING DIMENSION OF THE ARRAY X.
!*
!* N         INTEGER.
!*           N IS THE NUMBER OF ROWS OF THE MATRIX X.
!*
!* P         INTEGER.
!*           P IS THE NUMBER OF COLUMNS OF THE MATRIX X.
!*
!* LDU       INTEGER.
!*           LDU IS THE LEADING DIMENSION OF THE ARRAY U.
!*          (SEE BELOW).
!*
!* LDV       INTEGER.
!*           LDV IS THE LEADING DIMENSION OF THE ARRAY V.
!*           (SEE BELOW).
!*
!* JOB       INTEGER.
!*           JOB CONTROLS THE COMPUTATION OF THE SINGULAR
!*           VECTORS.  IT HAS THE DECIMAL EXPANSION AB
!*           WITH THE FOLLOWING MEANING
!*
!*           A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR VECTORS.
!*           A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS IN U.
!*           A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR VECTORS IN U.
!*           B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR VECTORS.
!*           B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS IN V.
!*
!* ON RETURN
!*
!* S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P).
!*           THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE SINGULAR
!*           VALUES OF X ARRANGED IN DESCENDING ORDER OF MAGNITUDE.
!*
!* E         DOUBLE PRECISION(P).
!*           E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
!*           DISCUSSION OF INFO FOR EXCEPTIONS.
!*
!* U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF
!*           JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2 THEN K.EQ.MIN(N,P).
!*           U CONTAINS THE MATRIX OF LEFT SINGULAR VECTORS.
!*           U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
!*           OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
!*           IN THE SUBROUTINE CALL.
!*
!* V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P.
!*           V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
!*           V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
!*           THEN V MAY BE IDENTIFIED WITH X IN THE
!*           SUBROUTINE CALL.
!*
!* INFO      INTEGER.
!*           THE SINGULAR VALUES (AND THEIR CORRESPONDING SINGULAR
!*           VECTORS) S(INFO+1),S(INFO+2),...,S(M) ARE CORRECT
!*           (HERE M=MIN(N,P)).  THUS IF INFO.EQ.0, ALL THE
!*           SINGULAR VALUES AND THEIR VECTORS ARE CORRECT.
!*           IN ANY EVENT, THE MATRIX B = TRANS(U)*X*V IS THE
!*           BIDIAGONAL MATRIX WITH THE ELEMENTS OF S ON ITS DIAGONAL
!*           AND THE ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
!*           IS THE TRANSPOSE OF U).  THUS THE SINGULAR VALUES
!*           OF X AND B ARE THE SAME.
!*
!* LINPACK. THIS VERSION DATED 03/19/79 .
!* G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
!* DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS.
!* EXTERNAL DROT
!* BLAS DAXPY,DDOT,DSCAL,DSWAP,DNRM2,DROTG
!* FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT

!* INTERNAL VARIABLES
integer :: iter, j, jobu, k, kase, kk, l, ll, lls,&
           lm1, lp1, ls, lu, m, maxit, mm, mm1, mp1,&
		   nct, nctp1, ncu, nrt, nrtp1
real(8) :: t, work(n)
real(8) :: b, c, cs, el, emm1, f, g, scale, shift, sl,&
           sm, sn, smm1, t1, test, ztest
logical :: wantu, wantv

!* SET THE MAXIMUM NUMBER OF ITERATIONS.
maxit = 30

!* DETERMINE WHAT IS TO BE COMPUTED.
wantu = .false.
wantv = .false.
jobu = mod(job,100)/10
ncu = n

if (jobu.gt.1) ncu = min(n,p)
if (jobu.ne.0) wantu = .true.
if (mod(job,10).ne.0) wantv = .true.

!* REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
!* IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.

info = 0
nct = min(n-1, p)
s(1:nct+1) = 0.0d0
nrt = max(0, min(p-2,n))
lu = max(nct,nrt)
if (lu.lt.1) go to 170
do  l = 1, lu
  lp1 = l + 1
  if (l.gt.nct) go to 20
  
  !* COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
  !* PLACE THE L-TH DIAGONAL IN S(L).

  s(l) = dsqrt(sum(x(l:n,l)**2.d0))
  if (s(l).eq.0.0d0) go to 10
  if (x(l,l).ne.0.0d0) s(l) = sign(s(l), x(l,l))
  x(l:n,l) = x(l:n,l) / s(l)
  x(l,l) = 1.0d0 + x(l,l)

  10 s(l) = -s(l)

  20 if (p.lt.lp1) go to 50
  do  j = lp1, p
    if (l.gt.nct) go to 30
    if (s(l).eq.0.0d0) go to 30
    
    !* APPLY THE TRANSFORMATION.
    t = -dot_product(x(l:n,l), x(l:n,j)) / x(l,l)
    x(l:n,j) = x(l:n,j) + t * x(l:n,l)
    
    !* PLACE THE L-TH ROW OF X INTO  E FOR THE
    !* SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
    30 e(j) = x(l,j)
  end do

  50 if (.not.wantu .or. l.gt.nct) go to 70
  
  !* PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK MULTIPLICATION.
  u(l:n,l) = x(l:n,l)
  70 if (l.gt.nrt) cycle
  
  !* COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
  !* L-TH SUPER-DIAGONAL IN E(L).
  e(l) = dsqrt(sum(e(lp1:p)**2.d0))
  if (e(l).eq.0.0d0) go to 80
  if (e(lp1).ne.0.0d0) e(l) = sign(e(l), e(lp1))
  e(lp1:lp1+p-l-1) = e(lp1:p) / e(l)
  e(lp1) = 1.0d0 + e(lp1)

  80 e(l) = -e(l)
  if (lp1.gt.n .or. e(l).eq.0.0d0) go to 120
  
  !* APPLY THE TRANSFORMATION.
  work(lp1:n) = 0.0D0
  do  j = lp1, p
    work(lp1:lp1+n-l-1) = work(lp1:lp1+n-l-1) + e(j) * x(lp1:lp1+n-l-1,j)
  end do
  do  j = lp1, p
    x(lp1:lp1+n-l-1,j) = x(lp1:lp1+n-l-1,j) - (e(j)/e(lp1)) * work(lp1:lp1+n-l-1)
  end do

  120 if (.not.wantv) cycle
  
  !* PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
  !* BACK MULTIPLICATION.
  v(lp1:p,l) = e(lp1:p)
end do

!* SET UP THE FINAL BIDIAGONAL MATRIX OF ORDER M.
170 m = min(p,n+1)
nctp1 = nct + 1
nrtp1 = nrt + 1
if (nct.lt.p) s(nctp1) = x(nctp1,nctp1)
if (n.lt.m) s(m) = 0.0d0
if (nrtp1.lt.m) e(nrtp1) = x(nrtp1,m)
e(m) = 0.0d0

!* IF REQUIRED, GENERATE U.
if (.not.wantu) go to 300
if (ncu.lt.nctp1) go to 200
do  j = nctp1, ncu
  u(1:n,j) = 0.0d0
  u(j,j) = 1.0d0
end do

200 do  ll = 1, nct
  l = nct - ll + 1
  if (s(l).eq.0.0d0) go to 250
  lp1 = l + 1
  if (ncu.lt.lp1) go to 220
  do  j = lp1, ncu
    t = -dot_product(u(l:n,l), u(l:n,j)) / u(l,l)
    u(l:n,j) = u(l:n,j) + t * u(l:n,l)
  end do

  220 u(l:n,l) = -u(l:n,l)
  u(l,l) = 1.0d0 + u(l,l)
  lm1 = l - 1
  if (lm1.lt.1) cycle
  u(1:lm1,l) = 0.0d0
  cycle

  250 u(1:n,l) = 0.0d0
  u(l,l) = 1.0d0
end do

!* IF IT IS REQUIRED, GENERATE V.
300 if (.not.wantv) go to 350
do  ll = 1, p
  l = p - ll + 1
  lp1 = l + 1
  if (l.gt.nrt) go to 320
  if (e(l).eq.0.0d0) go to 320
  do  j = lp1, p
    t = -dot_product(v(lp1:lp1+p-l-1,l), v(lp1:lp1+p-l-1,j)) / v(lp1,l)
    v(lp1:lp1+p-l-1,j) = v(lp1:lp1+p-l-1,j) + t * v(lp1:lp1+p-l-1,l)
  end do

  320 v(1:p,l) = 0.0d0
  v(l,l) = 1.0d0
end do

!* MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
350 mm = m
iter = 0

!* QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
!* ...EXIT
360 if (m.eq.0) go to 620

!* IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET FLAG AND RETURN.
if (iter.lt.maxit) go to 370
info = m
!* ...EXIT
go to 620
!* THIS SECTION OF THE PROGRAM INSPECTS FOR NEGLIGIBLE ELEMENTS
!* IN THE S AND E ARRAYS.  ON COMPLETION
!* THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
!* 
!* KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L < M
!* KASE = 2     IF S(L) IS NEGLIGIBLE AND L < M
!* KASE = 3     IF E(L-1) IS NEGLIGIBLE, L < M, AND
!*              S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
!* KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
370 do  ll = 1, m
  l = m - ll
!*...EXIT
  if (l.eq.0) exit
  test = dabs(s(l)) + dabs(s(l+1))
  ztest = test + dabs(e(l))
  if (ztest.ne.test) cycle
  e(l) = 0.0d0
!*...EXIT
  exit
end do

if (l.ne.m - 1) go to 410
kase = 4
go to 480

410 lp1 = l + 1
mp1 = m + 1
do  lls = lp1, mp1
  ls = m - lls + lp1
!*...EXIT
  if (ls.eq.l) exit
  test = 0.0d0
  if (ls.ne.m) test = test + dabs(e(ls))
  if (ls.ne.l + 1) test = test + dabs(e(ls-1))
  ztest = test + dabs(s(ls))
  if (ztest.ne.test) cycle
  s(ls) = 0.0d0
!*...EXIT
  exit
end do

if (ls.ne.l) go to 450
kase = 3
go to 480

450 if (ls.ne.m) go to 460
kase = 1
go to 480

460 kase = 2
l = ls
480 l = l + 1

!* PERFORM THE TASK INDICATED BY KASE.

select case (kase)
  case(1)
    go to 490
  case(2)
    go to 520
  case(3)
    go to 540
  case(4)
    go to 570
end select

!* DEFLATE NEGLIGIBLE S(M).

490 mm1 = m - 1
f = e(m-1)
e(m-1) = 0.0d0
do  kk = l, mm1
  k = mm1 - kk + l
  t1 = s(k)
  call drotg(t1,f,cs,sn)
  s(k) = t1
  if (k.eq.l) go to 500
  f = -sn*e(k-1)
  e(k-1) = cs*e(k-1)

  500 if (wantv) call drot1(p,v(1:,k),v(1:,m),cs,sn)
end do
go to 610

!* SPLIT AT NEGLIGIBLE S(L).

520 f = e(l-1)
e(l-1) = 0.0d0
do  k = l, m
  t1 = s(k)
  call drotg(t1, f, cs, sn)
  s(k) = t1
  f = -sn*e(k)
  e(k) = cs*e(k)
  if (wantu) call drot1(n, u(1:,k), u(1:,l-1), cs, sn)
end do
go to 610

!* PERFORM ONE QR STEP.
!* CALCULATE THE SHIFT.

540 scale = max(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)),dabs(s(l)),dabs(e(l)))
sm = s(m)/scale
smm1 = s(m-1)/scale
emm1 = e(m-1)/scale
sl = s(l)/scale
el = e(l)/scale
b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0D0
c = (sm*emm1)**2
shift = 0.0D0
if (b.eq.0.0d0 .and. c.eq.0.0d0) go to 550
shift = dsqrt(b**2+c)
if (b.lt.0.0d0) shift = -shift
shift = c/(b + shift)

550 f = (sl + sm)*(sl - sm) - shift
g = sl*el

!* CHASE ZEROS.
mm1 = m - 1
do  k = l, mm1
  call drotg(f,g,cs,sn)
  if (k.ne.l) e(k-1) = f
  f = cs*s(k) + sn*e(k)
  e(k) = cs*e(k) - sn*s(k)
  g = sn*s(k+1)
  s(k+1) = cs*s(k+1)
  if (wantv) call drot1(p,v(1:,k),v(1:,k+1),cs,sn)
  call drotg(f,g,cs,sn)
  s(k) = f
  f = cs*e(k) + sn*s(k+1)
  s(k+1) = -sn*e(k) + cs*s(k+1)
  g = sn*e(k+1)
  e(k+1) = cs*e(k+1)
  if (wantu.and.k.lt.n) call drot1(n,u(1:,k),u(1:,k+1),cs,sn)
end do
e(m-1) = f
iter = iter + 1
go to 610

!* CONVERGENCE.
!* MAKE THE SINGULAR VALUE  POSITIVE.
570 if (s(l).ge.0.0d0) go to 590
s(l) = -s(l)
if (wantv) v(1:p,l) = -v(1:p,l)
!* ORDER THE SINGULAR VALUE.
590 if (l.eq.mm) go to 600
!...exit
if (s(l).ge.s(l+1)) go to 600
t = s(l)
s(l) = s(l+1)
s(l+1) = t
if (wantv.and.l.lt.p) call dswap1(p, v(1:,l), v(1:,l+1))
if (wantu.and.l.lt.n) call dswap1(n, u(1:,l), u(1:,l+1))
l = l + 1
go to 590

600 iter = 0
m = m - 1

610 go to 360

620 return
end subroutine dsvdc
!*
!*===============================================
!* SVD_PRODUCT_TEST tests that A = U * S * V'.
!*===============================================
!*  M, N, the number of rows and columns.
!*
!*  A(M,N), the matrix whose singular value
!*  decomposition we are investigating.
!*
!*  U(M,M), S(M,N), V(N,N), the factors
!*  that form the singular value decomposition of A.
!*
subroutine svd_product_test(m,n,a,u,s,v)
  implicit none
  integer(4) :: i,j
  integer(4) :: m,n
  real(8) :: a(m,n)
  real(8) :: u(m,m),s(m,n),v(n,n),usv(m,n)
  real(8) :: a_norm,dif_norm
  real(8) :: mat_dif_fro,mat_norm_fro
  !*--------------------
  a_norm=mat_norm_fro(m,n,a)

  usv(1:m,1:n) = matmul ( u(1:m,1:m), &
                 matmul ( s(1:m,1:n), transpose ( v(1:n,1:n) ) ) )

!  call r8mat_print ( m, n, usv, '  The product U * S * V'':' )
  write(997,*) "=========================================================================="
  write(997,*) ' The SVD test'
  write(997,*)
  do i=1,m
    write(997,*) (usv(i,j),j=1,n)
  enddo

  dif_norm = mat_dif_fro( m, n, a, usv )

  write (997, '(a)' ) ' '
  write (997, '(a,g14.6)' ) '  Frobenius Norm of A, A_NORM = ', a_norm
  write (997, '(a)' ) ' '
  write (997, '(a)' ) '  ABSOLUTE ERROR for A = U*S*V'':'
  write (997, '(a,g14.6)' ) '  Frobenius norm of difference A-U*S*V'' = ', &
    dif_norm
  write (997, '(a)' ) ' '
  write (997, '(a)' ) '  RELATIVE ERROR for A = U*S*V'':'
  write (997, '(a,g14.6)' ) '  Ratio of DIF_NORM / A_NORM = ', &
    dif_norm / a_norm
  !*--------------------
end subroutine svd_product_test
!*
!*===============================================
!* MAT_NORM_FRO returns the Frobenius norm of an MAT.
!*===============================================
!*  The Frobenius norm is defined as
!*
!*  MAT_NORM_FRO = sqrt (sum(1<=I<=M)
!*                                     sum ( 1 <= j <= N ) A(I,J)^2 )
!*
!*  The matrix Frobenius norm is not derived from a vector norm, but
!*  is compatible with the vector L2 norm, so that:
!*
!*  vec_norm_l2 ( A * x ) <= mat_norm_fro ( A ) * vec_norm_l2 ( x ).
!*
!*  Parameters:
!*
!*  M, the number of rows in A; N, the number of columns in A.
!*
!*  A(M,N), the matrix whose Frobenius norm is desired.
!*  Output, MAT_NORM_FRO, the Frobenius norm of A.
!*
real(8) function mat_norm_fro(m,n,a)
  implicit none
  integer(4) :: m,n
  real(8) a(m,n)
  real ( kind = 8 ) r8mat_norm_fro
  !*--------------------
  mat_norm_fro=dsqrt(sum(a(1:m,1:n)**2.d0))
  return
  !*--------------------
end function mat_norm_fro
!*
!*===============================================
!* MAT_DIF_FRO, the Frobenius norm of the difference of two R8MAT's.
!*===============================================
!*  The Frobenius norm is defined as
!*
!*  MAT_DIF_FRO = sqrt (sum(1<=I<=M) 
!*                                  sum ( 1 <= j <= N ) ( A(I,J) - B(I,J) )^2 )
!*
!*   M, the number of rows; N, the number of columns.
!*
!*  Output, MAT_DIF_FRO, the Frobenius norm of
!*  the difference of A and B.
!*
real(8) function mat_dif_fro(m,n,a,b)
  implicit none
  integer(4) :: m,n
  real(8) :: a(m,n),b(m,n)
  !*--------------------
  mat_dif_fro=dsqrt(sum((a(1:m,1:n)-b(1:m,1:n))**2.d0))
  return
  !*--------------------
end function mat_dif_fro
!*
!*===================================================================
!*  PSEUDO INVERSE
!*===================================================================
subroutine pseudo_inverse(m,n,x,a_pseudo)
  implicit none
  integer :: i,j
  integer :: m,n
  integer :: info,Ener
  real(8) :: a_pseudo(n,m),x(m,n)
  real(8) :: s(m,n),u(m,m),v(n,n),sp(n,m)
  real(8) :: sing_val(m),e(m)
  real(8) :: at_pseudo(n,m)
  real(8) :: tol=1.d-6

  call dsvdc(x,m,n,sing_val,e,u,v,11,info)

  s=0.d0
  Ener=0
  do j=1,m
    if(sing_val(j).gt.tol) Ener=Ener+1
  enddo

  do j=1,Ener
    s(j,j)=sing_val(j)
  enddo

  sp(1:n,1:m)=0.d0
  do i=1,min(m,n)
    if (s(i,i).ne.0.d0)then
      sp(i,i)=1.d0/s(i,i)
    endif
  enddo

  a_pseudo(1:n,1:m) = matmul(v(1:n,1:n), &
    matmul(sp(1:n,1:m),transpose(u(1:m,1:m))))

  !at_pseudo(1:m,1:n) = matmul(u(1:m,1:m), &
  !  matmul(s(1:m,1:n),transpose(v(1:n,1:n))))
  !err=matmul(x(1:m,1:n)-at_pseudo(1:m,1:n))
  !E=matmul(x(1:m,1:n)-a_pseudo(1:n,1:m))


  !*
  return
end subroutine pseudo_inverse
!*
!*===================================================================
!*
!*===================================================================
subroutine podBmtx(m,n,P,minP,maxP,x)
  implicit none
  integer :: i,j
  integer :: n,m
  real(8) :: minP(m),maxP(m),P(m,n)
  real(8) :: x(m,n)
  !*
  do i=1,m
    minP(i)=minval(P(i,:))
    maxP(i)=maxval(P(i,:))
    do j=1,n
      x(i,j)=(P(i,j)-minP(i))/(maxP(i)-minP(i))
!      x(i,j)=(P(i,j))/(maxP(i)-minP(i))
    enddo
  enddo
  !*
end subroutine podBmtx
!*
!*===================================================================
!*  CONSTRUCT A SET OF RADIAL BASIS FUNCTION (SPLINE LINEAR)
!*===================================================================
subroutine construct_rbf(m,n,c,P,minP,maxP,phi)
  implicit none
  integer :: i,j
  integer :: m,n
  real(8) :: c
  real(8) :: x(m,n),phi(n,n),xi(m),xj(m)
  real(8) :: P(m,n),minP(m),maxP(m)
  !*
  call podBmtx(m,n,P,minP,maxP,x)
  !*
  do i=1,n
    xi(:)=x(:,i)
    do j=1,n
      xj(:)=x(:,j)-xi(:)
      phi(i,j)=dsqrt(dot_product(xj,xj)+c*c)
    enddo
  enddo
  !*
end subroutine construct_rbf
!*
!*===================================================================
!*  FORWARD RBF PROBLEM - COMPUTE VALUES
!*===================================================================
subroutine forward_problem(m,n,c,pX,P,vy)
  implicit none
  integer :: i,j
  integer :: n,m
  real(8) :: c
  real(8) :: vy(n)
  real(8) :: x(m,n),pX(m)
  real(8) :: P(m,n),minP(m),maxP(m)
  !*
  call podBmtx(m,n,P,minP,maxP,x)
  !*
  do i=1,n
    call solve_rbf(m,n,x,pX,vy,c)
  enddo
  !*
end subroutine forward_problem
!*
!*===================================================================
!*  INTERPOLATE A FUNCTION USING A PRE-CONSTRUCTED RADIAL BASIS
!*  FUNCTION COEFFICIENTS
!*===================================================================
subroutine solve_rbf(m,n,x,x1,xi,c)
  implicit none
  integer :: i,j
  integer :: n,m
  real(8) :: f1,c
  real(8) :: x(m,n),x1(m),xj(m),xi(n)
  !*  
  f1=0.d0
  !*
  do j=1,n
    xj(:)=x1(:)-x(:,j)
    f1=DSQRT(DOT_PRODUCT(xj,xj)+c*c)
    xi(j)=f1
  enddo
  !*
end subroutine solve_rbf
!*
!*===================================================================
!*  PRINT RESULTS
!*===================================================================
subroutine output(m,n,k,name_variable,x)
implicit none
  integer :: i,j,k
  character(len=*):: name_variable 
  integer :: m,n
  real(8) :: x(m,n)
  !*
  open(unit=k,file=trim(name_variable)//".dat")
  do i=1,m
    write(k,*) (x(i,j),j=1,n)
  end do
  close(unit=k)
end subroutine output
