program vulstraatgegevensaan

implicit none

integer     ::    i,j
integer, ALLOCATABLE :: aantalstrateninOSPMpunt(:)
real        ::   Xs,Ys
integer     :: aantalLijnbronnen,aantalOSPM,maxaantal
real, ALLOCATABLE  :: Xr(:),Yr(:),H(:),B(:)
real, ALLOCATABLE  :: L1(:),L2(:),theta(:),bijwelkelijn(:)
real, ALLOCATABLE  :: Xa(:),Ya(:),Xb(:),Yb(:),Bweg(:),Hweg(:),listexc(:,:)
integer, ALLOCATABLE :: kruisingstabel(:,:),kruisingstabelkort(:,:),temparray(:),nexc(:)
real :: dist,alpha,alpharad,m,noemer
integer :: aantalverwijderd
real :: tmpr
integer :: tmpi
character*6 :: tmpc

OPEN(UNIT=1,FILE='Lijnbronnen.txt')
READ(1,*) aantalLijnbronnen
write(*,*) aantalLijnbronnen
ALLOCATE(Xa(aantalLijnbronnen))
ALLOCATE(Xb(aantalLijnbronnen))
ALLOCATE(Ya(aantalLijnbronnen))
ALLOCATE(Yb(aantalLijnbronnen))
ALLOCATE(Bweg(aantalLijnbronnen))
ALLOCATE(Hweg(aantalLijnbronnen))

do i=1,aantalLijnbronnen
 READ(1,*) Xa(i),Ya(i),Xb(i),Yb(i),Bweg(i),Hweg(i)
enddo
write(*,*) Hweg(5896)
CLOSE(1)

OPEN(UNIT=1,FILE='straatgegevens.txt')
READ(1,*) aantalOSPM
ALLOCATE(aantalstrateninOSPMpunt(aantalOSPM))
ALLOCATE(Xr(aantalOSPM))
ALLOCATE(Yr(aantalOSPM))
ALLOCATE(H(aantalOSPM))
ALLOCATE(B(aantalOSPM))
ALLOCATE(L1(aantalOSPM))
ALLOCATE(L2(aantalOSPM))
ALLOCATE(theta(aantalOSPM))
ALLOCATE(bijwelkelijn(aantalOSPM))
ALLOCATE(kruisingstabel(aantalOSPM,100))
ALLOCATE(nexc(aantalOSPM))
ALLOCATE(listexc(aantalOSPM,36))
do i=1,aantalOSPM
 READ(1,*) Xr(i),Yr(i),H(i),B(i),L1(i),L2(i),theta(i),bijwelkelijn(i),nexc(i),listexc(i,:)
enddo
CLOSE(1)
i=50
kruisingstabel=0.

do i=1,aantalOSPM
if (mod(i,1000) .eq. 0) WRITE(*,*) i,aantalOSPM
if (theta(i) .ge. 180) then
 theta(i)=theta(i)-180
endif
if (theta(i) .ge. 180) then
 theta(i)=theta(i)-180
endif
if (theta(i) .lt. 0) then
 theta(i)=theta(i)+180
endif
if (theta(i) .lt. 0) then
 theta(i)=theta(i)+180
endif
if (theta(i) .lt. 0.1) then
 theta(i)=0.1
endif

if (L1(i) .lt. 30) then
L1(i)=30
endif
if (L2(i) .lt. 30) then
L2(i)=30
endif


! WRITE(*,*) i
 aantalstrateninOSPMpunt(i)=0
 do j=1,aantalLijnbronnen
  if ((Xa(j) .ne. Xb(j)) .and. (Ya(j) .ne. Yb(j)) .and. (theta(i) .ne. 0) .and. (theta(i) .ne. 90)) then
   alpha=90-theta(i)
   alpharad=alpha*3.141592654/180.
   m=tan(alpharad)
   noemer=(Yb(j)-Ya(j))/(Xb(j)-Xa(j))+1/m
   if (noemer .eq. 0) then
    goto 100
   endif
   Xs=(Yr(i)+Xr(i)/m-Yb(j)+(Yb(j)-Ya(j))/(Xb(j)-Xa(j))*Xb(j))/noemer
   Ys=-1/m*(Xs-Xr(i))+Yr(i)
 elseif ((Xa(j) .ne. Xb(j)) .and. (Ya(j) .ne. Yb(j)) .and. (theta(i) .eq. 0)) then
   Xs=Xr(i)
   Ys=(Yb(j)-Ya(j))/(Xb(j)-Xa(j))*(Xr(i)-Xb(j))+Xb(j)
  elseif ((Xa(j) .ne. Xb(j)) .and. (Ya(j) .ne. Yb(j)) .and. (theta(i) .eq. 0)) then
   Ys=Yr(i)
   Xs=(Xb(j)-Xa(j))/(Yb(j)-Ya(j))*(Yr(i)-Yb(j))+Xb(j)
  elseif ((Xa(j) .eq. Xb(j)) .and. (theta(i) .ne. 0) .and. (theta(i) .ne. 90)) then 
   alpha=90-theta(i)
   alpharad=alpha*3.141592654/180.
   m=tan(alpharad)
   Xs=Xa(j)
   Ys=-1/m*(Xa(j)-Xr(i))+Yr(i)
   !if (Xr(i) .eq. 153.5864407) then
   !WRITE(*,*) j,i,theta(i),alpha,alpharad,m,Xs,Ys
   !endif
  elseif ((Ya(j) .eq. Yb(j)) .and. (theta(i) .ne. 0) .and. (theta(i) .ne. 90)) then 
   alpha=90-theta(i)
   alpharad=alpha*3.141592654/180.
   m=tan(alpharad)
   Ys=Ya(j)
   Xs=-m*(Ya(j)-Yr(i))+Xr(i)
  elseif ((Xa(j) .eq. Xb(j)) .and. (theta(i) .eq. 90)) then
   Xs=Xa(j)
   Ys=Yr(i)
  elseif ((Ya(j) .eq. Yb(j)) .and. (theta(i) .eq. 0)) then
   Xs=Xr(i)
   Ys=Ya(j)
  else
   goto 100
  endif
  dist=sqrt((Xs-Xr(i))**2+(Ys-Yr(i))**2)
 
 if (Xa(j) .ne. Xb(j)) then
  if (dist .le. B(i)/1000. .and. (Xs-0.0001 .le. MAX(Xa(j),Xb(j))) .and. (Xs+0.0001 .gt. MIN(Xa(j),Xb(j))) .and. (Hweg(j) .le. H(i))) then
   aantalstrateninOSPMpunt(i)=aantalstrateninOSPMpunt(i)+1
!   WRITE(*,*) i,j,aantalstrateninOSPMpunt(i)
   kruisingstabel(i,aantalstrateninOSPMpunt(i))=j
  endif
 else
  if (dist .le. B(i)/1000. .and. (Ys-0.0001 .le. MAX(Ya(j),Yb(j))) .and. (Ys+0.0001 .gt. MIN(Ya(j),Yb(j))) .and. (Hweg(j) .le. H(i))) then
   aantalstrateninOSPMpunt(i)=aantalstrateninOSPMpunt(i)+1
!   WRITE(*,*) i,j,aantalstrateninOSPMpunt(i)
   kruisingstabel(i,aantalstrateninOSPMpunt(i))=j
  endif
 endif
  100 continue
 enddo
enddo

maxaantal=0
do i=1,aantalOSPM
if (aantalstrateninOSPMpunt(i) .gt. maxaantal) then
 maxaantal=aantalstrateninOSPMpunt(i)
endif
enddo

ALLOCATE(kruisingstabelkort(aantalOSPM,maxaantal))
ALLOCATE(temparray(maxaantal))

kruisingstabelkort=kruisingstabel(:,1:maxaantal)
WRITE(*,*) maxaantal
OPEN(UNIT=1,FILE='straatgegevens_config_exc.txt')
WRITE(1,*) aantalOSPM,maxaantal
do i=1,aantalOSPM
 temparray=kruisingstabel(i,1:maxaantal)
 WRITE(1,*) Xr(i),Yr(i),H(i),B(i),L1(i),L2(i),theta(i),bijwelkelijn(i),nexc(i),listexc(i,:),aantalstrateninOSPMpunt(i),temparray,'end'
enddo
CLOSE(1)

end program
