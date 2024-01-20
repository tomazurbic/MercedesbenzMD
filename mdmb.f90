module particles
   real*8, dimension(:,:), allocatable :: delci     ! positions
   real*8, dimension(:,:), allocatable :: vel     ! velocities
   real*8, dimension(:,:), allocatable :: acc     ! accelerations
   real*8 :: ekin,tc
end module particles

module parametri
   real*8 :: a1,a2  !MB LJ parameters
   real*8 :: s1,s2,rh,eh !MB parameters
   real*8 :: temp,g !temperature and density
   real*8 :: tstep,dt2,d2t  !time step, dt2=tstep/2, d2t=tstep**2/2
   integer :: tsteps,esteps !number of time steps for sampling and equalibration
   integer :: nseries !number of series
   integer :: N  !number of particles
   real*8 :: lbox ! box length, square box
   real*8 :: mass(3) !mass and rotational inertia
   real*8,parameter :: pi=4d0*datan(1d0)
   integer :: tevery !frequency of frame for movie
   real*8 :: hbc !HB cut off energy
   real*8 :: tau ! Brendensen thermostat tau parameter
   real*8 :: pri !desired presssure
   real*8 :: taub ! Brendensen barostat tau parameter
   integer :: snap  !store snapshots and enery vs time
end module parametri

module MC
   real*8 :: h1,h2,pp,v1,v2,vh
   real*8 :: u,p,ckot
   integer :: kor1
   integer :: m
   real*8 :: k1,ktemp
end module MC

module potencial
    real*8, dimension(:,:), allocatable :: pot,mpot,tla,mtla,d,d2
	real*8, dimension(:,:,:), allocatable :: rrr
	real*8, dimension(:,:,:), allocatable :: sila
end module potencial

module korel
    real*8, dimension(10000) :: gr
	real*8, dimension(0:1000) :: kot
	real*8 :: interv
	integer :: ngr
end module korel

module hbvezi
	real*8, dimension(10) :: stv
	real*8, dimension(2) :: hba
end module hbvezi

module vkorel
   real*8, dimension(:), allocatable :: vr,omega,shear,flux
   real*8, dimension(:,:), allocatable :: vel0    ! velocities
   real*8, dimension(:), allocatable :: shear0 !off diagonal shear element
   real*8, dimension(:), allocatable :: flux0 !flux in direction x
end module vkorel

program dinamika
 use parametri
 use particles
 use potencial
 use MC
 use korel
 use vkorel
 use hbvezi
 implicit none
 character tekst
 integer :: i,j
 integer,dimension(8) :: values1,values2
 real*8 :: t,dgr(10000),time1,time2,cp,kapa,alfa,dkot(0:1000)
 real*8 :: term(25),dterm(25),stemp,ftemp,dtemp,g0
 real*8, dimension(:), allocatable :: avr,avr2,aomega,aomega2,ash,ash2,afl,afl2
 real*8, dimension(:,:), allocatable :: ggr,termo,kkot
 character*4 ftm,ftm2

 call random_init(.true., .true.)
 call cpu_time(time1)
 call date_and_time(VALUES=values1)
 open(1002,file='md.log')
 open(223,file='rezultati')
 open(1011,file='!md_termodinamika')
 write(1002,*)'start date and time YYYY, MM, DD, UTC, HH, MIN, SS'
 write(1002,'(8i5)') values1
 open(11,file='input_md3_rose_MB')
      read(11,*)tekst
      read(11,*)tekst,rh
      read(11,*)tekst,eh
      read(11,*)tekst,s1
      read(11,*)tekst,s2
      read(11,*)tekst,a1,a2
      read(11,*)tekst
      read(11,*)tekst,stemp,ftemp,dtemp
      read(11,*)tekst,g
	  read(11,*)tekst,pri
      read(11,*)tekst
      read(11,*)tekst,tstep
      read(11,*)tekst,nseries
      read(11,*)tekst,tsteps
      read(11,*)tekst,esteps
      read(11,*)tekst,N
	  read(11,*)tekst,interv
	  read(11,*)tekst,ngr
	  read(11,*)tekst,tevery
	  read(11,*)tekst,hbc
	  read(11,*)tekst,tau
	  read(11,*)tekst,taub
	  read(11,*)tekst,snap
 close(11)
 dt2=tstep/2d0
 d2t=tstep**2/2d0
 mass(1)=1d0
 mass(2)=1d0
 mass(3)=3d0*0.1d0*0.35d0**2 !mass of arm, distance of arm from center)
 lbox=sqrt(dfloat(N)/g)
 g0=g
 
 temp=stemp-dtemp
do while (temp.ge.ftemp)
 temp=temp+dtemp
 lbox=sqrt(dfloat(N)/g0)
 write (ftm2,'(a1,I3.3)') 't',int(temp*1000)

 write(1002,*)'box1=',lbox
 write(1002,*)'N=',N
 write(1002,*)'temp=',temp
 write(1002,*)'start density=',g0
 write(1002,*)'pressure=',pri
 write(1002,*)'time step=',tstep
 write(1002,*)'molecular dynamics of 2D MB particles'
 write(1002,*)'masses and inertia=', mass
 allocate( delci(N,3) )
 allocate( vel(N,3) )
 allocate( acc(N,3) )
 allocate( pot(N,N) )
 allocate( tla(N,N) )
 allocate( mpot(N,N) )
 allocate( mtla(N,N) )
 allocate( rrr(N,N,2) )
 allocate( sila(N,N,3) )
 allocate( d(N,N) )
 allocate( d2(N,N) )
 allocate( vr(0:tsteps) )
 allocate( vel0(N,3) )
 allocate( shear0(N) )
 allocate( flux0(N) )
 allocate( omega(0:tsteps) )
 allocate( shear(0:tsteps) )
 allocate( flux(0:tsteps) )
 allocate( avr(0:tsteps) )
 allocate( avr2(0:tsteps) )
 allocate( aomega(0:tsteps) )
 allocate( aomega2(0:tsteps) )
 allocate( ash(0:tsteps) )
 allocate( ash2(0:tsteps) )
 allocate( afl(0:tsteps) )
 allocate( afl2(0:tsteps) )
 allocate( termo(25,nseries) )
 allocate( ggr(10000,nseries) )
 allocate( kkot(0:1000,nseries) )
 
 call initp()
 write(1002,*)'positions initialized' 
 if (snap.eq.1) call snapshot('init'//ftm2)
 call initv()
 write(1002,*)'velocities initialized'
 
 call ekvmd(ftm2)
 write(1002,*)'MD ekvalibration successful'
 write(1002,*)'h1/N=',h1/dfloat(N)	
 write(1002,*)'v1/N=',v1/dfloat(N),dfloat(N)/v1
 write(1002,*)'p=',pp,pp*dfloat(N)/v1*temp	
 if (snap.eq.1) call snapshot('eqmd'//ftm2)
 if (snap.eq.1) then
 open(11,file='box'//ftm2)
   write(11,'(2f8.2)')-lbox/2d0,-lbox/2d0
   write(11,'(2f8.2)')-lbox/2d0,lbox/2d0
   write(11,'(2f8.2)')lbox/2d0,lbox/2d0
   write(11,'(2f8.2)')lbox/2d0,-lbox/2d0
   write(11,'(2f8.2)')-lbox/2d0,-lbox/2d0
 close(11)
 endif
 write(223,*)'MD'
 
 avr=0d0
 avr2=0d0
 aomega=0d0
 aomega2=0d0
 ash=0d0
 ash2=0d0
 afl=0d0
 afl2=0d0
 do i=1,nseries
  call md(i) 
  g=dfloat(N)/v1
  cp=(h2-h1**2)/temp**2/dfloat(N)
  kapa=(v2-v1**2)/temp/v1
  alfa=(vh-v1*h1)/temp**2/v1
  write(1002,*) i,'-th run'
  write(1002,*)'h1/N=',h1/dfloat(N)	
  write(1002,*) h1,h2,v1,v2,vh
  write(1002,*)'density=',g
  write(1002,*)'heat capacity=',cp
  write(1002,*)'expansion=',alfa
  write(1002,*)'compress=',kapa	 
  write(1002,*)'p=',pp,pp*g*temp
  write(1002,*)'Ekin=',k1  
  write(1002,*)'T=',ktemp
  write(1002,*)'aver cos=',ckot
  write(1002,*)'aver num HB=',hba(1)
  termo(13,i)=ckot
  termo(14,i)=stv(1)
  termo(15,i)=stv(2)
  termo(16,i)=stv(3)
  termo(17,i)=stv(4)
  termo(18,i)=hba(1)
  termo(19,i)=hba(2)
  write(223,'(9f15.8)')temp,g,h1/dfloat(N),h2/dfloat(N)**2,cp,pp,pp*g*temp,k1,ktemp
  write (ftm,'(a1,I2.2)') 'd',i
  if (snap.eq.1) call snapshot(trim(ftm)//trim(ftm2)) 
  ggr(:,i)=gr(:)  
  kkot(:,i)=kot(:)
  termo(1,i)=h1/dfloat(N)
  termo(2,i)=g
  termo(3,i)=kapa
  termo(4,i)=alfa
  termo(5,i)=cp
  termo(6,i)=pp*g*temp	
  termo(7,i)=ktemp

  do j=0,tsteps
   avr(j)=avr(j)+vr(j)
   aomega(j)=aomega(j)+omega(j)
   ash(j)=ash(j)+shear(j)
   afl(j)=afl(j)+flux(j)
   avr2(j)=avr2(j)+vr(j)**2
   aomega2(j)=aomega2(j)+omega(j)**2
   ash2(j)=ash2(j)+shear(j)**2
   afl2(j)=afl2(j)+flux(j)**2
  enddo

 enddo
 term=0d0	  
 dterm=0d0
 open(102,file='md_data'//ftm2)
 do i=1,19
 do j=1,nseries
   term(i)=termo(i,j)+term(i)
 enddo
 term(i)=term(i)/dfloat(nseries)
 do j=1,nseries
  dterm(i)=dterm(i)+(term(i)-termo(i,j))**2
 enddo
 dterm(i)=dsqrt(dterm(i)/dfloat(nseries))
 write(102,*)term(i),dterm(i),nseries
 enddo
 close(102)
 write(1011,'(21f15.8)')temp,pri,term(1),term(2),term(3),term(4),term(5),term(6),term(7),term(8),term(9),term(10),term(11), &
 & term(12),term(13),term(14),term(15),term(16),term(17),term(18),term(19)
 
  open(11,file='md_koor'//ftm2)
   gr=0d0
   dgr=0d0
   do i=1,10000
    do j=1,nseries
      gr(i)=gr(i)+ggr(i,j)
    enddo
    gr(i)=gr(i)/dfloat(nseries)
    do j=1,nseries
     dgr(i)=dgr(i)+(ggr(i,j)-gr(i))**2
    enddo
    dgr(i)=dsqrt(dgr(i)/dfloat(nseries))
	if ((dfloat(i)-0.5)*interv .lt. lbox) write(11,'(3f15.5)') (dfloat(i)-0.5)*interv,gr(i),dgr(i)
   enddo
  close(11)
  open(11,file='md_kot'//ftm2)
   kot=0d0
   dkot=0d0
   do i=0,1000
    do j=1,nseries
      kot(i)=kot(i)+kkot(i,j)
    enddo
    kot(i)=kot(i)/dfloat(nseries)
    do j=1,nseries
     dkot(i)=dkot(i)+(kkot(i,j)-kot(i))**2
    enddo
    dkot(i)=dsqrt(dkot(i)/dfloat(nseries))
	write(11,'(3f15.5)') dfloat(i)*pi/5d2,kot(i),dkot(i)
   enddo
  close(11)
  open(11,file='md_akoor'//ftm2)
  avr=avr/dfloat(nseries)
  avr2=avr2/dfloat(nseries)
  aomega=aomega/dfloat(nseries)
  aomega2=aomega2/dfloat(nseries)
  ash=ash/dfloat(nseries)
  ash2=ash2/dfloat(nseries)
  afl=afl/dfloat(nseries)
  afl2=afl2/dfloat(nseries)
  do i=0,tsteps
   write(11,'(13f15.5)') dfloat(i)*tstep,avr(i),dsqrt(avr2(i)-avr(i)**2),aomega(i),  &
   & dsqrt(aomega2(i)-aomega(i)**2),avr(i)/avr(0),aomega(i)/aomega(0),ash(i), &
   & dsqrt(ash2(i)-ash(i)**2),ash(i)/ash(0),afl(i),dsqrt(afl2(i)-afl(i)**2),afl(i)/afl(0)
  enddo
  close(11)

 deallocate( delci )
 deallocate( vel )
 deallocate( acc )
 deallocate( pot )
 deallocate( tla )
 deallocate( mpot )
 deallocate( mtla )
 deallocate( rrr )
 deallocate( sila )
 deallocate( d )
 deallocate( d2 )
 deallocate( vr )
 deallocate( vel0 )
 deallocate( shear0 )
 deallocate( flux0 )
 deallocate( omega )
 deallocate( shear )
 deallocate( flux )
 deallocate( avr )
 deallocate( avr2 )
 deallocate( aomega )
 deallocate( aomega2 )
 deallocate( ash )
 deallocate( ash2 )
 deallocate( afl )
 deallocate( afl2 )
 deallocate( termo )
 deallocate( ggr )
 deallocate( kkot )
 write(1002,*)'----------------------------------------'
 enddo
 call cpu_time(time2)
 call date_and_time(VALUES=values2)
 write(1002,*)'CPU simulation time=',time2-time1
 write(1002,*)'start and finish date and time YYYY, MM, DD, UTC, HH, MIN, SS'
 write(1002,'(8i5)') values1
 write(1002,'(8i5)') values2
 time2=(values2(8)-values1(8))/1000.+values2(7)-values1(7)
 time2=time2+(values2(6)-values1(6))*60
 time2=time2+(values2(5)-values1(5))*60*60
 time2=time2+(values2(3)-values1(3))*60*60*24
 write(1002,*)'Real simulation time=',time2
 close(1002)
 close(223)
 close(1011)
end program dinamika
	  
!--------------------------------------------------------------------
!     function minimum image
      function image(ax,bx,ll)
      real*8 :: image
      real*8 :: ax,bx,ll
      real*8 :: cx
      cx=bx-ax
      cx=cx-ll*dnint(cx/ll)
      image=cx
      end function image
 
!------------------------------------------------------------------
! procedure calculates interaction betwen particles----------------
!------------------------------------------------------------------
 subroutine interakcija()
 use particles
 use MC
 use potencial
 use parametri
 integer :: i,j
 real*8 :: ljpot,image,p1,x,y,r,mbpot,r2
 real*8 :: fx,fy,m1,m2

  pot=0d0
  mpot=0d0
  tla=0d0
  mtla=0d0
  rrr=0d0  
  d=0d0
  d2=0d0
  sila=0d0
  do i=1,N
   do j=i,N
    if (i .eq. j) then 
     pot(i,i)=0d0
     tla(i,i)=0d0
	 rrr(i,i,:)=0d0
	 d(i,i)=0d0
    else                                 
     x=image(delci(j,1),delci(i,1),lbox)
     y=image(delci(j,2),delci(i,2),lbox)  
     rrr(i,j,1)=x
	 rrr(i,j,2)=y
     rrr(j,i,1)=-x
	 rrr(j,i,2)=-y
     r2=x**2+y**2
	 r=dsqrt(r2)
	 d2(i,j)=r2
	 d2(j,i)=d2(i,j)
	 d(i,j)=r
	 d(j,i)=d(i,j)
     pot(i,j)=ljpot(r2,a1,a2,p1)
     tla(i,j)=p1
     pot(j,i)=pot(i,j)
     tla(j,i)=tla(i,j)
     sila(i,j,1)=-tla(i,j)*rrr(i,j,1)/d2(i,j)
     sila(i,j,2)=-tla(i,j)*rrr(i,j,2)/d2(i,j) 
	 mpot(i,j)=mbpot(r,x,y,delci(i,3),delci(j,3),fx,fy,m1,m2,p1)
	 mtla(i,j)=p1
	 sila(i,j,1)=sila(i,j,1)-mtla(i,j)*rrr(i,j,1)/d2(i,j)+fx
     sila(i,j,2)=sila(i,j,2)-mtla(i,j)*rrr(i,j,2)/d2(i,j)+fy
	 sila(i,j,3)=sila(i,j,3)+m1
     mpot(j,i)=mpot(i,j)
     mtla(j,i)=mtla(i,j)	 
     sila(j,i,1)=-sila(i,j,1)
     sila(j,i,2)=-sila(i,j,2)
	 sila(j,i,3)=sila(j,i,3)+m2
    end if  
   enddo
  enddo 
  u=sum(pot)/2d0+sum(mpot)/2d0
  p=sum(tla)/2d0+sum(mtla)/2d0
 end subroutine interakcija
 
!-------------------------------------------------------------------
! function that returns the LJ potencial, r is in reduced units-----
!-------------------------------------------------------------------
function ljpot(r,a1,a2,x1)
 implicit none
 real*8 :: r,ljpot,x,a1,a2,x1
 x=a2**2/r
 x=x**2*x
 ljpot=4d0*a1*x*(x-1d0)
 x1=-48d0*a1*x*(x-0.5d0)
end function ljpot

!------------------------------------------------------------------------
! Gauss function G(x,a,s)
!------------------------------------------------------------------------
function gau(x,a,s)
 implicit none
 real*8 :: gau,a,s,x
  gau=exp(-(x-a)**2/2d0/s**2)
end function gau
!------------------------------------------------------------------------
! derivative of Gauss function G(x,a,s)
!------------------------------------------------------------------------
function dgau(x,a,s)
 implicit none
 real*8 :: gau,a,s,x,dgau
  dgau=-(x-a)/s**2*gau(x,a,s)
end function dgau


!----------------------------------------------------------------------
! function returns MB potential between MB particles------
!----------------------------------------------------------------------
function mbpot(r,x1,y1,f11,f21,fx,fy,m1,m2,xx)
use parametri
implicit none
 real*8 :: r,x1,y1,f11,f21,dgau,fx,fy,m1,m2,df,dy1,dy2
 real*8 :: mbpot,x,gau,y,z,v1(3),v2(3),xx,dy,dy3
 real*8 :: aa1,aa2,bb1,bb2,f1(3),f2(3),dv1(3),dv2(3)
 integer :: i,j
 f1(1)=f11+2d0*pi/3d0
 f1(2)=f11-2d0*pi/3d0
 f1(3)=f11
 f2(1)=f21+2d0*pi/3d0
 f2(2)=f21-2d0*pi/3d0
 f2(3)=f21
 if (abs(r-rh)/s1 .lt. 4d0) then
  aa1=x1/r
  aa2=y1/r
  mbpot=eh*gau(r,rh,s1)
  xx=eh*dgau(r,rh,s1)*r
  m1=mbpot
  m2=mbpot
  fx=mbpot
  fy=mbpot
  y=0d0
  dy=0d0
  dy1=0d0
  dy2=0d0
  dy3=0d0
  do i=1,3
   bb1=cos(f1(i))
   bb2=sin(f1(i))
   v1(i)=aa1*bb1+aa2*bb2
   dv1(i)=-aa1*bb2+aa2*bb1
  enddo
  do j=1,3
   bb1=cos(f2(j))
   bb2=sin(f2(j))
   v2(j)=aa1*bb1+aa2*bb2
   dv2(j)=-aa1*bb2+aa2*bb1
  enddo
  do i=1,3
   x=(v1(i)-1d0)
   do j=1,3
    z=(v2(j)+1d0)
    if ((abs(x/s2) .lt. 4d0) .and. (abs(z/s2) .lt. 4d0)) then 
	 y=y+gau(x,0d0,s2)*gau(z,0d0,s2)
	 dy=dy+dgau(x,0d0,s2)*dv1(i)*gau(z,0d0,s2)
	 dy3=dy3+gau(x,0d0,s2)*dv2(j)*dgau(z,0d0,s2)
	 dy1=dy1+dgau(x,0d0,s2)*dv1(i)*gau(z,0d0,s2)*y1/r**2
	 dy1=dy1+gau(x,0d0,s2)*dv2(j)*dgau(z,0d0,s2)*y1/r**2
	 dy2=dy2-dgau(x,0d0,s2)*dv1(i)*gau(z,0d0,s2)*x1/r**2
	 dy2=dy2-gau(x,0d0,s2)*dv2(j)*dgau(z,0d0,s2)*x1/r**2
	endif
   enddo
  enddo
  mbpot=mbpot*y
  xx=xx*y
  m1=-m1*dy
  fx=-fx*dy1
  fy=-fy*dy2
  m2=-m2*dy3
 else
  mbpot=0d0
  xx=0d0
  fx=0d0
  fy=0d0
  m1=0d0
  m2=0d0
 end if
end function mbpot  





!------------------------------------------------------------------
! procedure moves one particle and evaluate if we accept this move
!------------------------------------------------------------------
subroutine snapshot(tekst)
 use particles
 use parametri
 use potencial
 implicit none
 character(8) :: tekst
 integer :: i,j
 real*8 :: image,x,y
 open(11,file='slika_'//tekst)
 open(16,file='roke_'//tekst) 
 do i=1,N
  write(11,'(3f15.7)')delci(i,1),delci(i,2),a2/2d0
  write(16,'(2f8.2)')delci(i,1),delci(i,2)
  write(16,'(2f8.2)')delci(i,1)+0.5d0*dcos(delci(i,3)+3d0*pi/3d0),delci(i,2)+0.5d0*dsin(delci(i,3)+3d0*pi/3d0)
  write(16,'(2f8.2)')
  write(16,'(2f8.2)')delci(i,1),delci(i,2)
  write(16,'(2f8.2)')delci(i,1)+0.5d0*dcos(delci(i,3)-1d0*pi/3d0),delci(i,2)+0.5d0*dsin(delci(i,3)-1d0*pi/3d0)
  write(16,'(2f8.2)')
  write(16,'(2f8.2)')delci(i,1),delci(i,2)
  write(16,'(2f8.2)')delci(i,1)+0.5d0*dcos(delci(i,3)+1d0*pi/3d0),delci(i,2)+0.5d0*dsin(delci(i,3)+1d0*pi/3d0)
  write(16,'(2f8.2)')
 enddo
 close(11)
 close(16) 
 open(15,FILE='hb_vezi_'//tekst)
 do i=1,N
  do j=i,N
   if (mpot(i,j).lt.-hbc) then
    x=image(delci(i,1),delci(j,1),lbox)
    y=image(delci(i,2),delci(j,2),lbox)
    write(15,'(2f8.2)')delci(i,1),delci(i,2)
    write(15,'(2f8.2)')delci(i,1)+x,delci(i,2)+y	
    write(15,*)	
   endif
  enddo
 enddo
 close(15)	
 end subroutine snapshot
 
!------------------------------------------------------------------
! initialize positions  -------------------------------------------
!------------------------------------------------------------------
 subroutine initp()
 use particles
 use parametri
 implicit none
 integer :: i,j,k
 real*8 :: image
 real*8:: z(3)
 call random_number(z)
 delci(1,1)=(z(1)-0.5d0)*lbox
 delci(1,2)=(z(2)-0.5d0)*lbox
 delci(1,3)=z(3)*2d0*pi
 do i=2,N
   j=0
   do while (j .lt. 0.5)
    call random_number(z)
    delci(i,1)=(z(1)-0.5d0)*lbox
    delci(i,2)=(z(2)-0.5d0)*lbox
    delci(i,3)=z(3)*2d0*pi
    j=1
    do k=1,i-1
      if (image(delci(i,1),delci(k,1),lbox)**2+image(delci(i,2),delci(k,2),lbox)**2 .le. 0.4) j=0
    enddo
   enddo
 enddo   
 end subroutine initp
 
!------------------------------------------------------------------
! initialize velocities -------------------------------------------
!------------------------------------------------------------------
 subroutine initv()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 :: ran3,s,x1,x2,z1,z2,sigma(3),sumv(3)
 real*8:: z(4)
 sigma(:)=dsqrt(temp/mass(:))
 do i=1, N
  call random_number(z)
  vel(i,1)=sqrt(-2*log(z(1)))*cos(2*pi*z(2))*sigma(1)
  vel(i,2)=sqrt(-2*log(z(1)))*sin(2*pi*z(2))*sigma(2)
  vel(i,3)=sqrt(-2*log(z(3)))*cos(2*pi*z(4))*sigma(3)
 enddo
 sumv(1)=sum(vel(:,1))
 sumv(2)=sum(vel(:,2))
 sumv(3)=sum(vel(:,3))
 vel(:,1)=vel(:,1)-sumv(1)/dfloat(N)
 vel(:,2)=vel(:,2)-sumv(2)/dfloat(N)
 vel(:,3)=vel(:,3)-sumv(3)/dfloat(N)
 end subroutine initv
 
!-----------------------------------------------------------------------
! ekvalibration of particles by MD method ------------------------------
!-----------------------------------------------------------------------
 subroutine ekvmd(ftm)
 use parametri
 use particles
 use MC
 implicit none	
 integer :: i 
 real*8 :: t
 character(4) ftm
 if (snap.eq.1) open(10,FILE='ekvmd'//ftm)
 t=0d0
 h1=0d0
 pp=0d0
 v1=0d0
 m=0
 do i=1,esteps
 call forces()
 call move()
 call temperature1()
 call pressure2()
 call com_vel()
 t=t+tstep
    if (snap.eq.1) write(10,'(4f20.7)') t,tc,(dfloat(N)*tc-p/2.0d0)/lbox**2
    h1=h1+u+pri*lbox**2
	v1=v1+lbox**2
    pp=pp+p
    m=m+1
 enddo
 h1=h1/dfloat(m)
 v1=v1/dfloat(m)
 pp=pp/dfloat(m)
 if (snap.eq.1) write(10,*)'----------------------------------------'
 if (snap.eq.1) write(10,'(a10,2f20.7)') 'povprecje',h1,pp
 pp=1-pp/2d0/temp/dfloat(N)
 if (snap.eq.1) close(10)
 end subroutine ekvmd
 
!-----------------------------------------------------------------------
! forces calculations --------------------------------------------------
!-----------------------------------------------------------------------
 subroutine forces()
 use particles
 use potencial
 use parametri
 implicit none
 integer :: i
 real*8 :: nav
 call interakcija()
 do i=1,N                      
  acc(i,1)=(sum(sila(i,:,1)))/mass(1) 
  acc(i,2)=(sum(sila(i,:,2)))/mass(2)
  acc(i,3)=(sum(sila(i,:,3)))/mass(3) 
 enddo
 end subroutine forces
 
!-----------------------------------------------------------------------
! Verlet velocities - --------------------------------------------------
!-----------------------------------------------------------------------
 subroutine move()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 image
 do i=1,N
  vel(i,1)=vel(i,1)+dt2*acc(i,1)  !previous step
  vel(i,2)=vel(i,2)+dt2*acc(i,2)  !previous step
  vel(i,3)=vel(i,3)+dt2*acc(i,3)  !previous step
  ! current step
  delci(i,1)=delci(i,1)+vel(i,1)*tstep+acc(i,1)*d2t
  delci(i,2)=delci(i,2)+vel(i,2)*tstep+acc(i,2)*d2t
  delci(i,3)=delci(i,3)+vel(i,3)*tstep+acc(i,3)*d2t
  vel(i,1)=vel(i,1)+dt2*acc(i,1)
  vel(i,2)=vel(i,2)+dt2*acc(i,2)
  vel(i,3)=vel(i,3)+dt2*acc(i,3)
  delci(i,1)=image(0d0,delci(i,1),lbox)
  delci(i,2)=image(0d0,delci(i,2),lbox)
  delci(i,3)=image(0d0,delci(i,3),2*pi)
 enddo
 end subroutine move
 
!-----------------------------------------------------------------------
! velocities rescaling to maintain temeprature -------------------------
!-----------------------------------------------------------------------
 subroutine temperature1()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 rescale
 ekin=0d0
 do i=1,N
  ekin=ekin+(vel(i,1)**2*mass(1)+vel(i,2)**2*mass(2)+vel(i,3)**2*mass(3))/2d0
 enddo
 ekin=ekin/dfloat(N)
 tc=2d0*ekin/3d0
 rescale=dsqrt(tc/temp)
 vel=vel/rescale
 end subroutine temperature1
 
!-----------------------------------------------------------------------
! Berendsen thermostat--------------------------------------------------
!-----------------------------------------------------------------------
 subroutine temperature2()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 rescale
 ekin=0d0
 do i=1,N
  ekin=ekin+(vel(i,1)**2*mass(1)+vel(i,2)**2*mass(2)+vel(i,3)**2*mass(3))/2d0
 enddo
 ekin=ekin/dfloat(N)
 tc=2d0*ekin/3d0
 rescale=dsqrt(1d0+tstep/tau*(temp/tc-1d0))
 vel=vel*rescale
 end subroutine temperature2
 
!-----------------------------------------------------------------------
! canonical sampling through velocity rescaling (CSVR) thermostat-------
!-----------------------------------------------------------------------
 subroutine temperature3()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 rescale,z1(3*N),z2(3*N),norm(3*N)
 call random_number(z1)
 call random_number(z2)
 norm=sqrt(-2d0*log(z1))*cos(2*pi*z2)
 ekin=0d0
 do i=1,N
  ekin=ekin+(vel(i,1)**2*mass(1)+vel(i,2)**2*mass(2)+vel(i,3)**2*mass(3))/2d0
 enddo
 tc=2d0*ekin/3d0/dfloat(N)
! rescale=dsqrt(1d0+tstep/tau*(temp/tc-1d0))
 rescale=exp(-tstep/tau)+temp/2d0/ekin*(1d0-exp(-tstep/tau))*sum(norm**2)
 rescale=rescale+2d0*exp(-tstep/2d0/tau)*sqrt(temp/2d0/ekin*(1d0-exp(-tstep/tau)))*norm(1)
 vel=vel*sqrt(rescale)
 end subroutine temperature3
 
!-----------------------------------------------------------------------
! Remove COM velocity---------------------------------------------------
!-----------------------------------------------------------------------
 subroutine com_vel()
 use particles
 use parametri
 implicit none
 integer :: i
 real*8 comv(3),comm(3)
 comv=0d0
 comm=0d0
 do i=1,N
  comv(:)=comv(:)+vel(i,:)*mass(:)
  comm(:)=comm(:)+mass(:)
 enddo
 comv(:)=comv(:)/comm(:)
 do i=1,N
 vel(i,:)=vel(i,:)-comv(:)
 enddo
 end subroutine com_vel
 
 
!-----------------------------------------------------------------------
! Berendsen barostat--------------------------------------------------
!-----------------------------------------------------------------------
 subroutine pressure()
 use particles
 use parametri
 use MC
 implicit none
 integer :: i
 real*8 rescale
 rescale=dsqrt(1.0d0+tstep/taub*((dfloat(N)*tc-p/2.0d0)/lbox**2-pri))
 lbox=lbox*rescale
 delci=delci*rescale
 end subroutine pressure
 
 
!-----------------------------------------------------------------------
! SCR barostat--------------------------------------------------
!-----------------------------------------------------------------------
 subroutine pressure2()
 use particles
 use parametri
 use MC
 implicit none
 integer :: i
 real*8 rescale,z1(N),z2(N),norm(N),sumi
 call random_number(z1)
 call random_number(z2)
 norm=sqrt(-2d0*log(z1))*cos(2*pi*z2)
 sumi=sqrt(2d0*tc*tstep/lbox**2/taub)*norm(1)
 rescale=dsqrt(1.0d0+tstep/taub*((dfloat(N)*tc-p/2.0d0)/lbox**2-pri)+sumi)
 lbox=lbox*rescale
 delci=delci*rescale
 end subroutine pressure2

 
!-----------------------------------------------------------------------
! sampling of particles by MD method ------------------------------
!-----------------------------------------------------------------------
 subroutine md(ii)
 use parametri
 use particles
 use MC
 use korel
 use vkorel
 use potencial
 use hbvezi
 implicit none	
 integer :: i,ii,mm
 real*8 :: t,c,ei
 character*4 ftm,ftm2
 character*8 ftm3
 write (ftm,'(I3.3)') ii
 write (ftm2,'(a1,I3.3)') 't',int(temp*1000)
 if (snap.eq.1) open(10,FILE='md_sam'//trim(ftm)//trim(ftm2))
 t=0d0
 h1=0d0
 h2=0d0
 k1=0d0
 v1=0d0
 v2=0d0
 vh=0d0
 ktemp=0d0
 pp=0d0
 m=0
 stv=0d0
 hba=0d0
 gr=0d0
 vel0=vel
 ckot=0d0
 do i=1,N
  shear0(i)=mass(1)*vel(i,1)*vel(i,2)+delci(i,1)*mass(1)*acc(i,2)
  ei=0.5d0*mass(1)*(vel(i,1)**2+vel(i,2)**2)+(sum(pot(i,:))+sum(mpot(i,:)))/2d0
  flux0(i)=ei*vel(i,1)
  do mm=1,N
   flux0(i)=flux0(i)+(delci(i,1)-delci(mm,1))*(sila(i,mm,1)*vel(i,1)+sila(i,mm,2)*vel(i,2))
  enddo
 enddo
 call akoor(0)
 mm=0
 do i=1,tsteps
 call forces()
 call move()
 call temperature3()
 call pressure2()
 call com_vel()
 if (mod(i,ngr).eq.0) then
  call koor()
  call vezi()
  mm=mm+1
  write (ftm3,'(I2.2,I6.6)') ii,i/tevery
 endif
 if ((mod(i,tevery).eq.0).and.(ii.eq.1)) then
  write (ftm3,'(I2.2,I6.6)') ii,i/tevery
  call frame(ftm3//ftm2)
 endif
 t=t+tstep
    if (snap.eq.1) write(10,'(5f20.7)') t,tc,(dfloat(N)*tc-p/2.0d0)/lbox**2
    h1=h1+u+pri*lbox**2
	h2=h2+(u+pri*lbox**2)**2
	v1=v1+lbox**2
	v2=v2+lbox**4
	vh=vh+lbox**2*(u+pri*lbox**2)
    pp=pp+p
	k1=k1+ekin
	ktemp=ktemp+tc
    m=m+1
	call akoor(i)
 enddo
 h1=h1/dfloat(m)
 h2=h2/dfloat(m)
 v1=v1/dfloat(m)
 v2=v2/dfloat(m)
 vh=vh/dfloat(m)
 pp=pp/dfloat(m)
 k1=k1/dfloat(m)
 ckot=ckot/dfloat(N)/dfloat(mm)
 stv=stv/dfloat(N)/dfloat(mm)
 hba=hba/dfloat(N)/dfloat(mm)
 ktemp=ktemp/dfloat(m)
 if (snap.eq.1) write(10,*)'----------------------------------------'
 if (snap.eq.1) write(10,'(a10,2f20.7)') 'povprecje',h1,pp
 pp=1-pp/2d0/temp/dfloat(N)
 if (snap.eq.1) close(10)
 do i=1,10000
  c=pi*(2*i-1)*interv**2
  gr(i)=gr(i)/dfloat(N-1)/dfloat(mm)/dfloat(N)*v1/c
 enddo
 do i=0,1000
  kot(i)=kot(i)/dfloat(mm)/dfloat(N)/pi*500
 enddo
 end subroutine md
 
!-----------------------------------------------------------------------
! correlation of particles ---------------------------------------------
!-----------------------------------------------------------------------
 subroutine koor()
 use potencial
 use parametri
 use korel
 use MC
 use particles
 implicit none	
 integer :: i,j,k
 real*8 :: rr,image
  do i=1,N-1
  ckot=ckot+dcos(delci(i,3))
  rr=image(delci(i,3),0d0,2d0*pi)+pi
   k=int(rr/pi*500)
   if (k.lt.1000) then
     kot(k)=kot(k)+1d0
   endif  
  do j=i+1,N
   rr=d(i,j)
   k=int(rr/interv)
   if (k.lt.10000) then
     gr(k+1)=gr(k+1)+2d0
   endif
  enddo
  enddo
  ckot=ckot+dcos(delci(N,3))
  rr=image(delci(N,3),0d0,2d0*pi)+pi
   k=int(rr/pi*500)
   if (k.lt.1000) then
     kot(k)=kot(k)+1d0
   endif  
 end subroutine koor
 
!-----------------------------------------------------------------------
! correlation of particles ---------------------------------------------
!-----------------------------------------------------------------------
 subroutine akoor(i)
 use particles
 use parametri
 use vkorel
 use potencial
 implicit none	
 integer :: i,j,k
 real*8 :: sh,ei,fl
  do j=1,N
   omega(i)=omega(i)+vel(j,3)*vel0(j,3)
   vr(i)=vr(i)+vel(j,1)*vel0(j,1)+vel(j,2)*vel0(j,2)
   sh=mass(1)*vel(j,1)*vel(j,2)+delci(j,1)*mass(1)*acc(j,2)
   shear(i)=shear(i)+shear0(j)*sh
   ei=0.5d0*mass(1)*(vel(j,1)**2+vel(j,2)**2)+(sum(pot(j,:))+sum(mpot(j,:)))/2d0
   fl=ei*vel(j,1)
   do k=1,N
    fl=fl+(delci(j,1)-delci(k,1))*(sila(j,k,1)*vel(j,1)+sila(j,k,2)*vel(j,2))
   enddo
   flux(i)=flux(i)+flux0(j)*fl
  enddo
  omega(i)=omega(i)/dfloat(N)
  vr(i)=vr(i)/dfloat(N)
  shear(i)=shear(i)/dfloat(N)
  flux(i)=flux(i)/dfloat(N)
 end subroutine akoor
 
!------------------------------------------------------------------
! figures for movie
!------------------------------------------------------------------
 subroutine frame(tekst)
 use particles
 use parametri
 implicit none
 character(12) :: tekst
 integer :: i
 open(11,file='f'//tekst)
 do i=1,N
  write(11,'(2f8.2)')delci(i,1),delci(i,2)
 enddo
 close(11)
 end subroutine frame
 
 
!-------------------------------------------------------------------
! subroutine that evaluate hb statistics
      subroutine vezi()
	  use parametri
      use potencial
	  use hbvezi
      implicit none
      integer i,j,k


      do i=1,N
       k=0
       do j=1,N
         if (mpot(i,j).le.-hbc) k=k+1
       enddo
      if (k.le.9) stv(k+1)=stv(k+1)+1d0
      enddo
      do i=1,N
       do j=1,N
         if (mpot(i,j).le.-hbc) hba(1)=hba(1)+1
       enddo
      enddo
	  
      do i=1,N
       do j=1,N
         if (pot(i,j).le.-hbc) hba(2)=hba(2)+1
       enddo
      enddo
      return
      end subroutine vezi	
	  
!-------------------------------------------------------------------