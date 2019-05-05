program main
      implicit none 
      integer::i,j,n=108
      !real(kind=8)::fx(200),fy(200),fz(200) 
      !real(8)::v1,v2,f,k=0.0000001
      real(8),allocatable::a(:,:),rx(:),ry(:),rz(:),fx(:),fy(:),fz(:)

      allocate (a(n,n),rx(n),ry(n),rz(n),fx(n),fy(n),fz(n))
      
      open(1,file='argon.init')
     
      do  i=1,n
        read(1,*) (a(i,j),j=1,3)
        write(10,*)(a(i,j),j=1,3)
        rx(i) =a(i,1)
        ry(i) =a(i,2)
        rz(i) =a(i,3)
      end do
      call calc_force(n,rx,ry,rz,fx,fy,fz)
      
             
      
end program main

subroutine calc_force(n,rx,ry,rz,fx,fy,fz)
 implicit none
 integer::i,j,n
 real(8)::rx(200),ry(200),rz(200),boxl,sigsq,epslon,rxij,ryij,rzij,rijsq,sr2,sr6,sr12,rcutoffsq,rcutoff,sig,kb
 real(8)::fxij,fyij,fzij,fpr,fx(200),fy(200),fz(200)
 !real(8),allocatable::a(:,:)!,fx(:),fy(:),fz(:)
 n=108
 !do  i=1,n
  !write(15,*)rx(i),ry(i),rz(i)  
 !end do
 
 boxl =17.158d0
 sig = 3.405d0
 sigsq = sig * sig
 kb = 1.38E-23
 epslon = 120.0d0 !* kb

 rcutoff = 0.50d0 * boxl
 rcutoffsq = rcutoff * rcutoff

 do i=1,n
  rx(i) = rx(i) - (boxl * anint(rx(i)/boxl))
  ry(i) = ry(i) - (boxl * anint(ry(i)/boxl))
  rz(i) = rz(i) - (boxl * anint(rz(i)/boxl))
 end do

 do i=1,n
  fx(i)=0.0d0; fy(i)=0.0d0; fz(i)=0.0d0
 end do

 
 do i=1,n-1
  do j=i+1,n
   rxij = rx(i) - rx(j)
   ryij = ry(i) - ry(j)
   rzij = rz(i) - rz(j)
 
   rxij = rxij - (boxl * anint(rxij/boxl))
   ryij = ryij - (boxl * anint(ryij/boxl))
   rzij = rzij - (boxl * anint(rzij/boxl))
  
 
   rijsq = rxij**2 + ryij**2 + rzij**2
   sr2 = sigsq/rijsq
   sr6 = sr2 * sr2 * sr2
   sr12 = sr6 * sr6

   if(rijsq<rcutoffsq) then    
    fpr = (48.0d0 * epslon * (sr12 - (0.5d0*sr6)))/rijsq
    fxij = fpr * rxij
    fyij = fpr * ryij
    fzij = fpr * rzij

    fx(i) = fx(i) + fxij
    fy(i) = fy(i) + fyij     
    fz(i) = fz(i) + fzij
    fx(j) = fx(j) - fxij
    fy(j) = fy(j) - fyij     
    fz(j) = fz(j) - fzij
    
   end if
  end do  
  
   write(22,*)fx(i),fy(i),fz(i)              
      
 end do 
end subroutine calc_force
