program main
      implicit none 
      integer::i,j,n=108
      real(kind=8)::pot_lj 
      real(8)::v,v1,v2,f,k=0.0000001
      real(8),allocatable::a(:,:),rx(:),ry(:),rz(:)

      allocate (a(n,n),rx(n),ry(n),rz(n))
      
      open(1,file='argon.init')
     
     do  i=1,n
        read(1,*) (a(i,j),j=1,3)
        !write(10,*)(a(i,j),j=1,3)
        rx(i) =a(i,1)
        ry(i) =a(i,2)
        rz(i) =a(i,3)
     end do
     call calc_pot(n,rx,ry,rz,v)
     write(*,*)v
     rx(2)=rx(2)+k
     call calc_pot(n,rx,ry,rz,v1)
     write(*,*)v1
     rx(2)=rx(2)-(2*k)
     call calc_pot(n,rx,ry,rz,v2)
     write(*,*)v2
     f = ((v2 -v1)/(2*k))
     write(*,*)'force =',f
end program main
 

subroutine calc_pot(n,rx,ry,rz,pot_lj)
     implicit none 
     integer ::n
     real(kind=8) :: rx (200),ry(200),rz(200)
     real(kind=8) :: pot_lj
     integer::i,j 
     real(kind=8) :: boxl,sigsq,epslon
     real(kind=8) :: rxij,ryij,rzij,rijsq
     real(kind=8) :: sr2,sr6,sr12
     real(kind=8) :: rcutoffsq,rcutoff
     real(kind=8) :: sig,kb,a
     n=108 
     !do i= 1,n
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
        !write(3,*)rx(i),ry(i),rz(i)
     end do
     
     
     pot_lj=0.0d0
     do i=1,n-1
      do j=i+1,n
           rxij = rx(i) - rx(j)
           ryij = ry(i) - ry(j)
           rzij = rz(i) - rz(j)
           
           rxij = rxij - (boxl * anint(rxij/boxl))
           ryij = ryij - (boxl * anint(ryij/boxl))
           rzij = rzij - (boxl * anint(rzij/boxl))
           !write(18,*)rxij,ryij,rzij
           
           
           rijsq = rxij**2 + ryij**2 + rzij**2
           !write(19,*)rijsq
           sr2 = sigsq/rijsq
           sr6 = sr2 * sr2 * sr2
           sr12 = sr6 * sr6
        if(rijsq<rcutoffsq) then
          pot_lj= pot_lj + sr12 -sr6              
        endif
     end do
   end do 
    pot_lj = 4.0d0*epslon*pot_lj
end subroutine calc_pot
