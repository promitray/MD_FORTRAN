program main
      implicit none 
      integer::i,j,n=108,k,itime = 10000
      real(8)::v,pot_lj,KE,m,delta =0.005d0
      real(8),allocatable::a(:,:),rx(:),ry(:),rz(:),fx(:),fy(:),fz(:),vx(:),vy(:),vz(:)
      m=39.948

      allocate (a(n,n),rx(n),ry(n),rz(n),fx(n),fy(n),fz(n),vx(n),vy(n),vz(n))
      
      open(1,file='argon.init')
     
      do  i=1,n
        read(1,*) (a(i,j),j=1,3)
        rx(i) =a(i,1)
        ry(i) =a(i,2)
        rz(i) =a(i,3)
      end do
      call calc_force(n,rx,ry,rz,fx,fy,fz)
      call calc_pot(n,rx,ry,rz,v)
      do i=1,n
       vx(i)=0.0d0; vy(i)=0.0d0; vz(i)=0.0d0
      end do
      KE = 0.0d0
      k=0
      do i=1,n
        KE = KE + 0.5*m*((vx(i)*vx(i)) + (vy(i)*vy(i)) + (vz(i)*vz(i)))
        write(13,*)k*delta,i,rx(i),ry(i),rz(i)
        write(14,*)k*delta,i,vx(i),vy(i),vz(i)
      end do 
      write(47,*)k*delta,v,KE,v+KE

     !!main loop
      do k=1,itime
       KE=0.00d0
       do i=1,n
        vx(i) = vx(i) + 0.5 * (fx(i)/m) * delta
        vy(i) = vy(i) + 0.5 * (fy(i)/m) * delta
        vz(i) = vz(i) + 0.5 * (fz(i)/m) * delta
        
         
        rx(i) = rx(i) + (vx(i)*delta)
        ry(i) = ry(i) + (vy(i)*delta)
        rz(i) = rz(i) + (vz(i)*delta)
       end do
 
       call calc_force(n,rx,ry,rz,fx,fy,fz)
       do i=1,n
        vx(i) = vx(i) +  0.5 * (fx(i)/m) * delta
        vy(i) = vy(i) +  0.5 * (fy(i)/m) * delta
        vz(i) = vz(i) +  0.5 * (fz(i)/m) * delta
        
        write(13,*)k*delta,i,rx(i),ry(i),rz(i)
        write(14,*)k*delta,i,vx(i),vy(i),vz(i)
        KE= KE + 0.5 *m*((vx(i)*vx(i))+ (vy(i)*vy(i))+(vz(i)*vz(i)))
       end do
       call calc_pot(n,rx,ry,rz,v)
       write(47,*)k*delta,v,KE,v+KE
      end do 
 
           
      
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
  
   !write(22,*)fx(i),fy(i),fz(i)              
      
 end do 
end subroutine calc_force
