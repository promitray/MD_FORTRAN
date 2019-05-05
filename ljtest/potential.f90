program potential
 implicit none
 integer::i,j,n =108
 real(8)::rx(200),ry(200),rz(200),boxl,sigsq,epslon,rxij,ryij,rzij,rijsq,sr2,sr6,sr12,rcutoffsq,rcutoff,v,sig,kb,a

 dimension a(200,200)
 open(1,file='argon2.init')
 do  i=1,n
  read(1,*) (a(i,j),j=1,3)
  write(10,*)(a(i,j),j=1,3)
  rx(i) =a(i,1)
  ry(i) =a(i,2)
  rz(i) =a(i,3)
   
 end do

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

 v=0.0d0
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
    v = v + sr12 -sr6              
   endif
  end do
 end do 
 v = 4.0d0*epslon*v
 write(*,*)v

 end
