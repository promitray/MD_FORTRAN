program potential_calc
 implicit none
 real(8)::v(500),sig,sigsq,epslon,sr2,sr6,sr12,r
 integer::i
 sig = 3.405d0
 sigsq = sig * sig
 epslon = 120.0d0

 do i=20,150
   r = 0.10d0*i
   sr2 = (sigsq/(r**2))
   sr6 = sr2 * sr2 * sr2
   sr12 = sr6 * sr6
   v= 4*epslon*(sr12 -sr6)
   write(23,*)r,v
 end do 
end program potential_calc
