program potential
 implicit none
 integer::i,j,n =108
 real(kind=8)::boxl,sigsq,epslon,sr2,sr6,sr12,rcutoff,v,sig,kb
 real (kind=8) :: rij
 real (kind=8), allocatable :: r(:,:), rijvec(:)

 open(1,file='argon.init')


 allocate (r(n,3), rijvec(3)) 

 do  i=1,n
   read(1,*)r(i,:)
 end do

 boxl =17.158d0
 sig = 3.405d0
 sigsq = sig * sig
 kb = 1.38E-23
 epslon = 120.0d0 * kb

 rcutoff = 0.50d0 * boxl


 v=0.0d0
 do i=1,n-1
   do j=i+1,n
      rijvec(:) = r(i,:) - r(j,:)
      rijvec(:) = rijvec(:) - (boxl*anint(rijvec(:)/boxl))
      rij = dsqrt(dot_product(rijvec(:),rijvec(:)))
      !write(11,*)rijvec(:)
      !write(12,*)rij
            
      sr2 = (sigsq/(rij*rij))
      sr6 = sr2 * sr2 * sr2
      sr12 = sr6 * sr6
        if(rij<rcutoff) then
           v = v + sr12 -sr6
        endif

   end do 
 end do 



 v = 4.0d0*epslon*v
 write(*,*),v

end program potential
