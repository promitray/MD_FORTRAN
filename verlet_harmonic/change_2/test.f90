                IMPLICIT NONE 
                REAL(8)::x,vel,F,m,H
                INTEGER::ixx
                INTEGER,PARAMETER::itime = 100000
	        REAL,PARAMETER::xposeq = 0.00d0, val=1.000d0,delta=0.002d0
		
	        x = 0.00d0
		vel  = 1.00d0
                F= -val*(x-xposeq)
                m=1.00d0
                ixx=0
                
                write(16,*) ixx*delta,x
                write(18,*) ixx*delta,vel
                H = 0.5d0 *(vel**2) + 0.5d0 * val *(x**2)
                write(20,*) x,vel
                write(24,*) ixx*delta,H
              
		do ixx = 1, itime                
		  vel = vel + (0.5d0*(F/m)*delta)
                  x=x + (vel*delta)
                  F = -val*(x-xposeq)
                  vel = vel + (0.5d0*(F/m)*delta)     
                  write(16,*) ixx*delta,x
                  write(18,*) ixx*delta,vel
                  H = 0.5d0 *(vel**2) + 0.5d0 * val *(x**2)
                  write(20,*) x,vel
                  write(24,*) ixx*delta,H
		end do
		stop
		end
