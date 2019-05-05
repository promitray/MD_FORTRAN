                IMPLICIT NONE 
                REAL::x,vel,F,m,H
                INTEGER::ixx
                INTEGER,PARAMETER::itime = 10000
		REAL,PARAMETER::xposeq = 0.0, val=1.000,delta=0.02
		
		


		x = 0.00
		vel  = 1
                F= -val*(x-xposeq)
                m=1.00
                ixx=0
                
                !write(16,*) ixx*delta,x
                !write(18,*) ixx*delta,vel
                H = 0.5 *(vel**2) + 0.5 * val *(x**2)
                !write(20,*) x,vel
                write(28,*) ixx*delta,H
              


		
		do ixx = 1, itime 


               
		vel = vel + (0.5*(F/m)*delta)
                x=x + (vel*delta)
                F = -val*(x-xposeq)
                vel = vel + (0.5*(F/m)*delta)
        
           
		

		!write(16,*) ixx*delta,x
                !write(18,*) ixx*delta,vel
                H = 0.5 *(vel**2) + 0.5 * val *(x**2)
                !write(20,*) x,vel
                write(28,*) ixx*delta,H

		end do 


		stop
		end
