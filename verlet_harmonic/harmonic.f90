
		
		parameter(itime = 10000,delta=0.02)
		parameter(xposeq = 0.0, val=1.000)
		
		


		x = 0.00
		vel  = 1.00
                F= -val*(x-xposeq)
                m=1.00
                
                !write(16,*) 0,x
                !write(18,*) 0,vel
                H = 0.5 *(vel**2) + 0.5 * val *(x**2)
                !write(20,*) x,vel
                write(24,*) 0,H
              


		
		do ixx = 1, itime 


               
		vel = vel + (0.5*(F/m)*delta)
                x=x + (vel*delta)
                F = -val*(x-xposeq)
                vel = vel + (0.5*(F/m)*delta)
        
           
		

		!write(16,*) ixx*delta,x
                !write(18,*) ixx*delta,vel
                H = 0.5 *(vel**2) + 0.5 * val *(x**2)
                !write(20,*) x,vel
                write(24,*) ixx*delta,H

		end do 


		stop
		end
