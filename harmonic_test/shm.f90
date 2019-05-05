c		***** program for simple harmonic motion ****
		
		parameter(iobs = 1000, itime = 1000,delta=0.0002)
		parameter(xposeq = 0.0, val=1.000)
c		dimension xpos(0:itime*iobs) 

		
c		*** initial position ***

		xpos1 = 0.0
		xpos2 = 0.01

c                **************************

c		*** Verlet's dynamics ***

		do iyy =1,iobs

		do ixx = 1, itime 

c               **** force calculation ***
                force = - val*(xpos2-xposeq)
		xpos3 = 2.0*xpos2 - xpos1 + force*(delta**2)
		xpos1 = xpos2 
		xpos2 = xpos3
		end do

		write(12,*) iyy*itime*delta,xpos3

		end do 


		stop
		end
