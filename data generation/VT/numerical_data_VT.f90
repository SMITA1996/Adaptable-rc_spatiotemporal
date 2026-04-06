	implicit none
	integer           ::nt,i,j,kk,t0
	integer, parameter:: nx=40
	integer, parameter:: ny=40
	integer, parameter::step=100
	
	double precision, parameter::dt=1d-3
	double precision, parameter::dx=1d-1
	double precision, parameter::dy=1d-1
	
	double precision, parameter::a=10d-2
	double precision, parameter::b=1d0
	double precision, parameter::w=1d-1
	double precision, parameter::sigma=1d0	!w**2/2d0		!1.0
	double precision, parameter::D=1d-3
	double precision, parameter::twopi=8.0d00*datan(1.0d00)
	double precision, parameter::D0=125d-3       !noise amplitude
	double precision, parameter::tau=1d-1
	double precision, parameter::lambda=1d0/tau	
	double precision	    ::v(nx,ny),v2(nx,ny),u1(nx,ny),c,r5,r6,corr,E,h,cc
	double precision	    ::rxv,dfv,x,y,t,r1,r2,r3,r4,u,g,u3,suma,meand
	
	
	open(2000,file='data_mean.dat')

	open(200,file='pat.dat')
	
	
	
	
	call random_seed
	
	
	t0=0
	nt=1000                   
	c=10d-1
	
!====INITIAL CONDITION==================================
	do i=1,nx
    	do j=1,ny
        v(i,j)=2d-1					            
    	enddo
	enddo



	
!=========================================================================================================



	do j=2,ny-1
     	v(1,j)=v(nx-1,j)
	v(nx,j)=v(2,j)	
    	enddo


	do i=1,nx
        v(i,1)=v(i,ny-1)
	v(i,ny)=v(i,2)
    	enddo

	
	
	do i=1,nx				
	do j=1,ny       				
	call random_number (harvest=r1)
	call random_number (harvest=r2)	
	u1(i,j)=dsqrt(-2d0*D0*lambda*dlog(r1))*dcos(twopi*r2)							
	enddo
	enddo	

	do
	write(*,*)c,t0


	E=dexp(-lambda*dt)

	do kk=1,nt
	t=dt*float(kk)

        do i=2,nx-1
        do j=2,ny-1
                												
       !%%%%%%%%%%%%%%deterministic skeleton of model%%%%%%%%%%%%%%%%%%%%%
        g=v(i,j)
        rxv=a-(b*v(i,j))+ ((c*v(i,j)**2)/(1+(v(i,j)**2))) 	

        dfv=(D*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4.0*v(i,j)))/(dx**2)

        v2(i,j)= v(i,j)+dt*(rxv+dfv)+(dt)*u1(i,j)*v(i,j)



	call random_number (harvest=r5)
        call random_number (harvest=r6)
	h=dsqrt(-2d0*D0*lambda*(1d0-E**2)*dlog(r5))*dcos(twopi*r6)					
        u1(i,j)=u1(i,j)*E+h  
        
        
                   
        enddo
        enddo

!	 Boundary conditions:



	do j=2,ny-1
	v2(1,j)=v2(nx-1,j)
	v2(nx,j)=v2(2,j)	
	enddo


	do i=1,nx
	v2(i,1)=v2(i,ny-1)
	v2(i,ny)=v2(i,2)
	enddo


	  do i=1,nx
	  do j=1,ny
	  v(i,j)=v2(i,j)
	  v(i,j)=v2(i,j)
          enddo 
      	  enddo
   

   
	enddO   

	   
	   do i=1,nx	   
	   write(200,*)(v(i,j),j=1,ny)
	   enddo
	   	   
		
	  suma=0d0
          do i=1,nx
          do j=1,ny
          suma=suma+v(i,j)
	  enddo
	  enddo
	  
	  meand=suma/(nx*ny)		
	 
	  t0=t0+1
	  write(2000,*)c,t0,meand
		
	  c=c+1d-3
	   
	   if(c.gt.4d-1)exit
	   enddo






	stop
		
	end



