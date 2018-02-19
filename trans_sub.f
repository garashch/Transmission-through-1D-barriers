	program main
	implicit real*8 (a-h,o-z)
!
c------ SG 02/19/2018 the only purpose is to call 
c----- the subroutine computing transmission through 
c----- the Eckart, rectangular barrier and WKB for the latter
c----- needs potential parameters from 'pot.inc'
!
	write(*,*) 'enter energy as fraction of the barrier top'

	read(*,*) E

	call prob(E,Teck,Trec,Twkb)

	write(*,*) E, Teck, Trec, Twkb

	stop
	end

c-----------------------------------------------------------------

	subroutine prob(E,Teck,Trec,Twkb)
!-----------------------------------------------------------
!------ assume E is given as fraction of the barrier top V0
c-----------------------------------------------------------
c----------- S Garashchuk 02/19/2018 ------------------
c-------- Symmetric Eckart reactangular  BARRIER ------

c---- computes exact QM transmission probability ------
c---- through a symmetric rectangular barrier in 1D ---
c---- parameters are converted to atomic units  -------
c
c---- V0 is the barrier height ------------------------ 
c---- 1000 cm^{-1}/219474.62 = 0.004556335489 hartree -
c
c---- m is the particle mass  -------------------------
c---- am = 1836.15 a.m.u.; taken as the proton mass ---
c
c----  barrier width sqrt(md^2/m) is converted to bohr
c----  then we take d as half-width -------------------
c------------------------------------------------------
	implicit real*8 (a-h,o-z)

	real*8 md2

	include 'pot.inc'	! m, V0 and d are here 

	parameter(Eh_to_cm = 219474.62d0,bohr_to_ang = 0.529177d0) 

	parameter(au_mass = 9.10938356d-31, angstrom = 1d-10) 
c------------------------------------------------------

	pi = 4d0*atan(1d0)

	V0 = Vcm/Eh_to_cm
	
	E = E*V0+1d-16	! energy is given as fraction of V0

!	write(*,*) 'PARAMETERS ARE FROM pot.inc'
!	write(*,*)
!	write(*,*) 'MASS = ',am,' a.m.u'
!	write(*,*)
!	write(*,*) 'BARRIER TOP =',Vcm,' 1/cm'
!	write(*,*)
!	write(*,*) 'BARRIER m*d*d=',md2,' kg*m*m' 
!
	d =  dsqrt(md2/am/au_mass)/angstrom

!	write(*,*) 
!	write(*,*) 'FOR RECTANGULAR BARRIER WIDTH =', d,'Angstrom'
!	write(*,*) 

	d = d/bohr_to_ang

!	write(*,*) 'IN ATOMIC UNITS m,V0,d:'
!	write(*,*) 
!	write(*,*) am,v0,d

c-----	define the Eckart  parameter d;  V=V0/cosh(x/d)^2

	d = d/2d0 !to match the areas under the barrier int(V)=2*d*V0

!	write(*,*) 
!	write(*,*) 'FOR ECKART BARRIER V=V0/cosh(x/d)^2  d =', d,'bohr'
!	write(*,*) 
!	write(*,*) 'AREA under the barrier=',2d0*d*v0
!	write(*,*) 

	flox = dsqrt(abs(8d0*am*V0*d*d-1d0))*pi

        ak =  dsqrt(2*am*E)  	! wave vector before/after the barrier

	q = 2d0*pi*ak*d 	!Eckart 

	Teck = (cosh(q)-1d0)/(cosh(q)+cosh(flox)) !Eckart probability

	akap = dsqrt(2*am*abs(V0-E)) ! wave vector inside rect barrier

        x = 2d0*akap*d  	! argument of hyperbolic functions

        if(E.lt.V0) then       ! under the barrier

         eps = (akap/ak-ak/akap)/2d0

         Trec =  1d0/(cosh(x)**2+(eps*sinh(x))**2) ! exact probability

        else                   ! above the barrier

         eps = (akap/ak+ak/akap)/2d0

         Trec =  1d0/(cos(x)**2+(eps*sin(x))**2) ! exact probability

        endif

        Twkb = 1d0     ! WKB probability for rectangular barrier

        if(E.lt.V0) Twkb =  exp(-4d0*akap*d)

	write(*,*) E,Teck,Trec,Twkb

c------------------- plot the potential ------------------
!
!	open(10,file='pot_eckart')
!
!	np = 200
!
!	dx = 50d0/d/np
!
!	xmn = -np/2*dx
!
!	do i=0,np
!
!	 x =  xmn+i*dx
!
!	 pot = v0/cosh(x/d)**2 !THIS IS THE ECKART POTENTIAL
!	
!	 write(10,*) x,pot
!
!	enddo
!
!	close(10)
c------------------------------------------------------------
	return

	end

