	program main
c------------------------------------------------------
c----------- S Garashchuk 02/12/2018 ------------------
c-------- RECTANGULAR BARRIER -------------------------
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

	V0 = Vcm/Eh_to_cm

	write(*,*) 'PARAMETERS ARE FROM pot.inc'
	write(*,*)
	write(*,*) 'MASS = ',am,' a.m.u'
        write(*,*)
        write(*,*) 'BARRIER TOP =',Vcm,' 1/cm'
        write(*,*)
        write(*,*) 'BARRIER m*d*d=',md2,' kg*m*m'

	d =  dsqrt(md2/am/au_mass)/angstrom

	write(*,*) 
	write(*,*) 'RECTANGULAR BARRIER WIDTH =', d,'Angstrom'
	write(*,*) 

	d = d/bohr_to_ang

	write(*,*) 'IN ATOMIC UNITS m,V0,d:'
	write(*,*) 
	write(*,*) am,v0,d
	write(*,*) 
	write(*,*) 'Area under the barrier =', v0*d
	write(*,*)

	d = d/2d0 !half-width is used in the formulas

	write(*,*) 'emax/V0, number of points'

	read(*,*) emx,np

	de = emx*v0/np

	open(11,file='prob_rec') !open output file

	do i = 1,np

	 E = de*i+1d-16

	 akap = dsqrt(2*am*abs(V0-E)) ! wave vector in the barrier 

	 ak =  dsqrt(2*am*E)  ! wave vector before/after the barrier

	 x = 2d0*akap*d  ! argument of hyperbolic functions

	 if(E.lt.V0) then  	! under the barrier 

	  eps = (akap/ak-ak/akap)/2d0

	  Txqm =  1d0/(cosh(x)**2+(eps*sinh(x))**2) ! exact probability

	 else  			! above the barrier

	  eps = (akap/ak+ak/akap)/2d0

	  Txqm =  1d0/(cos(x)**2+(eps*sin(x))**2) ! exact probability

	 endif

	 Twkb = 1d0	! WKB probability

	 if(E.lt.V0) Twkb =  exp(-4d0*akap*d) 

	 write(11,*) E,txqm,twkb,twkb/txqm

	enddo

	close(11) 	!close output file

	stop

	end

