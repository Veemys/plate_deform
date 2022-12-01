program solveEquation
	implicit none

	integer, parameter								:: io = 69
	integer											:: n, max_iter, i, j
	double precision								:: E, Inert, thick, L, mu, H1, H2, Q1, Q2
	double precision								:: u, h, eps
	double precision, dimension(:), allocatable		:: x, b, w, g1, g2, g3
	double precision, dimension(:,:), allocatable	:: A

	! open(io, file = 'input.txt')
	! read(io,*) n						! number of nodes
	! read(io,*) E						! Young's modulus
	! read(io,*) delta					! 
	! read(io,*) mu						! viscosity
	! read(io,*) Q1						! 
	! read(io,*) Q2						! 
	! read(io,*) H1						! 
	! read(io,*) H2						! 
	! read(io,*) L						! length
	! close(io)
	
	E = 2e11
	thick = 0.01
	Inert = thick**(3.0 / 12.0)
	L = 0.8 			! Lenght of plate
	mu = 0.001 			! dynamic viscosity
	H1 = 0.04 			! height of upper channel
	H2 = 0.02 			! height of lower channel
	u = 0.008 			! mean velocity of water
	Q1 = u * H1 		! massflow in upper channer
	Q2 = u * H2 		! massflow in lower channel
	
	n = 10
	
	allocate(b(-1:n+2), x(-1:n+2), g1(-1:n+2), g2(-1:n+2), g3(-1:n+2), w(-1:n+2))
	allocate(A(-1:n+2,-1:n+2))
	
	write(*,*) 'E =', E
	write(*,*) 'thick =', thick
	write(*,*) 'Inert =', Inert
	write(*,*) 'L =', L
	write(*,*) 'mu =', mu
	write(*,*) 'H1 =', H1
	write(*,*) 'H2 =', H2
	write(*,*) 'u =', u
	write(*,*) 'Q1 =', Q1
	write(*,*) 'Q2 =', Q2
	write(*,*)
	
	max_iter = 100000
	eps = 1e-12
	
	call createGrid(n, L, x, h)
	call calcGCoeffs(x, mu, Q1, Q2, H1, H2, L, g1, g2, g3, n)
	
	do i = -1, n + 2
		write(*,*) g1(i), g2(i), g3(i)
	end do
	write(*,*)
	
	call formMatrixA(E, Inert, g1, g2, g3, h, n, A)
	! A = transpose(A)
	
	do i = -1, n + 2
		do j = -1, n + 2
			write(*,*) A(i,j)
		end do
		write(*,*)
	end do
	
	call formVectorB(n, b, g3)
	
	do i = -1, n + 2
		write(*,*) b(i)
	end do
	
	call simpleIteration(n, A, b, w, max_iter, eps)
	call writeToFile(n, x, w)

end

! create grid
subroutine createGrid(n, L, x, h)
	implicit none

	integer 								:: n, i
	double precision						:: L, h
	double precision, dimension(-1:n+2)		:: x

	h = L / (n - 1)

	x(-1) = - 2 * h
	do i = 0, n + 2

		x(i) = x(i-1) + h

	end do

end subroutine

! form matrix
subroutine formMatrixA(E, Inert, g1, g2, g3, h, n, A)
	implicit none

	integer										:: i, n
	double precision							:: E, Inert
	double precision							:: h
	double precision, dimension(-1:n+2)			:: g1, g2, g3, b
	double precision, dimension(-1:n+2,-1:n+2)	:: A

	A = 0.0

	! boundary condition	
	A(-1,-1) = 0.0
	A(-1,0) = 0.0
	A(-1,1) = 1.0
	A(-1,2) = 0.0
	A(-1,3) = 0.0
	b(-1) = 0.0
	
	A(0,-1) = 0.0
	A(0,0) = -1.0
	A(0,1) = 0.0
	A(0,2) = 1.0
	A(0,3) = 0.0
	b(0) = 0.0
	
	A(n+1,n-2) = 0.0
	A(n+1,n-1) = 1.0
	A(n+1,n) = -2.0
	A(n+1,n+1) = 1.0
	A(n+1,n+2) = 0.0
	b(n+1) = 0.0
	
	A(n+2,n-2) = -1.0
	A(n+2,n-1) = 2.0
	A(n+2,n) = 0.0
	A(n+2,n+1) = -2.0
	A(n+2,n+2) = 1.0
	b(n+2) = 0.0

	do i = 1, n

		A(i,i-2) = E * Inert / h**4
		A(i,i-1) = -4.0 * E * Inert / h**4 - g1(i) / h**2
		A(i,i) = 6.0 * E * Inert / h**4 + 2.0 * g1(i) / h**2 + g2(i)
		A(i,i+1) = -4.0 * E * Inert / h**4 - g1(i) / h**2
		A(i,i+2) = E * Inert / h**4

	end do

end subroutine

! form b vector
subroutine formVectorB(n, b, g3)
	implicit none

	integer									:: n
	double precision, dimension(-1:n+2)		:: b, g1, g2, g3

	b = g3

end subroutine

! calculation g's coeffs
subroutine calcGCoeffs(x, mu, Q1, Q2, H1, H2, L, g1, g2, g3, n)
	implicit none

	integer										:: i, n
	double precision							:: mu, Q1, Q2, H1, H2, L
	double precision, dimension(-1:n+2)			:: x, g1, g2, g3

	g1 = 0.0
	g2 = 0.0
	g3 = 0.0

	! do i = -1, n+2
	do i = 1, n
	
		g1(i) = 6.0 * mu * (Q1 / H1**2 + Q2 / H2**2) * (L - x(i))
		g2(i) = 36.0 * mu * (Q1 / H1**4 + Q2 / H2**4) * (L - x(i))
		g3(i) = -12.0 * mu * (Q1 / H1**3 - Q2 / H2**3) * (L - x(i))
		
	end do

end subroutine

! write
subroutine writeToFile(n, x, w)
	implicit none

	integer, parameter						:: io = 12
	integer									:: i, n
	double precision, dimension(n+2)		:: x, w

	open(io, file = "output.txt")

	do i = 1, n

		write(io,*) x(i), w(i)

	end do

end subroutine