! simple iteration method for solving SLAE
subroutine simpleIteration(n, A, b, x, max_iter, eps)
	implicit none

	integer										:: i, j, k, n, max_iter
	double precision							:: eps, resid
	double precision, dimension(-1:n+2)			:: x, x_old, b, beta
	double precision, dimension(-1:n+2,-1:n+2)	:: A, alpha

	do i = -1, n + 2
		do j = -1, n + 2
			if (i /= j) then
				alpha(i,j) = -A(i,j) / (A(i,i) + 1e-1)
			else
				alpha(i,j) = 0.0
			end if
		end do
		beta(i) = b(i) / (A(i,i) + 1e-1)
	end do

	x = b
	x_old = x
	resid = 1
	do k = 1, max_iter

		x = beta + matmul(alpha, x_old)

		resid = maxval(x - x_old)
		if (mod(k, 100) == 0) write(*,*) "res = ", resid
		if (resid < eps) exit

		x_old = x

	end do

end subroutine
