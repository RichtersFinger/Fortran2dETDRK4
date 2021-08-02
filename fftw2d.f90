
module FFTW2d
	use, intrinsic :: iso_c_binding
	implicit none
   	include 'fftw3.f03'

	real (kind = 8) :: PI = 4.0d0 * atan(1.0d0)
contains
	function dfft( n1, n2, rinput )
		integer (kind = 8) :: p_r2c
		integer :: n1, n2 ! samples
		real ( kind = 8 ) :: rinput(n1, n2)
		complex ( kind = 8 ) :: dfft(n1/2 + 1, n2), coutput(n1/2 + 1, n2)
		
		call dfftw_plan_dft_r2c_2d(p_r2c, n1, n2, rinput, coutput, FFTW_ESTIMATE )
		call dfftw_execute_dft_r2c( p_r2c, rinput / real(n1*n2), coutput)
		call dfftw_destroy_plan( p_r2c )
		
		dfft = coutput
	end function
	function idfft( n1, n2, cinput )
		integer (kind = 8) :: p_c2r
		integer :: n1, n2 ! samples
		real ( kind = 8 ) :: idfft(n1, n2), routput(n1, n2)
		complex ( kind = 8 ) :: cinput(n1/2 + 1, n2), temp(n1/2 + 1, n2)
		
		call dfftw_plan_dft_c2r_2d(p_c2r, n1, n2, temp, routput, FFTW_ESTIMATE )
		temp = cinput
		call dfftw_execute_dft_c2r( p_c2r, temp, routput)
		call dfftw_destroy_plan( p_c2r )
		
		idfft = routput
	end function
end module
