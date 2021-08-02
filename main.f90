! compile using : gfortran fftw2d.f90 ETDRK4.f90 main.f90 -lfftw3 -o run -O3
program run
	use FFTW2d
	use ETDRK4
	implicit none

	integer :: i, j, k, simmode
	real (kind = 8) :: m

	m = 20.0d0
	! example for Kuramoto-Sivashinsky equation (KSE) 
	call initETDRK4( 0.1d0, m * 2.0d0 * PI, m * 2.0d0 * PI, 256, 256, 'KSE', [0.0d0] )
	! example for Nikolaevskiy equation (NE)
	!call initETDRK4( 0.02d0, m * 2.0d0 * PI, m * 2.0d0 * PI, 256, 256, 'NE', [0.0d0, 0.1d0] )

	real_data = 0.0d0
	freq_data = 0.0d0

	call perturbation(0.1d0, 2411)
	do i = 1, 600
		call step( 10, .true. )

		call snapshot(real_data, 0)
		call snapshot(real_data, i)
		!call snapshot(idfft(2*nx, 2*ny, upscaledfield( freq_data, 2 )), i)

		write(*,*) i, time
	end do

	call ETDRK4_free()
end program
