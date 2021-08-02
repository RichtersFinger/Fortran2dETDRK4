! Copyright (c) 2021 Steffen Richters-Finger
! ETDRK4: 
! Cox, Steven M., and Paul C. Matthews. "Exponential time differencing for stiff systems." Journal of Computational Physics 176.2, 430-455 (2002).
! Kassam, Aly-Khan, and Lloyd N. Trefethen. "Fourth-order time-stepping for stiff PDEs." SIAM Journal on Scientific Computing 26.4 1214-1233 (2005).
module ETDRK4
	use FFTW2d
	use, intrinsic :: iso_c_binding
	implicit none
	
	type(C_PTR) :: p
	complex(C_DOUBLE_COMPLEX), pointer :: ctemp(:, :), freq_partialx(:, :), freq_partialy(:, :)
	real(C_DOUBLE), pointer :: rtemp(:, :), real_partialx(:, :), real_partialy(:, :)

	integer (kind = 8) :: fftw_plan_freq_partialx, fftw_plan_freq_partialy, fftw_plan_rtemp

	integer :: nx, ny
	real( kind = 8 ) :: dtime, Lx, Ly

	! model-specific parameters
	real( kind = 8 ), allocatable :: controlparam(:)
	integer :: modelid = -1

	integer, allocatable :: karrayx(:)
	integer, allocatable :: karrayy(:)
	integer, parameter :: Ncontourintegral = 64								
	real (kind = 8), allocatable :: qarray(:, :), derivx(:, :), derivy(:, :)
	complex (kind = 8), allocatable :: a0(:, :), a1(:, :), a2(:, :), a3(:, :), a4(:, :)
	real(kind = 8) :: time
	real(kind = 8), allocatable :: real_data(:, :)	
	complex(kind = 8), allocatable :: freq_data(:, :)

contains
	! initialize the module with a time step, system size information, number of samples, model identifier, and an array of control parameters
	subroutine initETDRK4( dtime0, Lx0, Ly0, nx0, ny0, modelid0, controlparam0 )
		real (kind = 8) :: dtime0, Lx0, Ly0
		character (*) :: modelid0
		real (kind = 8) :: controlparam0(:)
		integer :: ix, iy, nx0, ny0
		real(C_DOUBLE), pointer :: arr(:,:)

		! save parameters
		dtime = dtime0
		Lx = Lx0	
		Ly = Ly0	
		nx = nx0	
		ny = ny0
		select case (modelid0)
		case ('KSE') ! damped Kuramoto-Sivashinsky equation (KSE), parameter controls damping
			modelid = 0
		case ('NE') ! damped Nikolaevskiy equation (NE), parameters control damping + drive
			modelid = 1
		case default
			write(*,*) 'model abbreviation', modelid0, ' not recognized'
			modelid = -1
		end select
		allocate(controlparam(size(controlparam0(:))))
		controlparam = controlparam0
		
		time = 0.0d0
		
		! allocate memory
		allocate(real_data(nx, ny))
		allocate(freq_data(nx/2 + 1, ny))
		allocate(a0(nx/2 + 1, ny))
		allocate(a1(nx/2 + 1, ny))
		allocate(a2(nx/2 + 1, ny))
		allocate(a3(nx/2 + 1, ny))
		allocate(a4(nx/2 + 1, ny))
		allocate(karrayx(nx/2 + 1))
		allocate(karrayy(ny))
		allocate(qarray(nx/2 + 1, ny))
		allocate(derivx(nx/2 + 1, ny))
		allocate(derivy(nx/2 + 1, ny))	

		! allocate aligned memory (increased performance of FFTW)
  		p = fftw_alloc_real(int(nx * ny, C_SIZE_T))
  		call c_f_pointer(p, rtemp, [nx, ny])
  		p = fftw_alloc_real(int(nx * ny, C_SIZE_T))
  		call c_f_pointer(p, real_partialx, [nx, ny])
  		p = fftw_alloc_real(int(nx * ny, C_SIZE_T))
  		call c_f_pointer(p, real_partialy, [nx, ny])
		p = fftw_alloc_complex(int((nx/2 + 1) * ny, C_SIZE_T))
		call c_f_pointer(p, ctemp, [nx/2 + 1, ny])
		p = fftw_alloc_complex(int((nx/2 + 1) * ny, C_SIZE_T))
		call c_f_pointer(p, freq_partialx, [nx/2 + 1, ny])
		p = fftw_alloc_complex(int((nx/2 + 1) * ny, C_SIZE_T))
		call c_f_pointer(p, freq_partialy, [nx/2 + 1, ny])

		! make FFTW plans
		! note: fftw-c2r plans destroy input array data
		call dfftw_plan_dft_c2r_2d(fftw_plan_freq_partialx, &
					nx, ny, freq_partialx, real_partialx, FFTW_ESTIMATE )
		call dfftw_plan_dft_c2r_2d(fftw_plan_freq_partialy, &
					nx, ny, freq_partialy, real_partialy, FFTW_ESTIMATE )
		call dfftw_plan_dft_r2c_2d(fftw_plan_rtemp, nx, ny, rtemp, ctemp, FFTW_ESTIMATE )

		! prepare supporting arrays for integration
		do iy = 1, ny/2 + 1
 			karrayy(iy) = iy - 1
		end do
		do iy = ny/2 + 2, ny
			karrayy(iy) = -ny + iy - 1 
		end do
	
		do ix = 1, nx/2 + 1
			karrayx(ix) = ix - 1
		end do

		do iy = 1, ny
			derivx(:, iy) = 2.0d0 * PI / Lx * karrayx
		end do	
		do ix = 1, nx/2 + 1
			derivy(ix, :) = 2.0d0 * PI / Ly * karrayy
		end do	
		
		do iy = 1, ny
			do ix = 1, nx/2 + 1
				! qarray contains the linear part as in regular exponential time differencing
				select case (modelid)
				case (0) ! KSE
					qarray(ix, iy) = - controlparam(1) + ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
										+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)&
							- ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
								+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)**2.0d0
				case (1) ! NE
					qarray(ix, iy) = - controlparam(1) + (controlparam(2) - 1.0d0) * ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
										+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)&
							+ 2.0d0 * ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
								+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)**2.0d0&
							- ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
							+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)**3.0d0
				case default
					qarray(ix, iy) = ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
										+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)&
							- ((2.0d0 * PI / Lx * karrayx(ix))**2.0d0 &
								+ (2.0d0 * PI / Ly * karrayy(iy))**2.0d0)**2.0d0
				end select

				a0(ix, iy) = exp( dtime * qarray(ix, iy) / 2.0d0 )
				if (abs(dtime * qarray(ix, iy)) < 0.5d0) then
					! see: Kassam, Aly-Khan, and Lloyd N. Trefethen. "Fourth-order time-stepping for stiff PDEs." SIAM Journal on Scientific Computing 26.4 1214-1233 (2005).
					a1(ix, iy) = dtime * &!
							contourInta1( dtime * qarray(ix, iy) )
					a2(ix, iy) = dtime * &!
							contourInta2( dtime * qarray(ix, iy) )
					a3(ix, iy) = dtime * &!
							contourInta3( dtime * qarray(ix, iy) )
					a4(ix, iy) = dtime * &!
							contourInta4( dtime * qarray(ix, iy) )
				else
					a1(ix, iy) = (exp(dtime * qarray(ix, iy) / 2.0d0) - 1.0d0) / qarray(ix, iy)
					a2(ix, iy) = (-4.0d0 - dtime * qarray(ix, iy) + &
						exp(dtime * qarray(ix, iy)) * (4.0d0 - 3.0d0 * dtime * qarray(ix, iy) + (dtime * qarray(ix, iy))**2.0d0)) / &
						(qarray(ix, iy)**3.0d0 * dtime**2.0d0)
					a3(ix, iy) = 2.0d0 * (2.0d0 + dtime * qarray(ix, iy) + &
						exp(dtime * qarray(ix, iy)) * (-2.0d0 + dtime * qarray(ix, iy))) / &
						(qarray(ix, iy)**3.0d0 * dtime**2.0d0)
					a4(ix, iy) = (-4.0d0 - 3.0d0 * dtime * qarray(ix, iy) - (dtime * qarray(ix, iy))**2.0d0 + &
						exp(dtime * qarray(ix, iy)) * (4.0d0 - dtime * qarray(ix, iy))) / &
						(qarray(ix, iy)**3.0d0 * dtime**2.0d0)
				end if
			end do
		end do	

		! helper functions for initialization
		contains
		! evaluates contour-integral at z for coefficients a1
		real (kind = 8) function contourInta1( z )
			real (kind = 8) :: z, alpha
			complex (kind = 8) :: integral, t
		
			integer :: k

			! alpha = 0
			alpha = 0.0d0
			t = z + 1.0d0
			integral = (exp(t/2.0d0) - 1.0d0)/t &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			do k = 1, Ncontourintegral - 1
				alpha = k * PI / Ncontourintegral
				t = z + exp(alpha * (0.0d0, 1.0d0))
				integral = integral + & 
						2.0d0 * (exp(t/2.0d0) - 1.0d0)/t &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)
			end do
			! alpha = PI
			alpha = PI
			t = z - 1.0d0
			integral = integral + (exp(t/2.0d0) - 1.0d0)/t &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			contourInta1 = real(integral/ 2.0d0 / Ncontourintegral) 

		end function
		! evaluates contour-integral at z for coefficients a2
		real (kind = 8) function contourInta2( z )
			real (kind = 8) :: z, alpha
			complex (kind = 8) :: integral, t
		
			integer :: k

			! alpha = 0
			alpha = 0.0d0
			t = z + 1.0d0
			integral = (-4.0d0 - t + exp(t) * (4.0d0 - 3.0d0 * t + t * t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			do k = 1, Ncontourintegral - 1
				alpha = k * PI / Ncontourintegral
				t = z + exp(alpha * (0.0d0, 1.0d0))
				integral = integral + & 
						2.0d0 * (-4.0d0 - t + exp(t) * (4.0d0 - 3.0d0 * t + t * t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)
			end do
			! alpha = PI
			alpha = PI
			t = z - 1.0d0
			integral = integral + (-4.0d0 - t + exp(t) * (4.0d0 - 3.0d0 * t + t * t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			contourInta2 = real(integral/ 2.0d0 / Ncontourintegral) 

		end function
		! evaluates contour-integral at z for coefficients a3
		real (kind = 8) function contourInta3( z )
			real (kind = 8) :: z, alpha
			complex (kind = 8) :: integral, t
		
			integer :: k

			! alpha = 0
			alpha = 0.0d0
			t = z + 1.0d0
			integral = 2.0d0 * (2.0d0 + t + exp(t) * (-2.0d0 + t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			do k = 1, Ncontourintegral - 1
				alpha = k * PI / Ncontourintegral
				t = z + exp(alpha * (0.0d0, 1.0d0))
				integral = integral + & 
						2.0d0 * 2.0d0 * (2.0d0 + t + exp(t) * (-2.0d0 + t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)
			end do
			! alpha = PI
			alpha = PI
			t = z - 1.0d0
			integral = integral + 2.0d0 * (2.0d0 + t + exp(t) * (-2.0d0 + t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			contourInta3 = real(integral/ 2.0d0 / Ncontourintegral) 

		end function
		! evaluates contour-integral at z for coefficients a4
		real (kind = 8) function contourInta4( z )
			real (kind = 8) :: z, alpha
			complex (kind = 8) :: integral, t
		
			integer :: k

			! alpha = 0
			alpha = 0.0d0
			t = z + 1.0d0
			integral = (-4.0d0 - 3.0d0 * t - t * t + exp(t) * (4.0d0 - t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			do k = 1, Ncontourintegral - 1
				alpha = k * PI / Ncontourintegral
				t = z + exp(alpha * (0.0d0, 1.0d0))
				integral = integral + & 
						2.0d0 * (-4.0d0 - 3.0d0 * t - t * t + exp(t) * (4.0d0 - t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)
			end do
			! alpha = PI
			alpha = PI
			t = z - 1.0d0
			integral = integral + (-4.0d0 - 3.0d0 * t - t * t + exp(t) * (4.0d0 - t))/t**3.0d0 &
						* exp(alpha * (0.0d0, 1.0d0))/ (t - z)

			contourInta4 = real(integral/ 2.0d0 / Ncontourintegral) 

		end function
	end subroutine
	subroutine ETDRK4_free( )
		call fftw_free(p)
		call dfftw_destroy_plan( fftw_plan_freq_partialx )
		call dfftw_destroy_plan( fftw_plan_freq_partialy )
		call dfftw_destroy_plan( fftw_plan_rtemp )

		deallocate(real_data)
		deallocate(freq_data)
		deallocate(a0)
		deallocate(a1)
		deallocate(a2)
		deallocate(a3)
		deallocate(a4)
		deallocate(karrayx)
		deallocate(karrayy)
		deallocate(qarray)
		deallocate(derivx)
		deallocate(derivy)

		deallocate(controlparam)
	end subroutine
	! remove constant offset in current data arrays
	subroutine resetzeromode()
		freq_partialx = freq_data
		freq_partialx(1, 1) = 0.0d0
		call dfftw_execute_dft_c2r( fftw_plan_freq_partialx, freq_partialx, real_partialx)
		real_data = real_partialx
		freq_data(1, 1) = 0.0d0
	end subroutine
	! calculate the nonlinearity for the Fourier data "somefreq_data" at time "sometime"
	function getNonlinearity(somefreq_data, sometime)
		real (kind = 8) :: sometime
		complex (kind = 8), dimension(nx/2 + 1, ny) :: getNonlinearity
		complex (kind = 8), dimension(nx/2 + 1, ny), intent(in) :: somefreq_data

		select case (modelid)
		case (0) ! KSE
			freq_partialx = derivx * (0.0d0, 1.0d0) * somefreq_data
			freq_partialy = derivy * (0.0d0, 1.0d0) * somefreq_data
			call dfftw_execute_dft_c2r( fftw_plan_freq_partialx, freq_partialx, real_partialx)
			call dfftw_execute_dft_c2r( fftw_plan_freq_partialy, freq_partialy, real_partialy)
			rtemp = (real_partialx * real_partialx + real_partialy * real_partialy)			
			call dfftw_execute_dft_r2c( fftw_plan_rtemp, rtemp, ctemp)
			getNonlinearity = ctemp / real(nx * ny)
		case (1) ! NE
			freq_partialx = derivx * (0.0d0, 1.0d0) * somefreq_data
			freq_partialy = derivy * (0.0d0, 1.0d0) * somefreq_data
			call dfftw_execute_dft_c2r( fftw_plan_freq_partialx, freq_partialx, real_partialx)
			call dfftw_execute_dft_c2r( fftw_plan_freq_partialy, freq_partialy, real_partialy)
			rtemp = (real_partialx * real_partialx + real_partialy * real_partialy)			
			call dfftw_execute_dft_r2c( fftw_plan_rtemp, rtemp, ctemp)
			getNonlinearity = -0.5d0 * ctemp / real(nx * ny)
		case default
			getNonlinearity = 0.0d0
		end select

	end function	
	! nstep: execute nstep integration steps
	! breset: reset zero-mode after each integration step if desired like for the undamped KSE
	subroutine step( nstep, breset )
		real (kind = 8) :: time0
		complex (kind = 8), dimension(nx/2 + 1, ny) :: N1, N2, N3, N4
		complex (kind = 8), dimension(nx/2 + 1, ny) :: freqtemp, freqtemp2
		integer :: nstep, istep, jstep
		logical, optional :: breset

		time0 = time
		do istep = 1, nstep
			! calculate first interstep
			N1 = getNonlinearity(freq_data, time0 + real(istep-1)*dtime)
			freqtemp = a0 * freq_data + a1 * N1
			freqtemp2 = freqtemp ! backup for later interstep
		
			! calc second interstep
			N2 = getNonlinearity(freqtemp, time0 + (real(istep-1) + 0.5d0)*dtime)
			freqtemp = a0 * freq_data + a1 * N2
		
			! calc third interstep
			N3 = getNonlinearity(freqtemp, time0 + (real(istep-1) + 0.5d0)*dtime)
			freqtemp = a0 * freqtemp2 + a1 * (2.0d0 * N3 - N1)

			! do step
			N4 = getNonlinearity(freqtemp, time0 + (real(istep-1) + 1.0d0)*dtime)
			freq_data = a0 * a0 * freq_data + a2 * N1 + a3 * (N2 + N3) + a4 * N4
		
			! reset zero-mode if desired like for undamped KSE
			if (present(breset)) then
				if (breset) freq_data(1, 1) = 0.0d0
			end if

			! for real-valued field equation the corresponding complex-valued Fourier-field is symmetric
			! this is not reflected by the integration scheme and has to be restored
			! resymmetrize freq_data(1, :) and freq_data(nx/2 + 1, :)
			do jstep = 1, ny/2 - 1
				freq_data(1, ny + 1 - jstep) = dconjg(freq_data(1, jstep + 1))
				freq_data(nx/2 + 1, ny + 1 - jstep) = dconjg(freq_data(nx/2 + 1, jstep + 1))
			end do
			freq_data(1, 1) = real(freq_data(1, 1))
			freq_data(nx/2 + 1, 1) = real(freq_data(nx/2 + 1, 1))
			! resymmetrize using fft (easy and less error-prone but lazy and costly)
			!freq_partialx = freq_data
			!call dfftw_execute_dft_c2r( fftw_plan_freq_partialx, freq_partialx, real_partialx)
			!rtemp = real_partialx		
			!call dfftw_execute_dft_r2c( fftw_plan_rtemp, rtemp, ctemp)
			!freq_data = ctemp / real(nx * ny)
		end do
		! apply step to real field
		freq_partialx = freq_data
		call dfftw_execute_dft_c2r( fftw_plan_freq_partialx, freq_partialx, real_partialx)
		real_data = real_partialx
		time = time0 + real(nstep) * dtime

	end subroutine
	! use to read data from previous call to snapshot like "call snapshot( real_data, 0, 'backup.dat' )"
	! this requires the field to be written in a matrix format (like "snapshot" does)
	subroutine initfromfile( filename )
		integer :: i
		character(*) :: filename
		
		open (unit = 99, file = filename, status = 'old')

		do i = 1, ny
			read(99, *) real_data(:, i)
		end do

		freq_data = dfft(nx, ny, real_data)

		close(99)
	end subroutine
	! cinput is data in Fourier space; returns a complex array in fourier space with number of samples increased by factor
	function upscaledfield( cinput, factor )
		integer :: factor, ix, iy
		complex(kind = 8) :: upscaledfield(factor*nx/2 + 1, factor*ny), cinput(nx/2 + 1, ny)

		upscaledfield = 0.0d0
		do iy = 1, ny/2 + 1
			do ix = 1, nx/2 + 1
				upscaledfield(ix, iy) = cinput(ix, iy)
			end do
		end do
		do iy = ny, ny/2 + 2, -1
			do ix = 1, nx/2 + 1
				upscaledfield(ix, (factor - 1) * ny + iy) = cinput(ix, iy)
			end do
		end do
	end function
	! writes real data array rinput as matrix into file; can be used as initial state and loaded using initfromfile( filename )
	! ksnapshots: counter for standard output file format
	! filename0: optional parameter for output filepath
	subroutine snapshot( rinput, ksnapshots, filename0 )
		integer :: ksnapshots
		integer :: isnapshot, is
		character( len = 999 ) :: filename, dummy
		character( len = 8 ) :: out_format = '(I5.5)'
		real ( kind = 8 ) :: rinput( :, : )
		character( * ), optional :: filename0

		if (present(filename0)) then
			filename = filename0
		else
			! convert int to str
			write (dummy, out_format) ksnapshots
			! combine parts
			filename = 'data/out'//trim(dummy)//'.dat'
		end if
		
		open(unit = 110, file = trim(filename), status = 'replace')

		do is = 1, size(rinput(1, :))
			write(110, *) rinput(:, is)
		end do
		
		close(110)

	end subroutine
	! writes Fourier space data array cinput as matrix into file;
	! imode: identifier that determines the written data (real/imaginary part, absolute value, or phase)
	! ksnapshotsf: counter for standard output file format
	! filename0: optional parameter for output filepath
	subroutine snapshotf( cinput, imode, ksnapshotsf, filename0 )
		integer :: ksnapshotsf
		integer :: isnapshot, imode, is, js
		character( len = 999 ) :: filename, dummy
		character( len = 8 ) :: out_format = '(I5.5)'
		complex ( kind = 8 ) :: cinput( nx/2 + 1, ny )
		real ( kind = 8 ) :: routput( nx, ny )
		character( * ), optional :: filename0
		

		if (present(filename0)) then
			filename = filename0
		else
			! convert int to str
			write (dummy, out_format) ksnapshotsf
			! combine parts
			filename = 'data/out'//trim(dummy)//'.dat'
		end if

		open(unit = 110, file = trim(filename), status = 'replace')

		routput = 0.0d0

		select case( imode )
		case( 0 )
			do is = 0, nx/2 - 1
				do js = 0, ny/2
					routput(1 + is, ny/2 + js) = abs(cinput(nx/2 + 1 - is, js + 1))
					routput(nx/2 + 1 + is, ny/2 + js) = abs(cinput(is + 1, js + 1))
				end do
				do js = 1, ny/2 - 1
					routput(1 + is, js) = abs(cinput(nx/2 + 1 - is, ny/2 + js + 1))
					routput(nx/2 + 1 + is, js) = abs(cinput(is + 1, ny/2 + js + 1))
				end do
			end do
		case( 1 ) ! real
			do is = 0, nx/2 - 1
				do js = 0, ny/2
					routput(1 + is, ny/2 + js) = real(cinput(nx/2 + 1 - is, js + 1))
					routput(nx/2 + 1 + is, ny/2 + js) = real(cinput(is + 1, js + 1))
				end do
				do js = 1, ny/2 - 1
					routput(1 + is, js) = real(cinput(nx/2 + 1 - is, ny/2 + js + 1))
					routput(nx/2 + 1 + is, js) = real(cinput(is + 1, ny/2 + js + 1))
				end do
			end do
		case( 2 ) ! aimag
			do is = 0, nx/2 - 1
				do js = 0, ny/2
					routput(1 + is, ny/2 + js) = aimag(cinput(nx/2 + 1 - is, js + 1))
					routput(nx/2 + 1 + is, ny/2 + js) = aimag(cinput(is + 1, js + 1))
				end do
				do js = 1, ny/2 - 1
					routput(1 + is, js) = aimag(cinput(nx/2 + 1 - is, ny/2 + js + 1))
					routput(nx/2 + 1 + is, js) = aimag(cinput(is + 1, ny/2 + js + 1))
				end do
			end do
		case( 3 ) ! phase
			do is = 0, nx/2 - 1
				do js = 0, ny/2
					routput(1 + is, ny/2 + js) = atan2(aimag((cinput(nx/2 + 1 - is, js + 1))), &
										real((cinput(nx/2 + 1 - is, js + 1))))
					routput(nx/2 + 1 + is, ny/2 + js) = atan2(abs(cinput(is + 1, js + 1)), &
										real(cinput(is + 1, js + 1)))
				end do
				do js = 1, ny/2 - 1
					routput(1 + is, js) = atan2(aimag(cinput(nx/2 + 1 - is, ny/2 + js + 1)), &
										real(cinput(nx/2 + 1 - is, ny/2 + js + 1)))
					routput(nx/2 + 1 + is, js) = atan2(aimag(cinput(is + 1, ny/2 + js + 1)),&
										real(cinput(is + 1, ny/2 + js + 1)))
				end do
			end do
		end select
		do is = 1, ny 
			write(110, *) routput(:, is)
		end do
		close(110)
	end subroutine
	! adds a perturbation to the current data arrays
	! ramplitude: amplitude for perturbation
	! rseed: optional parameter to change seed to some specific value
	subroutine perturbation(ramplitude, rseed)
		real(kind = 8) :: ramplitude 
		integer :: irx, iry
		integer, optional :: rseed 

		if (present(rseed)) &
			call srand(rseed)

		do irx = 1, nx
			do iry = 1, ny
				real_data(irx, iry) = real_data(irx, iry) &
							+ 2.0d0 * ramplitude * (rand() - 0.5d0)
			end do
		end do
		freq_data = dfft(nx, ny, real_data)

	end subroutine
end module
