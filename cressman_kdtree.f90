 subroutine cressman(x, y, ob, xg, yg, roi, anal, mask, nobs, nx, ny)

   use kdtree2_module
   use time_kdtree

   implicit none
 
   real(8), intent(in),  dimension(nobs) :: x
   real(8), intent(in),  dimension(nobs) :: y
   real(4), intent(in),  dimension(nobs) :: ob
   real(8), intent(in),  dimension(nx)   :: xg
   real(8), intent(in),  dimension(ny)   :: yg
   real(4), intent(out), dimension(nx,ny) :: anal
   integer(4), intent(out), dimension(nx,ny) :: mask
   real(4), intent(in) :: roi
   integer, intent(in) :: nobs, nx, ny

! this is how you declare a tree in your main program

   real(kdkind), dimension(:,:), allocatable :: xy_obs
   real(kdkind), allocatable :: query(:)

   type(kdtree2), pointer :: tree

   type(kdtree2_result), allocatable :: results(:)
   integer, parameter :: nalloc = 10000
   integer            :: nfound

   integer n, i, j
   real(8) R2, w_sum, top, wk, rk2, t0, t1

   logical, parameter :: debug = .false.
 
   IF( debug ) THEN
      print *, 'Nobs:                   ', nobs
      print *, 'Nx/Ny:                  ', nx, ny
      print *, 'Maxval of anal before:  ', maxval(anal)
      print *, 'Minval of anal before:  ', minval(anal)
      print *, ''
      print *, 'Maxval of observations:  ', maxval(ob)
      print *, 'Minval of observations:  ', minval(ob)
      print *, 'Min/Max xob:    ', minval(x), maxval(x)
      print *, 'Min/Max yob:    ', minval(y), maxval(y)
      print *, 'Min/Max xgrid:  ', minval(xg), maxval(xg)
      print *, 'Min/Max ygrid:  ', minval(yg), maxval(yg)
      print *, 'Radius of influence:  ', roi
      print *, ''
   ENDIF

   R2 = roi**2

! Allocate coordinates

   allocate(xy_obs(2,nobs))
   allocate(results(nobs))
   allocate(query(2))

   results(:)%idx = -666

   DO n = 1,nobs
    xy_obs(1,n) = x(n)
    xy_obs(2,n) = y(n)
   ENDDO

!  call cpu_time(t0)

   tree => kdtree2_create(xy_obs,sort=.true.,rearrange=.true.)  ! this is how you create a tree.

!  call cpu_time(t1)
!  write (*,*) real(n)/real(t1-t0), ' points per second built for sort/rearranged tree.'
 
   DO j = 1,ny
    DO i = 1,nx

      query(1) = xg(i)
      query(2) = yg(j)

      call kdtree2_r_nearest(tp=tree, qv=query, r2=R2, nfound=nfound, nalloc=nalloc, results=results)

      w_sum = 0.0
      top   = 0.0  

      mask(i,j) = 0

      DO n = 1,nfound
        rk2 = results(n)%dis
        wk = max((R2-rk2) / (R2+rk2), 0.0)
        top = top + wk*ob(results(n)%idx)
        w_sum = w_sum + wk
      ENDDO
 
      IF (w_sum .ge. 0.01) THEN
        anal(i,j) = sngl(top/w_sum)
        mask(i,j) = 1
      ENDIF
 
    ENDDO
   ENDDO
 
   IF( debug ) THEN
      print *, 'Maxval of anal after:  ', maxval(anal)
      print *, 'Minval of anal after:  ', minval(anal)
   ENDIF

 end subroutine cressman
 
!======================================================================================================'
!
! Routine to cressman a set of obs to a 2D grid (could be a random grid).  Implemention
! uses the cKDTree algorithm in python to find the indices for each ob that is on the grid.
!
! Example code in python to create fields that are needed.
!=======================================================================================================
!# Begin python code
!
!import numpy, scipy.spatial
!import time
!
!# set up test grid
!
!x1d = numpy.arange(500) / 500.
!y1d = numpy.arange(600) / 600.
!z1d = numpy.arange(60)  / 60.
!
!y_array, z_array, x_array = numpy.meshgrid(y1d, z1d, x1d)
!combined_xyz_arrays = numpy.dstack([z_array.ravel(),y_array.ravel(),x_array.ravel()])[0]
!print 'XYZ ARRAY SHAPE', combined_xyz_arrays.shape
!
!obs = numpy.random.random(30).reshape(3,10)
!print 'Obs Shape:', obs.shape
!obs_list = list(obs.transpose())
!
!# Okay create the KDTree data structure
!
!mytree = scipy.spatial.cKDTree(combined_xyz_arrays)
!
!distances, indices1D = mytree.query(obs_list)
!
!indices3D = numpy.unravel_index(numpy.ravel(indices1D, y_array.size), y_array.shape)
!
!# these are the integer indices that you now pass into the fortran routine. They
!# are the un-raveled 3D index locations nearest the observation point in the 3D array
!
!k = indices3D[0]
!j = indices3D[1]
!i = indices3D[2]
!
!for n in numpy.arange(obs.shape[1]) :
!    print n, obs[0,n], obs[1,n], obs[2,n], z_array[k[n],j[n],i[n]], y_array[k[n],j[n],i[n]], x_array[k[n],j[n],i[n]]
! 
!# END python code
!======================================================================================================'
SUBROUTINE OBS_2_GRID2D(field, obs, xob, yob, xc, yc, ii, jj, hscale, missing, nobs, nx, ny)
                   
  implicit none

! Passed in variables

  real(kind=8), INTENT(OUT) :: field(ny,nx)      ! 2D analysis passed back to calling routine
      
  integer,      INTENT(IN)  :: nx, ny, nobs      ! grid dimensions
  real(kind=4), INTENT(IN)  :: xob(nobs)         ! x coords for each ob
  real(kind=4), INTENT(IN)  :: yob(nobs)         ! y coords for each ob
  real(kind=4), INTENT(IN)  :: obs(nobs)         ! obs 
  real(kind=4), INTENT(IN)  :: ii(nobs)          ! nearest index to the x-point on grid for ob
  real(kind=4), INTENT(IN)  :: jj(nobs)          ! nearest index to the x-point on grid for o 
  real(kind=4), INTENT(IN)  :: xc(ny,nx)         ! coordinates corresponding to WRF model grid locations
  real(kind=4), INTENT(IN)  :: yc(ny,nx)         ! coordinates corresponding to WRF model grid locations
  real(kind=4), INTENT(IN)  :: hscale
  real(kind=4), INTENT(IN)  :: missing
  
! Local variables

  integer i, j, i0, j0, n, i0m, i0p, j0m, j0p, idx, jdx       ! loop variables
  
  real(kind=8) dh, wgt
  real(kind=4), allocatable, dimension(:,:) :: sum, wgt_sum

! DEBUG
  logical, parameter :: debug = .false.
  
! Allocate local memory

  allocate(wgt_sum(ny, nx))
  allocate(sum(ny, nx))

! Initialize values

  field(:,:)   = missing  
  wgt_sum(:,:) = 0.0
  sum(:,:)     = 0.0

  IF( debug) THEN
    print *, "FORTRAN OBS_2_GRID2D:  dims  ", nx, ny, nobs
    print *, "FORTRAN OBS_2_GRID2D:  dx, hscale", xc(1,2) - xc(1,1), hscale
    print *, "FORTRAN OBS_2_GRID2D:  dy, hscale", yc(2,1) - yc(1,1), hscale
    print *, "FORTRAN OBS_2_GRID2D:    obs ", minval(obs), maxval(obs)
    print *, "FORTRAN OBS_2_GRID2D:  x-obs ", minval(xob), maxval(xob)
    print *, "FORTRAN OBS_2_GRID2D:  y-obs ", minval(yob), maxval(yob)
    print *, "FORTRAN OBS_2_GRID2D:  x-grid", minval(xc), maxval(xc)
    print *, "FORTRAN OBS_2_GRID2D:  y-grid", minval(yc), maxval(yc)
  ENDIF
     
  DO n = 1,nobs
  
    i0 = ii(n)
    j0 = jj(n)
    
    if(i0 .eq. 1) cycle
    if(j0 .eq. 1) cycle

    if(i0 .eq. nx) cycle
    if(j0 .eq. ny) cycle    

    idx = 1+nint(hscale/abs(xc(j0,i0) - xc(j0,i0-1)))
    jdx = 1+nint(hscale/abs(yc(j0,i0) - yc(j0-1,i0)))
    
    i0m = max(i0-idx,1)
    j0m = max(j0-jdx,1)
    
    i0p = min(i0+idx,nx)
    j0p = min(j0+jdx,ny)

    DO i = i0m, i0p
      DO j = j0m, j0p
        
          dh = ( (xc(j,i)-xob(n))**2 + (yc(j,i)-yob(n))**2) / hscale**2
          wgt = (1.0 - dh) / (1.0 + dh)
          
          IF (wgt .gt. 0.0) THEN
            sum(j,i) = sum(j,i) + wgt*obs(n)
            wgt_sum(j,i) = wgt_sum(j,i) + wgt
          ENDIF  
                    
      ENDDO    ! END J
    ENDDO     ! END I          
  ENDDO      ! END N
 
  WHERE( wgt_sum > 0.01 ) field = sum / wgt_sum
  
  deallocate(wgt_sum)
  deallocate(sum)

RETURN
END SUBROUTINE OBS_2_GRID2D
