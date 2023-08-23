program cam_test

  ! Read column(s) extracted from CAM's offline radiation driving dataset and produce fluxes.
  ! SW clear sky only.  MPI version.

  use mo_rte_kind,           only: wp, wl
  use mo_optical_props,      only: ty_optical_props, &
                                   ty_optical_props_arry, ty_optical_props_2str
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes,             only: ty_fluxes_broadband
  use mo_rte_sw,             only: rte_sw
  use mo_load_coefficients,  only: load_and_init
  use mo_rte_config,         only: rte_config_checks

  use netcdf,                only: nf90_open, NF90_WRITE, NF90_NOWRITE, NF90_NOERR, nf90_close
  use mo_simple_netcdf,      only: read_field, read_string, var_exists, get_dim_size, &
                                   write_field, create_dim, create_var
  use mpi

  implicit none

  ! Local Variables

  integer :: ierr, myid, numprocs

  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev
  real(wp), dimension(:,:),   allocatable :: p_lay_loc, t_lay_loc, p_lev_loc
  real(wp), dimension(:,:,:), allocatable :: gas_vmr
  real(wp), dimension(:,:,:), allocatable :: gas_vmr_loc

  ! Shortwave only
  real(wp), dimension(:),     allocatable :: mu0 ! cosine of solar zenith angle
  real(wp), dimension(:),     allocatable :: mu0_loc
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  real(wp), dimension(:),     allocatable :: asdir, asdif
  real(wp), dimension(:),     allocatable :: asdir_loc, asdif_loc

  ! Shortwave TOA flux
  real(wp) :: tsi
  real(wp), dimension(:,:), allocatable :: toa_flux

  ! Output variables
  real(wp), dimension(:,:), target, allocatable :: flux_up, flux_dn, flux_dir
  real(wp), dimension(:,:), target, allocatable :: flux_up_loc, flux_dn_loc, flux_dir_loc
  real(wp), dimension(:,:), target, allocatable :: sendbuf, rcvbuf

  ! Derived types from the RTE and RRTMGP libraries
  type(ty_gas_optics_rrtmgp)  :: k_dist
  type(ty_gas_concs)          :: gas_concs
  type(ty_optical_props_2str) :: atmos
  type(ty_fluxes_broadband)   :: fluxes


  integer :: nUserArgs = 0
  integer :: ncol, nlcol, nlay, nlev, nbnd, ngpt
  integer :: mxcols, mncols, ncolmx, nbinsmx
  integer :: i, igas, lcol, col_beg, col_end
  logical :: top_at_1

  integer, parameter :: ngas = 6           ! no data for 'co ' and 'n2 '
  character(len=3), dimension(ngas) :: &
     gas_names = ['h2o',       'co2',       'o3 ',       'n2o',       'ch4',       'o2 ']
  character(len=9), dimension(ngas) :: &
     cam_names = ['rad_Q    ', 'rad_CO2  ', 'rad_ozone', 'rad_N2O  ', 'rad_CH4  ', 'rad_O2   ']

  character(len=256) :: input_file, k_dist_file

  !======================================================================================
  ! Parse command line for file names
  nUserArgs = command_argument_count()
  if (nUserArgs /=  2) call stop_on_err("Need to supply input_file k_distribution_file.")
  call get_command_argument(1, input_file)
  call get_command_argument(2, k_dist_file)

  ! Read temperature, pressure, gas concentrations.
  ! All the data is read into each task before subsets of the columns are created
  ! for the SW flux calcs.
  call read_atmos(input_file, ncol, tsi, &
                  p_lay, t_lay, p_lev,   &
                  gas_vmr, mu0,          &
                  asdir, asdif)

  nlay = size(p_lay, 2)
  nlev = size(p_lev,2)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

  ! Divide the columns among bins (1 bin per task).  Max number of cols in each bin.
  mxcols = ceiling(real(ncol)/numprocs)

  ! Number of full bins
  nbinsmx = ncol/mxcols

  ! Number of columns in last bin.
  mncols = ncol - nbinsmx*mxcols

  ! Number of columns for this task
  nlcol = 0
  if (myid <= nbinsmx-1) then
     nlcol = mxcols
  else if (myid == nbinsmx) then
     nlcol = mncols
  end if

  ! Global begin and end indices for cols in local bin
  col_beg = myid*mxcols + 1
  col_end = min(col_beg+mxcols-1, ncol)

  ! If this task has no work then just finalize and stop this task, but
  ! allow the other tasks to complete...
  if (nlcol == 0) then
     call MPI_FINALIZE(ierr)
     stop
  end if

  ! allocate arrays for local task and distribute the data to them
  allocate( &
     p_lay_loc(nlcol,nlay), t_lay_loc(nlcol,nlay), p_lev_loc(nlcol,nlev), &
     gas_vmr_loc(nlcol,nlay,ngas), mu0_loc(nlcol), asdir_loc(nlcol), asdif_loc(nlcol))

  lcol = 0
  do i = col_beg, col_end
     lcol = lcol + 1
     p_lay_loc(lcol,:)     = p_lay(i,:)
     t_lay_loc(lcol,:)     = t_lay(i,:)
     p_lev_loc(lcol,:)     = p_lev(i,:)
     gas_vmr_loc(lcol,:,:) = gas_vmr(i,:,:)
     mu0_loc(lcol)         = mu0(i)
     asdir_loc(lcol)       = asdir(i)
     asdif_loc(lcol)       = asdif(i)
  end do

  ! Arrays containing all columns no longer needed.
  deallocate(p_lay, t_lay, p_lev, gas_vmr, mu0, asdir, asdif)

  ! Init object with gas vmr values for local columns.  The values from the
  ! local array are copied inside the gas_concs object.
  call stop_on_err(gas_concs%init(gas_names))
  do igas = 1, ngas
     call stop_on_err( &
        gas_concs%set_vmr(trim(gas_names(igas)), gas_vmr_loc(:,:,igas)))
  end do
  deallocate(gas_vmr_loc)

  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)

  call stop_on_err( &
     k_dist%set_tsi(tsi))! scales the TSI but does not change spectral distribution

  ! Problem sizes
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay_loc(1, 1) < p_lay_loc(1, nlay)

  ! Allocate arrays for the optical properties themselves.
  call stop_on_err( &
     atmos%alloc_2str(nlcol, nlay, k_dist))

  ! Boundary conditions
  allocate(toa_flux(nlcol, ngpt))
  allocate(sfc_alb_dir(nbnd, nlcol), sfc_alb_dif(nbnd, nlcol))

  ! CAM's albedos are not spectrally resolved, so set all bands to
  ! the same value.
  do i = 1, nlcol
     sfc_alb_dir(:,i) = asdir_loc(i)
     sfc_alb_dif(:,i) = asdif_loc(i)
  end do

  ! Fluxes
  allocate(flux_up_loc(nlcol,nlev), flux_dn_loc(nlcol,nlev))
  allocate(flux_dir_loc(nlcol,nlev))

  call rte_config_checks(logical(.false., wl))

  ! SW Solver
  fluxes%flux_up => flux_up_loc(:,:)
  fluxes%flux_dn => flux_dn_loc(:,:)
  fluxes%flux_dn_dir => flux_dir_loc(:,:)

  call stop_on_err( &
     k_dist%gas_optics(p_lay_loc, p_lev_loc, & ! in
                       t_lay_loc, gas_concs, & ! in
                       atmos, toa_flux))       ! out

  call stop_on_err( &
     rte_sw(atmos, top_at_1, &
            mu0_loc,   toa_flux, &
            sfc_alb_dir, sfc_alb_dif, &
            fluxes))
  
  ! Gather fluxes from all tasks.  Allocation for all columns only needs to be made
  ! on root task.
  ncolmx = mxcols * numprocs
  if (myid == 0) then
     allocate(flux_up(ncolmx,nlev), flux_dn(ncolmx,nlev), flux_dir(ncolmx,nlev))
     allocate(rcvbuf(nlev,ncolmx))
  end if

  ! Transpose the send arrays so the columns are sent in contiguous storage,
  allocate(sendbuf(nlev,nlcol))
  sendbuf = transpose(flux_up_loc)
  call MPI_GATHER( &
     sendbuf, nlcol*nlev, MPI_DOUBLE_PRECISION, rcvbuf, mxcols*nlev, &
     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (myid == 0) flux_up = transpose(rcvbuf)

  sendbuf = transpose(flux_dn_loc)
  call MPI_GATHER( &
     sendbuf, nlcol*nlev, MPI_DOUBLE_PRECISION, rcvbuf, mxcols*nlev, &
     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (myid == 0) flux_dn = transpose(rcvbuf)

  sendbuf = transpose(flux_dir_loc)
  call MPI_GATHER( &
     sendbuf, nlcol*nlev, MPI_DOUBLE_PRECISION, rcvbuf, mxcols*nlev, &
     MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  if (myid == 0) flux_dir = transpose(rcvbuf)

  ! only write the fluxes from rank 0:
  if (myid == 0) then

     call write_sw_fluxes(input_file, flux_up(:ncol,:), flux_dn(:ncol,:), flux_dir(:ncol,:))

  end if

  call MPI_FINALIZE(ierr)

!=========================================================================================
contains
!=========================================================================================

subroutine stop_on_err(error_msg)
  use iso_fortran_env, only : error_unit
  character(len=*), intent(in) :: error_msg

  if(error_msg /= "") then
    write (error_unit,*) trim(error_msg)
    write (error_unit,*) "rte_rrtmgp_test stopping"
    error stop 1
  end if
end subroutine stop_on_err

!=========================================================================================

subroutine read_atmos(fileName, ncol, tsi, p_lay, t_lay, p_lev, &
                      gas_vmr, mu0, asdir, asdif)

   character(len=*),                        intent(in ) :: fileName
   integer,                                 intent(out) :: ncol
   real(wp),                                intent(out) :: tsi
   real(wp), dimension(:,:),   allocatable, intent(out) :: p_lay, t_lay, p_lev
   real(wp), dimension(:,:,:), allocatable, intent(out) :: gas_vmr
   real(wp), dimension(:),     allocatable, intent(out) :: mu0, asdir, asdif
   ! -------------------
   integer :: igas
   integer :: ncid, nlay, nlev, pver, pverp
   real(wp), dimension(:,:,:,:), allocatable :: lay4d, lev4d
   real(wp), dimension(:,:,:), allocatable :: tmp3d
   real(wp), dimension(:,:), allocatable :: conc
   real(wp), dimension(:), allocatable :: a, b, alpha, beta
   ! ---------------------------------------------------------------------------

   if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

   ! The CAM dataset contains a hyperslab of columns extracted from a 4D array (lon,lat,lev,time).
   ! The variables in the NetCDF file have the shape (:,1,:,1).
   ! Allocate the state variables to include the "extra layer" that CAM uses for RRTMGP.
   ! This assumes that CAM's top is below 1 Pa.  Set extra layer values the same way
   ! that CAM does.
   ncol = get_dim_size(ncid, 'lon')
   pver = get_dim_size(ncid, 'lev')
   pverp = get_dim_size(ncid, 'ilev')
   nlay = pver + 1
   nlev = nlay + 1
   allocate(mu0(ncol), asdir(ncol), asdif(ncol))
   allocate(p_lev(ncol,nlev), p_lay(ncol,nlay), t_lay(ncol,nlay), conc(ncol,nlay))
   allocate(gas_vmr(ncol,nlay,ngas))
   allocate(a(ncol), b(ncol), alpha(ncol), beta(ncol))

   tsi = read_field(ncid, 'sol_tsi')

   tmp3d = read_field(ncid, 'rad_coszen', ncol, 1, 1)
   mu0 = tmp3d(:,1,1)

   tmp3d = read_field(ncid, 'rad_asdir', ncol, 1, 1)
   asdir = tmp3d(:,1,1)

   tmp3d = read_field(ncid, 'rad_asdif', ncol, 1, 1)
   asdif = tmp3d(:,1,1)

   lev4d = read_field(ncid, 'rad_pint', ncol, 1, pverp, 1)
   p_lev(:,2:nlev) = lev4d(:,1,:,1)
   ! top interface of extra layer set to 1.01 Pa
   p_lev(:,1) = 1.01_wp

   lay4d = read_field(ncid, 'rad_pmid', ncol, 1, pver, 1)
   p_lay(:,2:nlay) = lay4d(:,1,:,1)
   ! set top midpoint pressure the same way CAM does
   p_lay(:,1) = 0.5_wp * p_lev(:,2)

   lay4d = read_field(ncid, 'rad_temp', ncol, 1, pver, 1)
   t_lay(:,2:nlay) = lay4d(:,1,:,1)
   t_lay(:,1) = t_lay(:,2)

   ! factors used to compute the extra layer concentration of O3
   alpha = log(p_lev(:,2)/50._wp)
   beta = log(p_lay(:,2)/p_lev(:,2))/log(p_lay(:,2)/50._wp)
   a = ( (1._wp + alpha)*exp(-alpha) - 1._wp ) / alpha
   b = 1._wp - exp(-alpha)

   do igas = 1, ngas
      if(.not. var_exists(ncid, trim(cam_names(igas)))) &
         call stop_on_err("read_atmos: can't read concentration of " // trim(gas_names(igas)))
      lay4d = read_field(ncid, trim(cam_names(igas)), ncol, 1, pver, 1)
      conc(:,2:nlay) = lay4d(:,1,:,1)
      conc(:,1) = conc(:,2)

      ! convert conc to vmr
      select case (trim(cam_names(igas)))
      case ('rad_Q')
         ! specific humidity to mmr
         conc = conc / (1._wp - conc)
         ! mmr to vmr
         conc = conc * 1.607793_wp
      case ('rad_CO2')
         conc = conc * 0.658114_wp
      case ('rad_ozone')
         conc = conc * 0.603428_wp
         ! adjust "extra layer" value of ozone
         conc(:,1) = (conc(:,1) / (1._wp + beta)) * (a+b)
      case ('rad_N2O')
         conc = conc * 0.658090_wp
      case ('rad_CH4')
         conc = conc * 1.805423_wp
      case ('rad_O2')
         conc = conc * 0.905140_wp
      end select

      gas_vmr(:,:,igas) = conc

   end do

   ncid = nf90_close(ncid)

end subroutine read_atmos

!=========================================================================================

subroutine write_sw_fluxes(fileName, flux_up, flux_dn, flux_dir)

   character(len=*),         intent(in) :: fileName
   real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn, flux_dir
   ! -------------------
   integer :: ncid, ncol, nlev, nlon
   ! -------------------
   if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))

   ncol  = size(flux_up, dim=1)
   nlev  = size(flux_up, dim=2)

   ! add dimension for number of layer interfaces in RRTMGP grid
   call create_dim(ncid, 'plev_rad', nlev)

   ! check number of columns in CAM dataset
   nlon  = get_dim_size(ncid, 'lon')
   if (nlon /= ncol) then
      print*, 'ERROR: ncol, nlon= ', ncol, nlon
   end if

   call create_var(ncid, "sw_flux_up",  ["lon     ",  "plev_rad"], [ncol, nlev])
   call create_var(ncid, "sw_flux_dn",  ["lon     ",  "plev_rad"], [ncol, nlev])
   call create_var(ncid, "sw_flux_dir", ["lon     ",  "plev_rad"], [ncol, nlev])

   call stop_on_err(write_field(ncid, "sw_flux_up",  flux_up ))
   call stop_on_err(write_field(ncid, "sw_flux_dn",  flux_dn ))
   call stop_on_err(write_field(ncid, "sw_flux_dir", flux_dir))

   ncid = nf90_close(ncid)
end subroutine write_sw_fluxes

!=========================================================================================
end program cam_test
