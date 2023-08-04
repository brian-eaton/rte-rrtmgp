program cam_test

  ! Read column(s) extracted from CAM's offline radiation driving dataset and produce fluxes.
  ! SW clear sky only.

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

  implicit none

  ! Local Variables

  ! Arrays: dimensions (col, lay)
  real(wp), dimension(:,:),   allocatable :: p_lay, t_lay, p_lev

  ! Shortwave only
  real(wp), dimension(:),     allocatable :: mu0
  real(wp), dimension(:,:),   allocatable :: sfc_alb_dir, sfc_alb_dif ! First dimension is band
  !
  ! Source functions
  !   Shortwave
  real(wp), dimension(:,:), allocatable, save :: toa_flux

  ! Output variables
  real(wp), dimension(:,:), target, &
                            allocatable :: flux_up, flux_dn, flux_dir

  ! Derived types from the RTE and RRTMGP libraries
  type(ty_gas_optics_rrtmgp)  :: k_dist
  type(ty_gas_concs)          :: gas_concs, gas_concs_in
  type(ty_optical_props_2str) :: atmos
  type(ty_fluxes_broadband)   :: fluxes

  ! Inputs to RRTMGP
  logical :: top_at_1

  integer  :: ncol, nlay, nbnd, ngpt
  integer  :: igas
  integer  :: nUserArgs = 0

  integer, parameter :: ngas = 6           ! no data for 'co ' and 'n2 '
  character(len=3), dimension(ngas) &
                     :: gas_names = ['h2o', 'co2', 'o3 ', 'n2o', 'ch4', 'o2 ']
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
  !   Arrays are allocated as they are read
  call read_atmos(input_file, ncol,     &
                  p_lay, t_lay, p_lev,  &
                  gas_concs_in)

  nlay = size(p_lay, 2)

  call stop_on_err(gas_concs%init(gas_names))
  do igas = 1, ngas
    call vmr_2d_to_1d(gas_concs, gas_concs_in, gas_names(igas), size(p_lay, 1), nlay)
  end do

  ! load data into classes
  call load_and_init(k_dist, k_dist_file, gas_concs)

  ! Problem sizes
  nbnd = k_dist%get_nband()
  ngpt = k_dist%get_ngpt()
  top_at_1 = p_lay(1, 1) < p_lay(1, nlay)

  ! Allocate arrays for the optical properties themselves.
  call stop_on_err(atmos%alloc_2str( ncol, nlay, k_dist))

  !  Boundary conditions depending on whether the k-distribution being supplied
  allocate(toa_flux(ncol, ngpt))
  allocate(sfc_alb_dir(nbnd, ncol), sfc_alb_dif(nbnd, ncol), mu0(ncol))

  ! Ocean-ish values for no particular reason
  sfc_alb_dir = 0.06_wp
  sfc_alb_dif = 0.06_wp
  mu0 = .86_wp

  ! Fluxes
  allocate(flux_up(ncol,nlay+1), flux_dn(ncol,nlay+1))
  allocate(flux_dir(ncol,nlay+1))

  call rte_config_checks(logical(.false., wl))

  ! SW Solver
  fluxes%flux_up => flux_up(:,:)
  fluxes%flux_dn => flux_dn(:,:)
  fluxes%flux_dn_dir => flux_dir(:,:)

  call stop_on_err(k_dist%gas_optics(p_lay, p_lev, &
                                     t_lay,        &
                                     gas_concs,    &
                                     atmos,        &
                                     toa_flux))

  call stop_on_err(rte_sw(atmos, top_at_1, &
                          mu0,   toa_flux, &
                          sfc_alb_dir, sfc_alb_dif, &
                          fluxes))

  call write_sw_fluxes(input_file, flux_up, flux_dn, flux_dir)

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

subroutine vmr_2d_to_1d(gas_concs, gas_concs_in, name, sz1, sz2)
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_rte_kind,           only: wp

  type(ty_gas_concs), intent(in)    :: gas_concs_in
  type(ty_gas_concs), intent(inout) :: gas_concs
  character(len=*),   intent(in)    :: name
  integer,            intent(in)    :: sz1, sz2

  real(wp) :: tmp(sz1, sz2), tmp_col(sz2)

  call stop_on_err(gas_concs_in%get_vmr(name, tmp))

  tmp_col(:) = tmp(1, :)

  call stop_on_err(gas_concs%set_vmr(name, tmp_col))

end subroutine vmr_2d_to_1d

!=========================================================================================

subroutine read_atmos(fileName, ncol, p_lay, t_lay, p_lev, gas_concs)

   character(len=*),    intent(in   ) :: fileName
   integer,             intent(out  ) :: ncol
   real(wp), dimension(:,:), allocatable, &
                        intent(inout) :: p_lay, t_lay, p_lev
   type(ty_gas_concs), intent(inout) :: gas_concs
   ! -------------------
   integer :: igas
   integer :: ncid, nlay, nlev, nlat, nlon
   real(wp), dimension(:,:,:,:), allocatable :: lay4d, lev4d
   real(wp), dimension(:,:), allocatable :: conc
   ! ---------------------------------------------------------------------------

   if(nf90_open(trim(fileName), NF90_NOWRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("read_atmos: can't find file " // trim(fileName))

   ! CAM dataset has a single column extracted from a 4D array (lon,lat,lev,time).
   ! Profile variables have shape (1,1,:,1).
   ncol = 1
   nlay = get_dim_size(ncid, 'lev')
   nlev = get_dim_size(ncid, 'ilev')
   if(nlev /= nlay+1) call stop_on_err("read_atmos: nlev should be nlay+1")

   lay4d = read_field(ncid, 'rad_pmid', 1, 1, nlay, 1)
   p_lay = reshape(lay4d, [1, nlay])

   lay4d = read_field(ncid, 'rad_temp', 1, 1, nlay, 1)
   t_lay = reshape(lay4d, [1, nlay])

   lev4d = read_field(ncid, 'rad_pint', 1, 1, nlev, 1)
   p_lev = reshape(lev4d, [1, nlev])

   call stop_on_err(gas_concs%init(gas_names))
   do igas = 1, ngas
      if(.not. var_exists(ncid, trim(cam_names(igas)))) &
         call stop_on_err("read_atmos: can't read concentration of " // trim(gas_names(igas)))
      lay4d = read_field(ncid, trim(cam_names(igas)), 1, 1, nlay, 1)
      conc = reshape(lay4d, [1, nlay])

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
         case ('rad_N2O')
            conc = conc * 0.658090_wp
         case ('rad_CH4')
            conc = conc * 1.805423_wp
         case ('rad_O2')
            conc = conc * 0.905140_wp
         end select

      call stop_on_err(gas_concs%set_vmr(trim(gas_names(igas)), conc))
   end do

   ncid = nf90_close(ncid)

end subroutine read_atmos

!=========================================================================================

subroutine write_sw_fluxes(fileName, flux_up, flux_dn, flux_dir)
   character(len=*),         intent(in) :: fileName
   real(wp), dimension(:,:), intent(in) :: flux_up, flux_dn, flux_dir
   ! -------------------
   integer :: ncid, ncol, nlev, ilev
   ! -------------------
   if(nf90_open(trim(fileName), NF90_WRITE, ncid) /= NF90_NOERR) &
      call stop_on_err("write_fluxes: can't open file " // trim(fileName))
   !
   ! At present these dimension sizes aren't used
   !   We could certainly check the array sizes against these dimension sizes
   !
   ncol  = size(flux_up, dim=1)
   nlev  = size(flux_up, dim=2)
   ! number of layer interfaces in CAM dataset
   ilev  = get_dim_size(ncid, 'ilev')
   if (nlev /= ilev) then
      print*, 'ERROR: nlev, ilev= ', nlev, ilev
   end if
   call create_dim(ncid, "col_flx", ncol)

   call create_var(ncid, "sw_flux_up",  ["col_flx",  "ilev   "], [ncol, nlev])
   call create_var(ncid, "sw_flux_dn",  ["col_flx",  "ilev   "], [ncol, nlev])
   call create_var(ncid, "sw_flux_dir", ["col_flx",  "ilev   "], [ncol, nlev])

   call stop_on_err(write_field(ncid, "sw_flux_up",  flux_up ))
   call stop_on_err(write_field(ncid, "sw_flux_dn",  flux_dn ))
   call stop_on_err(write_field(ncid, "sw_flux_dir", flux_dir))

   ncid = nf90_close(ncid)
end subroutine write_sw_fluxes

!=========================================================================================
end program cam_test
