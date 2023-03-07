!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS) Coupler.
!*
!* FMS Coupler is free software: you can redistribute it and/or modify
!* it under the terms of the GNU Lesser General Public License as
!* published by the Free Software Foundation, either version 3 of the
!* License, or (at your option) any later version.
!*
!* FMS Coupler is distributed in the hope that it will be useful, but
!* WITHOUT ANY WARRANTY; without even the implied warranty of
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!* General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS Coupler.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atm_ice_wave_exchange_mod

  use mpp_mod,             only: mpp_pe, mpp_clock_id, mpp_clock_begin, mpp_clock_end, CLOCK_ROUTINE
  use mpp_domains_mod,     only: mpp_get_compute_domain
  use fms_mod,             only: clock_flag_default
  use constants_mod,       only: RADIUS
  use xgrid_mod,           only: xmap_type, setup_xmap, xgrid_count, stock_move, &
                                 put_to_xgrid, get_from_xgrid
  use time_manager_mod,    only: time_type
  use data_override_mod,   only: data_override
  use atmos_model_mod,     only: atmos_data_type
  use land_model_mod,      only: land_data_type
  use ice_model_mod,       only: ice_data_type
  use ocean_model_mod,     only: ocean_public_type
  use wave_model_mod,      only: wave_data_type, atmos_wave_boundary_type, ice_wave_boundary_type
  use stock_constants_mod, only: Lnd_stock, Ice_stock, ISTOCK_WATER, ISTOCK_SIDE
  use MOM_domains,         only: pass_vector, AGRID,pass_var,CGRID_NE

  implicit none
  private


  !---- exchange grid maps -----

  type(xmap_type), save :: xmap_atm_wav
  integer, save         :: n_xgrid_atm_wav
  type(xmap_type), save :: xmap_ice_wav
  integer, save         :: n_xgrid_ice_wav

  ! Exchange grid indices
  integer :: X2_GRID_ATM, X2_GRID_WAV

  public :: atm_wave_exchange_init, atm_to_wave, ice_wave_exchange_init, ice_to_wave

  integer :: cplClock, fluxLandIceClock
  logical :: do_runoff
  real    :: Dt_cpl
contains

  subroutine atm_wave_exchange_init(Atm, Wav, Atmos_wave_boundary)
    type(atmos_data_type),          intent(in)    :: Atm !< A derived data type to specify atmospheric boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(atmos_wave_boundary_type), intent(inout) :: atmos_wave_boundary !< A derived data type to specify properties
                                                                     !! passed from atmos to waves

    integer :: is, ie, js, je

    call setup_xmap(xmap_atm_wav, (/ 'ATM', 'WAV' /),       &
         (/ atm%Domain, Wav%Domain /),                    &
         "INPUT/grid_spec.nc", Atm%grid )
    ! exchange grid indices
    X2_GRID_ATM = 1; X2_GRID_WAV = 2;
    n_xgrid_atm_wav = max(xgrid_count(xmap_atm_wav),1)
    call mpp_get_compute_domain( Wav%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( atmos_wave_boundary%wavgrd_u10_mpp(is:ie,js:je,1) )
    allocate( atmos_wave_boundary%wavgrd_v10_mpp(is:ie,js:je,1) )

    atmos_wave_boundary%wavgrd_u10_mpp(:,:,:) = 0.0
    atmos_wave_boundary%wavgrd_v10_mpp(:,:,:) = 0.0

  end subroutine atm_wave_exchange_init

  subroutine ice_wave_exchange_init(Ice, Wav, Ice_wave_boundary)
    type(ice_data_type),          intent(in)      :: Ice !< A derived data type to specify ocean/ice boundary data
    type(wave_data_type),           intent(inout) :: Wav !< A derived data type to specify wave boundary data
    type(ice_wave_boundary_type), intent(inout)   :: Ice_wave_boundary !< A derived data type to specify properties
                                                                     !! passed from atmos to waves

    integer :: is, ie, js, je

    call setup_xmap(xmap_ice_wav, (/ 'WAV', 'OCN' /),       &
         (/ Wav%Domain, Ice%Domain /),                    &
         "INPUT/grid_spec.nc" )
    ! exchange grid indices
    n_xgrid_ice_wav = max(xgrid_count(xmap_ice_wav),1)
    call mpp_get_compute_domain( Wav%domain, is, ie, js, je )

    !allocate land_ice_boundary
    allocate( ice_wave_boundary%wavgrd_ucurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_ucurr_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%wavgrd_vcurr_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_vcurr_mpp(:,:,:) = 0.0

    allocate( wav%ustk0_mpp(is:ie,js:je) )
    wav%ustk0_mpp(:,:) = 0.0
    allocate( wav%vstk0_mpp(is:ie,js:je) )
    wav%vstk0_mpp(:,:) = 0.0
	
	!----------------shuoli202203---------------------
	allocate( wav%HS(is:ie,js:je) )
    wav%HS(:,:) = 0.0
	allocate( wav%K(is:ie,js:je) )
    wav%K(:,:) = 0.0
	allocate( wav%OMIGA(is:ie,js:je) )
    wav%OMIGA(:,:) = 0.0
	!----------------------------------------------

    allocate( wav%ustkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%ustkb_mpp(:,:,:) = 0.0
    allocate( wav%vstkb_mpp(is:ie,js:je,wav%num_stk_bands) )
    wav%vstkb_mpp(:,:,:) = 0.0

    ! This are a temporary and costly trick to make MPI work
    allocate( wav%glob_loc_X(is:ie,js:je) )
    wav%glob_loc_X(:,:) = 0
    allocate( wav%glob_loc_Y(is:ie,js:je) )
    wav%glob_loc_Y(:,:) = 0
    !----------------shuoli202111---------------------
    allocate( ice_wave_boundary%wavgrd_ithick_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_ithick_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%wavgrd_isize_mpp(is:ie,js:je,1) )
    ice_wave_boundary%wavgrd_isize_mpp(:,:,:) = 0.0
    
    allocate( wav%tauxice_mpp(is:ie,js:je) )
    wav%tauxice_mpp(:,:) = 0.0
    allocate( wav%tauyice_mpp(is:ie,js:je) )
    wav%tauyice_mpp(:,:) = 0.0
    !-------------------------------------------------
    call mpp_get_compute_domain( Ice%domain, is, ie, js, je )
    allocate( ice_wave_boundary%icegrd_ustk0_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_ustk0_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_vstk0_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_vstk0_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_ustkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_ustkb_mpp(:,:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_vstkb_mpp(is:ie,js:je,1,wav%num_stk_bands) )
    ice_wave_boundary%icegrd_vstkb_mpp(:,:,:,:) = 0.0
	!--------------------shuoli202203-------------------------------------
	allocate( ice_wave_boundary%icegrd_hs_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_hs_mpp(:,:,:) = 0.0
	allocate( ice_wave_boundary%icegrd_k_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_k_mpp(:,:,:) = 0.0
	allocate( ice_wave_boundary%icegrd_omiga_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_omiga_mpp(:,:,:) = 0.0
	!--------------------------------------------------------------
    !----------------shuoli202111---------------------
    allocate( ice_wave_boundary%icegrd_tauxice_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_tauxice_mpp(:,:,:) = 0.0
    allocate( ice_wave_boundary%icegrd_tauyice_mpp(is:ie,js:je,1) )
    ice_wave_boundary%icegrd_tauyice_mpp(:,:,:) = 0.0
	
	!-initiallize Ice%sCS%IST%fxw2i_str here, whole data_domain----
	allocate(Ice%sCS%IST%fxw2i_str(Ice%sCS%G%isdB:Ice%sCS%G%iedB,Ice%sCS%G%jsd:Ice%sCS%G%jed))
	Ice%sCS%IST%fxw2i_str(:,:) = 0.0
	allocate(Ice%sCS%IST%fyw2i_str(Ice%sCS%G%isd:Ice%sCS%G%ied,Ice%sCS%G%jsdB:Ice%sCS%G%jedB))
	Ice%sCS%IST%fyw2i_str(:,:) = 0.0
    !-------------------------------------------------

    return
  end subroutine ice_wave_exchange_init

  !> Does atmosphere TO wave operations (could do wave TO atmosphere operations too).
  subroutine atm_to_wave( Time, Atm, Wav, Atmos_wave_Boundary )
    type(time_type),                intent(in) :: Time !< Current time
    type(atmos_data_type),           intent(in) :: Atm
    type(wave_data_type),            intent(in) :: Wav
    type(atmos_wave_boundary_type), intent(inout):: Atmos_wave_Boundary

    real, dimension(n_xgrid_atm_wav) :: &
         ex_u_atm, &
         ex_v_atm

    integer :: remap_method

    remap_method = 1

    call put_to_xgrid (Atm%u_bot , 'ATM', ex_u_atm , xmap_atm_wav, remap_method=remap_method, complete=.false.)
    call put_to_xgrid (Atm%v_bot , 'ATM', ex_v_atm , xmap_atm_wav, remap_method=remap_method, complete=.true.)
    if (Wav%pe) then
       call get_from_xgrid(Atmos_Wave_Boundary%wavgrd_u10_mpp, 'WAV', ex_u_atm, xmap_atm_wav)
       call get_from_xgrid(Atmos_Wave_Boundary%wavgrd_v10_mpp, 'WAV', ex_v_atm, xmap_atm_wav)
    endif

  end subroutine atm_to_wave

  !> Does both ice TO wave and wave TO ice exchange grid operations.
  subroutine ice_to_wave( Time, Ice, Wav, Ice_wave_Boundary )
    type(time_type),                intent(in) :: Time !< Current time
    type(ice_data_type),            intent(inout) :: Ice !< The ice module container
    type(wave_data_type),           intent(in) :: Wav !< The wave module container
    type(ice_wave_boundary_type), intent(inout):: Ice_wave_Boundary !< The ice-wave boundary container

    real, dimension(n_xgrid_ice_wav) :: &
         ex_ucurr,  & ! Exchange grid x-current
         ex_vcurr,  & ! Exchange grid y-current
         ex_ustokes,& ! Exchange grid x-Stokes drift
         ex_vstokes,&   ! Exchange grid y-Stokes drift
         !----shuoli202111-------------------
         ex_ithick ,&   ! Exchange grid ice-to-wave ice thickness 
         ex_isize  ,&   ! Exchange grid ice-to-wave ice concentration
         ex_tauxice  ,&   ! Exchange grid wave-to-ice stress
         ex_tauyice  ,&   ! Exchange grid wave-to-ice peak direction
         !------------------------------------------    
		 !--------shuoli202203-------------------------------------
		 ex_hs,&   ! Exchange grid wave to ocean hs
		 ex_k,&   ! Exchange grid wave to ocean wave number
		 ex_omiga   ! Exchange grid wave to ocean peak angular freq.
		 !---------------------------------------------------------------
		 
    !----shuoli---
	real, dimension(:,:,:), allocatable :: ex_isize3d, ex_ithick3d, w2ipart_size, w2imH_ice
	real, dimension(:,:), allocatable :: ex_isize2d, ex_ithick2d, icefxw2i_str, icefyw2i_str
	!-------------
    integer :: remap_method ! Interpolation method (todo: list options)
    integer :: i_stk, i, j,is,ie,js,je,k
	!------------shuoli------
	real :: weights  ! A sum of the weights around a point.
    real :: I_wts    ! 1.0 / wts or 0 if wts is 0 [nondim].
	!--------------shuoli---
	!real, dimension(:,:), allocatable:: fxw2i  !< Zonal ice stress on ocean [Pa]
    !real, dimension(:,:), allocatable:: fyw2i  !< Meridional ice stress on ocean [Pa]
	!allocate(fxw2i(SZIB_(Ice%sCS%G),SZJ_(Ice%sCS%G)))
	!allocate(fyw2i(SZI_(Ice%sCS%G),SZJB_(Ice%sCS%G)))
	!fxw2i = 0.0
	!fyw2i = 0.0
	!--------------

    remap_method = 1
	
	!-----shuoli---
	call mpp_get_compute_domain( Ice%domain, is, ie, js, je )

	allocate(ex_isize3d(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed, 1))
	allocate(ex_ithick3d(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed, 1))
	ex_isize3d = 0.0
	ex_ithick3d = 0.0
	
	allocate(ex_isize2d(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed))
	allocate(ex_ithick2d(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed))
	ex_isize2d = 0.0
	ex_ithick2d = 0.0
	
	allocate(w2ipart_size(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed, 1:Ice%sCS%IG%CatIce))
	allocate(w2imH_ice(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed, 1:Ice%sCS%IG%CatIce))
	w2ipart_size = 0.0
	w2imH_ice = 0.0
	
	allocate(icefxw2i_str(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed))
	allocate(icefyw2i_str(Ice%sCS%G%isd:Ice%sCS%G%ied, Ice%sCS%G%jsd:Ice%sCS%G%jed))
    icefxw2i_str = 0.0
	icefyw2i_str = 0.0
	!---------------

    ! -> Put Ocean (ice) parameters onto exchange grid
    call put_to_xgrid (Ice%u_surf(:,:,:) , 'OCN', ex_ucurr , xmap_ice_wav)
    call put_to_xgrid (Ice%v_surf(:,:,:) , 'OCN', ex_vcurr , xmap_ice_wav)

    ! -> Only on wave-PEs, bring wave information off exchange grid
    if (Wav%pe) then
       call get_from_xgrid(Ice_Wave_Boundary%wavgrd_ucurr_mpp(:,:,1), 'WAV', ex_ucurr, xmap_ice_wav)
       call get_from_xgrid(Ice_Wave_Boundary%wavgrd_vcurr_mpp(:,:,1), 'WAV', ex_vcurr, xmap_ice_wav)
    endif

      ! -> Put Wave parameters onto exchange grid
    call put_to_xgrid (Wav%ustk0_mpp(:,:) , 'WAV', ex_ustokes , xmap_ice_wav)
    call put_to_xgrid (Wav%vstk0_mpp(:,:) , 'WAV', ex_vstokes , xmap_ice_wav)
    if (Ice%pe) then
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_ustk0_mpp(:,:,:), 'OCN', ex_ustokes, xmap_ice_wav)
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_vstk0_mpp(:,:,:), 'OCN', ex_vstokes, xmap_ice_wav)
    endif
    do i_stk = 1,wav%num_Stk_bands
      call put_to_xgrid (Wav%ustkb_mpp(:,:,i_stk) , 'WAV', ex_ustokes , xmap_ice_wav)
      call put_to_xgrid (Wav%vstkb_mpp(:,:,i_stk) , 'WAV', ex_vstokes , xmap_ice_wav)
      ! -> Only on ice-PEs, bring ice information off exchange grid
      if (Ice%pe) then
        call get_from_xgrid(Ice_Wave_Boundary%icegrd_ustkb_mpp(:,:,:,i_stk), 'OCN', ex_ustokes, xmap_ice_wav)
        call get_from_xgrid(Ice_Wave_Boundary%icegrd_vstkb_mpp(:,:,:,i_stk), 'OCN', ex_vstokes, xmap_ice_wav)
      endif

    enddo
	
	!%---------------------shuoli202203-------------------------------------
	! -> Wave parameters to ice to mom6, including wave number, hs, peak angular freq.
	! Hs,k,omiga
    call put_to_xgrid (Wav%HS(:,:) , 'WAV', ex_hs , xmap_ice_wav)
    call put_to_xgrid (Wav%K(:,:) , 'WAV', ex_k , xmap_ice_wav)
	call put_to_xgrid (Wav%OMIGA(:,:) , 'WAV', ex_omiga , xmap_ice_wav)
    if (Ice%pe) then
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_hs_mpp(:,:,:), 'OCN', ex_hs, xmap_ice_wav)
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_k_mpp(:,:,:), 'OCN', ex_k, xmap_ice_wav)
	  call get_from_xgrid(Ice_Wave_Boundary%icegrd_omiga_mpp(:,:,:), 'OCN', ex_omiga, xmap_ice_wav)
    endif
    !---------------------------------------------------------------------------------------
    !%-------------------shuoli202111----------------------------------------
       !--------ice to wave-----
    ! -> Put ice parameters onto exchange grid
	do k=1,Ice%sCS%IG%CatIce
	w2ipart_size(:,:,k) = Ice%sCS%IST%part_size(:,:,k)
	w2imH_ice(:,:,k) = Ice%sCS%IST%mH_ice(:,:,k)
	enddo
	
	!if (Ice%sCS%IST%valid_IST) then
     !  call pass_var(w2ipart_size, Ice%sCS%G%Domain)
	!   call pass_var(w2imH_ice, Ice%sCS%G%Domain, complete=.false.)
	!endif !---there's problem here
	
	ex_isize2d(:,:) = sum(w2ipart_size(:,:,:), DIM=3) ! sum fractional coverage over all thickness categories
    ex_ithick2d(:,:) = sum(w2imH_ice(:,:,:), DIM=3) * Ice%sCS%IG%H_to_kg_m2/917.0 
    ex_isize3d(:,:,1) = ex_isize2d(Ice%sCS%G%isd:Ice%sCS%G%ied,Ice%sCS%G%jsd:Ice%sCS%G%jed)
    ex_ithick3d(:,:,1) = ex_ithick2d(Ice%sCS%G%isd:Ice%sCS%G%ied,Ice%sCS%G%jsd:Ice%sCS%G%jed)
	!-------------
	!write(*,*)'computeD size, dataD size,GcomputeD', is, ie, js, je,Ice%sCS%G%isd,Ice%sCS%G%ied, Ice%sCS%G%jsd,Ice%sCS%G%jed,Ice%sCS%G%isc,Ice%sCS%G%iec,Ice%sCS%G%jsc,Ice%sCS%G%jec
    !---------
    call put_to_xgrid (ex_ithick3d(:,:,:) , 'OCN', ex_ithick , xmap_ice_wav)
    call put_to_xgrid (ex_isize3d(:,:,:) , 'OCN', ex_isize , xmap_ice_wav)

    ! -> Only on wave-PEs, bring ice information off exchange grid
    if (Wav%pe) then
       call get_from_xgrid(Ice_Wave_Boundary%wavgrd_ithick_mpp(:,:,1), 'WAV', ex_ithick, xmap_ice_wav)
       call get_from_xgrid(Ice_Wave_Boundary%wavgrd_isize_mpp(:,:,1), 'WAV', ex_isize, xmap_ice_wav)
    endif
	
	!----------
	!write(*,*)'from ice to wave ice con. ice thick: ',maxval(Ice_Wave_Boundary%wavgrd_ithick_mpp),maxval(Ice_Wave_Boundary%wavgrd_isize_mpp)
	!write(*,*)'part_size, mH_ice size', size(Ice%sCS%IST%part_size,DIM=1), size(Ice%sCS%IST%part_size,DIM=2), size(Ice%sCS%IST%part_size,DIM=3), &
	!size(Ice%sCS%IST%mH_ice,DIM=1), size(Ice%sCS%IST%mH_ice,DIM=2), size(Ice%sCS%IST%mH_ice,DIM=3),Ice%sCS%IG%CatIce,Ice%sCS%IG%NkIce,Ice%sCS%IG%NkSnow
	!--------
	
	
      !----------wave to ice-----------
      ! -> Put Wave parameters onto exchange grid
    call put_to_xgrid (Wav%tauxice_mpp(:,:) , 'WAV', ex_tauxice , xmap_ice_wav)
    call put_to_xgrid (Wav%tauyice_mpp(:,:) , 'WAV', ex_tauyice , xmap_ice_wav)
    
    if (Ice%pe) then!data on ice_compute_domain
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_tauxice_mpp(:,:,:), 'OCN', ex_tauxice, xmap_ice_wav)
      call get_from_xgrid(Ice_Wave_Boundary%icegrd_tauyice_mpp(:,:,:), 'OCN', ex_tauyice, xmap_ice_wav)
    endif
	

      icefxw2i_str(Ice%sCS%G%isc:Ice%sCS%G%iec,Ice%sCS%G%jsc:Ice%sCS%G%jec) = Ice_Wave_Boundary%icegrd_tauxice_mpp(:,:,1) * 1026.0
	  icefyw2i_str(Ice%sCS%G%isc:Ice%sCS%G%iec,Ice%sCS%G%jsc:Ice%sCS%G%jec) = Ice_Wave_Boundary%icegrd_tauyice_mpp(:,:,1) * 1026.0
	
	
	! set wave-to-ice stress, all introduced here
	call pass_vector(icefxw2i_str, icefyw2i_str, Ice%sCS%G%Domain, stagger=AGRID)
	
	do j=Ice%sCS%G%jsc-1,Ice%sCS%G%jec+1 ; do I=Ice%sCS%G%isc-1,Ice%sCS%G%iec
    weights = (Ice%sCS%G%areaT(i,j)*ex_isize2d(i,j) + Ice%sCS%G%areaT(i+1,j)*ex_isize2d(i+1,j))
	!weights_area = (Ice%sCS%G%areaT(i,j) + Ice%sCS%G%areaT(i+1,j))
    if (Ice%sCS%G%mask2dCu(I,j) * weights > 0.0) then ; I_wts = 1.0 / weights
      Ice%sCS%IST%fxw2i_str(I,j) = Ice%sCS%G%mask2dCu(I,j) * &
          (icefxw2i_str(i,j) * Ice%sCS%G%areaT(i,j)*ex_isize2d(i,j) + &
		  icefxw2i_str(i+1,j) * Ice%sCS%G%areaT(i+1,j)*ex_isize2d(i+1,j)) * I_wts ! 
	  !Ice%sCS%IST%fxw2i_str(I,j) = Ice%sCS%G%mask2dCu(I,j) * &
         ! (icefxw2i_str(i,j)  + icefxw2i_str(i+1,j) ) /2
         
	else
      Ice%sCS%IST%fxw2i_str(I,j) = 0.0
    end if
    enddo ; enddo
	
	do J=Ice%sCS%G%jsc-1,Ice%sCS%G%jec ; do i=Ice%sCS%G%isc-1,Ice%sCS%G%iec+1
    weights = (Ice%sCS%G%areaT(i,j)*ex_isize2d(i,j) + Ice%sCS%G%areaT(i,j+1)*ex_isize2d(i,j+1))
	!weights_area = (Ice%sCS%G%areaT(i,j) + Ice%sCS%G%areaT(i,j+1))
    if (Ice%sCS%G%mask2dCv(i,J) * weights > 0.0) then ; I_wts = 1.0 / weights
      Ice%sCS%IST%fyw2i_str(i,J) = Ice%sCS%G%mask2dCv(i,J) * &
          (icefyw2i_str(i,J) * Ice%sCS%G%areaT(i,J)*ex_isize2d(i,j) + &
		  icefyw2i_str(i,J+1) * Ice%sCS%G%areaT(i,J+1)*ex_isize2d(i,j+1)) * I_wts
	  !Ice%sCS%IST%fyw2i_str(i,J) = Ice%sCS%G%mask2dCv(i,J) * &
        !  (icefyw2i_str(i,J) + icefyw2i_str(i,J+1) ) /2
        
	else
      Ice%sCS%IST%fyw2i_str(i,J) = 0.0
    endif
    enddo ; enddo
	
	
	
	    if ( maxval(Ice%sCS%IST%fxw2i_str) >100 .or. maxval(Ice%sCS%IST%fyw2i_str) >100)then
	        write(*,*)'wave to ice stress x too large-->0 ',maxval(Ice%sCS%IST%fxw2i_str), maxval(Ice%sCS%IST%fyw2i_str)
	        !Ice%sCS%IST%fxw2i_str(i,j) = 0.0
			!Ice%sCS%IST%fyw2i_str(i,j) = 0.0
		end if
	
	
    !--------------------------------------------------------------

    return
  end subroutine ice_to_wave

end module atm_ice_wave_exchange_mod
