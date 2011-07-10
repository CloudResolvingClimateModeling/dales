!> \file modsurfdata.f90
!! Variable definitions and auxilary routines for the surface model

!>
!! Variable definitions and auxilary routines for surface model
!>
!! This routine should have no dependency on any other routine, save perhaps modglobal or modfields.
!!  \author Thijs Heus, MPI-M
!!  \todo Documentation
!!  \par Revision list
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!



module modsurfdata

! implicit none

SAVE
  integer :: isurf        = -1            !<   Flag for surface parametrization

  ! Soil properties

  ! Domain-uniform properties
  integer, parameter  :: ksoilmax = 4       !<  Number of soil layers [-]

  real              :: lambdasat          !<  heat conductivity saturated soil [W/m/K]
  real              :: Ke                 !<  Kersten number [-]

  real, allocatable :: zsoil  (:)         !<  Height of bottom soil layer from surface [m]
  real, allocatable :: zsoilc (:)        !<  Height of center soil layer from surface [m]
  real, allocatable :: dzsoil (:)         !<  Depth of soil layer [m]
  real, allocatable :: dzsoilh(:)         !<  Depth of soil layer between center of layers [m]

  ! Spatially varying properties
  real, allocatable :: lambda  (:,:,:)    !<  Heat conductivity soil layer [W/m/K]
  real, allocatable :: lambdah (:,:,:)    !<  Heat conductivity soil layer half levels [W/m/K]
  real, allocatable :: lambdas (:,:,:)    !<  Soil moisture diffusivity soil layer 
  real, allocatable :: lambdash(:,:,:)    !<  Soil moisture diffusivity soil half levels 
  real, allocatable :: gammas  (:,:,:)    !<  Soil moisture conductivity soil layer 
  real, allocatable :: gammash (:,:,:)    !<  Soil moisture conductivity soil half levels 
  real, allocatable :: Dh      (:,:,:)    !<  Heat diffusivity
  real, allocatable :: phiw    (:,:,:)    !<  Water content soil matrix [-]
  real, allocatable :: phiwm   (:,:,:)    !<  Water content soil matrix previous time step [-]
  real, allocatable :: phifrac (:,:,:)    !<  Relative water content per layer [-]
  real              :: phiwav  (ksoilmax)
  real, allocatable :: phitot  (:,:)      !<  Total soil water content [-]
  real, allocatable :: pCs     (:,:,:)    !<  Volumetric heat capacity [J/m3/K]
  real, allocatable :: rootf   (:,:,:)    !<  Root fraction per soil layer [-]
  real              :: rootfav (ksoilmax)
  real, allocatable :: tsoil   (:,:,:)    !<  Soil temperature [K]
  real, allocatable :: tsoilm  (:,:,:)    !<  Soil temperature previous time step [K]
  real              :: tsoilav (ksoilmax)
  real, allocatable :: tsoildeep (:,:)    !<  Soil temperature [K]
  real              :: tsoildeepav

  real, allocatable :: swdavn  (:,:,:)
  real, allocatable :: swuavn  (:,:,:)
  real, allocatable :: lwdavn  (:,:,:)
  real, allocatable :: lwuavn  (:,:,:)

  integer           :: nradtime  = 60

  ! Soil related constants [adapted from ECMWF]
  real, parameter   :: phi       = 0.472  !<  volumetric soil porosity [-]
  real, parameter   :: phifc     = 0.323  !<  volumetric moisture at field capacity [-]
  real, parameter   :: phiwp     = 0.171  !<  volumetric moisture at wilting point [-]

  real, parameter   :: pCm       = 2.19e6 !<  Volumetric soil heat capacity [J/m3/K]
  real, parameter   :: pCw       = 4.2e6  !<  Volumetric water heat capacity [J/m3/K]

  real, parameter   :: lambdadry = 0.190  !<  Heat conductivity dry soil [W/m/K]
  real, parameter   :: lambdasm  = 3.11   !<  Heat conductivity soil matrix [W/m/K]
  real, parameter   :: lambdaw   = 0.57   !<  Heat conductivity water [W/m/K]

  real, parameter   :: bc        = 6.04     !< Clapp and Hornberger non-dimensional exponent [-]
  real, parameter   :: gammasat  = 0.57e-6  !< Hydraulic conductivity at saturation [m s-1]
  real, parameter   :: psisat    = -0.388   !< Matrix potential at saturation [m]

  ! Land surface properties

  ! Surface properties
  real, allocatable :: z0m        (:,:) !<  Roughness length for momentum [m]
  real              :: z0mav    = -1
  real, allocatable :: z0h        (:,:) !<  Roughness length for heat [m]
  real              :: z0hav    = -1
  real, allocatable :: tskin      (:,:) !<  Skin temperature [K]
  real, allocatable :: tskinm     (:,:) !<  Skin temperature previous timestep [K]
  real, allocatable :: Wl         (:,:) !<  Liquid water reservoir [m]
  real              :: Wlav     = -1                                  
  real, parameter   :: Wmax     = 0.0002 !<  Maximum layer of liquid water on surface [m]
  real, allocatable :: Wlm        (:,:) !<  Liquid water reservoir previous timestep [m]
  real, allocatable :: qskin      (:,:) !<  Skin specific humidity [kg/kg]
  real, allocatable :: albedo     (:,:) !<  Surface albedo [-]
  real              :: albedoav = -1
  real, allocatable :: LAI        (:,:) !<  Leaf area index vegetation [-]
  real              :: LAIav    = -1
  real, allocatable :: cveg       (:,:) !<  Vegetation cover [-]
  real              :: cvegav   = -1
  real, allocatable :: cliq       (:,:) !<  Fraction of vegetated surface covered with liquid water [-]
  real, allocatable :: Cskin      (:,:) !<  Heat capacity skin layer [J]
  real              :: Cskinav  = -1
  real, allocatable :: lambdaskin (:,:) !<  Heat conductivity skin layer [W/m/K]
  real              :: lambdaskinav
  real              :: ps       = -1    !<  Surface pressure [Pa]

  ! Surface energy balance
  real, allocatable :: Qnet     (:,:)   !<  Net radiation [W/m2]
  real              :: Qnetav   = -1
  real, allocatable :: LE       (:,:)   !<  Latent heat flux [W/m2]
  real, allocatable :: H        (:,:)   !<  Sensible heat flux [W/m2]
  real, allocatable :: G0       (:,:)   !<  Ground heat flux [W/m2]
  real, allocatable :: ra       (:,:)   !<  Aerodynamic resistance [s/m]
  real, allocatable :: rs       (:,:)   !<  Composite resistance [s/m]
  real, allocatable :: rsveg    (:,:)   !<  Vegetation resistance [s/m]
  real, allocatable :: rssoil   (:,:)   !<  Soil evaporation resistance [s/m]
  real              :: rsisurf2 = 0.    !<  Vegetation resistance [s/m] if isurf2 is used
  real, allocatable :: rsmin    (:,:)   !<  Minimum vegetation resistance [s/m]
  real              :: rsminav = -1
  real, allocatable :: rssoilmin(:,:)   !<  Minimum soil evaporation resistance [s/m]
  real              :: rssoilminav = -1
  real, allocatable :: tendskin (:,:)   !<  Tendency of skin [W/m2]
  real, allocatable :: gD       (:,:)   !<  Response factor vegetation to vapor pressure deficit [-]
  real              :: gDav

  ! Turbulent exchange variables
  logical           :: lmostlocal  = .false.  !<  Switch to apply MOST locally to get local Obukhov length
  logical           :: lsmoothflux = .false.  !<  Create uniform sensible and latent heat flux over domain
  logical           :: lneutral    = .false.  !<  Disable stability corrections
  real, allocatable :: obl   (:,:)      !<  Obukhov length [m]
  real              :: oblav            !<  Spatially averaged obukhov length [m]
  real, allocatable :: Cm    (:,:)      !<  Drag coefficient for momentum [-]
  real, allocatable :: Cs    (:,:)      !<  Drag coefficient for scalars [-]
  real, allocatable :: ustar (:,:)      !<  Friction velocity [m/s]
  real, allocatable :: thlflux (:,:)    !<  Kinematic temperature flux [K m/s]
  real, allocatable :: qtflux  (:,:)    !<  Kinematic specific humidity flux [kg/kg m/s]
  real, allocatable :: svflux  (:,:,:)  !<  Kinematic scalar flux [- m/s]

  ! Surface gradients of prognostic variables
  real, allocatable :: dudz  (:,:)      !<  U-wind gradient in surface layer [1/s]
  real, allocatable :: dvdz  (:,:)      !<  V-wind gradient in surface layer [1/s]
  real, allocatable :: dqtdz (:,:)      !<  Specific humidity gradient in surface layer [kg/kg/m]
  real, allocatable :: dthldz(:,:)      !<  Liquid water potential temperature gradient in surface layer [K/m]

  ! Surface properties in case of prescribed conditions (previous isurf 2, 3 and 4)
  real              :: thls  = -1       !<  Surface liquid water potential temperature [K]
  real              :: qts              !<  Surface specific humidity [kg/kg]
  real              :: thvs             !<  Surface virtual temperature [K]
  real, allocatable :: svs   (:)        !<  Surface scalar concentration [-]
  real              :: z0    = -1       !<  Surface roughness length [m]

  ! prescribed surface fluxes
  real              :: ustin  = -1      !<  Prescribed friction velocity [m/s]
  real              :: wtsurf = -1      !<  Prescribed kinematic temperature flux [K m/s]
  real              :: wqsurf = -1      !<  Prescribed kinematic moisture flux [kg/kg m/s]
  real              :: wsvsurf(100) = 0 !<  Prescribed surface scalar(n) flux [- m/s]

  ! Heterogeneous surfaces
  integer, parameter:: max_lands                   = 10 !<  Amount of land types that can be defined
  integer, parameter:: mpatch                      = 16 !<  Maximum amount of patches that can be defined in each direction
  logical           :: lhetero                 = .false.!<  Switch to apply heterogeneous surfaces using surface.inp.xxx
  logical           :: loldtable               = .false.!<  Switch to use surface.inp.xxx instead of updated surface.<name>.inp.xxx
  integer           :: xpatches                    =  2 !<  Amount of patches in the x-direction
  integer           :: ypatches                    =  1 !<  Amount of patches in the y-direction
  integer           :: land_use(mpatch,mpatch)     =  0 !<  Indicator for the land type
  integer           :: landtype(max_lands)         = -1 !< Type nr of the land in surface.inp.xxx
  character(len=10),dimension(max_lands) :: landname = "none" !< Name of the land type in surface.inp.xxx
  real              :: z0mav_land(max_lands)       = -1 !< Roughness length per land type for momentum [m]
  real              :: z0hav_land(max_lands)       = -1 !< Roughness length per land type for heat [m]
  real, allocatable :: z0mav_patch(:,:)                 !< Rougness length per patch
  real, allocatable :: z0hav_patch(:,:)                 !< Rougness length per patch
  real              :: thls_land(max_lands)        = -1 !< Surface liquid water potential temperature [K]
  real, allocatable :: thls_patch(:,:)                  !< Surface liquid water potential temperature [K]
  real, allocatable :: qts_patch(:,:)                   !< Surface liquid water potential temperature [K]
  real, allocatable :: thvs_patch(:,:)                  !< Surface virtual potential temperature [K]
  real              :: ps_land(max_lands)          = -1 !< Surface pressure [Pa]
  real, allocatable :: ps_patch(:,:)                    !< Surface pressure [Pa]
  real              :: ustin_land(max_lands)       = -1 !< Prescribed friction velocity [m/s]
  real, allocatable :: ustin_patch(:,:)                 !< Prescribed friction velocity [m/s]
  real              :: wt_land(max_lands)          = -1 !< Prescribed kinematic temperature flux [K m/s]
  real              :: wq_land(max_lands)          = -1 !< Prescribed kinematic moisture flux [kg/kg m/s]
  real              :: wsv_land(100,max_lands)     = -1 !< Prescribed surface scalar(n) flux [- m/s]
  real, allocatable :: wt_patch(:,:)                    !< Prescribed kinematic temperature flux [K m/s]
  real, allocatable :: wq_patch(:,:)                    !< Prescribed kinematic moisture flux [kg/kg m/s]
  real, allocatable :: wsv_patch(:,:,:)                 !< Prescribed surface scalar(n) flux [- m/s]
  real, allocatable :: rsisurf2_patch(:,:)              !< Vegetation resistance [s/m] if isurf2 is used
  real              :: rsisurf2_land(max_lands)    = 0  !< Vegetation resistance [s/m] if isurf2 is used
  real, allocatable :: albedo_patch(:,:)                !< Albedo
  real              :: albedo_land(max_lands)      = -1 !< Albedo
  real, allocatable :: tsoil_patch(:,:,:)               !< Soil temperature [K]
  real, allocatable :: tsoildeep_patch(:,:)             !< Soil temperature [K]
  real, allocatable :: phiw_patch(:,:,:)                !< 
  real, allocatable :: rootf_patch(:,:,:)               !< Root fraction per soil layer [-]
  real, allocatable :: Cskin_patch(:,:)                 !< Heat capacity skin layer [J]
  real, allocatable :: lambdaskin_patch(:,:)            !< Heat conductivity skin layer [W/m/K]
  real, allocatable :: Qnet_patch(:,:)                  !< Net radiation [W/m2]
  real, allocatable :: cveg_patch(:,:)                  !< Vegetation cover [-]
  real, allocatable :: Wl_patch(:,:)                    !< Liquid water reservoir [m]
  real, allocatable :: rsmin_patch(:,:)                 !< Minimum vegetation resistance [s/m]
  real, allocatable :: LAI_patch(:,:)                   !< Leaf area index vegetation [-]
  real, allocatable :: gD_patch(:,:)                    !< Response factor vegetation to vapor pressure deficit [-]
  real              :: tsoil_land(ksoilmax,max_lands)=-1!< Soil temperature [K]
  real              :: tsoildeep_land(max_lands)   = -1 !< Soil temperature [K]
  real              :: phiw_land(ksoilmax,max_lands) =-1!< 
  real              :: rootf_land(ksoilmax,max_lands)=-1!< Root fraction per soil layer [-]
  real              :: Cskin_land(max_lands)       = -1 !< Heat capacity skin layer [J]
  real              :: lambdaskin_land(max_lands)  = -1 !< Heat conductivity skin layer [W/m/K]
  real              :: Qnet_land(max_lands)        = -1 !< Net radiation [W/m2]
  real              :: cveg_land(max_lands)        = -1 !< Vegetation cover [-]
  real              :: Wl_land(max_lands)          = -1 !< Liquid water reservoir [m]
  real              :: rsmin_land(max_lands)       = -1 !< Minimum vegetation resistance [s/m]
  real              :: LAI_land(max_lands)         = -1 !< Leaf area index vegetation [-]
  real              :: gD_land(max_lands)          = -1 !< Response factor vegetation to vapor pressure deficit [-]
  real, allocatable :: oblpatch(:,:)                    !<  Obukhov length [m]

end module modsurfdata
