# ESMF Regridding in LIS

This document decribes the code modifications made to introduce ESMF regridding in LIS.




## New Directory: `esmf_regrid`
This directory contains three Fortran modules:

- **LIS_create_gridMod.F90**: utility functions for creating ESMF grids
- **LIS_field_bundleMod.F90**: utility subroutine for manipulating ESMF bundles
- **LIS_ESMF_Regrid_Utils.F90**: utility subroutines for creating the ESMF routehandle, performing regridding operations (at the level of ESMF field or ESMF bundle)

ESMF supports three main types of grids: uniform, rectilinear, and curvilinear. 
The LIS_create_gridMod.F90 file has a ESMF utility function that creates a rectilinear grid. 
It takes as arguments (among others) the latitude and longitude grid points. 
This function can be used for both an ESMF uniform grid and an ESMF rectilinear grid. 
Its applies for the regular lat/lon grid and the Gaussian grid.

## Modifications in the core Directories

Three files were changed:

- **LIS_PRIV_rcMod.F90**: Added the logical variable do_esmfRegridding as member variable of the lisrcdec derived type.
- **LIS_readConfig.F90**: Added code statements to check if ESMF regridding is selected in the configuration file. This allows to set do_esmfRegridding which default value is .FALSE.. 
- **LIS_coreMod.F90**: Added the following variables in the lis_domain_type derived type:

      type(ESMF_FieldBundle)     :: nldas2_bundle
      type(ESMF_FieldBundle)     :: merra2_bundle
      type(ESMF_FieldBundle)     :: gdas_bundle
      type(ESMF_FieldBundle)     :: gdasT1534_bundle
      type(ESMF_STAGGERLOC)      :: staggerloc



## Modifications in Met Fotcing Directories
For this work, we implemented the ESMF regridding tool on the following met forcing:

| Forcing Name | Forcing Grid Type | Model Grid Type |
| --- | --- | --- |
| nldas2       | regular lat-lon   | regular lat-lon |
| merra2       | regular lat-lon   | regular lat-lon |
| gdas         | gaussian          | regular lat-lon |
| gdasT1543    | gaussian          | gaussian        |

To describe the code changes, we will focus on the **merra2** case only. The same steps apply to all the met forcings too.

### File: merra2_forcingMod.F90
I first added the following variables in the merra2_type_dec derived type:

     type(ESMF_FieldBundle)       :: forcing_bundle
     type(ESMF_RouteHandle)       :: routehandle
     type(ESMF_DynamicMask)       :: dynamicMask
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_STAGGERLOC)        :: staggerloc
     type(ESMF_RegridMethod_Flag) :: regridMethod
     real                         :: undefined_value 
 
 In addition, I included the following module variables:
 
     integer                      :: num_merra2_fields       ! number of available fields
     character(len=100)           :: list_merra2_fields(30)  ! list of name of fields
     character(len=30), parameter :: merra2_bundle_bname = "merra2_bundle_"
  
I then wrote the subroutines:

- **set_list_merra2_fields**: Set the number of fields to be regridded and set a unique name of each of them.
- **create_merra2_Forcing_ESMFbundle**: Create the MERRA2 forcing ESMF grid and bundle. This subroutine might be called several times depending on the integration date. Howver, it will be called once for any period when the MERRA2 forcing resolution does not change.
- **create_merra2_Model_ESMFbundle**: Create the model ESMF grid and ESMF bundle. This subroutine is called once as the model resolution does not change. To create the model ESMF grid, the longitude and latitude grid points are coming directly from the arrays LIS_domain(n)%glon and LIS_domain(n)%glat.
- **create_merra2_ESMFroutehandle**: Determine the ESMF routehandle needed for thr regridding between the MERRA2 forcing and the model.

In the subroutine **init_merra2**, included the lines:

    IF (LIS_rc%do_esmfRegridding) THEN
       CALL set_list_merra2_fields(findex)
    ENDIF
    
before looping over the domains. Inside the loop, I added:

       IF (LIS_rc%do_esmfRegridding) THEN
          merra2_struc(n)%undefined_value = LIS_rc%udef

          if ((LIS_rc%met_interp(findex)) .eq. "bilinear") then
             merra2_struc(n)%regridMethod = ESMF_REGRIDMETHOD_BILINEAR
          elseif(trim(LIS_rc%met_interp(findex)) .eq. "neighbor") then
             merra2_struc(n)%regridMethod = ESMF_REGRIDMETHOD_NEAREST_STOD
          elseif(trim(LIS_rc%met_interp(findex)) .eq. "conservative") then
             merra2_struc(n)%regridMethod =  ESMF_REGRIDMETHOD_CONSERVE
          endif

          merra2_struc(n)%staggerloc = ESMF_STAGGERLOC_CENTER
          LIS_domain(n)%staggerloc   = ESMF_STAGGERLOC_CENTER

          forcing_gridDesc(:) = 0
          forcing_gridDesc(2) = gridDesci(n, 2) ! num points along x
          forcing_gridDesc(3) = gridDesci(n, 3) ! num points along y
          forcing_gridDesc(4) = gridDesci(n, 4) ! lower lat
          forcing_gridDesc(5) = gridDesci(n, 5) ! lower lon
          forcing_gridDesc(7) = gridDesci(n, 7) ! upper lat
          forcing_gridDesc(8) = gridDesci(n, 8) ! upper lon
          forcing_gridDesc(9) = gridDesci(n, 9) ! x-grid size
          forcing_gridDesc(10)= gridDesci(n, 10) ! y-grid size

          CALL create_merra2_Forcing_ESMFbundle(n, forcing_gridDesc(1:10))
          CALL create_merra2_Model_ESMFbundle(n)
          CALL create_merra2_ESMFroutehandle(n)
       ELSE
          ! Setting up weights for Interpolation
          ...
       ENDIF

### File: read_merra2.F90

I wrote a subroutine **performESMFregrid_merra2** that perform the ESMF regridding by spacially interpolating a MERRA2 field to the LIS running domain. The declaration of the subroutine is:

       subroutine performESMFregrid_merra2(n, findex, month, input_var, var_index, merraforc)

      integer, intent(in)    :: n
      integer, intent(in)    :: findex
      integer, intent(in)    :: month
      real,    intent(in)    :: input_var(merra2_struc(n)%ncold, merra2_struc(n)%nrold, 24)
      integer, intent(in)    :: var_index
      real,    intent(inout) :: merraforc(merra2_struc(n)%nvars, 24, LIS_rc%lnc(n)*LIS_rc%lnr(n))
      
 The subroutine **read_merra2** was modified as:
 
      IF (LIS_rc%do_esmfRegridding) THEN
        CALL performESMFregrid_merra2(n, findex, month,  tair, 1, merraforc)
        CALL performESMFregrid_merra2(n, findex, month,  qair, 2, merraforc)
        CALL performESMFregrid_merra2(n, findex, month, uwind, 5, merraforc)
        CALL performESMFregrid_merra2(n, findex, month, vwind, 6, merraforc)
        CALL performESMFregrid_merra2(n, findex, month,    ps, 7, merraforc)
     ELSE
        call interp_merra2_var(n,findex,month,tair,  1, .false., merraforc)
        call interp_merra2_var(n,findex,month,qair,  2, .false., merraforc)
        call interp_merra2_var(n,findex,month,uwind, 5, .false., merraforc)
        call interp_merra2_var(n,findex,month,vwind, 6, .false., merraforc)
        call interp_merra2_var(n,findex,month,ps,    7, .false., merraforc)
     ENDIF


## Selecting ESMF REgrid Option
To use ESMF regridding, we only need the setting:

  **Regridding Tool: "withESMF"** 
        
in the configuration file at run time.
