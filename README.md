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
     
  


## Selecting ESMF REgrid Option
To use ESMF regridding, we only need the setting:

        **Regridding Tool: "withESMF"** 
        
in the configuration file at run time.
