# ESMF Regridding in LIS

This document decribes the code modifications made to introduce ESMF regridding in LIS. The goal is to regrid from forcing data into the model domain. Internally in LIS, the forcing data grid and the model domain grid are represented as one dimentional arrays. For this work, we will create new ESMF two-dimensional (2D) grids to facilitate the use of ESMF. The regridding process will then be performed on the two new grids before the regridded data are passed to the LIS model domain.

## General Description of Steps

The ESMF regridding implemented here is broken in two stages:

- **Stage 1**: Generation of an interpolation weight matrix that describes how points in the source grid contribute to points in the destination grid. This is done through the creation of an ESMF RouteHandle.
- **Stage 2**: Multiplication of values on the source grid by the interpolation weight matrix to produce values on the destination grid. 

Typically, ESMF regridding follows the steps:

      Create the source ESMF grid
      Create the source ESMF field
      Create the destination ESMF grid
      Create the destination ESMF field
      Create the ESMF RouteHandle (call of ESMF_FieldRegridStore())
      
      Loop of the reading of the data
              Read the source data
              Populate the source ESMF field
              Perform ESMF regridding (call of ESMF_FieldRegrid())

In ESMF, regridding is done at the level of the ESMF field level or at the ESMF bundle level. We chose to use the field case (though we also implemented at the bundle option) because we assume that the locations of missing values are dynamics (not known in advance). Our task is then to define ESMF fields (on both the source and destination grids) for all the variables to be regridded. To simplify the implementation, we group the source (forcing data) and destination (model domain) fields into source ESMF bundle and model ESMF bundle respectively.

**Options for Grid Types**

ESMF supports three <a href="https://www.dkrz.de/up/services/analysis/visualization/sw/ncl/docs/ncl-regridding-with-esmf?lang=en">types of grids</a>:

- Rectilinear grid: A 2-dimensional rectilinear grid have parallel grid axes, the x-axis values are monotonic increasing and perpendicular to the monotonic-increasing y-axis values. The x- and y-axis coordinates are 1-dimensional x(I) and y(I).
- Curvilinear grid: A curvilinear grid is characterized by curved coordinate lines. The x- and y-axis coordinates are 2-dimensional x(i,j) and y(i,j).
- Unstructured grid: An unstructured or irregular grid can have shapes such as triangle or tetrahedral in an irregular pattern defined by latitude and longitude vertices and the number of vertices for each cell. It also can represent irregular distributed points (point(x,y)). 

[fig_grid](https://slideplayer.com/slide/4799757/15/images/32/Grid+Types+Structured+Grids%3A+rectilinear+curvilinear+uniform+regular.jpg)

In this work, we focus on 2D regular lat-lon grid and gaussian grid. They can be seen as logically rectangular grids. We wrote a ESMF utility function that creates a ESMF rectilinear grid. The function takes as arguments (among other parameters) the longitude and latitude grid points that are predefined based on the type of grid.

**Options for Regridding Methods**

The _lis.config_ file contains the setting:

       Spatial interpolation method (met forcing):
       
 that determine the regridding method used to interpolate from the forcing data grid to the model grid. Its options are:
 
       bilinear, neighbor, consevative
       
 We use the same setting to determine which ESMF regridding option to choose:
 
- `ESMF_REGRIDMETHOD_BILINEAR`: Destination value is a linear combination of the source values in the cell which contains the destination point. The weights for the linear combination are based on the distance of destination point from each source value.
- `ESMF_REGRIDMETHOD_NEAREST_STOD`: Each destination point is mapped to the closest source point. A given source point may go to multiple destination points, but no destination point will receive input from more than one source point.
- `ESMF_REGRIDMETHOD_CONSERVE`: The main purpose of this method is to preserve the integral of the field between the source and destination. The value of a destination cell is calculated as the weighted sum of the values of the source cells that it overlaps. The weights are determined by the amount the source cell overlaps the destination cell. Needs corner coordinate values to be provided in the Grid. 

**[Dynamic Masking](http://esmf-cu.colorado.edu/esmf_releases/last_built/ESMF_refdoc/node5.html#RH:DynMask)**

Once a RouteHandle object is available, whether it was created with or without static masking, the associated regrid operation can further be masking during RouteHandle execution . This is called dynamic masking, because it can dynamically change between subsequent RouteHandle executions. The RouteHandle itself remains unchange during this process. The dynamic masking information is processed on the fly as the RouteHandle is applied. We use such approach with the assumption that locations of missing values are not know in advance.


## New Directory: `esmf_regrid`
This directory contains three Fortran modules:

- **LIS_create_gridMod.F90**: utility functions for creating ESMF grids
- **LIS_field_bundleMod.F90**: utility subroutine for manipulating ESMF bundles
- **LIS_ESMF_Regrid_Utils.F90**: utility subroutines for creating the ESMF RouteHandle, performing regridding operations (at the level of ESMF field or ESMF bundle)

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



## Modifications in Met Focing Directories
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
- **create_merra2_Forcing_ESMFbundle**: Create the MERRA2 forcing ESMF grid and bundle. This subroutine might be called several times depending on the integration date. However, it will be called once for any period when the MERRA2 forcing resolution does not change. Before creating the ESMF forcing data grid, the subroutine computes the latitude and longitude grid points using the gobal forcing data domain parameters (lat-lon corners, total number of grid points in each dimension). 
- **create_merra2_Model_ESMFbundle**: Create the model ESMF grid and ESMF bundle. This subroutine is called once as the model resolution does not change. To create the model ESMF grid, the longitude and latitude grid points are coming directly from the arrays LIS_domain(n)%glon and LIS_domain(n)%glat.
- **create_merra2_ESMFroutehandle**: Determine the ESMF RouteHandle needed for thr regridding between the MERRA2 forcing and the model.

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


## Selecting ESMF Regridding Option at Runtime
To use ESMF regridding, we only need the setting:

  **Regridding Tool: "withESMF"** 
        
in the configuration file at run time.
