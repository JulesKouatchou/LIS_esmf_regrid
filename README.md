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

In ESMF, regridding is done at the ESMF field level or at the ESMF bundle level. We chose to use the field case (though we also implemented the bundle option) because we assume that the locations of missing values are dynamics (not known in advance). Our task is then to define ESMF fields (on both the source and destination grids) for all the variables to be regridded. To simplify the implementation, we group the source (forcing data) and destination (model domain) fields into source ESMF bundle and model ESMF bundle respectively. The ESMF bundles are place holders to easily move data around. 

**Options for Grid Types**

ESMF supports three <a href="https://www.dkrz.de/up/services/analysis/visualization/sw/ncl/docs/ncl-regridding-with-esmf?lang=en">types of grids</a>:

- Rectilinear grid: A 2-dimensional rectilinear grid have parallel grid axes, the x-axis values are monotonic increasing and perpendicular to the monotonic-increasing y-axis values. The x- and y-axis coordinates are 1-dimensional x(I) and y(I).
- Curvilinear grid: A curvilinear grid is characterized by curved coordinate lines. The x- and y-axis coordinates are 2-dimensional x(i,j) and y(i,j).
- Unstructured grid: An unstructured or irregular grid can have shapes such as triangle or tetrahedral in an irregular pattern defined by latitude and longitude vertices and the number of vertices for each cell. It also can represent irregular distributed points (point(x,y)). 

![fig_grid](https://slideplayer.com/slide/4799757/15/images/32/Grid+Types+Structured+Grids%3A+rectilinear+curvilinear+uniform+regular.jpg)

In this work, we focus on 2D regular lat-lon grid and gaussian grid. They can be seen as logically rectangular grids. We wrote a ESMF utility function that creates a ESMF rectilinear grid. The function takes as arguments (among other parameters) the longitude and latitude grid points at centers (that are predefined based on the type of grid) and at corners. We wrote routines that determines the grid points at corners based on the ones at centers.

In the process of creating the grid, it is important to indicate the coordinate system to use. It is used to indicate to other users the type of the coordinates, and also to control how the coordinates are interpreted in regridding methods. Three options are available:

   - `ESMF_COORDSYS_CART`: Cartesian coordinate system. The Cartesian coordinates are mapped to the grid coordinate dimensions in the following order: x, y, z.
   - `ESMF_COORDSYS_SPH_DEG`: Spherical coordinates in degrees. The spherical coordinates are mapped to the grid coordinate dimensions in the following order: longitude, latitude, radius. This is the option we will use here.
   - `ESMF_COORDSYS_SPH_RAD`: Spherical coordinates in radians. The spherical coordinates are mapped to the grid coordinate dimensions in the following order: longitude, latitude, radius.

**Options for Regridding Methods**

The _lis.config_ file contains the setting:

       Spatial interpolation method (met forcing):
       
 that determine the regridding method used to interpolate from the forcing data grid to the model grid. Its options are:
 
       bilinear, neighbor, budget-bilinear (for consevative)
       
 We use the same setting to determine which ESMF regridding option to choose:
 
- `ESMF_REGRIDMETHOD_BILINEAR`: Destination value is a linear combination of the source values in the cell which contains the destination point. The weights for the linear combination are based on the distance of destination point from each source value.
- `ESMF_REGRIDMETHOD_NEAREST_STOD`: Each destination point is mapped to the closest source point. A given source point may go to multiple destination points, but no destination point will receive input from more than one source point.
- `ESMF_REGRIDMETHOD_CONSERVE`: The main purpose of this method is to preserve the integral of the field between the source and destination. The value of a destination cell is calculated as the weighted sum of the values of the source cells that it overlaps. The weights are determined by the amount the source cell overlaps the destination cell. Needs corner coordinate values to be provided in the Grid. 

Depending on the regridding method considered, the tool will use one of the two lines types:

- `ESMF_LINETYPE_CART`: Specifies that the line between two points follows a straight path through the 3D Cartesian space in which the sphere is embedded. Distances are measured along this 3D Cartesian line. We can only use the bilinear and nearest neighbor methods here.
- `ESMF_LINETYPE_GREAT_CIRCLE`: Specifies that the line between two points follows a great circle path (the shortest path between two points on a sphere) along the sphere surface. We can only use the bilinear and conservative methods here.


We wrote a configuarable standalone code that performs the regridding of two synthetic fields. The source and destination fields can be either on the Gaussian or lat/lon grid. We observed that the settings (coordinate system and line type) on the table below were best (in terms of meeting the regridding accuracy) to carry out regridding:

| | Nearest neighbor | Bilinear | Conservative |
| --- | --- | --- | --- |
| **lat-lon to lat-lon** | ESMF_COORDSYS_CART | ESMF_COORDSYS_CART | ESMF_COORDSYS_SPH_DEG |
|  (nldas2, merra2)      | ESMF_LINETYPE_CART | ESMF_LINETYPE_CART | ESMF_LINETYPE_GREAT_CIRCLE |
| **Gaussian to lat-lon**  | ESMF_COORDSYS_CART | ESMF_COORDSYS_SPH_DEG | ESMF_COORDSYS_SPH_DEG |
|    (gdas)               | ESMF_LINETYPE_CART | ESMF_LINETYPE_GREAT_CIRCLE | ESMF_LINETYPE_GREAT_CIRCLE |
| **Gaussian to Gaussian** | ESMF_COORDSYS_CART |  ?        | ? |
|    (gdasT1534)       | ESMF_LINETYPE_CART |   ?        | ?  |

The above table gives us some guidelines on what type of grids we need to create as  function of the regridding method selected. In LIS, when the "budget-bilinear" regridding method is seleected, one set of fields is regridded using the bilinear method and another one using the conservative method. This causes us a challenge when the `nldas2` and `merra2` forcing data are used. In such a case, we might need to create two ESMF grids and compute two sets of weights in order to handle the two regridding methods.

**[Dynamic Masking](http://esmf-cu.colorado.edu/esmf_releases/last_built/ESMF_refdoc/node5.html#RH:DynMask)**

Once a RouteHandle object is available, whether it was created with or without static masking, the associated regrid operation can further be masking during RouteHandle execution . This is called dynamic masking, because it can dynamically change between subsequent RouteHandle executions. The RouteHandle itself remains unchange during this process. The dynamic masking information is processed on the fly as the RouteHandle is applied. We use such approach with the assumption that locations of missing values are not know in advance.


## New Directory: `esmf_regrid`
This directory contains three Fortran modules:

- **LIS_create_gridMod.F90**: utility functions for creating ESMF grids
- **LIS_field_bundleMod.F90**: utility subroutine for manipulating ESMF bundles and fields.
- **LIS_ESMF_Regrid_Utils.F90**: utility subroutines for creating the ESMF RouteHandle, performing regridding operations (at the level of ESMF field or ESMF bundle)

ESMF supports three main types of grids: uniform, rectilinear, and curvilinear. 
The `LIS_create_gridMod.F90` file has a ESMF utility function that creates a rectilinear grid. 
It takes as arguments (among others) the latitude and longitude grid points at the centers and 
latitude and longitude grid points at the corners. It also takes a logical flag (`periodic` optional argument)
to state if the grid being created is periodic or not. 
The grid will have the following features:

- The indexing is local: the local domain always start at (1,1).
- The grid points at centers and corners are attached to the grid coordinates.
- It can be used in the contect of conservative, bilinear and nearest neighbor regridding.

This function can be used for both an ESMF uniform grid and an ESMF rectilinear grid. 
It applies for the regular lat/lon grid and the Gaussian grid.

The grid points at the centers are provided indirectly by the LIS model. We wrote subroutines to determine
the latitude and longitude grid points at corners given the ones at the centers.

## Modifications in the `core` Directory

Two files were changed:

- **LIS_PRIV_rcMod.F90**: Added the logical variable `do_esmfRegridding` as member variable of the lisrcdec derived type.
- **LIS_readConfig.F90**: Added code statements to check if ESMF regridding is selected in the configuration file. This allows to set `do_esmfRegridding` which default value is .FALSE.. 


## Modifications in Met Forcing Directories
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

     type(ESMF_RouteHandle)       :: routehandle_bilinear
     type(ESMF_DynamicMask)       :: dynamicMask_bilinear
     type(ESMF_RegridMethod_Flag) :: regridMethod_bilinear = ESMF_REGRIDMETHOD_BILINEAR
     type(ESMF_RouteHandle)       :: routehandle_conserve
     type(ESMF_DynamicMask)       :: dynamicMask_conserve
     type(ESMF_RegridMethod_Flag) :: regridMethod_conserve = ESMF_REGRIDMETHOD_CONSERVE
     type(ESMF_RouteHandle)       :: routehandle_neighbor
     type(ESMF_DynamicMask)       :: dynamicMask_neighbor
     type(ESMF_RegridMethod_Flag) :: regridMethod_neighbor = ESMF_REGRIDMETHOD_NEAREST_STOD
     type(ESMF_TypeKind_Flag)     :: type_kind = ESMF_TYPEKIND_R4
     type(ESMF_Grid)              :: forcing_grid
     type(ESMF_Grid)              :: model_grid   
     type(ESMF_Field)             :: forcing_field  
     type(ESMF_Field)             :: model_field
     type(ESMF_Grid)              :: forcing_gridCS
     type(ESMF_Grid)              :: model_gridCS
     type(ESMF_Field)             :: forcing_fieldCS
     type(ESMF_Field)             :: model_fieldCS
     real                         :: undefined_value
 
  
I then wrote the subroutines:

- **create_merra2_Forcing_ESMFobjects**: Create the MERRA2 forcing ESMF grids and fields. This subroutine might be called several times if the resolution of the forcing data changes over time. However, it will be called once for any period when the MERRA2 forcing resolution does not change. Before creating the ESMF forcing data grid, the subroutine computes the lat-lon grid points at centers and corners using the gobal forcing data domain parameters (lat-lon corners, total number of grid points in each dimension). 
- **create_merra2_Model_ESMFobjects**: Create the model ESMF grids and ESMF fields. This subroutine is called once as the model resolution does not change. To create the model ESMF grids, the lat-lon grid points at centers are coming directly from the arrays LIS_domain(n)%glon and LIS_domain(n)%glat. The lat-lon points at corners are computes before the ESMF grid creations.
- **create_merra2_ESMFroutehandle**: Determine the ESMF RouteHandle(s) needed for the regridding between the MERRA2 forcing and the model.

In the subroutine **init_merra2**, included the lines inside the loops:

       IF (LIS_rc%do_esmfRegridding) THEN
          merra2_struc(n)%undefined_value = LIS_rc%udef

          forcing_gridDesc(:) = 0
          forcing_gridDesc(2) = gridDesci(n, 2) ! num points along x
          forcing_gridDesc(3) = gridDesci(n, 3) ! num points along y
          forcing_gridDesc(4) = gridDesci(n, 4) ! lower lat
          forcing_gridDesc(5) = gridDesci(n, 5) ! lower lon
          forcing_gridDesc(7) = gridDesci(n, 7) ! upper lat
          forcing_gridDesc(8) = gridDesci(n, 8) ! upper lon
          forcing_gridDesc(9) = gridDesci(n, 9) ! x-grid size
          forcing_gridDesc(10)= gridDesci(n, 10) ! y-grid size

          CALL create_merra2_Forcing_ESMFobjects(n, forcing_gridDesc(1:10))
          CALL create_merra2_Model_ESMFobjects(n)
          CALL create_merra2_ESMFroutehandle(n, findex)
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
        CALL performESMFregrid_merra2(n, findex, month,  tair, 1, .FALSE., merraforc)
        CALL performESMFregrid_merra2(n, findex, month,  qair, 2, .FALSE., merraforc)
        CALL performESMFregrid_merra2(n, findex, month, uwind, 5, .FALSE., merraforc)
        CALL performESMFregrid_merra2(n, findex, month, vwind, 6, .FALSE., merraforc)
        CALL performESMFregrid_merra2(n, findex, month,    ps, 7, .FALSE., merraforc)
     ELSE
        call interp_merra2_var(n, findex, month, tair,  1, .false., merraforc)
        call interp_merra2_var(n, findex, month, qair,  2, .false., merraforc)
        call interp_merra2_var(n, findex, month, uwind, 5, .false., merraforc)
        call interp_merra2_var(n, findex, month, vwind, 6, .false., merraforc)
        call interp_merra2_var(n, findex, month, ps,    7, .false., merraforc)
     ENDIF

and

     IF (LIS_rc%do_esmfRegridding) THEN
        CALL performESMFregrid_merra2(n, findex, month,  prectot, 8, .TRUE., merraforc)
     ELSE
        call interp_merra2_var(n, findex, month, prectot,  8, .true., merraforc)
     ENDIF


## GDAS Special Case
With the GDAS forcing data, the resulotion of the dataset changes over time. It is possible that for the same experiment, GDAS data will have various resolutions. The GDAS ESMF bundle and grid are therefore not created during intializations but each time there is a new set of data. The model ESMF bundle and grid are still created once at the initialization stage.

## Selecting ESMF Regridding Option at Runtime
To use ESMF regridding, we only need the setting:

   **Regridding Tool: "withESMF"** 
        
in the configuration file at run time.

## Initial Tests
We ran a series of of experiments (under various configurations) to test the efeectiveness of the ESMF regridding. All the model integrations are just for a couple of days and produced hourly outputs.
We use as reference the results obtained using the original version of the LIS code and we plot a time series of the errors (absolute maximum, absolote mean and rms).

|           | neighbor       | bilinear       | budget-bilinear |
| ---       | ---            | ---            | --- |
| nldas2    | Good (default) | N/A            | N/A |
| merra2    | issue with `Wind_f_tavg` | Good (default) |  issue with  `CRainf_f_tavg` and `Snow_f_tavg`  |
| gdas      | Good           | Good (default) | Good but issue with `Rainf_f_inst` |
| gdasT1543 | Good (default) |  issue with  `Rainf_f_inst` and `SWdown_f_inst`      |     |

|           |         | neighbor | bilinear | budget-bilinear |
| ---       | ---     | ---      | ---      | --- |
| nldas2 (4x4)| Regular | 43.49  |          |     |
|           | ESMF    |   33.21  |          |     |
| merra2 (1x2) | Regular | 126.33   | 104.56   | 105.11 |
|           | ESMF    | 148.88   | 132.90   | 144.00 |
| gdas (2x2) | Regular |  10.89 |  10.73   | 9.49 |
|           | ESMF    |    4.01 |   4.54  |  8.96 |
| gdasT1543 (4x4) | Regular | 2193.28 |          |      |
|           | ESMF    | 2208.37 |          |      |


## Obtaining the Code

The code is in a local Git repository on `discover` and can be obtained using the command:

       git clone /discover/nobackup/jkouatch/LIS_PROJECT/sourceCode/LISF
       
The following files were modified:

      make/plugins.py
      core/LIS_readConfig.F90
      core/LIS_coreMod.F90
      metforcing/gdas/get_gdas.F90
      metforcing/gdas/gdas_forcingMod.F90
      metforcing/gdas/read_gdas.F90
      metforcing/gdasT1534/gdasT1534_forcingMod.F90
      metforcing/gdasT1534/read_gdasT1534.F90
      metforcing/merra2/merra2_forcingMod.F90
      metforcing/merra2/read_merra2.F90
      metforcing/nldas2/nldas2_forcingMod.F90
      metforcing/nldas2/read_nldas2a.F90
      
A new directory, `esmf_regrid`, was created and it contains the files:

      esmf_regrid/LIS_create_gridMod.F90
      esmf_regrid/LIS_ESMF_Regrid_Utils.F90
      esmf_regrid/LIS_field_bundleMod.F90
      metforcing/nldas2/read_nldas2b.F90
      
      
      




