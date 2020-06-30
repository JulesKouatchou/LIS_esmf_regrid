      program test_gaussian_latlon

      use ESMF
   
      implicit none

      integer :: LIS_localPet, npets, rc, localrc, userrc
      character(len=ESMF_MAXSTR) :: cname1, cname2, cplname
      type(ESMF_VM):: vm

      type(ESMF_STAGGERLOC) :: stagger_loc = ESMF_STAGGERLOC_CENTER
      type(ESMF_Grid) :: gdas_grid, model_grid
      type(ESMF_FIELD) :: model_field, gdas_field
      character(len=20) :: gridname = "gdas grid"
      real             :: forcing_gridDesc(10)
      real             :: model_gridDesc(10)
      Integer :: x_procs = 1
      Integer :: y_procs = 1
      logical :: LIS_masterproc
      integer :: ftc(2), ftlb(2), ftub(2), i, j
      real, pointer :: farray2dd(:,:), farray2d_mod(:,:)
      type(ESMF_ArraySpec)         :: arrayspec
      type(ESMF_RouteHandle)       :: routehandle
      type(ESMF_DynamicMask)       :: dynamicMask
      type(ESMF_RegridMethod_Flag) :: regridMethod = ESMF_REGRIDMETHOD_BILINEAR
      real                         :: undefined_value = -999.0
      integer :: srcTerm

      real, pointer    :: lat_points(:), slat(:), lat_weights(:)
      real, pointer    :: lon_points(:)
      real             :: dy, dx, min_lon, min_lat, max_lat
      integer          :: num_lons, num_lats, ic
      real, parameter  :: pi=3.14159265358979
      real, parameter  :: dpr=180.0/pi
      integer, pointer :: factorIndexList(:,:)
      real(ESMF_KIND_R8), pointer :: factorList(:)

!EOP
!------------------------------------------------------------------------------
      !--------------
      ! ESMF Settings
      !--------------
      ! Initialize framework and get back default global VM
      call ESMF_Initialize(vm=vm, defaultlogfilename="FieldRegridSTest.Log", &
                          logkindflag=ESMF_LOGKIND_MULTI, rc=localrc)
      call check_error(localrc)

      call ESMF_LogSet (flush=.true.)
      ! Get number of PETs we are running with
      call ESMF_VMGet(vm, petCount=npets, localPet=LIS_localPet, rc=localrc)
      call check_error(localrc)

      LIS_masterproc = .FALSE.
      IF (LIS_localPet == 0 ) LIS_masterproc = .TRUE.

      ! Check for correct number of PETs
      if ( npets < x_procs*y_procs ) then
         call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
      endif

      !---------------------------------------
      ! Create the forcing data grid and field
      !---------------------------------------
      forcing_gridDesc = 0
      forcing_gridDesc(1) = 4
      forcing_gridDesc(2) = 512
      forcing_gridDesc(3) = 256
      forcing_gridDesc(4) = 89.463
      forcing_gridDesc(5) = 0
      forcing_gridDesc(6) = 128
      forcing_gridDesc(7) = -89.463
      forcing_gridDesc(8) = -0.703125
      forcing_gridDesc(9) =  0.703125
      forcing_gridDesc(10) = 128

      num_lons = forcing_gridDesc(2)
      num_lats = forcing_gridDesc(3)
      min_lat  = forcing_gridDesc(4)
      min_lon  = forcing_gridDesc(5)
      max_lat  = forcing_gridDesc(7)
      dx       = forcing_gridDesc(9)

      ALLOCATE(lon_points(num_lons))
      do ic = 1, num_lons
         lon_points(ic) = (ic-1)*dx + min_lon
         ! values must be between -180 to 180
         IF (lon_points(ic) > 180.0) lon_points(ic) = lon_points(ic) - 360.0
      enddo

      ! Determine the global Gaussian latitude grid points

      ALLOCATE(lat_points(num_lats))
      ALLOCATE(slat(num_lats))
      ALLOCATE(lat_weights(num_lats))

      call gausslat(num_lats, slat, lat_weights)
      do ic = 1, num_lats
         lat_points(ic) = dpr*asin(slat(ic))
         !lat_points(ic) = dpr*asin(slat(num_lats-ic+1))
      enddo

      gdas_grid = create_rectilinear_grid(lon_points, lat_points, &
                               "GDAS Grid", x_procs, y_procs,  &
                                staggerloc = stagger_loc)

      DEALLOCATE(lon_points, lat_points)
      DEALLOCATE(slat, lat_weights)

      call ESMF_GridWriteVTK(gdas_grid, stagger_loc, "forcing_gaussian", rc=rc)

      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R4)

      gdas_field = ESMF_FieldCreate(gdas_grid, arrayspec, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "ForcingPressure", rc=rc)
      call check_error(rc)

      call ESMF_FieldGet(gdas_field, localDe=0, farrayPtr=farray2dd, &
                   totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
        print"(a14,7i4)", "ForcingGrid: ", LIS_localPet,ftlb,ftub,ftc
      call check_error(rc)

      farray2dd(:,:) = LIS_localPet ! 100.00
      !farray2dd(ftlb(1):ftlb(1)+1,ftub(2)-3:ftub(2)) = undefined_value

      print*,"Forcing Pressure: ", LIS_localPet, minval(farray2dd), maxval(farray2dd)

      !--------------------------------
      ! Create the model grid and field
      !--------------------------------
      model_gridDesc = 0.0
      model_gridDesc(2) = 45.00   ! LIS_rc%gridDesc(n,2) ! Global num points along x
      model_gridDesc(3) = 33.00   ! LIS_rc%gridDesc(n,3) ! Global num points along y
      model_gridDesc(4) = 32.875 ! lower lat
      model_gridDesc(5) = -104.875 ! lower lon
      model_gridDesc(7) =  40.875 ! 7) ! upper lat
      model_gridDesc(8) = -93.875 ! 8) ! upper lon
      model_gridDesc(9) = 0.25 ! x-grid size
      model_gridDesc(10)= 0.25 ! y-grid size

      num_lons = model_gridDesc(2)
      num_lats = model_gridDesc(3)
      min_lat  = model_gridDesc(4)
      min_lon  = model_gridDesc(5)
      max_lat  = model_gridDesc(7)
      dx       = model_gridDesc(9)
      dy       = model_gridDesc(10)

      ALLOCATE(lon_points(num_lons))
      do ic = 1, num_lons
         lon_points(ic) = (ic-1)*dx + min_lon
         ! values must be between -180 to 180
         IF (lon_points(ic) > 180.0) lon_points(ic) = lon_points(ic) - 360.0
      enddo

      ALLOCATE(lat_points(num_lats))
      do ic = 1, num_lats     
         lat_points(ic) = (ic-1)*dy + min_lat
      enddo

      model_grid = create_rectilinear_grid(lon_points, lat_points, &
                               "Model Grid", x_procs, y_procs,  &
                                staggerloc = stagger_loc)
      print*,"Create model grid: ", LIS_localPet

      DEALLOCATE(lon_points, lat_points)

      call ESMF_GridWriteVTK(model_grid, stagger_loc, "model_latlon", rc=rc)

      call ESMF_ArraySpecSet(arrayspec, rank=2, typekind=ESMF_TYPEKIND_R4)

      model_field = ESMF_FieldCreate(model_grid, arrayspec, &
                              totalLWidth=(/0,0/), totalUWidth=(/0,0/), &
                              name = "ModelPressure", rc=rc)
      call check_error(rc)
      print*,"Create model field: ", LIS_localPet

      call ESMF_FieldGet(model_field, localDe=0, farrayPtr=farray2d_mod, &
          totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
      call check_error(rc)
      farray2d_mod = -99999.0
      print"(a20,2f15.4)", "Model Data Before: ", minval(farray2d_mod), maxval(farray2d_mod)

      !----------------------------
      ! Create the ESMF RouteHandle
      !----------------------------
      srcTerm = 0
      call ESMF_FieldRegridStore(srcField          = gdas_field, &
                                 dstField          = model_field, &
                                 srcTermProcessing = srcTerm, &
                                 srcMaskValues     = [0], &
                                 unmappedAction    = ESMF_UNMAPPEDACTION_IGNORE, &
                                 !lineType          = ESMF_LINETYPE_GREAT_CIRCLE, &
                                 factorList=factorList, factorIndexList=factorIndexList, &
                                 routeHandle       = routehandle, &
                                 regridmethod      = regridMethod, rc=rc)
      call check_error(rc)
      print*, LIS_localPet, "  ---> Done with ESMF_FieldRegridStore "

      print*," factorList:      ", LIS_localPet, size(factorList)
      print*," factorIndexList: ", LIS_localPet, size(factorIndexList, 1), size(factorIndexList, 2)

      call ESMF_DynamicMaskSetR4R8R4(dynamicMask, &
                                     simpleDynMaskProc, &
                                     dynamicSrcMaskValue=undefined_value, rc=rc)
      call check_error(rc)
      print*, LIS_localPet, "  ---> Done with ESMF_DynamicMaskSetR4R8R4 "

      !----------------------------
      ! Perform the ESMF Regridding
      !----------------------------
    
       call ESMF_FieldRegrid(srcField      = gdas_field,        &
                             dstField      = model_field,        &
                             routehandle   = routehandle,     &
                             dynamicMask   = dynamicMask,     &
                             zeroregion    = ESMF_REGION_SELECT, &
                             termorderflag = ESMF_TERMORDER_SRCSEQ, &
                             checkflag     = .TRUE., &
                             rc=rc)
      call check_error(rc)
      print*, LIS_localPet, "  ---> Done with ESMF_FieldRegrid "

      call ESMF_FieldGet(model_field, localDe=0, farrayPtr=farray2d_mod, &
          totalLBound=ftlb, totalUBound=ftub, totalCount=ftc, rc=rc)
      call check_error(rc)
      print*, LIS_localPet, "******* Model Data: ", minval(farray2d_mod), maxval(farray2d_mod)


!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------
!BOP
   function create_rectilinear_grid(lon_points, lat_points, grid_name, Nx, Ny, &
                                 staggerloc) result(new_grid)
      integer,          intent(in)    :: Nx ! number of PEs along longitudes
      integer,          intent(in)    :: Ny ! number of PEs along latitudes
      character(len=*), intent(in)    :: grid_name
      real,             intent(in)    :: lat_points(:)
      real,             intent(in)    :: lon_points(:)
      type(ESMF_STAGGERLOC), optional :: staggerloc
!
! !RETURNED VALUE:
      type(ESMF_Grid)                 :: new_grid
!
! !DESCRIPTION:
! Create a generic ESMF 2D rectilinear grid given arbitrary latitude and
! longitudes grid points.
!
! !LOCAL VARIABLES:
      integer                         :: i, j, num_lats, num_lons
      integer                         :: tlb(1), tub(1)
      integer                         :: rc, STATUS
      real(ESMF_KIND_R8), pointer     :: coordX(:)
      real(ESMF_KIND_R8), pointer     :: coordY(:)
      type(ESMF_STAGGERLOC)           :: staggerloc_
      integer                         :: countsPerDEDimX(Nx)
      integer                         :: countsPerDEDimY(Ny)
      CHARACTER(len=ESMF_MAXSTR)      :: IAm = "create_rectilinear_grid"
!EOP
!------------------------------------------------------------------------------
!BOC
      rc = ESMF_SUCCESS

      !write(LIS_logunit,*) '[INFO] Create the ESMF Rectilinear Grid: '// TRIM(grid_name)

      ! Determine the number of grid points to distributed to processors
      !-----------------------------------------------------------------
      num_lons = SIZE(lon_points)
      call decomposeDim(num_lons, countsPerDEDimX, Nx )

      num_lats = SIZE(lat_points)
      call decomposeDim(num_lats, countsPerDEDimY, Ny )

      ! Create the ESMF Grid
      !---------------------
      new_grid = ESMF_GridCreateNoPeriDim( &
                          ! Define an irregular distribution
                          countsPerDEDim1 = countsPerDEDimX, &
                          countsPerDEDim2 = countsPerDEDimY, &
                          ! Specify mapping of coords dim to Grid dim
                          coordDep1       = (/1/), & ! 1st coord is 1D and depends on 1st Grid dim
                          coordDep2       = (/2/), & ! 2nd coord is 2D and depends on 2nd Grid dim
                          indexflag       = ESMF_INDEX_GLOBAL, &
                          coordSys        = ESMF_COORDSYS_CART, & ! use cartesian coordinates
                          name            = TRIM(grid_name), rc=rc)

      call check_error(rc)

      !-------------------------------------------------------------------
      ! Allocate coordinate storage and associate it with the center
      ! stagger location.  Since no coordinate values are specified in
      ! this call no coordinate values are set yet.
      !-------------------------------------------------------------------
      staggerloc_ = ESMF_STAGGERLOC_CENTER
      if (PRESENT(staggerloc)) staggerloc_ = staggerloc

      call ESMF_GridAddCoord(new_grid, staggerloc=staggerloc_, rc=rc)
      call check_error(rc)

      !-------------------------------------------------------------------
      ! Get the pointer to the first coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 1, &
                                       farrayPtr           = coordX, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call check_error(rc)

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the first dimension.
      !-------------------------------------------------------------------
      do i = tlb(1), tub(1)
         coordX(i) = lon_points(i)    ! longitude grid points
      enddo

      !-------------------------------------------------------------------
      ! Get the pointer to the second coordinate array and the bounds
      ! of its global indices on the local DE.
      !-------------------------------------------------------------------
      call ESMF_GridGetCoord(new_grid, localDE             = 0, &
                                       staggerLoc          = staggerloc_, &
                                       coordDim            = 2, &
                                       farrayPtr           = coordY, &
                                       computationalLBound = tlb, &
                                       computationalUBound = tub, rc=rc)
      call check_error(rc)

      !-------------------------------------------------------------------
      ! Calculate and set coordinates in the second dimension.
      !-------------------------------------------------------------------
      do j = tlb(1), tub(1)
         coordY(j) = lat_points(j)    ! latitude  grid points
      enddo

      !write(LIS_logunit,*) '[INFO] Done creating the ESMF Rectilinear Grid: '// TRIM(grid_name)

   end function create_rectilinear_grid
!EOC
!------------------------------------------------------------------------------
!BOP
   subroutine simpleDynMaskProc(dynamicMaskList, dynamicSrcMaskValue, &
      dynamicDstMaskValue, rc)
      type(ESMF_DynamicMaskElementR4R8R4), pointer        :: dynamicMaskList(:)
      real(ESMF_KIND_R4),            intent(in), optional :: dynamicSrcMaskValue
      real(ESMF_KIND_R4),            intent(in), optional :: dynamicDstMaskValue
      integer,                       intent(out)          :: rc
      integer :: i, j
      real(ESMF_KIND_R8)  :: renorm
!EOP
!-------------------------------------------------------------------------
!BOC
      if (associated(dynamicMaskList)) then
         do i=1, size(dynamicMaskList)
            dynamicMaskList(i)%dstElement = 0.d0 ! set to zero
            renorm = 0.d0 ! reset
            do j=1, size(dynamicMaskList(i)%factor)
               if (dynamicSrcMaskValue /= dynamicMaskList(i)%srcElement(j)) then
                  dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement &
                         + dynamicMaskList(i)%factor(j) &
                         * dynamicMaskList(i)%srcElement(j)
                  renorm = renorm + dynamicMaskList(i)%factor(j)
               endif
            enddo
            if (renorm > 0.d0) then
               dynamicMaskList(i)%dstElement = dynamicMaskList(i)%dstElement / renorm
            else if (present(dynamicSrcMaskValue)) then
               dynamicMaskList(i)%dstElement = dynamicSrcMaskValue
            else
               rc = ESMF_RC_ARG_BAD  ! error detected
               return
            endif
         enddo
      endif
      ! return successfully
      rc = ESMF_SUCCESS
    end subroutine

!EOC
!-------------------------------------------------------------------------
!BOP

    subroutine runESMF_Regridding(srcField, dstField, routehandle, &
                                       dynamicMask, rc)

      type(ESMF_RouteHandle), intent(inOut) :: routehandle
      type(ESMF_Field),       intent(in) :: srcField
      type(ESMF_Field),       intent(inOut) :: dstField
      type(ESMF_DynamicMask), intent(inOut) :: dynamicMask
      integer,                intent(out)   :: rc

!EOP
!-------------------------------------------------------------------------
!BOC
       call ESMF_FieldRegrid(srcField      = srcField,        &
                             dstField      = dstField,        &
                             routehandle   = routehandle,     &
                             dynamicMask   = dynamicMask,     &
                             zeroregion    = ESMF_REGION_SELECT, &
                             termorderflag = ESMF_TERMORDER_SRCSEQ, &
                             checkflag   = .TRUE., rc=rc)
      call check_error(rc)
       !call LIS_verify(rc, 'ESMF_FieldRegrid failed')

    end subroutine runESMF_Regridding
!EOC
!-------------------------------------------------------------------------
!BOP
      subroutine decomposeDim(dim_world, dim_array, NDEs )
!
! !INPUT PARAMETERS:
      integer, intent(in)  :: dim_world ! total number of grid points
      integer, intent(in)  :: NDEs      ! number of DEs
!
! !OUTPUT PARAMETERS:
      integer, intent(out) :: dim_array(0:NDEs-1) 
!
! !DESCRIPTION:
! For a given number of grid points and a number of available processors,
! this subroutines determines the number of grid points assigned to each
! processor.
!
! !LOCAL VARIABLES:
      integer ::   n, im, rm
!EOP
!------------------------------------------------------------------------------
!BOC
      im = dim_world/NDEs
      rm = dim_world-NDEs*im
      do n=0,NDEs-1
                      dim_array(n) = im
      if( n.le.rm-1 ) dim_array(n) = im+1
      enddo
      end subroutine decomposeDim
!EOC
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: gaussian_comp_lats
! \label{gaussian_comp_lats}
!
! !INTERFACE:
subroutine gaussian_comp_lats(jmax, gaussian_lat_array)
! !USES:
  !use LIS_logMod , only : LIS_logunit

   implicit none
! !ARGUMENTS:
   integer, intent(in)                :: jmax
   real, dimension(jmax), intent(out) :: gaussian_lat_array

!
! !DESCRIPTION: 
!  This routine computes the latitudes of the global quasi-regular Gaussian 
!  grid corresponding to the total number of latitude circles given by jmax.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \end{description}
!EOP

   real*8, parameter :: pi=3.14159265358979

   real, allocatable, dimension(:) :: slat, wlat

   real*8  :: dpr
   integer :: j

   dpr =180./pi

   allocate(slat(jmax))
   allocate(wlat(jmax))

   call gausslat(jmax,slat,wlat)

   do j = 1, jmax
      gaussian_lat_array(j) = dpr*asin(slat(jmax-j+1))
   enddo

   deallocate(slat)
   deallocate(wlat)

   IF (LIS_masterproc) write(*,fmt=*) 'gaussian lats:', gaussian_lat_array
   !write(unit=LIS_logunit,fmt=*) 'gaussian lats:', gaussian_lat_array
end subroutine gaussian_comp_lats


!BOP
! !ROUTINE: gaussian_find_row
! \label{gaussian_find_row}
                                                
! !INTERFACE:
function gaussian_find_row(jmax, gaussian_lat_array, lat)
! !USES:
   !use LIS_logMod , only : LIS_logunit, LIS_endrun

   implicit none
! !ARGUMENTS:
   integer, intent(in)               :: jmax
   real, dimension(jmax), intent(in) :: gaussian_lat_array
   real, intent(in)                  :: lat
!
! !DESCRIPTION: 
!  This function computes the row number within the global quasi-regular 
!  Gaussian grid corresponding to the latitude given by lat.
!
!   \begin{description}
!   \item [jmax]
!      the total number of latitude circles
!   \item [gaussian\_lat\_array]
!      array of computed latitudes, in degrees
!   \item [lat]
!      latitude to search for
!   \end{description}
!EOP

   integer :: gaussian_find_row

   real    :: eps
   integer :: r

   eps = 180./(2*jmax)

   gaussian_find_row = -9999
   do r = 1, jmax
      if ( abs(gaussian_lat_array(r) - lat) < eps ) then
      !if ( gaussian_latitudes(1,r) == lat ) then
         gaussian_find_row = r
      endif
   enddo

   if ( gaussian_find_row == -9999 ) then
      STOP 9999
      !write(LIS_logunit,fmt=*) '[ERR] gaussian_find_row -- '// &
      !                         'Could not find lat',lat
      !call LIS_endrun
   endif

end function gaussian_find_row

!BOP  
!        
! !ROUTINE : gausslat
! \label{gausslat}
!
! !REVISION HISTORY:
!   04-16-92 Mark Iredell; Initial Specification
!   10-20-97 Mark Iredell; Increased precision
!   05-14-02 Urzula Jambor; Reduced limit of eps from e-12 to e-7
!     
! !INTERFACE:
subroutine gausslat(jmax,slat,wlat)
  implicit none
! !ARGUMENTS:
  integer       :: jmax
  real          :: slat(jmax)
  real          :: wlat(jmax)
! !DESCRIPTION:
!   This subroutine computes gaussian latitudes
!   Computes cosines of colatitude and gaussian weights
!   on the gaussian latitudes.  the gaussian latitudes are at
!   the zeroes of the legendre polynomial of the given order.
!
!  The arguments are:
!  \begin{description}
!    \item[jmax]
!     input number of latitudes
!    \item[slat]
!     cosines of colatitude
!    \item[wlat]
!     gaussian weights
!  \end{description}
!EOP
  real, parameter :: pi=3.14159265358979
  real, parameter :: eps=1.e-7
  integer, parameter :: jz=50
  real :: c
  integer:: jh, jhe, n, j
  real :: spmax, sp, r
  real :: pk(jmax/2),pkm1(jmax/2),pkm2(jmax/2)
  real :: bz(jz)
  data bz        / 2.4048255577,  5.5200781103, &
       8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679, &
       21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684, &
       33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132, &
       46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550, &
       58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299, &
       71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
       84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, &
       96.6052679510, 99.7468198587, 102.888374254, 106.029930916, &
       109.171489649, 112.313050280, 115.454612653, 118.596176630, &
       121.737742088, 124.879308913, 128.020877005, 131.162446275, &
       134.304016638, 137.445588020, 140.587160352, 143.728733573, &
       146.870307625, 150.011882457, 153.153458019, 156.295034268 /

  c=(1.-(2./pi)**2)*0.25
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  jh=jmax/2
  jhe=(jmax+1)/2
  r=1./sqrt((jmax+0.5)**2+c)
  do j=1,min(jh,jz)
     slat(j)=cos(bz(j)*r)
  enddo
  do j=jz+1,jh
     slat(j)=cos((bz(jz)+(j-jz)*pi)*r)
  enddo
  spmax=1.
  do while(spmax.gt.eps)
     spmax=0.
     do j=1,jh
        pkm1(j)=1.
        pk(j)=slat(j)
     enddo
     do n=2,jmax
        do j=1,jh
           pkm2(j)=pkm1(j)
           pkm1(j)=pk(j)
           pk(j)=((2*n-1)*slat(j)*pkm1(j)-(n-1)*pkm2(j))/n
        enddo
     enddo
     do j=1,jh
        sp=pk(j)*(1.-slat(j)**2)/(jmax*(pkm1(j)-slat(j)*pk(j)))
        slat(j)=slat(j)-sp
        spmax=max(spmax,abs(sp))
     enddo
  enddo
  do j=1,jh
     wlat(j)=(2.*(1.-slat(j)**2))/(jmax*pkm1(j))**2
     slat(jmax+1-j)=-slat(j)
     wlat(jmax+1-j)=wlat(j)
  enddo
  if(jhe.gt.jh) then
     slat(jhe)=0.
     wlat(jhe)=2./jmax**2
     do n=2,jmax,2
        wlat(jhe)=wlat(jhe)*n**2/(n-1)**2
     enddo
  endif
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return
end subroutine gausslat
!------------------------------------------------------------------------------
   subroutine check_error(userstatus)
       integer :: status, userstatus
       if (ESMF_LogFoundError(rcToCheck=userstatus, &
                              msg=ESMF_LOGERR_PASSTHRU, &
                              line=__LINE__, file=__FILE__, &
                              rcToReturn=status)) &
                              call ESMF_Finalize(endflag=ESMF_END_ABORT)
    end subroutine check_error


      end program test_gaussian_latlon
