! ******************************************************************************************************************************** !
! biogem_data_netCDF.f90
! BioGEochemical Model
! DATA LOADING/SAVING ROUTINES
! ******************************************************************************************************************************** !


MODULE biogem_data_netCDF


  USE gem_netcdf
  USE biogem_lib
  USE biogem_box
  IMPLICIT NONE
  SAVE


CONTAINS

  
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !
  ! netCDF INITIALIZATION ROUTINES
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! SAVE NETCDF RESTART DATA
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_data_netCDF_ncrstsave(dum_name,dum_yr,dum_iou)
    ! -------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    character(LEN=*),INTENT(IN)::dum_name                      !
    REAL,INTENT(in)::dum_yr                                    !
    INTEGER,INTENT(OUT)::dum_iou                               !
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    integer::io,is,l
    integer::loc_ntrec,loc_iou
    integer::loc_id_lonm,loc_id_latm,loc_id_lon_e,loc_id_lat_e
    integer::loc_id_zt,loc_id_zt_e
    integer,dimension(1:1)::loc_it_1
    integer,dimension(1:2)::loc_it_2
    integer,dimension(1:3)::loc_it_3
    character(127)::loc_title,loc_timunit
    character(7)::loc_string_year
    real::loc_c0,loc_c1
    real,dimension(0:n_i)::loc_lon_e
    real,dimension(0:n_j)::loc_lat_e
    real,dimension(0:n_k)::loc_zt_e
    real,dimension(n_i,n_j,n_k)::loc_ijk,loc_ijk_mask
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    loc_c0 = 0.0
    loc_c1 = 1.0
    ! -------------------------------------------------------- !
    ! WRITE TO FILE
    ! -------------------------------------------------------- !
    ! -------------------------------------------------------- ! open file
    call sub_opennew(dum_name,loc_iou)
    ! -------------------------------------------------------- ! start definitions
    call sub_redef(loc_iou)
    ! -------------------------------------------------------- ! set global attributes
    loc_string_year = fun_conv_num_char_n(8,int(dum_yr))
    loc_title = 'BIOGEM restart @ year '//loc_string_year
    call sub_putglobal(loc_iou,dum_name,loc_title,string_ncrunid,loc_timunit)
    ! -------------------------------------------------------- ! define dimensions
    call sub_defdim ('lon',loc_iou,n_i,loc_id_lonm)
    call sub_defdim ('lat',loc_iou,n_j,loc_id_latm)
    call sub_defdim ('lon_edges',loc_iou,n_i+1,loc_id_lon_e)
    call sub_defdim ('lat_edges',loc_iou,n_j+1,loc_id_lat_e)
    call sub_defdim ('zt',loc_iou, n_k,loc_id_zt)
    call sub_defdim ('zt_edges',loc_iou,n_k+1,loc_id_zt_e)
    ! -------------------------------------------------------- ! define 1d data (t)
    loc_it_1(1) = loc_id_lonm
    call sub_defvar ('lon',loc_iou,1,loc_it_1,loc_c0,loc_c0,'X','D','longitude of the t grid','longitude','degrees_east')
    loc_it_1(1) = loc_id_latm
    call sub_defvar ('lat',loc_iou,1,loc_it_1,loc_c0,loc_c0,'Y','D','latitude of the t grid','latitude','degrees_north')
    loc_it_1(1) = loc_id_lon_e
    call sub_defvar ('lon_edges',loc_iou,1,loc_it_1,loc_c0,loc_c0,' ','D','longitude of t grid edges',' ','degrees')
    loc_it_1(1) = loc_id_lat_e
    call sub_defvar ('lat_edges',loc_iou,1,loc_it_1,loc_c0,loc_c0,' ','D','latitude of t grid edges',' ','degrees')
    loc_it_1(1) = loc_id_zt
    call sub_defvar ('zt',loc_iou,1,loc_it_1,loc_c0,loc_c0,'Z','D','depth of z grid',' ','cm')
    loc_it_1(1) = loc_id_zt_e
    call sub_defvar ('zt_edges',loc_iou,1,loc_it_1,loc_c0,loc_c0,' ','D','depth of z grid edges',' ','m')
    loc_it_2(1) = loc_id_lonm
    loc_it_2(2) = loc_id_latm
    loc_it_3(1) = loc_id_lonm
    loc_it_3(2) = loc_id_latm
    loc_it_3(3) = loc_id_zt
    ! -------------------------------------------------------- ! define (3D) tracer variables -- dissolved
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       call sub_defvar('ocn_'//trim(string_ocn(io)),loc_iou,3,loc_it_3,loc_c0,loc_c0,' ','F', &
            & string_longname_ocn(io),'Ocean tracer - '//trim(string_ocn(io)),' ')
    end do
    ! -------------------------------------------------------- ! define (3D) tracer variables -- particulate
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       call sub_defvar('bio_part_'//trim(string_sed(is)),loc_iou,3,loc_it_3,loc_c0,loc_c0,' ','F', &
            & string_longname_sed(is),'Particulate tracer - '//trim(string_sed(is)),' ')
    end do
    ! -------------------------------------------------------- ! define (3D) tracer variables -- [H+]
    call sub_defvar('carb_'//trim(string_carb(ic_H)),loc_iou,3,loc_it_3,loc_c0,loc_c0,' ','F', &
         & string_longname_carb(ic_H),'Carbonate chemsitry - '//trim(string_carb(ic_H)),' ')
    ! -------------------------------------------------------- ! end definitions
    call sub_enddef (loc_iou)
    call sub_sync(loc_iou)
    ! -------------------------------------------------------- !
    loc_ntrec = 1
    ! -------------------------------------------------------- ! write 1D variables
    call sub_putvar1d ('lon',loc_iou,n_i,loc_ntrec,n_i,phys_ocn(ipo_lon,:,1,n_k),loc_c1,loc_c0)
    call edge_maker (1,loc_lon_e,phys_ocn(ipo_lon,:,1,n_k),phys_ocn(ipo_lone,:,1,n_k),phys_ocn(ipo_dlon,:,1,n_k),n_i)
    call sub_putvar1d ('lon_edges',loc_iou,n_i+1,loc_ntrec,n_i+1,loc_lon_e,loc_c1,loc_c0)
    call sub_putvar1d ('lat',loc_iou,n_j,loc_ntrec,n_j,phys_ocn(ipo_lat,1,:,n_k),loc_c1,loc_c0)
    call edge_maker (1,loc_lat_e,phys_ocn(ipo_lat,1,:,n_k),phys_ocn(ipo_latn,1,:,n_k),phys_ocn(ipo_dlat,1,:,n_k),n_j)
    call sub_putvar1d ('lat_edges',loc_iou, n_j+1,loc_ntrec,n_j+1,loc_lat_e,loc_c1,loc_c0)
    call sub_putvar1d ('zt',loc_iou,n_k,loc_ntrec,n_k,phys_ocn(ipo_Dmid,1,1,n_k:1:-1),loc_c1,loc_c0)
    call edge_maker (1,loc_zt_e,phys_ocn(ipo_Dmid,1,1,n_k:1:-1),phys_ocn(ipo_Dbot,1,1,n_k:1:-1), &
         & phys_ocn(ipo_dD,1,1,n_k:1:-1),n_k)
    loc_zt_e(0)=0.0
    call sub_putvar1d ('zt_edges',loc_iou,n_k+1,loc_ntrec,n_k+1,loc_zt_e,loc_c1,loc_c0)
    ! -------------------------------------------------------- ! write (3D) tracer variables -- dissolved
    loc_ijk_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ijk(:,:,:) = ocn(io,:,:,:)
       call sub_putvar3d('ocn_'//trim(string_ocn(io)),loc_iou,n_i,n_j,n_k,loc_ntrec, &
            & loc_ijk(:,:,n_k:1:-1),loc_ijk_mask(:,:,n_k:1:-1))
    end do
    ! -------------------------------------------------------- ! write (3D) tracer variables -- particulate
    loc_ijk_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ijk(:,:,:) = bio_part(is,:,:,:)
       call sub_putvar3d('bio_part_'//trim(string_sed(is)),loc_iou,n_i,n_j,n_k,loc_ntrec, &
            & loc_ijk(:,:,n_k:1:-1),loc_ijk_mask(:,:,n_k:1:-1))
    end do
    ! -------------------------------------------------------- ! write (3D) tracer variables -- [H+]
    loc_ijk_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    loc_ijk(:,:,:) = carb(ic_H,:,:,:)
    call sub_putvar3d('carb_'//trim(string_carb(ic_H)),loc_iou,n_i,n_j,n_k,loc_ntrec, &
         & loc_ijk(:,:,n_k:1:-1),loc_ijk_mask(:,:,n_k:1:-1))    
    ! -------------------------------------------------------- ! close file and return IOU
    call sub_closefile(loc_iou)
    dum_iou = loc_iou
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  END SUBROUTINE sub_data_netCDF_ncrstsave
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INITIALIZE netCDF
  SUBROUTINE sub_init_netcdf(dum_name,dum_iou,dum_dd)
    !-----------------------------------------------------------------------
    !       dummy arguments
    !-----------------------------------------------------------------------
    character(LEN=*),INTENT(IN) :: dum_name
    INTEGER,INTENT(OUT):: dum_iou
    INTEGER,INTENT(IN)::dum_dd
    !-----------------------------------------------------------------------
    !       define local variables
    !-----------------------------------------------------------------------
    character(255) :: loc_title, loc_timunit
    real           :: loc_c0, loc_c1
    integer        :: loc_it(6), loc_id_time, loc_id_lonm, loc_id_latp, loc_id_ztp
    integer        :: loc_id_lonps, loc_id_latps
    integer        :: loc_id_latm, loc_id_zt, loc_id_lon_e, loc_id_xu, loc_id_yu
    integer        :: loc_id_lat_e, loc_id_zt_e, loc_id_xu_e, loc_id_yu_e
    integer        :: loc_id_misc
    integer::loc_id_latp_e,loc_id_ztp_e
    integer::loc_id_lonps_e,loc_id_latps_e
    !-----------------------------------------------------------------------
    !       initialize local variables
    !-----------------------------------------------------------------------
    loc_c0 = 0.
    loc_c1 = 1.
    !-----------------------------------------------------------------------
    !       open file
    !-----------------------------------------------------------------------
    call sub_opennew (dum_name, dum_iou)
    !-----------------------------------------------------------------------
    !       start definitions
    !-----------------------------------------------------------------------
    call sub_redef (dum_iou)
    !-----------------------------------------------------------------------
    !       set global attributes
    !-----------------------------------------------------------------------
    loc_title = 'Time averaged integrals'
    write (loc_timunit,'(a)') 'Year mid-point'
    call sub_putglobal (dum_iou, dum_name, loc_title, string_ncrunid, loc_timunit)
    !-----------------------------------------------------------------------
    !       define dimensions
    !-----------------------------------------------------------------------
    call sub_defdim ('time', dum_iou, const_integer_zero, loc_id_time)
    call sub_defdim ('xu', dum_iou, n_i, loc_id_xu)
    call sub_defdim ('lon', dum_iou, n_i, loc_id_lonm)
    call sub_defdim ('lat', dum_iou, n_j, loc_id_latm)
    call sub_defdim ('zt', dum_iou, n_k, loc_id_zt)
    call sub_defdim ('yu', dum_iou, n_j, loc_id_yu)
    call sub_defdim ('lon_edges', dum_iou, n_i+1, loc_id_lon_e)
    call sub_defdim ('lat_edges', dum_iou, n_j+1, loc_id_lat_e)
    call sub_defdim ('zt_edges', dum_iou, n_k+1, loc_id_zt_e)
    call sub_defdim ('xu_edges', dum_iou, n_i+1, loc_id_xu_e)
    call sub_defdim ('yu_edges', dum_iou, n_j+1, loc_id_yu_e)
    call sub_defdim ('lat_moc', dum_iou, n_j+1, loc_id_latp)
    call sub_defdim ('zt_moc', dum_iou, n_k+1, loc_id_ztp)
    call sub_defdim ('lat_moc_edges', dum_iou, n_j+2, loc_id_latp_e)
    call sub_defdim ('zt_moc_edges', dum_iou, n_k+2, loc_id_ztp_e)
    call sub_defdim ('lon_psi', dum_iou, n_i, loc_id_lonps)
    call sub_defdim ('lat_psi', dum_iou, n_j+1, loc_id_latps)
    call sub_defdim ('lon_psi_edges', dum_iou, n_i+1, loc_id_lonps_e)
    call sub_defdim ('lat_psi_edges', dum_iou, n_j+2, loc_id_latps_e)
    call sub_defdim ('para', dum_iou, const_integer_one, loc_id_misc)
    !-----------------------------------------------------------------------
    !       define 1d data (t)
    !-----------------------------------------------------------------------
    loc_it(1) = loc_id_time
    call sub_defvar ('time', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'T', 'D' &
         &, 'Year', 'time', trim(loc_timunit))
    call sub_defvar ('year', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'T', 'F','year', ' ',' ')
    !-----------------------------------------------------------------------
    !       define 1d data (x, y or z)
    !-----------------------------------------------------------------------
    loc_it(1) = loc_id_lonm
    call sub_defvar ('lon', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of the t grid', 'longitude', 'degrees_east')
    loc_it(1) = loc_id_latm
    call sub_defvar ('lat', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of the t grid', 'latitude', 'degrees_north')
    loc_it(1) = loc_id_zt
    call sub_defvar ('zt', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Z', 'D' , &
         &'z-level mid depth', 'depth', 'm')
    loc_it(1) = loc_id_xu
    call sub_defvar ('xu', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of the u grid', 'longitude', 'degrees_east')
    loc_it(1) = loc_id_yu
    call sub_defvar ('yu', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of the u grid', 'latitude', 'degrees_north')
    loc_it(1) = loc_id_lon_e
    call sub_defvar ('lon_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of t grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_lat_e
    call sub_defvar ('lat_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of t grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_zt_e
    call sub_defvar ('zt_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Z', 'D' , &
         &'depth of t grid edges', ' ', 'm')
    loc_it(1) = loc_id_xu_e
    call sub_defvar ('xu_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
         &'longitude of u grid edges', ' ', 'degrees')
    loc_it(1) = loc_id_yu_e
    call sub_defvar ('yu_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
         &'latitude of u grid edges', ' ', 'degrees')
    SELECT CASE (dum_dd)
    CASE (2)
       ! MOC
       loc_it(1) = loc_id_latp
       call sub_defvar ('lat_moc', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
            &'latitude of moc grid', 'latitude', 'degrees_north')
       loc_it(1) = loc_id_ztp
       call sub_defvar ('zt_moc', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Z', 'D' , &
            &'depth of moc grid', 'depth', 'm')
       call sub_putatttext ('zt_moc', dum_iou, 'positive', 'down')
       loc_it(1) = loc_id_latp_e
       call sub_defvar ('lat_moc_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
            &'latitude of moc grid edges', 'latitude', 'degrees_north')
       loc_it(1) = loc_id_ztp_e
       call sub_defvar ('zt_moc_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Z', 'D' , &
            &'depth of moc grid edges', 'depth', 'm')
       ! PSI
       loc_it(1) = loc_id_lonps
       call sub_defvar ('lon_psi', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
            &'longitude of psi grid', 'longitude', 'degrees_east')
       loc_it(1) = loc_id_latps
       call sub_defvar ('lat_psi', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
            &'latitude of psi grid', 'latitude', 'degrees_north')
       loc_it(1) = loc_id_lonps_e
       call sub_defvar ('lon_psi_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'X', 'D' , &
            &'longitude of psi grid edges', 'longitude', 'degrees_east')
       loc_it(1) = loc_id_latps_e
       call sub_defvar ('lat_psi_edges', dum_iou, 1, loc_it(1), loc_c0, loc_c0, 'Y', 'D' , &
            &'latitude of psi grid edges', 'latitude', 'degrees_north')
    end select
    !-----------------------------------------------------------------------
    !       define basic 2d/3d data (x,y)
    !-----------------------------------------------------------------------
    SELECT CASE (dum_dd)
    CASE (2)
       loc_it(1) = loc_id_lonm
       loc_it(2) = loc_id_latm
       call sub_defvar('2Dgrid_level',dum_iou,2,loc_it(1:2),loc_c0,100.0,' ','I', &
            & 'grid definition','model_level_number','n/a')
       call sub_defvar ('2Dgrid_mask_landsea', dum_iou, 2, loc_it(1:2), loc_c0, 100.0, ' ', 'F', &
            &'land-sea mask', ' ' ,'n/a')
       call sub_defvar ('2Dgrid_topo', dum_iou, 2, loc_it(1:2), loc_c0, 100000., ' ', 'F', &
            &'ocean depth ', ' ' ,'m')
!!$       call sub_defvar ('axes_area', dum_iou, 2, loc_it(1:2), loc_c0, 0.5099044E+15, ' ', 'F', &
!!$            &'grid area ', ' ' ,'m2')
    end select
    SELECT CASE (dum_dd)
    CASE (3,4)
       loc_it(1) = loc_id_lonm
       loc_it(2) = loc_id_latm
       loc_it(3) = loc_id_zt
       call sub_defvar('2Dgrid_level',dum_iou,2,loc_it(1:2),loc_c0,100.0,' ','I', &
            & 'grid definition','model_level_number','n/a')
       call sub_defvar ('2Dgrid_mask_landsea', dum_iou, 2, loc_it(1:2), loc_c0, 100.0, ' ', 'F', &
            &'land-sea mask', ' ' ,'n/a')
       call sub_defvar ('2Dgrid_topo', dum_iou, 2, loc_it(1:2), loc_c0, 100000., ' ', 'F', &
            &'ocean depth ', ' ' ,'m')
       call sub_defvar('3Dgrid_mask_seafloor',dum_iou,3,loc_it(1:3),loc_c0, 1.,' ','F', &
            & 'ocean mask',' ','n/a')
!!$       call sub_defvar ('axes_area', dum_iou, 2, loc_it(1:2), loc_c0, 0.5099044E+15, ' ', 'F', &
!!$            &'grid area ', ' ' ,'m2')
    end select
    !-----------------------------------------------------------------------
    call sub_enddef (dum_iou)
    call sub_closefile (dum_iou)
    !-----------------------------------------------------------------------
  END SUBROUTINE sub_init_netcdf
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_save_netcdf(dum_yr,dum_dd)
    !-----------------------------------------------------------------------
    !       dummy arguments
    !-----------------------------------------------------------------------
    INTEGER,INTENT(IN)::dum_dd
    REAL,INTENT(in):: dum_yr
    !-----------------------------------------------------------------------
    !       local variables
    !-----------------------------------------------------------------------
    character(255) :: loc_name
    real           :: loc_c0, loc_c1
    integer        :: i, j, k, loc_i, loc_iou, loc_ntrec
    real,dimension(n_i,n_j) :: loc_mask_ALL,loc_data
    real,dimension(n_i+1) :: loc_lon_e, loc_xu_e
    real,dimension(n_j+1) :: loc_lat_e, loc_yu_e
    real,dimension(0:n_k) :: loc_zt_e, loc_help
    REAL,DIMENSION(0:n_k+1)::loc_grid_dz, loc_tmp_k
    real,dimension(0:n_j+1)::loc_lat_moc_e
    REAL,DIMENSION(0:n_k+1)::loc_zt_moc_e
    logical :: loc_defined
    !-----------------------------------------------------------------------
    !       initialize local variables
    !-----------------------------------------------------------------------
    loc_c0 = 0.
    loc_c1 = 1.
    !
    loc_mask_ALL(:,:)  = 1.0
    loc_data(:,:)      = 0.0
    loc_lon_e(:)       = 0.0
    loc_xu_e(:)        = 0.0
    loc_lat_e(:)       = 0.0
    loc_yu_e(:)        = 0.0
    loc_zt_e(:)        = 0.0
    loc_help(:)        = 0.0
    loc_grid_dz(:)     = 0.0
    loc_tmp_k(:)       = 0.0
    loc_lat_moc_e(:)   = 0.0
    loc_zt_moc_e(:)    = 0.0
    ! set local netCDF variables
    SELECT CASE (dum_dd)
    CASE (-1)
       loc_name = string_ncrst
       loc_iou = ncrst_iou
       loc_ntrec = ncrst_ntrec
    CASE (2)
       loc_name = string_ncout2d
       loc_iou = ncout2d_iou
       loc_ntrec = ncout2d_ntrec
    CASE (3)
       loc_name = string_ncout3d
       loc_iou = ncout3d_iou
       loc_ntrec = ncout3d_ntrec
    CASE (4)
       loc_name = string_ncout3dsig
       loc_iou = ncout3dsig_iou
       loc_ntrec = ncout3dsig_ntrec
    CASE DEFAULT
       CALL sub_report_error( &
            & 'biogem_data_netCDF','sub_save_netcdf', &
            & 'illegal netCDF dimension', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    end select
    ! open file and get latest record number
    loc_defined = .true.
    loc_i = 0
    if (loc_ntrec .eq. 0) then
       loc_defined = .false.
       loc_i = 1
    end if
    call sub_opennext (loc_name, dum_yr, loc_i, loc_ntrec, loc_iou)
    ! write time data
    call sub_putvars  ('time', loc_iou, loc_ntrec, dum_yr, loc_c1, loc_c0)
    call sub_putvarIs  ('year', loc_iou, loc_ntrec, nint(dum_yr), loc_c1, loc_c0)
    if(.not. loc_defined) then
       ! write 1d data: x
       call sub_putvar1d ('lon', loc_iou, n_i, loc_ntrec, n_i, &
            & phys_ocn(ipo_lon,:,1,1), loc_c1, loc_c0)
       call edge_maker (1, loc_lon_e, phys_ocn(ipo_lon,:,1,1), &
            & phys_ocn(ipo_lone,:,1,1), phys_ocn(ipo_dlon,:,1,1), n_i)
       call sub_putvar1d ('lon_edges', loc_iou, n_i+1, loc_ntrec, n_i+1, &
            & loc_lon_e, loc_c1, loc_c0)
       call sub_putvar1d ('xu', loc_iou, n_i, loc_ntrec, n_i, &
            & loc_lon_e(1:n_i), loc_c1, loc_c0)
       call edge_maker (2, loc_xu_e, phys_ocn(ipo_lon,:,1,1), &
            & phys_ocn(ipo_lone,:,1,1), phys_ocn(ipo_dlon,:,1,1), n_i)
       call sub_putvar1d ('xu_edges', loc_iou, n_i+1, loc_ntrec, n_i+1, &
            & loc_xu_e, loc_c1, loc_c0)
       ! write 1d data: y
       call sub_putvar1d ('lat', loc_iou, n_j, loc_ntrec, n_j, &
            & phys_ocn(ipo_lat,1,:,1), loc_c1, loc_c0)
       call edge_maker (1, loc_lat_e, phys_ocn(ipo_lat,1,:,1), &
            & phys_ocn(ipo_latn,1,:,1), phys_ocn(ipo_dlat,1,:,1), n_j)
       call sub_putvar1d ('lat_edges', loc_iou, n_j+1, loc_ntrec, n_j+1, &
            & loc_lat_e, loc_c1, loc_c0)
       call sub_putvar1d ('yu', loc_iou, n_j, loc_ntrec, n_j, &
            & loc_lat_e(1:n_j), loc_c1, loc_c0)
       call edge_maker (2, loc_yu_e, phys_ocn(ipo_lat,1,:,1), &
            & phys_ocn(ipo_latn,1,:,1), phys_ocn(ipo_dlat,1,:,1), n_j)
       call sub_putvar1d ('yu_edges', loc_iou, n_j+1, loc_ntrec, n_j+1, &
            & loc_yu_e, loc_c1, loc_c0)
       ! write 1d data: z
       call sub_putvar1d ('zt', loc_iou, n_k, loc_ntrec, n_k, &
            & phys_ocn(ipo_Dmid,1,1,n_k:1:-1), loc_c1, loc_c0)
       loc_help(1:n_k) = phys_ocn(ipo_dD,1,1,n_k:1:-1)
       call edge_maker (1, loc_zt_e, phys_ocn(ipo_Dmid,1,1,n_k:1:-1), &
            & phys_ocn(ipo_Dbot,1,1,n_k:1:-1), loc_help , n_k)
       loc_zt_e(0)=0.0
       call sub_putvar1d ('zt_edges', loc_iou, n_k+1, loc_ntrec, n_k+1, &
            & loc_zt_e, loc_c1, loc_c0)
       !
       SELECT CASE (dum_dd)
       CASE (2)
          ! MOC
          call sub_putvar1d ('lat_moc', loc_iou, n_j+1, loc_ntrec, n_j+1, &
               & (180.0/const_pi) * ASIN(goldstein_sv(:)), loc_c1, loc_c0)
          loc_grid_dz(1:n_k) = goldstein_dz(:)
          DO k=n_k,0,-1
             loc_tmp_k(n_k-k) = SUM(goldstein_dsc * loc_grid_dz(k+1:n_k+1))
          ENDDO
          call sub_putvar1d ('zt_moc', loc_iou, n_k+1, loc_ntrec, n_k+1, &
               & loc_tmp_k, loc_c1, loc_c0)
          loc_lat_moc_e(0) = phys_ocn(ipo_lat,1,1,1) - (phys_ocn(ipo_latn,1,1,1) + 90.0)
          loc_lat_moc_e(1:n_j) = phys_ocn(ipo_lat,1,1:n_j,1)
          loc_lat_moc_e(n_j+1) = phys_ocn(ipo_lat,1,n_j,1) + (phys_ocn(ipo_latn,1,n_j,1) - phys_ocn(ipo_latn,1,n_j-1,1))
          call sub_putvar1d ('lat_moc_edges', loc_iou, n_j+2, loc_ntrec, n_j+2, &
               & loc_lat_moc_e(:), loc_c1, loc_c0)
          loc_zt_moc_e(0) = 0.0
          loc_zt_moc_e(1:n_k) = phys_ocn(ipo_Dmid,1,1,n_k:1:-1)
          loc_zt_moc_e(n_k+1) = phys_ocn(ipo_Dbot,1,1,1)
          call sub_putvar1d ('zt_moc_edges', loc_iou, n_k+2, loc_ntrec, n_k+2, &
               & loc_zt_moc_e(:), loc_c1, loc_c0)
          ! PSI
          call sub_putvar1d('lon_psi',loc_iou,n_i,loc_ntrec,n_i,loc_lon_e(2:n_i+1),loc_c1,loc_c0)
          call sub_putvar1d('lat_psi',loc_iou,n_j+1,loc_ntrec,n_j+1,loc_lat_e(1:n_j+1),loc_c1,loc_c0)
          call sub_putvar1d('lon_psi_edges',loc_iou,n_i+1,loc_ntrec,n_i+1,loc_xu_e(:),loc_c1,loc_c0)
          call sub_putvar1d('lat_psi_edges',loc_iou,n_j+2,loc_ntrec,n_j+2,loc_lat_moc_e(:),loc_c1,loc_c0)
       end select
       ! set maximum ocean depth
       do i=1,n_i
          do j=1,n_j
             if(phys_ocn(ipo_mask_ocn,i,j,n_k) == 1.0) then
                loc_data(i,j) = phys_ocn(ipo_Dbot,i,j,goldstein_k1(i,j))
             else
                loc_data(i,j) = 0.0
             end if
          end do
       end do
       ! write 2D grid data
       SELECT CASE (dum_dd)
       CASE (2)
          call sub_putvar2dI ('2Dgrid_level', loc_iou, n_i, n_j, loc_ntrec, &
               & goldstein_k1(:,:))
          call sub_putvar2d ('2Dgrid_mask_landsea', loc_iou, n_i, n_j, loc_ntrec, &
               & phys_ocn(ipo_mask_ocn,:,:,n_k), phys_ocn(ipo_mask_ocn,:,:,n_k))
          call sub_putvar2d ('2Dgrid_topo', loc_iou, n_i, n_j, loc_ntrec, &
               & loc_data(:,:), phys_ocn(ipo_mask_ocn,:,:,n_k))
!!$          call sub_putvar2d ('axes_area', loc_iou, n_i, n_j, loc_ntrec, &
!!$               & phys_ocnatm(ipoa_A,:,:), loc_mask_ALL)
       end select
       ! write 3D grid data
       SELECT CASE (dum_dd)
       CASE (3,4)
          call sub_putvar2dI ('2Dgrid_level', loc_iou, n_i, n_j, loc_ntrec, &
               & goldstein_k1(:,:))
          call sub_putvar2d ('2Dgrid_mask_landsea', loc_iou, n_i, n_j, loc_ntrec, &
               & phys_ocn(ipo_mask_ocn,:,:,n_k), phys_ocn(ipo_mask_ocn,:,:,n_k))
          call sub_putvar2d ('2Dgrid_topo', loc_iou, n_i, n_j, loc_ntrec, &
               & loc_data(:,:), phys_ocn(ipo_mask_ocn,:,:,n_k))
          call sub_putvar3d('3Dgrid_mask_seafloor',loc_iou,n_i,n_j,n_k,loc_ntrec, &
               & phys_ocn(ipo_mask_ocn,:,:,n_k:1:-1),phys_ocn(ipo_mask_ocn,:,:,n_k:1:-1))
!!$          call sub_putvar2d ('axes_area', loc_iou, n_i, n_j, loc_ntrec, &
!!$               & phys_ocn(ipo_A,:,:,n_k), phys_ocn(ipo_mask_ocn,:,:,n_k))
       end select
    end if
    !-----------------------------------------------------------------------
    ! update record number
    SELECT CASE (dum_dd)
    CASE (-1)
       ncrst_ntrec = loc_ntrec
    CASE (2)
       ncout2d_ntrec = loc_ntrec
    CASE (3)
       ncout3d_ntrec = loc_ntrec
    CASE (4)
       ncout3dsig_ntrec = loc_ntrec
    end select
    !
    call sub_sync(loc_iou)
    !-----------------------------------------------------------------------
  END SUBROUTINE sub_save_netcdf
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_init_netcdf_TM(dum_wet_grid)
    !-----------------------------------------------------------------------
    !       dummy arguments
    !-----------------------------------------------------------------------
    integer::dum_wet_grid
    !-----------------------------------------------------------------------
    !       define local variables
    !-----------------------------------------------------------------------
    integer::ncid,status
    integer::coo_dimid,index_dimid,var_id

    ! open file
    status=nf90_create(TRIM(par_outdir_name)//'transport_matrix_COO.nc',nf90_clobber,ncid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    ! dimensions
    status=nf90_def_dim(ncid,'coo',nf90_unlimited,coo_dimid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_dim(ncid,'index',dum_wet_grid,index_dimid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    ! define variables
    status=nf90_def_var(ncid,'coo_val',nf90_double,(/coo_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'coo_row',nf90_int,(/coo_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'coo_col',nf90_int,(/coo_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'coo_avg_n',nf90_byte,(/coo_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'index_i',nf90_int,(/index_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'index_j',nf90_int,(/index_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    status=nf90_def_var(ncid,'index_k',nf90_int,(/index_dimid/),var_id)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    !end definition
    status=nf90_enddef(ncid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    ! close file
    status=nf90_close(ncid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

  END SUBROUTINE sub_init_netcdf_TM
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_save_netcdf_TM(dum_TM_flag,dum_start,dum_val,dum_col,dum_row,dum_avg_n,dum_i,dum_j,dum_k)

    integer::dum_start
    real::dum_val
    integer::dum_TM_flag
    integer::dum_col,dum_row,dum_avg_n,dum_i,dum_j,dum_k

    integer::ncid,status,nc_record_count,dimid
    character(len=100)::name

    ! open file
    status=nf90_open(TRIM(par_outdir_name)//'transport_matrix_COO.nc',nf90_write,ncid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

    select case (dum_TM_flag)
    case(1) ! write TM data

       ! coo is unlimited dimension that we have to append to
       ! find the current length of coo:
       status=nf90_inq_dimid(ncid, 'coo', dimid)
       if(status /= nf90_NoErr) print*,trim(nf90_strerror(status)//',coo')

       status=nf90_inquire_dimension(ncid, dimid, name, nc_record_count)
       if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

       nc_record_count=nc_record_count+1 ! start count for appending data

       call sub_putvars ('coo_val', ncid, nc_record_count, dum_val,1.0,0.0)

       call sub_putvarIs ('coo_col', ncid, nc_record_count, dum_col,1.0,0.0)

       call sub_putvarIs ('coo_row', ncid, nc_record_count, dum_row,1.0,0.0)

       call sub_putvarIs ('coo_avg_n', ncid, nc_record_count, dum_avg_n,1.0,0.0)

    case(0) ! write index data

       call sub_putvarIs ('index_i', ncid, dum_start, dum_i,1.0,0.0)

       call sub_putvarIs ('index_j', ncid, dum_start, dum_j,1.0,0.0)

       call sub_putvarIs ('index_k', ncid, dum_start, dum_k,1.0,0.0)
    end select

    ! close file
    status=nf90_close(ncid)
    if(status /= nf90_NoErr) print*,trim(nf90_strerror(status))

  END SUBROUTINE sub_save_netcdf_TM
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! TIME-SERIES 3-D FIELDS
  SUBROUTINE sub_save_netcdf_3d_sig()
    !-----------------------------------------------------------------------
    !       DEFINE LOCAL VARIABLES
    !-----------------------------------------------------------------------
    INTEGER::l,i,j,k,loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    real::loc_tot,loc_frac,loc_standard
    !-----------------------------------------------------------------------
    !       INITIALIZE LOCAL VARIABLES
    !-----------------------------------------------------------------------
    loc_iou   = ncout3dsig_iou
    loc_ntrec = ncout3dsig_ntrec
    loc_mask  = phys_ocn(ipo_mask_ocn,:,:,:)
    !----------------------------------------------------------------
    !       SAVE OCEAN TRACER FIELD
    !----------------------------------------------------------------
    DO l=1,n_l_ocn
       loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                SELECT CASE (ocn_type(l2io(l)))
                CASE (0)
                   if (l == io2l(io_T)) then
                      loc_ijk(i,j,k) = int_misc_3D_sig(l,i,j,k)/int_t_sig - const_zeroC
                   else
                      loc_ijk(i,j,k) = int_misc_3D_sig(l,i,j,k)/int_t_sig
                   end if
                CASE (1)
                   loc_ijk(i,j,k) = int_misc_3D_sig(l,i,j,k)/int_t_sig
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_misc_3D_sig(io2l(ocn_dep(l2io(l))),i,j,k)/int_t_sig
                   loc_frac = int_misc_3D_sig(l,i,j,k)/int_t_sig
                   loc_standard = const_standards(ocn_type(l2io(l)))
                   loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                END SELECT
             end do
          end do
       end do
       SELECT CASE (ocn_type(l2io(l)))
       CASE (0,1,n_itype_min:n_itype_max)
          call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn_tname(l)), &
               & trim(string_ocn_tlname(l)),trim(string_ocn_unit(l)),ocn_mima(l,1),ocn_mima(l,2))
          call sub_putvar3d_g('ocn_'//trim(string_ocn(l2io(l))),loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       END SELECT
    END DO
    !----------------------------------------------------------------
  END SUBROUTINE sub_save_netcdf_3d_sig
  ! ****************************************************************************************************************************** !
  
  
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !
  ! 3-D FIELDS
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! *** save time-slice data ***
  SUBROUTINE sub_save_netcdf_3d(dum_t)
    ! dummy arguments
    REAL,INTENT(in)::dum_t
    ! ---------------------------------------------------------------- !
    ! reservoir FIELDS
    ! ---------------------------------------------------------------- !
    if (ctrl_save_basic_reservoirs) CALL sub_3d_save_reservoirs_basic()
    if (ctrl_save_advanced_reservoirs) CALL sub_3d_save_reservoirs_advanced()
    ! ---------------------------------------------------------------- !
    ! geochemsitry FIELDS
    ! ---------------------------------------------------------------- !
    if (ctrl_save_basic_geochemistry) CALL sub_3d_save_geochemistry_basic()
    if (ctrl_save_advanced_geochemistry) CALL sub_3d_save_geochemistry_advanced()
    ! ---------------------------------------------------------------- !
    ! biological pump FIELDS
    ! ---------------------------------------------------------------- !
    if (ctrl_save_basic_biologicalpump) CALL sub_3d_save_biologicalpump_basic()
    if (ctrl_save_advanced_biologicalpump) CALL sub_3d_save_biologicalpump_advanced()
    ! ---------------------------------------------------------------- !
    ! proxy-related FIELDS
    ! ---------------------------------------------------------------- !
    if (ctrl_save_basic_proxies) CALL sub_3d_save_proxies_basic()
    if (ctrl_save_advanced_proxies) CALL sub_3d_save_proxies_advanced()
    !----------------------------------------------------------------
    ! HIDDEN
    ! ---------------------------------------------------------------- ! grid
    if (ctrl_save_hidden_grid) call sub_3d_save_hidden_grid
    ! ---------------------------------------------------------------- ! climate
    if (ctrl_save_hidden_climate) call sub_3d_save_hidden_climate
    ! ---------------------------------------------------------------- ! preformed tracer related FIELDS
    if (ctrl_save_hidden_preformedtracers) CALL sub_3d_save_hidden_preformedtracers()
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_save_netcdf_3d
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_reservoirs_basic()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    INTEGER::l,i,j,k,io,is,ip,ic,icc,loc_iou,loc_ntrec
    integer::id
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask,loc_sed_mask
    real::loc_tot,loc_frac,loc_standard
    real::loc_d13C,loc_d14C
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask        = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! OCEAN TRACER PROPERTIES
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                SELECT CASE (ocn_type(io))
                CASE (0)
                   if (io == io_T) then
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice - const_zeroC
                   else
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                   end if
                CASE (1)
                   loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,k)/int_t_timeslice
                   loc_frac = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                   loc_standard = const_standards(ocn_type(io))
                   loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                END SELECT
             end do
          end do
       end do
       SELECT CASE (ocn_type(io))
       CASE (0)
          ! T,S
          call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn_tname(l)), &
               & trim(string_ocn_tlname(l)),trim(string_ocn_unit(l)),ocn_mima(l,1),ocn_mima(l,2))
          call sub_putvar3d_g('ocn_'//trim(string_ocn(io)),loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       CASE (1)
          ! bulk tracers
          call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn_tname(l)), &
               & trim(string_ocn_tlname(l)),trim(string_ocn_unit(l)),ocn_mima(l,1),ocn_mima(l,2))
          call sub_putvar3d_g('ocn_'//trim(string_ocn(io)),loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       CASE (n_itype_min:n_itype_max)
          ! isotopes
          If (ctrl_save_basic_proxies) then
             call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn_tname(l)), &
                  & trim(string_ocn_tlname(l)),trim(string_ocn_unit(l)),ocn_mima(l,1),ocn_mima(l,2))
             call sub_putvar3d_g('ocn_'//trim(string_ocn(io)),loc_iou,n_i,n_j,n_k, &
                  & loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end if
       END SELECT
    END DO
    ! radiocarbon in DELTA notation
    IF (ocn_select(io_DIC_13C) .AND. ocn_select(io_DIC_14C)) THEN
       loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice
                loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,k)/int_t_timeslice
                loc_standard = const_standards(ocn_type(io_DIC_13C))
                loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_frac = int_ocn_timeslice(io_DIC_14C,i,j,k)/int_t_timeslice
                loc_standard = const_standards(ocn_type(io_DIC_14C))
                loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_ijk(i,j,k) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
             end do
          end do
       end do
       loc_unitsname = 'o/oo'
       If (ctrl_save_basic_proxies) then
          call sub_adddef_netcdf(loc_iou,4,'ocn_DIC_D14C', &
               & ' oceanic D14C (big delta)',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('ocn_DIC_D14C',loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end if
    end if
    ! radiocarbon AGE
    ! NOTE: assuming the values already in loc_ijk (above)
    !       BUT need to (re)calculate atmospheric D14C ...
    IF (ocn_select(io_DIC_13C) .AND. ocn_select(io_DIC_14C)) THEN
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             ! first calculate D14C for the atmopshere
             loc_tot  = int_sfcatm1_timeslice(ia_pCO2,i,j)/int_t_timeslice
             loc_frac = int_sfcatm1_timeslice(ia_pCO2_13C,i,j)/int_t_timeslice
             loc_standard = const_standards(atm_type(ia_pCO2_13C))
             loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
             loc_frac = int_sfcatm1_timeslice(ia_pCO2_14C,i,j)/int_t_timeslice
             loc_standard = const_standards(atm_type(ia_pCO2_14C))
             loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
             loc_ij(i,j) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
             ! now convert to radiocarbon age
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = fun_convert_D14Ctoage(loc_ijk(i,j,k),loc_ij(i,j))
             end do
          end do
       end do
       loc_unitsname = 'years'
       If (ctrl_save_basic_proxies) then
          call sub_adddef_netcdf(loc_iou,4,'ocn_DIC_D14C_age', &
               & ' oceanic D14C age',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('ocn_DIC_D14C_age',loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end if
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_reservoirs_basic
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_reservoirs_advanced()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,is,l,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    real::loc_ocn_mean_S
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! SALINITY-NORMALIZED OCEAN TRACER FIELD
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_extra) then
       loc_ocn_mean_S = SUM(int_ocn_timeslice(io_S,:,:,:)*phys_ocn(ipo_M,:,:,:))/SUM(phys_ocn(ipo_M,:,:,:))
       DO l=3,n_l_ocn
          io = conv_iselected_io(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   SELECT CASE (ocn_type(io))
                   CASE (0,1)
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)* &
                           & (loc_ocn_mean_S/int_ocn_timeslice(io_S,i,j,k))/int_t_timeslice
                      loc_unitsname = 'mol kg-1'
                   END SELECT
                end DO
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (0,1)
             call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn(io))//'_Snorm', &
                  & trim(string_ocn(io))//' normalized by salinity',trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar3d_g('ocn_'//trim(string_ocn(io))//'_Snorm',loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          END SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! SAVE TRACER INVENTORIES
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_extra) then
       DO l=3,n_l_ocn
          io = conv_iselected_io(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   SELECT CASE (ocn_type(io))
                   CASE (1,n_itype_min:n_itype_max)
                      loc_ijk(i,j,k) = phys_ocn(ipo_M,i,j,k)*int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                      loc_unitsname = 'mol'
                   END SELECT
                end DO
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (1)
             ! bulk tracers
             If (ctrl_save_basic_reservoirs) then
                call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn(io))//'_tot', &
                     & trim(string_ocn(io))//' volume integrated inventory',trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar3d_g('ocn_'//trim(string_ocn(io))//'_tot',loc_iou, &
                     & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
             end if
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou,4,'ocn_'//trim(string_ocn(io))//'_tot', &
                     & trim(string_ocn(io))//' volume integrated inventory',trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar3d_g('ocn_'//trim(string_ocn(io))//'_tot',loc_iou, &
                     & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
             end if
          END SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_reservoirs_advanced
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_geochemistry_basic
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,is,l,id,ic,icc
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! CARBONATE CHEMISTRY FIELD
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       DO ic=1,n_carb
          loc_ijk(:,:,:) = int_carb_timeslice(ic,:,:,:)/int_t_timeslice
          SELECT CASE (ic)
          CASE (ic_pHsws)
             loc_unitsname = 'pH(sws)'                   
          CASE (ic_ohm_cal,ic_ohm_arg,ic_RF0)     
             loc_unitsname = 'n/a' 
          case default
             loc_unitsname = 'mol kg-1'
          end SELECT
          SELECT CASE (ic)
          CASE (ic_pHsws,ic_conc_CO2,ic_conc_CO3,ic_conc_HCO3,ic_ohm_cal,ic_ohm_arg)
             call sub_adddef_netcdf(loc_iou,4,'carbchem_'//trim(string_carb(ic)), &
                  & 'carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar3d_g('carbchem_'//trim(string_carb(ic)),loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! GEOCHEMICAL DIAGNOSTICS -- PRECIP
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol kg-1 yr-1'
    DO id=1,n_diag_precip
       loc_ijk(:,:,:) = int_diag_precip_timeslice(id,:,:,:)/int_t_timeslice
       if (ctrl_save_n_diag_precip(id)) then
          call sub_adddef_netcdf(loc_iou,4,'react_'//trim(string_diag_precip(id)), &
               & 'production rate - '//trim(string_diag_precip(id)),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('react_'//trim(string_diag_precip(id)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end if
    end DO
    ! ---------------------------------------------------------------- !
    ! GEOCHEMICAL DIAGNOSTICS -- REACTION
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol kg-1 yr-1'
    DO id=1,n_diag_react
       loc_ijk(:,:,:) = int_diag_react_timeslice(id,:,:,:)/int_t_timeslice
       if (ctrl_save_n_diag_react(id)) then
          call sub_adddef_netcdf(loc_iou,4,'react_'//trim(string_diag_react(id)), &
               & 'reaction rate - '//trim(string_diag_react(id)),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('react_'//trim(string_diag_react(id)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end if
    end DO
    ! ---------------------------------------------------------------- !
    ! GEOCHEMICAL DIAGNOSTICS -- REDOX -- aqueous reactions
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_redox) then
       loc_unitsname = 'mol kg-1 yr-1'
       DO id=1,n_diag_redox_aq
          loc_ijk(:,:,:) = int_diag_redox_timeslice(id,:,:,:)/int_t_timeslice
          call sub_adddef_netcdf(loc_iou,4,'react_'//trim(string_diag_redox(id)), &
               & 'redox transformation rate - '//trim(string_diag_redox(id)),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('react_'//trim(string_diag_redox(id)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end DO
    end If
    !-----------------------------------------------------------------------
    !       Fe speciation
    !-----------------------------------------------------------------------
    IF ( ocn_select(io_Fe) .AND. ocn_select(io_FeL) ) THEN
       loc_unitsname = 'nmol kg-1'
       loc_ijk(:,:,:) = 1.0E9*(int_ocn_timeslice(io_Fe,:,:,:) + int_ocn_timeslice(io_FeL,:,:,:))/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,4,'geochem_Fe_FeT','Total dissolved iron concentration', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('geochem_Fe_FeT',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       loc_ijk(:,:,:) = 1.0E9*(int_ocn_timeslice(io_L,:,:,:) + int_ocn_timeslice(io_FeL,:,:,:))/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,4,'geochem_Fe_LT','Total Fe-binding ligand concentration', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('geochem_Fe_LT',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    !----------------------------------------------------------------
    !       Fe SPECIATION [complete]
    !----------------------------------------------------------------
    If ( ocn_select(io_Fe) .OR. ocn_select(io_TDFe) ) then
       DO id=1,n_diag_iron
          select case (id)
          CASE (idiag_iron_TFe3pct)
             loc_unitsname = 'percent (%)'
          case default
             loc_unitsname = 'mol kg-1'
          end SELECT
          loc_ijk(:,:,:) = int_diag_iron_timeslice(id,:,:,:)/int_t_timeslice
          call sub_adddef_netcdf(loc_iou,4,'geochem_Fe_'//trim(string_diag_iron(id)), &
               & 'water-column Fe speciation - '//trim(string_diag_iron(id)), &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('geochem_Fe_'//trim(string_diag_iron(id)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end DO
       !----------------------------------------------------------------
       !       Fe speciation [isotopes]
       !----------------------------------------------------------------
       If (ctrl_save_basic_proxies) then
          loc_ijk(:,:,:) = const_real_null
          loc_standard   = const_standards(ocn_type(io_Fe_56Fe))
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   if (ocn_select(io_Fe_56Fe) .AND. ocn_select(io_FeL_56Fe)) then
                      loc_tot  = (int_ocn_timeslice(io_Fe,i,j,k) + int_ocn_timeslice(io_FeL,i,j,k))/int_t_timeslice
                      loc_frac = (int_ocn_timeslice(io_Fe_56Fe,i,j,k) + int_ocn_timeslice(io_FeL_56Fe,i,j,k))/int_t_timeslice
                      loc_save = .true.
                   elseif (ocn_select(io_Fe_56Fe) .AND. ocn_select(io_Fe2_56Fe)) then
                      loc_tot  = (int_ocn_timeslice(io_Fe,i,j,k) + int_ocn_timeslice(io_Fe2,i,j,k))/int_t_timeslice
                      loc_frac = (int_ocn_timeslice(io_Fe_56Fe,i,j,k) + int_ocn_timeslice(io_Fe2_56Fe,i,j,k))/int_t_timeslice
                      loc_save = .true.
                   elseif (ocn_select(io_TDFe_56Fe)) then
                      loc_tot  = int_ocn_timeslice(io_TDFe,i,j,k)/int_t_timeslice
                      loc_frac = int_ocn_timeslice(io_TDFe_56Fe,i,j,k)/int_t_timeslice
                      loc_save = .true.
                   else
                      loc_save = .false.
                   end if
                   if (loc_save) loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                END DO
             END DO
          END DO
          if (loc_save) then
             loc_unitsname  = 'o/oo'
             loc_min        = -1000.0
             loc_max        = +1000.0
             loc_name = 'geochem_Fe_TDFe_d56Fe'
             call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'water-column Fe speciation - d56Fe (total dissolved iron)', &
                  & trim(loc_unitsname),loc_min,loc_max)
             call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end if
       end If
    end if
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_geochemistry_basic
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_geochemistry_advanced
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,is,l,id,ic,icc
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! CARBONATE CHEMISTRY FIELD
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       DO ic=1,n_carb
          loc_ijk(:,:,:) = int_carb_timeslice(ic,:,:,:)/int_t_timeslice
          SELECT CASE (ic)
          CASE (ic_pHsws)
             loc_unitsname = 'pH(sws)'           
          CASE (ic_H,ic_pH_n,ic_err)  
             loc_unitsname = 'n/a' 
          case default
             loc_unitsname = 'mol kg-1'
          end SELECT
          SELECT CASE (ic)
          CASE (ic_H,ic_pH_n,ic_err,ic_pHsws)
             call sub_adddef_netcdf(loc_iou,4,'carbchem_'//trim(string_carb(ic)), &
                  & 'carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar3d_g('carbchem_'//trim(string_carb(ic)),loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end SELECT
       END DO
       ! save accumulated carbchem error occurrence array
       loc_ijk(:,:,:) = diag_carb_errsum(:,:,:)
       loc_unitsname = 'n/a'
       loc_shortname = 'carbchem_errsum'
       loc_longname  = 'carbonate chemistry -- '//'Summed occurrence of failure of the pH calculation.'
       call sub_adddef_netcdf(loc_iou,4,trim(loc_shortname),trim(loc_longname),' ',const_real_zero,const_real_zero)
       call sub_putvar3d_g(trim(loc_shortname),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! ---------------------------------------------------------------- !
    ! CARBONATE CHEMISTRY CONSTANTS (YAWN)
    ! ---------------------------------------------------------------- !
    DO icc=1,n_carbconst
       loc_ijk(:,:,:) = int_carbconst_timeslice(icc,:,:,:)/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,4, 'carbchem_const_'//trim(string_carbconst(icc)), &
            & 'carbonate chemistry dissociation constants - '//trim(string_carbconst(icc)),' ',const_real_zero,const_real_zero)
       call sub_putvar3d_g('carbchem_const_'//trim(string_carbconst(icc)),loc_iou, &
            & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    END DO
    ! ---------------------------------------------------------------- !
    ! ocean carbonate system isotopic properties
    ! ---------------------------------------------------------------- !
    If (ctrl_save_basic_proxies) then
       if (ocn_select(io_DIC_13C)) then
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   loc_ijk(i,j,k) = fun_calc_isotope_delta( &
                        & int_carb_timeslice(ic_conc_CO2,i,j,k), &
                        & int_carbisor_timeslice(ici_CO2_r13C,i,j,k)*int_carb_timeslice(ic_conc_CO2,i,j,k), &
                        & const_standards(11), &
                        & .FALSE., &
                        & const_real_null &
                        & )
                END DO
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,4,'carbchem_d13C_CO2','carbonate chemistry properties - '//'d13C of CO2(aq)','o/oo', &
               & const_real_zero,const_real_zero)
          call sub_putvar3d_g('carbchem_d13C_CO2',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   loc_ijk(i,j,k) = fun_calc_isotope_delta( &
                        & int_carb_timeslice(ic_conc_HCO3,i,j,k), &
                        & int_carbisor_timeslice(ici_HCO3_r13C,i,j,k)*int_carb_timeslice(ic_conc_HCO3,i,j,k), &
                        & const_standards(11), &
                        & .FALSE., &
                        & const_real_null &
                        & )
                END DO
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,4,'carbchem_d13C_HCO3','carbonate chemistry properties - '//'d13C of HCO3-','o/oo', &
               & const_real_zero,const_real_zero)
          call sub_putvar3d_g('carbchem_d13C_HCO3',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   loc_ijk(i,j,k) = fun_calc_isotope_delta( &
                        & int_carb_timeslice(ic_conc_CO3,i,j,k), &
                        & int_carbisor_timeslice(ici_CO3_r13C,i,j,k)*int_carb_timeslice(ic_conc_CO3,i,j,k), &
                        & const_standards(11), &
                        & .FALSE., &
                        & const_real_null &
                        & )
                END DO
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,4,'carbchem_d13C_CO32','carbonate chemistry properties - '//'d13C of CO32-','o/oo', &
               & const_real_zero,const_real_zero)
          call sub_putvar3d_g('carbchem_d13C_CO32',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end if
    end if
    ! ---------------------------------------------------------------- !
    ! GEOCHEMICAL DIAGNOSTICS -- REDOX -- POM-aqueous reactions
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_redox) then
       loc_unitsname = 'mol kg-1 yr-1'
       DO id=n_diag_redox_aq+1,n_diag_redox
          loc_ijk(:,:,:) = int_diag_redox_timeslice(id,:,:,:)/int_t_timeslice
          call sub_adddef_netcdf(loc_iou,4,'redox_'//trim(string_diag_redox(id)), &
               & 'redox transformation rate - '//trim(string_diag_redox(id)),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar3d_g('redox_'//trim(string_diag_redox(id)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end DO
    end If
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_geochemistry_advanced
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_biologicalpump_basic
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,is,l,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ijk(:,:,:) = const_real_zero
       !-------------------------------------------------------------- ! flux density
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_abio, &
                     & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                     & par_sed_type_scavenged)
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                   loc_unitsname = 'mol m-2 yr-1'
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                   loc_frac = int_bio_settle_timeslice(is,i,j,k)*phys_ocn(ipo_rA,i,j,k)/int_t_timeslice
                   loc_standard = const_standards(sed_type(is))
                   loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                   loc_unitsname = 'o/oo'
                CASE (par_sed_type_frac)
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)
                   loc_unitsname = 'n/a'
                end SELECT
             end do
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_abio, &
            & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
            & par_sed_type_scavenged,par_sed_type_frac)
          ! bulk tracers
          call sub_adddef_netcdf(loc_iou,4,'biop_fsinking_'//trim(string_sed(is)), &
               & 'particulate flux (density) - '//trim(string_sed(is)),loc_unitsname,const_real_zero,const_real_zero)
          call sub_putvar3d_g('biop_fsinking_'//trim(string_sed(is)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       CASE (n_itype_min:n_itype_max)
          ! isotopes
          If (ctrl_save_basic_proxies) then
             call sub_adddef_netcdf(loc_iou,4,'biop_fsinking_'//trim(string_sed(is)), &
                  & 'particulate flux (density) - '//trim(string_sed(is)),loc_unitsname,const_real_zero,const_real_zero)
             call sub_putvar3d_g('biop_fsinking_'//trim(string_sed(is)),loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end if
       end SELECT
    end do
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_biologicalpump_basic
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_biologicalpump_advanced
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,is,l,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! REMINERALIZATION FIELD
    ! ---------------------------------------------------------------- !
    DO l=3,n_l_ocn
       io = conv_iselected_io(l)
       is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
       if (is == 0) then
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   SELECT CASE (ocn_type(io))
                   CASE (0,1)
                      loc_ijk(i,j,k) = int_bio_remin_timeslice(io,i,j,k)/int_t_timeslice
                      loc_unitsname = 'mol yr-1'
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_bio_remin_timeslice(ocn_dep(io),i,j,k)/int_t_timeslice
                      loc_frac = int_bio_remin_timeslice(io,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(ocn_type(io))
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                      loc_unitsname = 'o/oo'
                   END SELECT
                end do
             end do
          end do
          SELECT CASE (ocn_type(io))
          CASE (1)
             ! bulk tracers
             If (ctrl_save_basic_reservoirs) then
                call sub_adddef_netcdf(loc_iou,4,'biop_remin_'//trim(string_ocn(io))//'_remin', &
                     & 'remineralization flux - '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar3d_g('biop_remin_'//trim(string_ocn(io))//'_remin',loc_iou, &
                     & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
             end if
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou,4,'biop_remin_'//trim(string_ocn(io))//'_remin', &
                     & 'remineralization flux - '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar3d_g('biop_remin_'//trim(string_ocn(io))//'_remin',loc_iou, &
                     & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
             end if
          END SELECT
       end if
    END DO
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ijk(:,:,:) = const_real_zero
       !-------------------------------------------------------------- ! total flux (per grid cell)
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_abio, &
                     & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                     & par_sed_type_scavenged)
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)/int_t_timeslice
                   loc_unitsname = 'mol yr-1'
                CASE (par_sed_type_frac)
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is,i,j,k)
                   loc_unitsname = 'n/a'
                end SELECT
             end do
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_abio, &
            & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
            & par_sed_type_scavenged, &
            & par_sed_type_frac)
          call sub_adddef_netcdf(loc_iou,4,'biop_fsinking_tot_'//trim(string_sed(is)), &
               & 'particulate flux (total) - '//trim(string_sed(is)),loc_unitsname,const_real_zero,const_real_zero)
          call sub_putvar3d_g('biop_fsinking_tot_'//trim(string_sed(is)), loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end SELECT
       !-------------------------------------------------------------- ! normalized flux
       DO i=1,n_i
          DO j=1,n_j
             if (int_bio_settle_timeslice(is,i,j,n_k) > 0.0) then
                loc_ijk(i,j,1:n_k-1) = int_bio_settle_timeslice(is,i,j,1:n_k-1)/int_bio_settle_timeslice(is,i,j,n_k)
                loc_ijk(i,j,n_k)     = 1.0
             else
                loc_ijk(i,j,:)       = 0.0
             end if
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_det)
          call sub_adddef_netcdf(loc_iou,4,'biop_fsinking_norm_'//trim(string_sed(is)), &
               & 'export-normalized particulate flux - '//trim(string_sed(is)),'n/a',const_real_zero,const_real_zero)
          call sub_putvar3d_g('biop_fsinking_norm_'//trim(string_sed(is)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end SELECT
       !
    END DO
    !----------------------------------------------------------------- ! CaCO3:POC 'rain ratio'
    IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'n/a'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (int_bio_settle_timeslice(is_POC,i,j,k) > const_real_nullsmall) then
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is_CaCO3,i,j,k)/int_bio_settle_timeslice(is_POC,i,j,k)
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_ratio_CaCO3toPOC','CaCO3 to POC rain ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_ratio_CaCO3toPOC',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    !----------------------------------------------------------------- ! opal:POC 'rain ratio'
    IF (sed_select(is_opal) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'n/a'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (int_bio_settle_timeslice(is_POC,i,j,k) > const_real_nullsmall) then
                   loc_ijk(i,j,k) = int_bio_settle_timeslice(is_opal,i,j,k)/int_bio_settle_timeslice(is_POC,i,j,k)
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_ratio_opaltoPOC','opal to POC rain ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_ratio_opaltoPOC',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    !----------------------------------------------------------------- ! Particulate flux fractions
    IF (sed_select(is_POC_frac2)) THEN
       loc_unitsname = 'n/a'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_bio_settle_timeslice(is_POC_frac2,i,j,k)/real(int_t_timeslice_count)
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_frac_POC2','POC fraction #2', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_frac_POC2',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    IF (sed_select(is_CaCO3_frac2)) THEN
       loc_unitsname = 'n/a'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_bio_settle_timeslice(is_CaCO3_frac2,i,j,k)/real(int_t_timeslice_count)
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_frac_CaCO32','CaCO3 fraction #2', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_frac_CaCO32',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    IF (sed_select(is_opal_frac2)) THEN
       loc_unitsname = 'n/a'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_bio_settle_timeslice(is_opal_frac2,i,j,k)/real(int_t_timeslice_count)
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_frac_opal2','opal fraction #2', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_frac_opal2',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! PARTICULATE CONCENTRATION FIELD
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_extra) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_ijk(:,:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_abio, &
                        & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                        & par_sed_type_scavenged)
                      loc_ijk(i,j,k) = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice
                      loc_unitsname = 'mol kg-1'
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_bio_part_timeslice(sed_dep(is),i,j,k)/int_t_timeslice
                      loc_frac = int_bio_part_timeslice(is,i,j,k)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                      loc_unitsname = 'o/oo'
                   CASE (par_sed_type_frac)
                      loc_ijk(i,j,k) = int_bio_part_timeslice(is,i,j,k)
                      loc_unitsname = 'n/a'
                   END SELECT
                end do
             end do
          end do
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged,par_sed_type_frac)
             call sub_adddef_netcdf(loc_iou,4,'biop_conc_'//trim(string_sed(is)), &
                  & 'particulate density - '//trim(string_sed(is)),loc_unitsname,const_real_zero,const_real_zero)
             call sub_putvar3d_g('biop_conc_'//trim(string_sed(is)),loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou,4,'biop_conc_'//trim(string_sed(is)), &
                     & 'particulate density - '//trim(string_sed(is)),loc_unitsname,const_real_zero,const_real_zero)
                call sub_putvar3d_g('biop_conc_'//trim(string_sed(is)),loc_iou, &
                     & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
             end if
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! N-star
    ! ---------------------------------------------------------------- !
    IF (ocn_select(io_PO4) .AND. ocn_select(io_NO3)) THEN
       loc_unitsname = 'umol kg-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) =                                                           &
                     & 1.0E6*                                                              &
                     & (                                                                   &
                     &   int_ocn_timeslice(io_NO3,i,j,k) -                                 &
                     &   bio_part_red(is_POP,is_PON,i,j)*int_ocn_timeslice(io_PO4,i,j,k) + &
                     &   par_bio_Nstar_offset                                              &
                     & )                                                                   &
                     & /int_t_timeslice
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_misc_Nstar','N-star', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_misc_Nstar',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! P-star
    ! ---------------------------------------------------------------- !
    IF (ocn_select(io_PO4) .AND. ocn_select(io_NO3)) THEN
       loc_unitsname = 'umol kg-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (bio_part_red(is_POP,is_PON,i,j) > const_real_nullsmall) then
                   loc_ijk(i,j,k) =                                                           &
                        & 1.0E6*                                                              &
                        & (                                                                   &
                        &   int_ocn_timeslice(io_PO4,i,j,k) -                                 &
                        &   int_ocn_timeslice(io_NO3,i,j,k)/bio_part_red(is_POP,is_PON,i,j)   &
                        & )                                                                   &
                        & /int_t_timeslice
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_misc_Pstar','P-star', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_misc_Pstar',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! DINex
    ! ---------------------------------------------------------------- !
    IF (ocn_select(io_PO4) .AND. ocn_select(io_NO3) .AND. ocn_select(io_NH4)) THEN
       loc_unitsname = 'umol kg-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) =                                                           &
                     & 1.0E6*                                                              &
                     & (                                                                   &
                     &   int_ocn_timeslice(io_NO3,i,j,k)+int_ocn_timeslice(io_NH4,i,j,k) - &
                     &   par_bio_red_POP_PON*int_ocn_timeslice(io_PO4,i,j,k)               &
                     & )                                                                   &
                     & /int_t_timeslice
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'biop_misc_DINex','DIN excess', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('biop_misc_DINex',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_biologicalpump_advanced
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_proxies_basic()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! Cd trace metal ratios (in sea-water)
    ! ---------------------------------------------------------------- !
    IF ( ocn_select(io_Cd) .AND. ocn_select(io_Ca) ) THEN
       loc_unitsname = 'nmol kg-1 (mmol kg-1)-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (int_ocn_timeslice(io_Ca,i,j,k) > const_real_nullsmall) then
                   loc_ijk(i,j,k) = 1.0E6*int_ocn_timeslice(io_Cd,i,j,k)/int_ocn_timeslice(io_Ca,i,j,k)
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_ratio_CdtoCa','Cd:Ca trace metal ratio (ocean)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_r_CdtoCa',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    IF (ocn_select(io_Cd) .AND. ocn_select(io_PO4)) THEN
       loc_unitsname = 'nmol kg-1 (umol kg-1)-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (int_ocn_timeslice(io_PO4,i,j,k) > const_real_nullsmall) then
                   loc_ijk(i,j,k) = 1.0E3*int_ocn_timeslice(io_Cd,i,j,k)/int_ocn_timeslice(io_PO4,i,j,k)
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_ratio_CdtoPO4','Cd:PO4 trace metal ratio (ocean)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_r_CdtoPO4',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    IF (ocn_select(io_Cd) .AND. ocn_select(io_PO4) .AND. ocn_select(io_Ca)) THEN
       loc_unitsname = 'umol kg-1 (mmol kg-1)-1'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                if (int_ocn_timeslice(io_Ca,i,j,k) > const_real_nullsmall) then
                   loc_ijk(i,j,k) = 1.0E3*int_ocn_timeslice(io_PO4,i,j,k)/int_ocn_timeslice(io_Ca,i,j,k)
                end if
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_ratio_PO4toCa','PO4:Ca ratio (ocean)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_r_PO4toCa',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! Os isotope ratios
    ! ---------------------------------------------------------------- !
    If ( ocn_select(io_Os) .AND. ocn_select(io_Os_187Os) ) then
       loc_unitsname = ''
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_ocn_timeslice(io_Os_187Os,i,j,k)/int_ocn_timeslice(io_Os_188Os,i,j,k)
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_187Osr188Os', &
            & 'water-column Os isotope ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_187Osr188Os',loc_iou, &
            & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end If
    ! ---------------------------------------------------------------- !
    ! Sr 87/86 isotope ratios
    ! ---------------------------------------------------------------- !
    If ( ocn_select(io_Sr) .AND. ocn_select(io_Sr_87Sr) ) then
       loc_unitsname = ''
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_ocn_timeslice(io_Sr_87Sr,i,j,k)/(ocn(io_Sr,i,j,k)-ocn(io_Sr_87Sr,i,j,k)-ocn(io_Sr_88Sr,i,j,k))
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_87Srr86Sr', &
            & 'water-column Sr isotope ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_87Srr86Sr',loc_iou, &
            & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end If
    ! ---------------------------------------------------------------- !
    ! Sr d88Sr isotope ratios
    ! ---------------------------------------------------------------- !
    If ( ocn_select(io_Sr) .AND. ocn_select(io_Sr_88Sr) ) then
       loc_unitsname = ''
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = fun_calc_isotope_deltaR(ocn(io_Sr,i,j,k)-ocn(io_Sr_87Sr,i,j,k)-ocn(io_Sr_88Sr,i,j,k), &
                     & ocn(io_Sr_88Sr,i,j,k),const_standardsR(ocn_type(io_Sr_88Sr)),const_real_null)
             END DO
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_d88Sr', &
            & 'water-column Sr isotope ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_d88Sr',loc_iou, &
            & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end If
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_proxies_basic
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_proxies_advanced()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! ocean carbonate system isotopic properties -- proxies
    ! ---------------------------------------------------------------- !
    If (ocn_select(io_DIC_13C)) then
       ! Schmittner LA1 d13C proxy
       loc_unitsname = 'o/oo'
       loc_ijk(:,:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = int_ocn_timeslice(io_DIC,i,j,k)
                loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,k)
                loc_standard = const_standards(ocn_type(io_DIC_13C))
                loc_ijk(i,j,k) = const_d13C_LA1_a + &
                     & const_d13C_LA1_b*fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null) + &
                     & const_d13C_LA1_c*1.0E6*int_carb_timeslice(ic_conc_CO3,i,j,k) + &
                     & const_d13C_LA1_d*int_phys_ocn_timeslice(ipo_Dmid,i,j,k)/int_t_timeslice
             end DO
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,4,'proxy_d13C_LA1','Schmittner d13C (LA1) proxy','o/oo',const_real_zero,const_real_zero)
       call sub_putvar3d_g('proxy_d13C_LA1',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end If
    ! ---------------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_proxies_advanced
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_hidden_grid()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    INTEGER::l,i,j,k,io,is,ip,ic,icc,loc_iou,loc_ntrec
    integer::id
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask,loc_sed_mask
    real::loc_tot,loc_frac,loc_standard
    real::loc_d13C,loc_d14C
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask        = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! phys_ocn -- OCEAN GRID
    ! ---------------------------------------------------------------- !
    if (ctrl_save_hidden_extra) then
       DO ip=1,n_phys_ocn
          loc_ijk(:,:,:) = int_phys_ocn_timeslice(ip,:,:,:)/int_t_timeslice
          call sub_adddef_netcdf(loc_iou,4,'3Dgrid_'//trim(string_phys_ocn(ip)), &
               & 'ocean physics & grid - '//trim(string_phys_ocn(ip)),' ',const_real_zero,const_real_zero)
          call sub_putvar3d_g('3Dgrid_'//trim(string_phys_ocn(ip)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       END DO
    else
       DO ip=1,n_phys_ocn
          loc_ijk(:,:,:) = int_phys_ocn_timeslice(ip,:,:,:)/int_t_timeslice
          SELECT CASE (ip)
          CASE (ipo_dD,ipo_A,ipo_V,ipo_M) 
             call sub_adddef_netcdf(loc_iou,4,'3Dgrid_'//trim(string_phys_ocn(ip)), &
                  & 'ocean grid - '//trim(string_phys_ocn(ip)),' ',const_real_zero,const_real_zero)
             call sub_putvar3d_g('3Dgrid_'//trim(string_phys_ocn(ip)),loc_iou, &
                  & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_hidden_grid
  ! ****************************************************************************************************************************** !
  

  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_hidden_climate()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    INTEGER::l,i,j,k,io,is,ip,ic,icc,loc_iou,loc_ntrec
    integer::id
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask,loc_sed_mask
    real::loc_tot,loc_frac,loc_standard
    real::loc_d13C,loc_d14C
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask        = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! OCEAN PROPERTIES
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ijk(:,:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                SELECT CASE (ocn_type(io))
                CASE (0)
                   if (io == io_T) then
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice - const_zeroC
                   else
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                   end if
                END SELECT
             end do
          end do
       end do
       SELECT CASE (ocn_type(io))
       CASE (0)
          ! T,S
          call sub_adddef_netcdf(loc_iou,4,'climate_'//trim(string_ocn_tname(l)), &
               & 'ocean climatology - '//trim(string_ocn_tlname(l)),trim(string_ocn_unit(l)),ocn_mima(l,1),ocn_mima(l,2))
          call sub_putvar3d_g('climate_'//trim(string_ocn(io)),loc_iou,n_i,n_j,n_k, &
               & loc_ntrec,loc_ijk(:,:,:),loc_mask)
       END SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! phys_ocn -- OCEAN 'CLIMATE'
    ! ---------------------------------------------------------------- !
    DO ip=1,n_phys_ocn
       loc_ijk(:,:,:) = int_phys_ocn_timeslice(ip,:,:,:)/int_t_timeslice
       SELECT CASE (ip)
       CASE (ipo_rho,ipo_gu,ipo_gv,ipo_gw) 
          call sub_adddef_netcdf(loc_iou,4,'climate_'//trim(string_phys_ocn(ip)), &
               & 'ocean climatology - '//trim(string_phys_ocn(ip)),' ',const_real_zero,const_real_zero)
          call sub_putvar3d_g('climate_'//trim(string_phys_ocn(ip)),loc_iou, &
               & n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
       end SELECT
    END DO
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_hidden_climate
  ! ****************************************************************************************************************************** !
  
  
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_3d_save_hidden_preformedtracers()
    ! ---------------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    integer::i,j,k,io,id
    INTEGER::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j,n_k)::loc_ijk,loc_mask
    CHARACTER(len=255)::loc_name
    CHARACTER(len=255)::loc_unitsname
    real::loc_tot,loc_frac,loc_standard
    real::loc_min,loc_max
    logical::loc_save
    ! ---------------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout3d_iou
    loc_ntrec = ncout3d_ntrec
    loc_ij(:,:)     = 0.0
    loc_ijk(:,:,:)  = 0.0
    loc_mask(:,:,:) = phys_ocn(ipo_mask_ocn,:,:,:)
    ! ---------------------------------------------------------------- !
    ! COLOR AGE TRACER
    ! ---------------------------------------------------------------- !
    if ( ctrl_force_ocn_age .AND. ocn_select(io_colr) ) then
       loc_unitsname = 'yrs'
       loc_ijk(:,:,:) = int_ocn_timeslice(io_colr,:,:,:)/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,4,'diag_age','color tracer; ventilation age','(yrs)',const_real_zero,const_real_zero)
       call sub_putvar3d_g('diag_age',loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! ---------------------------------------------------------------- !
    ! Preformed nutrients and things
    ! ---------------------------------------------------------------- !
    ! NOTE: save (loc_save == .true.) unless Csoft and not ctrl_bio_remin_redox_save
    do io=io_col0,io_col9
       if (ocn_select(io)) then
          loc_ijk(:,:,:) = const_real_null
          loc_unitsname = '???'
          loc_save = .true.
          loc_min = const_real_zero
          loc_max = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                DO k=goldstein_k1(i,j),n_k
                   select case (io)
                   CASE (io_col0:io_col6)
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                      loc_unitsname  = 'mol kg-1'
                   CASE (io_col7)
                      if (ocn_select(io_col0)) then
                         loc_tot        = int_ocn_timeslice(io_col0,i,j,k)/int_t_timeslice
                         loc_frac       = int_ocn_timeslice(io_col7,i,j,k)/int_t_timeslice
                         loc_standard   = const_standards(ocn_type(io_DIC_13C))
                         loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                         loc_unitsname  = 'o/oo'
                         loc_min        = -1000.0
                         loc_max        = +1000.0
                      end if
                   CASE (io_col8)
                      if (ocn_select(io_DIC_14C)) then
                         loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                         loc_unitsname  = 'yrs'
                      elseif (ocn_select(io_col9)) then
                         loc_tot        = int_ocn_timeslice(io_col9,i,j,k)/int_t_timeslice
                         loc_frac       = int_ocn_timeslice(io_col8,i,j,k)/int_t_timeslice
                         loc_standard   = const_standards(ocn_type(io_DIC_13C))
                         loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                         loc_unitsname  = 'o/oo'
                         loc_min        = -1000.0
                         loc_max        = +1000.0
                      end if
                   CASE (io_col9)
                      loc_ijk(i,j,k) = int_ocn_timeslice(io,i,j,k)/int_t_timeslice
                      loc_unitsname  = 'mol kg-1'
                   case default
                      ! NOTHING DOING
                   end select
                END DO
             END DO
          END DO
          loc_name = 'diag_pre_NULL_'//fun_conv_num_char_n(2,io)
          select case (io)
          CASE (io_col0)
             if (ocn_select(io_DIC))     loc_name = 'diag_pre_'//trim(string_ocn(io_DIC))
          CASE (io_col1)
             if (ocn_select(io_ALK))     loc_name = 'diag_pre_'//trim(string_ocn(io_ALK))
          CASE (io_col2)
             if (ocn_select(io_O2))      loc_name = 'diag_pre_'//trim(string_ocn(io_O2))
          CASE (io_col3)
             if (ocn_select(io_PO4))     loc_name = 'diag_pre_'//trim(string_ocn(io_PO4))
          CASE (io_col4)
             if (ocn_select(io_NO3))     loc_name = 'diag_pre_'//trim(string_ocn(io_NO3))
          CASE (io_col5)
             if (ocn_select(io_Fe))      loc_name = 'diag_pre_'//trim(string_ocn(io_Fe))
             if (ocn_select(io_TDFe))    loc_name = 'diag_pre_'//trim(string_ocn(io_TDFe))
          CASE (io_col6)
             if (ocn_select(io_SiO2))    loc_name = 'diag_pre_'//trim(string_ocn(io_SiO2))
          CASE (io_col7)
             if (ocn_select(io_DIC_13C)) loc_name = 'diag_pre_'//trim(string_ocn(io_DIC_13C))
             if (.NOT. ocn_select(io_col0)) loc_save = .false.
          CASE (io_col8)
             if (ocn_select(io_DIC_14C)) then
                loc_name = 'diag_pre_d14C_age'
             else
                loc_name = 'diag_pre_Csoft_d13C'
                if (.NOT. ctrl_bio_remin_redox_save) loc_save = .false.
                if (.NOT. ocn_select(io_col9)) loc_save = .false.
             end if
          CASE (io_col9)
             if (ocn_select(io_DIC)) loc_name = 'diag_pre_Csoft'
             if (.NOT. ctrl_bio_remin_redox_save) loc_save = .false.
          end select
          ! write to netCDF
          if (loc_save) then
             call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Preformed tracer',trim(loc_unitsname),loc_min,loc_max)
             call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
          end if
       end if
    end do
    ! ---------------------------------------------------------------- !
    ! Related preformed metrics! :o)
    ! (1)  diag_reg_O2sat -- O2(sat)
    ! (2)  diag_reg_AOU -- AOU (from O2(sat) and O2)
    ! (3)  diag_reg_AOU_P -- AOU-regenP (regen-P from AOU) (assuming Redfield)
    ! (4)  diag_reg_AOU_P_C -- AOU-regenP-regenC (regenC from regenP from AOU) (assuming Redfield)
    ! (5)  diag_reg_AOU_P_C_d13C -- AOU-regenP-regenC d13C (AOU-regenP-regenC d13C in proportion to AOU-regenP-regenC/DIC)
    !      (assuming POC d13C)
    ! (6)  diag_reg_OU -- OU (from O2(pre) and actual O2)
    ! (7)  diag_reg_P -- regenP (from PO4(pre) and actual PO4)
    ! (8)  diag_reg_P_C -- regenP-regenC (regenC from regenP) (assuming Redfield)
    ! (9)  diag_reg_P_C_d13C -- regenP-regenC d13C (regenP-regenC d13C in proportion to regenP-regenC/DIC)
    !      (assuming POC d13C)
    ! (10) diag_reg_C -- regenC (Csoft (real regenC))
    ! (11) diag_reg_C_d13C -- regenC d13C (regenC d13C in proportion to regenC/DIC) (assuming POC d13C)
    ! (12) diag_reg_d13C -- regend13C (Csoft d13C in proportion to Csoft/DIC)
    ! (13) diag_reg_DIC -- DIC (from PreC + regenC)
    ! (14) diag_reg_DIC_d13C -- DIC d13C (from PreC d13C + regenC d13C)
    ! (15) Ccarb (DIC minus (DIC(pre) plus Csoft))
    ! (16) Ccarb d13C (DIC d13C minus (d13C(pre) plus Csoft d13C)) in respective proportions ...
    ! ---------------------------------------------------------------- !
    ! NOTE: all is deliberately and tediously re-done behind ctrl_bio_preformed and then behind seperate 'ifs' 
    !       to make the code clearer and mroe transparent(?)
    ! (1) O2(sat)
    if (ocn_select(io_O2)) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_O2sat'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                ! calculate solubility constant 
                ! NOTE: that each time-step, only the surface values are calculated (hence need for full 3D uopdate here)
                loc_ijk(i,j,k) = fun_calc_solconst(ia_pO2, &
                     & int_ocn_timeslice(io_T,i,j,k)/int_t_timeslice, &
                     & int_ocn_timeslice(io_S,i,j,k)/int_t_timeslice, &
                     & int_phys_ocn_timeslice(ipo_rho,i,j,k)/int_t_timeslice &
                     & )
                ! calculate O2(sat) (scale to pO2)
                loc_ijk(i,j,k) = (int_sfcatm1_timeslice(ia_pO2,i,j)/int_t_timeslice)*loc_ijk(i,j,k)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (2) AOU (from O2(sat) and O2)
    ! NOTE: utilize the fact that loc_ijk already contains the value of O2(sat)
    ! NOTE: loop rather than manipulate the full arrays directly to preserve the null values for non ocean points
    !       (doing null minus null might be messy, and we don't like messy here)
    if (ocn_select(io_O2)) then
       loc_name       = 'diag_reg_AOU'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = loc_ijk(i,j,k) - int_ocn_timeslice(io_O2,i,j,k)/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (3) AOU-regenP (regen P from AOU)
    ! NOTE: cascade the values of loc_ijk ... containing AOU ... onwards! 
    ! NOTE: default: par_bio_red_POP_PO2=-138.0 -- invert sign as AOU is positive
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    if (ocn_select(io_O2)) then
       loc_name       = 'diag_reg_AOU_P'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = loc_ijk(i,j,k)/(-1.0*par_bio_red_POP_PO2)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (4) AOU-regenC (regen C from regen P from AOU)
    ! NOTE: cascade the values of loc_ijk ... containing regen P ... from AOU ... onwards!
    ! NOTE: default: par_bio_red_POP_POC=106.0
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    if (ocn_select(io_O2)) then
       loc_name       = 'diag_reg_AOU_P_C'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = par_bio_red_POP_POC*loc_ijk(i,j,k)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (5) AOU-regen13C (regen 13C regen C from from regen P from AOU)
    ! NOTE: cascade the values of loc_ijk ... containing regen C from regen P from AOU ... onwards! 
    !       To the moon!
    ! NOTE: take surface export POC d13C, then 'dilute' isotopic signature in the proportion of regen C over DIC
    !       => need to calculate d13C from POC and POC_13C
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    ! NOTE: no need to change units of int_bio_settle_timeslice
    if (ocn_select(io_O2) .AND. ocn_select(io_DIC_13C)) then
       loc_name       = 'diag_reg_AOU_P_C_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = int_bio_settle_timeslice(is_POC,i,j,n_k)
                loc_frac = int_bio_settle_timeslice(is_POC_13C,i,j,n_k)
                loc_standard = const_standards(sed_type(is_POC_13C))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_ijk(i,j,k) = loc_ij(i,j)*loc_ijk(i,j,k)/(int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (6) real OU (preformed O2 minus actual O2)
    if (ocn_select(io_O2) .AND. ocn_select(io_col2)) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_OU'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = (int_ocn_timeslice(io_col2,i,j,k) - int_ocn_timeslice(io_O2,i,j,k))/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (7) real regen P (actual PO4 minus preformed P)
    if (ocn_select(io_PO4) .AND. ocn_select(io_col3)) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_P'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = (int_ocn_timeslice(io_PO4,i,j,k) - int_ocn_timeslice(io_col3,i,j,k))/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (8) regenP-regenC (regen C from regen P)
    ! NOTE: cascade the values of loc_ijk ... containing regen P ... onwards!
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    if (ocn_select(io_PO4) .AND. ocn_select(io_col3)) then
       loc_name       = 'diag_reg_P_C'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = par_bio_red_POP_POC*loc_ijk(i,j,k)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (9) regenP-regenC regen d13C (regen 13C from regen C from regen P)
    ! NOTE: cascade the values of loc_ijk ... containing regen C from regen P ... onwards!
    ! NOTE: loc_ij already contains the values of POC d13C!
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    if (ocn_select(io_PO4) .AND. ocn_select(io_col3) .AND. ocn_select(io_DIC_13C)) then
       loc_name       = 'diag_reg_P_C_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = loc_ij(i,j)*loc_ijk(i,j,k)/(int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (10) Csoft regenC for completeness
    if (ocn_select(io_col9) .AND. ctrl_bio_remin_redox_save) then
       loc_name       = 'diag_reg_C'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_ocn_timeslice(io_col9,i,j,k)/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (11) regenC regen d13C (regen 13C from regen C)
    ! NOTE: cascade the values of loc_ijk ... containing regen C ... onwards!
    ! NOTE: loc_ij already contains the values of POC d13C!
    ! NOTE: int_t_timeslice already taken account of in loc_ijk
    if (ocn_select(io_PO4) .AND. ocn_select(io_col3) .AND. ocn_select(io_DIC_13C)) then
       loc_name       = 'diag_reg_C_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = loc_ij(i,j)*loc_ijk(i,j,k)/(int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (12) regen d13C (diag_reg_Csoft_d13C in proportion to Csoft/DIC)
    ! NOTE: need to calculate regen d13C from regen C and regen 13C
    ! NOTE: remember D14C needs to be NOT selected
    ! NOTE: int_t_timeslices cancel out
    ! NOTE: no need to change units of regen C and regen 13C
    if ( &
         & ocn_select(io_DIC_13C) .AND. (.NOT. ocn_select(io_DIC_14C)) .AND. ctrl_bio_remin_redox_save .AND. &
         & ocn_select(io_col8) .AND. ocn_select(io_col9) &
         & ) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = int_ocn_timeslice(io_col9,i,j,k)
                loc_frac = int_ocn_timeslice(io_col8,i,j,k)
                loc_standard = const_standards(ocn_type(io_DIC_13C))
                loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_ijk(i,j,k) = loc_ijk(i,j,k)*int_ocn_timeslice(io_col9,i,j,k)/int_ocn_timeslice(io_DIC,i,j,k)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (13) DIC (from PreC + regenC)
    if (ocn_select(io_col0) .AND. ocn_select(io_col9) .AND. ctrl_bio_remin_redox_save) then
       loc_name       = 'diag_reg_DIC'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = (int_ocn_timeslice(io_col0,i,j,k) + int_ocn_timeslice(io_col9,i,j,k))/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (14) DIC d13C (from PreC d13C + regenC d13C)
    ! NOTE: int_t_timeslices cancel out
    if ( &
         & ocn_select(io_DIC_13C) .AND. (.NOT. ocn_select(io_DIC_14C)) .AND. ctrl_bio_remin_redox_save .AND. &
         & ocn_select(io_col0) .AND. ocn_select(io_col9) .AND. ocn_select(io_col7) .AND. ocn_select(io_col8) &
         & ) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_DIC_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = (int_ocn_timeslice(io_col0,i,j,k) + int_ocn_timeslice(io_col9,i,j,k))
                loc_frac = (int_ocn_timeslice(io_col7,i,j,k) + int_ocn_timeslice(io_col8,i,j,k))
                loc_standard = const_standards(ocn_type(io_DIC_13C))
                loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (15) Ccarb (DIC minus (DIC(pre) plus Csoft))
    if (ocn_select(io_col0) .AND. ocn_select(io_col9) .AND. ctrl_bio_remin_redox_save) then
       loc_name       = 'diag_reg_CaCO3'
       loc_unitsname  = 'mol kg-1'
       loc_min        = -1.0
       loc_max        = +1.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_ijk(i,j,k) = int_ocn_timeslice(io_DIC,i,j,k)/int_t_timeslice - &
                     & (int_ocn_timeslice(io_col0,i,j,k) + int_ocn_timeslice(io_col9,i,j,k))/int_t_timeslice
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! (16) Ccarb d13C (DIC d13C minus d13C(pre) minus Csoft d13C in respective proportions ...
    ! NOTE: int_t_timeslices cancel out
    if ( &
         & ocn_select(io_DIC_13C) .AND. (.NOT. ocn_select(io_DIC_14C)) .AND. ctrl_bio_remin_redox_save .AND. &
         & ocn_select(io_col0) .AND. ocn_select(io_col9) .AND. ocn_select(io_col7) .AND. ocn_select(io_col8) &
         & ) then
       loc_ijk(:,:,:) = const_real_null
       loc_name       = 'diag_reg_CaCO3_d13C'
       loc_unitsname  = 'o/oo'
       loc_min        = -1000.0
       loc_max        = +1000.0
       DO i=1,n_i
          DO j=1,n_j
             DO k=goldstein_k1(i,j),n_k
                loc_tot  = int_ocn_timeslice(io_DIC,i,j,k) - &
                     & (int_ocn_timeslice(io_col0,i,j,k) + int_ocn_timeslice(io_col9,i,j,k))
                loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,k) - &
                     & (int_ocn_timeslice(io_col7,i,j,k) + int_ocn_timeslice(io_col8,i,j,k))
                loc_standard = const_standards(ocn_type(io_DIC_13C))
                loc_ijk(i,j,k) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_ijk(i,j,k) = loc_ijk(i,j,k)*loc_tot/int_ocn_timeslice(io_DIC,i,j,k)
             END DO
          END DO
       END DO
       ! write to netCDF
       call sub_adddef_netcdf(loc_iou,4,trim(loc_name),'Regenerated tracer',trim(loc_unitsname),loc_min,loc_max)
       call sub_putvar3d_g(trim(loc_name),loc_iou,n_i,n_j,n_k,loc_ntrec,loc_ijk(:,:,:),loc_mask)
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_3d_save_hidden_preformedtracers
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !
  ! 2-D FIELDS
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !
 

  ! ****************************************************************************************************************************** !
  ! *** save time-slice data ***
  SUBROUTINE sub_save_netcdf_2d(dum_dtyr)
    ! ---------------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! ---------------------------------------------------------------- !
    real,intent(in)::dum_dtyr
    ! ---------------------------------------------------------------- !
    ! reservoirs
    ! ---------------------------------------------------------------- !
    If (ctrl_save_basic_reservoirs) CALL sub_2d_save_reservoirs_basic()
    If (ctrl_save_advanced_reservoirs) CALL sub_2d_save_reservoirs_advanced()
    ! ---------------------------------------------------------------- !
    ! geochemsitry
    ! ---------------------------------------------------------------- !
    If (ctrl_save_basic_geochemistry) CALL sub_2d_save_geochemistry_basic()
    If (ctrl_save_advanced_geochemistry) CALL sub_2d_save_geochemistry_advanced()
    ! ---------------------------------------------------------------- !
    ! biological pump
    ! ---------------------------------------------------------------- !
    If (ctrl_save_basic_biologicalpump) CALL sub_2d_save_biologicalpump_basic()
    If (ctrl_save_advanced_biologicalpump) CALL sub_2d_save_biologicalpump_advanced()
    ! ---------------------------------------------------------------- !
    ! proxies
    ! ---------------------------------------------------------------- !
    If (ctrl_save_basic_proxies) CALL sub_2d_save_proxies_basic()
    If (ctrl_save_advanced_proxies) CALL sub_2d_save_proxies_advanced()
    ! ---------------------------------------------------------------- !
    ! hidden
    ! ---------------------------------------------------------------- !
    !----------------------------------------------------------------- ! grid
    If (ctrl_save_hidden_grid) CALL sub_2d_save_hidden_grid()
    !----------------------------------------------------------------- ! climate
    If (ctrl_save_hidden_climate) CALL sub_2d_save_hidden_climate(dum_dtyr)
    !----------------------------------------------------------------- ! GENIE exercises
    If (ctrl_save_hidden_fossilfuelco2) CALL sub_2d_save_hidden_fossilfuelco2()
    ! ---------------------------------------------------------------- ! interface (fluxes)
    if (ctrl_save_hidden_interfacefluxes) CALL sub_2d_save_hidden_interfacefluxes()
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_save_netcdf_2d
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save geochemical composition data ***
  SUBROUTINE sub_2d_save_reservoirs_basic()  
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,is,io,id
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    real::loc_tot,loc_frac,loc_standard
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask_sur,loc_mask_sed
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0   
    loc_mask_sur(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_sed(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
    ! ---------------------------------------------------------------- !
    ! ocean surface tracer properties
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                SELECT CASE (ocn_type(io))
                CASE (0)
                   if (io == io_T) then
                      loc_ij(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice - const_zeroC
                   else
                      loc_ij(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice
                   end if
                CASE (1)
                   loc_ij(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,n_k)
                   loc_frac = int_ocn_timeslice(io,i,j,n_k)
                   loc_standard = const_standards(ocn_type(io))
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                end SELECT
             end IF
          end DO
       end DO
       SELECT CASE (ocn_type(io))
       CASE (0)
          If (io == io_T) then
             loc_unitsname = 'degrees C'
          else
             loc_unitsname = 'o/oo'
          end If
       CASE (1)
          loc_unitsname = 'mol kg-1'
       case (n_itype_min:n_itype_max)
          loc_unitsname = 'o/oo'
       end SELECT
       SELECT CASE (ocn_type(io))
       CASE (0)
          ! T,S
             call sub_adddef_netcdf(loc_iou, 3, 'ocn_seasur_'//trim(string_ocn(io)), &
                  & 'surface-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('ocn_seasur_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
       CASE (1)
          ! bulk tracers
             call sub_adddef_netcdf(loc_iou, 3, 'ocn_seasur_'//trim(string_ocn(io)), &
                  & 'surface-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('ocn_seasur_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
       CASE (n_itype_min:n_itype_max)
          ! isotopes
          If (ctrl_save_basic_proxies) then
             call sub_adddef_netcdf(loc_iou, 3, 'ocn_seasur_'//trim(string_ocn(io)), &
                  & 'surface-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('ocn_seasur_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
          end if
       end SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! seafloor tracer properties
    ! ---------------------------------------------------------------- !
    if (ctrl_save_hidden_seafloor) then
       DO l=1,n_l_ocn
          io = conv_iselected_io(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   SELECT CASE (ocn_type(io))
                   CASE (0)
                      if (io == io_T) then
                         loc_ij(i,j) = int_ocn_timeslice(io,i,j,loc_k1)/int_t_timeslice - const_zeroC
                      else
                         loc_ij(i,j) = int_ocn_timeslice(io,i,j,loc_k1)/int_t_timeslice
                      end if
                   CASE (1)
                      loc_ij(i,j) = int_ocn_timeslice(io,i,j,loc_k1)/int_t_timeslice
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,loc_k1)
                      loc_frac = int_ocn_timeslice(io,i,j,loc_k1)
                      loc_standard = const_standards(ocn_type(io))
                      loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                   end SELECT
                end IF
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (0)
             If (io == io_T) then
                loc_unitsname = 'degrees C'
             else
                loc_unitsname = 'o/oo'
             end If
          CASE (1)
             loc_unitsname = 'mol kg-1'
          case (n_itype_min:n_itype_max)
             loc_unitsname = 'o/oo'
          end SELECT
          SELECT CASE (ocn_type(io))
          CASE (0)
             ! T,S
             call sub_adddef_netcdf(loc_iou, 3, 'ocn_seafloor_'//trim(string_ocn(io)), &
                  & 'bottom-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('ocn_seafloor_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
          CASE (1)
             ! bulk tracers
             call sub_adddef_netcdf(loc_iou, 3, 'ocn_seafloor_'//trim(string_ocn(io)), &
                  & 'bottom-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('ocn_seafloor_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou, 3, 'ocn_seafloor_'//trim(string_ocn(io)), &
                     & 'bottom-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar2d('ocn_seafloor_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
             end if
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! save core-top sediment composition data
    ! ---------------------------------------------------------------- !
    If (flag_sedgem) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio)
             loc_unitsname = 'wt%'
          CASE (par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det,par_sed_type_scavenged)
             loc_unitsname = 'ppm'
          CASE (n_itype_min:n_itype_max)
             loc_unitsname = 'o/oo'
          CASE (par_sed_type_age)
             loc_unitsname = 'years'
          end SELECT
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged,par_sed_type_age)
             call sub_adddef_netcdf(loc_iou,3,'sed_seafloor_'//trim(string_sed(is)), &
                  & 'sediment core-top '//trim(string_sed(is)),trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('sed_seafloor_'//trim(string_sed(is)),loc_iou,n_i,n_j, &
                  & loc_ntrec,int_sfcsed1_timeslice(is,:,:)/int_t_timeslice,loc_mask_sed)
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou,3,'sed_seafloor_'//trim(string_sed(is)), &
                     & 'sediment core-top '//trim(string_sed(is)),trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar2d('sed_seafloor_'//trim(string_sed(is)),loc_iou,n_i,n_j, &
                     & loc_ntrec,int_sfcsed1_timeslice(is,:,:)/int_t_timeslice,loc_mask_sed)
             end if
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_reservoirs_basic
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save geochemical composition data ***
  SUBROUTINE sub_2d_save_reservoirs_advanced()  
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,ia,io
    integer::loc_iou,loc_ntrec
    real::loc_tot,loc_frac,loc_standard
    real::loc_d13C,loc_d14C
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask_sur,loc_mask_sur_ALL
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_mask_sur(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_sur_ALL(:,:) = const_real_one
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED TRACER INVENTORIES
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol m-2'
    DO l=3,n_l_ocn
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = &
                     & sum(phys_ocn(ipo_M,i,j,:)*int_ocn_timeslice(io,i,j,:))* &
                     & phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
             end If
          end DO
       end DO
       SELECT CASE (ocn_type(io))
       CASE (1)
          call sub_adddef_netcdf(loc_iou,3,'ocn_wcint_'//trim(string_ocn(io)), &
               & trim(string_ocn(io))//' water-column integrated tracer inventory', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('ocn_wcint_'//trim(string_ocn(io)),loc_iou, &
               & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_sur)
       END SELECT
    end DO
    ! ---------------------------------------------------------------- !
    ! save atmosphere tracer field
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_atm
       ia = conv_iselected_ia(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             SELECT CASE (atm_type(ia))
             CASE (0,1)
                loc_ij(i,j) = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice
             case (n_itype_min:n_itype_max)
                loc_tot  = int_sfcatm1_timeslice(atm_dep(ia),i,j)/int_t_timeslice
                loc_frac = int_sfcatm1_timeslice(ia,i,j)/int_t_timeslice
                loc_standard = const_standards(atm_type(ia))
                loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
             END SELECT
          end DO
       end DO
       SELECT CASE (atm_type(ia))
       CASE (0,1)
          ! bulk atmosphere (no-one really needs to see this!)
             call sub_adddef_netcdf(loc_iou,3,'atm_'//trim(string_atm_tname(l)), &
                  & trim(string_atm_tlname(l)),trim(string_atm_unit(l)),atm_mima(l,1),atm_mima(l,2))
             call sub_putvar2d('atm_'//trim(string_atm(ia)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur_ALL)
       CASE (n_itype_min:n_itype_max)
          ! bulk atmosphere isotopes (no-one really needs to see this!)
          If (ctrl_save_basic_proxies) then
             call sub_adddef_netcdf(loc_iou,3,'atm_'//trim(string_atm_tname(l)), &
                  & trim(string_atm_tlname(l)),trim(string_atm_unit(l)),atm_mima(l,1),atm_mima(l,2))
             call sub_putvar2d('atm_'//trim(string_atm(ia)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur_ALL)
          end if
       END SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! save atmospheric 14C data in dumb-ass D14C units
    ! ---------------------------------------------------------------- !
    ! NOTE: saved as time-series (no-one really needs to see 2D)
    If (ctrl_save_basic_proxies) then
       IF (atm_select(ia_pCO2_13C) .AND. atm_select(ia_pCO2_14C)) THEN
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_tot  = int_sfcatm1_timeslice(ia_pCO2,i,j)/int_t_timeslice
                loc_frac = int_sfcatm1_timeslice(ia_pCO2_13C,i,j)/int_t_timeslice
                loc_standard = const_standards(atm_type(ia_pCO2_13C))
                loc_d13C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_frac = int_sfcatm1_timeslice(ia_pCO2_14C,i,j)/int_t_timeslice
                loc_standard = const_standards(atm_type(ia_pCO2_14C))
                loc_d14C = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                loc_ij(i,j) = fun_convert_delta14CtoD14C(loc_d13C,loc_d14C)
                loc_unitsname = 'o/oo'
             end do
          end do
          call sub_adddef_netcdf(loc_iou,3,'atm_pCO2_D14C', &
               & ' atmospheric pCO2 D14C',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('atm_pCO2_D14C',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur_ALL)
       end if
    end if
    !--------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_reservoirs_advanced
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save basic geochemsitry ***
  SUBROUTINE sub_2d_save_geochemistry_basic()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::l,i,j,ia,io,is
    integer::ib,id,ip,ic
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    real::loc_tot,loc_frac,loc_standard
    real::loc_d13C,loc_d14C
    real,DIMENSION(n_i,n_j)::loc_ij,loc_ij_tot,loc_mask_surf,loc_mask_surf_ALL,loc_sed_mask
    real,DIMENSION(n_i,n_j)::loc_ij_OC,loc_ij_NC,loc_ij_SC
    real,DIMENSION(n_i,n_j)::loc_ij_ON,loc_ij_NN,loc_ij_SN
    real,DIMENSION(n_i,n_j)::loc_ij_OP,loc_ij_NP,loc_ij_SP
    CHARACTER(len=31)::loc_string     !
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:)            = 0.0
    loc_mask_surf(:,:)     = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_surf_ALL(:,:) = 1.0
    loc_sed_mask(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
    ! ---------------------------------------------------------------- !
    ! surface ocean carbonate chemsitry
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       DO ic=1,n_carb
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   loc_ij(i,j) = int_carb_timeslice(ic,i,j,n_k)/int_t_timeslice
                end IF
             end DO
          end DO
          SELECT CASE (ic)
          CASE (ic_pHsws)
             loc_unitsname = 'pH(sws)'                   
          CASE (ic_ohm_cal,ic_ohm_arg,ic_RF0)     
             loc_unitsname = 'n/a' 
          case default
             loc_unitsname = 'mol kg-1'
          end SELECT
          SELECT CASE (ic)
          CASE (ic_pHsws,ic_conc_CO2,ic_conc_CO3,ic_conc_HCO3,ic_ohm_cal,ic_ohm_arg,ic_RF0)
             call sub_adddef_netcdf(loc_iou,3,'carbchem_seasur_'//trim(string_carb(ic)), &
                  & 'Ocean surface carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('carbchem_seasur_'//trim(string_carb(ic)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
          end select
       END DO
    end if
    ! ---------------------------------------------------------------- !
    ! seafloor carbonate chemsitry
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       if (ctrl_save_hidden_seafloor) then
          DO ic=1,n_carb
             loc_ij(:,:) = const_real_zero
             DO i=1,n_i
                DO j=1,n_j
                   loc_k1 = goldstein_k1(i,j)
                   IF (n_k >= loc_k1) THEN
                      loc_ij(i,j) = int_carb_timeslice(ic,i,j,loc_k1)/int_t_timeslice
                   end IF
                end DO
             end DO
             SELECT CASE (ic)
             CASE (ic_pHsws)
                loc_unitsname = 'pH(sws)'                
             CASE (ic_ohm_cal,ic_ohm_arg)
                loc_unitsname = 'n/a'
             case default
                loc_unitsname = 'mol kg-1'
             end SELECT
             SELECT CASE (ic)
             CASE (ic_pHsws,ic_ohm_cal,ic_ohm_arg,ic_dCO3_cal,ic_dCO3_arg)
                call sub_adddef_netcdf(loc_iou,3,'carbchem_seafloor_'//trim(string_carb(ic)), &
                     & 'Benthic carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                     & trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar2d('carbchem_seafloor_'//trim(string_carb(ic)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
             end select
          END DO
       end if
    end if
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED REDOX TRANSFORMATION ANALYSIS
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_redox) then
       ! NOTE: divide by conv_sed_ocn_x(io_y,is_POz) to convert from oxidant consumed to POC remineralized
       !       sign switch (for consumption flux) is done by dividing by conv_sed_ocn_x(io_y,is_POz)
       ! NOTE: POM only (and hence ocean interior)
       ! NOTE: scale int_diag_redox_timeslice by area and integration time
       loc_ij_OC(:,:)  = const_real_zero
       loc_ij_NC(:,:)  = const_real_zero
       loc_ij_SC(:,:)  = const_real_zero
       loc_ij_ON(:,:)  = const_real_zero
       loc_ij_NN(:,:)  = const_real_zero
       loc_ij_SN(:,:)  = const_real_zero
       loc_ij_OP(:,:)  = const_real_zero
       loc_ij_NP(:,:)  = const_real_zero
       loc_ij_SP(:,:)  = const_real_zero
       ! -------------------------------------------------------- ! (1) oxic remineralization
       if (ocn_select(io_O2)) then
          loc_ij_tot(:,:) = const_real_zero
          loc_unitsname = 'mol C+N+P m-2 yr-1'
          loc_shortname = 'redox_wcint_remin_O2_POM'
          loc_longname  = 'total water-column integrated remineraliaztion of POM by O2'
          DO i=1,n_i
             DO j=1,n_j
                If (goldstein_k1(i,j) <= n_k) then
                   if (sed_select(is_POC) .AND. (abs(conv_sed_ocn_O(io_O2,is_POC)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POC_dO2'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_OC(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_O(io_O2,is_POC)
                      loc_ij_OC(i,j) = loc_ij_OC(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_PON) .AND. (abs(conv_sed_ocn_O(io_O2,is_PON)) > const_real_nullsmall)) then
                      loc_string = 'reminP_PON_dO2'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_ON(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_O(io_O2,is_PON)
                      loc_ij_ON(i,j) = loc_ij_ON(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_POP) .AND. (abs(conv_sed_ocn_O(io_O2,is_POP)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POP_dO2'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_OP(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_O(io_O2,is_POP)
                      loc_ij_OP(i,j) = loc_ij_OP(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                end if
             end DO
          end DO
          loc_ij_tot(:,:) = loc_ij_OC(:,:) + loc_ij_ON(:,:) + loc_ij_OP(:,:)
          call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
               & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij_tot(:,:),loc_mask_surf)
       end if
       ! -------------------------------------------------------- ! (2) denitrification
       if (ocn_select(io_NO3)) then
          loc_unitsname = 'mol C+N+P m-2 yr-1'
          loc_shortname = 'redox_wcint_remin_NO3_POM' 
          loc_longname  = 'total water-column integrated remineraliaztion of POM by NO3'
          DO i=1,n_i
             DO j=1,n_j
                If (goldstein_k1(i,j) <= n_k) then
                   if (sed_select(is_POC) .AND. (abs(conv_sed_ocn_N(io_NO3,is_POC)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POC_dNO3'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_NC(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_N(io_NO3,is_POC)
                      loc_ij_NC(i,j) = loc_ij_NC(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_PON) .AND. (abs(conv_sed_ocn_N(io_NO3,is_PON)) > const_real_nullsmall)) then
                      loc_string = 'reminP_PON_dNO3'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_NN(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_N(io_NO3,is_PON)
                      loc_ij_NN(i,j) = loc_ij_NN(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_POP) .AND. (abs(conv_sed_ocn_N(io_NO3,is_POP)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POP_dNO3'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_NP(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_N(io_NO3,is_POP)
                      loc_ij_NP(i,j) = loc_ij_NP(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                end if
             end DO
          end DO
          loc_ij_tot(:,:) = loc_ij_NC(:,:) + loc_ij_NN(:,:) + loc_ij_NP(:,:)
          call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
               & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij_tot(:,:),loc_mask_surf)
       end if
       ! -------------------------------------------------------- ! (3) sulphate reduction
       if (ocn_select(io_SO4)) then
          loc_unitsname = 'mol C+N+P m-2 yr-1'
          loc_shortname = 'redox_wcint_remin_SO4_POM' 
          loc_longname  = 'total water-column integrated remineraliaztion of POM by SO4'  
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                If (goldstein_k1(i,j) <= n_k) then
                   if (sed_select(is_POC) .AND. (abs(conv_sed_ocn_S(io_SO4,is_POC)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POC_dSO4'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_SC(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_S(io_SO4,is_POC)
                      loc_ij_SC(i,j) = loc_ij_SC(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_PON) .AND. (abs(conv_sed_ocn_S(io_SO4,is_PON)) > const_real_nullsmall)) then
                      loc_string = 'reminP_PON_dSO4'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_SN(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_S(io_SO4,is_PON)
                      loc_ij_SN(i,j) = loc_ij_SN(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                   if (sed_select(is_POP) .AND. (abs(conv_sed_ocn_S(io_SO4,is_POP)) > const_real_nullsmall)) then
                      loc_string = 'reminP_POP_dSO4'
                      id = fun_find_str_i(trim(loc_string),string_diag_redox)
                      loc_ij_SP(i,j) = sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))/conv_sed_ocn_S(io_SO4,is_POP)
                      loc_ij_SP(i,j) = loc_ij_SP(i,j)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   end if
                end if
             end DO
          end DO
          loc_ij_tot(:,:) = loc_ij_SC(:,:) + loc_ij_SN(:,:) + loc_ij_SP(:,:)
          call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
               & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij_tot(:,:),loc_mask_surf)
       end if
       ! -------------------------------------------------------- ! (4) totals
       loc_unitsname = 'mol C m-2 yr-1'
       loc_shortname = 'redox_wcint_ALL_POC'
       loc_longname  = 'total water-column integrated remineralization of POC by all oxidants'
       loc_ij_tot(:,:) = loc_ij_OC(:,:) + loc_ij_NC(:,:) + loc_ij_SC(:,:)
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij_tot(:,:),loc_mask_surf)
       loc_unitsname = 'mol C+N+P m-2 yr-1'
       loc_shortname = 'redox_wcint_ALL_POM'
       loc_longname  = 'total water-column integrated remineralization of POM by all oxidants'
       loc_ij_tot(:,:) = &
            & loc_ij_OC(:,:) + loc_ij_ON(:,:) + loc_ij_OP(:,:) + &
            & loc_ij_NC(:,:) + loc_ij_NN(:,:) + loc_ij_NP(:,:) + &
            & loc_ij_SC(:,:) + loc_ij_SN(:,:) + loc_ij_SP(:,:)
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij_tot(:,:),loc_mask_surf)
       ! -------------------------------------------------------- ! (5) percentatges -- denitrification
       loc_unitsname = '% C'
       loc_shortname = 'redox_wcint_pctden_POC'
       loc_longname  = '% of total water-column POC remineralized via denitrification' 
       loc_ij_tot(:,:) = loc_ij_OC(:,:) + loc_ij_NC(:,:) + loc_ij_SC(:,:)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             If (loc_ij_tot(i,j) > const_real_nullsmall) then
                loc_ij(i,j) = 100.0*loc_ij_NC(i,j)/loc_ij_tot(i,j)
             end if
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       loc_unitsname = '% C+N+P'
       loc_shortname = 'redox_wcint_pctden_POM'
       loc_longname  = '% of total water-column POM remineralized via denitrification' 
       loc_ij_tot(:,:) = &
            & loc_ij_OC(:,:) + loc_ij_ON(:,:) + loc_ij_OP(:,:) + &
            & loc_ij_NC(:,:) + loc_ij_NN(:,:) + loc_ij_NP(:,:) + &
            & loc_ij_SC(:,:) + loc_ij_SN(:,:) + loc_ij_SP(:,:)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             If (loc_ij_tot(i,j) > const_real_nullsmall) then
                loc_ij(i,j) = 100.0*(loc_ij_NC(i,j)+loc_ij_NN(i,j)+loc_ij_NP(i,j))/loc_ij_tot(i,j)
             end if
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! -------------------------------------------------------- ! (5) percentatges -- sulphate reduction
       loc_unitsname = '% C'
       loc_shortname = 'redox_wcint_pctsul_POC'
       loc_longname  = '% of total water-column POC remineralized via sulphate reduction' 
       loc_ij_tot(:,:) = loc_ij_OC(:,:) + loc_ij_NC(:,:) + loc_ij_SC(:,:)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             If (loc_ij_tot(i,j) > const_real_nullsmall) then
                loc_ij(i,j) = 100.0*loc_ij_SC(i,j)/loc_ij_tot(i,j)
             end if
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       loc_unitsname = '% C+N+P'
       loc_shortname = 'redox_wcint_pctsul_POM'
       loc_longname  = '% of total water-column POM remineralized via sulphate reduction' 
       loc_ij_tot(:,:) = &
            & loc_ij_OC(:,:) + loc_ij_ON(:,:) + loc_ij_OP(:,:) + &
            & loc_ij_NC(:,:) + loc_ij_NN(:,:) + loc_ij_NP(:,:) + &
            & loc_ij_SC(:,:) + loc_ij_SN(:,:) + loc_ij_SP(:,:)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             If (loc_ij_tot(i,j) > const_real_nullsmall) then
                loc_ij(i,j) = 100.0*(loc_ij_SC(i,j)+loc_ij_SN(i,j)+loc_ij_SP(i,j))/loc_ij_tot(i,j)
             end if
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,''//trim(loc_shortname), &
            & trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end If
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_geochemistry_basic
  ! ****************************************************************************************************************************** !
  
  
  ! ****************************************************************************************************************************** !
  ! *** save geochemical composition data ***
  SUBROUTINE sub_2d_save_geochemistry_advanced()  
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j
    integer::id,ic
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask_surf
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0   
    loc_mask_surf(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    ! ---------------------------------------------------------------- !
    ! ocean surface carbonate chemistry data
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       DO ic=1,n_carb
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   loc_ij(i,j) = int_carb_timeslice(ic,i,j,n_k)/int_t_timeslice
                end IF
             end DO
          end DO
          SELECT CASE (ic)
          CASE (ic_pHsws)
             loc_unitsname = 'pH(sws)'           
          CASE (ic_ohm_cal,ic_ohm_arg,ic_RF0,ic_RdDICdALK,ic_RdfCO2dDIC,ic_RdALKdDIC,ic_pH_n,ic_err)  
             loc_unitsname = 'n/a' 
          case default
             loc_unitsname = 'mol kg-1'
          end SELECT
          SELECT CASE (ic)
          CASE (ic_RdDICdALK,ic_RdfCO2dDIC,ic_RdALKdDIC,ic_pH_n,ic_err)  
             call sub_adddef_netcdf(loc_iou,3,'carbchem_seasur_'//trim(string_carb(ic)), &
                  & 'Ocean surface carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('carbchem_seasur_'//trim(string_carb(ic)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
          end SELECT
       END DO
       ! save accumulated carbchem error occurrence array
       loc_ij(:,:) = diag_carb_errsum(:,:,n_k)
       loc_shortname = 'carbchem_seasur_errsum'
       loc_longname  = 'Ocean surface carbonate chemistry -- '//'Summed occurrence of failure of the pH calculation.'
       loc_unitsname = 'n/a'
       call sub_adddef_netcdf(loc_iou,3,trim(loc_shortname),trim(loc_longname), &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! seafloor carbonate chemistry
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       if (ctrl_save_hidden_seafloor) then
          DO ic=1,n_carb
             loc_ij(:,:) = const_real_zero
             DO i=1,n_i
                DO j=1,n_j
                   loc_k1 = goldstein_k1(i,j)
                   IF (n_k >= loc_k1) THEN
                      loc_ij(i,j) = int_carb_timeslice(ic,i,j,loc_k1)/int_t_timeslice
                   end IF
                end DO
             end DO
             SELECT CASE (ic)
             CASE (ic_pHsws)
                loc_unitsname = 'pH(sws)'                
             CASE (ic_ohm_cal,ic_ohm_arg,ic_pH_n,ic_err)
                loc_unitsname = 'n/a'
             case default
                loc_unitsname = 'mol kg-1'
             end SELECT
             SELECT CASE (ic)
             CASE (ic_conc_CO2,ic_conc_CO3,ic_conc_HCO3,ic_pH_n,ic_err)
                call sub_adddef_netcdf(loc_iou,3,'carbchem_seafloor_'//trim(string_carb(ic)), &
                     & 'Benthic carbonate chemistry -- '//trim(string_longname_carb(ic)), &
                     & trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar2d('carbchem_seafloor_'//trim(string_carb(ic)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
             end select
          END DO
       end if
    end if
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED PRODUCTION RATE
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol m-2 yr-1'
    DO id=1,n_diag_precip
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = &
                     & sum(phys_ocn(ipo_M,i,j,:)*int_diag_precip_timeslice(id,i,j,:))* &
                     & phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
             end If
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,'geochem_wcint_'//trim(string_diag_precip(id)), &
            & 'water-column integrated production rate - '//trim(string_diag_precip(id)), &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('geochem_wcint_'//trim(string_diag_precip(id)),loc_iou, &
            & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end DO
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED REACTION RATE
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol m-2 yr-1'
    DO id=1,n_diag_react
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = &
                     & sum(phys_ocn(ipo_M,i,j,:)*int_diag_react_timeslice(id,i,j,:))* &
                     & phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
             end If
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,'geochem_wcint_'//trim(string_diag_react(id)), &
            & 'water-column integrated reaction rate - '//trim(string_diag_react(id)), &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('geochem_wcint_'//trim(string_diag_react(id)),loc_iou, &
            & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end DO
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED REDOX TRANSFORMATIONS
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol m-2 yr-1'
    DO id=1,n_diag_redox
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = &
                     & sum(phys_ocn(ipo_M,i,j,:)*int_diag_redox_timeslice(id,i,j,:))* &
                     & phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
             end If
          end DO
       end DO
       call sub_adddef_netcdf(loc_iou,3,'redox_wcint_'//trim(string_diag_redox(id)), &
            & 'water-column integrated redox transformation rate - '//trim(string_diag_redox(id)), &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('redox_wcint_'//trim(string_diag_redox(id)),loc_iou, &
            & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end DO
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_geochemistry_advanced
  ! ****************************************************************************************************************************** !

    
  ! ****************************************************************************************************************************** !
  ! *** save nutrient data ***
  SUBROUTINE sub_2d_save_biologicalpump_basic()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,io,is,ib,m
    integer::loc_k1
    integer::loc_m,loc_tot_m
    integer::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask,loc_mask_surf
    real::loc_tot,loc_frac,loc_standard
    real,DIMENSION(n_sed,n_i,n_j)::loc_isij
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0    
    loc_mask_surf(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES -- surface export -- as flux density
    ! ---------------------------------------------------------------- !
    ! NOTE: bio_settle is in units of mol per time interval (so divide by grid cell area for per m2)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged)
                   loc_ij(i,j) = int_bio_settle_timeslice(is,i,j,n_k)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   loc_unitsname = 'mol m-2 yr-1'
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,n_k)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   loc_frac = int_bio_settle_timeslice(is,i,j,n_k)*phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
                   loc_standard = const_standards(sed_type(is))
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                   loc_unitsname = 'o/oo'
                end SELECT
             end if
          end do
       end do
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged)
          ! bulk tracers
          If (ctrl_save_basic_reservoirs) then
             call sub_adddef_netcdf(loc_iou,3,'biop_seasur_f_'//trim(string_sed(is)), &
                  & 'particulate biological export flux density - '//trim(string_sed(is)), &
                  & loc_unitsname,const_real_zero,const_real_zero)
             call sub_putvar2d('biop_seasur_f_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
          end if
       case (n_itype_min:n_itype_max)
          ! isotopes
          If (ctrl_save_basic_proxies) then
             call sub_adddef_netcdf(loc_iou,3,'biop_seasur_f_'//trim(string_sed(is)), &
                  & 'particulate biological export mean isotopic composition - '//trim(string_sed(is)), &
                  & loc_unitsname,const_real_zero,const_real_zero)
             call sub_putvar2d('biop_seasur_f_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
          end if
       end SELECT
    end do
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES -- seafloor -- flux density
    ! ---------------------------------------------------------------- !
    ! NOTE: bio_settle is in units of mol per time interval (so divide by grid cell area for per m2)
    If (ctrl_save_hidden_seafloor) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged)
                      loc_ij(i,j) = int_bio_settle_timeslice(is,i,j,loc_k1)*phys_ocn(ipo_rA,i,j,loc_k1)/int_t_timeslice
                      loc_unitsname = 'mol m-2 yr-1'
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_bio_settle_timeslice(sed_dep(is),i,j,loc_k1)*phys_ocn(ipo_rA,i,j,loc_k1)/int_t_timeslice
                      loc_frac = int_bio_settle_timeslice(is,i,j,loc_k1)*phys_ocn(ipo_rA,i,j,loc_k1)/int_t_timeslice
                      loc_standard = const_standards(sed_type(is))
                      loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                      loc_unitsname = 'o/oo'
                   end SELECT
                end if
             end do
          end do
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_scavenged)
             ! bulk tracers
             call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_f_'//trim(string_sed(is)), &
                  & 'particulate sediment rain flux density - '//trim(string_sed(is)), &
                  & loc_unitsname,const_real_zero,const_real_zero)
             call sub_putvar2d('biop_seafloor_f_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
          case (n_itype_min:n_itype_max)
             ! isotopes
             If (ctrl_save_basic_proxies) then
                call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_f_'//trim(string_sed(is)), &
                     & 'particulate sediment rain mean isotopic composition - '//trim(string_sed(is)), &
                     & loc_unitsname,const_real_zero,const_real_zero)
                call sub_putvar2d('biop_seafloor_f_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
             end if
          end SELECT
       end do
    end if
    ! ---------------------------------------------------------------- !
    ! CaCO3/POC 'rain ratio' -- seasurface
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'n/a'
    IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POC,i,j,n_k) > const_real_nullsmall) then
                   loc_ij(i,j) = int_bio_settle_timeslice(is_CaCO3,i,j,n_k)/int_bio_settle_timeslice(is_POC,i,j,n_k)
                else
                   loc_ij(i,j) = 0.0
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_CaCO3toPOC','CaCO3/POC surface ocean export rain ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_CaCO3toPOC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! CaCO3/POC 'rain ratio' -- seafloor
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_seafloor) then
       loc_unitsname = 'n/a'
       IF (sed_select(is_CaCO3) .AND. sed_select(is_POC)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   if (int_bio_settle_timeslice(is_POC,i,j,loc_k1) > const_real_nullsmall) then
                      loc_ij(i,j) = int_bio_settle_timeslice(is_CaCO3,i,j,loc_k1)/int_bio_settle_timeslice(is_POC,i,j,loc_k1)
                   else
                      loc_ij(i,j) = 0.0
                   end if
                end IF
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_r_CaCO3toPOC','CaCO3/POC seafloor rain ratio', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('biop_seafloor_r_CaCO3toPOC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end if
    end if
    ! ---------------------------------------------------------------- !
    ! opal/POC surface ocean export 'rain ratio'
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'n/a'
    IF (sed_select(is_opal) .AND. sed_select(is_POC)) THEN
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POC,i,j,n_k) > const_real_nullsmall) then
                   loc_ij(i,j) = int_bio_settle_timeslice(is_opal,i,j,n_k)/int_bio_settle_timeslice(is_POC,i,j,n_k)
                else
                   loc_ij(i,j) = 0.0
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_opaltoPOC','opal/POC surface ocean export rain ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_opaltoPOC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                if (int_bio_settle_timeslice(is_POC,i,j,loc_k1) > const_real_nullsmall) then
                   loc_ij(i,j) = int_bio_settle_timeslice(is_opal,i,j,loc_k1)/int_bio_settle_timeslice(is_POC,i,j,loc_k1)
                else
                   loc_ij(i,j) = 0.0
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_r_opaltoPOC','opal/POC seafloor rain ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seafloor_r_opaltoPOC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! C/P export cellular quotient ratio 
    ! ---------------------------------------------------------------- !  
    IF (sed_select(is_POP) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'n/a'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POP,i,j,n_k) > const_rns) then
                   loc_ij(i,j) = int_bio_settle_timeslice(is_POC,i,j,n_k)/int_bio_settle_timeslice(is_POP,i,j,n_k)
                else
                   loc_ij(i,j) = 0.0
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_CtoP','C/P ratio of surface ocean POM export', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_CtoP',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! C/Fe export cellular quotient ratio
    ! ---------------------------------------------------------------- !
    IF (sed_select(is_POFe) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'n/a'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POFe,i,j,n_k) > const_real_nullsmall) then
                   loc_ij(i,j) = int_bio_settle_timeslice(is_POC,i,j,n_k)/int_bio_settle_timeslice(is_POFe,i,j,n_k)
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_CtoFe','C/Fe ratio of surface ocean POM export', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_CtoFe',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end IF
    ! ---------------------------------------------------------------- !
    ! nutrient availablity diagnostics
    ! ---------------------------------------------------------------- !
    if (ocn_select(io_PO4) .AND. ocn_select(io_SiO2)) then
       loc_unitsname = 'n/a'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                loc_ij(i,j) = int_ocn_timeslice(io_SiO2,i,j,n_k)/int_t_timeslice - &
                     & par_bio_red_POP_PON*int_ocn_timeslice(io_PO4,i,j,n_k)/int_t_timeslice
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_SiSTAR','Si Star (using PO4 and assumed N:P)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_SiSTAR',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! Fe diagnostics
    ! ---------------------------------------------------------------- !
    If (ctrl_data_save_slice_sur .AND. (ocn_select(io_Fe) .OR. ocn_select(io_TDFe))) then
       ! total aeolian Fe flux
       loc_unitsname = 'mg Fe m-2 yr-1'
       loc_ij(:,:) = conv_mol_mmol*par_det_Fe_frac*conv_det_mol_g* &
            & (int_phys_ocn_timeslice(ipo_rA,:,:,n_k)*int_bio_settle_timeslice(is_det,:,:,n_k))/(int_t_timeslice**2)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fluxden','Total aeolian iron flux density to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fluxden',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! solubility (%)
       loc_unitsname = '%'
       loc_ij(:,:) = 100.0*int_phys_ocnatm_timeslice(ipoa_solFe,:,:)/int_t_timeslice
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_sol','Aeolian iron solubility', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_sol',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! nutrient limitation
    ! ---------------------------------------------------------------- !
    ! NOTE: -1.0 equates to dominance of PO4 limitation, +1.0 to dominance of Fe limitation
    !       a value of ~0.0 represents ~equal limitation
    ! for reference: loc_kPO4 = loc_PO4/(loc_PO4 + par_bio_c0_PO4), loc_kFe = loc_FeT/(loc_FeT + par_bio_c0_Fe)
    !                i.e. 1.0 == no limitation, 0.0 == complete limitation
    if ( ocn_select(io_PO4) .AND. (ocn_select(io_Fe) .OR. ocn_select(io_TDFe)) ) then
       loc_unitsname = 'n/a'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                loc_ij(i,j) = &
                     & (int_diag_bio_timeslice(idiag_bio_kPO4,i,j) - int_diag_bio_timeslice(idiag_bio_kFe,i,j)) &
                     & / &
                     & (int_diag_bio_timeslice(idiag_bio_kPO4,i,j) + int_diag_bio_timeslice(idiag_bio_kFe,i,j))
             end if
          end do
       end do
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_FevsPO4lim','occurrence of Fe (+1.0) vs. PO4 (-1.0) limitation', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_FevsPO4lim',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! ECOGEM diagnostics
    ! ---------------------------------------------------------------- !
    if (flag_ecogem) then
       ! save POM
       DO l=1,n_l_sed
          loc_ij(:,:) = 0.0
          is = conv_iselected_is(l)
          if ( (is == is_POC) .OR. (sed_type(is) == par_sed_type_POM)) then
             loc_ij(:,:) = int_diag_ecogem_part(is,:,:)/int_t_timeslice
             call sub_adddef_netcdf(loc_iou,3,'eco_f_POC_'//trim(string_sed(is)), &
                  & 'ECOGEM particulate organic matter production - '//trim(string_sed(is)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('eco_f_POC_'//trim(string_sed(is)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_isij(is,:,:),loc_mask_surf)
          end if
       end do
       ! save DOM
       ! NOTE: only save the io tracers that convert (non zero) to POM (i.e. DOM)
       DO l=3,n_l_ocn
          loc_ij(:,:) = 0.0
          io = conv_iselected_io(l)
          loc_tot_m = conv_DOM_POM_i(0,io)
          do loc_m=1,loc_tot_m
             loc_ij(:,:) = int_diag_ecogem_remin(io,:,:)/int_t_timeslice
             call sub_adddef_netcdf(loc_iou,3,'eco_f_DOC_'//trim(string_ocn(io)), &
                  & 'ECOGEM dissolved matter production - '//trim(string_ocn(io)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('eco_f_DOC_'//trim(string_ocn(io)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
          end do
       end DO
       ! calculate POM equivalnt of DOM
       loc_isij(:,:,:) = 0.0
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                DO l=3,n_l_ocn
                   io = conv_iselected_io(l)
                   loc_tot_m = conv_DOM_POM_i(0,io)
                   do loc_m=1,loc_tot_m
                      is = conv_DOM_POM_i(loc_m,io)
                      loc_isij(is,i,j) = loc_isij(is,i,j) + conv_DOM_POM(is,io)*int_diag_ecogem_remin(io,i,j)
                   end do
                end do
             end If
          end DO
       end DO
       ! calculate DOM ratio (replace values in same local array) and save as netCDF
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          if ( (is == is_POC) .OR. (sed_type(is) == par_sed_type_POM)) then
             DO i=1,n_i
                DO j=1,n_j
                   If (goldstein_k1(i,j) <= n_k) then
                      if ((loc_isij(is,i,j)+int_diag_ecogem_part(is,i,j)) > const_real_nullsmall) then
                         loc_isij(is,i,j) = loc_isij(is,i,j)/(loc_isij(is,i,j)+int_diag_ecogem_part(is,i,j))
                      else
                         loc_isij(is,i,j) = 0.0
                      end if
                   end If
                end DO
             end DO
             call sub_adddef_netcdf(loc_iou,3,'eco_r_DOMfract_'//trim(string_sed(is)), &
                  & 'ECOGEM dissolved matter production fraction - '//trim(string_sed(is)), &
                  & trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('eco_r_DOMfract_'//trim(string_sed(is)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_isij(is,:,:),loc_mask_surf)
          end if
       end do
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_biologicalpump_basic
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save nutrient data ***
  SUBROUTINE sub_2d_save_biologicalpump_advanced()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,io,is,ib
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask,loc_mask_surf,loc_sed_mask
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0   
    loc_mask_surf(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_sed_mask(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES -- surface export -- total flux
    ! ---------------------------------------------------------------- !
    ! NOTE: bio_settle is in units of mol per time interval
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                SELECT CASE (sed_type(is))
                CASE (par_sed_type_bio,par_sed_type_abio, &
                     & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                     & par_sed_type_scavenged)
                   loc_ij(i,j) = int_bio_settle_timeslice(is,i,j,n_k)/int_t_timeslice
                end SELECT
             end if
          end do
       end do
       loc_unitsname = 'mol yr-1'
       SELECT CASE (sed_type(is))
       CASE (par_sed_type_bio,par_sed_type_abio, &
            & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
            & par_sed_type_scavenged)
          loc_unitsname = 'mol yr-1'
          call sub_adddef_netcdf(loc_iou,3,'biop_seasur_ftot_'//trim(string_sed(is)), &
               & 'total biological particulate export flux - '//trim(string_sed(is)), &
               & loc_unitsname,const_real_zero,const_real_zero)
          call sub_putvar2d('biop_seasur_ftot_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
       end SELECT
    end do
    ! ---------------------------------------------------------------- !
    ! PARTICULATE FLUXES -- seafloor -- total flux
    ! ---------------------------------------------------------------- !
    ! NOTE: bio_settle is in units of mol per time interval
    if (ctrl_save_hidden_seafloor) then
       DO l=1,n_l_sed
          is = conv_iselected_is(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_abio, &
                        & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                        & par_sed_type_scavenged)
                      loc_ij(i,j) = int_bio_settle_timeslice(is,i,j,loc_k1)/int_t_timeslice
                   end SELECT
                end if
             end do
          end do
          SELECT CASE (sed_type(is))
          CASE (par_sed_type_bio,par_sed_type_abio, &
               & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
               & par_sed_type_scavenged)
             loc_unitsname = 'mol yr-1'
             If (ctrl_save_advanced_biologicalpump) then
                call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_ftot_'//trim(string_sed(is)), &
                     & 'total particulate rain flux to the sediments - '//trim(string_sed(is)), &
                     & loc_unitsname,const_real_zero,const_real_zero)
                call sub_putvar2d('biop_seafloor_ftot_'//trim(string_sed(is)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
             end if
          end SELECT
       end do
    end if
    ! ---------------------------------------------------------------- !
    ! particulate flux fractions -- seafloor
    ! ---------------------------------------------------------------- !
    If (ctrl_save_hidden_seafloor) then
       loc_unitsname = 'n/a'
       IF (sed_select(is_POC_frac2)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   loc_ij(i,j) = int_bio_settle_timeslice(is_POC_frac2,i,j,loc_k1)/real(int_t_timeslice_count)
                END if
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_r_POCfrac2', &
               & 'POC fraction #2',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('biop_seafloor_r_POCfrac2',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_sed_mask)
       end IF
       IF (sed_select(is_CaCO3_frac2)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   loc_ij(i,j) = int_bio_settle_timeslice(is_CaCO3_frac2,i,j,loc_k1)/real(int_t_timeslice_count)
                END if
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_r_CaCO3frac2', &
               & 'CaCO3 fraction #2',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('biop_seafloor_r_CaCO3frac2',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_sed_mask)
       end IF
       IF (sed_select(is_opal_frac2)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   loc_ij(i,j) = int_bio_settle_timeslice(is_opal_frac2,i,j,loc_k1)/real(int_t_timeslice_count)
                END if
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'biop_seafloor_r_opalfrac2', &
               & 'opal fraction #2',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('biop_seafloor_r_opalfrac2',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_sed_mask)
       end IF
    end IF
    ! ---------------------------------------------------------------- !
    ! P/C export cellular quotient ratio
    ! ---------------------------------------------------------------- !
    IF (sed_select(is_POP) .AND. sed_select(is_POC)) THEN
       loc_unitsname = 'o/oo'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POC,i,j,n_k) > const_rns) then
                   loc_ij(i,j) = 1.0E3*int_bio_settle_timeslice(is_POP,i,j,n_k)/int_bio_settle_timeslice(is_POC,i,j,n_k)
                else
                   loc_ij(i,j) = const_real_null
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_PtoC','P/C ratio of ocean surface POM export (in units of per mil)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_PtoC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end IF
    ! ---------------------------------------------------------------- !
    ! Fe/C export cellular quotient ratio
    ! ---------------------------------------------------------------- !
    IF (sed_select(is_POFe) .AND. sed_select(is_POC)) THEN
       loc_unitsname = '10^3 o/oo'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             IF (n_k >= goldstein_k1(i,j)) THEN
                if (int_bio_settle_timeslice(is_POC,i,j,n_k) > const_real_nullsmall) then
                   loc_ij(i,j) = 1.0E6*int_bio_settle_timeslice(is_POFe,i,j,n_k)/int_bio_settle_timeslice(is_POC,i,j,n_k)
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'biop_seasur_r_FetoC','average POM export Fe/C cellular ratio', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_seasur_r_FetoC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end IF
    ! ---------------------------------------------------------------- !
    ! Fe diagnostics
    ! ---------------------------------------------------------------- !
    If (ocn_select(io_Fe) .OR. ocn_select(io_TDFe)) then
       ! total aeolian Fe flux (mass)
       ! NOTE: new calculation of loc_ij(:,:)
       loc_unitsname = 'mg Fe m-2 yr-1'
       loc_ij(:,:) = conv_mol_mmol*par_det_Fe_frac*conv_det_mol_g* &
            & (int_phys_ocn_timeslice(ipo_rA,:,:,n_k)*int_bio_settle_timeslice(is_det,:,:,n_k))/(int_t_timeslice**2)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fFetot_g','Total aeolian iron flux density to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fFetot_g',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! total aeolian Fe flux (moles)
       ! NOTE: based on preceeding step calculation of loc_ij(:,:)
       loc_unitsname = 'mmol Fe m-2 yr-1'
       loc_ij(:,:) = conv_Fe_g_mol*loc_ij(:,:)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fFetot_mol','Total aeolian iron flux density to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fFetot_mol',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! solulablized aeolian Fe flux
       ! NOTE: based on preceeding step calculation of loc_ij(:,:)
       loc_unitsname = 'umol Fe m-2 yr-1'
       loc_ij(:,:) = conv_mol_umol*conv_mmol_mol*(int_phys_ocnatm_timeslice(ipoa_solFe,:,:)/int_t_timeslice)*loc_ij(:,:)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fFe_mol','Solulablized aeolian iron flux density to surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fFe_mol',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! solulablized aeolian Fe flux (2)
       ! NOTE: based on preceeding step calculation of loc_ij(:,:)
       loc_unitsname = 'mol Fe yr-1'
       loc_ij(:,:) = conv_umol_mol*(int_phys_ocn_timeslice(ipo_A,:,:,n_k)/int_t_timeslice)*loc_ij(:,:)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fFe','Solulablized aeolian iron flux to surface grid points', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fFe',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! particulate Fe loss
       ! NOTE: new calculation of loc_ij(:,:)
       loc_unitsname = 'umol Fe m-2 yr-1'
       loc_ij(:,:) = conv_mol_umol* &
            & int_phys_ocn_timeslice(ipo_rA,:,:,n_k)*int_bio_settle_timeslice(is_POFe,:,:,n_k)/(int_t_timeslice**2)
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fpartFe','Particulate organic matter iron loss from surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fpartFe',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! total scavenged Fe loss
       loc_unitsname = 'umol Fe m-2 yr-1'
       loc_ij(:,:) = (conv_mol_umol*int_phys_ocn_timeslice(ipo_rA,:,:,n_k)/(int_t_timeslice**2))* &
            & ( &
            &   int_bio_settle_timeslice(is_POM_Fe,:,:,n_k) &
            & )
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fscavFetot','Total scavenged Fe loss from surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fscavFetot',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       ! total scavenged FeOOH loss
       loc_unitsname = 'umol Fe m-2 yr-1'
       loc_ij(:,:) = (conv_mol_umol*int_phys_ocn_timeslice(ipo_rA,:,:,n_k)/(int_t_timeslice**2))* &
            & ( &
            &   int_bio_settle_timeslice(is_POM_FeOOH,:,:,n_k)   + &
            &   int_bio_settle_timeslice(is_CaCO3_FeOOH,:,:,n_k) + &
            &   int_bio_settle_timeslice(is_opal_FeOOH,:,:,n_k)  + &
            &   int_bio_settle_timeslice(is_det_FeOOH,:,:,n_k)     &
            & )
       call sub_adddef_netcdf(loc_iou,3,'biop_iron_fscavFeOOHtot','Total scavenged FeOOH loss from surface', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_iron_fscavFeOOHtot',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! ---------------------------------------------------------------- !
    ! Biological productivity controls
    ! ---------------------------------------------------------------- !
    DO ib=1,n_diag_bio
       select case (ib)
       CASE (idiag_bio_dPO4,idiag_bio_dPO4_1,idiag_bio_dPO4_2)
          loc_unitsname = 'mol kg-1 yr-1'
       CASE (idiag_bio_N2fixation,idiag_bio_NH4assim)
          loc_unitsname = 'mol kg-1 yr-1'
       CASE (idiag_bio_DOMlifetime)
          loc_unitsname = 'yr'
       case default
          loc_unitsname = 'n/a'
       end select
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = int_diag_bio_timeslice(ib,i,j)/int_t_timeslice
             end If
          end DO
       end DO
       select case (ib)
       CASE (idiag_bio_CaCO3toPOC_nsp,idiag_bio_opaltoPOC_sp,idiag_bio_fspPOC)
          ! correct for the number of sub-slices to create an average
          loc_ij(:,:) = loc_ij(:,:)/real(int_t_timeslice_count)
       end select
       call sub_adddef_netcdf(loc_iou,3,'biop_bio_'//trim(string_diag_bio(ib)), &
            & 'biological productivity control - '//trim(string_diag_bio(ib)), &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('biop_bio_'//trim(string_diag_bio(ib)),loc_iou, &
            & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end DO
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_biologicalpump_advanced
  ! ****************************************************************************************************************************** !
 
  
  ! ****************************************************************************************************************************** !
  ! *** save proxy data ***
  SUBROUTINE sub_2d_save_proxies_basic()
    ! ---------------------------------------------------------------- !
    ! definelocal variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,io
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask,loc_mask_sur
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:)           = 0.0
    loc_mask_sur(:,:)     = phys_ocnatm(ipoa_mask_ocn,:,:)
    ! ---------------------------------------------------------------- !
    ! I/Ca
    ! ---------------------------------------------------------------- !
    If (ocn_select(io_IO3)) then
       ! I/Ca proxy model
       ! NOTE: from 'I/Ca evidence for upper ocean deoxygenation during the PETM' [Zhou te al., 2014]; 10.1002/2014PA002702!
       ! KD is calculated as [I/Ca]/[IO3-], with I/Ca == umol/mol and [IO3] == umol l-1
       ! KD = -0.16T + 13.65 where T is in oC
       ! => converting to mol/mol and umol kg-1:
       ! I/Ca = (1.0E-3*conv_m3_kg)*[IO3-] * KD
       ! (1) surface ocean I/Ca
       loc_ij(:,:) = const_real_zero
       loc_unitsname = 'mol/mol'
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                loc_ij(i,j) =(1.0E-3*conv_m3_kg)*int_ocn_timeslice(io_IO3,i,j,n_k)/int_t_timeslice * &
                     & 1.0E6*(-0.16*(int_ocn_timeslice(io_T,i,j,n_k)/int_t_timeslice - const_zeroC) + 13.65)                
             END if
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'proxy_seasur_ICa','ocean surface I/Ca', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('proxy_seasur_ICa',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_sur)
       ! (2) benthic I/Ca
       loc_ij(:,:) = const_real_zero
       loc_unitsname = 'mol/mol'
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                loc_ij(i,j) =(1.0E-3*conv_m3_kg)*int_ocn_timeslice(io_IO3,i,j,loc_k1)/int_t_timeslice * &
                     & 1.0E6*(-0.16*(int_ocn_timeslice(io_T,i,j,loc_k1)/int_t_timeslice - const_zeroC) + 13.65) 
             END if
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'proxy_seafloor_ICa','seafloor I/Ca', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('proxy_seafloor_ICa',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_sur)
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_proxies_basic
  ! ****************************************************************************************************************************** !
 
  
  ! ****************************************************************************************************************************** !
  ! *** save proxy data ***
  SUBROUTINE sub_2d_save_proxies_advanced()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,io,is
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    CHARACTER(len=255)::loc_unitsname
    real,DIMENSION(n_i,n_j)::loc_ij,loc_ij_1,loc_ij_2,loc_mask,loc_mask_surf,loc_sed_mask
    real::loc_tot,loc_frac,loc_standard
    real::loc_tot1,loc_frac1,loc_tot2,loc_frac2
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0   
    loc_mask_surf(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_sed_mask(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k)
    ! ---------------------------------------------------------------- !
    ! save planktic-benthic difference
    ! ---------------------------------------------------------------- !
    ! NOTE: exclude dissolved organic matter tracers
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
       if (is == 0) then
          loc_ij(:,:) = const_real_zero
          loc_ij_2(:,:) = const_real_zero
          loc_ij_2(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                loc_k1 = goldstein_k1(i,j)
                IF (n_k >= loc_k1) THEN
                   SELECT CASE (ocn_type(io))
                   CASE (0)
                      loc_ij_1(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice
                      loc_ij_2(i,j) = int_ocn_timeslice(io,i,j,loc_k1)/int_t_timeslice
                   CASE (1)
                      loc_ij_1(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice
                      loc_ij_2(i,j) = int_ocn_timeslice(io,i,j,loc_k1)/int_t_timeslice
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,n_k)
                      loc_frac = int_ocn_timeslice(io,i,j,n_k)
                      loc_standard = const_standards(ocn_type(io))
                      loc_ij_1(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.,const_real_null)
                      loc_tot  = int_ocn_timeslice(ocn_dep(io),i,j,loc_k1)
                      loc_frac = int_ocn_timeslice(io,i,j,loc_k1)
                      loc_standard = const_standards(ocn_type(io))
                      loc_ij_2(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.,const_real_null)
                   end SELECT
                end IF
             end DO
          end DO
          loc_ij(:,:) = loc_ij_1(:,:) - loc_ij_2(:,:)
          SELECT CASE (ocn_type(io))
          CASE (0)
             If (io == io_T) loc_unitsname = 'degrees C'
             If (io == io_S) loc_unitsname = 'o/oo'
          CASE (1)
             loc_unitsname = 'mol kg-1'
          case (n_itype_min:n_itype_max)
             loc_unitsname = 'o/oo'
          end SELECT
          SELECT CASE (ocn_type(io))
          CASE (0,1,n_itype_min:n_itype_max)
             call sub_adddef_netcdf(loc_iou, 3,'proxy_D_'//trim(string_ocn(io)), &
                  & 'planktic-benthic difference '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('proxy_D_'//trim(string_ocn(io)),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_sed_mask)
          end SELECT
       end if
    END DO
    ! ---------------------------------------------------------------- !
    ! Cd particulate surface ocean export trace metal ratios
    ! ---------------------------------------------------------------- !
    IF (ocn_select(io_Cd)) THEN
       loc_unitsname = 'nmol kg-1 (umol kg-1)-1'
       IF (sed_select(is_POCd) .AND. sed_select(is_POC)) THEN
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   if (int_bio_settle_timeslice(is_POC,i,j,n_k) > const_real_nullsmall) then
                      loc_ij(i,j) = 1.0E3*int_bio_settle_timeslice(is_POCd,i,j,n_k)/int_bio_settle_timeslice(is_POC,i,j,n_k)
                   end if
                end IF
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'proxy_seasur_r_CdtoPOC','Cd to C organic matter export ratio', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('proxy_seasur_r_CdtoPOC',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end IF
       IF (sed_select(is_POCd) .AND. sed_select(is_POP)) THEN
          loc_unitsname = 'nmol kg-1 (umol kg-1)-1'
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   if (int_bio_settle_timeslice(is_POP,i,j,n_k) > const_real_nullsmall) then
                      loc_ij(i,j) = 1.0E3*int_bio_settle_timeslice(is_POCd,i,j,n_k)/int_bio_settle_timeslice(is_POP,i,j,n_k)
                   end if
                end IF
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'proxy_seasur_r_CdtoPOP','Cd to P organic matter export ratio', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('proxy_seasur_r_CdtoPOP',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end IF
       IF (sed_select(is_CdCO3) .AND. sed_select(is_CaCO3)) THEN
          loc_unitsname = 'nmol kg-1 (mmol kg-1)-1'
          loc_ij(:,:) = const_real_null
          DO i=1,n_i
             DO j=1,n_j
                IF (n_k >= goldstein_k1(i,j)) THEN
                   if (int_bio_settle_timeslice(is_CaCO3,i,j,n_k) > const_real_nullsmall) then
                      loc_ij(i,j) = 1.0E6*int_bio_settle_timeslice(is_CdCO3,i,j,n_k)/int_bio_settle_timeslice(is_CaCO3,i,j,n_k)
                   end if
                end IF
             END DO
          END DO
          call sub_adddef_netcdf(loc_iou,3,'proxy_seasur_r_CdtoCaCO3','Cd:Ca trace metal ratio (carbonate)', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('proxy_seasur_r_CdtoCaCO3',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
       end IF
    end IF
    ! ---------------------------------------------------------------- !
    ! save overlying Cd trace metal ratios
    ! ---------------------------------------------------------------- !
    IF (ocn_select(io_Cd) .AND. ocn_select(io_Ca)) THEN
       loc_unitsname = 'nmol kg-1 (mmol kg-1)-1'
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                if (int_ocn_timeslice(io_Ca,i,j,loc_k1) > const_real_nullsmall) then
                   loc_ij(i,j) = 1.0E6*int_ocn_timeslice(io_Cd,i,j,loc_k1)/int_ocn_timeslice(io_Ca,i,j,loc_k1)
                end if
             end IF
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'proxy_seasur_r_CdtoCa','Cd:Ca trace metal ratio (ocean)', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('proxy_seasur_r_CdtoCa',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_sed_mask)
    end IF
    ! ---------------------------------------------------------------- !
    ! MISC proxies
    ! ---------------------------------------------------------------- !
    ! DpH
    if (ocn_select(io_DIC) .AND. ocn_select(io_ALK)) then
       loc_unitsname = 'pH units (SWS)'
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                loc_ij(i,j) = (-LOG10(int_carb_timeslice(ic_H,i,j,n_k)/int_t_timeslice)) - &
                     & (-LOG10(int_carb_timeslice(ic_H,i,j,loc_k1)/int_t_timeslice))
             END if
          END DO
       END DO
       call sub_adddef_netcdf(loc_iou,3,'proxy_D_pH','surface-benthic DpH', &
            & trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('proxy_D_pH',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    end if
    ! Dd13C (DIC)
    loc_unitsname = 'o/oo'
    loc_ij(:,:) = const_real_zero
    DO i=1,n_i
       DO j=1,n_j
          loc_k1 = goldstein_k1(i,j)
          IF (n_k >= loc_k1) THEN
             loc_tot1     = int_ocn_timeslice(io_DIC,i,j,n_k)
             loc_frac1    = int_ocn_timeslice(io_DIC_13C,i,j,n_k)
             loc_tot2     = int_ocn_timeslice(io_DIC,i,j,loc_k1)
             loc_frac2    = int_ocn_timeslice(io_DIC_13C,i,j,loc_k1)
             loc_standard = const_standards(ocn_type(io_DIC_13C))
             loc_ij(i,j) = fun_calc_isotope_delta(loc_tot1,loc_frac1,loc_standard,.FALSE.,const_real_null) - &
                  & fun_calc_isotope_delta(loc_tot2,loc_frac2,loc_standard,.FALSE.,const_real_null)
          END if
       END DO
    END DO
    call sub_adddef_netcdf(loc_iou,3,'proxy_D_d13C','surface-benthic Dd13C', &
         & trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('proxy_D_d13C',loc_iou,n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
    ! benthic d13C of DIC
    loc_ij(:,:) = const_real_zero
    loc_unitsname = 'o/oo'
    DO i=1,n_i
       DO j=1,n_j
          loc_k1 = goldstein_k1(i,j)
          IF (n_k >= loc_k1) THEN
             loc_tot  = int_ocn_timeslice(io_DIC,i,j,loc_k1)
             loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,loc_k1)
             loc_standard = const_standards(ocn_type(io_DIC_13C))
             loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
          end IF
       end DO
    end DO
    call sub_adddef_netcdf(loc_iou,3,'proxy_seafloor_DIC_d13C', &
         & 'bottom-water DIC d13C',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('proxy_seafloor_DIC_d13C',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
    ! benthic d13C of HCO3-
    loc_ij(:,:) = const_real_zero
    loc_unitsname = 'o/oo'
    DO i=1,n_i
       DO j=1,n_j
          loc_k1 = goldstein_k1(i,j)
          IF (n_k >= loc_k1) THEN
             loc_tot  = int_carb_timeslice(ic_conc_HCO3,i,j,loc_k1)
             loc_frac = int_carbisor_timeslice(ici_HCO3_r13C,i,j,loc_k1)*int_carb_timeslice(ic_conc_HCO3,i,j,loc_k1)
             loc_standard = const_standards(ocn_type(io_DIC_13C))
             loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
          end IF
       end DO
    end DO
    call sub_adddef_netcdf(loc_iou,3,'proxy_seafloor_HCO3_d13C', &
         & 'bottom-water HCO3- d13C',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('proxy_seafloor_HCO3_d13C',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
    ! benthic Schmittner d13C LAI
    loc_ij(:,:) = const_real_zero
    loc_unitsname = 'o/oo'
    DO i=1,n_i
       DO j=1,n_j
          loc_k1 = goldstein_k1(i,j)
          IF (n_k >= loc_k1) THEN
             loc_tot  = int_ocn_timeslice(io_DIC,i,j,loc_k1)
             loc_frac = int_ocn_timeslice(io_DIC_13C,i,j,loc_k1)
             loc_standard = const_standards(ocn_type(io_DIC_13C))
             loc_ij(i,j) = const_d13C_LA1_a + &
                  & const_d13C_LA1_b*fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null) + &
                  & const_d13C_LA1_c*1.0E6*int_carb_timeslice(ic_conc_CO3,i,j,loc_k1) + &
                  & const_d13C_LA1_d*int_phys_ocn_timeslice(ipo_Dbot,i,j,loc_k1)/int_t_timeslice
          end IF
       end DO
    end DO
    call sub_adddef_netcdf(loc_iou,3,'proxy_seafloor_LA1', &
         & 'bottom-water Schmittner d13C (LA1)',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('proxy_seafloor_LA1',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_surf)
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_proxies_advanced
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! *** save time-slice data ***
  SUBROUTINE sub_2d_save_hidden_grid()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::ip
    integer::loc_iou,loc_ntrec
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask_surf_ALL
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou                = ncout2d_iou
    loc_ntrec              = ncout2d_ntrec
    loc_mask_surf_ALL(:,:) = 1.0
    ! ---------------------------------------------------------------- !
    ! phys_atm -- ATMOSPHERE PHYSICS AND CLIMATE
    ! ---------------------------------------------------------------- !
    if (ctrl_save_hidden_extra) then
       loc_ij(:,:) = const_real_zero
       DO ip=1,n_phys_ocnatm
          loc_ij(:,:) = int_phys_ocnatm_timeslice(ip,:,:)/int_t_timeslice
          call sub_adddef_netcdf(loc_iou,3,'2Dgrid_'//trim(string_phys_ocnatm(ip)), &
               & 'atmosphere physics and grid - '//trim(string_phys_ocnatm(ip)),' ',const_real_zero,const_real_zero)
          call sub_putvar2d('2Dgrid_'//trim(string_phys_ocnatm(ip)),loc_iou, &
               & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf_ALL)
       END DO
    else
       loc_ij(:,:) = const_real_zero
       DO ip=1,n_phys_ocnatm
          SELECT CASE (ip)
          CASE (ipoa_A) 
             loc_ij(:,:) = int_phys_ocnatm_timeslice(ip,:,:)/int_t_timeslice
             call sub_adddef_netcdf(loc_iou,3,'2Dgrid_'//trim(string_phys_ocnatm(ip)), &
                  & 'atmosphere grid - '//trim(string_phys_ocnatm(ip)),' ',const_real_zero,const_real_zero)
             call sub_putvar2d('2Dgrid_'//trim(string_phys_ocnatm(ip)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf_ALL)
          end SELECT
       END DO
    end if
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_hidden_grid
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save time-slice data ***
  SUBROUTINE sub_2d_save_hidden_climate(dum_dtyr)
    ! ---------------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! ---------------------------------------------------------------- !
    real,intent(in)::dum_dtyr
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    integer::i,j,io,ia,ip,l
    integer::loc_k1
    integer::loc_iou,loc_ntrec
    REAL::loc_scale
    real,DIMENSION(n_i,n_j)::loc_ij
    real,DIMENSION(n_i,n_j)::loc_mask_sur,loc_mask_sur_ALL
    real,DIMENSION(0:n_j,0:n_k)::loc_mask_opsi,loc_tmp_jk
    real,DIMENSION(1:n_i,0:n_j)::loc_mask_psi,loc_tmp_ij
    CHARACTER(len=255)::loc_unitsname,loc_shortname,loc_longname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_mask_sur(:,:)     = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_sur_ALL(:,:) = const_real_one
    loc_tmp_jk(:,:)    = 0.0
    loc_tmp_ij(:,:)    = 0.0
    ! ---------------------------------------------------------------- !
    ! SET LOCAL VARIABLES
    ! ---------------------------------------------------------------- !
    loc_ntrec = ncout2d_ntrec
    loc_iou   = ncout2d_iou
    loc_scale = goldstein_dsc*goldstein_usc*const_rEarth*1.0E-6
    ! ---------------------------------------------------------------- !
    ! phys_atm -- 'CLIMATE'
    ! ---------------------------------------------------------------- !
    ! (2) surface wind speed
    loc_unitsname = 'm s-1'
    call sub_adddef_netcdf(loc_iou,3,'climate_wspeed','windspeed',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('climate_wspeed',loc_iou,n_i,n_j,loc_ntrec, &
         & int_phys_ocnatm_timeslice(ipoa_wspeed,:,:)/int_t_timeslice,loc_mask_sur)
    ! (3) fractional sea-ice cover
    loc_unitsname = 'n/a'
    call sub_adddef_netcdf(loc_iou,3,'climate_seaice','sea-ice cover (%)',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('climate_seaice',loc_iou,n_i,n_j,loc_ntrec, &
         & 100.0*int_phys_ocnatm_timeslice(ipoa_seaice,:,:)/int_t_timeslice,loc_mask_sur)
    ! (4) sea-ice thickness
    loc_unitsname = 'm'
    call sub_adddef_netcdf (loc_iou,3,'climate_seaice_th','sea-ice thickness', &
         & trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('climate_seaice_th',loc_iou,n_i,n_j,loc_ntrec, &
         & int_phys_ocnatm_timeslice(ipoa_seaice_th,:,:)/int_t_timeslice,loc_mask_sur)
    ! (5) incident sw radiation
    loc_unitsname = 'W m-2'
    call sub_adddef_netcdf(loc_iou,3,'climate_fxsw','incident sw radiation',trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('climate_fxsw',loc_iou,n_i,n_j,loc_ntrec, &
         & int_phys_ocnatm_timeslice(ipoa_fxsw,:,:)/int_t_timeslice,loc_mask_sur_ALL)
    ! (6) convective 'cost' (need to un-do the time-step weighting)
    loc_unitsname = 'yr-1'
    loc_longname = 'convective cost (column integrated adjustments per year)'
    call sub_adddef_netcdf(loc_iou,3,'climate_cost',trim(loc_longname),trim(loc_unitsname),const_real_zero,const_real_zero)
    call sub_putvar2d('climate_cost',loc_iou,n_i,n_j,loc_ntrec, &
         & int_phys_ocnatm_timeslice(ipoa_cost,:,:)/int_t_timeslice/dum_dtyr,loc_mask_sur)
    if (ctrl_save_hidden_extra) then
       ! (7) solar forcing
       loc_unitsname = 'W m-2'
       call sub_adddef_netcdf(loc_iou,3,'climate_solfor','solar forcing',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('climate_solfor',loc_iou,n_i,n_j,loc_ntrec, &
            & int_phys_ocnatm_timeslice(ipoa_solfor,:,:)/int_t_timeslice,loc_mask_sur_ALL)
       ! (8) wind stress
       loc_unitsname = 'N/m-2'
       call sub_adddef_netcdf (loc_iou,3,'climate_tau_u','wind stress (u)',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('climate_tau_u',loc_iou,n_i,n_j,loc_ntrec, &
            & int_phys_ocnatm_timeslice(ipoa_tau_u,:,:)/int_t_timeslice,loc_mask_sur)
       loc_unitsname = 'N/m-2'
       call sub_adddef_netcdf(loc_iou,3,'climate_tau_v','wind stress (v)',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('climate_tau_v',loc_iou,n_i,n_j,loc_ntrec, &
            & int_phys_ocnatm_timeslice(ipoa_tau_v,:,:)/int_t_timeslice,loc_mask_sur)
       ! (9) air-sea gas exchange coefficient
       if (atm_select(ia_pCO2)) then
          loc_unitsname = 'mol m-2 yr-1 uatm-1'
          call sub_adddef_netcdf(loc_iou,3,'climate_KCO2','air-sea gasx coef.',trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('climate_KCO2',loc_iou,n_i,n_j,loc_ntrec, &
               & int_phys_ocnatm_timeslice(ipoa_KCO2,:,:)/int_t_timeslice,loc_mask_sur)
       end if
       ! (10) MLD
       loc_unitsname = 'm'
       call sub_adddef_netcdf(loc_iou,3,'climate_MLD','mixed layer depth',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('climate_MLD',loc_iou,n_i,n_j,loc_ntrec, &
            & int_phys_ocnatm_timeslice(ipoa_mld,:,:)/int_t_timeslice,loc_mask_sur)
       loc_unitsname = 'n/a'
       call sub_adddef_netcdf(loc_iou,3,'climate_MLD_k','mixed layer level',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('climate_MLD_k',loc_iou,n_i,n_j,loc_ntrec, &
            & int_phys_ocnatm_timeslice(ipoa_mld_k,:,:)/int_t_timeslice,loc_mask_sur)
    end if
    ! ---------------------------------------------------------------- !
    ! atmosphere climate fields
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_atm
       ia = conv_iselected_ia(l)
       SELECT CASE (atm_type(ia))
       CASE (0)
          loc_ij(:,:) = int_sfcatm1_timeslice(ia,:,:)/int_t_timeslice
       END SELECT
       if (ia == ia_T) loc_shortname = 'climate_SAT'
       if (ia == ia_q) loc_shortname = 'climate_humidity'
       SELECT CASE (atm_type(ia))
       CASE (0)
          call sub_adddef_netcdf(loc_iou,3,trim(loc_shortname), &
               & trim(string_atm_tlname(l)),trim(string_atm_unit(l)),atm_mima(l,1),atm_mima(l,2))
          call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur_ALL)
       END SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! ocean surface climate fields
    ! ---------------------------------------------------------------- !
    DO l=1,n_l_ocn
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_zero
       DO i=1,n_i
          DO j=1,n_j
             loc_k1 = goldstein_k1(i,j)
             IF (n_k >= loc_k1) THEN
                SELECT CASE (ocn_type(io))
                CASE (0)
                   if (io == io_T) then
                      loc_ij(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice - const_zeroC
                   else
                      loc_ij(i,j) = int_ocn_timeslice(io,i,j,n_k)/int_t_timeslice
                   end if
                end SELECT
             end IF
          end DO
       end DO
       if (io == io_T) loc_shortname = 'climate_SST'
       if (io == io_S) loc_shortname = 'climate_SSS'
       if (io == io_T) loc_unitsname = 'degrees C'
       if (io == io_S) loc_unitsname = 'PSU'
       SELECT CASE (ocn_type(io))
       CASE (0)
          ! T,S
             call sub_adddef_netcdf(loc_iou, 3, trim(loc_shortname), &
                  & 'surface-water '//trim(string_ocn(io)), trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d(trim(loc_shortname),loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
       end SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! MOC
    ! ---------------------------------------------------------------- !
    ! NOTE: flip the data in the vertical to match the axes ...
    ! global
    loc_mask_opsi(:,:)  = const_real_one
    loc_tmp_jk(:,:) = loc_scale*int_opsi_timeslice(:,:)/int_t_timeslice
    loc_tmp_jk(:,n_k:0:-1) = loc_tmp_jk(:,0:n_k:1)
    where(abs(loc_tmp_jk) < const_real_nullsmall)
       loc_mask_opsi = const_real_zero
    endwhere
    call sub_adddef_netcdf_moc(loc_iou,'climate_opsi','Global streamfunction','Sv',const_real_zero,const_real_zero)
    call sub_putvar2d('climate_opsi',loc_iou,n_j+1,n_k+1,loc_ntrec,loc_tmp_jk,loc_mask_opsi)
    ! Atlantic & Pacific -- modern topos only
    select case (fname_topo)
    case ('worbe2', 'worjh2', 'worjh4', 'worlg2', 'worlg4', 'wv2jh2', 'wv3jh2', 'worri4', 'p_worbe2', 'p_worjh2')
       ! Atlantic
       loc_tmp_jk(:,:) = loc_scale*int_opsia_timeslice(:,:)/int_t_timeslice
       loc_tmp_jk(:,n_k:0:-1) = loc_tmp_jk(:,0:n_k:1)
       loc_mask_opsi = const_real_one
       where(abs(loc_tmp_jk) < const_real_nullsmall)
          loc_mask_opsi = const_real_zero
       endwhere
       call sub_adddef_netcdf_moc(loc_iou,'climate_opsia','Atlantic streamfunction','Sv',const_real_zero,const_real_zero)
       call sub_putvar2d('climate_opsia',loc_iou,n_j+1,n_k+1,loc_ntrec,loc_tmp_jk,loc_mask_opsi)
       ! Pacific
       loc_tmp_jk(:,:) = loc_scale*int_opsip_timeslice(:,:)/int_t_timeslice
       loc_tmp_jk(:,n_k:0:-1) = loc_tmp_jk(:,0:n_k:1)
       loc_mask_opsi = const_real_one
       where(abs(loc_tmp_jk) < const_real_nullsmall)
          loc_mask_opsi = const_real_zero
       endwhere
       call sub_adddef_netcdf_moc(loc_iou,'climate_opsip','Pacific streamfunction','Sv',const_real_zero,const_real_zero)
       call sub_putvar2d('climate_opsip',loc_iou,n_j+1,n_k+1,loc_ntrec,loc_tmp_jk,loc_mask_opsi)
    end select
    ! ---------------------------------------------------------------- !
    ! PSI
    ! ---------------------------------------------------------------- !
    loc_tmp_ij(:,:)   = int_psi_timeslice(1:n_i,0:n_j)/int_t_timeslice
    loc_mask_psi(:,:) = const_real_one
    call sub_adddef_netcdf_psi(loc_iou,'climate_psi','Barotropic streamfunction','Sv',const_real_zero,const_real_zero)
    call sub_putvar2d('climate_psi',loc_iou,n_i,n_j+1,loc_ntrec,loc_tmp_ij,loc_mask_sur)
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_hidden_climate
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save fossil fuel CO2 related data ***
  SUBROUTINE sub_2d_save_hidden_fossilfuelco2()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,ia,io
    integer::loc_iou,loc_ntrec
    CHARACTER(len=255)::loc_unitsname
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask,loc_mask_sur
    real::loc_tot,loc_frac,loc_standard
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:)       = 0.0
    loc_mask          = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_sur(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    ! ---------------------------------------------------------------- !
    ! WATER-COLUMN INTEGRATED CO2 INVENTORIES
    ! ---------------------------------------------------------------- !
    loc_unitsname = 'mol m-2'
    DO l=3,n_l_ocn
       io = conv_iselected_io(l)
       loc_ij(:,:) = const_real_null
       DO i=1,n_i
          DO j=1,n_j
             If (goldstein_k1(i,j) <= n_k) then
                loc_ij(i,j) = &
                     & sum(phys_ocn(ipo_M,i,j,:)*int_ocn_timeslice(io,i,j,:))* &
                     & phys_ocn(ipo_rA,i,j,n_k)/int_t_timeslice
             end If
          end DO
       end DO
       SELECT CASE (io)
       CASE (io_DIC)
          call sub_adddef_netcdf(loc_iou,3,'CO2_wcint_'//trim(string_ocn(io)), &
               & trim(string_ocn(io))//' water-column integrated tracer inventory', &
               & trim(loc_unitsname),const_real_zero,const_real_zero)
          call sub_putvar2d('CO2_wcint_'//trim(string_ocn(io)),loc_iou, &
               & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_sur)
       END SELECT
    end DO
    ! ---------------------------------------------------------------- !
    ! save air-sea delta pCO2 data
    ! ---------------------------------------------------------------- !
    IF (atm_select(ia_pCO2)) THEN
       loc_ij(:,:) = phys_ocn(ipo_mask_ocn,:,:,n_k) * &
            & (int_carb_timeslice(ic_fug_CO2,:,:,n_k) - int_sfcatm1_timeslice(ia_pCO2,:,:))/int_t_timeslice
       loc_unitsname = 'atm'
       call sub_adddef_netcdf(loc_iou,3,'CO2_Dairsea','air-sea pCO2 diff',trim(loc_unitsname),const_real_zero,const_real_zero)
       call sub_putvar2d('CO2_Dairsea',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask_sur)
    end if
    ! ---------------------------------------------------------------- !
    ! save derived CO2 flux data
    ! ---------------------------------------------------------------- !
    loc_ij(:,:) = int_diag_airsea_timeslice(ia_pCO2,:,:)/int_t_timeslice
    loc_unitsname = 'mol yr-1'
    call sub_adddef_netcdf(                                                    &
         & loc_iou,3,'CO2_ftotgasex','pCO2 net sea->air gas exchange flux per grid point', &
         & trim(loc_unitsname),const_real_zero,const_real_zero                 &
         & )
    call sub_putvar2d('CO2_ftotgasex',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask)
    loc_ij(:,:) = (int_diag_airsea_timeslice(ia_pCO2,:,:)/int_t_timeslice) * &
         & (5.0/phys_ocnatm(ipoa_dlon,:,:))*(4.0/phys_ocnatm(ipoa_dlat,:,:))
    loc_unitsname = 'mol (5x4)-1 yr-1'
    call sub_adddef_netcdf(                                                                &
         & loc_iou,3,'CO2_fgridgasex','pCO2 net sea->air gas exchange flux per 5 x 4 grid', &
         & trim(loc_unitsname),const_real_zero,const_real_zero                             &
         & )
    call sub_putvar2d('CO2_fgridgasex',loc_iou,n_i,n_j,loc_ntrec,loc_ij,loc_mask)
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_hidden_fossilfuelco2
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! *** save interface flux data ***
  SUBROUTINE sub_2d_save_hidden_interfacefluxes()
    ! ---------------------------------------------------------------- !
    ! define local variables
    ! ---------------------------------------------------------------- !
    INTEGER::i,j,l,ia,io,is
    integer::loc_iou,loc_ntrec
    real::loc_tot,loc_frac,loc_standard
    real,DIMENSION(n_i,n_j)::loc_ij,loc_mask,loc_mask_surf,loc_sed_mask
    CHARACTER(len=255)::loc_unitsname
    ! ---------------------------------------------------------------- !
    ! initialize local variables
    ! ---------------------------------------------------------------- !
    loc_iou   = ncout2d_iou
    loc_ntrec = ncout2d_ntrec
    loc_ij(:,:) = 0.0
    loc_mask = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_mask_surf(:,:) = phys_ocnatm(ipoa_mask_ocn,:,:)
    loc_sed_mask(:,:)  = phys_ocn(ipo_mask_ocn,:,:,n_k)
    ! ---------------------------------------------------------------- !
    ! save flux density data
    ! ---------------------------------------------------------------- !
    ! NOTE: use atmospheric grid point physics array to avoid the zero
    !       area values of dry grid points in the (ocean) physics array
    ! NOTE: a positive value of the array represents net ocean to atmosphere transfer
    DO l=3,n_l_atm
       ia = conv_iselected_ia(l)
       SELECT CASE (atm_type(ia))
       CASE (1)
          loc_ij(:,:) = (int_diag_airsea_timeslice(ia,:,:)/phys_ocnatm(ipoa_A,:,:))/int_t_timeslice
          loc_unitsname = 'mol m-2 yr-1'
          call sub_adddef_netcdf(                                                                                     &
               & loc_iou,3,'interf_focnatm_'//trim(string_atm(ia)),trim(string_atm(ia))//                             &
               & ': net sea->air gas exchange flux density',                                                          &
               & trim(loc_unitsname),const_real_zero,const_real_zero                                                  &
               & )
          call sub_putvar2d ('interf_focnatm_'//trim(string_atm(ia)),loc_iou,n_i,n_j, &
               & loc_ntrec,loc_ij,loc_mask)
       end SELECT
    END DO
    ! ---------------------------------------------------------------- !
    ! sediment-ocean exchange data
    ! ---------------------------------------------------------------- !
    If (flag_sedgem) then
       !-------------------------------------------------------------- ! save ocn->sed interface flux data
       If (ctrl_data_save_slice_focnsed .OR. ctrl_data_save_slice_diag_proxy) then
          DO l=1,n_l_sed
             is = conv_iselected_is(l)
             loc_ij(:,:) = const_real_zero
             DO i=1,n_i
                DO j=1,n_j
                   SELECT CASE (sed_type(is))
                   CASE (par_sed_type_bio,par_sed_type_abio, &
                        & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                        & par_sed_type_scavenged)
                      loc_ij(i,j) = int_focnsed_timeslice(is,i,j)
                      loc_unitsname = 'mol yr-1'
                   CASE (par_sed_type_age)
                      if (int_focnsed_timeslice(sed_dep(is),i,j) > 0.0) then
                         loc_ij(i,j) = int_focnsed_timeslice(is,i,j)/int_focnsed_timeslice(sed_dep(is),i,j)
                         loc_unitsname = 'years'
                      end if
                   case (n_itype_min:n_itype_max)
                      loc_tot  = int_focnsed_timeslice(sed_dep(is),i,j)
                      loc_frac = int_focnsed_timeslice(is,i,j)
                      loc_standard = const_standards(sed_type(is))
                      loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.,const_real_null)
                      loc_unitsname = 'o/oo'
                   end SELECT
                end DO
             end DO
             SELECT CASE (sed_type(is))
             CASE (par_sed_type_bio,par_sed_type_abio, &
                  & par_sed_type_POM,par_sed_type_CaCO3,par_sed_type_opal,par_sed_type_det, &
                  & par_sed_type_scavenged,par_sed_type_age,n_itype_min:n_itype_max)
                call sub_adddef_netcdf(loc_iou,3,'interf_focnsed_'//trim(string_sed(is)), &
                     & trim(string_sed(is))//' ocean->sediment flux',trim(loc_unitsname),const_real_zero,const_real_zero)
                call sub_putvar2d('interf_focnsed_'//trim(string_sed(is)),loc_iou,n_i,n_j, &
                     & loc_ntrec, loc_ij, loc_sed_mask)
             end SELECT
          END DO
       end if
       !-------------------------------------------------------------- ! save sed->ocn interface flux data
       ! NOTE: exclude dissolved organic matter tracers
       If (ctrl_data_save_slice_fsedocn) then
          DO l=3,n_l_ocn
             io = conv_iselected_io(l)
             is = maxval(maxloc(abs(conv_DOM_POM(:,io))))-1
             if (is == 0) then
                loc_ij(:,:) = const_real_zero
                DO i=1,n_i
                   DO j=1,n_j
                      SELECT CASE (ocn_type(io))
                      CASE (1)
                         loc_ij(i,j) = int_fsedocn_timeslice(io,i,j)
                         loc_unitsname = 'mol yr-1'
                      case (n_itype_min:n_itype_max)
                         loc_tot  = int_fsedocn_timeslice(ocn_dep(io),i,j)
                         loc_frac = int_fsedocn_timeslice(io,i,j)
                         loc_standard = const_standards(ocn_type(io))
                         loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.TRUE.,const_real_null)
                         loc_unitsname = 'o/oo'
                      end SELECT
                   end DO
                end DO
                SELECT CASE (ocn_type(io))
                CASE (1,n_itype_min:n_itype_max)
                   call sub_adddef_netcdf(loc_iou,3,'interf_fsedocn_'//trim(string_ocn(io)), &
                        & trim(string_ocn(io))//' sediment->ocean flux',trim(loc_unitsname),const_real_zero,const_real_zero)
                   call sub_putvar2d('interf_fsedocn_'//trim(string_ocn(io)),loc_iou,n_i,n_j, &
                        & loc_ntrec,loc_ij,loc_sed_mask)
                end SELECT
             end if
          END DO
       end if
    end if
    ! ---------------------------------------------------------------- !
    ! weathering diagnostics data
    ! ---------------------------------------------------------------- !
    If (flag_rokgem) then
       DO l=3,n_l_ocn
          io = conv_iselected_io(l)
          loc_ij(:,:) = const_real_zero
          DO i=1,n_i
             DO j=1,n_j
                SELECT CASE (ocn_type(io))
                CASE (1)
                   loc_ij(i,j) = int_diag_weather_timeslice(io,i,j)
                   loc_unitsname = 'mol yr-1'
                case (n_itype_min:n_itype_max)
                   loc_tot  = int_diag_weather_timeslice(ocn_dep(io),i,j)
                   loc_frac = int_diag_weather_timeslice(io,i,j)
                   loc_standard = const_standards(ocn_type(io))
                   loc_ij(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_real_null)
                   loc_unitsname = 'o/oo'
                end SELECT
             end DO
          end DO
          SELECT CASE (ocn_type(io))
          CASE (1)
             ! bulk tracers
             call sub_adddef_netcdf(loc_iou,3,'interf_fweather_'//trim(string_ocn(io)), &
                  & 'weathering flux - '//trim(string_ocn(io)),trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('interf_fweather_'//trim(string_ocn(io)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
          CASE (n_itype_min:n_itype_max)
             ! isotopes
             call sub_adddef_netcdf(loc_iou,3,'interf_fweather_'//trim(string_ocn(io)), &
                  & 'weathering flux - '//trim(string_ocn(io)),trim(loc_unitsname),const_real_zero,const_real_zero)
             call sub_putvar2d('interf_fweather_'//trim(string_ocn(io)),loc_iou, &
                  & n_i,n_j,loc_ntrec,loc_ij(:,:),loc_mask_surf)
          end SELECT
       END DO
    end If
    ! ---------------------------------------------------------------- !
  END SUBROUTINE sub_2d_save_hidden_interfacefluxes
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! ****************************************************************************************************************************** !


END MODULE biogem_data_netCDF
