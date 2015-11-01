!not needded outside PM
module pmeshpar
    double precision , save              :: PI, PI2, PI4,DT

    double precision, allocatable,save   :: velx_pm(:,:,:), vely_pm(:,:,:),velz_pm(:,:,:),qx_pm(:,:,:), qy_pm(:,:,:)
    double precision, allocatable,save   :: velphix_pm(:,:,:),velphiy_pm(:,:,:),velphiz_pm(:,:,:)
    double precision,allocatable, save   :: SOL_pm(:,:,:,:),SOL_0_pm(:,:,:,:)
    double precision,allocatable,save    :: source_bound(:,:),x_s(:,:),y_s(:,:),z_s(:,:),d_s(:),cos_s(:),sin_s(:)
    double precision,allocatable,save    :: source_bound_lev(:,:,:),xs_lev(:,:,:),ys_lev(:,:,:),zs_lev(:,:,:),ds_lev(:,:,:)
    integer         ,allocatable,save    :: nbound_lev(:)
    double precision,allocatable,save    :: Psiz_pm_0(:,:,:),Psiz_pm_f(:,:,:)
    double precision, allocatable,save   :: Cont_pm(:,:)
    integer                              :: levmax,npar_cell,ND

    integer, save                        :: nbound,ndumcell,NVR_CFD_sa,IDVPM

end module pmeshpar

module pmgrid
    double precision,allocatable,target,save       :: velvrx_pm(:,:,:),velvry_pm(:,:,:),velvrz_pm(:,:,:)
    double precision,allocatable, target,save       :: RHS_pm(:,:,:,:)

    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DD, DXpm, DYpm, DZpm,DXpm2,DYpm2,DZpm2,DVpm,EPSVOL
    integer ,save                        :: NXpm, NYpm, NZpm, NXpm_par,NYpm_par,NZpm_par    
    integer ,save                        :: NXs_bl(10),NYs_bl(10),NXf_bl(10),NYf_bl(10),NZs_bl(10),NZf_bl(10),NBlocks

end module pmgrid
!     QP
!   -->1 Vorticity X
!   -->2 Vorticity Y
!   -->3 Vorticity Z
!   -->4 Dilatation
!   -->5 Pseudopressure
!   -->6 Mass
!   -->7 Volume
module parvar
    integer                          ::NVR
    double precision, pointer , save :: XP(:,:),QP(:,:)
    double precision, pointer , save :: UP(:,:),GP(:,:)
    integer,allocatable,save            :: NVR_projtype(:)

    double precision,allocatable :: XP_CFD_Sa(:,:),UT_CFD_Sa(:),UN_CFD_Sa(:),&
        DS_CFD_Sa(:)  ,Vo_CFD_Sa(:)
end module parvar

module output
    character * 15          :: filename1,filename2,filename3,filename4
end module output
