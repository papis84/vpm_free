Module test_mod
  double precision,allocatable,target:: XPR(:,:),QPR(:,:),UPR(:,:),GPR(:,:),XPO(:,:),QPO(:,:),&
                                        XP_in(:,:),QP_in(:,:),XSOUR(:,:),QSOUR(:,:),XP_all(:,:),QP_all(:,:)
  double precision,pointer:: velx(:,:,:), vely(:,:,:), velz(:,:,:)
  double precision,pointer:: RHS_pm_in(:,:,:,:),RHS_pm_out(:,:,:,:)
  integer ,allocatable    :: qflag(:)
  integer :: NVR_ext,NVR_sources,NVR_sources_init,NVR_ext_init
end Module test_mod

