Module test_mod
  double precision,allocatable,target:: XPR(:,:),QPR(:,:),UPR(:,:),GPR(:,:)
  double precision,pointer:: velx(:,:,:), vely(:,:,:), velz(:,:,:)
  double precision,pointer:: RHS_pm_in(:,:,:,:)
  integer :: NVR_ext
end Module test_mod

