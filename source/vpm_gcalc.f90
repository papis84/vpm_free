Subroutine calc_velocity_serial_3d
    use vpm_vars         
    use pmeshpar
    use parvar
    use pmgrid
    use openmpth
    Implicit None
    double precision ::  dphidx, dphidy,dphidz,dpsidx(3), dpsidy(3),dpsidz(3)
    double precision ::  wdudx, wdvdy, wdwdz, velxp, velyp, velzp, velxm, velym, velzm
    double precision ::  upi,umi,upj,umj,upk,umk
    double precision ::  vpi,vmi,vpj,vmj,vpk,vmk
    double precision ::  wpi,wmi,wpj,wmj,wpk,wmk
    integer          :: i, j, k


    DXpm2 = 2 * DXpm
    DYpm2 = 2 * DYpm
    DZpm2 = 2 * DZpm

    velvrx_pm=0
    velvry_pm=0
    velvrz_pm=0
    !$omp parallel private(i,j,k,dphidx,dphidy,dphidz,dpsidx,dpsidy,dpsidz) num_threads(OMPTHREADS)
    !$omp do
    do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
        do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
           do i = NXs_bl(1) + 1, NXf_bl(1) - 1

                !--> dpsi(x,y,z)d(xyz)
                dpsidx(1:3) = (SOL_pm(1:3,i + 1, j, k)  - SOL_pm(1:3,i - 1, j, k)) / DXpm2
                dpsidy(1:3) = (SOL_pm(1:3,i, j + 1, k)  - SOL_pm(1:3,i, j - 1, k)) / DYpm2
                dpsidz(1:3) = (SOL_pm(1:3,i, j, k + 1)  - SOL_pm(1:3,i, j, k - 1)) / DZpm2
                ! U = grad x psi
                velvrx_pm(i, j, k)  = +  dpsidy(3) - dpsidz(2)
                velvry_pm(i, j, k)  = - (dpsidx(3) - dpsidz(1))
                velvrz_pm(i, j, k)  = +  dpsidx(2) - dpsidy(1)
            enddo
        enddo
    enddo
    !$omp enddo
    !$omp endparallel
  !Sol of vorticity is no longer need thus we use it for storing deformation
    SOL_pm=0.d0
  !REMEMBER VORTICITY CARRIED IS -OMEGA and the quantity transferes is -OMEGA thus
  !deformation = - omega*\gradu 
  !$omp parallel private(i,j,k,velxp,velxm,velyp,velym,velzp,velzm,wdudx,wdvdy,wdwdz,&
  !$omp                  upi,upj,upk,vpi,vpj,vpk,wpi,wpj,wpk,umi,umj,umk,vmi,vmj,vmk,&
  !$omp                  wmi,wmj,wmk) num_threads(OMPTHREADS)
  !$omp do
    do k = NZs_bl(1) + 2, NZf_bl(1) - 2
        do j = NYs_bl(1) + 2, NYf_bl(1) - 2 
           do i =  NXs_bl(1) + 2, NXf_bl(1) -2 
                velxp = velvrx_pm(i + 1, j, k) 
                velxm = velvrx_pm(i - 1, j, k) 

                velyp = velvry_pm(i,j + 1, k)
                velym = velvry_pm(i,j - 1, k)

                velzp = velvrz_pm(i,j, k + 1) 
                velzm = velvrz_pm(i,j, k - 1)

                upi = velvrx_pm(i+1,j,k) 
                upj = velvrx_pm(i,j+1,k) 
                upk = velvrx_pm(i,j,k+1) 

                vpi = velvry_pm(i+1,j,k) 
                vpj = velvry_pm(i,j+1,k) 
                vpk = velvry_pm(i,j,k+1) 

                wpi = velvrz_pm(i+1,j,k) 
                wpj = velvrz_pm(i,j+1,k) 
                wpk = velvrz_pm(i,j,k+1) 

                umi = velvrx_pm(i-1,j,k) 
                umj = velvrx_pm(i,j-1,k) 
                umk = velvrx_pm(i,j,k-1) 

                vmi = velvry_pm(i-1,j,k) 
                vmj = velvry_pm(i,j-1,k) 
                vmk = velvry_pm(i,j,k-1) 

                wmi = velvrz_pm(i-1,j,k) 
                wmj = velvrz_pm(i,j-1,k) 
                wmk = velvrz_pm(i,j,k-1) 

                wdudx  = (RHS_pm(1,i + 1,j,k) * upi - (RHS_pm(1,i - 1,j,k)) * umi)/DXpm2
                wdvdy  = (RHS_pm(2,i,j + 1,k) * upj - (RHS_pm(2,i,j - 1,k)) * umj)/DYpm2
                wdwdz  = (RHS_pm(3,i,j,k + 1) * upk - (RHS_pm(3,i,j,k - 1)) * umk)/DZpm2

                SOL_pm(1,i,j,k) =  wdudx + wdvdy + wdwdz

                ! Wy * (thu/thx + thv/thy + thw/thz)
                wdudx  = (RHS_pm(1,i + 1,j,k) * vpi - (RHS_pm(1,i - 1,j,k)) * vmi)/DXpm2
                wdvdy  = (RHS_pm(2,i,j + 1,k) * vpj - (RHS_pm(2,i,j - 1,k)) * vmj)/DYpm2
                wdwdz  = (RHS_pm(3,i,j,k + 1) * vpk - (RHS_pm(3,i,j,k - 1)) * vmk)/DZpm2

                SOL_pm(2,i,j,k) =  wdudx + wdvdy + wdwdz

                ! Wy * (thu/thx + thv/thy + thw/thz)
                wdudx  = (RHS_pm(1,i + 1,j,k) * wpi  - (RHS_pm(1,i - 1,j,k)) * wmi)/DXpm2
                wdvdy  = (RHS_pm(2,i,j + 1,k) * wpj  - (RHS_pm(2,i,j - 1,k)) * wmj)/DYpm2
                wdwdz  = (RHS_pm(3,i,j,k + 1) * wpk  - (RHS_pm(3,i,j,k - 1)) * wmk)/DZpm2

                SOL_pm(3,i,j,k) =  wdudx + wdvdy + wdwdz
            enddo
        enddo
    enddo
   !$omp enddo
   !$omp endparallel

End Subroutine calc_velocity_serial_3d
