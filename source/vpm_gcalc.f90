Subroutine calc_velocity_serial_3d(idcalc)
    use vpm_vars         
    use pmeshpar
    use parvar
    use pmgrid
    use openmpth
    Implicit None
    integer, intent(in) :: idcalc
    double precision ::  dphidx, dphidy,dphidz,dpsidx(3), dpsidy(3),dpsidz(3)
    double precision ::  wdudx, wdvdy, wdwdz, velxp, velyp, velzp, velxm, velym, velzm
    double precision ::  upi,umi,upj,umj,upk,umk
    double precision ::  vpi,vmi,vpj,vmj,vpk,vmk
    double precision ::  wpi,wmi,wpj,wmj,wpk,wmk
    integer          :: i, j, k

  if (idcalc.ge.0) then 
    DXpm2 = 2 * DXpm
    DYpm2 = 2 * DYpm
    DZpm2 = 2 * DZpm
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
                velvrx_pm(i, j, k)  = velvrx_pm(i, j, k) +  dpsidy(3) - dpsidz(2)  
                velvry_pm(i, j, k)  = velvry_pm(i, j, k) - (dpsidx(3) - dpsidz(1)) 
                velvrz_pm(i, j, k)  = velvrz_pm(i, j, k) +  dpsidx(2) - dpsidy(1)   

               !velvrx_pm(i, j, k)  =  +  dpsidy(3) - dpsidz(2)
               !velvry_pm(i, j, k)  =  - (dpsidx(3) - dpsidz(1))
               !velvrz_pm(i, j, k)  =  +  dpsidx(2) - dpsidy(1)
            enddo
        enddo
    enddo
   !$omp enddo
   !$omp endparallel

 if (neqpm.eq.4) then 
    !$omp parallel private(i,j,k,dphidx,dphidy,dphidz,dpsidx,dpsidy,dpsidz) num_threads(OMPTHREADS)
    !$omp do
    do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
        do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
           do i = NXs_bl(1) + 1, NXf_bl(1) - 1

                !--> dpsi(x,y,z)d(xyz)
                dpsidx(1) = (SOL_pm(4,i + 1, j, k)  - SOL_pm(4,i - 1, j, k)) / DXpm2
                dpsidy(1) = (SOL_pm(4,i, j + 1, k)  - SOL_pm(4,i, j - 1, k)) / DYpm2
                dpsidz(1) = (SOL_pm(4,i, j, k + 1)  - SOL_pm(4,i, j, k - 1)) / DZpm2
                ! U = grad x psi
                velvrx_pm(i, j, k)  = velvrx_pm(i, j, k) + dpsidx(1)
                velvry_pm(i, j, k)  = velvry_pm(i, j, k) + dpsidy(1)
                velvrz_pm(i, j, k)  = velvrz_pm(i, j, k) + dpsidz(1)
            enddo
        enddo
    enddo
    !$omp enddo
    !$omp endparallel
 endif
endif!
  !Sol of vorticity is no longer need thus we use it for storing deformation
   if (idcalc.eq.0) return
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
 !DEFORMATION WITH A MINUS BECAUSE WE HAVE STORED MINUS VORTICITY
                wdudx  = -(RHS_pm(1,i + 1,j,k) * upi - (RHS_pm(1,i - 1,j,k)) * umi)/DXpm2
                wdvdy  = -(RHS_pm(2,i,j + 1,k) * upj - (RHS_pm(2,i,j - 1,k)) * umj)/DYpm2
                wdwdz  = -(RHS_pm(3,i,j,k + 1) * upk - (RHS_pm(3,i,j,k - 1)) * umk)/DZpm2

                SOL_pm(1,i,j,k) =  wdudx + wdvdy + wdwdz

                ! Wy * (thu/thx + thv/thy + thw/thz)
                wdudx  = -(RHS_pm(1,i + 1,j,k) * vpi - (RHS_pm(1,i - 1,j,k)) * vmi)/DXpm2
                wdvdy  = -(RHS_pm(2,i,j + 1,k) * vpj - (RHS_pm(2,i,j - 1,k)) * vmj)/DYpm2
                wdwdz  = -(RHS_pm(3,i,j,k + 1) * vpk - (RHS_pm(3,i,j,k - 1)) * vmk)/DZpm2

                SOL_pm(2,i,j,k) =  wdudx + wdvdy + wdwdz

                ! Wy * (thu/thx + thv/thy + thw/thz)
                wdudx  = -(RHS_pm(1,i + 1,j,k) * wpi  - (RHS_pm(1,i - 1,j,k)) * wmi)/DXpm2
                wdvdy  = -(RHS_pm(2,i,j + 1,k) * wpj  - (RHS_pm(2,i,j - 1,k)) * wmj)/DYpm2
                wdwdz  = -(RHS_pm(3,i,j,k + 1) * wpk  - (RHS_pm(3,i,j,k - 1)) * wmk)/DZpm2

                SOL_pm(3,i,j,k) =  wdudx + wdvdy + wdwdz
            enddo
        enddo
    enddo
   !$omp enddo
   !$omp endparallel

End Subroutine calc_velocity_serial_3d

Subroutine diffuse_vort_3d
    use vpm_vars         
    use pmeshpar
    use parvar
    use pmgrid
    use openmpth
    Implicit None
    double precision ::  dwxdx, dwydy,dwzdz,VIS
    integer          :: i, j, k


    DXpm2 =  DXpm**2
    DYpm2 =  DYpm**2 
    DZpm2 =  DZpm**2
    Sol_pm=0.d0
    !$omp parallel private(i,j,k,dwxdx,dwydy,dwzdz,VIS) num_threads(OMPTHREADS)
    !$omp do
    do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
        do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
           do i = NXs_bl(1) + 1, NXf_bl(1) - 1
                if (neqpm.eq.3) then 
                  VIS =  NI
                else 
                  VIS = RHS_pm(4,i,j,k) + NI
                endif
                !write(*,*) VIS,neqpm
                !--> Remember that RHS = -w 
                dwxdx = (RHS_pm(1, i + 1, j, k)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i - 1, j, k)) / DXpm2
                dwydy = (RHS_pm(1, i, j + 1, k)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i, j - 1, k)) / DYpm2
                dwzdz = (RHS_pm(1, i, j, k + 1)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                SOL_pm(1,i,j,k) = -VIS*(dwxdx+dwydy+dwzdz) ! because RHS=-w

                dwxdx = (RHS_pm(2, i + 1, j, k)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i - 1, j, k)) / DXpm2
                dwydy = (RHS_pm(2, i, j + 1, k)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i, j - 1, k)) / DYpm2
                dwzdz = (RHS_pm(2, i, j, k + 1)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                SOL_pm(2,i,j,k) = -VIS*(dwxdx+dwydy+dwzdz) ! because RHS=-w

                dwxdx = (RHS_pm(3, i + 1, j, k)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i - 1, j, k)) / DXpm2
                dwydy = (RHS_pm(3, i, j + 1, k)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i, j - 1, k)) / DYpm2
                dwzdz = (RHS_pm(3, i, j, k + 1)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                SOL_pm(3,i,j,k) = -VIS*(dwxdx+dwydy+dwzdz) ! because RHS=-w
            enddo
        enddo
    enddo
    !$omp enddo
    !$omp endparallel

   !!$omp parallel private(i,j,k) num_threads(OMPTHREADS)
   !!$omp do
   !do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
   !    do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
   !       do i = NXs_bl(1) + 1, NXf_bl(1) - 1

   !            !--> Remember that RHS = -w 
   !            RHS_pm(1:3, i, j, k) = RHS_pm(1:3, i, j, k) - NI * SOL_pm(1,i,j,k) 
   !        enddo
   !    enddo
   !enddo
   !!$omp enddo
   !!$omp endparallel
End Subroutine diffuse_vort_3d 


Subroutine calc_antidiffusion
use vpm_vars         
use pmeshpar
use parvar
use pmgrid
use openmpth
Implicit None
double precision ::  dwxdx, dwydy,dwzdz,Ct
integer          :: i, j, k
double precision,allocatable :: laplvort(:,:,:,:)

    allocate (laplvort(3,NXpm,NYpm,NZpm))
    laplvort=0.d0
    DXpm2 =  DXpm**2
    DYpm2 =  DYpm**2 
    DZpm2 =  DZpm**2
    Sol_pm=0.d0
    Ct = 6.8 * DXpm**2 /4
    do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
        do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
           do i = NXs_bl(1) + 1, NXf_bl(1) - 1
                !write(*,*) VIS,neqpm
                !--> Remember that RHS = -w 
                dwxdx =-(RHS_pm(1, i + 1, j, k)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i - 1, j, k)) / DXpm2
                dwydy =-(RHS_pm(1, i, j + 1, k)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i, j - 1, k)) / DYpm2
                dwzdz =-(RHS_pm(1, i, j, k + 1)  - 2 * RHS_pm(1, i, j, k) &
                       + RHS_pm(1, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                laplvort(1,i,j,k)=dwxdx+dwydy+dwzdz

                dwxdx =-(RHS_pm(2, i + 1, j, k)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i - 1, j, k)) / DXpm2
                dwydy =-(RHS_pm(2, i, j + 1, k)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i, j - 1, k)) / DYpm2
                dwzdz =-(RHS_pm(2, i, j, k + 1)  - 2 * RHS_pm(2, i, j, k) &
                       + RHS_pm(2, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                laplvort(2,i,j,k)=dwxdx+dwydy+dwzdz

                dwxdx =-(RHS_pm(3, i + 1, j, k)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i - 1, j, k)) / DXpm2
                dwydy =-(RHS_pm(3, i, j + 1, k)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i, j - 1, k)) / DYpm2
                dwzdz =-(RHS_pm(3, i, j, k + 1)  - 2 * RHS_pm(3, i, j, k) &
                       + RHS_pm(3, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                laplvort(3,i,j,k)=dwxdx+dwydy+dwzdz
            enddo
        enddo
    enddo

    do k = NZs_bl(1) + 1, NZf_bl(1)- 1 
        do j = NYs_bl(1) + 1, NYf_bl(1 )- 1
           do i = NXs_bl(1) + 1, NXf_bl(1) - 1
                !write(*,*) VIS,neqpm
                !Minus because of (-w) has been included in laplvort
                dwxdx = (laplvort(1, i + 1, j, k)  - 2 * laplvort(1, i, j, k) &
                       + laplvort(1, i - 1, j, k)) / DXpm2
                dwydy = (laplvort(1, i, j + 1, k)  - 2 * laplvort(1, i, j, k) &
                       + laplvort(1, i, j - 1, k)) / DYpm2
                dwzdz = (laplvort(1, i, j, k + 1)  - 2 * laplvort(1, i, j, k) &
                       + laplvort(1, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                Sol_pm(1,i,j,k)=Sol_pm(1,i,j,k) + Ct*dwxdx+dwydy+dwzdz

                dwxdx = (laplvort(2, i + 1, j, k)  - 2 * laplvort(2, i, j, k) &
                       + laplvort(2, i - 1, j, k)) / DXpm2
                dwydy = (laplvort(2, i, j + 1, k)  - 2 * laplvort(2, i, j, k) &
                       + laplvort(2, i, j - 1, k)) / DYpm2
                dwzdz = (laplvort(2, i, j, k + 1)  - 2 * laplvort(2, i, j, k) &
                       + laplvort(2, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                Sol_pm(2,i,j,k)=Sol_pm(2,i,j,k) + Ct*dwxdx+dwydy+dwzdz

                dwxdx = (laplvort(3, i + 1, j, k)  - 2 * laplvort(3, i, j, k) &
                       + laplvort(3, i - 1, j, k)) / DXpm2
                dwydy = (laplvort(3, i, j + 1, k)  - 2 * laplvort(3, i, j, k) &
                       + laplvort(3, i, j - 1, k)) / DYpm2
                dwzdz = (laplvort(3, i, j, k + 1)  - 2 * laplvort(3, i, j, k) &
                       + RHS_pm(3, i, j, k - 1)) / DZpm2
                ! U = grad x psi
                Sol_pm(3,i,j,k)=Sol_pm(3,i,j,k) + Ct*dwxdx+dwydy+dwzdz
            enddo
        enddo
    enddo

End Subroutine calc_antidiffusion
