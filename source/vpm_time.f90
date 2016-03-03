Subroutine convect_first_order(Xbound,Dpm,NN,NN_Bl)
    use vpm_vars        
    use pmeshpar
    use parvar
    use pmgrid
    use MPI
    use pmlib
    use projlib
    use yapslib
    double precision,intent(in)  :: Xbound(6),Dpm(3)
    integer,intent(in)           :: NN(3),NN_bl(6)
    double precision,allocatable :: QINF(:)
    integer         ,allocatable :: ieq(:)
    integer                      :: i

   
     do i= 1,NVR
         XP(1:3,i) = XP(1:3,i)  + UP(1:3,i) * DT_c
         QP(1:3,i) = QP(1:3,i)  + GP(1:3,i) * DT_c
     enddo
End Subroutine convect_first_order

!---------------------------------------------------------------------------!
!-> Subroutine back_to_particles                                            !
!   This subroutine interpolates PM grid values back to particles at the    !
!   positions they ARE.Convections takes place afterwards.                  !
!   Input :                                                                 !
!          itype (1,2) defines what value to interpolate to the particles   !
!---------------------------------------------------------------------------!
Subroutine back_to_particles_3D(SOL_pm,RHS_pm,XP,QP,UP,GP,&
                                velvrx_pm,velvry_pm,velvrz_pm,&
                                Xbound,Dpm,NN,NN_bl,NVR,neqpm,iproj,itype,NVRM)

    use openmpth
    Implicit None
    integer, intent(in)             :: NN(3),NVR,iproj,NN_bl(6),neqpm,NVRM,itype
    double precision, intent(in)    :: RHS_pm(neqpm+1,NN(1),NN(2),NN(3))
    double precision, intent(in)    :: SOL_pm(neqpm,NN(1),NN(2),NN(3))
    double precision, intent(inout) :: QP(neqpm+1,NVRM),XP(3,NVRM),UP(3,NVRM),GP(3,NVRM)
!   double precision, intent(inout) :: QP(:,:),XP(:,:),UP(:,:),GP(:,:)
    double precision, intent(in)    :: velvrx_pm(NN(1),NN(2),NN(3)),velvry_pm(NN(1),NN(2),NN(3))
    double precision, intent(in)    :: velvrz_pm(NN(1),NN(2),NN(3))

    double precision, intent(in)    :: Xbound(6),Dpm(3)

    double precision :: fx, fy, fz, f, x, y, z, projection_fun,Gloc(3)
    integer          :: inode, jnode, knode, i, j , k, nv,ips,ipf,ivortx,ivorty,DVpm



    if(iproj.eq.2) then
        ips =  0
        ipf =  1
    else if (iproj.eq.3) then
        ips =  1
        ipf =  2
    else if (iproj.eq.4) then
        ips =  1
        ipf =  2
    endif
    DVpm = Dpm(1)*Dpm(2)*Dpm(3)
    if(itype.eq.1) then 
    !Itype==1 normal back to part
       !!$omp parallel private(nv,inode,jnode,knode,ivortx,ivorty,x,y,z,fx,fy,fz,f,i,j,k) num_threads(OMPTHREADS)
       !!$omp do
        do  nv = 1, NVR
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1,nv) - XBound(1)) / Dpm(1)) + 1
            jnode = int((XP(2,nv) - XBound(2)) / Dpm(2)) + 1
            knode = int((XP(3,nv) - XBound(3)) / Dpm(3)) + 1
            !--We search the 4 nodes close to the particles
            Gloc(1:3)=0.d0
            do k = knode - ips, knode + ipf
               do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x  = (XP(1,nv)- XBound(1) - (i-1) * Dpm(1)) / Dpm(1)
                        fx = projection_fun(iproj,x)

                        y  = (XP(2,nv)- XBound(2) - (j-1) * Dpm(2)) / Dpm(2)
                        fy = projection_fun(iproj,y)

                        z  = (XP(3,nv)- XBound(3) - (k-1) * Dpm(3)) / Dpm(3)
                        fz = projection_fun(iproj,z)



                        f = fx * fy * fz
                      
                       !QP(1,nv)    = QP(1,nv)    +  f * RHS_pm(1,i,j,k)
                       !QP(2,nv)    = QP(2,nv)    +  f * RHS_pm(2,i,j,k)
                       !QP(3,nv)    = QP(3,nv)    +  f * RHS_pm(3,i,j,k)
                        UP(1,nv)    = UP(1,nv)    +  f * (velvrx_pm(i,j,k))
                        UP(2,nv)    = UP(2,nv)    +  f * (velvry_pm(i,j,k))
                        UP(3,nv)    = UP(3,nv)    +  f * (velvrz_pm(i,j,k))

                        Gloc(1:3)    = Gloc(1:3)   +  f * SOL_pm(1:3,i,j,k) 

                    enddo
                enddo
            enddo

         !  QP(1:3,nv)  = QP(1:3,nv) * DVpm!QP(neqpm+1,nv)
            GP(1:3,nv)  =  Gloc(1:3)*QP(neqpm+1,nv)
        enddo
       !!$omp enddo
       !!$omp endparallel
      endif

      if (itype.eq.2) then 
       !diffusion
        GP(1:3,:)=0.d0
      !Itype==1 normal back to part
        !$omp parallel private(nv,inode,jnode,knode,ivortx,ivorty,x,y,z,fx,fy,fz,f,i,j,k) num_threads(OMPTHREADS)
        !$omp do
        do  nv = 1, NVR
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1,nv) - XBound(1)) / Dpm(1)) + 1
            jnode = int((XP(2,nv) - XBound(2)) / Dpm(2)) + 1
            knode = int((XP(3,nv) - XBound(3)) / Dpm(3)) + 1
            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
               do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x  = (XP(1,nv)- XBound(1) - (i-1) * Dpm(1)) / Dpm(1)
                        fx = projection_fun(iproj,x)

                        y  = (XP(2,nv)- XBound(2) - (j-1) * Dpm(2)) / Dpm(2)
                        fy = projection_fun(iproj,y)

                        z  = (XP(3,nv)- XBound(3) - (k-1) * Dpm(3)) / Dpm(3)
                        fz = projection_fun(iproj,z)



                        f = fx * fy * fz
                      
                        GP(1,nv)    = GP(1,nv)    +  f * Sol_pm(1,i,j,k)
                        GP(2,nv)    = GP(2,nv)    +  f * Sol_pm(2,i,j,k)
                        GP(3,nv)    = GP(3,nv)    +  f * Sol_pm(3,i,j,k)

                    enddo
                enddo
            enddo

            GP(1,nv)  = GP(1,nv) * QP(neqpm+1,nv)
            GP(2,nv)  = GP(2,nv) * QP(neqpm+1,nv)
            GP(3,nv)  = GP(3,nv) * QP(neqpm+1,nv)
        enddo
        !$omp enddo
        !$omp endparallel
      endif
End Subroutine back_to_particles_3D


