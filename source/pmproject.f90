Module projlib

    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DXpm, DYpm, DZpm,DVpm
    double precision, save               :: EPSVOL
    integer ,save                        :: NXpm, NYpm, NZpm,NXs,NXf,NYs,NYf,NZs,NZf
    integer ,save                        :: IDVPM,ND



    private :: XMIN_pm,XMAX_pm,YMIN_pm,YMAX_pm,ZMIN_pm,ZMAX_pm,DXpm,DYpm,DZpm,NXpm,NYpm,NZpm,ND,EPSVOL,DVpm,IDVPM
    private :: NXs,NXf,NYs,NYf,NZs,NZf
contains
    Subroutine projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
        double precision, intent(in) :: Xbound(6),Dpm(3),EPSVOL
        integer         , intent(in) :: NN(3),IDVPM,ND,NN_bl(6)


        XMIN_pm=Xbound(1);YMIN_pm=Xbound(2);ZMIN_pm=Xbound(3)
        XMAX_pm=Xbound(4);YMAX_pm=Xbound(5);ZMAX_pm=Xbound(6)
        DXpm   =Dpm(1)   ;DYpm   = Dpm(2)  ;DZpm   =Dpm(3)

        NXpm = NN(1);NYpm=NN(2);NZpm=NN(3)
        NXs  = NN_bl(1);NYs=NN_bl(2);NZs=NN_bl(3)
        NXf  = NN_bl(4);NYf=NN_bl(5);NZf=NN_bl(6)
        if (ND.eq.2) then
            DVpm = DXpm*DYpm
        else
            DVpm = DXpm*DYpm*DZPm
        endif

    End Subroutine projlibinit
    ! --------------------------------------------------------------------------!
    !-->Subroutine project_particles_3D                                          !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    Subroutine project_particles_3D(Qproj,Qpar,QpX,Qprojtype,ipar,isize,ieq,neq,QINF,iparsize)
        !use pmgrid

        Implicit None
        integer,intent(in) :: ipar,isize,ieq(neq),iparsize
        double precision ,intent(out):: Qproj(isize,NXpm,NYpm,NZpm)
        double precision ,intent(in) :: Qpar(isize,iparsize),QpX(3,iparsize),QINF(neq)
        integer,intent(in) :: Qprojtype(iparsize)
        double precision   :: fx, fy, fz, f, x, y, z, projection_fun
        integer            :: inode, jnode, knode, i, j , k, nv, itype,nbj,nb,neq,ips,ipf

        !-->Projection function (TSC)
        Qproj = 0.d0
        do  nv = 1, ipar
            itype = Qprojtype(nv)
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((QpX(1,nv) - XMIN_pm) / DXpm) + 1
            jnode = int((QpX(2,nv) - YMIN_pm) / DYpm) + 1
            knode = int((QpX(3,nv) - ZMIN_pm) / DZpm) + 1
            if(itype.eq.2) then
                ips =  0
                ipf =  1
            else 
                ips =  1
                ipf =  2
            endif

            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
                do j = jnode - ips, jnode  + ipf
                    do i = inode - ips, inode  + ipf

                        x  = (QpX(1,nv)- XMIN_pm - (i-1) * DXpm) / DXpm
                        fx = projection_fun(itype,x)

                        y  = (QpX(2,nv)- YMIN_pm - (j-1) * DYpm) / DYpm
                        fy = projection_fun(itype,y)

                        z  = (QpX(3,nv)- ZMIN_pm - (k-1) * DZpm) / DZpm
                        fz = projection_fun(itype,z)

                        f = fx * fy * fz
                        if (k.lt.1) print *, QPX(1:3,nv),nv,ipar,ZMIN_pm,DZpm,nv
                        Qproj(ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k) +&
                            f * (QPar(ieq(1:neq-1),nv))!-QINF(1:neq-1)*QPar(ieq(neq),nv))
                      
                        Qproj(ieq(neq),i,j,k) = Qproj(neq,i,j,k) +&
                            f * (QPar(ieq(neq),nv)-QINF(neq))
                    enddo
                enddo
            enddo
        enddo
        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions

    End Subroutine project_particles_3D

    !--------------------------------------------------------------------------!
    !-->Subroutine project_vol2d                                               !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    Subroutine project_vol3d(Qproj,isize,ieq,neq,iflag)
        !use pmgrid

        Implicit None
        integer,intent(in) :: isize,ieq(neq),iflag
        double precision   :: Qproj(isize,NXpm,NYpm,NZpm)
        double precision   :: fx, fy, fz, f, x, y, z, projection_fun,EPSVOLt
        integer            :: inode, jnode, knode, i, j , k, nv, itype,nbj,nb,neq,IDVPMt

        EPSVOLt=EPSVOL
        IDVPMt=IDVPM
        if(iflag.eq.1) IDVPMt=1
        !  if(iflag.eq.1) EPSVOLt=0.5
        if (IDVPMt.eq.0) then
            do k =1, NZpm
                do j =1, NYpm
                    do i =1, NXpm
                        if (i.lt.NXs.or.i.gt.NXf.or.j.lt.NYs.or.j.ge.NYf.or.k.lt.NZs.or.k.gt.NZf) then
                            Qproj (ieq(1:neq-1),i,j,k) = 0.d0
                            Qproj (ieq(neq),i,j,k) = DVpm
                            cycle
                        endif
                            Qproj (ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k)/(Qproj(ieq(neq),i,j,k))
                    enddo
                enddo
            enddo
        else if (IDVPMt.eq.1) then
            do k =1, NZpm
                do j =1, NYpm
                    do i =1, NXpm
                      ! if (i.lt.NXs.or.i.gt.NXf.or.j.lt.NYs.or.j.ge.NYf.or.k.lt.NZs.or.k.gt.NZf) then
                      !     Qproj (ieq(1:neq-1),i,j,k) = 0.d0
                      !     cycle
                      ! endif
                            Qproj (ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k)/ DVpm
                            Qproj (ieq(neq),i,j,k)=DVpm
                    enddo
                enddo
            enddo
        else
            write(*,*) 'WRONG IDVPM'
            STOP
        endif
    End Subroutine project_vol3d
    !--------------------------------------------------------------------------!
    !-->Subroutine project_particles                                           !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    Subroutine project_particles_2D(Qproj,Qpar,QpX,Qprojtype,ipar,isize,ieq,neq,QINF)
        !use pmgrid

        Implicit None
        integer,intent(in) :: ipar,isize,ieq(neq)
        double precision ,intent(out):: Qproj(isize,NXpm,NYpm,NZpm)
        double precision ,allocatable:: Qprojpriv(:,:,:,:)
        double precision ,intent(in) :: Qpar(isize,ipar),QpX(3,ipar),QINF(neq)
        integer,intent(in) :: Qprojtype(ipar)
        double precision   :: fx, fy, fz, f, x, y, z, projection_fun
        integer            :: inode, jnode, knode, i, j , k, nv, itype,nbj,nb,neq,ips,ipf

        integer            ::omp_get_max_threads,omp_get_num_threads
        !-->Projection function (TSC)
        Qproj = 0.d0
        !!$omp parallel private(nv,Qprojpriv,jnode,inode,itype,ips,ipf,i,j,k,f,fx,fy,x,y)
        !!allocate(Qprojpriv(NXpm,NYpm,NZpm,isize))
        !!Qprojpriv=0.d0
        !!$omp do
        do  nv = 1, ipar
            itype = Qprojtype(nv)
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((QpX(1,nv) - XMIN_pm) / DXpm) + 1
            jnode = int((QpX(2,nv) - YMIN_pm) / DYpm) + 1
            ! knode = int(ZVR(nv) / DZpm) 3D

            !    enddo

            !--We search the 4 nodes close to the particles
            !    do k = knode - 1, knode + 2 3D
            if(itype.eq.2) then
                ips =  0
                ipf =  1
            else
                ips =  1
                ipf =  2
            endif
           k = 1 !nbj
            do j = jnode - ips, jnode  + ipf
                do i = inode - ips, inode  + ipf

                    x  = (QpX(1,nv)- XMIN_pm - (i-1) * DXpm) / DXpm
                    fx = projection_fun(itype,x)

                    y  = (QpX(2,nv)- YMIN_pm - (j-1) * DYpm) / DYpm
                    fy = projection_fun(itype,y)

                    ! z  = (ZVR(nv) - (k-1) * DZpm) / DZpm 3D
                    ! fz = projection_fun(itype,z) * z     3D


                    f = fx * fy !* fz 3D
                    Qproj(ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k) +&
                        f * (QPar(ieq(1:neq-1),nv)-QINF(1:neq-1)*QPar(ieq(neq),nv))
                    Qproj(ieq(neq),i,j,k) = Qproj(ieq(neq),i,j,k) +&
                        f * (QPar(ieq(neq),nv)-QINF(neq))
                enddo
            enddo
        enddo
        !!$omp end do
        !!$omp critical
        ! !Qproj=Qproj+Qprojpriv
        !!$omp end critical
        !!deallocate(Qprojpriv)
        !!$omp end parallel

        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions

    End Subroutine project_particles_2D

    Subroutine project_particles_2D_vol(Qproj,Qpar,QpX,Qprojtype,ipar,isize,ieq,neq,QINF)
        !use pmgrid

        Implicit None
        integer,intent(in) :: ipar,isize,ieq(neq)
        double precision ,intent(out):: Qproj(NXpm,NYpm,NZpm)
        double precision ,allocatable:: Qprojpriv(:,:,:,:)
        double precision ,intent(in) :: Qpar(isize,ipar),QpX(3,ipar),QINF(neq)
        integer,intent(in) :: Qprojtype(ipar)
        double precision   :: fx, fy, fz, f, x, y, z, projection_fun
        integer            :: inode, jnode, knode, i, j , k, nv, itype,nbj,nb,neq,ips,ipf

        integer            ::omp_get_max_threads,omp_get_num_threads
        !-->Projection function (TSC)
        Qproj = 0.d0
        !!$omp parallel private(nv,Qprojpriv,jnode,inode,itype,ips,ipf,i,j,k,f,fx,fy,x,y)
        !!allocate(Qprojpriv(NXpm,NYpm,NZpm,isize))
        !!Qprojpriv=0.d0
        !!$omp do
        do  nv = 1, ipar
            itype = Qprojtype(nv)
            if (abs(Qpar(1,nv)).lt.0.5*DVpm) cycle 
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((QpX(1,nv) - XMIN_pm) / DXpm) + 1
            jnode = int((QpX(2,nv) - YMIN_pm) / DYpm) + 1
            ! knode = int(ZVR(nv) / DZpm) 3D

            !    enddo

            !--We search the 4 nodes close to the particles
            !    do k = knode - 1, knode + 2 3D
            if(itype.eq.2) then
                ips =  0
                ipf =  1
            else
                ips =  1
                ipf =  2
            endif
           k = 1 !nbj
            do j = jnode - ips, jnode  + ipf
                do i = inode - ips, inode  + ipf

                    x  = (QpX(1,nv)- XMIN_pm - (i-1) * DXpm) / DXpm
                    fx = projection_fun(itype,x)

                    y  = (QpX(2,nv)- YMIN_pm - (j-1) * DYpm) / DYpm
                    fy = projection_fun(itype,y)

                    ! z  = (ZVR(nv) - (k-1) * DZpm) / DZpm 3D
                    ! fz = projection_fun(itype,z) * z     3D


                    f = fx * fy !* fz 3D
                    Qproj(i,j,k) = Qproj(i,j,k) +&
                        f * (QPar(ieq(neq),nv))
                enddo
            enddo
        enddo
        !!$omp end do
        !!$omp critical
        ! !Qproj=Qproj+Qprojpriv
        !!$omp end critical
        !!deallocate(Qprojpriv)
        !!$omp end parallel

        !--After the projection all the values will be divided by Vol_pm to take into account the volume
        !  effect in the interpolation.We extend volume and density values to no-particle regions

    End Subroutine project_particles_2D_vol
    !--------------------------------------------------------------------------!
    !-->Subroutine project_vol2d                                               !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
    Subroutine project_vol2d(Qproj,isize,ieq,neq,iflag)
        !use pmgrid

        Implicit None
        integer,intent(in) :: isize,ieq(neq),iflag
        double precision   :: Qproj(isize,NXpm,NYpm,NZpm)
        double precision   :: fx, fy, fz, f, x, y, z, projection_fun,EPSVOLt
        integer            :: inode, jnode, knode, i, j , k, nv, itype,nbj,nb,neq,IDVPMt

        EPSVOLt=EPSVOL
        IDVPMt=IDVPM
        if(iflag.eq.1) IDVPMt=1
        !  if(iflag.eq.1) EPSVOLt=0.5
        if (IDVPMt.eq.0) then
          !$omp parallel private(i,j,k)
            k = 1 !nbj
          !$omp do
            do j =1, NYpm
                do i =1, NXpm
                    if (i.lt.NXs.or.i.gt.NXf.or.j.lt.NYs.or.j.ge.NYf) then
                        Qproj (ieq(1:neq-1),i,j,k) = 0.d0
                        cycle
                    endif
                   !if (Qproj(i,j,k,ieq(neq)).gt.EPSVOLt*DVpm) then
                        Qproj (ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k)/(Qproj(ieq(neq),i,j,k))
                   !else
                   !    Qproj (i,j,k,ieq(1:neq-1)) = 0.d0
                   !    Qproj (i,j,k,ieq(neq)) =DVpm
                   !endif
                enddo
            enddo
            !$omp end do
            !$omp end parallel
        else if (IDVPMt.eq.1) then
          !$omp parallel private(i,j,k)
            k = 1 !nbj
          !$omp do
            do j =1, NYpm
                do i =1, NXpm
                    if (i.lt.NXs.or.i.gt.NXf.or.j.lt.NYs.or.j.ge.NYf) then
                        Qproj (ieq(neq),i,j,k)=DVpm
                        Qproj (ieq(1:neq-1),i,j,k) = 0.d0
                        cycle
                    endif
                   !if (Qproj(i,j,k,ieq(neq)).gt.EPSVOLt*DVpm) then
                        Qproj (ieq(1:neq-1),i,j,k) = Qproj(ieq(1:neq-1),i,j,k)/ DVpm
                        Qproj (ieq(neq),i,j,k)=DVpm
                   !endif
                enddo
            enddo
          !$omp end do
          !$omp end parallel
        else
            write(*,*) 'WRONG IDVPM'
            STOP
        endif
    End Subroutine project_vol2d

    !--------------------------------------------------------------------------!
    !-->Subroutine project_particles_3D                                          !
    !   This subroutine projects particle values on the PM grid                !
    !   The values projected are :                                             !
    !      - Mass  -> becomes Density on the grid                              !
    !      - Volume                                                            !
    !      - Dilatation                                                        !
    !      - Phi                                                               !
    !      - PsiX , PsiY                                                       !
    !--------------------------------------------------------------------------!
End module projlib
!-----------------------------------------------------------------------------!
!-->Function projection_fun                                                   !
!   This subroutine defines the projection functions and finds the projected  !
!   value of the Particle Value on the PMgrid                                 !
!   Input :                                                                   !
!      itype       : type of projection function                              !
!      x           : position of projection
!-----------------------------------------------------------------------------!
double precision function projection_fun(itype, x)
    Implicit None
    double precision, intent(in) :: x
    integer, intent(in)          :: itype

    double precision             :: xabs

    xabs = abs(x)
    if(itype.eq.2) then
        if(xabs.gt.1) then
            projection_fun = 0.d0
        else  if(xabs.le.1) then
            projection_fun = 1.d0 - xabs
        else
            write(*,*) xabs
        endif
    else if (itype.eq.3) then
        !--Triangular-Shaped Cloud function
        if(xabs.gt.3.d0/2.d0)                   projection_fun = 0.d0
        if(xabs.ge.0.5d0.and.xabs.le.3.d0/2.d0) projection_fun = 0.5d0 * (-xabs + 3.d0/2.d0) ** 2
        if(xabs.ge.0.d0.and.xabs.lt.0.5d0)      projection_fun = 0.5d0 * ( xabs + 3.d0/2.d0) ** 2 - 3.d0/2.d0 * (xabs + 0.5d0) ** 2
    else if (itype.eq.4) then
        if(xabs.gt.2.d0)                        projection_fun = 0.d0
        if(xabs.ge.1.d0.and.xabs.le.2.d0)       projection_fun = 0.5d0 * (2.d0 - xabs)**2 * (1.d0 - xabs)
        if(xabs.lt.1.d0)                        projection_fun = 1 - 2.5d0 * xabs ** 2  + 1.5d0 * xabs ** 3
    else if (itype.eq.5) then 
        if(xabs.le.0.5d0)                       projection_fun = 1.d0 - xabs**2
        if(xabs.gt.0.5d0.and.xabs.le.1.5d0)     projection_fun = 0.5d0 * (1.d0 - xabs) * (2.d0 -xabs)
        if(xabs.gt.1.5d0)                       projection_fun = 0.d0
    else if (itype.eq.6) then 
        if(xabs.le.1.d0)                        projection_fun = 1.d0 - xabs**2
        if(xabs.gt.1.d0.and.xabs.le.2.d0)       projection_fun = 0.5d0 * (1.d0 - xabs) * (2.d0 -xabs)
        if(xabs.gt.2.d0)                        projection_fun = 0.d0
        
    else
        write(*,* ) 'No such projection function',itype
        STOP
    endif


End Function projection_fun
