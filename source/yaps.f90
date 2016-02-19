!This library solves the poisson problem using domain decomposition method
Module yapslib
    use projlib
    use pmlib
    use MPI

    private

    double precision,allocatable      :: SOL_pm_coarse(:,:,:,:),RHS_pm_coarse(:,:,:,:),SOL_pm_sample(:,:,:,:,:)
    double precision,allocatable      :: SOL_pm_sumsample(:,:,:,:)

    double precision             :: DXpm,DYpm,DZpm,DXpm_c,DYpm_c,DZpm_c
    integer                      :: NXpm,NYpm,NZpm,NXpm_c,NYpm_c,NZPm_c
    integer                      :: ibctyp,idvpm,epsvol,ND,iproj,ndumcell,npmsize
    double precision             :: XMIN_pm,YMIN_pm,ZMIN_pm,XMAX_pm,YMAX_pm,ZMAX_pm
    double precision             :: XMIN_pm_c,YMIN_pm_c,ZMIN_pm_c,XMAX_pm_c,YMAX_pm_c,ZMAX_pm_c
    double precision             :: MACH
    double precision,allocatable :: Xbound_bl(:,:)
    integer                      :: BLOCKS,NXB,NYB,NBI,NBJ,NBK,NB,i,j,k,NXbl,NYbl,NN(3),NN_bl(6)
    integer,allocatable          :: NNbl(:,:),NNbl_bl(:,:),NN_coarse_map(:,:),map_nodes(:,:,:,:),nnb(:)
    double precision             :: projection_fun,fx,fy,fz,f,xc,yc,zc,X(3),addsol,starttime,endtime
    integer                      :: ic,jc,kc,inode,jnode,knode
    integer                      :: i_nb,j_nb,k_nb
    integer                      :: NXs,NYs,NZs,NXf,NYf,NZf,ib,jb,kb,ibj,jbj,kbj,ixs,ixf,jxs,jxf,izs,izf
    integer                      :: nc,NN_map(6),isubtrackt,node,neq,isizex,isizey,isizez,nbc


    double precision,pointer                 :: SOL_pm_bl(:,:,:,:), RHS_pm_bl(:,:,:,:),QP(:,:),XP(:,:)

    double precision,allocatable             :: BBound(:,:,:)


    integer  :: status(MPI_STATUS_SIZE),source,ierr,my_rank,np,mat4,mat5,dest
    character *25                :: outfil1,outfil2

    public  :: yaps2d,yaps3d
contains
    include 'yaps2d.f90'
    include 'yaps3d.f90'
    Subroutine mpimat2(mat2,nsize1,nsize2)
        use MPI

        integer ierr
        integer::typelist(2)
        integer ::imat(2),mat(2),start(2)
        integer ::istart
        integer ::orig1,orig2,nsize1,nsize2,mat2
        !allocate(struct%AS_ij(nsize,nsize))

        mat(1) = nsize1
        mat(2) = nsize2
        imat = mat
        start(1) = 0
        start(2) = 0

        typelist(1)=MPI_DOUBLE_PRECISION
            typelist(2)=MPI_DOUBLE_PRECISION

                !write (*,*) nsize
                call MPI_TYPE_CREATE_SUBARRAY(2,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat2,ierr)
                call MPI_TYPE_COMMIT(mat2,ierr)

            end subroutine mpimat2

            Subroutine mpimat4(mat4,nsize1,nsize2,nsize3,nsize4)
                use MPI

                integer, intent(in)::nsize1,nsize2,nsize3,nsize4
                integer ierr
                integer ::imat(4),mat(4),start(4)
                integer ::mat4
                !allocate(struct%AS_ij(nsize,nsize))
                imat(1)=nsize1
                imat(2)=nsize2
                imat(3)=nsize3
                imat(4)=nsize4
                mat(1)=nsize1
                mat(2)=nsize2
                mat(3)=nsize3
                mat(4)=nsize4
                start(1)=0
                start(2)=0
                start(3)=0
                start(4)=0

                !write (*,*) nsize
                call MPI_TYPE_CREATE_SUBARRAY(4,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat4,ierr)
                call MPI_TYPE_COMMIT(mat4,ierr)

            End subroutine mpimat4


            Subroutine mpimat5(mat5,nsize1,nsize2,nsize3,nsize4,nsize5)
                use MPI

                integer, intent(in)::nsize1,nsize2,nsize3,nsize4,nsize5
                integer ierr
                integer ::imat(5),mat(5),start(5)
                integer ::mat5
                !allocate(struct%AS_ij(nsize,nsize))
                imat(1)=nsize1
                imat(2)=nsize2
                imat(3)=nsize3
                imat(4)=nsize4
                imat(5)=nsize5
                mat(1)=nsize1
                mat(2)=nsize2
                mat(3)=nsize3
                mat(4)=nsize4
                mat(5)=nsize5
                start(1)=0
                start(2)=0
                start(3)=0
                start(4)=0
                start(5)=0

                !write (*,*) nsize
                call MPI_TYPE_CREATE_SUBARRAY(5,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat5,ierr)
                call MPI_TYPE_COMMIT(mat5,ierr)

            End subroutine mpimat5

        End Module yapslib

        Subroutine writesol_bl(RHS_pm,SOL_pm,Dpm,outfil,Xbound,NN_bl,NN)

            character*25      :: outfil
            integer           :: i,j,k,NN_bl(6),NN(3)
            double precision  :: RHS_pm(7,NN(1),NN(2),NN(3))
            double precision  :: SOL_pm(7,NN(1),NN(2),NN(3))
            double precision  :: Dpm(3),Xbound(6)
            double precision  :: XPM,YPM,ZPM,velx,vely,velz


            open(1,file=outfil)
            WRITE(1,'(a150)')'VARIABLES = "X" "Y" "U" "V" "VORTZ""DIL""PSI""PHI"'
            WRITE(1,*)'ZONE T='//outfil//' I=',NN_bl(4)-NN_bl(1)+1,' J=',NN_bl(5) - NN_bl(2) + 1,' F=POINT'
            do j=NN_bl(2),NN_bl(5)
                do i=NN_bl(1),NN_bl(4)
                    ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
                    ! do j=1,NYpm
                    !   do i=1,NXpm
                    XPM=Xbound(1)+(I-1)*Dpm(1)
                    YPM=Xbound(2)+(J-1)*Dpm(2)
                    velx=0
                    vely=0
                    !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                    !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)
                    !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                    !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)

                    WRITE(1,'(10(e28.17,1x))')XPM,YPM,velx,vely,-RHS_pm(1,I,J,1),RHS_pm(2,I,J,1),SOL_pm(1,I,J,1),SOL_pm(2,I,J,1)

                enddo
            enddo
            close(1)
            !   ---FOR PLOTTING PURPOSES ONLY
            call system('preplot '//outfil//' >/dev/null')
            call system('rm '//outfil)



        End Subroutine writesol_bl

        Subroutine writesol_bl_3d(outfil,Dpm,Xbound,NN_bl,NN)

            character*25      :: outfil
            integer           :: i,j,k,NN_bl(6),NN(3)
            double precision  :: Dpm(3),Xbound(6)
            double precision  :: XPM,YPM,ZPM

            open(1,file=outfil)
            WRITE(1,'(a250)')'VARIABLES = "X" "Y" "Z"'
            WRITE(1,*)'ZONE T='//outfil//' I=',NN(1),' J=',NN(2),' K=',NN(3),&
           ' F=POINT'
            do k=1,NN(3)
               do j=1,NN(2)
                   do i=1,NN(1)
                       ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
                       ! do j=1,NYpm
                       !   do i=1,NXpm
                       XPM=Xbound(1)+(I-1)*Dpm(1)
                       YPM=Xbound(2)+(J-1)*Dpm(2)
                       ZPM=Xbound(3)+(K-1)*Dpm(3)
                       !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                       !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)
                       !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                       !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)
               
                      WRITE(1,'(3(e28.17,1x))')XPM,YPM,ZPM  
                   enddo
               enddo
            enddo
            close(1)
            !   ---FOR PLOTTING PURPOSES ONLY
            call system('preplot '//outfil//' >/dev/null')
            call system('rm '//outfil)



        End Subroutine writesol_bl_3d

        Subroutine writegrow(RHS_pm,SOL_pm,Dpm,outfil,Xbound,NN_bl,NN)

            character*25      :: outfil
            integer           :: i,j,k,NN_bl(6),NN(3)
            double precision  :: RHS_pm(3,NN(1),NN(2),NN(3))
            double precision  :: SOL_pm(3,NN(1),NN(2),NN(3))
            double precision  :: Dpm(3),Xbound(6)
            double precision  :: XPM,YPM,ZPM,velx,vely,velz


            open(1,file=outfil)
            WRITE(1,'(a150)')'VARIABLES = "X" "Y""VORTZ""DIL""PSI""PHI"'
            WRITE(1,*)'ZONE T='//outfil//' I=',NN(1),' J=',NN(2),' F=POINT'
            do j=1,NN(2)
                do i=1,NN(1)
                    ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
                    ! do j=1,NYpm
                    !   do i=1,NXpm
                    XPM=Xbound(1)+(I-1)*Dpm(1)
                    YPM=Xbound(2)+(J-1)*Dpm(2)
                    !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                    !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)

                    WRITE(1,'(10(e28.17,1x))')XPM,YPM,-RHS_pm(1,I,J,1),RHS_pm(2,I,J,1),SOL_pm(1,I,J,1),SOL_pm(2,I,J,1)

                enddo
            enddo
            close(1)
            !   ---FOR PLOTTING PURPOSES ONLY
            call system('preplot '//outfil//' >/dev/null')
            call system('rm '//outfil)



        End Subroutine writegrow


      
        Subroutine writegrow_3d(RHS_pm,SOL_pm,Dpm,outfil,Xbound,NN_bl,NN)

            character*25      :: outfil
            integer           :: i,j,k,NN_bl(6),NN(3)
            double precision  :: RHS_pm(3,NN(1),NN(2),NN(3))
            double precision  :: SOL_pm(3,NN(1),NN(2),NN(3))
            double precision  :: Dpm(3),Xbound(6)
            double precision  :: XPM,YPM,ZPM,velx,vely,velz


            open(1,file=outfil)
            WRITE(1,'(a150)')'VARIABLES = "X" "Y" "Z" "VORTZ""DIL""PSI""PHI"'
            WRITE(1,*)'ZONE T='//outfil//' I=',NN(1),' J=',NN(2),' K=',NN(3),&
           ' F=POINT'
            do k=1,NN(3)
             do j=1,NN(2)
                do i=1,NN(1)
                    ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
                    ! do j=1,NYpm
                    !   do i=1,NXpm
                    XPM=Xbound(1)+(I-1)*Dpm(1)
                    YPM=Xbound(2)+(J-1)*Dpm(2)
                    ZPM=Xbound(3)+(K-1)*Dpm(3)
                    !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                    !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)

                    WRITE(1,'(11(e28.17,1x))')XPM,YPM,ZPM,-RHS_pm(1,I,J,K),-RHS_pm(2,I,J,K),SOL_pm(1,I,J,K),SOL_pm(2,I,J,K)

                enddo
            enddo
            enddo
            close(1)
            !   ---FOR PLOTTING PURPOSES ONLY
            call system('preplot '//outfil//' >/dev/null')
            call system('rm '//outfil)



        End Subroutine writegrow_3d
