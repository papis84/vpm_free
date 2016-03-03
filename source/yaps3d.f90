    Subroutine yaps3d(DSOL_pm,DRHS_pm,Xbound_bl,Xbound_coarse,Dpm_fine,Dpm_coarse,NNbl,NNbl_bl,&
            NN_coarse,NN_bl_coarse,ND,BLOCKS,ibctyp,neqs,neqf,nc,NBI,NBJ,NBK,nb_i,nb_j,nb_k,ireturn,&
            iyntree,ilevmax,npmsize)

        use projlib
        use pmlib
        use MPI
        Implicit None
        integer,intent(in)            :: ibctyp,neqs,neqf,nc,ireturn,iyntree,ilevmax,npmsize
        integer,intent(in)            :: ND,BLOCKS,NNbl(3,BLOCKS),NNBl_bl(6,BLOCKS)
        integer,intent(in)            :: NN_coarse(3),NN_bl_coarse(6),nb_i,nb_j,nb_k,NBI,NBJ,NBK

        double precision ,intent(in)            :: Xbound_bl(6,BLOCKS),Xbound_coarse(6)
        double precision ,intent(in)            :: Dpm_fine(3),Dpm_coarse(3)
        double precision,intent(inout),target   :: DSOL_pm(:,:,:,:), DRHS_pm(:,:,:,:)

        double precision,allocatable :: SOL_pm_tmp(:,:,:,:),RHS_pm_tmp(:,:,:,:)
        double precision             :: Xbound_tmp(6)
        integer                      :: NN_tmp(3),NN_bl_tmp(6),iynbc,iface12,iface34,iface56,ibound,itree,lmax

        !Assign variables
        nullify(SOL_pm_bl,RHS_pm_bl)
        SOL_pm_bl=>DSOL_pm;RHS_pm_bl=>DRHS_pm
        !normally QP,XP not needed
        nullify(QP)
        nullify(XP)

        call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

        !The equations to be solved should be from 1 to neqf
        neq = neqf
        nb = my_rank+1
        Xbound_tmp(1:6) = Xbound_bl(1:6,nb)
        NN_tmp(1:3)     = NNbl(1:3,nb)
        NN_bl_tmp(1:6)  = NNbl_bl(1:6,nb)

        if (my_rank.eq.0)starttime = MPI_WTIME()
        !initialise solution vector
        SOL_pm_bl=0.d0
        iynbc=1!infinite domain bc's
        itree=iyntree !tree algorithm for sources
        lmax=ilevmax!maximum level
        ! 1 is the Nblocks not needed needs fix
        call pmesh(SOL_pm_bl,RHS_pm_bl,QP,XP,&
                   Xbound_tmp,Dpm_fine,NN_tmp,NN_bl_tmp,ND,1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)

        !write(outfil1,'(a9,i2.2)') 'blockgrow',nb
        !call  writesol_bl_3d(outfil1,Dpm_fine,Xbound_tmp,NN_bl_tmp,NN_tmp)
        !---Block definitions

        !Define coarse pm from which values will be interpolated for the final solve
        NXpm_c = NN_coarse(1);NYpm_c=NN_coarse(2);NZpm_c=NN_coarse(3)



          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pmesh',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
        allocate(SOL_pm_coarse(npmsize,NXpm_c,NYpm_c,NZpm_c),RHS_pm_coarse(npmsize,NXpm_c,NYpm_c,NZpm_c))


        !call definevort(RHS_pm_coarse,MACH,Xbound_coarse,Dpm_coarse,NN_coarse,NN_bl_coarse,iproj)


        !SOL_pm_coarse=0.d0

        !call pmesh(SOL_pm_coarse,RHS_pm_coarse,QP,XP,velphix_coarse,velphiy_coarse,velphiz_coarse,&
            !                        velvrx_coarse,velvry_coarse,velvrz_coarse,Xbound_coarse,DPm_coarse,NN_coarse,NN_bl_coarse,ND,1,ibctyp,1,1,1,1)

        !if (my_rank.eq.0) then
            ! write(*,*) 'Coarse block dimensions'
            ! write(*,*) NN_bl_coarse(1:6)
            !outfil2='coarse'
            !call writesol(RHS_pm_coarse,SOL_pm_coarse,velphix_coarse,velphiy_coarse,velvrx_coarse,velvry_coarse,Dpm_coarse,outfil2,Xbound_coarse,NN_bl_coarse,NN_coarse)DU91-W2-250
        !endif



        !Variable map_nodes if 1 then the  block nb contains inode,jnode of coarse grid   
        allocate(map_nodes(NXpm_c,NYpm_c,NZpm_c,BLOCKS),nnb(BLOCKS))
        allocate(SOL_pm_sample(npmsize,NXpm_c,NYpm_c,NZpm_c,BLOCKS))
        !NN_coarse_map matches the the boundaries of the fine domains with the boundaries at he coarse domains
        allocate (NN_coarse_map(6,BLOCKS))
        NN_coarse_map=-1 !to help debbuging...
        SOL_pm_sample=0.d0
        map_nodes=0

        nb     = my_rank+1 !(j_nb-1)*NBB + i_nb
       
        !This is done by ALL blocks for ALL blocks.The main reason is NN_coarse_map.An alternate would be
        !broadcasting this information

        if(my_rank.eq.0)starttime = MPI_WTIME()
        do nbc=1,BLOCKS
            NXs    = 1!NNbl_bl(1,nb)
            NXf    = NNbl(1,nbc)!NNbl_bl(4,nb)

            NYs = 1
            NYf = NNbl(2,nbc)

            NZs = 1
            NZf = NNbl(3,nbc)

           !NXs    = NNbl_bl(1,nbc)
           !NXf    = NNbl_bl(4,nbc)

           !NYs = NNbl_bl(2,nbc)
           !NYf = NNbl_bl(5,nbc)

           !NZs = NNbl_bl(3,nbc)
           !NZf = NNbl_bl(6,nbc)
            !---------------------------------------------------------------------------------
            do k=NZs,NZf,nc
                 do i=NXs,NXf,nc
                      do j=NYs,NYf,nc
                          X(1) = Xbound_bl(1,nbc) + (i-1)*Dpm_fine(1)
                          X(2) = Xbound_bl(2,nbc) + (j-1)*Dpm_fine(2)
                          X(3) = Xbound_bl(3,nbc) + (k-1)*Dpm_fine(3)
                          inode = int((X(1) - Xbound_coarse(1))/ Dpm_coarse(1)) +1
                          jnode = int((X(2) - Xbound_coarse(2))/ Dpm_coarse(2)) +1
                          knode = int((X(3) - Xbound_coarse(3))/ Dpm_coarse(3)) +1
                          if (i.eq.NXs) NN_coarse_map(1,nbc) = inode
                          if (i.eq.NXf) NN_coarse_map(4,nbc) = inode
                          if (j.eq.NYs) NN_coarse_map(2,nbc) = jnode
                          if (j.eq.NYf) NN_coarse_map(5,nbc) = jnode
                          if (k.eq.NZs) NN_coarse_map(3,nbc) = knode
                          if (k.eq.NZf) NN_coarse_map(6,nbc) = knode
                          !each rank fills the SOL_pm_sample matrix(Sampling of fine solution)
                          if(nbc.eq.nb) then 
                            SOL_pm_sample(neqs:neqf,inode,jnode,knode,nb) =SOL_pm_bl(neqs:neqf,i,j,k)
                          endif
                          map_nodes(inode,jnode,knode,nbc)= 1
                      enddo
                  enddo
                enddo
           enddo

          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pm_sample1',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        if(my_rank.eq.0)starttime = MPI_WTIME()
        !--BCAST Sol_pm_sample
        nb=my_rank+1
        if (my_rank.eq.0) then
            allocate(SOL_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3)))
            allocate(RHS_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3)))
            SOL_pm_tmp(:,:,:,:)= SOL_pm_sample(:,:,:,:,nb)
            NN_map(:)=NN_coarse_map(:,nb)
            !calculate the RHS of the poisson problem using the sampled values
            call calc_laplacian_coarse_3d(SOL_pm_tmp,RHS_pm_tmp,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_map,neqs,neqf)
            RHS_pm_coarse=0
            RHS_pm_coarse=RHS_pm_coarse + RHS_pm_tmp
            !0 proccesor gathers the rhs using sampled values from all cpu's
            do nbc=2,BLOCKS
                call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3))
                source = nbc-1
                call MPI_RECV(SOL_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(RHS_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                SOL_pm_sample(:,:,:,:,nbc)=SOL_pm_tmp(:,:,:,:)
                RHS_pm_coarse=RHS_pm_coarse + RHS_pm_tmp
                call MPI_TYPE_FREE(mat4,ierr)
            enddo
            deallocate(SOL_pm_tmp,RHS_pm_tmp)
        else
            allocate(SOL_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3)))
            allocate(RHS_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3)))
            SOL_pm_tmp(:,:,:,:)= SOL_pm_sample(:,:,:,:,nb)
            NN_map(:)=NN_coarse_map(:,nb)
            call calc_laplacian_coarse_3d(SOL_pm_tmp,RHS_pm_tmp,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_map,neqs,neqf)
            dest=0
            call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3))
            call MPI_SEND(SOL_pm_tmp,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(RHS_pm_tmp,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
            call MPI_TYPE_FREE(mat4,ierr)
            deallocate(SOL_pm_tmp,RHS_pm_tmp)
        endif
        !We broadcast SOL_pm_sample(which is the addition of Sampled solution from proccesors to all cpu's
          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pm_sample11',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
        if(my_rank.eq.0)starttime = MPI_WTIME()
        call mpimat5(mat5,npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3),BLOCKS)
        call MPI_BCAST(SOL_pm_sample,1,mat5,0,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat5,ierr)

        !We broadcast RHS based on the sampled solutions of all ranks
        call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),NN_coarse(3))
        call MPI_BCAST(RHS_pm_coarse,1,mat4,0,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat4,ierr)

          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pm_bcast',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        !add Sampling blocks
       !allocate (SOL_pm_sumsample(NN_coarse(1),NN_coarse(2),NN_coarse(3),7))
       !SOL_pm_sumsample=0.d0
       !do nbj=1,BLOCKS
       !    SOL_pm_sumsample(:,:,:,:)=SOL_pm_sumsample(:,:,:,:)+SOL_pm_sample(:,:,:,:,nbj)
       !enddo
        !call calc_laplacian(SOL_pm_sumsample,RHS_pm_coarse,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_bl_coarse)

        !if(my_rank.eq.0)then
        !    outfil2='coarse_sumsample_grow'
        !   call  writesol_bl_3d(outfil2,Dpm_coarse,Xbound_coarse,NN_bl_coarse,NN_coarse)
        !endif

        ! if (my_rank.eq.1)starttime = MPI_WTIME()
        ! if(my_rank.eq.1)endtime = MPI_WTIME()
        ! if(my_rank.eq.1) write(*,*)'Laplace=',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        !if (my_rank.eq.1)starttime = MPI_WTIME()
        SOL_pm_coarse=0.d0
        iynbc=1
        itree=0!iyntree
        lmax=ilevmax
        if(my_rank.eq.0)starttime = MPI_WTIME()
        if (itree.eq.0) levmax=1
        call pmesh(SOL_pm_coarse,RHS_pm_coarse,QP,XP,&
             Xbound_coarse,DPm_coarse,NN_coarse,NN_bl_coarse,ND,1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)
        ! if(my_rank.eq.1)endtime = MPI_WTIME()
        ! if(my_rank.eq.1) write(*,*)'Poisson Coarse=',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pmesh_coarse',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
      ! if (my_rank.eq.0) then
      !      outfil2='coarse_laplacian'
      !      call writegrow_3d(RHS_pm_coarse,SOL_pm_coarse,Dpm_coarse,outfil2,Xbound_coarse,NN_bl_coarse,NN_coarse)
      ! endif
        ! Interpolate to all blocks
        
        ! At this point we calculate boundary conditions for each block and solve it 
        if(my_rank.eq.0)starttime = MPI_WTIME()
        nb  = my_rank+1
        NXs = NNbl_bl(1,nb)
        NXf = NNbl_bl(4,nb)

        NYs = NNbl_bl(2,nb)
        NYf = NNbl_bl(5,nb)

        NZs = NNbl_bl(3,nb)
        NZf = NNbl_bl(6,nb)

        isizex  = NNbl(1,nb)
        isizey  = NNbl(2,nb)
        isizez  = NNbl(3,nb)
        !NODS COUNT FROM FACE1,FACE2-->FACE3,FACE4--->FACE5,FACE6
        !iface12 corresponds to planes normal x faces (1-2)
        iface12 = isizey * isizez
        !iface34 corresponds to planes normal y faces (3-4)
        iface34 = isizex * isizez
        !iface56 corresponds to planes normal z faces (5-6)
        iface56 = isizex * isizey
        ibound = 2*(iface12 + iface34 + iface56)
        !BBounds are the matrix which contain the boundary conditions of each blokc
        allocate(BBound(ibound,neq,BLOCKS))
        call mapnodes_bl
        !---------------------------------------------------------------------------------
        !Xs boundary ,9 point stencil
        !face 1

          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pm_bc1',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
          if(my_rank.eq.0)starttime = MPI_WTIME()
        i=NXs
        do k = NZs, NZf
          do j = NYs, NYf
            node   = 0
            node=node + (k-1)*isizey + j
            call interp_stencil
          enddo
        enddo
        !face 2
        i=NXf
        do k = NZs, NZf
          do j = NYs, NYf
            node   = iface12
            node=node + (k-1)*isizey + j
            call interp_stencil
          enddo
        enddo

        !face 3
        j=NYs
        do k = NZs, NZf
          do i =NXs+1, NXf -1
            node   = 2*iface12
            node=node + (k-1)*isizex + i
            call interp_stencil
          enddo
        enddo

        !face 4
        j=NYf
        do k = NZs, NZf
          do i =NXs+1, NXf -1
            node   = 2*iface12 + iface34
            node=node + (k-1)*isizex + i
            call interp_stencil
          enddo
        enddo

        !face 5
        k=NZs
        do j = NYs+1, NYf-1
          do i =NXs+1, NXf -1
            node   = 2*iface12 + 2*iface34
            node=node + (j-1)*isizex + i
            call interp_stencil
          enddo
        enddo

        k=NZf
        do j = NYs+1, NYf-1
          do i =NXs+1, NXf -1
            node   = 2*iface12 + 2*iface34 +iface56
            node=node + (j-1)*isizex + i
            call interp_stencil
          enddo
        enddo
        !--Finished setting coarse boundary conditions
        !---------------------------------------------------------------------------------

        Xbound_tmp(1:6) = Xbound_bl(1:6,nb)
        NN_tmp(1:3)     = NNbl(1:3,nb)
        NN_bl_tmp(1:6)  = NNbl_bl(1:6,nb)



        !in case we want to add something external in bc's
        deallocate(NN_coarse_map,SOL_pm_coarse,RHS_pm_coarse,map_nodes,nnb,SOL_pm_sample,BBound)
        if (ireturn.eq.1) then
            nullify(SOL_pm_bl,RHS_pm_bl)
            return
        endif
          if(my_rank.eq.0)endtime = MPI_WTIME()
          if(my_rank.eq.0) write(199,*)'pm_bc2',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
        ! write(*,*) 'Solving for Block',nb
        !iynbc=0 means that the bc's of the poisson solver are already defined
        iynbc=0
        itree=0
        lmax=0
        if (my_rank.eq.0) starttime=MPI_WTIME()
        call pmesh(SOL_pm_bl,RHS_pm_bl,QP,XP,&
                   Xbound_tmp,Dpm_fine,NN_tmp,NN_bl_tmp,ND,1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)
        if(my_rank.eq.0)endtime = MPI_WTIME()
        if(my_rank.eq.0) write(199,*)'pmesh_final',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        ! write(outfil1,'(a5,i2.2)') 'block',nb
        ! call writesol_bl_3d(RHS_pm_bl,SOL_pm_bl,Dpm_fine,outfil1,Xbound_tmp,NN_bl_tmp,NN_tmp)

        !allocate(SOL_pm_er(NN_fine(1),NN_fine(2),1,7))
        !allocate(velvrx_tmp(NXpm,NYpm,NZpm) ,velvry_tmp(NXpm,NYpm,NZpm),velvrz_tmp(NXpm,NYpm,NZpm))
        ! write(*,*) NXf-NXs+1
        ! ixs=(i_nb-1)*(NXf-NXs)+NN_bl_fine(1)
        ! jxs=(j_nb-1)*(NYf-NXs)+NN_bl_fine(2)

        ! ixf= ixs + (NXf-NXs+1)-1
        ! jxf= jxs + (NYf-NYs+1)-1
        ! write(*,*) 'BLOCKmap',ixs,ixf,jxs,jxf
        ! SOL_pm_er(ixs:ixf,jxs:jxf,1,1)= SOL_pm(NN_bl(1):NN_bl(4),NN_bl(2):NN_bl(5),1,1)
        ! velvrx_tmp(ixs:ixf,jxs:jxf,1) = velvrx(NN_bl(1):NN_bl(4),NN_bl(2):NN_bl(5),1)
        ! velvry_tmp(ixs:ixf,jxs:jxf,1) = velvry(NN_bl(1):NN_bl(4),NN_bl(2):NN_bl(5),1)


        !SOL_pm_er  = abs(SOL_pm_er -SOL_pm_fine)
        !velphix_fine=0
        !velphiy_fine=0
        !velvrx_tmp = abs(velvrx_tmp - velvrx_fine)
        !velvry_tmp = abs(velvry_tmp - velvry_fine)
        !write(outfil1,'(a5,i2.2)') 'error',nc
        !call writesol(RHS_pm_fine,SOL_pm_er,velphix_fine,velphiy_fine,velvrx_tmp,velvry_tmp,Dpm_fine,outfil1,Xbound_fine,NN_bl_fine,NN_fine)




    contains

        Subroutine mapnodes_bl

            integer :: i_nb,j_nb,icnb,jcnb,inode1,jnode1,nbc,nb_plus,nb_minus,i_check,j_check,ik,jk,kk
            integer :: nod,ibound,isizex,isizey,isizez,i_plus,j_plus,k_plus,i_minus,j_minus,k_minus
            integer :: itarget ,jtarget, ktarget, idir,nodstart
            integer :: iface12 ,iface34 ,iface56,m
            double precision :: X1(3)
            double precision,allocatable :: Bound_sol(:,:)
            integer  :: status(MPI_STATUS_SIZE),source,ierr,mat2,dest

            NXs = NNbl_bl(1,nb)
            NXf = NNbl_bl(4,nb)

            NYs = NNbl_bl(2,nb)
            NYf = NNbl_bl(5,nb)

            NZs = NNbl_bl(3,nb)
            NZf = NNbl_bl(6,nb)

            nod=0
            nb = my_rank+1
            isizex  = NNbl(1,nb)
            isizey  = NNbl(2,nb)
            isizez  = NNbl(3,nb)
            !NODS COUNT FROM FACE1,FACE2-->FACE3,FACE4--->FACE5,FACE6
            !iface12 corresponds to planes normal x faces (1-2)
            iface12 = isizey * isizez
            !iface34 corresponds to planes normal y faces (3-4)
            iface34 = isizex * isizez
            !iface56 corresponds to planes normal z faces (5-6)
            iface56 = isizex * isizey
            !total nodes of faces of each block
            ibound = 2*(iface12 + iface34 + iface56)
           !ibound = isizex*isizey*isizez

            allocate(Bound_sol(ibound,neq))
            BBound=0.d0
            call mpimat2(mat2,ibound,neq)
            !(ATTENTION BOX EDGES ARE COMING ALWAYS FROM FACES 1-2 and the remaining from 3-4)
            !NODS COUNT FROM FACE1,FACE2-->FACE3,FACE4--->FACE5,FACE6--Corners are at 1-2
            ! at each face first x,second y,third z (node convention numbering)
            !Bound_sol is with respect to the receiving.
            !SEND TO THE RIGHT BLOCKS
            !RIGHT SWEEP in I direction,node blocks are only done once.
            do idir =1,17
                ik=0;jk=0;kk=0;Bound_sol=0
                if (idir.eq.1) then
                    !SENDING TO I+1
                    !sending from face 2 of block my_rank to face 12 of i+1 block
                    i= NXf
                    nodstart=0
                    do k=NZs,NZf
                       do j=NYs,NYf
                          ktarget = k
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-3 to face 1-3
                    j=NYs
                    nodstart = 2 * iface12
                    do k=NZs,NZf
                       do i=NXf,NNbl(1,nb)
                          ktarget = k
                          itarget = i - NXf + NXs  !to start from NXs + 1(NXs before)
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-4 to face 1-4
                    j=NYf
                    nodstart = 2 * iface12 + iface34
                    do k=NZs,NZf
                       do i=NXf,NNbl(1,nb)
                          ktarget = k
                          itarget = i - NXf + NXs  !to start from NXs + 1(NXs before)
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-5 to face 1-5
                    k=NZs
                    nodstart =  2 * iface12 + 2 * iface34
                    do j=NYs,NYf
                       do i=NXf,NNbl(1,nb)
                          jtarget = j
                          itarget = i -  NXf  + NXs  !to start from NXs + 1(NXs before)
                          nod     = nodstart+ (jtarget - 1) * isizex + itarget
 
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-6 to face 1-6
                    k=NZf
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,NYf
                       do i=NXf,NNbl(1,nb)
                          jtarget = j
                          itarget = i - NXf + NXs !to start from NXs + 1(NXs before)
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                         
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1
                endif
                if (idir.eq.2)  then
                    !SENDING TO J+1
                    !sending from face 4 of block my_rank to face 2 of j+1 block
                    nodstart=2*iface12
                    j= NYf
                    do k=NZs,NZf
                       do i=NXs,NXf
                          ktarget = k
                          itarget = i 
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 1-4 to face 1-3
                    i=NXs
                    nodstart = 0
                    do k=NZs,NZf
                       do j=NYf,NNbl(2,nb)
                          ktarget = k
                          jtarget = j - NYf + NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 2-4 to face 2-3
                    i=NXf
                    nodstart = iface12
                    do k=NZs,NZf
                       do j=NYf,NNbl(2,nb)
                          ktarget = k
                          jtarget = j - NYf + NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 5-4 to face 5-3
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,NXf
                          jtarget = j - NYf + NYs
                          itarget = i
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 6-4 to face 6-3
                    k=NZf
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,NXf
                          jtarget = j - NYf + NYs
                          itarget = i
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    jk=1
                endif
                if (idir.eq.3) then
                    !SENDING TO K + 1
                    !sending from face 6 of block my_rank to face 5 of j+1 block
                    
                    k= NZf
                    nodstart = 2 * iface12 + 2 * iface34
                    do j=NYs,NYf
                       do i=NXs,NXf
                          jtarget = j
                          itarget = i 
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 1-6 to face 1-5
                    i=NXs
                    nodstart = 0
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,NYf
                          ktarget = k - NZf + NZs
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 2-6 to face 2-5
                    i=NXf
                    nodstart = iface12
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,NYf
                          ktarget = k - NZf + NZs
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 3-6 to face 3-5
                    j=NYs
                    nodstart = 2 * iface12
                    do k=NZf,NNbl(3,nb)
                       do i=NXs,NXf
                          ktarget = k - NZf + NZs
                          itarget = i
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                     

                    !from face 4-6 to face 4-5
                    j=NYf
                    nodstart = 2 * iface12 + iface34
                    do k=NZf,NNbl(3,nb)
                       do i=NXs,NXf
                          ktarget = k - NZf + NZs
                          itarget = i
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    kk=1
                endif
                if (idir.eq.4) then
                    !SENDING TO I+1,J+1
                    !sending from face 2 to face 1(/ diagonal) 
                    i= NXf
                    nodstart = 0
                    do k=NZs,NZf
                       do j=NYf,NNbl(2,nb)
                          ktarget = k
                          jtarget = j - NYf + NYs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-4 to face 1-3
                    j= NYf
                    nodstart = 2 * iface12
                    do k=NZs,NZf
                       do i=NXf,NNbl(1,nb)
                          ktarget = k
                          itarget = i - NXf + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 2-4(on side 5) to  face 1-3 (on side5)
                    k=NZs
                    nodstart = 2 * iface12 + 2* iface34
                    do j=NYf,NNbl(2,nb)
                       do i=NXf,NNbl(1,nb)
                          jtarget = j - NYf  + NYs 
                          itarget = i - NXf  + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 2-4(on side 6) to  face 1-3 (on side6)
                    k=NZf
                    nodstart = 2 * iface12 + 2* iface34 + iface56
                    do j=NYf,NNbl(2,nb)
                       do i=NXf,NNbl(1,nb)
                          jtarget = j - NYf  + NYs
                          itarget = i - NXf  + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=1;kk=0
                endif
                if (idir.eq.5) then 
                    !SENDING TO I+1,J-1
                    !sending from face 2 to face 1(\ diagonal) 
                    i= NXf
                    nodstart = 0
                    do k=NZs,NZf
                       do j=NYs,1,-1
                          ktarget = k
                          jtarget = j - NYs +NYf!Lower diagonal starts at NYf
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-3 to face 1-4
                    j= NYs
                    nodstart = 2 * iface12 + iface34
                    do k=NZs,NZf
                       do i=NXf,NNbl(1,nb)
                          ktarget = k
                          itarget = i - NXf + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 2-3 to face 1-4 side5
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34
                    do j=NYs,1,-1
                       do i=NXf,NNbl(1,nb)
                          jtarget = j - NYs + NYf  !Lower diagonal starts at NYf
                          itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 2-3 to face 1-4 side6
                    k=NZf
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,1,-1
                       do i=NXf,NNbl(1,nb)
                          jtarget = j - NYs + NYf !Lower diagonal starts at NYf
                          itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=-1;kk=0
                endif
                if (idir.eq.6) then
                    !sending to I+1,K+1
                    !sending from face 2 to face 1(/ diagonal in XZ plane) 
                    i= NXf
                    nodstart = 0
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,NYf
                          ktarget = k - NZf + NZs
                          jtarget = j              
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-6 to face 1-5
                    k= NZf
                    nodstart = 2 * iface12 + 2 * iface34
                    do j=NYs,NYf
                       do i=NXf,NNbl(1,nb)
                          jtarget = j
                          itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 2-6(on side 3) to  face 1-5 (on side3)
                    j=NYs
                    nodstart = 2 * iface12 
                    do k=NZf,NNbl(3,nb)
                       do i=NXf,NNbl(1,nb)
                          ktarget = k - NZf  + NZs 
                          itarget = i - NXf  + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 2-4(on side 6) to  face 1-3 (on side6)
                    j=NYf
                    nodstart = 2 * iface12 +  iface34
                    do k=NZf,NNbl(3,nb)
                       do i=NXf,NNbl(1,nb)
                          ktarget = k - NZf  + NZs 
                          itarget = i - NXf  + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=0;kk=1
                endif
                if (idir.eq.7) then 
                    !SENDING TO I+1,K-1
                    !sending from face 2 to face 1(\ diagonal in XZ plane) 
                    i= NXf
                    nodstart = 0
                    do k=NZs,1,-1
                       do j=NYs,NYf
                          ktarget = k - NZs + NZf
                          jtarget = j 
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 2-5 to face 1-6
                    k= NZs
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,NYf
                       do i=NXf,NNbl(1,nb)
                          jtarget = j
                          itarget = i - NXf + NXs !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 2-5 to face 1-6 side3
                    j=NYs
                    nodstart = 2 * iface12 
                    do k=NZs,1,-1
                       do i=NXf,NNbl(1,nb)
                          ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                          itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 2-5 to face 1-6 side4
                    j=NYf
                    nodstart = 2 * iface12 +  iface34 
                    do k=NZs,1,-1
                       do i=NXf,NNbl(1,nb)
                          ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                          itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=0;kk=-1
                endif
                if (idir.eq.8) then
                    !sending to J+1,K+1
                    !sending from face 4 to face 3(/ diagonal in YZ plane) 
                    j= NYf
                    nodstart = 2 * iface12
                    do k=NZf,NNbl(3,nb)
                       do i=NXs,NXf
                          ktarget = k - NZf + NZs
                          itarget = i              
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 4-6 to face 3-5
                    k= NZf
                    nodstart = 2 * iface12 + 2 * iface34
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,NXf
                          jtarget = j - NYf + NYs
                          itarget = i   !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 4-6(on side 1) to  face 3-5 (on side1)
                    i=NXs
                    nodstart = 0
                    do k=NZf,NNbl(3,nb)
                       do j=NYf,NNbl(2,nb)
                          ktarget = k - NZf  + NZs 
                          jtarget = j - NYf  + NYs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                   !from face 4-6(on side 2) to  face 3-5 (on side2)
                    i=NXf
                    nodstart = iface12
                    do k=NZf,NNbl(3,nb)
                       do j=NYf,NNbl(2,nb)
                          ktarget = k - NZf  + NZs 
                          jtarget = j - NYf  + NYs !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=0;jk=1;kk=1
                endif
                if (idir.eq.9) then 
                    !SENDING TO J+1,K-1
                    !sending from face 4 to face 3(\ diagonal in YZ plane) 
                    j= NYf
                    nodstart = 2 * iface12
                    do k=NZs,1,-1
                       do i=NXs,NXf
                          ktarget = k - NZs + NZf
                          itarget = i 
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 4-5 to face 3-6
                    k= NZs
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,NXf
                          jtarget = j - NYf + NYs
                          itarget = i   !Upper diagonal starts at NYs
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 4-5 to face 3-6 side1
                    i=NXs
                    nodstart = 0
                    do k=NZs,1,-1
                       do j=NYf,NNbl(2,nb)
                          ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                          jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 2-5 to face 1-6 side4
                    i=NXf
                    nodstart = iface12
                    do k=NZs,1,-1
                       do j=NYf,NNbl(2,nb)
                          ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                          jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=0;jk=1;kk=-1
                endif
                if (idir.eq.10) then
                   !node overlapping boxes sending to  I+1,J+1,K+1
                   !face 2 to face 1
                   i=NXf
                   nodstart=0
                   do k=NZf,NNbl(3,nb)
                      do j=NYf,NNbl(2,nb)
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 4 to face 3
                   j=NYf
                   nodstart=2 * iface12
                   do k=NZf,NNbl(3,nb)
                      do i=NXf,NNbl(1,nb)
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 6 to face 5
                   k=NZf
                   nodstart=2 * iface12 + 2 * iface34
                   do j=NYf,NNbl(2,nb)
                      do i=NXf,NNbl(1,nb)
                         jtarget = j - NYf + NYs  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=1;jk=1;kk=1
                endif
                if (idir.eq.11) then
                   !node overlapping boxes sending to  I+1,J+1,K-1
                   !face 2 to face 1
                   i=NXf
                   nodstart=0
                   do k=NZs,1,-1
                      do j=NYf,NNbl(2,nb)
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 4 to face 3
                   j=NYf
                   nodstart=2 * iface12
                   do k=NZs,1,-1
                      do i=NXf,NNbl(1,nb)
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 5 to face 6
                   k=NZs
                   nodstart=2 * iface12 + 2 * iface34 + iface56
                   do j=NYf,NNbl(2,nb)
                      do i=NXf,NNbl(1,nb)
                         jtarget = j - NYf + NYs  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=1;jk=1;kk=-1
                endif
                if (idir.eq.12) then
                   !node overlapping boxes sending to  I+1,J-1,K-1
                   !face 2 to face 1
                   i=NXf
                   nodstart=0
                   do k=NZs,1,-1
                      do j=NYs,1,-1
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         jtarget = j - NYs + NYf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 3 to face 4
                   j=NYs
                   nodstart=2 * iface12  + iface34
                   do k=NZs,1,-1
                      do i=NXf,NNbl(1,nb)
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 5 to face 6
                   k=NZs
                   nodstart=2 * iface12 + 2 * iface34 + iface56
                   do j=NYs,1,-1
                      do i=NXf,NNbl(1,nb)
                         jtarget = j - NYs + NYf  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=1;jk=-1;kk=-1
                endif
                if (idir.eq.13) then
                   !node overlapping boxes sending to  I+1,J-1,K+1
                   !face 2 to face 1
                   i=NXf
                   nodstart=0
                   do k=NZf,NNbl(3,nb)
                      do j=NYs,1,-1
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         jtarget = j - NYs + NYf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 3 to face 4
                   j=NYs
                   nodstart=2 * iface12  + iface34
                   do k=NZf,NNbl(3,nb)
                      do i=NXf,NNbl(1,nb)
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 6 to face 5
                   k=NZf
                   nodstart=2 * iface12 + 2 * iface34
                   do j=NYs,1,-1
                      do i=NXf,NNbl(1,nb)
                         jtarget = j - NYs + NYf  !Lower diagonal starts at NYf
                         itarget = i - NXf + NXs  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=1;jk=-1;kk=1
                endif
                if (idir.eq.14) then
                   !node overlapping boxes sending to  I-1,J+1,K+1
                   !face 1 to face 2
                   i=NXs
                   nodstart=iface12
                   do k=NZf,NNbl(3,nb)
                      do j=NYf,NNbl(2,nb)
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 4 to face 3
                   j=NYf
                   nodstart=2 * iface12
                   do k=NZf,NNbl(3,nb)
                      do i=NXs,1,-1
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 6 to face 5
                   k=NZf
                   nodstart=2 * iface12 + 2 * iface34
                   do j=NYf,NNbl(2,nb)
                      do i=NXs,1,-1
                         jtarget = j - NYf + NYs  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=-1;jk=1;kk=1
                endif
                if (idir.eq.15) then
                   !node overlapping boxes sending to  I-1,J+1,K-1
                   !face 1 to face 2
                   i=NXs
                   nodstart=iface12
                   do k=NZs,1,-1
                      do j=NYf,NNbl(2,nb)
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         jtarget = j - NYf + NYs  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 4 to face 3
                   j=NYf
                   nodstart=2 * iface12
                   do k=NZs,1,-1
                      do i=NXs,1,-1
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 5 to face 6
                   k=NZs
                   nodstart=2 * iface12 + 2 * iface34 + iface56
                   do j=NYf,NNbl(2,nb)
                      do i=NXs,1,-1
                         jtarget = j - NYf + NYs  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=-1;jk=1;kk=-1
                endif
                if (idir.eq.16) then
                   !node overlapping boxes sending to  I-1,J-1,K-1
                   !face 1 to face 2
                   i=NXs
                   nodstart=iface12
                   do k=NZs,1,-1
                      do j=NYs,1,-1
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         jtarget = j - NYs + NYf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 3 to face 4
                   j=NYs
                   nodstart=2 * iface12  + iface34
                   do k=NZs,1,-1
                      do i=NXs,1,-1
                         ktarget = k - NZs + NZf  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 5 to face 6
                   k=NZs
                   nodstart=2 * iface12 + 2 * iface34 + iface56
                   do j=NYs,1,-1
                      do i=NXs,1,-1
                         jtarget = j - NYs + NYf  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=-1;jk=-1;kk=-1
                endif
                if (idir.eq.17) then
                   !node overlapping boxes sending to  I-1,J-1,K+1
                   !face 1 to face 2
                   i=NXs
                   nodstart=iface12
                   do k=NZf,NNbl(3,nb)
                      do j=NYs,1,-1
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         jtarget = j - NYs + NYf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizey + jtarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 3 to face 4
                   j=NYs
                   nodstart=2 * iface12  + iface34
                   do k=NZf,NNbl(3,nb)
                      do i=NXs,1,-1
                         ktarget = k - NZf + NZs  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (ktarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   !face 6 to face 5
                   k=NZf
                   nodstart=2 * iface12 + 2 * iface34
                   do j=NYs,1,-1
                      do i=NXs,1,-1
                         jtarget = j - NYs + NYf  !Lower diagonal starts at NYf
                         itarget = i - NXs + NXf  !Upper diagonal starts at NYs
                         nod     = nodstart + (jtarget - 1) * isizex + itarget
                         Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                      enddo
                   enddo
                   ik=-1;jk=-1;kk=1
                endif

                i_plus  = nb_i + ik
                j_plus  = nb_j + jk !because in equation (j-1)
                k_plus  = nb_k + kk
                nb_plus = (k_plus - 1) * NBJ * NBI + (j_plus - 1) * NBI + i_plus
                if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1.and.&
                    k_plus.le.NBK.and.k_plus.ge.1) then
                    dest = nb_plus-1 ! because nb = rank-1
                    !write(*,*) nb-1,'sending to ',dest,idir
                    call MPI_SEND(Bound_sol,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
                    i_minus  = nb_i - ik
                    j_minus  = nb_j - jk  
                    k_minus  = nb_k - kk
                    nb_minus = (k_minus - 1) * NBJ * NBI + (j_minus - 1) * NBI + i_minus
                    if (i_minus.ge.1.and.j_minus.ge.1.and.j_minus.le.NBJ.and.i_minus.le.NBI.and.&
                        k_minus.ge.1.and.k_minus.le.NBK)  then
                        source=nb_minus-1
                        !  write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_minus) = BBound(:,:,nb_minus) + Bound_sol(:,:)
                    endif
                else
                    i_minus  = nb_i - ik
                    j_minus  = nb_j - jk  
                    k_minus  = nb_k - kk
                    nb_minus = (k_minus - 1) * NBJ * NBI + (j_minus - 1) * NBI + i_minus

                    if (i_minus.ge.1.and.j_minus.ge.1.and.j_minus.le.NBJ.and.i_minus.le.NBI.and.&
                        k_minus.ge.1.and.k_minus.le.NBK)  then
                        source=nb_minus-1
                        !  write(*,*) nb-1,' received from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_minus) = BBound(:,:,nb_minus) + Bound_sol(:,:)
                    endif
                endif
            enddo
            !LEFT SWEEP in I

          !write(*,*) 'RIGHT SWEEP OK'
            do idir =1,9
                ik=0;jk=0;kk=0;Bound_sol=0
                if (idir.eq.1) then
                    !SENDING TO I-1

                    !sending from face 1 of block my_rank to face 2 of i+1 block
                    i= NXs
                    nodstart = iface12
                    do k=NZs,NZf
                       do j=NYs,NYf
                          ktarget = k
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-3 to face 2-3
                    j=NYs
                    nodstart = 2 * iface12
                    do k=NZs,NZf
                       do i=NXs,1,-1
                          ktarget = k
                          itarget = i - NXs + NXf  !to start from NXs + 1(NXs before)
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-4 to face 2-4
                    j=NYf
                    nodstart = 2 * iface12 + iface34
                    do k=NZs,NZf
                       do i=NXs,1,-1
                          ktarget = k
                          itarget = i - NXs + NXf !to start from NXs + 1(NXs before)
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-5 to face 2-5
                    k=NZs
                    nodstart =  2 * iface12 + 2 * iface34
                    do j=NYs,NYf
                       do i=NXs,1,-1
                          jtarget = j
                          itarget = i - NXs + NXf  !to start from NXs + 1(NXs before)
                          nod     = nodstart + (jtarget - 1) * isizex + itarget
 
                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-6 to face 2-6
                    k=NZf
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,NYf
                       do i=NXs,1,-1
                          jtarget = j
                          itarget = i - NXs + NXf !to start from NXs + 1(NXs before)
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1
                endif
                if (idir.eq.2)  then
                    !SENDING TO J-1
                    !sending from face 3 of block my_rank to face 4 of j-1 block
                    j= NYs
                    nodstart=2*iface12 + iface34
                    do k=NZs,NZf
                       do i=NXs,NXf
                          ktarget = k
                          itarget = i 
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 1-3 to face 1-4
                    i=NXs
                    nodstart = 0
                    do k=NZs,NZf
                       do j = NYs ,1,-1
                          ktarget = k
                          jtarget = j - NYs + NYf
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 2-3 to face 2-4
                    i=NXf
                    nodstart =  iface12
                    do k=NZs,NZf
                       do j = NYs,1,-1
                          ktarget = k
                          jtarget = j - NYs + NYf
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 5-3 to face 5-4
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34
                    do j = NYs,1,-1
                       do i = NXs, NXf-1
                          jtarget = j - NYs  + NYf
                          itarget = i
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 6-3 to face 6-4
                    k=NZf
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j = NYs,1,-1
                       do i=NXs,NXf
                          jtarget = j - NYs + NYf
                          itarget = i
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    jk=1
                endif
                if (idir.eq.3) then
                    !SENDING TO K - 1
                    !sending from face 5 of block my_rank to face 6 of k-1 block
                    
                    k= NZs
                    nodstart= 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,NYf
                       do i=NXs,NXf
                          jtarget = j
                          itarget = i 
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 1-5 to face 1-6
                    i=NXs
                    nodstart = 0
                    do k=NZs,1,-1
                       do j=NYs,NYf
                          ktarget = k - NZs + NZf
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    !from face 2-5 to face 2-6
                    i=NXf
                    nodstart = iface12
                    do k = NZs, 1, -1
                       do j = NYs, NYf
                          ktarget = k - NZs + NZf
                          jtarget = j
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face 3-5 to face 3-6
                    j=NYs
                    nodstart = 2 * iface12
                    do k = NZs, 1, -1
                       do i=NXs,NXf 
                          ktarget = k - NZs + NZf 
                          itarget = i
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                     

                    !from face 4-5 to face 4-6
                    j=NYf
                    nodstart = 2 * iface12 + iface34
                    do k = NZs , 1, -1
                       do i = NXs, NXf
                          ktarget = k - NZs + NZf 
                          itarget = i
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    kk=1
                endif
                if (idir.eq.4) then
                    !SENDING TO I-1,J-1 (left /)
                    !sending from face 1 of block my_rank to face 2 of i+1 block
                    i = NXs
                    nodstart = iface12
                    do k=NZs,NZf
                       do j=NYs,1,-1
                          ktarget = k
                          jtarget = j - NYs + NYf
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-3 to face 2-4
                    j=NYs
                    nodstart = 2 * iface12 +iface34
                    do k=NZs,NZf
                       do i=NXs,1,-1
                          ktarget = k
                          itarget = i - NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(1-3) to face(2-4 on side 5
                    k=NZs
                    nodstart = 2 * iface12 + 2*iface34
                    do j=NYs,1,-1
                       do i=NXs,1,-1
                          jtarget = j - NYs + NYf  !to start from NXf - 1 and backwards
                          itarget = i - NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(1-3) to face(2-4 on side 6
                    k=NZf
                    nodstart = 2 * iface12 + 2*iface34 +iface56
                    do j=NYs,1,-1
                       do i=NXs,1,-1
                          jtarget = j - NYs + NYf  !to start from NXf - 1 and backwards
                          itarget = i - NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=1;kk=0
                 endif

                if (idir.eq.5) then
                    !SENDING TO I-1,J+1 (left \)
                    !sending from face 1 of block my_rank to face 2 of i+1 block
                    i= NXs
                    nodstart = iface12
                    do k=NZs,NZf
                       do j=NYf,NNbl(2,nb)
                          ktarget = k
                          jtarget = j - NYf + NYs
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-4 to face 2-3
                    j=NYf
                    nodstart = 2 * iface12
                    do k=NZs,NZf
                       do i=NXs,1,-1
                          ktarget = k
                          itarget = i - NXs + NXf !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-4 to face 2-3 on side 5
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34 
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,1,-1
                          jtarget = j - NYf + NYs 
                          itarget = i - NXs + NXf   !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 1-4 to face 2-3 on side 6
                    k=NZf
                    nodstart = 2 * iface12 + 2*iface34 +iface56
                    do j=NYf,NNbl(2,nb)
                       do i=NXs,1,-1
                          jtarget = j - NYf + NYs 
                          itarget = i - NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=-1;kk=0
                 endif
                if (idir.eq.6) then
                    !SENDING TO I-1,K-1 (left / in XZ plane)
                    !sending from face 1 of block my_rank to face 2 of i+1 block
                    i = NXs
                    nodstart = iface12
                    do k=NZs,1,-1
                       do j=NYs,NYf
                          ktarget = k - NZs + NZf
                          jtarget = j 
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-5 to face 2-6
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,NYf
                       do i=NXs,1,-1
                          jtarget = j
                          itarget = i - NXs +NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(1-5) to face(2-6 on side 3
                    j=NYs
                    nodstart = 2 * iface12 
                    do k=NZs,1,-1
                       do i=NXs,1,-1
                          ktarget = k - NZs +NZf !to start from NXf - 1 and backwards
                          itarget = i - NXs +NXf !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(1-5) to face(2-6 on side 4
                    j=NYf
                    nodstart = 2 * iface12 + iface34 
                    do k=NZs,1,-1
                       do i=NXs,1,-1
                          ktarget = k -  NZs + NZf  !to start from NXf - 1 and backwards
                          itarget = i -  NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=0;kk=1
                 endif

                if (idir.eq.7) then
                    !SENDING TO I-1,K+1 (left \ in XZ plane)
                    !sending from face 1 of block my_rank to face 2 of i+1 block
                    i= NXs
                    nodstart = iface12
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,NYf
                          ktarget = k - NZf + NZs
                          jtarget = j 
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-6 to face 2-5
                    k=NZf
                    nodstart = 2 * iface12 + 2* iface34
                    do j=NYs,NYf
                       do i=NXs,1,-1
                          jtarget = j
                          itarget = i - NXs + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 1-6 to face 2-5 on side 3
                    j=NYs
                    nodstart = 2 * iface12
                    do k=NZf,NNbl(3,nb)
                       do i=NXs,1,-1
                          ktarget = k - NZf  + NZs 
                          itarget = i - NXs  + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 1-6 to face 2-5 on side 4
                    j=NYf
                    nodstart = 2 * iface12 + iface34
                    do k=NZf ,NNbl(3,nb)
                       do i=NXs,1,-1
                          ktarget = k - NZf  + NZs 
                          itarget = i - NXs  + NXf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=1;jk=0;kk=-1
                 endif
                if (idir.eq.8) then
                    !SENDING TO J-1,K-1 (left / in YZ plane)
                    !sending from face 3 of block my_rank to face 4 of i+1 block
                    j = NYs
                    nodstart = 2*iface12 + iface34
                    do k=NZs,1,-1
                       do i=NXs,NXf
                          ktarget = k - NZs + NZf
                          itarget = i 
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 3-5 to face 4-6
                    k=NZs
                    nodstart = 2 * iface12 + 2 * iface34 + iface56
                    do j=NYs,1,-1
                       do i=NXs,NXf
                          jtarget = j - NYs + NYf
                          itarget = i 
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(3-5) to face(4-6 on side 1
                    i=NXs
                    nodstart = 0
                    do k=NZs,1,-1
                       do j=NYs,1,-1
                          ktarget = k - NZs +NZf !to start from NXf - 1 and backwards
                          jtarget = j - NYs +NYf !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    !from face(3-5) to face(4-6 on side 2
                    i=NXf
                    nodstart = iface12
                    do k=NZs,1,-1
                       do j=NYs,1,-1
                          ktarget = k - NZs +NZf !to start from NXf - 1 and backwards
                          jtarget = j - NYs +NYf !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=0;jk=1;kk=1
                 endif

                if (idir.eq.9) then
                    !SENDING TO J-1,K+1 (left \ in XZ plane)
                    !sending from face 3 of block my_rank to face 4 of i+1 block
                    j= NYs
                    nodstart = 2 * iface12 + iface34
                    do k=NZf,NNbl(3,nb)
                       do i=NXs,NXf
                          ktarget = k - NZf + NZs
                          itarget = i
                          nod     = nodstart + (ktarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 3-6 to face 4-5
                    k=NZf
                    nodstart = 2 * iface12 + 2* iface34
                    do j=NYs,1,-1
                       do i=NXs,NXf
                          jtarget = j - NYs + NYf
                          itarget = i                                                    
                          nod     = nodstart + (jtarget - 1) * isizex + itarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ! from face 3-6 to face 4-5 on side 1
                    i=NXs
                    nodstart = 0
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,1,-1
                          ktarget = k - NZf  + NZs 
                          jtarget = j - NYs  + NYf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo

                    ! from face 1-6 to face 2-5 on side 4
                    i=NXf
                    nodstart = iface12 
                    do k=NZf,NNbl(3,nb)
                       do j=NYs,1,-1
                          ktarget = k - NZf  + NZs 
                          jtarget = j - NYs  + NYf  !to start from NXf - 1 and backwards
                          nod     = nodstart + (ktarget - 1) * isizey + jtarget

                          Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,k)
                       enddo
                    enddo
                    ik=0;jk=1;kk=-1
                 endif


                !fix Nodes at corners we want to put all contributions(from y nodes,corners node exist twice)
                ! in (1,isizex,isizex+1,2*isizex)

                i_minus = nb_i - ik
                j_minus = nb_j - jk  
                k_minus = nb_k - kk  
                nb_minus = (k_minus - 1) * NBJ * NBI + (j_minus - 1) * NBI + i_minus
                if (i_minus.ge.1.and.j_minus.ge.1.and.i_minus.le.NBI.and.j_minus.le.NBJ.and.&
                    k_minus.ge.1.and.k_minus.le.NBK) then
                    dest = nb_minus-1 ! because nb = rank-1
                    !write(*,*) nb-1,'sending to ',dest
                    call MPI_SEND(Bound_sol,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
                    i_plus  = nb_i + ik
                    j_plus  = nb_j + jk !because in equation (j-1)
                    k_plus  = nb_k + kk
                    nb_plus = (k_plus - 1) * NBJ * NBI + (j_plus - 1) * NBI + i_plus
                    if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1.and.&
                        k_plus.le.NBK.and.k_plus.ge.1)  then
                        source=nb_plus-1
                        !    write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_plus) = BBound(:,:,nb_plus) + Bound_sol(:,:)
                    endif
                else
                    i_plus  = nb_i + ik
                    j_plus  = nb_j + jk !because in equation (j-1)
                    k_plus  = nb_k + kk
                    nb_plus = (k_plus - 1) * NBJ * NBI + (j_plus - 1) * NBI + i_plus
                    if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1.and.&
                        k_plus.le.NBK.and.k_plus.ge.1)  then
                        source=nb_plus-1
                        !    write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_plus) = BBound(:,:,nb_plus) + Bound_sol(:,:)
                    endif
                endif

            enddo

            call MPI_TYPE_FREE(mat2,ierr)




        End Subroutine mapnodes_bl

!-----
        Subroutine interp_stencil
            double precision :: projection_fun,addlocal,add(neq),add_sample(neq)
            integer          :: i_nb,j_nb,k_nb

            iproj=3
            X(1) = Xbound_bl(1,nb) + (i-1)*Dpm_fine(1)
            X(2) = Xbound_bl(2,nb) + (j-1)*Dpm_fine(2)
            X(3) = Xbound_bl(3,nb) + (k-1)*Dpm_fine(3)

            !   SOL_pm_bl(i,j,1,1,nb)=SOL_pm_fine(inode,jnode,1,1)
            !inode =int(((X(1) - Xbound_coarse(1)) / Dpm_coarse(1))) + 1
            !jnode =int(((X(2) - Xbound_coarse(2)) / Dpm_coarse(2))) + 1
            inode =int(nint((X(1) - Xbound_coarse(1)) / Dpm_coarse(1))) + 1
            jnode =int(nint((X(2) - Xbound_coarse(2)) / Dpm_coarse(2))) + 1
            knode =int(nint((X(3) - Xbound_coarse(3)) / Dpm_coarse(3))) + 1
            !--We search the 4 nodes close to the particles
            nnb=0
            addlocal=0
            do nbc = 1,BLOCKS
                if(map_nodes(inode-1,jnode-1,knode-1,nbc).ne.1.or.&
                        map_nodes(inode+1,jnode+1,knode+1,nbc).ne.1)then
                    nnb(nbc)=1
                endif
            enddo
            
            do kc = knode-1, knode + 1
               do jc = jnode-1 , jnode + 1
                   do ic = inode-1 , inode + 1
                       add=0.d0
                       add_sample=0.d0
                       do nbc=1,BLOCKS
                           add_sample(neqs:neqf) = SOL_pm_sample(neqs:neqf,ic,jc,kc,nbc)+add_sample(neqs:neqf)
                           if(nnb(nbc).eq.1) then
                               !add fine local contribution when available
                               add(neqs:neqf)=SOL_pm_sample(neqs:neqf,ic,jc,kc,nbc) + add(neqs:neqf)
                           endif
                       enddo
                       xc  = (X(1)- Xbound_coarse(1) - (ic-1) * Dpm_coarse(1)) / Dpm_coarse(1)
                       fx = projection_fun(iproj,xc)
               
                       yc  = (X(2)- Xbound_coarse(2) - (jc-1) * Dpm_coarse(2)) / Dpm_coarse(2)
                       fy = projection_fun(iproj,yc)
               
                       zc  = (X(3)- Xbound_coarse(3) - (kc-1) * Dpm_coarse(3)) / Dpm_coarse(3)
                       fz = projection_fun(iproj,zc)

                       f = fx * fy * fz
                       !SOL_pm_coa
                       SOL_pm_bl(neqs:neqf,i,j,k) = SOL_pm_bl(neqs:neqf,i,j,k)+ &
                           f * (SOL_pm_coarse(neqs:neqf,ic,jc,kc)-add_sample(neqs:neqf)+&
                           add(neqs:neqf))
                   enddo
               enddo
            enddo
!          write(*,*) 'ok',node
            !add local data in the blocks that sampled data was interpolated
            do k_nb = nb_k-1, nb_k +1
               do j_nb = nb_j-1, nb_j+1
                   do i_nb = nb_i-1, nb_i+1
                       if (i_nb.lt.1.or.i_nb.gt.NBI.or.j_nb.lt.1.or.j_nb.gt.NBJ.or.k_nb.lt.1.or.k_nb.gt.NBK)cycle
                       nbc = (k_nb-1)*NBI*NBJ + (j_nb-1)*NBI + i_nb
                       if (nnb(nbc).eq.1.or.nbc.eq.nb) cycle
                       SOL_pm_bl(neqs:neqf,i,j,k)= SOL_pm_bl(neqs:neqf,i,j,k) +&
                                                    BBound(node,neqs:neqf ,nbc)
                   enddo
               enddo
            enddo


        End Subroutine interp_stencil


    End Subroutine yaps3d

   Subroutine calc_laplacian_coarse_3d(SOL_pm,RHS_pm,NN,NN_bl,Dpm,NN_map,neqs,neqf)
       Implicit none

       integer,intent(in) :: NN(3),NN_bl(6),NN_map(6),neqs,neqf
       double precision,intent(out)   :: RHS_pm(neqf,NN(1),NN(2),NN(3))
       double precision,intent(in)    :: SOL_pm(neqf,NN(1),NN(2),NN(3)),Dpm(3)
       integer                        :: i,j,k



       RHS_pm=0.d0
       do k=NN_map(3)+1,NN_map(6)-1
          do j=NN_map(2)+1,NN_map(5)-1
             do i=NN_map(1)+1,NN_map(4)-1
                RHS_pm(neqs:neqf,i,j,k)=&
               (SOL_pm(neqs:neqf,i+1,j,k) - 2*SOL_pm(neqs:neqf,i,j,k) +&
                SOL_pm(neqs:neqf,i-1,j,k))/Dpm(1)**2 + &
               (SOL_pm(neqs:neqf,i,j+1,k) - 2*SOL_pm(neqs:neqf,i,j,k) +&
                SOL_pm(neqs:neqf,i,j-1,k))/Dpm(2)**2 + &
               (SOL_pm(neqs:neqf,i,j,k+1) - 2*SOL_pm(neqs:neqf,i,j,k) +&
                SOL_pm(neqs:neqf,i,j,k-1))/Dpm(3)**2  
     
              enddo
          enddo
       enddo


        End Subroutine calc_laplacian_coarse_3d

