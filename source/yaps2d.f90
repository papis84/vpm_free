    Subroutine yaps2d(DSOL_pm,DRHS_pm,Xbound_bl,Xbound_coarse,Dpm_fine,Dpm_coarse,NNbl,NNbl_bl,&
            NN_coarse,NN_bl_coarse,ND,BLOCKS,ibctyp,neqs,neqf,nc,NBI,NBJ,nb_i,nb_j,ireturn,&
            iyntree,ilevmax,npmsize)

        use projlib
        use pmlib
        use MPI
        Implicit None
        integer,intent(in)            :: ibctyp,neqs,neqf,nc,ireturn,iyntree,ilevmax,npmsize
        integer,intent(in)            :: ND,BLOCKS,NNbl(3,BLOCKS),NNBl_bl(6,BLOCKS)
        integer,intent(in)            :: NN_coarse(3),NN_bl_coarse(6),nb_i,nb_j,NBI,NBJ

        double precision ,intent(in)            :: Xbound_bl(6,BLOCKS),Xbound_coarse(6)
        double precision ,intent(in)            :: Dpm_fine(3),Dpm_coarse(3)
        double precision,intent(inout),target   :: DSOL_pm(:,:,:,:), DRHS_pm(:,:,:,:)

        double precision,allocatable :: SOL_pm_tmp(:,:,:,:),RHS_pm_tmp(:,:,:,:)
        double precision             :: Xbound_tmp(6)
        integer                      :: NN_tmp(3),NN_bl_tmp(6),iynbc,itree,lmax

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
        iynbc=1 !infinite domain bc's
        itree=iyntree !tree algorithm for sources
        lmax=ilevmax!maximum level
        ! 1 is the Nblocks not needed needs fix
        call pmesh(SOL_pm_bl,RHS_pm_bl,QP,XP,&
                   Xbound_tmp,Dpm_fine,NN_tmp,NN_bl_tmp,ND,1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)

        !write(outfil1,'(a9,i2.2)') 'blockgrow',nb
        !call writegrow(RHS_pm_bl,SOL_pm_bl,Dpm_fine,outfil1,Xbound_tmp,NN_bl_tmp,NN_tmp)
        !---Block definitions


        !4efine coarse pm from which values will be interpolated for the final solve
        NXpm_c = NN_coarse(1);NYpm_c=NN_coarse(2);NZpm_c=NN_coarse(3)



        allocate(SOL_pm_coarse(npmsize,NXpm_c,NYpm_c,NZpm_c),RHS_pm_coarse(npmsize,NXpm_c,NYpm_c,NZpm_c))


        !call definevort(RHS_pm_coarse,MACH,Xbound_coarse,Dpm_coarse,NN_coarse,NN_bl_coarse,iproj)


        !SOL_pm_coarse=0.d0

        !call pmesh(SOL_pm_coarse,RHS_pm_coarse,QP,XP,velphix_coarse,velphiy_coarse,velphiz_coarse,&
            !                        velvrx_coarse,velvry_coarse,velvrz_coarse,Xbound_coarse,DPm_coarse,NN_coarse,NN_bl_coarse,ND,1,ibctyp,1,1,1,1)

        !if (my_rank.eq.0) then
            ! write(*,*) 'Coarse block dimensions'
            ! write(*,*) NN_bl_coarse(1:6)
            !outfil2='coarse'
            !call writesol(RHS_pm_coarse,SOL_pm_coarse,velphix_coarse,velphiy_coarse,velvrx_coarse,velvry_coarse,Dpm_coarse,outfil2,Xbound_coarse,NN_bl_coarse,NN_coarse)
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
        isizey = NNbl_bl(4,nb) - NNbl_bl(1,nb)+ 1
        isizex = NNbl_bl(5,nb) - NNbl_bl(2,nb)+ 1
        !BBounds are the matrix which contain the boundary conditions of each blokc
        allocate(BBound(2*(isizex+isizey),neq,BLOCKS))

        !This is done by ALL blocks for ALL blocks.The main reason is NN_coarse_map.An alternate would be
        !broadcasting this information
        do nbc=1,BLOCKS
            NXs    = 1!NNbl_bl(1,nb)
            NXf    = NNbl(1,nbc)!NNbl_bl(4,nb)

            NYs = 1
            NYf = NNbl(2,nbc)

            NZs = 1
            NZf = NNbl(3,nbc)
            !---------------------------------------------------------------------------------
              do i=NXs,NXf,nc
                   do j=NYs,NYf,nc
                       X(1) = Xbound_bl(1,nbc) + (i-1)*Dpm_fine(1)
                       X(2) = Xbound_bl(2,nbc) + (j-1)*Dpm_fine(2)
                       X(3) = Xbound_bl(3,nbc) + (k-1)*Dpm_fine(3)
                       inode = int(nint((X(1) - Xbound_coarse(1))/ Dpm_coarse(1))) +1
                       jnode = int(nint((X(2) - Xbound_coarse(2))/ Dpm_coarse(2))) +1
                       if (i.eq.NXs) NN_coarse_map(1,nbc) = inode
                       if (i.eq.NXf) NN_coarse_map(4,nbc) = inode
                       if (j.eq.NYs) NN_coarse_map(2,nbc) = jnode
                       if (j.eq.NYf) NN_coarse_map(5,nbc) = jnode
                       !each rank fills the SOL_pm_sample matrix(Sampling of fine solution)
                       if(nbc.eq.nb) SOL_pm_sample(neqs:neqf,inode,jnode,1,nb) =SOL_pm_bl(neqs:neqf,i,j,1)
                       map_nodes(inode,jnode,1,nbc)= 1
                   enddo
               enddo
           enddo

        !--BCAST Sol_pm_sample
        nb=my_rank+1
        if (my_rank.eq.0) then
            allocate(SOL_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),1))
            allocate(RHS_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),1))
            SOL_pm_tmp(:,:,:,:)= SOL_pm_sample(:,:,:,:,nb)
            NN_map(:)=NN_coarse_map(:,nb)
            !calculate the RHS of the poisson problem using the sampled values
            call calc_laplacian_coarse(SOL_pm_tmp,RHS_pm_tmp,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_map,neqs,neqf)
            RHS_pm_coarse=0
            RHS_pm_coarse=RHS_pm_coarse + RHS_pm_tmp
            !0 proccesor gathers the rhs using sampled values from all cpu's
            do nbc=2,BLOCKS
                call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),1)
                source = nbc-1
                call MPI_RECV(SOL_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                call MPI_RECV(RHS_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                SOL_pm_sample(:,:,:,:,nbc)=SOL_pm_tmp(:,:,:,:)
                RHS_pm_coarse=RHS_pm_coarse + RHS_pm_tmp
                call MPI_TYPE_FREE(mat4,ierr)
            enddo
            deallocate(SOL_pm_tmp,RHS_pm_tmp)
        else
            allocate(SOL_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),1))
            allocate(RHS_pm_tmp(npmsize,NN_coarse(1),NN_coarse(2),1))
            SOL_pm_tmp(:,:,:,:)= SOL_pm_sample(:,:,:,:,nb)
            NN_map(:)=NN_coarse_map(:,nb)
            call calc_laplacian_coarse(SOL_pm_tmp,RHS_pm_tmp,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_map,neqs,neqf)
            dest=0
            call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),1)
            call MPI_SEND(SOL_pm_tmp,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
            call MPI_SEND(RHS_pm_tmp,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
            call MPI_TYPE_FREE(mat4,ierr)
            deallocate(SOL_pm_tmp,RHS_pm_tmp)
        endif
        !We broadcast SOL_pm_sample(which is the addition of Sampled solution from proccesors to all cpu's
        call mpimat5(mat5,npmsize,NN_coarse(1),NN_coarse(2),1,BLOCKS)
        call MPI_BCAST(SOL_pm_sample,1,mat5,0,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat5,ierr)

        !We broadcast RHS based on the sampled solutions of all ranks
        call mpimat4(mat4,npmsize,NN_coarse(1),NN_coarse(2),1)
        call MPI_BCAST(RHS_pm_coarse,1,mat4,0,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat4,ierr)


        !add Sampling blocks
       !allocate (SOL_pm_sumsample(NN_coarse(1),NN_coarse(2),NN_coarse(3),7))
       !SOL_pm_sumsample=0.d0
       !do nbj=1,BLOCKS
       !    SOL_pm_sumsample(:,:,:,:)=SOL_pm_sumsample(:,:,:,:)+SOL_pm_sample(:,:,:,:,nbj)
       !enddo
        !call calc_laplacian(SOL_pm_sumsample,RHS_pm_coarse,NN_coarse,NN_bl_coarse,Dpm_coarse,NN_bl_coarse)

        !if(my_rank.eq.0)then
            !outfil2='coarse_sumsample_grow'
            !call writegrow(RHS_pm_coarse,SOL_pm_sumsample,Dpm_coarse,outfil2,Xbound_coarse,NN_bl_coarse,NN_coarse)
        !endif

        ! if (my_rank.eq.1)starttime = MPI_WTIME()
        ! if(my_rank.eq.1)endtime = MPI_WTIME()
        ! if(my_rank.eq.1) write(*,*)'Laplace=',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        !if (my_rank.eq.1)starttime = MPI_WTIME()
        SOL_pm_coarse=0.d0
        iynbc=1
        itree= iyntree
        lmax=ilevmax
        if (itree.eq.0) lmax=1
        if (my_rank.eq.0) write(199,*) 'coarse'
        call pmesh(SOL_pm_coarse,RHS_pm_coarse,QP,XP,Xbound_coarse,DPm_coarse,NN_coarse,NN_bl_coarse,ND,&
            1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)
        ! if(my_rank.eq.1)endtime = MPI_WTIME()
        ! if(my_rank.eq.1) write(*,*)'Poisson Coarse=',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'

        call mapnodes_bl
       !if (my_rank.eq.0) then
       !     outfil2='coarse_laplacian'
       !     call writegrow(RHS_pm_coarse,SOL_pm_coarse,Dpm_coarse,outfil2,Xbound_coarse,NN_bl_coarse,NN_coarse)
       !endif
        ! Interpolate to all blocks
        
        ! At this point we calculate boundary conditions for each block and solve it 
        nb  = my_rank+1
        NXs = NNbl_bl(1,nb)
        NXf = NNbl_bl(4,nb)

        NYs = NNbl_bl(2,nb)
        NYf = NNbl_bl(5,nb)

        isizex = NXf - NXs + 1
        isizey = NYf - NYs + 1
        !---------------------------------------------------------------------------------
        !Xs boundary ,9 point stencil
        i=NXs
        node   = 0
        do j=NYs,NYf
            node=node + 1
            call interp_stencil
        enddo

        i=NXf
        node = isizey
        do j=NYs,NYf
            node = node + 1
            call interp_stencil
        enddo
        j=NYs
        node =2*isizey+1
        do i=NXs+1,NXf-1
            node = node + 1
            call interp_stencil
        enddo

        j=NYf
        node =2*isizey + isizex+1
        do i=NXs+1,NXf-1
            node = node + 1
            call interp_stencil
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
        ! write(*,*) 'Solving for Block',nb
        !iynbc=0 means that the bc's of the poisson solver are already defined
        iynbc=0
        itree=0
        lmax=0
        call pmesh(SOL_pm_bl,RHS_pm_bl,QP,XP,&
                   Xbound_tmp,Dpm_fine,NN_tmp,NN_bl_tmp,ND,1,ibctyp,neqs,neqf,iynbc,0,itree,lmax)

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        if(my_rank.eq.0)endtime = MPI_WTIME()
        if(my_rank.eq.0) write(199,*)'Parallel Poiss=',int((endtime-starttime)/60),'m',mod(endtime-starttime,60.d0),'s'
         !write(outfil1,'(a5,i2.2)') 'block',nb
         !call writesol_bl(RHS_pm_bl,SOL_pm_bl,Dpm_fine,outfil1,Xbound_tmp,NN_bl_tmp,NN_tmp)


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

            integer :: i_nb,j_nb,icnb,jcnb,inode1,jnode1,nbc,nb_plus,nb_minus,i_check,j_check,ik,jk
            integer :: nod,ibound,isizex,isizey,i_plus,j_plus,i_minus,j_minus
            double precision :: X1(3)
            double precision,allocatable :: Bound_sol(:,:)
            integer  :: status(MPI_STATUS_SIZE),source,ierr,mat2,dest

            NXs = NNbl_bl(1,nb)
            NXf = NNbl_bl(4,nb)
            NYs = NNbl_bl(2,nb)
            NYf = NNbl_bl(5,nb)
            nod=0
            nb = my_rank+1
            isizey = NYf-NYs + 1
            isizex = NXf-NXs + 1
            ibound = 2*(isizey + isizex)
            allocate(Bound_sol(ibound,neq))
            BBound=0.d0
            call mpimat2(mat2,ibound,neq)

            !SEND TO THE RIGHT BLOCKS
            !RIGHT SWEEP in I direction
            do k =1,4
                ik=0;jk=0;Bound_sol=0
                if (k.eq.1) then
                    i=NXf
                    nod=0!it's receiving from 1 to isizey
                    do j=NYs,NYf
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo                                               
                                                                        
                    j=NYs                                               
                    !because of corner nodes(isizey+1)                  
                    nod = 2*isizey + 1                                  
                    do i=NXf+1,NNbl(1,nb)                               
                        nod = nod + 1                                   
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo                                               
                                                                        
                    j=NYf                                               
                    nod = 2*isizey + isizex + 1                         
                    do i=NXf+1,NNbl(1,nb)                               
                        nod = nod + 1                                   
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo                                               
                    ik=1                                                
                endif                                                   
                if (k.eq.2)  then                                       
                    j= NYf                                              
                    !it's receiving from at NYs so 2*isizey             
                    nod=2*isizey+1                                      
                    do i=NXs+1,NXf-1                                    
                        nod = nod + 1                                   
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo                                               
                                                                        
                    i= NXs                                              
                    nod =0                                              
                    do j=NYf,NNbl(2,nb)                                 
                        nod = nod + 1                                   
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo                                               
                                                                        
                    i= NXf                                              
                    nod =isizey                                         
                    do j=NYf,NNbl(2,nb)                                 
                        nod = nod + 1                                   
                        Bound_sol(nod,neqs:neqf ) = SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    jk=1;ik=0
                endif
                if (k.eq.3) then
                    !DIAGONAL(/)
                    i=NXf
                    nod=0!the first part of the diagonal is at NXs so zero
                    do j=NYf,NNbl(2,nb)
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    nod = 2*isizey+1!the second diagonal part is at NYs
                    j=NYf
                    do i=NXf+1,NNbl(1,nb)
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    ik=1;jk=1
                endif
                if (k.eq.4) then
                    !DIAGONAL(\)
                    i=NXf
                    nod=isizey+1!It is at Nxs but from isize and down(WATCH nod=nod-1)
                    do j=NYs,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    j=NYs
                    nod = 2*isizey + isizex+1 !Because it is at NYf b
                    do i=NXf+1,NNbl(1,nb)
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    ik=1;jk=-1
                endif

                i_plus = nb_i + ik
                j_plus = nb_j-1 + jk +  1!because in equation (j-1)
                nb_plus = (j_plus-1)*NBI + i_plus

                if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1) then
                    dest = nb_plus-1 ! because nb = rank-1
                    !write(*,*) nb-1,'sending to ',dest
                    call MPI_SEND(Bound_sol,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
                    i_minus = nb_i - ik
                    j_minus = nb_j - 1 - jk  + 1
                    nb_minus = (j_minus - 1) * NBI + i_minus
                    if (i_minus.ge.1.and.j_minus.ge.1.and.j_minus.le.NBJ.and.i_minus.le.NBI)  then
                        source=nb_minus-1
                        !  write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_minus) = BBound(:,:,nb_minus) + Bound_sol(:,:)
                    endif
                else
                    i_minus = nb_i - ik
                    j_minus = nb_j - 1 - jk  + 1
                    nb_minus = (j_minus - 1) * NBI + i_minus

                    if (i_minus.ge.1.and.j_minus.ge.1.and.j_minus.le.NBJ.and.i_minus.le.NBI)  then
                        source=nb_minus-1
                        !  write(*,*) nb-1,' received from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_minus) = BBound(:,:,nb_minus) + Bound_sol(:,:)
                    endif
                endif
            enddo
            !LEFT SWEEP in I


            do k =1,4
                ik=0;jk=0;Bound_sol=0
                if (k.eq.1) then
                    i= NXs
                    nod=isizey!because it is at the start of NXf
                    do j=NYs,NYf
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    j=NYs
                    nod =2*isizey+isizex
                    do i=NXs-1,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    j=NYf
                    nod =2*isizey+2*isizex
                    do i=NXs-1,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    ik=1;jk=0;
                endif
                if (k.eq.2)  then
                    j= NYs
                    nod=2*isizey + isizex+1!!because it at the start of NYf
                    do i=NXs+1,NXf-1
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    i=NXs
                    nod =isizey+1
                    do j=NYs,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo

                    i=NXf
                    nod =2*isizey+1
                    do j=NYs,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    jk=1;ik=0

                endif
                if (k.eq.3) then
                    !DIAGONAL(/)
                    i=NXs
                    nod=2*isizey+1!first part at NXf but for 2*isizey
                    do j=NYs,1,-1
                        nod = nod -1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    !NXs,NYs from above
                    j=NYs
                    nod = 2*isizey+2*isizex!same
                    do i=NXs-1,1,-1!-1 because it counst differentely
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    ik=1;jk=1
                endif
                if (k.eq.4) then
                    !DIAGONAL(\)
                    i=NXs
                    nod = isizey
                    do j=NYf,NNbl(2,nb)
                        nod = nod + 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    j=NYf
                    nod = 2*isizey + isizex
                    !NXs,NYs from above
                    do i=NXs-1,1,-1
                        nod = nod - 1
                        Bound_sol(nod,neqs:neqf ) =  SOL_pm_bl(neqs:neqf,i,j,1)
                    enddo
                    ik=1;jk=-1
                endif

                !fix Nodes at corners we want to put all contributions(from y nodes,corners node exist twice)
                ! in (1,isizey,isizey+1,2*isizey)

                i_minus = nb_i - ik
                j_minus = nb_j - 1 - jk  + 1
                nb_minus = (j_minus - 1) * NBI + i_minus
                if (i_minus.ge.1.and.j_minus.ge.1.and.i_minus.le.NBI.and.j_minus.le.NBJ) then
                    dest = nb_minus-1 ! because nb = rank-1
                    !write(*,*) nb-1,'sending to ',dest
                    call MPI_SEND(Bound_sol,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
                    i_plus = nb_i + ik
                    j_plus = nb_j-1 + jk +  1!because in equation (j-1)
                    nb_plus = (j_plus-1)*NBI + i_plus
                    if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1)  then
                        source=nb_plus-1
                        !    write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_plus) = BBound(:,:,nb_plus) + Bound_sol(:,:)
                    endif
                else
                    i_plus = nb_i + ik
                    j_plus = nb_j-1 + jk +  1!because in equation (j-1)
                    nb_plus = (j_plus-1)*NBI + i_plus
                    if (i_plus.le.NBI.and.j_plus.le.NBJ.and.i_plus.ge.1.and.j_plus.ge.1)  then
                        source=nb_plus-1
                        !    write(*,*) nb-1,'receiving  from',source
                        call MPI_RECV(Bound_sol,1,mat2,source,1,MPI_COMM_WORLD,status,ierr)
                        BBound(:,:,nb_plus) = BBound(:,:,nb_plus) + Bound_sol(:,:)
                    endif
                endif

            enddo

            call MPI_TYPE_FREE(mat2,ierr)
            ! if (my_rank.eq.0) then
            ! do nbj=1,BLOCKS
            !     write(*,*) nbj,BBound(isizey,1,nbj),BBound(2*isizey+isizex+1,1,nbj)
            ! enddo
            ! endif
            ! BBound(1,1:neq,1:BLOCKS)        = BBound(1,1:neq,1:BLOCKS) + &
                !                                   BBound(2*isizey+1,1:neq,1:BLOCKS)

            ! BBound(isizey,1:neq,1:BLOCKS)   = BBound(isizey,1:neq,1:BLOCKS) + &
                !                                   BBound(2*isizey+isizex+1,1:neq,1:BLOCKS)

            ! BBound(isizey+1,1:neq,1:BLOCKS) = BBound(isizey+1,1:neq,1:BLOCKS) + &
                !                                   BBound(2*isizey + isizex,1:neq,1:BLOCKS)

            ! BBound(2*isizey,1:neq,1:BLOCKS) = BBound(2*isizey,1:neq,1:BLOCKS) + &
                !                                   BBound(2*isizey + 2*isizex,1:neq,1:BLOCKS)

        End Subroutine mapnodes_bl

        Subroutine interp_stencil
            double precision :: projection_fun,addlocal,X1(2),add(neq),add_sample(neq)
            integer          :: inode1,jnode1,i_nb,j_nb

            iproj=3
            X(1) = Xbound_bl(1,nb) + (i-1)*Dpm_fine(1)
            X(2) = Xbound_bl(2,nb) + (j-1)*Dpm_fine(2)

            !   SOL_pm_bl(i,j,1,1,nb)=SOL_pm_fine(inode,jnode,1,1)
            !inode =int(((X(1) - Xbound_coarse(1)) / Dpm_coarse(1))) + 1
            !jnode =int(((X(2) - Xbound_coarse(2)) / Dpm_coarse(2))) + 1
            inode =int(nint((X(1) - Xbound_coarse(1)) / Dpm_coarse(1))) + 1
            jnode =int(nint((X(2) - Xbound_coarse(2)) / Dpm_coarse(2))) + 1
            !--We search the 4 nodes close to the particles
            k=1
            nnb=0
            addlocal=0
            do nbc = 1,BLOCKS
                if(map_nodes(inode-1,jnode-1,1,nbc).ne.1.or.&
                        map_nodes(inode+1,jnode+1,1,nbc).ne.1)then
                    nnb(nbc)=1
                endif
            enddo
            do jc = jnode-1 , jnode + 1
                do ic = inode-1 , inode + 1
                    add=0.d0
                    add_sample=0.d0
                    do nbc=1,BLOCKS
                        add_sample(neqs:neqf) = SOL_pm_sample(neqs:neqf,ic,jc,1,nbc)+add_sample(neqs:neqf)
                        if(nnb(nbc).eq.1) then
                            !add fine local contribution when available
                            add(neqs:neqf)=SOL_pm_sample(neqs:neqf,ic,jc,1,nbc) + add(neqs:neqf)
                        endif
                    enddo
                    xc  = (X(1)- Xbound_coarse(1) - (ic-1) * Dpm_coarse(1)) / Dpm_coarse(1)
                    fx = projection_fun(iproj,xc)

                    yc  = (X(2)- Xbound_coarse(2) - (jc-1) * Dpm_coarse(2)) / Dpm_coarse(2)
                    fy = projection_fun(iproj,yc)

                    f = fx * fy
                    !SOL_pm_coa
                    SOL_pm_bl(neqs:neqf,i,j,1) = SOL_pm_bl(neqs:neqf,i,j,1)+ &
                        f * (SOL_pm_coarse(neqs:neqf,ic,jc,1)-add_sample(neqs:neqf)+&
                        add(neqs:neqf))
                enddo
            enddo

            !add local data in the blocks that sampled data was interpolated
            do j_nb = nb_j-1, nb_j+1
                do i_nb = nb_i-1, nb_i+1
                    if (i_nb.lt.1.or.i_nb.gt.NBI.or.j_nb.lt.1.or.j_nb.gt.NBJ)   cycle
                    nbc = (j_nb-1)*NBI + i_nb
                    if (nnb(nbc).eq.1.or.nbc.eq.nb) cycle
                    SOL_pm_bl(neqs:neqf,i,j,1)=SOL_pm_bl(neqs:neqf,i,j,1) + BBound(node,neqs:neqf ,nbc)
                enddo
            enddo


        End Subroutine interp_stencil


    End Subroutine yaps2d

    Subroutine calc_laplacian_coarse(SOL_pm,RHS_pm,NN,NN_bl,Dpm,NN_map,neqs,neqf)
        Implicit none

        integer,intent(in) :: NN(3),NN_bl(6),NN_map(6),neqs,neqf
        double precision,intent(out)   :: RHS_pm(neqf,NN(1),NN(2),NN(3))
        double precision,intent(in)    :: SOL_pm(neqf,NN(1),NN(2),NN(3)),Dpm(3)
        integer                        :: i,j



        RHS_pm=0.d0
        do i=NN_map(1)+1,NN_map(4)-1
            do j=NN_map(2)+1,NN_map(5)-1
                RHS_pm(neqs:neqf,i,j,1)=(SOL_pm(neqs:neqf,i+1,j,1) - 2*SOL_pm(neqs:neqf,i,j,1) +&
                                         SOL_pm(neqs:neqf,i-1,j,1))/Dpm(1)**2 + &
                                        (SOL_pm(neqs:neqf,i,j+1,1) - 2*SOL_pm(neqs:neqf,i,j,1) +&
                                         SOL_pm(neqs:neqf,i,j-1,1))/Dpm(2)**2



            enddo
        enddo

    End Subroutine calc_laplacian_coarse
