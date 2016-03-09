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

Subroutine mpimat4int(mat4int,nsize1,nsize2,nsize3,nsize4)
    use MPI

    integer, intent(in)::nsize1,nsize2,nsize3,nsize4
    integer ierr
    integer ::imat(4),mat(4),start(4)
    integer ::mat4int
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

    call MPI_TYPE_CREATE_SUBARRAY(4,imat,mat,start,MPI_ORDER_FORTRAN,MPI_INTEGER,mat4int,ierr)
    call MPI_TYPE_COMMIT(mat4int,ierr)

End subroutine mpimat4int

Subroutine mpimat3_pm(mat3,nsize1,nsize2,nsize3)
    use MPI

    integer, intent(in)::nsize1,nsize2,nsize3
    integer ierr
    integer ::imat(3),mat(3),start(3)
    integer ::mat3
    !allocate(struct%AS_ij(nsize,nsize))
    imat(1)=nsize1
    imat(2)=nsize2
    imat(3)=nsize3
    mat(1)=nsize1
    mat(2)=nsize2
    mat(3)=nsize3
    start(1)=0
    start(2)=0
    start(3)=0

    call MPI_TYPE_CREATE_SUBARRAY(3,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat3,ierr)
    call MPI_TYPE_COMMIT(mat3,ierr)

End subroutine mpimat3_pm


Subroutine rhsbcast(RHS_pm,NN,neq)
    use MPI

    integer,intent(in)          ::NN(3),neq
    double precision,intent(inout) ::RHS_pm(neq,NN(1),NN(2),NN(3))
    integer :: my_rank,np,ierr,mat4


    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

    !---------------------------------------------
    call mpimat4(mat4,neq,NN(1),NN(2),NN(3))
    call MPI_BCAST(RHS_pm,1,mat4,0,MPI_COMM_WORLD,ierr)
    call MPI_TYPE_FREE(mat4,ierr)
    !-----------------------------------
End Subroutine rhsbcast

Subroutine rhsscat(BLOCKS,NN_tmp,NNbl,NNbl_bl,NN_bl,nb_i,nb_j,RHS_pm_bl)
    use vpm_vars        
    use pmgrid
    use pmeshpar
    use MPI
    Implicit None
    integer,intent(in) ::BLOCKS,NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS),nb_i,nb_j,NN_bl(6),NN_tmp(3)
    double precision,intent(out) ::RHS_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
    integer :: my_rank,ierr,source,dest,status(MPI_STATUS_SIZE),mat4
    integer :: ixs,jxs,ixf,jxf,nb,NXs,NXf,NYs,NYf,j,NN(3),nbs

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)

    RHS_pm_bl=0
    nb  = my_rank + 1
    NN(1:3)=NNbl(1:3,nb)
    NXs = NNbl_bl(1,nb)
    NXf = NNbl_bl(4,nb)

    NYs = NNbl_bl(2,nb)
    NYf = NNbl_bl(5,nb)
    ixs=(nb_i-1)*(NXf-NXs)+NN_bl(1)
    jxs=(nb_j-1)*(NYf-NYs)+NN_bl(2)

    ixf= ixs + (NXf-NXs+1)-1
    jxf= jxs + (NYf-NYs+1)-1

    RHS_pm_bl(1:neqpm,NXs:NXf,NYs:NYf,1) = RHS_pm(1:neqpm,ixs:ixf,jxs:jxf,1)
    if (nb_i.gt.1) RHS_pm_bl(:,NXs,:,:)=0.d0
    if (nb_j.gt.1) RHS_pm_bl(:,:,NYs,:)=0.d0


End Subroutine rhsscat

Subroutine solget(BLOCKS,NBI,NBJ,NN_tmp,NNbl,NNbl_bl,NN_bl,SOL_pm_bl)
    use vpm_vars        
    use pmgrid
    use pmeshpar
    use MPI
    Implicit None
    integer,intent(in) ::BLOCKS,NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS),NBI,NBJ,NN_bl(6),NN_tmp(3)
    double precision,intent(in)  :: SOL_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
    double precision,allocatable :: SOL_pm_tmp(:,:,:,:)
    integer :: my_rank,ierr,source,dest,status(MPI_STATUS_SIZE),mat4
    integer :: ixs,jxs,ixf,jxf,nb,NXs,NXf,NYs,NYf,j,NN_block(3),i,nbs

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)


    nb  = my_rank + 1
    !NN(1:3)=NN_tmp(1:3)
    !do nbs=1,BLOCKS
    !  write(*,*) my_rank,NN_tmp
    !enddo
    if (my_rank.eq.0) then
        j=1
        i=1
        NXs = NNbl_bl(1,1)
        NXf = NNbl_bl(4,1)

        NYs = NNbl_bl(2,1)
        NYf = NNbl_bl(5,1)
        ixs=(i-1)*(NXf-NXs)+NN_bl(1)
        jxs=(j-1)*(NYf-NYs)+NN_bl(2)

        ixf= ixs + (NXf-NXs+1)-1
        jxf= jxs + (NYf-NYs+1)-1

        SOL_pm(1:neqpm,ixs:ixf,jxs:jxf,1)= SOL_pm_bl(1:neqpm,NXs:NXf,NYs:NYf,1)
        !-->Assign
        do j=1,NBJ
            do i=1,NBI
                nbs = (j-1)*NBI + i
                if (nbs.eq.1)cycle
                NN_block(1:3)=NNbl(1:3,nbs)
                allocate(SOL_pm_tmp(neqpm,NN_block(1),NN_block(2),NN_block(3)))
                NXs = NNbl_bl(1,nbs)
                NXf = NNbl_bl(4,nbs)

                NYs = NNbl_bl(2,nbs)
                NYf = NNbl_bl(5,nbs)
                ixs=(i-1)*(NXf-NXs)+NN_bl(1)
                jxs=(j-1)*(NYf-NYs)+NN_bl(2)

                ixf= ixs + (NXf-NXs+1)-1
                jxf= jxs + (NYf-NYs+1)-1
                call mpimat4(mat4,neqpm,NN_block(1),NN_block(2),NN_block(3))
                source = nbs-1
                call MPI_RECV(SOL_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                SOL_pm(1:neqpm,ixs:ixf,jxs:jxf,1)= SOL_pm_tmp(1:neqpm,NXs:NXf,NYs:NYf,1)
                call MPI_TYPE_FREE(mat4,ierr)
                deallocate(SOL_pm_tmp)
            enddo
        enddo
    else
        dest=0
        call mpimat4(mat4,neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
        call MPI_SEND(SOL_pm_bl,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat4,ierr)
    endif

End Subroutine solget

Subroutine rhsscat_3d(BLOCKS,NN_tmp,NNbl,NNbl_bl,NN_bl,nb_i,nb_j,nb_k,RHS_pm_bl)
    use vpm_vars          
    use pmgrid
    use pmeshpar
    use MPI
    Implicit None
    integer,intent(in) ::BLOCKS,NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS),nb_i,nb_j,nb_k,NN_bl(6),NN_tmp(3)
    double precision,intent(out) ::RHS_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
    integer :: my_rank,ierr,source,dest,status(MPI_STATUS_SIZE),mat4
    integer :: ixs,jxs,kxs,ixf,jxf,kxf,nb,NXs,NXf,NYs,NYf,NZs,NZf,j,NN(3),nbs

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)

    RHS_pm_bl=0
    nb  = my_rank + 1
    NN(1:3)=NNbl(1:3,nb)

    NXs = NNbl_bl(1,nb)
    NXf = NNbl_bl(4,nb)

    NYs = NNbl_bl(2,nb)
    NYf = NNbl_bl(5,nb)

    NZs = NNbl_bl(3,nb)
    NZf = NNbl_bl(6,nb)
    ixs=(nb_i-1)*(NXf-NXs) +NN_bl(1)
    jxs=(nb_j-1)*(NYf-NYs) +NN_bl(2)
    kxs=(nb_k-1)*(NZf-NZs) +NN_bl(3)

    ixf= ixs + (NXf-NXs+1)-1
    jxf= jxs + (NYf-NYs+1)-1
    kxf= kxs + (NZf-NZs+1)-1
 
    RHS_pm_bl(1:neqpm,NXs:NXf,NYs:NYf,NZs:NZf) = RHS_pm(1:neqpm,ixs:ixf,jxs:jxf,kxs:kxf)
    if (nb_i.gt.1) RHS_pm_bl(:,NXs,:,:)=0.d0
    if (nb_j.gt.1) RHS_pm_bl(:,:,NYs,:)=0.d0
    if (nb_k.gt.1) RHS_pm_bl(:,:,:,NZs)=0.d0


End Subroutine rhsscat_3d

Subroutine solget_3d(BLOCKS,NBI,NBJ,NBK,NN_tmp,NNbl,NNbl_bl,NN_bl,SOL_pm_bl)
    use vpm_vars         
    use pmgrid
    use pmeshpar
    use MPI
    Implicit None
    integer,intent(in)           ::BLOCKS,NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS),NBI,NBJ,NBK,NN_bl(6),NN_tmp(3)
    double precision,intent(in)  :: SOL_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
    double precision,allocatable :: SOL_pm_tmp(:,:,:,:)
    integer :: my_rank,ierr,source,dest,status(MPI_STATUS_SIZE),mat4
    integer :: ixs,jxs,kxs,ixf,jxf,kxf,nb,NXs,NXf,NYs,NYf,NZs,NZf,j,k,NN_block(3),i,nbs

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)


    nb  = my_rank + 1
    !NN(1:3)=NN_tmp(1:3)
    !do nbs=1,BLOCKS
    !  write(*,*) my_rank,NN_tmp
    !enddo
    if (my_rank.eq.0) then
        j=1
        i=1
        k=1
        NXs = NNbl_bl(1,nb)
        NXf = NNbl_bl(4,nb)

        NYs = NNbl_bl(2,nb)
        NYf = NNbl_bl(5,nb)

        NZs = NNbl_bl(3,nb)
        NZf = NNbl_bl(6,nb)

        ixs=(i-1)*(NXf-NXs)+NN_bl(1)
        jxs=(j-1)*(NYf-NYs)+NN_bl(2)
        kxs=(k-1)*(NZf-NZs)+NN_bl(3)

        ixf= ixs + (NXf-NXs+1)-1
        jxf= jxs + (NYf-NYs+1)-1
        kxf= kxs + (NZf-NZs+1)-1

        SOL_pm(1:neqpm,ixs:ixf,jxs:jxf,kxs:kxf)= SOL_pm_bl(1:neqpm,NXs:NXf,NYs:NYf,NZs:NZf)
       !write(*,*) maxval(SOL_pm_bl),minval(SOL_pm_bl)
        !-->Assign
        do k=1,NBK
           do j=1,NBJ
               do i=1,NBI
                   nbs = (k-1)*NBI*NBJ+(j-1)*NBI + i
                   if (nbs.eq.1)cycle
                   NN_block(1:3)=NNbl(1:3,nbs)
                   allocate(SOL_pm_tmp(neqpm,NN_block(1),NN_block(2),NN_block(3)))
                   NXs = NNbl_bl(1,nbs)
                   NXf = NNbl_bl(4,nbs)
           
                   NYs = NNbl_bl(2,nbs)
                   NYf = NNbl_bl(5,nbs)

                   NZs = NNbl_bl(3,nbs)
                   NZf = NNbl_bl(6,nbs)

                   ixs=(i-1)*(NXf-NXs)+NN_bl(1)
                   jxs=(j-1)*(NYf-NYs)+NN_bl(2)
                   kxs=(k-1)*(NZf-NZs)+NN_bl(3)

                   ixs=(i-1)*(NXf-NXs)+NN_bl(1)
                   jxs=(j-1)*(NYf-NYs)+NN_bl(2)
                   kxs=(k-1)*(NZf-NZs)+NN_bl(3)
           
                   ixf= ixs + (NXf-NXs+1)-1
                   jxf= jxs + (NYf-NYs+1)-1
                   kxf= kxs + (NZf-NZs+1)-1
                   
                   call mpimat4(mat4,neqpm,NN_block(1),NN_block(2),NN_block(3))
                   source = nbs-1
                   call MPI_RECV(SOL_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
                   SOL_pm(1:neqpm,ixs:ixf,jxs:jxf,kxs:kxf)= SOL_pm_tmp(1:neqpm,NXs:NXf,NYs:NYf,NZs:NZf)
                   call MPI_TYPE_FREE(mat4,ierr)
                   deallocate(SOL_pm_tmp)
               enddo
           enddo
        enddo
    else
        dest=0
        call mpimat4(mat4,neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3))
        call MPI_SEND(SOL_pm_bl,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat4,ierr)
    endif

End Subroutine solget_3d


Subroutine velbcast_3d
    use vpm_vars        
    use pmgrid
    use pmeshpar
    use MPI
    integer :: my_rank,np,ierr,mat3

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

    !---------------------------------------------
    call mpimat3_pm(mat3,NXpm,NYpm,NZpm)
    call MPI_BCAST(velvrx_pm,1,mat3,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(velvry_pm,1,mat3,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(velvrz_pm,1,mat3,0,MPI_COMM_WORLD,ierr)
    call MPI_TYPE_FREE(mat3,ierr)
    !--------------------------------------------
End Subroutine velbcast_3d

Subroutine particles_scat
    use vpm_vars        
    use pmgrid
    use pmeshpar
    use parvar
    use MPI
    
    Implicit None
    integer :: my_rank,np,ierr,i
    integer :: dest,NVR_pr,NVR_r,source,mat2
    integer :: status(MPI_STATUS_SIZE)
    character *15 :: filname,filname2
    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

    !---------------------------------------------
    if (my_rank.eq.0) then 
        XP_scatt(1:3,1:NVR_p)= XP(1:3,1:NVR_p)
        QP_scatt(1:neqpm+1,1:NVR_p)= QP(1:neqpm+1,1:NVR_p)
        NVR_pr = NVR_p
        NVR_r  = NVR/np
        if (NVR_r.gt.0)  then 
            do i=2,np
               dest = i-1
               call mpimat2_pm(mat2,3,NVR_size,3,NVR_r,NVR_pr)!NVR_pr because counting starts from 0
               call MPI_SEND(XP,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
               call MPI_TYPE_FREE(mat2,ierr)
           
               call mpimat2_pm(mat2,neqpm+1,NVR_size,neqpm+1,NVR_r,NVR_pr)
               call MPI_SEND(QP,1,mat2,dest,1,MPI_COMM_WORLD,ierr)
               call MPI_TYPE_FREE(mat2,ierr)
               NVR_pr=NVR_pr+NVR_r
            enddo
        endif
    else
       if (NVR_p.gt.0) then 
          call mpimat2_pm(mat2,3,NVR_p,3,NVR_p,0)
          call MPI_RECV(XP_scatt,1,mat2,0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_TYPE_FREE(mat2,ierr)
          
          call mpimat2_pm(mat2,neqpm+1,NVR_p,neqpm+1,NVR_p,0)
          call MPI_RECV(QP_scatt,1,mat2,0,1,MPI_COMM_WORLD,status,ierr)
          call MPI_TYPE_FREE(mat2,ierr)
       endif
    endif

   !write(filname,'(i1)') my_rank+1
   !open(15,file=filname)
   !write(15,*) 'VARIABLES="X" "Y" "Z"'
   !do i = 1,NVR_p
   !   write(15,'(7(e28.17,1x))') XP_scatt(1:3,i)!,QP_scatt(1:neqpm+1,i)
   !enddo
End Subroutine particles_scat

Subroutine particles_gath
    use vpm_vars        
    use pmgrid
    use pmeshpar
    use parvar
    use MPI
    
    Implicit None
    integer :: my_rank,np,ierr,i
    integer :: dest,NVR_pr,NVR_r,source,mat2
    integer :: status(MPI_STATUS_SIZE)
    character *15 :: filname,filname2
    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

    !---------------------------------------------
    if (my_rank.eq.0) then 
        XP(1:3,1:NVR_p)= XP_scatt(1:3,1:NVR_p)
        QP(1:neqpm+1,1:NVR_p)= QP_scatt(1:neqpm+1,1:NVR_p)
        UP(1:3,1:NVR_p)= UP_scatt(1:3,1:NVR_p)
        GP(1:3,1:NVR_p)= GP_scatt(1:3,1:NVR_p)
        deallocate(QP_scatt,XP_scatt,UP_scatt,GP_scatt)
        NVR_pr = NVR_p
        NVR_r  = NVR/np
        allocate(XP_scatt(3,NVR_r),QP_scatt(neqpm+1,NVR_r),UP_scatt(3,NVR_r),GP_scatt(3,NVR_r))
        if (NVR_r.gt.0)  then 
            do i=2,np
               dest = i-1
               call mpimat2_pm(mat2,3,NVR_r,3,NVR_r,0)!NVR_pr because counting starts from 0
               call MPI_RECV(XP_scatt,1,mat2,dest,1,MPI_COMM_WORLD,status,ierr)
               call MPI_TYPE_FREE(mat2,ierr)
           
               call mpimat2_pm(mat2,neqpm+1,NVR_r,neqpm+1,NVR_r,0)
               call MPI_RECV(QP_scatt,1,mat2,dest,1,MPI_COMM_WORLD,status,ierr)
               call MPI_TYPE_FREE(mat2,ierr)

               call mpimat2_pm(mat2,3,NVR_r,3,NVR_r,0)!NVR_pr because counting starts from 0
               call MPI_RECV(UP_scatt,1,mat2,dest,1,MPI_COMM_WORLD,status,ierr)
               call MPI_TYPE_FREE(mat2,ierr)

               call mpimat2_pm(mat2,3,NVR_r,3,NVR_r,0)!NVR_pr because counting starts from 0
               call MPI_RECV(GP_scatt,1,mat2,dest,1,MPI_COMM_WORLD,status,ierr)
               call MPI_TYPE_FREE(mat2,ierr)

               XP(1:3,NVR_pr+1:NVR_pr+NVR_r)= XP_scatt(1:3,1:NVR_r)
               QP(1:neqpm+1,NVR_pr+1:NVR_pr+NVR_r)= QP_scatt(1:neqpm+1,1:NVR_r)
               UP(1:3,NVR_pr+1:NVR_pr+NVR_r)= UP_scatt(1:3,1:NVR_r)
               GP(1:3,NVR_pr+1:NVR_pr+NVR_r)= GP_scatt(1:3,1:NVR_r)

               NVR_pr=NVR_pr+NVR_r
            enddo
        endif
    else
       if (NVR_p.gt.0) then 
          call mpimat2_pm(mat2,3,NVR_p,3,NVR_p,0)
          call MPI_SEND(XP_scatt,1,mat2,0,1,MPI_COMM_WORLD,ierr)
          call MPI_TYPE_FREE(mat2,ierr)
          
          call mpimat2_pm(mat2,neqpm+1,NVR_p,neqpm+1,NVR_p,0)
          call MPI_SEND(QP_scatt,1,mat2,0,1,MPI_COMM_WORLD,ierr)
          call MPI_TYPE_FREE(mat2,ierr)

          call mpimat2_pm(mat2,3,NVR_p,3,NVR_p,0)
          call MPI_SEND(UP_scatt,1,mat2,0,1,MPI_COMM_WORLD,ierr)
          call MPI_TYPE_FREE(mat2,ierr)
          
          call mpimat2_pm(mat2,3,NVR_p,3,NVR_p,0)
          call MPI_SEND(GP_scatt,1,mat2,0,1,MPI_COMM_WORLD,ierr)
          call MPI_TYPE_FREE(mat2,ierr)
       endif
    endif

   !write(filname,'(i1)') my_rank+1
   !open(15,file=filname)
   !write(15,*) 'VARIABLES="X" "Y" "Z"'
   !do i = 1,NVR_p
   !   write(15,'(7(e28.17,1x))') XP_scatt(1:3,i)!,QP_scatt(1:neqpm+1,i)
   !enddo
End Subroutine particles_gath

Subroutine proj_gath(NN)
    use vpm_vars       
    use pmgrid
    use pmeshpar
    use parvar
    use MPI
    Implicit None
    integer,intent(in) :: NN(3)
    integer :: my_rank,np,ierr,i
    integer :: dest,NVR_pr,NVR_r,source,mat4
    integer :: status(MPI_STATUS_SIZE)
    double precision,allocatable:: RHS_pm_tmp(:,:,:,:)
    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)

    

    if (my_rank.eq.0) then 
      allocate (RHS_pm_tmp(neqpm+1,NN(1),NN(2),NN(3)))
      do i=2,np
         source =i-1
         call mpimat4(mat4,neqpm+1,NN(1),NN(2),NN(3))
         call MPI_RECV(RHS_pm_tmp,1,mat4,source,1,MPI_COMM_WORLD,status,ierr)
         call MPI_TYPE_FREE(mat4,ierr)
         RHS_pm=RHS_pm + RHS_pm_tmp
      enddo
      deallocate(RHS_pm_tmp)
    else
        dest=0
        call mpimat4(mat4,neqpm+1,NN(1),NN(2),NN(3))
        call MPI_SEND(RHS_pm,1,mat4,dest,1,MPI_COMM_WORLD,ierr)
        call MPI_TYPE_FREE(mat4,ierr)
    endif
End Subroutine 

Subroutine mpimat2_pm(mat2,orig1,orig2,nsize1,nsize2,istart)
use MPI

Implicit None
integer ierr
integer::typelist(2)
integer ::imat(2),mat(2),start(2)
integer ::istart
integer ::orig1,orig2,nsize1,nsize2,mat2
!allocate(struct%AS_ij(nsize,nsize))
imat(1)= orig1
imat(2)= orig2

mat(1) = nsize1
mat(2) = nsize2

start(1) = 0
start(2) = istart

typelist(1)=MPI_DOUBLE_PRECISION
typelist(2)=MPI_DOUBLE_PRECISION

!write (*,*) nsize
call MPI_TYPE_CREATE_SUBARRAY(2,imat,mat,start,MPI_ORDER_FORTRAN,MPI_DOUBLE_PRECISION,mat2,ierr)
call MPI_TYPE_COMMIT(mat2,ierr)

End Subroutine mpimat2_pm

