 Module vpm_vars
    double precision, allocatable :: XP_scatt(:,:),QP_scatt(:,:),UP_scatt(:,:),GP_scatt(:,:)
    double precision              :: DT_c,V_ref, NI
    integer ,allocatable          :: NVR_projscatt(:)
    integer                       :: interf_iproj,ncell_rem,iynhj

    integer                       :: ncoarse,nparcell1d
    integer                       :: neqpm,NVR_p,NVR_size,iwrite,NTIME_pm


    integer,save                       :: IPMWRITE,mrem,idefine,iynslice
    integer,save                       :: IPMWSTART(10),IPMWSTEPS(10)   
 End Module vpm_vars
 
 Module vpm_size
    double precision,save :: Xbound(6),Dpm(3),Xbound0(6),Dpm0(3)
    integer,save          :: NN_bl(6),NN(3),NN0_bl(6),NN0(3)
    integer,save          :: NXs0_bl(10),NYs0_bl(10),NXf0_bl(10),NYf0_bl(10),NZs0_bl(10),NZf0_bl(10)

    double precision,allocatable,save  :: Xbound_bl(:,:)
    integer         ,allocatable,save  :: NNbl_bl(:,:),NNbl(:,:)


    double precision,save              :: Xbound_tmp(6),Xbound_coarse(6),Dpm_coarse(3)
    integer,save                       :: NN_tmp(3),NN_bl_tmp(6),NN_coarse(3),NN_bl_coarse(6)
    integer,save                       :: nb_i,nb_j,nb_k,NBB,NXbl,NYbl,NZbl,BLOCKS,NXB,NYB,NZB,ndumcell_coarse ,ndumcell_bl
    double precision              :: starttime,endtime,st,et,ct
    integer,save                       :: iynbc,iret,NBI,NBJ,NBK,NVR_out_thres,NREMESH,ntorder,&
                                     iyntree,ilevmax,itree,nsize_out(3),ibctyp,NWRITE

    
 End Module vpm_size

Module test_app
    double precision,allocatable    :: analytic_sol(:,:,:,:)
End module test_app

Module openmpth
   integer                       ::OMPTHREADS
End Module openmpth

Module vpm_lib
  use vpm_vars
  use vpm_size
  use openmpth

 !private ::  starttime,endtime,st,et,ct
 !private ::  nb_i,nb_j,nb_k,NBB,NXbl,NYbl,NZbl,BLOCKS,NXB,NYB,NZB,ndumcell_coarse ,ndumcell_bl
 !private :: II,iynbc,iret,NBI,NBJ,NBK,NVR_out_thres,NREMESH,ntorder,&
 !                                   iyntree,ilevmax,itree,nsize_out,ibctyp
 
contains

 Subroutine vpm(XP_in,QP_in,UP_in,GP_in,NVR_in,neqpm_in,WhatToDo,&
                RHS_pm_in,Velx,Vely,Velz,NTIME_in,NI_in,NVRM_in)
    use vpm_vars
    use vpm_size
    use pmeshpar
    use parvar
    use pmgrid
    use MPI
    use pmlib
    use projlib
    use yapslib
    use openmpth

    Implicit None 
    
    double precision,intent(inout),target  :: XP_in(:,:),QP_in(:,:),UP_in(:,:),GP_in(:,:)
    double precision,intent(inout),pointer :: RHS_pm_in(:,:,:,:)
    double precision,intent(inout),pointer :: velx(:,:,:),vely(:,:,:),velz(:,:,:)
    integer,         intent(inout)         :: NVR_in
    integer,         intent(in)            :: neqpm_in,WhatToDo,NVRM_in,NTIME_in
    double precision, intent(in)           :: NI_in
    double precision :: a,XPM,YPM,totmass,totvor,MACH,error,pr
    integer          :: i , j, k , nb, NXs,&
                         NYs, NXf, NYf, NZs,NZf,ivel,ij,iter,IPM,ialloc
    integer,allocatable :: ieq(:)

    integer          :: ijk,npar,itype,nel,iflag,inf,iface,rk_step
    integer          :: ierr,my_rank,mat4,mat3,np,ICONV,omp_get_num_threads,omp_get_max_threads
    double precision :: xi,yi,ksi1,ksi2,th1,th2,w1,w2,xl,yl,XO(3),r
    double precision              :: errorCFD,errorf
    double precision ,allocatable :: QINF(:)
    double precision,allocatable  :: SOL_pm_bl(:,:,:,:),RHS_pm_bl(:,:,:,:)
   !double precision              :: Xbound_tmp(6)
   !integer                       :: NN_tmp(3),NN_bl_tmp(6),
    integer                       :: II,itypeb
         
    character *50                 :: outfil

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
    II=1;iynhj=0
    NVR_size=NVRM_in
    NTIME_pm=NTIME_in
    if (II.eq.1) then
        call omp_set_num_threads(1)
       !call mkl_set_num_threads(1)
       !call mkl_domain_set_num_threads(1)
    endif
! if (    PI   = 4.d0 * atan(1.d0)
    PI2  = 2.d0 * PI
    PI4  = 4.d0 * PI
    neqpm=neqpm_in
    NI   =  NI_in
  if (my_rank.eq.0) then 
   nullify(QP)
   nullify(XP)
   nullify(UP)
   nullify(GP)
  QP=>QP_in; XP=>XP_in
  UP=>UP_in; GP=>GP_in
  NVR = NVR_in
   
  else 
   nullify(QP)
   nullify(XP)
   nullify(UP)
   nullify(GP)
  !nullify(velvrx_pm); nullify(velvry_pm); nullify(velvrz_pm)
  endif
  

  nullify (RHS_pm_in)

  nb = my_rank + 1
  
 if (my_rank.eq.0)write(*,*) achar(27)//'[1;31m REALTIME ',NTIME_pm,'WhatToDo = ',WhatTodo,achar(27)//'[0m'
  if (WhatToDo.eq.0) then 
      ND=3
   !  if (NTIME_pm.eq.0) then 
      open(1,file='pm.inp')
      read(1,*) DXpm,DYpm,DZpm
      read(1,*) interf_iproj
      read(1,*) ibctyp
      read(1,*) IDVPM
      read(1,*) EPSVOL
      read(1,*) ncoarse
      read(1,*) NBI,NBJ,NBK
      read(1,*) NREMESH, ncell_rem 
      read(1,*) iyntree,ilevmax
      read(1,*) OMPTHREADS
      read(1,*) idefine
      read(1,*) iynslice
      read(1,*) IPMWRITE
      if(IPMWRITE.gt.10) stop  !maximume value of writes equal to 10
      if(IPMWRITE.GT.0) then
        do i=1,IPMWRITE !max value 10
      read(1,*) IPMWSTART(i), IPMWSTEPS(i)
        enddo
      endif
      close(1)
   !  endif
      call define_sizes
     !if (my_rank.eq.0) then 
          if (allocated(velvrx_pm)) then 
              deallocate(velvrx_pm,velvry_pm,velvrz_pm)
              allocate (velvrx_pm(NXpm,NYpm,NZpm),velvry_pm(NXpm,NYpm,NZpm),velvrz_pm(Nxpm,NYpm,NZpm))
              velvrx_pm=0.d0;velvry_pm=0.d0;velvrz_pm=0.d0
          else 
              allocate (velvrx_pm(NXpm,NYpm,NZpm),velvry_pm(NXpm,NYpm,NZpm),velvrz_pm(Nxpm,NYpm,NZpm))
              velvrx_pm=0.d0;velvry_pm=0.d0;velvrz_pm=0.d0
          endif      

         nullify(velx);nullify(vely);nullify(velz)
         velx=>velvrx_pm; vely=>velvry_pm;velz=>velvrz_pm
         iwrite=0
     !endif
      return
  endif
  call MPI_BCAST(NVR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (NVR.eq.0) return
  NN_tmp(1:3)     = NNbl(1:3,nb)
  NN_bl_tmp(1:6)  = NNbl_bl(1:6,nb)
  allocate(SOL_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3)),RHS_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3)))
  SOL_pm_bl=0.d0;RHS_pm_bl=0.d0

!allocate RHS_pm for parallel projection
  if (allocated(RHS_pm)) then 
      deallocate(RHS_pm)
      allocate(RHS_pm(neqpm+1,NXpm,NYpm,NZpm))
      RHS_pm=0.d0
  else 
      allocate(RHS_pm(neqpm+1,NXpm,NYpm,NZpm))
      RHS_pm=0.d0
  endif
  
! if(my_rank.eq.0) then
     if (allocated(SOL_pm)) then 
         deallocate(SOL_pm)
         allocate(SOL_pm(neqpm,NXpm,NYpm,NZpm))
         SOL_pm=0.d0
     else
         allocate(SOL_pm(neqpm,NXpm,NYpm,NZpm))
         SOL_pm=0.d0
     endif
! endif

  call project_particles 
 !if (my_rank.eq.0) then 
 !    if(mod(NTIME,1).eq.0) call writesol(NTIME-1)
 !endif
  if(my_rank.eq.0)then 
     st=MPI_WTIME()
    !call diffuse_vort_3d
  endif 

if (WhatTodo.lt.4) then 
  call rhsbcast(RHS_pm,NN,neqpm+1)
        !------------------------------
     
  IF (II.ne.0) then
         nb = my_rank + 1
         NN_tmp(1:3)     =  NNbl(1:3,nb)
         call rhsscat_3d(BLOCKS,NN_tmp,NNbl,NNbl_bl,NN_bl,nb_i,nb_j,nb_k,RHS_pm_bl)
  ENDIF
  if(my_rank.eq.0)then 
            et=MPI_WTIME()
            write(*,*) 'RHSBCAST',int((et-st)/60),'m',mod(et-st,60.d0),'s',maxval(RHS_pm_bl(neqpm,:,:,:))
  endif

   if (my_rank.eq.0) st=MPI_WTIME()
   call pmesh_solve
   if (my_rank.eq.0)then
       et= MPI_WTIME()
       write(*,*) 'PMSOLVE',int((et-st)/60),'m',mod(et-st,60.d0),'s'
   endif


  if (WhatToDo.eq.1) then 
     if (my_rank.eq.0)then 
         velvrx_pm=0.d0;velvry_pm=0.d0;velvrz_pm=0.d0

         !call convect_first_order(Xbound,Dpm,NN,NN_bl)
         call calc_velocity_serial_3d(0)
      !  if(mod(NTIME,NWRITE).eq.0) call writesol(NTIME)

         !if (ND.eq.3) then 
         ! call hill_error(NN,NN_bl,Xbound,Dpm,SOL_pm,velvrx_pm,velvry_pm,velvrz_pm)
         ! call writesol(NTIME)
          !stop
         !endif

      !  RHS_pm_in=>RHS_pm
         deallocate(SOL_pm)
         deallocate(RHS_pm)
     endif
     deallocate(SOL_pm_bl,RHS_pm_bl)
     return
  endif

  if (WhatToDo.eq.2) then 
     if (my_rank.eq.0)then 
         !call convect_first_order(Xbound,Dpm,NN,NN_bl)
         st=MPI_WTIME()
         velvrx_pm=0;velvry_pm=0;velvrz_pm=0
         call calc_velocity_serial_3d(1)
         et= MPI_WTIME()
         write(*,*) 'fd',int((et-st)/60),'m',mod(et-st,60.d0),'s'
   !     call calc_antidiffusion
       ! itypeb=1
       ! call back_to_particles_3D(SOL_pm,RHS_pm,XP,QP,UP,GP,&
       !                           velvrx_pm,velvry_pm,velvrz_pm,&
       !                           Xbound,Dpm,NN,NN_bl,NVR,neqpm,interf_iproj,itypeb,NVR_size)
  !      if(mod(NTIME_pm,20).eq.0.or.NTIME_pm.eq.1) call writesol
         if(mod(NTIME_pm,100).eq.0.or.NTIME_pm.eq.1) call writesol
         call writeline

    endif
        itypeb=1
        call back_to_particles_par
         
           iwrite=0
     !if(IPMWRITE.GT.0) then
     !  do i=1,IPMWRITE
     ! if(NTIME.ge.IPMWSTART(i).and.NTIME.le.(IPMWSTART(i)+IPMWSTEPS(i))) call writesol(NTIME)
     !  enddo
     !endif
         !if (ND.eq.3) then 
         ! call hill_error(NN,NN_bl,Xbound,Dpm,SOL_pm,velvrx_pm,velvry_pm,velvrz_pm)
         ! call writesol(NTIME)
          !stop
         !endif

      !  RHS_pm_in=>RHS_pm
         deallocate(SOL_pm)
         deallocate(RHS_pm)
         deallocate(velvrx_pm,velvry_pm,velvrz_pm)
     deallocate(SOL_pm_bl,RHS_pm_bl)
     return
  endif
endif

  if (WhatToDo.eq.4) then 
     if (my_rank.eq.0)then 


         call calc_velocity_serial_3d(-1)
    !    itypeb=1!normal back to particles
    !    call back_to_particles_3D(SOL_pm,RHS_pm,XP,QP,UP,GP,&
    !                              velvrx_pm,velvry_pm,velvrz_pm,&
    !                              Xbound,Dpm,NN,NN_bl,NVR,neqpm,interf_iproj,itypeb,NVR_size)

        if(IPMWRITE.GT.0) then
          do i=1,IPMWRITE
          if(NTIME_pm.ge.IPMWSTART(i).and.NTIME_pm.le.(IPMWSTART(i)+IPMWSTEPS(i))) call writesol
          enddo
        endif
        if(iynslice.eq.1)call writesolXavatar
      endif
        itypeb=1
        call back_to_particles_par
         
        iwrite=0
      !  RHS_pm_in=>RHS_pm
      !  velx=>velvrx_pm; vely=>velvry_pm;velz=>velvrz_pm
         deallocate(SOL_pm)
         deallocate(RHS_pm)
         deallocate(SOL_pm_bl,RHS_pm_bl)
     return
  endif
  
 if (WhatToDo.eq.5) then 
     if (my_rank.eq.0) then 
 !diffusion stores -NI*grad^2 w * Vol in GP(1,:)
         
         call diffuse_vort_3d
       ! itypeb=2!back to particles the diffused vorticity
       ! call back_to_particles_3D(SOL_pm,RHS_pm,XP,QP,UP,GP,&
       !                           velvrx_pm,velvry_pm,velvrz_pm,&
       !                           Xbound,Dpm,NN,NN_bl,NVR,neqpm,interf_iproj,itypeb,NVR_size)
      endif
        itypeb=2
        call back_to_particles_par
       deallocate(SOL_pm,RHS_pm)
       deallocate(SOL_pm_bl,RHS_pm_bl)
     return

 endif 
contains
     Subroutine pmesh_solve
            !Yaps or Serial Pmesh
            IF(II.eq.0) then
                IF(my_rank.eq.0) then
                        write(*,*) 'Solving_pm'
                        SOL_pm(1:neqpm,:,:,:)=0.0
                        itree=0
                        iynbc=1!for infinite domain bc's
                        call pmesh(SOL_pm,RHS_pm,QP,XP,Xbound,DPm,NN,NN_bl,ND,Nblocks,ibctyp,1,3,&
                                   iynbc,NVR,itree,ilevmax)
                      ! call calc_velocity_serial_3d(1)
                        write(*,*) 'PM_solved'
                ENDIF
                !--------------------------------------------
                    call velbcast_3d
            ELSE
                iret =0
                call yaps3d(SOL_pm_bl,RHS_pm_bl,Xbound_bl,Xbound_coarse,Dpm,Dpm_coarse,NNbl,NNbl_bl,&
                    NN_coarse,NN_bl_coarse,ND,BLOCKS,ibctyp,1,neqpm,ncoarse,NBI,NBJ,NBK,nb_i,nb_j,nb_k,&
                    iret,iyntree,ilevmax,neqpm)
                
                nb              =  my_rank + 1
                NN_tmp(1:3)     =  NNbl(1:3,nb)
                NN_bl_tmp(1:6)  =  NNbl_bl(1:6,nb)
                Xbound_tmp(1:6) =  Xbound_bl(1:6,nb)
                call solget_3d(BLOCKS,NBI,NBJ,NBK,NN_tmp,NNbl,NNbl_bl,NN_bl,SOL_pm_bl)
               !if (my_rank.eq.0) call calc_velocity_serial_3d(1)
               !call velbcast_3d
            ENDIF
            !--------------------------------------------

     End Subroutine pmesh_solve

     
     Subroutine project_particles
            call MPI_BCAST(NVR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            NVR_p=NVR/np
            if (my_rank.eq.0) NVR_p = NVR_p + mod(NVR,np)
            allocate (XP_scatt(3,NVR_p),QP_scatt(neqpm+1,NVR_p),NVR_projscatt(NVR_p))
            NVR_projscatt=interf_iproj
            call particles_scat
            call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
            if(my_rank.eq.0) call omp_set_num_threads(1)
                st=MPI_WTIME()
                allocate(ieq(neqpm+1),QINF(neqpm+1))
                QINF=0.d0
                do i=1,neqpm+1
                   ieq(i)=i
                enddo
                call project_particles_3D(RHS_pm,QP_scatt,XP_scatt,NVR_projscatt,NVR_p,neqpm+1,ieq,neqpm+1,QINF,NVR_p)

                call proj_gath(NN)
                if (my_rank.eq.0) then 
                  !allocate(ieq(neqpm+1),QINF(neqpm+1))
                  !QINF=0.d0
                  !do i=1,neqpm+1
                  !   ieq(i)=i
                  !enddo
                  !allocate (NVR_projscatt(NVR_size))
                  !NVR_projscatt=interf_iproj
                  !call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
                  !call project_particles_3D(RHS_pm,QP,XP,NVR_projscatt,NVR,neqpm+1,ieq,neqpm+1,QINF,NVR_size)
                    call omp_set_num_threads(OMPTHREADS)
                    call project_vol3d(RHS_pm,neqpm+1,ieq,neqpm+1,IDVPM)
                    et=MPI_WTIME()
                    write(*,*) 'Proj',int((et-st)/60),'m',mod(et-st,60.d0),'s'
                    call omp_set_num_threads(1)
              ! deallocate(ieq,QINF)
              ! deallocate(NVR_projscatt)
                endif
                deallocate(ieq,QINF)
                deallocate (XP_scatt,QP_scatt,NVR_projscatt)

     End Subroutine project_particles 

     Subroutine back_to_particles_par
       if (my_rank.eq.0) st=MPI_WTIME()
         call rhsbcast(RHS_pm,NN,neqpm+1)
         call rhsbcast(SOL_pm,NN,neqpm)
         call velbcast_3d
         call MPI_BCAST(NVR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         NVR_p=NVR/np
         if (my_rank.eq.0) NVR_p = NVR_p + mod(NVR,np)
         allocate (XP_scatt(3,NVR_p),QP_scatt(neqpm+1,NVR_p),&
                   UP_scatt(3,NVR_p),GP_scatt(3,NVR_p))
      
         UP_scatt=0.d0;GP_scatt=0.d0
         call particles_scat
         if (my_rank.eq.0)then
             et= MPI_WTIME()
               write(*,*) 'backcomm',int((et-st)/60),'m',mod(et-st,60.d0),'s'
         endif
        if (my_rank.eq.0) st=MPI_WTIME()
         call back_to_particles_3D(SOL_pm,RHS_pm,XP_scatt,QP_scatt,UP_scatt,GP_scatt,&
                                   velvrx_pm,velvry_pm,velvrz_pm,&
                                   Xbound,Dpm,NN,NN_bl,NVR_p,neqpm,interf_iproj,itypeb,NVR_p)
          call particles_gath
         deallocate (XP_scatt,QP_scatt,UP_scatt,GP_scatt)
         if (my_rank.eq.0)then
             et= MPI_WTIME()
               write(*,*) 'back',int((et-st)/60),'m',mod(et-st,60.d0),'s'
         endif

     End Subroutine back_to_particles_par
End Subroutine vpm

Subroutine define_sizes
  use vpm_vars        
  use vpm_size
  use pmlib
  use pmeshpar
  use parvar
  use pmgrid
  use MPI
  integer    :: nsiz(3),nsiz_bl(3)

  integer    :: i,j,k,np,my_rank,ierr,nb
 !double precision              :: Xbound_tmp(6)
 !integer                       :: NN_tmp(3),NN_bl_tmp(6)
 
    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)


            BLOCKS=np
!-------First Change Dpm so that the numbers of cells divides
!-------by nsize i.e with NBI,NBJ,ncoarse,levmax depending on the criterion
!------that's why ndumcell=0
            if(my_rank.eq.0) then 
               if ((NTIME_pm.eq.0.and.idefine.eq.1).or.idefine.eq.0) then 
                 XMIN_pm=minval(XP(1,1:NVR)) - interf_iproj*DXpm 
                 YMIN_pm=minval(XP(2,1:NVR)) - interf_iproj*DYpm 
                 ZMIN_pm=minval(XP(3,1:NVR)) - interf_iproj*DZpm 
                
                 XMAX_pm=maxval(XP(1,1:NVR)) + interf_iproj*DXpm 
                 YMAX_pm=maxval(XP(2,1:NVR)) + interf_iproj*DYpm 
                 ZMAX_pm=maxval(XP(3,1:NVR)) + interf_iproj*DZpm 
              !else
              !  XMAX_pm=maxval(XP(1,1:NVR)) + interf_iproj*DXpm 
               endif
            endif
   
            call MPI_BCAST(XMIN_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(YMIN_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(ZMIN_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

            call MPI_BCAST(XMAX_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(YMAX_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(ZMAX_pm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
            Xbound(1)=XMIN_pm;Xbound(2)=YMIN_pm;Xbound(3)=ZMIN_pm
            Xbound(4)=XMAX_pm;Xbound(5)=YMAX_pm;Xbound(6)=ZMAX_pm
            Dpm   (1)=DXpm   ;Dpm(2)   =DYpm   ;Dpm(3)   =DZpm
            ndumcell=0
            nsiz(1) = NBI * ncoarse
            nsiz(2) = NBJ * ncoarse
            nsiz(3) = NBK * ncoarse
             if ((NTIME_pm.eq.0.and.idefine.eq.1).or.idefine.eq.0) then
              call  definepm(3,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
           !else
           !  call  definepm(4,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
            endif               
            XMIN_pm=Xbound(1);YMIN_pm=Xbound(2);ZMIN_pm=Xbound(3)
            XMAX_pm=Xbound(4);YMAX_PM=Xbound(5);ZMAX_pm=Xbound(6)
            NXpm=NN(1);NYpm=NN(2);NZpm=NN(3)
            NXs_bl(1)=NN_bl(1);NYs_bl(1)=NN_bl(2);NZs_bl(1)=NN_bl(3)
            NXf_bl(1)=NN_bl(4);NYf_bl(1)=NN_bl(5);NZf_bl(1)=NN_bl(6)
            DXpm = Dpm(1); DYpm = Dpm(2); DZpm = Dpm(3)
            if (ND.eq.2) then
                DVpm=DXpm*DYpm
                DXpm2=2*DXpm
                DYpm2=2*DYpm
            else
                DVpm=DXpm*DYpm*DZpm
                DXpm2=2*DXpm
                DYpm2=2*DYpm
                DZpm2=2*DZpm
            endif
            if (my_rank.eq.0) then 
                write(*,*) 'New Dpm(1),Dpm(2),Dpm(3)'
                write(*,*) Dpm(1),Dpm(2),Dpm(3)
                write(*,*) NN
                write(*,*) Xbound
             endif
!define block grid so the are divided by ncoarse and ilevmax
!so as to have the coarse information at the boundaries exactly.
            ndumcell_bl= ncoarse
            nsiz_bl=ncoarse!ndumcell_coarse!*2*2**ilevmax
              if (.not. allocated(XBound_bl)) then 
                  allocate (Xbound_bl(6,BLOCKS),NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS))
              endif
                 NXB   = int(nint(((Xbound (4)-Xbound (1))/Dpm (1))))
                 NYB   = int(nint(((Xbound (5)-Xbound (2))/Dpm (2))))
                 NZB   = int(nint(((Xbound (6)-Xbound (3))/Dpm (3))))
                 NXbl  = NXB/NBI
                 NYbl  = NYB/NBJ
                 NZbl  = NZB/NBK
                 if (my_rank.eq.0) write(*,*) 'Check sizes',NXbl,float(NXB)/NBI,NYbl,float(NYB)/NBJ,float(NZB)/NBK,ncoarse
                 if (my_rank.eq.0) write(*,*) 'Check sizes coarse',float(NXbl)/ncoarse,float(NYbl)/ncoarse,float(NZbl)/ncoarse
                 do k=1,NBK
                    do j=1,NBJ
                        do i=1,NBI
                            nb  = (k-1)*NBJ*NBI + (j-1)*NBI + i
                            Xbound_bl(1,nb)=Xbound (1) + (i-1)*(NXbl)*Dpm (1)
                            Xbound_bl(4,nb)=Xbound (1) + (i)  *(NXbl)*Dpm (1)
      
                            Xbound_bl(2,nb)=Xbound (2) + (j-1)*(NYbl)*Dpm (2)
                            Xbound_bl(5,nb)=Xbound (2) + (j)  *(NYbl)*Dpm (2)
                    
                    
                            Xbound_bl(3,nb)=Xbound (3) + (k-1)*(NZbl)*Dpm (3)
                            Xbound_bl(6,nb)=Xbound (3) + (k)  *(NZbl)*Dpm (3)
                            !write(*,*) 'Block',nb
                            !write(*,*) Xbound_bl(1,nb),Xbound_bl(4,nb),Xbound_bl(2,nb),Xbound_bl(5,nb)
                            Xbound_tmp(1:6) = Xbound_bl(1:6,nb)
                            call definepm(1,Xbound_tmp,Dpm ,ND,ndumcell_bl,nsiz_bl,NN_tmp,NN_bl_tmp)
                            Xbound_bl(1:6,nb) = Xbound_tmp(1:6)
                            NNbl(1:3,nb)      = NN_tmp(1:3)
                            NNbl_bl(1:6,nb)   = NN_bl_tmp(1:6)
                            if (nb.eq.my_rank+1) then
                                nb_i=i
                                nb_j=j
                                nb_k=k
                            endif
                        enddo
                    enddo
                 enddo
        
            !--B
            nb = my_rank+1
!define coarse grid must cover block grids
            
            Xbound(1)=XMIN_pm!minval(Xbound_bl(1,:))
            Xbound(2)=YMIN_pm!minval(Xbound_bl(2,:))
            Xbound(3)=ZMIN_pm!minval(Xbound_bl(3,:))
            Xbound(4)=XMAX_pm!maxval(Xbound_bl(4,:))
            Xbound(5)=YMAX_pm!maxval(Xbound_bl(5,:))
            Xbound(6)=ZMAX_pm!maxval(Xbound_bl(6,:))
            Xbound_coarse = Xbound
            Dpm_coarse    = ncoarse*Dpm
            ndumcell_coarse = 4!2**ilevmax
            nsiz_bl=2**ilevmax
            call definepm(1,Xbound_coarse,Dpm_coarse,ND,ndumcell_coarse,nsiz_bl,NN_coarse,NN_bl_coarse)
!add to dummy cells to the grid globally used for remeshing purposes mainly
           
           !ndumcell=4
           !Xbound(1)=XMIN_pm;Xbound(2)=YMIN_pm;Xbound(3)=ZMIN_pm
           !Xbound(4)=XMAX_pm;Xbound(5)=YMAX_pm;Xbound(6)=ZMAX_pm
           !Dpm   (1)=DXpm   ;Dpm(2)   =DYpm   ;Dpm(3)   =DZpm
           !call  definepm(1,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
           !XMIN_pm=Xbound(1);YMIN_pm=Xbound(2);ZMIN_pm=Xbound(3)
           !XMAX_pm=Xbound(4);YMAX_PM=Xbound(5);ZMAX_pm=Xbound(6)
           !NXpm=NN(1);NYpm=NN(2);NZpm=NN(3)
           !NXs_bl(1)=NN_bl(1);NYs_bl(1)=NN_bl(2);NZs_bl(1)=NN_bl(3)
           !NXf_bl(1)=NN_bl(4);NYf_bl(1)=NN_bl(5);NZf_bl(1)=NN_bl(6)
           !print *,'final mesh',NN
            !-----
        
            
     
     


End Subroutine define_sizes

Subroutine writesol
    use vpm_vars        
    use pmeshpar
    use parvar
    use pmgrid
    use MPI

    
    character*50        :: filout
    integer           :: i,j,k
    double precision  :: XPM,YPM,ZPM,velocx,velocy,velocz
    
      ! if(iwrite.ne.0) return
        write(filout,'(i5.5,a)') NTIME_pm,'solution.dat'
        open(1,file=filout)
        WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
        WRITE(1,*)'ZONE I=',NXf_bl(1)-NXs_bl(1)+1,' J=',NYf_bl(1) - NYs_bl(1) + 1,&
            ' K=',NZf_bl(1) - NZs_bl(1) +1 ,' F=POINT'
        do k=NZs_bl(1),NZf_bl(1)
            do j=NYs_bl(1),NYf_bl(1)
                do i=NXs_bl(1),NXf_bl(1)
                    ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
                    ! do j=1,NYpm
                    !   do i=1,NXpm
                    XPM=XMIN_pm+(I-1)*DXpm
                    YPM=YMIN_pm+(J-1)*DYpm
                    ZPM=ZMIN_pm+(K-1)*DZpm
                    velocx = VelvrX_pm(i,j,k)
                    velocy = VelvrY_pm(i,j,k)
                    velocz = VelvrZ_pm(i,j,k)

                    WRITE(1,'(16(e28.17,1x))')XPM,YPM,ZPM,velocx,velocy,velocz,-RHS_pm(1,I,J,K),&
                                             -RHS_pm(2,I,J,K),&
                                             -RHS_pm(3,I,J,K)!,RHS_pm(4,I,J,K),SOL_pm(1,I,J,K),SOL_pm(2,I,J,K), SOL_pm(3,I,J,K)


                enddo
            enddo
        enddo
        close(1)
    iwrite=1
    !   ---FOR PLOTTING PURPOSES ONLY
    call system('~/bin/preplot '//filout//' >/dev/null')
 
  ! call system('rm '//filout)

  return
    write(filout,'(i5.5,a)') NTIME_pm,'flowfield.dat'
    open(1,file=filout,form='unformatted')
    write(1) NXs_bl(1),NXf_bl(1)
    write(1) NYs_bl(1),NYf_bl(1)
    write(1) NZs_bl(1),NZf_bl(1)
    write(1) DXpm,DYpm,DZpm
    write(1) XMIN_pm,YMIN_pm,ZMIN_pm
    do k=NZs_bl(1),NZf_bl(1)
        do j=NYs_bl(1),NYf_bl(1)
            do i=NXs_bl(1),NXf_bl(1)
                XPM=XMIN_pm+(I-1)*DXpm
                YPM=YMIN_pm+(J-1)*DYpm
                ZPM=ZMIN_pm+(K-1)*DZpm
                velocx = VelvrX_pm(i,j,k)
                velocy = VelvrY_pm(i,j,k)
                velocz = VelvrZ_pm(i,j,k)

                WRITE(1)       XPM,YPM,ZPM,velocx,velocy,velocz,-RHS_pm(1,I,J,K),&
                               -RHS_pm(2,I,J,K),&
                               -RHS_pm(3,I,J,K)

            enddo
        enddo
    enddo

    close(1)

End Subroutine writesol

Subroutine writesolXavatar

! AVATAR T2.4 specific: rotor + turbulent particles inflow [PM]

    use vpm_vars        
    use pmeshpar
    use parvar
    use pmgrid
    use MPI

    character*50      :: filout
    integer           :: i,j,k
    double precision  :: XPM,YPM,ZPM,velocx,velocy,velocz
    integer,dimension(15) :: NX_AVA
    integer               :: NXPOS_AVA_512_128_128, ii
    
     NXPOS_AVA_512_128_128=15
     NX_AVA(1)=2
     NX_AVA(2)=15
     NX_AVA(3)=40
     NX_AVA(4)=65
     NX_AVA(5)=90
     NX_AVA(6)=103
     NX_AVA(7)=128
     NX_AVA(8)=154
     NX_AVA(9)=205
     NX_AVA(10)=256
     NX_AVA(11)=308
     NX_AVA(12)=365
     NX_AVA(13)=417
     NX_AVA(14)=470
     NX_AVA(15)=512
     write(filout,'(i5.5,a)') NTIME_pm,'solX.dat'
     open(1,file=filout)

     WRITE(1,'(a100)')'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'! "Vol" "PSIX" "PSIY" "PSIZ"'
     do ii=1,NXPOS_AVA_512_128_128
        i=NX_AVA(ii) !/2  PM: DX=DY=DZ=8m
        WRITE(1,'(a11,i3,a8,i4,a7,i4,a7,i4,a11)') 'ZONE T= "',i,'", I=',1,&
                                             ', J=',NYf_bl(1)-NYs_bl(1)+1, &
                                             ', K=',NZf_bl(1)-NZs_bl(1)+1,  ', F=POINT'
        do k=NZs_bl(1),NZf_bl(1)
        do j=NYs_bl(1),NYf_bl(1)
             XPM=XMIN_pm+(I-1)*DXpm
             YPM=YMIN_pm+(J-1)*DYpm
             ZPM=ZMIN_pm+(K-1)*DZpm
             velocx = VelvrX_pm(i,j,k)
             velocy = VelvrY_pm(i,j,k)
             velocz = VelvrZ_pm(i,j,k)
             WRITE(1,'(16(e28.17,1x))')XPM,YPM,ZPM,velocx,velocy,velocz,-RHS_pm(1,I,J,K),&
                                      -RHS_pm(2,I,J,K),&
                                      -RHS_pm(3,I,J,K)!,RHS_pm(4,I,J,K),SOL_pm(1,I,J,K),SOL_pm(2,I,J,K), SOL_pm(3,I,J,K)
        enddo !j
        enddo !k
     enddo !ii
     close(1)
    call system('gzip '//filout)

  return

End Subroutine writesolXavatar

Subroutine writeline
    use vpm_vars        
    use pmeshpar
    use parvar
    use pmgrid
    use MPI

    
    character*50        :: filout
    integer           :: i,j,k,jmat(9),kmat(9),NNJ,NNK,il
    double precision  :: XPM,YPM,ZPM,velocx,velocy,velocz


    NNJ= NYf_bl(1) - NYs_bl(1)+1
    NNK= NZf_bl(1) - NZs_bl(1)+1
    jmat(1)=0.5 *NNJ;kmat(1)=0.5 *NNK
    jmat(2)=0.25*NNJ;kmat(2)=0.25*NNK
    jmat(3)=0.25*NNJ;kmat(3)=0.75*NNK
    jmat(4)=0.75*NNJ;kmat(4)=0.25*NNK
    jmat(5)=0.75*NNJ;kmat(5)=0.75*NNK
    jmat(6)=0.5*NNJ ;kmat(6)=0.25*NNK
    jmat(7)=0.5*NNJ;kmat(7)=0.75*NNK
    jmat(8)=0.25*NNJ;kmat(8)=0.5*NNK
    jmat(9)=0.75*NNJ;kmat(9)=0.5*NNK
    do i= 1, 9 
       j=jmat(i);k=kmat(i)
       write(filout,'(i2.2,a)') i,'hist.bin'
       open(1,file=filout,access='APPEND',form='unformatted')
       if (NTIME_pm.eq.1) then
          rewind(1)
          write(1) XMIN_pm,YMIN_pm,ZMIN_pm
          write(1) DXpm,DYpm,DZpm
          write(1) NXs_bl(1),NXf_bl(1),j,k
       endif
       WRITE(1) NTIME_pm,(velvrx_pm(il,j,k),il=NXs_bl(1),NXf_bl(1)),&
                        (velvry_pm(il,j,k),il=NXs_bl(1),NXf_bl(1)),&
                        (velvrz_pm(il,j,k),il=NXs_bl(1),NXf_bl(1))
                          
                              

       close(1)
    enddo


End Subroutine writeline

End Module vpm_lib
