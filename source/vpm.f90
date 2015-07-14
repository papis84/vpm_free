 Module vpm_vars
    double precision, allocatable :: XP_scatt(:,:),QP_scatt(:,:)
    double precision              :: DT_c,V_ref
    integer ,allocatable          :: NVR_projscatt(:)
    integer                       :: interf_iproj,ncell_rem,iynhj

    integer                       :: ncoarse,nparcell1d
    integer                       :: neqpm,NVR_p

 End Module vpm_vars
 
 Module vpm_size
    double precision :: Xbound(6),Dpm(3),Xbound0(6),Dpm0(3)
    integer          :: NN_bl(6),NN(3),NN0_bl(6),NN0(3)
    integer          :: NXs0_bl(10),NYs0_bl(10),NXf0_bl(10),NYf0_bl(10),NZs0_bl(10),NZf0_bl(10)

    double precision,allocatable  :: Xbound_bl(:,:)
    integer         ,allocatable  :: NNbl_bl(:,:),NNbl(:,:)


    double precision              :: Xbound_tmp(6),Xbound_coarse(6),Dpm_coarse(3)
    integer                       :: NN_tmp(3),NN_bl_tmp(6),NN_coarse(3),NN_bl_coarse(6)
    integer                       :: nb_i,nb_j,nb_k,NBB,NXbl,NYbl,NZbl,BLOCKS,NXB,NYB,NZB,ndumcell_coarse ,ndumcell_bl
    double precision              :: starttime,endtime,st,et,ct
    integer                       :: II,iynbc,iret,NBI,NBJ,NBK,NVR_out_thres,NREMESH,ntorder,&
                                     iyntree,ilevmax,itree,nsize_out(3),ibctyp

 End Module vpm_size

 Module openmpth
   integer                       ::OMPTHREADS
 End Module openmpth
Module vpm_lib
  use vpm_vars
  use vpm_size
  use openmpth
  
contains

Subroutine vpm_init(XP_in,QP_in,UP_in,GP_in,NVR_in,neqpm_in)
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
    double precision,intent(inout),target :: XP_in(:,:),QP_in(:,:),UP_in(:,:),GP_in(:,:)
    integer,         intent(inout)        :: NVR_in
    integer,         intent(in)           :: neqpm_in
    integer          :: ierr,my_rank,np

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
    II=1;iynhj=0
    if (II.eq.1) then
        call omp_set_num_threads(1)
       !call mkl_set_num_threads(1)
       !call mkl_domain_set_num_threads(1)
    endif
    PI   = 4.d0 * atan(1.d0)
    PI2  = 2.d0 * PI
    PI4  = 4.d0 * PI
    DT_c=0.2d0
  neqpm=neqpm_in
  if (my_rank.eq.0) then 
  QP=>QP_in; XP=>XP_in
  UP=>UP_in; GP=>GP_in
  NVR = NVR_in
  else 
   nullify(QP)
   nullify(XP)
   nullify(UP)
   nullify(GP)
  endif

  ND=3
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
  close(1)

End Subroutine vpm_init

 Subroutine vpm(XP_in,QP_in,UP_in,GP_in,NVR_in,neqpm_in,NTIME)
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
    
    double precision,intent(inout),target :: XP_in(:,:),QP_in(:,:),UP_in(:,:),GP_in(:,:)
    integer,         intent(inout)        :: NVR_in
    integer,         intent(in)           :: neqpm_in,NTIME
    double precision :: a,XPM,YPM,velx,vely,totmass,totvor,MACH,error,pr
    integer          :: i , j, k , nb, NXs,&
                         NYs, NXf, NYf, NZs,NZf,ivel,ij,iter,IPM,ialloc
    integer,allocatable :: ieq(:)

    integer          :: ijk,npar,itype,nel,iflag,inf,iface,rk_step
    integer          :: ierr,my_rank,mat4,mat3,np,ICONV,omp_get_num_threads,omp_get_max_threads
    double precision :: xi,yi,ksi1,ksi2,th1,th2,w1,w2,xl,yl,XO(3),r
    double precision              :: errorCFD,errorf
    double precision ,allocatable :: QINF(:)
    double precision,allocatable  :: SOL_pm_bl(:,:,:,:),RHS_pm_bl(:,:,:,:)
    double precision              :: Xbound_tmp(6)
    integer                       :: NN_tmp(3),NN_bl_tmp(6)
         
    character *50                 :: outfil

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
   
! if (NTIME.eq.1) call define_sizes(NTIME)
  nb = my_rank + 1
  NN_tmp(1:3)     = NNbl(1:3,nb)
  NN_bl_tmp(1:6)  = NNbl_bl(1:6,nb)
  allocate(SOL_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3)),RHS_pm_bl(neqpm,NN_tmp(1),NN_tmp(2),NN_tmp(3)))
  allocate(RHS_pm(4,NXpm,NYpm,NZpm))
  if (my_rank.eq.0)write(*,*) achar(27)//'[1;31m REALTIME',NTIME,achar(27)//'[0m'

  if(my_rank.eq.0) then
     allocate(SOL_pm(neqpm,NXpm,NYpm,NZpm))
  endif

  call project_particles 
  if(my_rank.eq.0)st=MPI_WTIME()
  call rhsbcast(RHS_pm,NN,4)
        !------------------------------
  allocate(velvrx_pm(NXpm,NYpm,NZpm),velvry_pm(NXpm,NYpm,NZpm),velvrz_pm(Nxpm,NYpm,NZpm))
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

  if (my_rank.eq.0) then 
      if(mod(NTIME,10).eq.0) call writesol(NTIME)
      call convect_first_order(Xbound,Dpm,NN,NN_bl)
  endif

  deallocate(Velvrx_pm,Velvry_pm,Velvrz_pm)
  deallocate(RHS_pm,SOL_pm_bl,RHS_pm_bl)
  if (my_rank.eq.0) deallocate(SOL_pm)
  
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
                        call calc_velocity_serial_3d
                        write(*,*) 'PM_solved'
                ENDIF
                !--------------------------------------------
                    call velbcast_3d
            ELSE
                iret =0
                call yaps3d(SOL_pm_bl,RHS_pm_bl,Xbound_bl,Xbound_coarse,Dpm,Dpm_coarse,NNbl,NNbl_bl,&
                    NN_coarse,NN_bl_coarse,ND,BLOCKS,ibctyp,1,3,ncoarse,NBI,NBJ,NBK,nb_i,nb_j,nb_k,&
                    iret,iyntree,ilevmax,neqpm)
                
                nb              =  my_rank + 1
                NN_tmp(1:3)     =  NNbl(1:3,nb)
                NN_bl_tmp(1:6)  =  NNbl_bl(1:6,nb)
                Xbound_tmp(1:6) =  Xbound_bl(1:6,nb)

                call solget_3d(BLOCKS,NBI,NBJ,NBK,NN_tmp,NNbl,NNbl_bl,NN_bl,SOL_pm_bl)
                if (my_rank.eq.0) call calc_velocity_serial_3d
               !call velbcast_3d
            ENDIF
            !--------------------------------------------

     End Subroutine pmesh_solve

     
     Subroutine project_particles
            call MPI_BCAST(NVR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
            NVR_p=NVR/np
            if (my_rank.eq.0) NVR_p = NVR_p + mod(NVR,np)
            allocate (XP_scatt(3,NVR_p),QP_scatt(4,NVR_p),NVR_projscatt(NVR_p))
            NVR_projscatt=interf_iproj
            call particles_scat
            call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
            if(my_rank.eq.0) call omp_set_num_threads(1)
                st=MPI_WTIME()
                allocate(ieq(4),QINF(4))
                QINF(1)=0.d0;QINF(2)=0.d0;QINF(3)=0.d0;QINF(4)=0!DVpm 
                ieq(1)=1;ieq(2)=2;ieq(3)=3;ieq(4)=4
                call project_particles_3D(RHS_pm,QP_scatt,XP_scatt,NVR_projscatt,NVR_p,4,ieq,4,QINF)
                call proj_gath(NN)
                if (my_rank.eq.0) then 
                    RHS_pm(4,:,:,:)= RHS_pm(4,:,:,:) !+ DVpm
                    call omp_set_num_threads(OMPTHREADS)
                    call project_vol3d(RHS_pm,4,ieq,4,IDVPM)
                    et=MPI_WTIME()
                    write(*,*) 'Proj',int((et-st)/60),'m',mod(et-st,60.d0),'s'
                    call omp_set_num_threads(1)
                endif
                deallocate(ieq,QINF)
                deallocate (XP_scatt,QP_scatt,NVR_projscatt)

     End Subroutine project_particles 

End Subroutine vpm

Subroutine define_sizes(NTIME)
  use vpm_vars        
  use vpm_size
  use pmlib
  use pmeshpar
  use parvar
  use pmgrid
  use MPI
  integer, intent(in) ::NTIME
  integer    :: nsiz(3),nsiz_bl(3)

  integer    :: i,j,k,np,my_rank,ierr,nb
  double precision              :: Xbound_tmp(6)
  integer                       :: NN_tmp(3),NN_bl_tmp(6)
 
    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)


            BLOCKS=np
!-------First Change Dpm so that the numbers of cells divides
!-------by nsize i.e with NBI,NBJ,ncoarse,levmax depending on the criterion
!------that's why ndumcell=0
            if(my_rank.eq.0) then 
               XMIN_pm=minval(XP(1,1:NVR)) - 2*interf_iproj*DXpm 
               YMIN_pm=minval(XP(2,1:NVR)) - 2*interf_iproj*DYpm 
               ZMIN_pm=minval(XP(3,1:NVR)) - 2*interf_iproj*DZpm 
              
               XMAX_pm=maxval(XP(1,1:NVR)) + 2*interf_iproj*DXpm 
               YMAX_pm=maxval(XP(2,1:NVR)) + 2*interf_iproj*DYpm 
               ZMAX_pm=maxval(XP(3,1:NVR)) + 2*interf_iproj*DZpm 
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
            if (iyntree.eq.0) ilevmax=1
            nsiz(1) = max((2**(ilevmax)*NBI),NBI*ncoarse)
            nsiz(2) = max((2**(ilevmax)*NBJ),NBJ*ncoarse)
            nsiz(3) = max((2**(ilevmax)*NBK),NBK*ncoarse)
            call  definepm(2,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
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
             endif
!define block grid so the are divided by ncoarse and ilevmax
!so as to have the coarse information at the boundaries exactly.
            ndumcell_bl= 4
            nsiz_bl=max(2**ilevmax,ncoarse)!ndumcell_coarse!*2*2**ilevmax
            if (NTIME.eq.1) allocate (Xbound_bl(6,BLOCKS),NNbl(3,BLOCKS),NNbl_bl(6,BLOCKS))
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
                            call definepm(3,Xbound_tmp,Dpm ,ND,ndumcell_bl,nsiz_bl,NN_tmp,NN_bl_tmp)
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
            Xbound_tmp(1:6) = Xbound_bl(1:6,nb)

!define coarse grid must cover block grids
            Xbound(1:3)=minval(Xbound_bl(1:3,:))
            Xbound(4:6)=maxval(Xbound_bl(4:6,:))
            Xbound_coarse = Xbound
            Dpm_coarse    = ncoarse*Dpm
            ndumcell_coarse = 4!2**ilevmax
            nsiz_bl=2**(ilevmax-1)!ndumcell_coarse!*2*2**ilevmax
            call definepm(1,Xbound_coarse,Dpm_coarse,ND,ndumcell_coarse,nsiz,NN_coarse,NN_bl_coarse)
             
!add to dummy cells to the grid globally used for remeshing purposes mainly
            ndumcell=4
            Xbound(1)=XMIN_pm;Xbound(2)=YMIN_pm;Xbound(3)=ZMIN_pm
            Xbound(4)=XMAX_pm;Xbound(5)=YMAX_pm;Xbound(6)=ZMAX_pm
            Dpm   (1)=DXpm   ;Dpm(2)   =DYpm   ;Dpm(3)   =DZpm
            call  definepm(1,Xbound,Dpm,ND,ndumcell,nsiz,NN,NN_bl)
            XMIN_pm=Xbound(1);YMIN_pm=Xbound(2);ZMIN_pm=Xbound(3)
            XMAX_pm=Xbound(4);YMAX_PM=Xbound(5);ZMAX_pm=Xbound(6)
            NXpm=NN(1);NYpm=NN(2);NZpm=NN(3)
            NXs_bl(1)=NN_bl(1);NYs_bl(1)=NN_bl(2);NZs_bl(1)=NN_bl(3)
            NXf_bl(1)=NN_bl(4);NYf_bl(1)=NN_bl(5);NZf_bl(1)=NN_bl(6)
            !-----
        
            
     
     


End Subroutine define_sizes

Subroutine writesol(NTIME)
    use vpm_vars        
    use pmeshpar
    use parvar
    use pmgrid
    use MPI

    
    character*50        :: filout
    integer           :: i,j,k,NTIME
    double precision  :: XPM,YPM,ZPM,velx,vely,velz
    
         
        write(filout,'(i5.5,a)') NTIME,'solution.dat'
        open(1,file=filout)
        WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ" "Vol" "PSIX" "PSIY" "PSIZ"'
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
                    !velx = VelphiX_pm(i,j,1) + VelvrX_pm(i,j,1)
                    !vely = VelphiY_pm(i,j,1) + VelvrY_pm(i,j,1)
                    velx = VelvrX_pm(i,j,k)
                    vely = VelvrY_pm(i,j,k)
                    velz = VelvrZ_pm(i,j,k)

                    WRITE(1,'(16(e28.17,1x))')XPM,YPM,ZPM,velx,vely,velz,-RHS_pm(1,I,J,K),-RHS_pm(2,I,J,K),&
                                             -RHS_pm(3,I,J,K),RHS_pm(4,I,J,K),SOL_pm(1,I,J,K),SOL_pm(2,I,J,K),&
                                              SOL_pm(3,I,J,K)

                enddo
            enddo
        enddo
        close(1)

    !   ---FOR PLOTTING PURPOSES ONLY
    call system('~/bin/preplot '//filout//' >/dev/null')
 
    call system('rm '//filout)



End Subroutine writesol

End Module vpm_lib
