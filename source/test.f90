
Program test_pm
use vpm_lib
use test_mod
use pmgrid
use MPI
double precision :: Vref,NI_in,DT_in,RMETM,OMET,OG,FACDEF,T,XMIN,XMAX,UINF(3)
double precision,allocatable ::velsavex(:,:,:)
double precision,allocatable ::XPDUM(:,:),QPDUM(:,:)
integer          :: Noutput, NDimoutput,NVR_turb
integer :: my_rank,np,ierr,i,neq,j

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 NI_in=-0.1
 DT_in=0.761905/2.  !=dx/U=8/10.5   !1000
 DT_in=3.76     !=dx/U=8/10.5   !1000

 if (my_rank.eq.0) then 
    open(1,file='particles.bin',form='unformatted')
    read(1) NVR_ext
    read(1) Vref
    allocate(XPR(3,NVR_ext),QPR(4,NVR_ext))
    write(*,*) 'NVR=',NVR_ext,Vref
    do i=1,NVR_ext
       read(1) XPR(1,i),XPR(2,i),XPR(3,i),QPR(1,i),QPR(2,i),QPR(3,i)
    enddo

 QPR(1:3,:) = -QPR(1:3,:) * Vref
!QPR(1:3,:)=0; QPR(1,:) =-0.1d0*Vref
 QPR(4,:) =Vref
 RMETM=0.001
endif
mrem=1
neq=3
UINF=0;UINF(1)=10.5
 !-Iwhattodo
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,0,NI_in,NVR_ext)
if (my_rank.eq.0) then
    allocate (RHS_pm_out(neqpm+1,NXpm,NYpm,NZpm))
    RHS_pm_out=0.d0
endif
 call remesh_particles_3d(1)
if (my_rank.eq.0) then 
   allocate(QPO(1:4,NVR_ext),XPO(1:3,NVR_ext),Qflag(NVR_ext))
   allocate(QP_in(1:4,NVR_ext),XP_in(1:3,NVR_ext))
   Qflag=0
   QPO=QPR
   XPO=XPR
   NVR_turb=NVR_ext
   XMIN = XMIN_pm +(NXs_bl(1)  - 1)*DXpm
   XMAX = XMIN_pm +(NXf_bl(1)  - 1)*DXpm
   XPO(1,:)=XPR(1,:)-(XMAX-XMIN)
   QP_in=QPO;XP_in=XPO
   call writepar(0,XPO,NVR_ext)
  !NVR_ext=NYpm*NZpm
endif
!call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,0,NI_in,NVR_ext)
 call remesh_particles_3d(1)
T=0
do i=1,100
!get velocities and deformations
 T = DT_in 

if (my_rank.eq.0) then 
   allocate(UPR(3,NVR_ext),GPR(3,NVR_ext))
   UPR=0;GPR=0
endif
call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,2,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
if (my_rank.eq.0) then 
     do j= 1,NVR_ext
         XPR(1:3,j) = XPR(1:3,j)  + (UPR(1:3,j)+UINF(1:3)) * DT_in

         FACDEF = 1.
        !OMET   = sqrt ( QPR(1,j)**2 + QPR(2,j)**2 + QPR(3,j)**2 )
        !OG     = sqrt ( GPR(1,j)**2 + GPR(2,j)**2 + GPR(3,j)**2 )
        ! if (OG.ne.0.) then                                                      !RMETM
        !    if (OMET.gt.0.001)  FACDEF = OMET*MIN(RMETM,DT_in*OG/OMET)/OG/DT_in  !RMETM
        ! endif                                                                   !RMETM
 
         QPR(1:3,j) = QPR(1:3,j)  - FACDEF * GPR(1:3,j) * DT_in !minus beacuse defromation is negative
         
     enddo
     deallocate (UPR,GPR)
     call find_par_in(T,UINF(1),NVR_turb)
     call find_par_out
endif
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 if (mod(i,1).eq.0) call remesh_particles_3d(1)
!call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,5,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
!if (my_rank.eq.0) then 
!   write(*,*) maxval(GPR(:,:))
!    do j= 1,NVR_ext
!        QPR(1:3,j) = QPR(1:3,j)  -GPR(1:3,j) * DT_in
!    enddo
!endif
!get velocities and deformation
enddo

!if (my_rank.eq.0) then 
!allocate (QPDUM(7,NVR_ext),XPDUM(3,NVr_ext))
!XPDUM =XPR
!QPDUM(1:3,:)=QPR(1:3,:)
!QPDUM(4:6,:)=QPR(1:3,:)
!QPDUM(7,:)  =QPR(4,:/)
!endif
! call vpm(XPDUM,QPDUM,UPR,GPR,NVR_ext,6,1,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
!  if (mod(i,1).eq.0) call remesh_particles_3d(1)
 call MPI_FINALIZE(ierr)
End Program test_pm

 Subroutine writepar(NTIME,XPR,NVR)
    Implicit None
     integer,intent(in) :: NTIME,NVR
     double precision,intent(in):: XPR(3,NVR)
     integer ::i
    character*80 :: filout1
    write(filout1,'(i5.5,a)') NTIME,'vr.dat'
    open(10,file=filout1)
    WRITE(10,*)'VARIABLES = "X" "Y" "Z" '
    do  i=1, NVR
        write(10,'(3(F20.10,1x))') XPR(1,i), XPR(2,i),XPR(3,i)
    enddo
    call system('~/bin/preplot '&
        //filout1//' >/dev/null')
    call system('rm '//filout1)
    close(10)

 End Subroutine writepar

 Subroutine find_par_out
   use vpm_lib
   use test_mod
   use pmgrid
   use projlib
   use MPI
   Implicit None
   integer  :: i,NVR_pm,NVR_in,NVR_out,NVR_out_max,j,k
   integer         ,allocatable :: NVR_projout(:),ieq(:)
   double precision             :: XMAX,XMIN,YMAX,YMIN,ZMIN,ZMAX,EPSX,EPSY,EPSZ,XPM,YPM,ZPM
   double precision             :: xc,yc,zc,dx,dy,dz
   double precision,allocatable :: XP_out(:,:),QP_out(:,:),XP_tmp(:,:),QP_tmp(:,:),QINF(:)
   character *50                ::filout

    !an eps for arithmetic reasons
    EPSX = 0.01*DXpm  
    EPSY = 0.01*DYpm
    EPSZ = 0.01*DZpm
    XMAX = XMIN_pm + (NXf_bl(1)-(interf_iproj/2 + 0) - 1)*DXpm - EPSX
    XMIN = XMIN_pm + (NXs_bl(1)+(interf_iproj/2 + 0) - 1)*DXpm + EPSX
    YMAX = YMIN_pm + (NYf_bl(1)-(interf_iproj/2 + 0) - 1)*DYpm - EPSY
    YMIN = YMIN_pm + (NYs_bl(1)+(interf_iproj/2 + 0) - 1)*DYpm + EPSY
    ZMAX = ZMIN_pm + (NZf_bl(1)-(interf_iproj/2 + 0) - 1)*DZpm - EPSZ
    ZMIN = ZMIN_pm + (NZs_bl(1)+(interf_iproj/2 + 0) - 1)*DZpm + EPSZ
    xc = 0.5d0*(XMAX + XMIN)
    yc = 0.5d0*(YMAX + YMIN)
    zc = 0.5d0*(ZMAX + ZMIN)
    dx  = XMAX - XMIN
    dy  = YMAX - YMIN
    dz  = ZMAX - ZMIN
    NVR_out=0
     
  ! do i=1,NVR_ext
  !    if (XPR(1,i).lt.XMIN.or.XPR(1,i).gt.XMAX) then 
  !         XPR(1,i) = XPR(1,i)  - int( (XPR(1,i) - xc)/(0.5d0*dx))*dx
  !       !NVR_out = NVR_out + 1 
  !       !XP(1,i) = XPR(1,i)  -2*(XPR(1,i)-xc)
  !    endif
  !    if (XPR(2,i).lt.YMIN.or.XPR(2,i).gt.YMAX) then 
  !        XPR(2,i) = XPR(2,i)  - int( (XPR(2,i) - yc)/(0.5d0*dy))*dy
  !       !XPR(2,i) = XPR(2,i)  -2*(XPR(2,i)-yc)
  !       !NVR_out = NVR_out + 1 
  !    endif
  !    if (XPR(3,i).lt.ZMIN.or.XPR(3,i).gt.ZMAX) then 
  !        XPR(3,i) = XPR(3,i)  - int( (XPR(3,i) - zc)/(0.5d0*dz))*dz
  !       !XPR(3,i) = XPR(3,i)  -2*(XPR(3,i)-zc)
  !      !NVR_out = NVR_out + 1 
  !    endif
  ! enddo
  !return
   NVR_out_max= (2*NXpm*NYpm+2*NYpm*NZpm+2*NZpm*NXpm)*3*mrem**2
    allocate(XP_out(1:3,NVR_out_max),QP_out(1:4,NVR_out_max),NVR_projout(NVR_out_max))
   allocate(XP_tmp(1:3,NVR_ext),QP_tmp(1:4,NVR_ext))
   NVR_projout=2!interf_iproj
   NVR_out=0
   NVR_in =0
   do i=1,NVR_ext
      if(XPR(2,i).lt.YMIN.or.XPR(2,i).gt.YMAX&
         .or.XPR(3,i).lt.ZMIN.or.XPR(3,i).gt.ZMAX) then 
         NVR_out = NVR_out + 1 
         XP_out(1:3,NVR_out) = XPR(1:3,i)
         QP_out(1:4,NVR_out) = QPR(1:4,i)
      else
         NVR_in = NVR_in + 1 
         XP_tmp(1:3,NVR_in) = XPR(1:3,i)
         QP_tmp(1:4,NVR_in) = QPR(1:4,i)
      endif
   enddo
 !if (NVR_out.eq.0) RHS_pm_out=0.d0
   write(*,*) 'Particles out',NVR_out
  deallocate(XPR,QPR)
  allocate(XPR(3,NVR_in),QPR(4,NVR_in))
  XPR(1:3,1:NVR_in)=XP_tmp(1:3,1:NVR_in);QPR(1:4,1:NVR_in)=QP_tmp(1:4,1:NVR_in)
  NVR_ext = NVR_in
 !deallocate(XP_tmp,QP_tmp)
 !allocate(ieq(neqpm+1),QINF(neqpm+1))
 !QINF=0.d0
 !do i=1,neqpm+1
 !   ieq(i)=i
 !enddo
 !call project_particles_3D(RHS_pm_out,QP_out,XP_out,NVR_projout,NVR_out,neqpm+1,ieq,neqpm+1,QINF,NVR_out_max)
 !RHS_pm_out(neqpm+1,:,:,:)=DVpm
 !write(*,*) maxval(abs(QP_out(1,:))),maxval(abs(QPR(1,:)))
 !call project_vol3d(RHS_pm_out,neqpm+1,ieq,neqpm+1,1)
 !write(*,*) 'out',NVR_out,NVR_in,maxval(abs(RHS_pm_out(1,:,:,:)))
  deallocate(XP_out,QP_out)
 

 !      write(filout,'(a)') 'solout.dat'
 !      open(1,file=filout)
 !      WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "VORTX" "VORTY" "VORTZ" '
 !      WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,&
 !          ' K=',NZpm,' F=POINT'
 !      do k=1,NZpm
 !          do j=1,NYpm
 !              do i=1,NXpm
 !                  ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
 !                  ! do j=1,NYpm
 !                  !   do i=1,NXpm
 !                  XPM=XMIN_pm+(I-1)*DXpm
 !                  YPM=YMIN_pm+(J-1)*DYpm
 !                  ZPM=ZMIN_pm+(K-1)*DZpm

 !                  WRITE(1,'(6(e28.17,1x))')XPM,YPM,ZPM,-RHS_pm_out(1,I,J,K),-RHS_pm_out(2,I,J,K),&
 !                                           -RHS_pm_out(3,I,J,K)


 !              enddo
 !          enddo
 !      enddo
 !      close(1)
 !  !   ---FOR PLOTTING PURPOSES ONLY
 !  call system('~/bin/preplot '//filout//' >/dev/null')
 
!   call system('rm '//filout)
End Subroutine find_par_out

Subroutine find_par_in(T_in,U,NVR_turb)
   use vpm_lib
   use test_mod
   use pmgrid
   use projlib
   use MPI
   Implicit None
   integer ,intent(in) ::NVR_turb
   double precision,intent(in):: T_in,U
   double precision           :: XO,XMIN
   integer                    :: NVR_in
   integer  :: i,NVR_pm,NVR_in_max,j,k
   double precision,allocatable :: XP_out(:,:),QP_out(:,:),XP_tmp(:,:),QP_tmp(:,:),QINF(:)
   XMIN = XMIN_pm +(NXs_bl(1) + 2 - 1)*DXpm
   XO = XMIN
   if (minval(Qflag).eq.1) then 
       Qflag=0
       QPO=QP_in;XPO=XP_in
       write(155,*) 'Reseting inflow Particles'
   endif
  NVR_in_max= (2*NXpm*NYpm+2*NYpm*NZpm+2*NZpm*NXpm)*3*mrem**2
  allocate(XP_tmp(3,NVR_ext+NVR_in_max),QP_tmp(4,NVR_ext+NVR_in_max))
  XP_tmp(1:3,1:NVR_ext)=XPR(1:3,1:NVR_ext)
  QP_tmp(1:4,1:NVR_ext)=QPR(1:4,1:NVR_ext)
  NVR_in=0
  do i=1,NVR_turb
     XPO(1,i)=XPO(1,i) + U*T_in
     if (XPO(1,i).gt.XO.and.qflag(i).eq.0) then 
         NVR_in=NVR_in+1
         Qflag(i)=1
         XP_tmp(1:3,NVR_ext+NVR_in)=XPO(1:3,i)
         QP_tmp(1:4,NVR_ext+NVR_in)=QPO(1:4,i)
     endif
  enddo
 deallocate (XPR,QPR)
 NVR_ext=NVR_ext+NVR_in
 allocate(XPR(1:3,NVR_ext),QPR(1:4,NVR_ext))
 XPR=XP_tmp
 QPR=QP_tmp
 deallocate(XP_tmp,QP_tmp)


  write(*,*) 'Particles in',NVR_in



End Subroutine find_par_in

