
Program test_pm
use vpm_lib
use test_mod
use pmgrid
use MPI
double precision :: Vref,NI_in,DT_in,RMETM,OMET,OG,FACDEF
double precision,allocatable ::velsavex(:,:,:)
double precision,allocatable ::XPDUM(:,:),QPDUM(:,:)
integer          :: Noutput, NDimoutput
integer :: my_rank,np,ierr,i,neq,j

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 NI_in=-0.1
 DT_in=0.761905/2.  !=dx/U=8/10.5   !1000
 DT_in=0.761905     !=dx/U=8/10.5   !1000

 if (my_rank.eq.0) then 
    open(1,file='particles.bin')
    read(1,*) NVR_ext
    read(1,*) Vref
    allocate(XPR(3,NVR_ext),QPR(4,NVR_ext))
    allocate(UPR(3,NVR_ext),GPR(3,NVR_ext))
    write(*,*) 'NVR=',NVR_ext,Vref
    do i=1,NVR_ext
       read(1,*) XPR(1,i),XPR(2,i),XPR(3,i),QPR(1,i),QPR(2,i),QPR(3,i)
    enddo

 QPR(1:3,:) = -QPR(1:3,:) * Vref
 QPR(4,:) =Vref
 RMETM=0.001
endif
neq=3
 !-Iwhattodo
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,0,NI_in,NVR_ext)
 call remesh_particles_3d(1)

do i=1,1000
!get velocities and deformations
!call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,2,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 !if (mod(i,1).eq.0) call remesh_particles_3d(1)
if (my_rank.eq.0) then 
    write(*,*) maxval(UPR(:,:))
     do j= 1,NVR_ext
         XPR(1:3,j) = XPR(1:3,j)  + UPR(1:3,j) * DT_in

         FACDEF = 1.
        !OMET   = sqrt ( QPR(1,j)**2 + QPR(2,j)**2 + QPR(3,j)**2 )
        !OG     = sqrt ( GPR(1,j)**2 + GPR(2,j)**2 + GPR(3,j)**2 )
        ! if (OG.ne.0.) then                                                      !RMETM
        !    if (OMET.gt.0.001)  FACDEF = OMET*MIN(RMETM,DT_in*OG/OMET)/OG/DT_in  !RMETM
        ! endif                                                                   !RMETM
 
         QPR(1:3,j) = QPR(1:3,j)  - FACDEF * GPR(1:3,j) * DT_in !minus beacuse defromation is negative
         
     enddo
  !  call writepar(i,XPR,NVR_ext)
 GPR=0
endif
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,5,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
if (my_rank.eq.0) then 
    write(*,*) maxval(GPR(:,:))
     do j= 1,NVR_ext
         QPR(1:3,j) = QPR(1:3,j)  -GPR(1:3,j) * DT_in
     enddo
endif
  if (mod(i,1).eq.0) call remesh_particles_3d(1)
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
