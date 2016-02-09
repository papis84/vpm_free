
Program test_pm
use vpm_lib
use test_mod
use pmgrid
use MPI
double precision :: Vref,NI_in,DT_in
double precision,allocatable ::velsavex(:,:,:)
double precision,allocatable ::XPDUM(:,:),QPDUM(:,:)
integer          :: Noutput, NDimoutput
integer :: my_rank,np,ierr,i,neq,j

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 NI_in=1e-06
 DT_in=0.761905  !=dx/U=8/10.5   !1000

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
endif
neq=3
 !-Iwhattodo
!call remesh_particles_3d(1)

do i=1,1000
!get velocities and deformations
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,0,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 call vpm(XPR,QPR,UPR,GPR,NVR_ext,neq,2,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
 !if (mod(i,1).eq.0) call remesh_particles_3d(1)
if (my_rank.eq.0) then 
     do j= 1,NVR_ext
         XPR(1:3,j) = XPR(1:3,j)  + UPR(1:3,j) * DT_in
         QPR(1:3,j) = QPR(1:3,j)  + GPR(1:3,j) * DT_in
     enddo
endif
!get velocities and deformation
enddo

!if (my_rank.eq.0) then 
!allocate (QPDUM(7,NVR_ext),XPDUM(3,NVr_ext))
!XPDUM =XPR
!QPDUM(1:3,:)=QPR(1:3,:)
!QPDUM(4:6,:)=QPR(1:3,:)
!QPDUM(7,:)  =QPR(4,:)
!endif
! call vpm(XPDUM,QPDUM,UPR,GPR,NVR_ext,6,1,RHS_pm_in,velx,vely,velz,i,NI_in,NVR_ext)
!  if (mod(i,1).eq.0) call remesh_particles_3d(1)
 call MPI_FINALIZE(ierr)
End Program test_pm
