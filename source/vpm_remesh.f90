!Subroutine remesh_particles
!This Subroutine remeshes particles  on the pm grid (4 particles per cell)
!---------------------------------------------------------------------------------------------------

Subroutine remesh_particles_3d(iflag)
 use vpm_vars
 use pmeshpar
 use pmgrid
 use parvar
 use projlib
 use openmpth
 use mpi
 use test_mod
  Implicit None
 
  integer,intent(in)  :: iflag
! double precision,intent(inout):: XP_in(:,:),QP_in(:,:)
  double precision    :: X(8), Y(8),Z(8),Vol,XMIN_vr,YMIN_vr,DXvr,DYvr,ANG,dens1,dens2,Mach,DZvr,ZMIN_vr
  double precision, allocatable:: XC(:), YC(:), ZC(:)
  double precision, allocatable,target::XP_tmp(:,:),QP_tmp(:,:)
  integer             :: i ,j ,k, NXpm1, NYpm1, NZpm1, ncell, npar,ndumc,nv,inode,jnode,knode,itype
  integer             :: nxstart, nxfin, nystart, nyfin, nzstart, nzfin,ndum_rem,nnod,nc,iis,jjs,kks,iif,jjf,kkf
  integer             :: iis2,jjs2,kks2,iif2,jjf2,kkf2
  double precision    :: fx,fy,f,projection_fun
  integer,allocatable :: ieq(:)
  double precision    :: Xbound(6), Dpm(3),wmag
  double precision    :: w1,w2,r1,r2,core,radi,th,xx,yy
  double precision,allocatable :: QINF(:),rhsper(:,:,:,:)
  integer             :: my_rank,ierr,np,NN(3),NN_bl(6)

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
    XMIN_vr = XMIN_pm
    YMIN_vr = YMIN_pm
    ZMIN_vr = ZMIN_pm
    DXvr    = DXpm
    DYvr    = DYpm
    DZvr    = DZpm
    Ndumc   = 1
    npar_cell = 1
    
    Dpm(1)= DXpm; Dpm(2)= DYpm; Dpm(3)= DZpm
    Xbound(1) = XMIN_pm;Xbound(2) = YMIN_pm ;Xbound(3)= ZMIN_pm 
    Xbound(4) = XMAX_pm;Xbound(5) = YMAX_pm ;Xbound(6)= ZMAX_pm 
!->PM grid is orthogonal (Volume of PM cell

!--------------------------------------------------------------------------------!
!-->The loops starts from 2 because we need cells that DO NOT contain particles  !
!-->Total Number of Cells                                                        !
!--------------------------------------------------------------------------------!
 NN(1) = NXpm; NN(2)= NYpm; NN(3)= NZpm
 NN_bl(1) = NXs_bl(1);NN_bl(2)= NYs_bl(1);NN_bl(3)=NZs_bl(1) 
 NN_bl(4) = NXf_bl(1);NN_bl(5)= NYf_bl(1);NN_bl(6)=NZf_bl(1) 

 NN = NN *  mrem
 NN_bl = NN_bl *mrem
 Dpm(1) = (Xbound(4)-Xbound(1))/(NN(1)-1)
 Dpm(2) = (Xbound(5)-Xbound(2))/(NN(2)-1)
 Dpm(3) = (Xbound(6)-Xbound(3))/(NN(3)-1)

DVpm = Dpm(1)*Dpm(2)*Dpm(3)
if (my_rank.eq.0) then 
NVR=NVR_ext
 XP=>XPR
 QP=>QPR
endif
if(iflag.eq.1) then
       if (allocated(RHS_pm)) then 
           deallocate(RHS_pm)
           allocate(RHS_pm(neqpm+1,NN(1),NN(2),NN(3)))
       else
           allocate(RHS_pm(neqpm+1,NN(1),NN(2),NN(3)))
       endif
       call MPI_BCAST(NVR,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       NVR_p=NVR/np
       if (my_rank.eq.0) NVR_p = NVR_p + mod(NVR,np)
       allocate (XP_scatt(3,NVR_p),QP_scatt(neqpm+1,NVR_p),NVR_projscatt(NVR_p))
       NVR_projscatt=interf_iproj
       call particles_scat
       call  projlibinit(Xbound,Dpm,NN,NN_bl,EPSVOL,IDVPM,ND)
       allocate(ieq(neqpm+1),QINF(neqpm+1))
       QINF=0.d0
       do i=1,neqpm+1
          ieq(i)=i
       enddo
       call project_particles_3D(RHS_pm,QP_scatt,XP_scatt,NVR_projscatt,NVR_p,neqpm+1,ieq,neqpm+1,QINF,NVR_p)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call proj_gath(NN)
       deallocate(XP_scatt,QP_scatt,NVR_projscatt)

       if (my_rank.eq.0) then 
           call omp_set_num_threads(OMPTHREADS)
           RHS_pm(neqpm+1,:,:,:)=DVpm
           call project_vol3d(RHS_pm,neqpm+1,ieq,neqpm+1,IDVPM)
           call omp_set_num_threads(1)
          !call hill_assign(NN,NN_bl,Xbound,Dpm,RHS_pm,neqpm+1)

       endif

       deallocate(ieq,QINF)
endif

if (my_rank.eq.0) then 
    deallocate(XPR,QPR)
    ncell =ncell_rem
    ndum_rem=2
    nnod  =1 !if ncell gt 1 particles IN cell else in nodes
    if (ncell.eq.1) nnod=0
    allocate (XC(ncell), YC(ncell), ZC(ncell))
    if (ncell.eq.1) then 
        nxfin   = NN_bl(4)  -interf_iproj/2 
        nxstart = NN_bl(1)  +interf_iproj/2 
        nyfin   = NN_bl(5)  -interf_iproj/2 
        nystart = NN_bl(2)  +interf_iproj/2  
        nzfin   = NN_bl(6)  -interf_iproj/2 
        nzstart = NN_bl(3)  +interf_iproj/2   
    else
        nxfin   = NN_bl(4) -1 
        nxstart = NN_bl(1) 
        nyfin   = NN_bl(5) -1 
        nystart = NN_bl(2) 
        nzfin   = NN_bl(6) -1 
        nzstart = NN_bl(3) 
    endif

    NXpm1   = nxfin - nxstart + 1 
    NYpm1   = nyfin - nystart + 1
    NZpm1   = nzfin - nzstart + 1
    NVR       =  NXpm1*NYpm1*NZpm1*ncell
    allocate(XP_tmp(3,NVR),QP_tmp(neqpm+1,NVR))
    XP=>XP_tmp
    QP=>QP_tmp
    XP=0
    QP=0
    npar = 0
    V_ref = 1.d0/float(ncell)*DVpm
!!$omp parallel private(i,j,k,npar,X,Y,Z) num_threads(OMPTHREADS)
!!$omp do
   do k = nzstart, nzfin
    do j = nystart, nyfin
     do i = nxstart, nxfin
! !-> Get PM cell nodes (orthogonal structured grid    
         X(1) = XMIN_pm + Dpm(1) * (i-1)
         X(2) = XMIN_pm + Dpm(1) * (i)
         X(3) = XMIN_pm + Dpm(1) * (i)
         X(4) = XMIN_pm + Dpm(1) * (i-1)
         X(5) = XMIN_pm + Dpm(1) * (i-1)
         X(6) = XMIN_pm + Dpm(1) * (i)
         X(7) = XMIN_pm + Dpm(1) * (i)
         X(8) = XMIN_pm + Dpm(1) * (i-1)
              
         Y(1) = YMIN_pm + Dpm(2) * (j-1)
         Y(2) = YMIN_pm + Dpm(2) * (j-1)
         Y(3) = YMIN_pm + Dpm(2) * (j)
         Y(4) = YMIN_pm + Dpm(2) * (j)
         Y(5) = YMIN_pm + Dpm(2) * (j-1)
         Y(6) = YMIN_pm + Dpm(2) * (j-1)
         Y(7) = YMIN_pm + Dpm(2) * (j)
         Y(8) = YMIN_pm + Dpm(2) * (j)
              
         Z(1) = ZMIN_pm + Dpm(3) * (k-1)
         Z(2) = ZMIN_pm + Dpm(3) * (k-1)
         Z(3) = ZMIN_pm + Dpm(3) * (k-1)
         Z(4) = ZMIN_pm + Dpm(3) * (k-1)
         Z(5) = ZMIN_pm + Dpm(3) * (k)
         Z(6) = ZMIN_pm + Dpm(3) * (k)
         Z(7) = ZMIN_pm + Dpm(3) * (k)
         Z(8) = ZMIN_pm + Dpm(3) * (k)

       !npar = ((k-nzstart)*NXpm1*NYpm1 +(j-nystart)*NXpm1 + i-nxstart)*ncell

        if (ncell.gt.1) then 
           call cell3d_interp_euler(X, XC,ncell,2)
           call cell3d_interp_euler(Y, YC,ncell,2)
           call cell3d_interp_euler(Z, ZC,ncell,2)
           do nc = 1, ncell
              npar = npar+1
              XP(1,npar) = XC(nc)
              XP(2,npar) = YC(nc)
              XP(3,npar) = ZC(nc)
              QP(neqpm+1,npar) =1.d0/float(ncell)*DVpm
           enddo

        else
            wmag = sqrt(RHS_pm(1,i,j,k)**2 +RHS_pm(2,i,j,k)**2+RHS_pm(3,i,j,k)**2)
            if (wmag.lt.1e-09) cycle
            npar = npar + 1
            XP(1,npar)= X(1)
            XP(2,npar)= Y(1)
            XP(3,npar)= Z(1)
                 
            QP(1:neqpm,npar)= RHS_pm(1:neqpm,i,j,k)* DVpm
            QP(neqpm+1,npar)= DVpm
        endif
     enddo
    enddo
   enddo
 ! !$omp enddo
 ! !$omp endparallel
    NVR=npar
    NVR_ext=NVR
    NVR_size=NVR
    allocate(XPR(3,NVR),QPR(neqpm+1,NVR))
    XPR(1:3,1:NVR)=XP_tmp(1:3,1:NVR)
    QPR(1:neqpm+1,1:NVR)=QP_tmp(1:neqpm+1,1:NVR)
    deallocate(XP_tmp,QP_tmp)
    XP=>XPR
    QP=>QPR
   if (ncell.gt.1) call back_to_particles_3D_rem(RHS_pm,XP,QP,Xbound,Dpm,NN,NVR,4)
   if (iflag.eq.0) deallocate(RHS_pm)

write(*,*) 'remesh complete',NVR,npar,DVpm,maxval(QPR(neqpm,:))
endif
!call back_to_particles_2D(4)

 ! open(1,file='vr.dat')
 !  WRITE(1,*)'VARIABLES = "X" "Y" "Z"'
 !  do  i=1, NVR
 !       write(1,'(2(F20.10,1x))') XP(i,1), XP(i,2),XP(i,3)
 !  enddo
!!-----FOR PLOTTING PURPOSES ONLY
 !     call system('preplot vr.dat>/dev/null')
 !     call system('rm vr.dat')
!!----FOR PLOTTING PURPOSES ONLY
 ! close(1)
deallocate(RHS_pm)
End Subroutine remesh_particles_3d


!---------------------------------------------------------------------------!
!-> Subroutine back_to_particles                                            !
!   This subroutine interpolates PM grid values back to particles at the    !
!   positions they ARE.Convections takes place afterwards.                  !
!   Input :                                                                 !
!          itype (1,2) defines what value to interpolate to the particles   !
!---------------------------------------------------------------------------!
Subroutine back_to_particles_3D_rem(RHS_pm,XP,QP,Xbound,Dpm,NN,NVR,iproj)

    use openmpth
    Implicit None
    integer, intent(in)             :: NN(3),NVR,iproj
    double precision, intent(in)    :: RHS_pm(4,NN(1),NN(2),NN(3)) 
    double precision, intent(inout) :: QP(4,NVR),XP(3,NVR)

    double precision, intent(in)    :: Xbound(6),Dpm(3)

    double precision :: fx, fy, fz, f, x, y, z, projection_fun
    integer          :: inode, jnode, knode, i, j , k, nv,ips,ipf,ivortx,ivorty,DVpm



    if(iproj.eq.2) then
        ips =  0
        ipf =  1
    else if (iproj.eq.3) then
        ips =  1
        ipf =  2
    else if (iproj.eq.4) then
        ips =  1
        ipf =  2
    endif
        QP(:,1:3)=0
        do  nv = 1, NVR
            !-->Find the cell/node  the  particle belongs for X and Y and Z direction.
            inode = int((XP(1,nv) - XBound(1)) / Dpm(1)) + 1
            jnode = int((XP(2,nv) - XBound(2)) / Dpm(2)) + 1
            knode = int((XP(3,nv) - XBound(3)) / Dpm(3)) + 1

            !--We search the 4 nodes close to the particles
            do k = knode - ips, knode + ipf
                do j = jnode - ips, jnode + ipf
                    do i = inode - ips, inode + ipf
                        x  = (XP(1,nv)- XBound(1) - (i-1) * Dpm(1)) / Dpm(1)
                        fx = projection_fun(iproj,x)

                        y  = (XP(2,nv)- XBound(2) - (j-1) * Dpm(2)) / Dpm(2)
                        fy = projection_fun(iproj,y)

                        z  = (XP(3,nv)- XBound(3) - (k-1) * Dpm(3)) / Dpm(3)
                        fz = projection_fun(iproj,z)

                        f = fx * fy * fz
                        QP(1:3,nv)    = QP(1:3,nv)    +  f * RHS_pm(1:3,i,j,k)

                    enddo
                enddo
            enddo
            QP(1:3,nv) = QP(1:3,nv) * QP(4,nv)
        enddo

End Subroutine back_to_particles_3D_rem

!--------------------------------------------------------------------------------
!>@function
! SUBROUTINE    cell3d_interp_euler    
!>
!>@author Papis
!>
!>@brief       
!>Subroutine cell3d_interp_euler creates 4 or more particles per cell using ksi ita
!>coordinates
!REVISION HISTORY
!> 17/7/2013 - Initial Version
!> TODO_dd
!>
!>@param [in]  F is the value at the global coordinates
!>@param [out] FC is the value at global coordinates of the interpolated value
!--------------------------------------------------------------------------------

Subroutine cell3d_interp_euler(F, FC, N, M)
  Implicit None
  
  integer         , intent(in)  :: N, M
  double precision, intent(in)  :: F(8)
  double precision, intent(out) :: FC(N)

  double precision              :: KSIC(8) , HTAC(8),ZETAC(8),KSI(N), HTA(N), ZETA(N) 

  integer                       :: i,j
    FC=0
!-->Define KSI,HTA corners
    KSIC(1) = -1.d0
    KSIC(2) =  1.d0
    KSIC(3) =  1.d0
    KSIC(4) = -1.d0
    KSIC(5) = -1.d0
    KSIC(6) =  1.d0
    KSIC(7) =  1.d0
    KSIC(8) = -1.d0

    HTAC(1) = -1.d0
    HTAC(2) = -1.d0
    HTAC(3) =  1.d0
    HTAC(4) =  1.d0
    HTAC(5) = -1.d0
    HTAC(6) = -1.d0
    HTAC(7) =  1.d0
    HTAC(8) =  1.d0

    ZETAC(1) = -1.d0
    ZETAC(2) = -1.d0
    ZETAC(3) = -1.d0
    ZETAC(4) = -1.d0
    ZETAC(5) =  1.d0
    ZETAC(6) =  1.d0
    ZETAC(7) =  1.d0
    ZETAC(8) =  1.d0
    
    call get_ksi_ita_pos_3d(N,M,KSIC,HTAC,ZETAC,KSI,HTA,ZETA)

    do i = 1, 8 !cell nodes
       do j= 1, N 
          FC(j) = FC(j) + F(i) * (1.d0 + KSI(j)*KSIC(i)) * (1.d0 + HTA(j)*HTAC(i)) *(1.d0+ ZETA(j)*ZETAC(i))  
       enddo
    enddo
      
    FC    = 0.125d0 * FC !1/8
 End Subroutine cell3d_interp_euler

!--------------------------------------------------------------------------------
!>@function
! SUBROUTINE    get_ksi_ita_pos        
!>
!>@author Papis
!>
!>@brief       
!>Subroutine get_ksi_ita_pos  depending on the defined number of particles(must be perfect square)
!REVISION HISTORY
!> 22/7/2013 - Initial Version
!> TODO_dd
!>
!>@param [in]  N is the number of particles                             
!>@param [in]  KSIC(4),HTAC(4) is the corner coordinates in the KSI,HTA                         
!>@param [out] KSI(2*N),HTA(2*N) local position                                   
!--------------------------------------------------------------------------------
 Subroutine get_ksi_ita_pos_3d(N,M,KSIC,HTAC,ZETAC,KSI,HTA,ZETA)
 Implicit None

 integer, intent(in)           :: N,M
 double precision, intent(in)  :: KSIC(8), HTAC(8),ZETAC(8)
 double precision, intent(out) :: KSI(N), HTA(N),ZETA(N)
 double precision              :: DKSI, DHTA,DZETA
 integer                       :: i,j,k,nod, L
  
!--> find position minus ksi minus ita quadrant and then by symmerty * 4

 KSI=0.d0
 HTA=0.d0
 
 DKSI = dabs(2.d0/float(M))
 DHTA = dabs(2.d0/float(M))
 DZETA = dabs(2.d0/float(M))
 do k = 1, M 
    do i = 1, M
       do j = 1, M
          nod = (k-1)* M * M + (j-1)*M + i 
          KSI(nod)  = KSIC(1)  + (i-1./2.) * DKSI
          HTA(nod)  = HTAC(1)  + (j-1./2.) * DHTA
          ZETA(nod) = ZETAC(1) + (k-1./2.) * DZETA
       enddo
    enddo
 enddo
 
 End Subroutine get_ksi_ita_pos_3d
