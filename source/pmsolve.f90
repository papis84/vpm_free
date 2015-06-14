!-----------------------------------------------------------------------!
!-> Subroutine solve_phiz                                                !
!   This Subroutines calls the fft library to solve for Phi poisson     !
!   in all the points of Particle mesh.Dirichlet Boundary Cond. are used!
!-----------------------------------------------------------------------!
Subroutine solve_eq(NXs,NXf,NYs,NYf,neq)
    use MKL_DFTI
    Implicit None
    integer,intent(in)   :: NXs,NXf,NYs,NYf,neq
    Integer              :: i,j,NWORK,INFO,NX,NY,nbj
    integer              :: ipar(128),stat
    integer              :: NN,nod
    double precision, allocatable ::dpar(:)
    double precision     :: XPM,YPM,XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,pertrb
    double precision,allocatable::SOL_pm2(:,:),WORK(:)
    double precision,allocatable::f(:),bd_ax(:),bd_bx(:),bd_ay(:),bd_by(:)
    type(DFTI_DESCRIPTOR)    :: xhandle

    !--> Assignment of Boundary Values
    ipar=0
    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    NX  = NXf - NXs + 1
    NY  = NYf - NYs + 1
    !-->Set Right hand Side term (except for boundary conditions)

    NN= NX*NY
    allocate(f(NN))
    do j=1,NY
        do i= 1,NX
            nod = (j-1)*NX + i
            f(nod) = -RHS_pm(neq,NXs+i-1,NYs+j-1, 1)
        enddo
    enddo

    !-->Set Boundary Conditions
    !---> XMIN,XMAX

    !    Psiz_pm2(1:5,:)  = 0.d0
    !    Psiz_pm2(NX:NX-5,:) = 0.d0
    !    Psiz_pm2(:,1:5)  = 0.d0
    !    Psiz_pm2(:,NY:NY-5) = 0.d0

    NN=NY
    allocate(bd_ax(NN),bd_bx(NN))
    do j =1,NY
        bd_ax(j)      =  SOL_pm(neq,NXs,j + NYs -1, 1)
        bd_bx(j)      =  SOL_pm(neq,NXf,j + NYs -1, 1)
    enddo
    !---> YMIN,YMAX
    NN=NX
    allocate(bd_ay(NN),bd_by(NN))

    do i=1,NX
        bd_ay(i)  = SOL_pm(neq,i + NXs-1,NYs , 1)
        bd_by(i)  = SOL_pm(neq,i + NXs-1,NYf , 1)
    enddo

    !--SOLVE PHI ON PMESH
    NWORK = 4 * (NY-1) + (13 + int(log(2.*((NY-1)+1))))*(NX+1)
    allocate (WORK(NWORK))
    allocate (dpar(int(5*(NX-1)/2) + 7))
    call d_init_Helmholtz_2D(XminCalc,XmaxCalc,YminCalc,YmaxCalc,NX-1,NY-1,'DDDD',0.d0,ipar,dpar,stat)
    call d_commit_Helmholtz_2D(f,bd_ax,bd_bx,bd_ay,bd_by,xhandle,ipar,dpar,stat)
    call d_Helmholtz_2D(f,bd_ax,bd_bx,bd_ay,bd_by,xhandle,ipar,dpar,stat)
    call free_Helmholtz_2D(xhandle,ipar,stat)

    do j=1,NY
        do i=1,NX
            nod = (j-1)*NX + i
            SOL_pm(neq,NXs+i-1,NYs+j-1, 1)=f(nod)
        enddo
    enddo

    !    OPEN(UNIT=1,FILE=filename2//char(48+nb))
    !  ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
    !    DO J=NYs,NYf
    !      DO I=NXs,NXf
    !        XPM=XMIN_pm+(I-1)*DXpm
    !        YPM=YMIN_pm+(J-1)*DYpm
    !        WRITE(1,'(5(F20.16,1x))')XPM,YPM,SOL_pm(I,J,nb,neq),RHS_pm(I,J,nb,neq),Cont_pm(I,j)
    !      ENDDO
    !    ENDDO
    !    CLOSE(1)

    !    OPEN(UNIT=1,FILE=filename4//char(48+nb))
    !  ! WRITE(1,*)'ZONE I=',NXpm,' J=',NYpm,' F=POINT'
    !    DO J=NYs,NYf
    !      DO I=NXs,NXf
    !        XPM=XMIN_pm+(I-1)*DXpm
    !        YPM=YMIN_pm+(J-1)*DYpm
    !        WRITE(1,'(3(F20.16,1x))')XPM,YPM,abs(SOL_pm(I,J,nb,neq)-Cont_pm(I,j))/abs(Cont_pm(I,J))
    !      ENDDO
    !        WRITE(1,*)
    !    ENDDO
    !    DO I=NXs,NXf
    !     DO J=NYs,NYf
    !        XPM=XMIN_pm+(I-1)*DXpm
    !        YPM=YMIN_pm+(J-1)*DYpm
    !        WRITE(1,'(3(F20.16,1x))')XPM,YPM,abs(SOL_pm(I,J,nb,neq)-Cont_pm(I,j))/abs(Cont_pm(I,J))
    !      ENDDO
    !        WRITE(1,*)
    !    ENDDO
    !    CLOSE(1)
End Subroutine solve_eq!_i

Subroutine solve_eq_0(NXs,NXf,NYs,NYf,neq)
    use MKL_DFTI
    Implicit None
    integer,intent(in)   :: NXs,NXf,NYs,NYf,neq
    Integer              :: i,j,NWORK,INFO,NX,NY,nbj
    integer              :: ipar(128),stat
    integer              :: NN,nod
    double precision, allocatable ::dpar(:)
    double precision     :: XPM,YPM,XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,pertrb
    double precision,allocatable::SOL_pm2(:,:),WORK(:)
    double precision,allocatable::f(:),bd_ax(:),bd_bx(:),bd_ay(:),bd_by(:)
    type(DFTI_DESCRIPTOR)    :: xhandle

    !--> Assignment of Boundary Values
    ipar=0
    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    NX  = NXf - NXs + 1
    NY  = NYf - NYs + 1
    !-->Set Right hand Side term (except for boundary conditions)

    NN= NX*NY
    allocate(f(NN))
    do j=1,NY
        do i= 1,NX
            nod = (j-1)*NX + i
            f(nod) = -RHS_pm(neq,NXs+i-1,NYs+j-1, 1)
        enddo
    enddo

    !-->Set Boundary Conditions
    !---> XMIN,XMAX

    !    Psiz_pm2(1:5,:)  = 0.d0
    !    Psiz_pm2(NX:NX-5,:) = 0.d0
    !    Psiz_pm2(:,1:5)  = 0.d0
    !    Psiz_pm2(:,NY:NY-5) = 0.d0

    NN=NY
    allocate(bd_ax(NN),bd_bx(NN))
    bd_ax=0.d0;bd_bx=0.d0;!  SOL_pm(NXs,j + NYs -1,nb,neq)
    !---> YMIN,YMAX
    NN=NX
    allocate(bd_ay(NN),bd_by(NN))
    bd_ay =0.d0;bd_by =0.d0! SOL_pm(i + NXs-1,NYf ,nb,neq)

    !--SOLVE PHI ON PMESH
    allocate (dpar(int(5*(NX-1)/2) + 7 ))
    call d_init_Helmholtz_2D(XminCalc,XmaxCalc,YminCalc,YmaxCalc,NX-1,NY-1,'DDDD',0.d0,ipar,dpar,stat)
    call d_commit_Helmholtz_2D(f,bd_ax,bd_bx,bd_ay,bd_by,xhandle,ipar,dpar,stat)
    call d_Helmholtz_2D(f,bd_ax,bd_bx,bd_ay,bd_by,xhandle,ipar,dpar,stat)
    call free_Helmholtz_2D(xhandle,ipar,stat)

    !call hwscrt(XminCalc,XmaxCalc,NX-1,1,0,0,YminCalc,YmaxCalc,NY-1,1,0,0,  &
        !                   0,SOL_pm2,NX,pertrb,INFO,WORK)


    !  SOL_pm(NXs:NXf,NYs:NYf,nb,neq)=SOL_pm2(1:NX,1:NY)
    do j=1,NY
        do i=1,NX
            nod = (j-1)*NX + i
            SOL_0_pm(neq,NXs+i-1,NYs+j-1, 1)=f(nod)
        enddo
    enddo

End Subroutine solve_eq_0
!-----------------------------------------------------------------------!
!-> Subroutine solve_phiz                                                !
!   This Subroutines calls the fft library to solve for Phi poisson     !
!   in all the points of Particle mesh.Dirichlet Boundary Cond. are used!
!-----------------------------------------------------------------------!
Subroutine solve_eq_3d_i(NXs,NXf,NYs,NYf,NZs,NZf,neq)
    Implicit None
    integer,intent(in) :: NXs,NXf,NYs,NYf,NZs,NZf,neq
    Integer            :: i,j,k,NWORK,INFO,NX,NY,NZ,nbj,LPEROD,MPEROD,NPEROD,IERROR
    double precision   :: XPM,YPM,XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,ZminCalc,ZmaxCalc,CX,CY,CZ
    double precision,allocatable::SOL_pm2(:,:,:),Apois(:),Bpois(:),Cpois(:)


    !--> Assignment of Boundary Values

    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    ZminCalc = ZMIN_pm +  (NZs-1) * DZpm
    ZmaxCalc = ZMIN_pm +  (NZf-1) * DZpm

    NX  = NXf - NXs + 1 -2
    NY  = NYf - NYs + 1 -2
    NZ  = NZf - NZs + 1 -2
    allocate(SOL_pm2(NX,NY,NZ))
    allocate(Apois(NZ),Bpois(NZ),Cpois(NZ))

    CX = 1.d0/DXpm**2
    CY = 1.d0/DYpm**2
    CZ = 1.d0/DZpm**2
    Apois(1:NZ) = CZ
    Bpois(1:NZ) = -2.d0 * CZ
    Cpois(1:NZ) = CZ
    Apois(1)=0.d0
    Cpois(NZ)=0.d0
    !-->Set Right hand Side term (except for boundary conditions)
    SOL_pm2(1:NX,1:NY,1:NZ)=RHS_pm(neq,NXs+1:NXf-1,NYs+1:NYf-1,NZs+1:NZf-1)
    !-->Set Boundary Conditions
    !---> XMIN,XMAX

    !    Psiz_pm2(1:5,:)  = 0.d0
    !    Psiz_pm2(NX:NX-5,:) = 0.d0
    !    Psiz_pm2(:,1:5)  = 0.d0
    !    Psiz_pm2(:,NY:NY-5) = 0.d0
    do k =1,NZ
        do j =1,NY
            SOL_pm2(1, j,k) =SOL_pm2(1,j,k)  -CX * SOL_pm(NXs,j + NYs -1+1,k + NZs -1+1,neq)
            SOL_pm2(NX,j,k) =SOL_pm2(NX,j,k) -CX * SOL_pm(NXf,j + NYs -1+1,k + NZs -1+1,neq)
        enddo
    enddo
    !---> YMIN,YMAX
    do k =1,NZ
        do i =1,NX
            SOL_pm2(i,1,k)  = SOL_pm2(i,1,k) - CY * SOL_pm(i + NXs-1+1,NYs ,k + NZs -1+1,neq)
            SOL_pm2(i,NY,k) = SOL_pm2(i,NY,k)- CY * SOL_pm(i + NXs-1+1,NYf ,k + NZs -1+1,neq)
        enddo
    enddo

    do j =1,NY
        do i =1,NX
            SOL_pm2(i,j,1)  = SOL_pm2(i,j,1)  - CZ * SOL_pm(i + NXs-1+1,j + NYs - 1+1 ,NZs,neq)
            SOL_pm2(i,j,NZ) = SOL_pm2(i,j,NZ) - CZ * SOL_pm(i + NXs-1+1,j + NYs - 1+1 ,NZf,neq)
        enddo
    enddo



    call POIS3D(1, NX, CX, 1, NY, CY, 1, NZ, Apois, Bpois, Cpois,&
        NX,NY, Sol_pm2, IERROR)
    SOL_pm(neq,NXs+1:NXf-1,NYs+1:NYf-1,NZs+1:NZf-1)=SOL_pm2(1:NX,1:NY,1:NZ)
    if(IERROR.ne.0) write(*,*)'WRONG SOLUTION',ierror
End Subroutine solve_eq_3d_i



Subroutine solve_eq_0_3d_i(NXs,NXf,NYs,NYf,NZs,NZf,neq)
    Implicit None
    integer,intent(in) :: NXs,NXf,NYs,NYf,NZs,NZf,neq
    Integer            :: i,j,k,NWORK,INFO,NX,NY,NZ,nbj,LPEROD,MPEROD,NPEROD,IERROR
    double precision   :: XPM,YPM,XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,ZminCalc,ZmaxCalc,CX,CY,CZ
    double precision,allocatable::SOL_pm2(:,:,:),Apois(:),Bpois(:),Cpois(:)


    !--> Assignment of Boundary Values

    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    ZminCalc = ZMIN_pm +  (NZs-1) * DZpm
    ZmaxCalc = ZMIN_pm +  (NZf-1) * DZpm

    NX  = NXf - NXs + 1 -2
    NY  = NYf - NYs + 1 -2
    NZ  = NZf - NZs + 1 -2
    allocate(SOL_pm2(NX,NY,NZ))
    allocate(Apois(NZ),Bpois(NZ),Cpois(NZ))

    CX = 1.d0/DXpm**2
    CY = 1.d0/DYpm**2
    CZ = 1.d0/DZpm**2
    Apois(1:NZ) = CZ
    Bpois(1:NZ) = -2.d0 * CZ
    Cpois(1:NZ) = CZ
    Apois(1)=0.d0
    Cpois(NZ)=0.d0
    !-->Set Right hand Side term (except for boundary conditions)

    SOL_pm2(1:NX,1:NY,1:NZ)=RHS_pm(neq,NXs+1:NXf-1,NYs+1:NYf-1,NZs+1:NZf-1)
    !-->Set Boundary Conditions
    !---> XMIN,XMAX

    !    Psiz_pm2(1:5,:)  = 0.d0
    !    Psiz_pm2(NX:NX-5,:) = 0.d0
    !    Psiz_pm2(:,1:5)  = 0.d0
    !    Psiz_pm2(:,NY:NY-5) = 0.d0



    call POIS3D(1, NX, CX, 1, NY, CY, 1, NZ, Apois, Bpois, Cpois,&
        NX,NY, Sol_pm2, IERROR)

    SOL_0_pm(neq,NXs+1:NXf-1,NYs+1:NYf-1,NZs+1:NZf-1)=SOL_pm2(1:NX,1:NY,1:NZ)


    if(IERROR.ne.0) write(*,*)'WRONG SOLUTION',ierror
End Subroutine solve_eq_0_3d_i

Subroutine solve_eq_3d(NXs,NXf,NYs,NYf,NZs,NZf,neq)
    use MKL_DFTI
    Implicit None
    integer,intent(in) :: NXs,NXf,NYs,NYf,NZs,NZf,neq
    Integer            :: i,j,k,NWORK,INFO,NX,NY,NZ,nbj,LPEROD,MPEROD,NPEROD,IERROR
    double precision   :: XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,ZminCalc,ZmaxCalc
    integer              :: ipar(128),stat
    integer              :: NN,nod
    double precision, allocatable ::dpar(:)
    double precision,allocatable::f(:),bd_ax(:),bd_bx(:),bd_ay(:),bd_by(:),bd_az(:),bd_bz(:)
    type(DFTI_DESCRIPTOR)    :: xhandle,yhandle

   
    !--> Assignment of Boundary Values
    ipar=0
    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    ZminCalc = ZMIN_pm +  (NZs-1) * DZpm
    ZmaxCalc = ZMIN_pm +  (NZf-1) * DZpm

    NX  = NXf - NXs + 1!-2
    NY  = NYf - NYs + 1!-2
    NZ  = NZf - NZs + 1!-2
    NN=NX*NY*NZ
    allocate(f(NN))
    do k=1,NZ
        do j=1,NY
            do i= 1,NX
                nod = (k-1)*NX*NY + (j-1)*NX + i
                f(nod) = -RHS_pm(neq,NXs+i-1,NYs+j-1,NZs+k-1)
            enddo
        enddo
    enddo
    NN=NY*NZ
    allocate(bd_ax(NN),bd_bx(NN))
    do k=1,NZ
        do j =1,NY
            nod = (k-1)*NY+j
            bd_ax(nod)      =  SOL_pm(neq,NXs,j + NYs -1,k+NZs-1)
            bd_bx(nod)      =  SOL_pm(neq,NXf,j + NYs -1,k+NZs-1)
        enddo
    enddo
    !---> YMIN,YMAX
    NN=NX*NZ
    allocate(bd_ay(NN),bd_by(NN))
    do k=1,NZ
        do i=1,NX
            nod = (k-1)*NX+i
            bd_ay(nod)  = SOL_pm(neq,i + NXs-1,NYs ,k+NZs-1)
            bd_by(nod)  = SOL_pm(neq,i + NXs-1,NYf ,k+NZs-1)
        enddo
    enddo

    NN=NX*NY
    allocate(bd_az(NN),bd_bz(NN))
    do j=1,NY
        do i=1,NX
            nod = (j-1)*NX+i
            bd_az(nod)  = SOL_pm(neq,i + NXs-1,j + NYs -1,NZs)
            bd_bz(nod)  = SOL_pm(neq,i + NXs-1,j + NYs -1,NZf)
        enddo
    enddo
    allocate (dpar(int(5*(NX-1+NY-1)/2) + 9 ))
    call d_init_Helmholtz_3D(XminCalc,XmaxCalc,YminCalc,YmaxCalc,ZminCalc,ZmaxCalc,NX-1,NY-1,NZ-1,'DDDDDD',0.d0,ipar,dpar,stat)
    call d_commit_Helmholtz_3D(f,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,xhandle,yhandle,ipar,dpar,stat)
    call d_Helmholtz_3D(f,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,xhandle,yhandle,ipar,dpar,stat)
    call free_Helmholtz_3D(xhandle,yhandle,ipar,stat)

    do k=1,NZ
        do j=1,NY
            do i=1,NX
                nod = (k-1)*NX*NY + (j-1)*NX + i
                SOL_pm(neq,NXs+i-1,NYs+j-1,NZs+k-1)=f(nod)
            enddo
        enddo
    enddo

End Subroutine solve_eq_3d

Subroutine solve_eq_0_3d(NXs,NXf,NYs,NYf,NZs,NZf,neq)
    use MKL_DFTI
    Implicit None
    integer,intent(in) :: NXs,NXf,NYs,NYf,NZs,NZf,neq
    Integer            :: i,j,k,NWORK,INFO,NX,NY,NZ,nbj,LPEROD,MPEROD,NPEROD,IERROR
    double precision   :: XMinCalc,XmaxCalc,YMinCalc,YmaxCalc,ZminCalc,ZmaxCalc
    integer              :: ipar(128),stat
    integer              :: NN,nod
    double precision, allocatable ::dpar(:)
    double precision,allocatable::f(:),bd_ax(:),bd_bx(:),bd_ay(:),bd_by(:),bd_az(:),bd_bz(:)
    type(DFTI_DESCRIPTOR)    :: xhandle,yhandle


    !--> Assignment of Boundary Values
    ipar=0
    XminCalc = XMIN_pm +  (NXs-1) * DXpm
    XmaxCalc = XMIN_pm +  (NXf-1) * DXpm

    YminCalc = YMIN_pm +  (NYs-1) * DYpm
    YmaxCalc = YMIN_pm +  (NYf-1) * DYpm

    ZminCalc = ZMIN_pm +  (NZs-1) * DZpm
    ZmaxCalc = ZMIN_pm +  (NZf-1) * DZpm

    NX  = NXf - NXs + 1!-2
    NY  = NYf - NYs + 1!-2
    NZ  = NZf - NZs + 1!-2

    NN=NX*NY*NZ
    allocate(f(NN))

    do k=1,NZ
        do j=1,NY
            do i= 1,NX
                nod = (k-1)*NX*NY + (j-1)*NX + i
                f(nod) = -RHS_pm(neq,NXs+i-1,NYs+j-1,NZs+k-1)
            enddo
        enddo
    enddo

    NN=NY*NZ
    allocate(bd_ax(NN),bd_bx(NN))
    bd_ax=0;bd_bx=0;
    !---> YMIN,YMAX
    NN=NX*NZ
    allocate(bd_ay(NN),bd_by(NN))
    bd_ay=0;bd_by=0;
    NN=NX*NY
    allocate(bd_az(NN),bd_bz(NN))
    bd_az=0;bd_bz=0;
    allocate (dpar(int(5*(NX-1+NY-1)/2) + 9))
    call d_init_Helmholtz_3D(XminCalc,XmaxCalc,YminCalc,YmaxCalc,ZminCalc,ZmaxCalc,NX-1,NY-1,NZ-1,'DDDDDD',0.d0,ipar,dpar,stat)
    call d_commit_Helmholtz_3D(f,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,xhandle,yhandle,ipar,dpar,stat)
    call d_Helmholtz_3D(f,bd_ax,bd_bx,bd_ay,bd_by,bd_az,bd_bz,xhandle,yhandle,ipar,dpar,stat)
    call free_Helmholtz_3D(xhandle,yhandle,ipar,stat)

    do k=1,NZ
        do j=1,NY
            do i=1,NX
                nod = (k-1)*NX*NY + (j-1)*NX + i
                SOL_0_pm(neq,NXs+i-1,NYs+j-1,NZs+k-1)=f(nod)
            enddo
        enddo
    enddo

End Subroutine solve_eq_0_3d
