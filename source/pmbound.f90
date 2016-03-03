!--------------------------------------------------------------------------------
!>@file 
!>@brief Calculations of boundary conditions using particles or sources at domain
!> boundaries
!------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
!-->Subroutine Bounds2d                                                      !
!   This subroutine calculates the boundary conditions for the Solver on PM  !
!   The boundary conditions change for the 2d case                           !
!   The equation solved div(grad)F = P  needs the exact values of F in PM    !
!   boundaries.In 2d:                                                        !
!   For one particle : F = P * (-lnr / (2pi)  )                              !
!   The boundary values is the sum for all particles at each i,j             !
!----------------------------------------------------------------------------!
Subroutine Bounds2d(itype,NXs,NXf,NYs,NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf,neqs,neqf
    integer           :: iconst, jconst, iplane,i,j
    double precision  :: X,Y
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)

    !-->In case of infinite domain bc's(sources are used),In case 
    !-->XMIN
    iconst = Nxs
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_s(iplane,iconst,NYs,NYf,neqs,neqf)
    endif
    !-->XMAX
    iconst = NXf
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        iplane=-iplane
        call calc_boundinf_2d_s(iplane,iconst,NYs,NYf,neqs,neqf)
    endif

    !-->YMIN
    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    jconst = NYs
    iplane = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_s(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

    !-->YMAX
    jconst = NYf
    iplane  = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        iplane=-iplane
        call calc_boundinf_2d_s(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

End Subroutine Bounds2d

!----------------------------------------------------------------------------!
!-->Subroutine Bounds3d                                                      !
!   This subroutine calculates the boundary conditions for the Solver on PM  !
!   The boundary conditions change for the 3d case                           !
!   The equation solved div(grad)F = P  needs the exact values of F in PM    !
!   boundaries.In 2d:                                                        !
!   For one particle : F = P * (-lnr / (2pi)  )                                 !
!   The boundary values is the sum for all particles at each i,j             !
!----------------------------------------------------------------------------!
Subroutine Bounds3d(itype,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer           :: iconst, jconst, kconst, iplane
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)
    !-->iplane is the plane of calculation of the bc's (i.e. for X=const a Y plane is defined)
    !-->iconst defines the poisition of the plane to be calculated
    !-->N*s,N*f is the nodes on the plane to be calculated
    !-->neqs,neqf is the bc's for more than one equations

    !-->XMIN
    iconst = Nxs
    iplane  =  2
    if (itype.eq.1) then
        call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then 
        call calc_boundinf_3d_s(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    endif
    !-->XMAX
    iconst = NXf
    iplane  = 2
    if (itype.eq.1) then
       call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
       call calc_boundinf_3d(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.3) then 
       iplane  = -iplane
       call calc_boundinf_3d_s(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    endif

    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    !-->YMIN
    jconst = NYs
    iplane =  1
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
         call calc_boundinf_3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then 
         call calc_boundinf_3d_s(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    endif
    !-->YMAX
    jconst = NYf
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if(itype.eq.3) then 
        iplane  =  -iplane
        call calc_boundinf_3d_s(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    endif

    !-->ZMIN
    kconst = NZs
    iplane = 3
    if (itype.eq.1) then
        call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then 
    endif
    !-->ZMAX
    kconst = NZf
    if (itype.eq.1) then
       call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
       call calc_boundinf_3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then
       iplane  =  -iplane
       call calc_boundinf_3d_s(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    endif



End Subroutine Bounds3d
!-----------------------------------------------------------------------------!
!->Subroutine calc_bound2d                                                    !
!  This subroutine calculates Biot-Savart Law for 2D applications             !
!  Input:                                                                     !
!        iplane : axis of plane                                               !
!        iconst : Calculation side in nodal value                             !
!-----------------------------------------------------------------------------!
Subroutine calc_bound2d(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r
    integer             :: i, j,nv

    !calculate bc's of all particles on the specified plane defined at iconst
    !-->Y=constant plane
    if (iplane.eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        do nv = 1, NVR
            YR = XP(2,nv) - Y

            do i = Ns, Nf

                X  = XMIN_pm + (i - 1)   * DXpm
                XR = XP(1,nv) - X
                r = sqrt(XR**2 + YR**2)
                !-->Remember Bio Savart :  logr/2pr for grad phi = delta
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1)  + QP(neqs:neqf,nv) * log(r)/pi2
            enddo
        enddo
        !-->X=constant plane
    else if (iplane.eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

        do nv = 1, NVR
            XR = XP(1,nv) - X

            do j = Ns, Nf !Because corner is calculated TWICE

                Y = YMIN_pm + (j - 1) * DYpm
                YR = XP(2,nv) - Y
                r = sqrt(XR**2 + YR**2)
                !  write(*,*) r,pi2
                !-->Remember Bio Savart :  logr/2pi
                SOL_pm(neqs:neqf,iconst, j, 1) = SOL_pm(neqs:neqf,iconst, j, 1)   + QP(neqs:neqf,nv) * log(r)/pi2
            enddo
        enddo

    endif



End Subroutine calc_bound2d

!-----------------------------------------------------------------------------!
!->Subroutine calc_bound3d                                                    !
!  This subroutine calculates Biot-Savart Law for 2D applications             !
!  Input:                                                                     !
!        iplane : axis of plane                                               !
!        iconst : Calculation side in nodal value                             !
!-----------------------------------------------------------------------------!
Subroutine calc_bound3d(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR, Z, ZR, r
    integer             :: i, j, k, nv

    !calculate bc's of all particles on the specified plane defined at iconst
    !-->Y=constant plane
    if (iplane.eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm

        do nv = 1, NVR
            YR = XP(2,nv) - Y
            do k = Ns2,Nf2
                do i = Ns, Nf
                    X  = XMIN_pm + (i - 1) * DXpm
                    XR = XP(1,nv) - X
                    Z  = ZMIN_pm + (k-1) * DZpm
                    ZR = XP(3,nv) - Z
                    r = sqrt(XR**2 + YR**2+ ZR**2)

                    !-->Remember Bio Savart : - 1/(PI4*R) for gradphi= delta
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k)  -&
                                                     QP(neqs:neqf,nv) / (PI4 * r)
                enddo
            enddo
        enddo
        !-->X=constant plane
    else if (iplane.eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

        do nv = 1, NVR
            XR = XP(1,nv) - X
            do k = Ns2,Nf2
                do j = Ns, Nf !Because corner is calculated TWICE

                    Y = YMIN_pm + (j - 1) * DYpm
                    YR = XP(2,nv) - Y
                    Z  = ZMIN_pm + (k-1) * DZpm
                    ZR = XP(3,nv) - Z
                    r = sqrt(XR**2 + YR**2 + ZR**2)

                    SOL_pm(neqs:neqf,iconst, j, k) = SOL_pm(neqs:neqf,iconst, j, k)   - &
                                                     QP(neqs:neqf,nv) / (PI4 * r)
                enddo
            enddo
        enddo

    else if (iplane.eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm

        do nv = 1, NVR
            ZR = XP(3,nv) - Z
            do j = Ns2, Nf2 !Because corner is calculated TWICE
                do i = Ns,Nf

                    X  = XMIN_pm + (i - 1) * DXpm
                    XR = XP(1,nv) - X
                    Y  = YMIN_pm + (j - 1) * DYpm
                    YR = XP(2,nv) - Y
                    r  = sqrt(XR**2 + YR**2 + ZR**2)

                    SOL_pm(neqs:neqf,i, j,iconst) = SOL_pm(neqs:neqf,i, j, iconst)   - QP(neqs:neqf,nv) / (PI4 * r)
                enddo
            enddo
        enddo
    endif



End Subroutine calc_bound3d
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!   Same as particles
!-------------------------------------------------------------------------------!

Subroutine calc_boundinf_2d(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r,DS
    integer             :: i, j,nv

    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        do nv = 1, nbound
            YR = 0.5*(y_s(1,nv)+y_s(2,nv)) - Y

            do i = Ns,Nf
                X  = XMIN_pm + (i - 1)   * DXpm
                XR = 0.5d0 * (x_s(1,nv)+x_s(2,nv)) -  X

                r = sqrt(XR**2 + YR**2)
                DS = d_s(nv)

                !-->Remember Bio Savart :  logr/2pr for grad phi = delta
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + &
                                                 source_bound(neqs:neqf,nv) * log(r) * DS / pi2

            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

        do nv = 1, nbound
            XR = 0.5* (x_s(1,nv)+x_s(2,nv)) - X

            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                YR = 0.5*(y_s(1,nv)+y_s(2,nv)) - Y

                r = sqrt(XR**2 + YR**2)
                DS = d_s(nv)

                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + &
                                                  source_bound(neqs:neqf,nv) * log(r) * DS / pi2

            enddo
        enddo

    endif


End Subroutine calc_boundinf_2d
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!   Same as particles
!-------------------------------------------------------------------------------!

Subroutine calc_boundinf_2d_s(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, greenint,cosb,sinb,DS
    integer             :: i, j,nv

    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        cosb=-1.d0*sign(1,iplane)
        sinb= 0.d0
        do nv = 1, nbound
            XR = x_s(1,nv)
            YR = y_s(1,nv)
            do i = Ns,Nf
                X  = XMIN_pm + (i - 1)   * DXpm

                DS = d_s(nv)

                call PHIELS(X,Y,XR,YR,DS,cosb,sinb,PI2,greenint)
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + &
                                                 source_bound(neqs:neqf,nv) * greenint 
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        cosb=0.d0
        sinb=1.d0*sign(1,iplane)

        do nv = 1, nbound
            XR = x_s(1,nv)
            YR = y_s(1,nv)
            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                DS = d_s(nv)

                call PHIELS(X,Y,XR,YR,DS,cosb,sinb,PI2,greenint)
                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + &
                                                  source_bound(neqs:neqf,nv) * greenint 

            enddo
        enddo

    endif


End Subroutine calc_boundinf_2d_s

!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!-------------------------------------------------------------------------------!
Subroutine calc_boundinf_3d(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS
    integer             :: i, j, k,nv
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm

        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound
             YR = 0.25*(y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv)) - Y
            do k = Ns2,Nf2
                do i = Ns,Nf

                    Z  = ZMIN_pm + (k - 1)   * DZpm
                     ZR = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) -  Z

                    X  = XMIN_pm + (i - 1)   * DXpm
                     XR = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) -  X

                    r = sqrt(XR**2 + YR**2 + ZR**2)
                    DS = d_s(nv)
                    !Green function -1/(4PIR)
                    if (r.gt.1d-05) then 
                        SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) - &
                                                     source_bound(neqs:neqf,nv) * DS / (PI4 * r) 
                    else 

                        SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + &
                         source_bound(neqs:neqf,nv) *  (DXpm+DZpm)*log((sqrt(2.d0)-1)/(sqrt(2.d0)+1))/PI4
                    endif

                enddo
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

        do nv = 1, nbound
              XR = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) -  X
             !XR = x_s(1,nv) -  X

            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    ZR = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) -  Z
                   !ZR = z_s(1,nv) -  Z
                    Y  = YMIN_pm + (j - 1) * DYpm
                    YR = 0.25*(y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv)) - Y
                   !YR =y_s(1,nv)- Y

                    r = sqrt(XR**2 + YR**2 + ZR**2)
                    DS = d_s(nv)

                    if (r.gt.1d-05) then 
                        SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k) - &
                                                      source_bound(neqs:neqf,nv) * DS / (PI4 * r)
                    else 
                        SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k) + &
                        source_bound(neqs:neqf,nv) * (DYpm + DZpm) *log((sqrt(2.d0)-1)/(sqrt(2.d0)+1))/PI4
                    endif

                enddo
            enddo
        enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm

        do nv = 1, nbound
            ZR = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) -  Z
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1)   * DXpm
                    XR = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) -  X
                    Y  = YMIN_pm + (j - 1) * DYpm
                    YR = 0.25*(y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv)) - Y

                    r = sqrt(XR**2 + YR**2 + ZR**2)
                    DS = d_s(nv)

                    if (r.gt.1d-05) then 
                        SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) - &
                                                      source_bound(neqs:neqf,nv) * DS / (PI4 * r)
                    else 
                         SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + &
                         source_bound(neqs:neqf,nv) *  (DXpm+DYpm)*log((sqrt(2.d0)-1)/(sqrt(2.d0)+1))/PI4
                    endif

                enddo
            enddo
        enddo

    endif


End Subroutine calc_boundinf_3d
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!-------------------------------------------------------------------------------!
Subroutine calc_boundinf_3d_s(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS
    integer             :: i, j, k,nv

    double precision    :: XO(3),RG(3),E1(3),E2(3),E3(3),S(4),T(4),SINB(4),COSB(4),D(4),&
                           AREA,DIAG,EPSS,FIS
    integer             :: ISING,NSIDE,si
                    
    !-->Y=constant plane
    ISING=0
    NSIDE=0
    EPSS =1d-14
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DZpm**2)
        E1    = 0.d0; E1(1) =  1.d0 * si 
        E2    = 0.d0; E2(3) =  1.d0 
        E3    = 0.d0; E3(2) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DZpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DXpm;D(3)= DZpm; D(4)=DXpm
        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound
           RG(1) = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) 
           RG(2) = 0.25d0 * (y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv))
           RG(3) = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) 
            do k = Ns2,Nf2
                do i = Ns,Nf
                   X  = XMIN_pm + (i - 1)   * DXpm
                   Z  = ZMIN_pm + (k - 1)   * DZpm
                   XO(1) = X; XO(2) = Y; XO(3) = Z
                   call FSOUR_A4 (XO,RG,E1,E2,E3,&
                                  S,T,D,SINB,COSB,&
                                  DIAG,AREA,PI4,NSIDE,EPSS,ISING,FIS)


                   !Green function -1/(4PIR)
                   SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + &
                                                    FIS* source_bound(neqs:neqf,nv)
                enddo
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        si = sign(1,iplane)
        DIAG = sqrt(DYpm**2+DZpm**2)
        E1    = 0.d0; E1(3) =  1.d0 
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(1) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DYpm * DZpm
        S(1)=-0.5 * DYpm ; S(2)=-0.5 * DYpm ; S(3)=0.5d0 *DYpm ; S(4)= 0.5d0 * DYpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DYpm;D(3)= DZpm; D(4)=DYpm
        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound
           RG(1) = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) 
           RG(2) = 0.25d0 * (y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv))
           RG(3) = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) 
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    Y  = YMIN_pm + (j - 1) * DYpm
                    XO(1) = X; XO(2) = Y; XO(3) = Z
                    !Green function -1/(4PIR)
                    call FSOUR_A4 (XO,RG,E1,E2,E3,&
                                   S,T,D,SINB,COSB,&
                                   DIAG,AREA,PI4,NSIDE,EPSS,ISING,FIS)

                   SOL_pm(neqs:neqf,iconst, j, k) = SOL_pm(neqs:neqf,iconst, j, k) + &
                                                    FIS* source_bound(neqs:neqf,nv)



                enddo
            enddo
        enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DYpm**2)
        E1    = 0.d0; E1(1) =  1.d0  
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(3) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DYpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DYpm ; T(2)= 0.5 * DYpm ; T(3)=0.5d0 *DYpm ; T(4)=-0.5d0 * DYpm 
        D(1)=DYpm; D(2)=DXpm;D(3)= DYpm; D(4)=DXpm
        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound
           RG(1) = 0.25d0 * (x_s(1,nv)+x_s(2,nv)+x_s(3,nv)+x_s(4,nv)) 
           RG(2) = 0.25d0 * (y_s(1,nv)+y_s(2,nv)+y_s(3,nv)+y_s(4,nv))
           RG(3) = 0.25d0 * (z_s(1,nv)+z_s(2,nv)+z_s(3,nv)+z_s(4,nv)) 
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1) * DXpm
                    Y  = YMIN_pm + (j - 1) * DYpm
                    XO(1) = X; XO(2) = Y; XO(3) = Z
                    !Green function -1/(4PIR)
                    call FSOUR_A4 (XO,RG,E1,E2,E3,&
                                   S,T,D,SINB,COSB,&
                                   DIAG,AREA,PI4,NSIDE,EPSS,ISING,FIS)

                    SOL_pm(neqs:neqf, i, j, iconst) = SOL_pm(neqs:neqf, i, j, iconst) + &
                                                    FIS* source_bound(neqs:neqf,nv)



                enddo
            enddo
        enddo

    endif


End Subroutine calc_boundinf_3d_s

Subroutine Bounds2d_lev(itype,NXs,NXf,NYs,NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf,neqs,neqf
    integer           :: iconst, jconst, iplane,i,j
    double precision  :: X,Y
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)

    !-->In case of infinite domain bc's(sources are used),In case 
    !-->XMIN
    iconst = Nxs
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_lev_s(iplane,iconst,NYs,NYf,neqs,neqf)
    endif
    !-->XMAX
    iconst = NXf
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_2d_lev_s(iplane,iconst,NYs,NYf,neqs,neqf)
    endif

    !-->YMIN
    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    jconst = NYs
    iplane = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_lev_s(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

    !-->YMAX
    jconst = NYf
    iplane  = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_2d_lev_s(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

End Subroutine Bounds2d_lev

Subroutine Bounds2d_lev_new(itype,NXs,NXf,NYs,NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf,neqs,neqf
    integer           :: iconst, jconst, iplane,i,j
    double precision  :: X,Y
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)

    !-->In case of infinite domain bc's(sources are used),In case 
    !-->XMIN
    iconst = Nxs
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev_new(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_lev_snew(iplane,iconst,NYs,NYf,neqs,neqf)
    endif
    !-->XMAX
    iconst = NXf
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound2d(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev_new(iplane,iconst,NYs,NYf,neqs,neqf)
       !call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_2d_lev_snew(iplane,iconst,NYs,NYf,neqs,neqf)
       !call calc_boundinf_2d_lev(iplane,iconst,NYs,NYf,neqs,neqf)
    endif

    !-->YMIN
    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    jconst = NYs
    iplane = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev_new(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_2d_lev_snew(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

    !-->YMAX
    jconst = NYf
    iplane  = 1
    if (itype.eq.1) then
        call calc_bound2d(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_2d_lev_new(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_2d_lev_snew(iplane,jconst,NXs+1,NXf-1,neqs,neqf)
    endif

End Subroutine Bounds2d_lev_new
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!   Same as particles
!-------------------------------------------------------------------------------!

Subroutine calc_boundinf_2d_lev(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(neqf)
    integer             :: i, j,nv
    integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm

        do nv = 1, nbound_lev(0)
            do i = Ns,Nf
                X         = XMIN_pm + (i - 1)   * DXpm
                SOURCE    = 0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d(nv,0,leafstart,branch,X,Y,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + SOURCE(neqs:neqf)/ pi2
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        do nv = 1, nbound_lev(0)
            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                SOURCE    =0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d(nv,0,leafstart,branch,X,Y,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + SOURCE(neqs:neqf)/ pi2

            enddo
        enddo

    endif


End Subroutine calc_boundinf_2d_lev

Recursive Subroutine tree_calc_2d(nv,nlev,leafstart,branch,X,Y,SOURCE,neqs,neqf)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: X , Y
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                :: YR, XR, r, DS

   leafstart = leafstart + 2**max(0,(nlev-1))
   leaffin   = leafstart + 2**(nlev) -1
   newlev    = nlev +1 
   leafs     = leafstart + (branch -1) * 2
   leaff     = leafstart +  1
   if(nlev.eq.0) leaff=1
   ibranch   = 0
   do nleaf = leafs,leaff
      ibranch = ibranch + 1
      YR     = ys_lev(nv,0,nleaf) - Y
      XR     = xs_lev(nv,0,nleaf) - X
      r      = sqrt(XR**2 + YR**2)
      DS = ds_lev(nv,0,nleaf)
      if (r.lt.10*DS.and.nlev.lt.levmax) then 
          call tree_calc_2d(nv,nlev + 1,leafstart,ibranch,X,Y,SOURCE,neqs,neqf)
      else 
          SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nv,neqs:neqf,0,nleaf) * DS * log(r)
      endif
   enddo

End Subroutine tree_calc_2d


Subroutine calc_boundinf_2d_lev_new(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(neqf)
    integer             :: i, j,nv
    integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
      ! do nv = 1, nbound_lev(levmax)
         nv=1
            do i = Ns,Nf
                X         = XMIN_pm + (i - 1)   * DXpm
                SOURCE    = 0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_new(nv,levmax,leafstart,branch,X,Y,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + SOURCE(neqs:neqf)/ pi2
            enddo
      ! enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
      ! do nv = 1, nbound_lev(levmax)
           nv=1
            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                SOURCE    =0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_new(nv,levmax,leafstart,branch,X,Y,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + SOURCE(neqs:neqf)/ pi2

            enddo
      ! enddo

    endif


End Subroutine calc_boundinf_2d_lev_new

Recursive Subroutine tree_calc_2d_new(nv,nlev,leafstart,branch,X,Y,SOURCE,neqs,neqf)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: X , Y
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                :: YR, XR, r, DS

   leaffin   = min(2*nv,nbound_lev(nlev))!in case of ilev_t.eq.-1!leafstart + 2**(nlev) -1
   leafstart = 2*nv -1
   if (nlev.eq.levmax) then 
      leafstart=1
      leaffin = nbound_lev(nlev)
   endif
   newlev    = nlev -1 
   ibranch   = 0
 ! write(*,*) nlev,leafstart,leaffin
 ! read(*,*) 
   if (leafstart.gt.nbound_lev(nlev)) then 
        leafstart=leaffin
   endif
   do nleaf = leafstart,leaffin
      ibranch = ibranch + 1
      YR     = ys_lev(nleaf,nlev,2) - Y
      XR     = xs_lev(nleaf,nlev,2) - X
      r      = sqrt(XR**2 + YR**2)
      DS = ds_lev(nleaf,nlev,2)
      if ((r.lt.10*DS.and.nlev.gt.0).or.ilev_t(nleaf,nlev).eq.-1) then 
          call tree_calc_2d_new(nleaf,newlev,leafstart,ibranch,X,Y,SOURCE,neqs,neqf)
      else 
          SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nleaf,neqs:neqf,nlev,2) * DS * log(r)
      endif
   enddo
End Subroutine tree_calc_2d_new
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!   Same as particles
!-------------------------------------------------------------------------------!

Subroutine calc_boundinf_2d_lev_s(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r, a,b,ra,rb,greenint,cosb,sinb,DS,SOURCE(neqf)
    integer             :: i, j,nv
    integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        cosb=-1.d0*sign(1,iplane)
        sinb= 0.d0
        do nv = 1, nbound_lev(0)
            do i = Ns,Nf
                X         = XMIN_pm + (i - 1)   * DXpm
                SOURCE    = 0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_s(nv,0,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + SOURCE(neqs:neqf)
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        cosb=0.d0
        sinb=1.d0*sign(1,iplane)
        do nv = 1, nbound_lev(0)
            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                SOURCE    =0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_s(nv,0,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + SOURCE(neqs:neqf)

            enddo
        enddo

    endif


End Subroutine calc_boundinf_2d_lev_s

Recursive Subroutine tree_calc_2d_s(nv,nlev,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
  Implicit None
  integer, intent(in) :: nlev,nv,neqs,neqf,branch
  integer, intent(inout) :: leafstart
  double precision, intent(in)  :: X , Y, cosb, sinb,PI2
  double precision, intent(inout) :: SOURCE(neqf)
  
  integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
  double precision                :: YR, XR, r, DS,greenint

  leafstart = leafstart + 2**max(0,(nlev-1))
  leaffin   = leafstart + 2**(nlev) -1
  newlev    = nlev +1 
  leafs     = leafstart + (branch -1) * 2
  leaff     = leafstart +  1
  if(nlev.eq.0) leaff=1
  ibranch   = 0
  do nleaf = leafs,leaff
     ibranch = ibranch + 1
     DS = ds_lev(nv,0,nleaf)
     YR     = ys_lev(nv,0,nleaf) 
     XR     = xs_lev(nv,0,nleaf) 
     r      = sqrt((XR - X)**2 + (YR - Y)**2)
     if (r.lt.4*sqrt(DS).and.nlev.lt.levmax) then 
         call tree_calc_2d_s(nv,nlev + 1,leafstart,ibranch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
     else 
       if (nlev.lt.levmax) then
           SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nv,neqs:neqf,0,nleaf)& 
                                                     * DS * log(r)/PI2
       else  
           call PHIELS(X,Y,XR,YR,DS,cosb,sinb,PI2,greenint)
           SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nv,neqs:neqf,0,nleaf)&
                                                     * greenint
       endif 
     endif
  enddo

End Subroutine tree_calc_2d_s

!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!   Same as particles
!-------------------------------------------------------------------------------!

Subroutine calc_boundinf_2d_lev_snew(iplane,iconst,Ns,Nf,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,neqs,neqf

    double precision    :: X, Y, XR, YR, r, a,b,ra,rb,greenint,cosb,sinb,DS,SOURCE(neqf)
    integer             :: i, j,nv
    integer             :: leafstart, leaffin, lev, nlev, nleaf, branch
    !calculate bc's of all sources on the specified plane defined at iconst
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        cosb=-1.d0*sign(1,iplane)
        sinb= 0.d0
            nv=1
            do i = Ns,Nf
                X         = XMIN_pm + (i - 1)   * DXpm
                SOURCE    = 0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_snew(nv,levmax,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,i, iconst, 1) = SOL_pm(neqs:neqf,i, iconst, 1) + SOURCE(neqs:neqf)
            enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        cosb=0.d0
        sinb=1.d0*sign(1,iplane)
            nv=1
            do j = Ns,Nf !Because corners are calculated twice
                Y  = YMIN_pm + (j - 1) * DYpm
                SOURCE    =0.d0
                leafstart = 0
                branch    = 1
                call tree_calc_2d_snew(nv,levmax,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
                SOL_pm(neqs:neqf,iconst, j , 1) = SOL_pm(neqs:neqf,iconst,j, 1) + SOURCE(neqs:neqf)

            enddo

    endif


End Subroutine calc_boundinf_2d_lev_snew

Recursive Subroutine tree_calc_2d_snew(nv,nlev,leafstart,branch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
  Implicit None
  integer, intent(in) :: nlev,nv,neqs,neqf,branch
  integer, intent(inout) :: leafstart
  double precision, intent(in)  :: X , Y, cosb, sinb,PI2
  double precision, intent(inout) :: SOURCE(neqf)
  
  integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
  double precision                :: YR, XR, r, DS,greenint


   leaffin   = min(2*nv,nbound_lev(nlev))!in case of ilev_t.eq.-1!leafstart + 2**(nlev) -1
   leafstart = 2*nv -1
   if (nlev.eq.levmax) then 
      leafstart=1
      leaffin = nbound_lev(nlev)
   endif
   newlev    = nlev -1 
   ibranch   = 0
 ! write(*,*) nlev,leafstart,leaffin
 ! read(*,*) 

   if (leafstart.gt.nbound_lev(nlev)) then 
        leafstart=leaffin
   endif
  do nleaf = leafstart,leaffin
     ibranch = ibranch + 1
     DS = ds_lev(nleaf,nlev,2)
     YR     = ys_lev(nleaf,nlev,2) 
     XR     = xs_lev(nleaf,nlev,2) 
     r      = sqrt((XR - X)**2 + (YR - Y)**2)
     if (r.lt.10*sqrt(DS).and.nlev.gt.0) then 
         call tree_calc_2d_snew(nleaf,newlev,leafstart,ibranch,X,Y,cosb,sinb,PI2,SOURCE,neqs,neqf)
     else 
       if (nlev.gt.0) then
           SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nleaf,neqs:neqf,nlev,2)& 
                                                     * DS * log(r)/PI2
       else  
           call PHIELS(X,Y,XR,YR,DS,cosb,sinb,PI2,greenint)
           SOURCE(neqs : neqf) = SOURCE(neqs : neqf) + source_bound_lev(nleaf,neqs:neqf,nlev,2)&
                                                     * greenint
       endif 
     endif
  enddo

End Subroutine tree_calc_2d_snew

!----------------------------------------------------------------------------!
!-->Subroutine Bounds3d                                                      !
!   This subroutine calculates the boundary conditions for the Solver on PM  !
!   The boundary conditions change for the 3d case                           !
!   The equation solved div(grad)F = P  needs the exact values of F in PM    !
!   boundaries.In 2d:                                                        !
!   For one particle : F = P * (-lnr / (2pi)  )                                 !
!   The boundary values is the sum for all particles at each i,j             !
!----------------------------------------------------------------------------!
Subroutine Bounds3d_lev(itype,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer           :: iconst, jconst, kconst, iplane
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)
    !-->iplane is the plane of calculation of the bc's (i.e. for X=const a Y plane is defined)
    !-->iconst defines the poisition of the plane to be calculated
    !-->N*s,N*f is the nodes on the plane to be calculated
    !-->neqs,neqf is the bc's for more than one equations

    !-->XMIN
    iconst = Nxs
    iplane  = 2
    if (itype.eq.1) then
        call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_s(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    endif
    !-->XMAX
    iconst = NXf
    if (itype.eq.1) then
        call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.3) then
        iplane = -iplane
        call calc_boundinf_3d_lev_s(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    endif

    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    !-->YMIN
    jconst = NYs
    iplane = 1
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_s(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    endif
    !-->YMAX
    jconst = NYf
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        iplane =-iplane
        call calc_boundinf_3d_lev_s(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    endif

    !-->ZMIN
    kconst = NZs
    iplane = 3
    if (itype.eq.1) then
        call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_s(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    endif
    !-->ZMAX
    kconst = NZf
    if (itype.eq.1) then
        call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_3d_lev_s(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    endif



End Subroutine Bounds3d_lev
!----------------------------------------------------------------------------!
!-->Subroutine Bounds3d                                                      !
!   This subroutine calculates the boundary conditions for the Solver on PM  !
!   The boundary conditions change for the 3d case                           !
!   The equation solved div(grad)F = P  needs the exact values of F in PM    !
!   boundaries.In 2d:                                                        !
!   For one particle : F = P * (-lnr / (2pi)  )                                 !
!   The boundary values is the sum for all particles at each i,j             !
!----------------------------------------------------------------------------!
Subroutine Bounds3d_lev_new(itype,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: itype, NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer           :: iconst, jconst, kconst, iplane,ip,jp,kp,i,j,k
    double precision  :: XMIN,YMIN,ZMIN,XMIN0,YMIN0,ZMIN0,DDX,DDY,DDZ
    !-->Calculate boundary conditions for each boundary (XMIN,XMAX,YMIN,YMAX)
    !-->iplane is the plane of calculation of the bc's (i.e. for X=const a Y plane is defined)
    !-->iconst defines the poisition of the plane to be calculated
    !-->N*s,N*f is the nodes on the plane to be calculated
    !-->neqs,neqf is the bc's for more than one equations

    !-->XMIN
    iconst = Nxs
    iplane  = 2
    DDX= XMIN_pm + (NXf_bl(1)-NXs_bl(1))*DXpm
    DDY= YMIN_pm + (NYf_bl(1)-NYs_bl(1))*DYpm
    DDZ= ZMIN_pm + (NZf_bl(1)-NZs_bl(1))*DZpm
    ip=1;jp=3;kp=3
    XMIN0 = XMIN_pm -(nint(ip/2.)-1)*DDX
    YMIN0 = YMIN_pm -(nint(jp/2.)-1)*DDY
    ZMIN0 = ZMIN_pm -(nint(kp/2.)-1)*DDZ
    if (itype.eq.1) then
        call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_snew(iplane ,iconst ,NYs ,NYf ,NZs,NZf,neqs,neqf)
    else if (itype.eq.1002) then 
         
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,iconst,NYs,NYf,NZs,NZf,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif
    !-->XMAX
    iconst = NXf
    if (itype.eq.1) then
        call calc_bound3d(iplane ,iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.3) then
        iplane = -iplane
        call calc_boundinf_3d_lev_snew(iplane, iconst, NYs, NYf, NZs, NZf,neqs,neqf)
    else if (itype.eq.1002) then 
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,iconst,NYs,NYf,NZs,NZf,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif

    !We use Nxs + 1,NXf - 1 For corners since they already calculated
    !-->YMIN
    jconst = NYs
    iplane = 1
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_snew(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.1002) then 
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif
    !-->YMAX
    jconst = NYf
    if (itype.eq.1) then
        call calc_bound3d(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.3) then
        iplane =-iplane
        call calc_boundinf_3d_lev_snew(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf)
    else if (itype.eq.1002) then 
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,jconst,NXs+1,NXf-1,NZs,NZf,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif

    !-->ZMIN
    kconst = NZs
    iplane = 3
    if (itype.eq.1) then
        call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then
        call calc_boundinf_3d_lev_snew(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.1002) then 
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,kconst,NXs+1,NXf-1,NZs+1,NZf-1,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif
    !-->ZMAX
    kconst = NZf
    if (itype.eq.1) then
        call calc_bound3d(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.2) then
        call calc_boundinf_3d_lev_new(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.3) then
        iplane  = -iplane
        call calc_boundinf_3d_lev_snew(iplane,kconst,NXs+1,NXf-1,NYs+1,NYf-1,neqs,neqf)
    else if (itype.eq.1002) then 
        do k=1,kp
           do j= 1,jp
              do i= 1,ip
                 XMIN = XMIN0 + (i-1)*DDX
                 YMIN = YMIN0 + (j-1)*DDY
                 ZMIN = ZMIN0 + (k-1)*DDZ
                 call calc_boundinf_3d_lev_newper(iplane,kconst,NXs+1,NXf-1,NZs+1,NZf-1,neqs,neqf,XMIN,YMIN,ZMIN)
              enddo
           enddo
        enddo
    endif



End Subroutine Bounds3d_lev_new
!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!-------------------------------------------------------------------------------!
Subroutine calc_boundinf_3d_lev(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(1:neqf)
    integer             :: i, j, k,nv
    integer             :: leafstart,branch
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound_lev(0)
            do k = Ns2,Nf2
                do i = Ns,Nf
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    X  = XMIN_pm + (i - 1)   * DXpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d(nv,0,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + SOURCE(neqs:neqf)/PI4 
                enddo
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

        do nv = 1, nbound_lev(0)
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    Y  = YMIN_pm + (j - 1) * DYpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d(nv,0,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k)  + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
        enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm

        do nv = 1, nbound_lev(0)
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1)   * DXpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1
                    call tree_calc_3d(nv,0,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
        enddo

    endif


End Subroutine calc_boundinf_3d_lev

Recursive Subroutine tree_calc_3d(nv,nlev,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: X , Y, Z
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                :: XR, YR, ZR, r, DS

   leafstart = leafstart + 4**max(0,(nlev-1))
   leaffin   = leafstart + 4**(nlev) -1
   newlev    = nlev +1 
   leafs     = leafstart + (branch -1) * 4
   leaff     = leafstart +  1
   if(nlev.eq.0) leaff=1
   ibranch   = 0
   do nleaf = leafs,leaff
      ibranch = ibranch + 1
      XR     = xs_lev(nv,0,nleaf) - X
      YR     = ys_lev(nv,0,nleaf) - Y
      ZR     = zs_lev(nv,0,nleaf) - Z
      r      = sqrt(XR**2 + YR**2 + ZR**2)
      DS = ds_lev(nv,0,nleaf)
             !if(nv.eq.1) then 
             !   write(*,*)'--------------------'
             !   write(*,*)'nlev=', nlev
             !   write(*,*) DS,r  
             !   write(*,*)xs_lev(nv,0,nleaf),ys_lev(nv,0,nleaf),zs_lev(nv,0,nleaf)
             !   write(*,*) X,Y,Z
             !   write(*,*)'--------------------'
             !endif
      if (r.lt.10*sqrt(DS).and.nlev.lt.levmax) then 
          call tree_calc_3d(nv,newlev,leafstart,ibranch,X,Y,Z,SOURCE,neqs,neqf)
      else 
          !Green function  -1/(PI4*R)
          if (r.gt.1d-05) then
              SOURCE(neqs : neqf) = SOURCE(neqs : neqf) - source_bound_lev(nv,neqs:neqf,0,nleaf) * DS/r
          else
              SOURCE(neqs : neqf) = SOURCE(neqs : neqf) +&
              source_bound_lev(nv,neqs:neqf,0,nleaf) *  2* sqrt(DS)*log((sqrt(2.d0)-1)/(sqrt(2.d0)+1))
          endif

      endif
   enddo

End Subroutine tree_calc_3d

Subroutine calc_boundinf_3d_lev_new(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(1:neqf)
    integer             :: i, j, k,nv
    integer             :: leafstart,branch
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        !calculate bc's of all sources on the specified plane defined at iconst
      ! do nv = 1, nbound_lev(0)
            nv=1
            do k = Ns2,Nf2
                do i = Ns,Nf
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    X  = XMIN_pm + (i - 1)   * DXpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + SOURCE(neqs:neqf)/PI4 
                enddo
            enddo
     !  enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm

     !  do nv = 1, nbound_lev(0)
            nv=1
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    Y  = YMIN_pm + (j - 1) * DYpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k)  + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
     !  enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm

       !do nv = 1, nbound_lev(0)
            nv=1
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1)   * DXpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1
                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
       !enddo

    endif


End Subroutine calc_boundinf_3d_lev_new

Recursive Subroutine tree_calc_3d_new(nv,nlev,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: X , Y, Z
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                :: XR, YR, ZR, r, DS,ss(neqf)
   leaffin   = min(4*nv,nbound_lev(nlev))!in case of ilev_t.eq.-1!leafstart + 2**(nlev) -1
   leafstart = 4*nv -3
   if (nlev.eq.levmax) then 
      leafstart=1
      leaffin = nbound_lev(nlev)
   endif
   newlev    = nlev -1 
   ibranch   = 0
 ! write(*,*) nlev,leafstart,leaffin
 ! read(*,*) 

   do nleaf = leafstart,leaffin
      ibranch = ibranch + 1
      XR     = xs_lev(nleaf,nlev,2) - X
      YR     = ys_lev(nleaf,nlev,2) - Y
      ZR     = zs_lev(nleaf,nlev,2) - Z
      r      = sqrt(XR**2 + YR**2 + ZR**2)
      DS = ds_lev(nleaf,nlev,2)
             !if(nv.eq.1) then 
             !   write(*,*)'--------------------'
             !   write(*,*)'nlev=', nlev
             !   write(*,*) DS,r  
             !   write(*,*)xs_lev(nv,0,nleaf),ys_lev(nv,0,nleaf),zs_lev(nv,0,nleaf)
             !   write(*,*) X,Y,Z
             !   write(*,*)'--------------------'
             !endif
      if (r.lt.10*sqrt(DS).and.nlev.gt.0) then 
          call tree_calc_3d_new(nleaf,newlev,leafstart,ibranch,X,Y,Z,SOURCE,neqs,neqf)
      else 
          !Green function  -1/(PI4*R)
          if (r.gt.1d-05) then
           !ss =source_bound_lev(nleaf,neqs:neqf,nlev,2)

           !if (ss(1).ge.0.or.ss.le.0.eqv.(.FALSE.)) then
           ! write(*,*) 'NaN error',ss
           !  STOP
           !endif

              SOURCE(neqs : neqf) = SOURCE(neqs : neqf) - source_bound_lev(nleaf,neqs:neqf,nlev,2) * DS/r
          else
              SOURCE(neqs : neqf) = SOURCE(neqs : neqf) +&
              source_bound_lev(nleaf,neqs:neqf,nlev,2) *  2* sqrt(DS)*log((sqrt(2.d0)-1)/(sqrt(2.d0)+1))
          endif

      endif
   enddo

End Subroutine tree_calc_3d_new

!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!-------------------------------------------------------------------------------!
Subroutine calc_boundinf_3d_lev_s(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(1:neqf)
    integer             :: i, j, k,nv
    integer             :: leafstart,branch

    double precision    :: XO(3),RG(3),E1(3),E2(3),E3(3),S(4),T(4),SINB(4),COSB(4),D(4),&
                           AREA,DIAG,EPSS,FIS
    integer             :: ISING,NSIDE,si
                    
    !-->Y=constant plane
    ISING=0
    NSIDE=0
    EPSS =1d-14
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DZpm**2)
        E1    = 0.d0; E1(1) =  1.d0 * si 
        E2    = 0.d0; E2(3) =  1.d0 
        E3    = 0.d0; E3(2) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DZpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DXpm;D(3)= DZpm; D(4)=DXpm

        !calculate bc's of all sources on the specified plane defined at iconst
        do nv = 1, nbound_lev(0)
            do k = Ns2,Nf2
                do i = Ns,Nf
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    X  = XMIN_pm + (i - 1)   * DXpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_s(nv,0,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + SOURCE(neqs:neqf) 
                enddo
            enddo
        enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        si = sign(1,iplane)
        DIAG = sqrt(DYpm**2+DZpm**2)
        E1    = 0.d0; E1(3) =  1.d0 
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(1) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DYpm * DZpm
        S(1)=-0.5 * DYpm ; S(2)=-0.5 * DYpm ; S(3)=0.5d0 *DYpm ; S(4)= 0.5d0 * DYpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DYpm;D(3)= DZpm; D(4)=DYpm

        do nv = 1, nbound_lev(0)
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_s(nv,0,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
                    SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k)  + SOURCE(neqs:neqf)

                enddo
            enddo
        enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm
        Z = ZMIN_pm + (iconst - 1) * DZpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DYpm**2)
        E1    = 0.d0; E1(1) =  1.d0  
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(3) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DYpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DYpm ; T(2)= 0.5 * DYpm ; T(3)=0.5d0 *DYpm ; T(4)=-0.5d0 * DYpm 
        D(1)=DYpm; D(2)=DXpm;D(3)= DYpm; D(4)=DXpm

        do nv = 1, nbound_lev(0)
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1)   * DXpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1
                    call tree_calc_3d_s(nv,0,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)


                    SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + SOURCE(neqs:neqf)

                enddo
            enddo
        enddo

    endif


End Subroutine calc_boundinf_3d_lev_s

Recursive Subroutine tree_calc_3d_s(nv,nlev,leafstart,branch,XO,SOURCE,neqs,neqf,&
                     DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: DIAG,E1(3),E2(3),E3(3),PI4
   double precision, intent(in)  :: COSB(4),SINB(4),S(4),T(4),D(4)
   double precision, intent(in)  :: EPSS
   integer         , intent(in)  :: ISING,NSIDE
   double precision, intent(in)  :: XO(3)
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                ::r, DS, RG(3),FIS,RATIO

   leafstart = leafstart + 4**max(0,(nlev-1))
   leaffin   = leafstart + 4**(nlev) -1
   newlev    = nlev +1 
   leafs     = leafstart + (branch -1) * 4
   leaff     = leafstart +  1
   if(nlev.eq.0) leaff=1
   ibranch   = 0
   do nleaf = leafs,leaff
      ibranch = ibranch + 1
      RG(1) = xs_lev(nv,0,nleaf) 
      RG(2) = ys_lev(nv,0,nleaf) 
      RG(3) = zs_lev(nv,0,nleaf) 
      r      = sqrt((XO(1)-RG(1))**2 + (XO(2)-RG(2))**2 + (XO(3)-RG(3))**2)
      DS     = ds_lev(nv,0,nleaf)
      RATIO  = r/sqrt(DS)
      if (RATIO.lt.10.and.nlev.lt.levmax) then 
          call tree_calc_3d_s(nv,newlev,leafstart,ibranch,XO,SOURCE,neqs,neqf,&
                            DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
      else 
       !write(*,*) newlev,ratio,r,sqrt(DS)
        if (nlev.ne.levmax) then 
          SOURCE = SOURCE  - source_bound_lev(nv,neqs:neqf,0,nleaf) * DS/(PI4*r)
        else
          call FSOUR_A4 (XO,RG,E1,E2,E3,&
                         S,T,D,SINB,COSB,&
                         DIAG,DS,PI4,NSIDE,EPSS,ISING,FIS)
          SOURCE = SOURCE + source_bound_lev(nv,neqs:neqf,0,nleaf) * FIS
        endif
      endif

   enddo

End Subroutine tree_calc_3d_s

!-------------------------------------------------------------------------------!
!-> Subroutine calc_boundinf                                                    !
!   This Subroutine calculates boundary conditions for the sources              !
!-------------------------------------------------------------------------------!
Subroutine calc_boundinf_3d_lev_snew(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf

    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(1:neqf)
    integer             :: i, j, k,nv
    integer             :: leafstart,branch

    double precision    :: XO(3),RG(3),E1(3),E2(3),E3(3),S(4),T(4),SINB(4),COSB(4),D(4),&
                           AREA,DIAG,EPSS,FIS
    integer             :: ISING,NSIDE,si
                    
    !-->Y=constant plane
    ISING=0
    NSIDE=0
    EPSS =1d-14
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN_pm + (iconst - 1) * DYpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DZpm**2)
        E1    = 0.d0; E1(1) =  1.d0 * si 
        E2    = 0.d0; E2(3) =  1.d0 
        E3    = 0.d0; E3(2) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DZpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DXpm;D(3)= DZpm; D(4)=DXpm

        !calculate bc's of all sources on the specified plane defined at iconst
            nv=1
            do k = Ns2,Nf2
                do i = Ns,Nf
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    X  = XMIN_pm + (i - 1)   * DXpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_snew(nv,levmax,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + SOURCE(neqs:neqf) 
                enddo
            enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN_pm + (iconst - 1) * DXpm
        si = sign(1,iplane)
        DIAG = sqrt(DYpm**2+DZpm**2)
        E1    = 0.d0; E1(3) =  1.d0 
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(1) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DYpm * DZpm
        S(1)=-0.5 * DYpm ; S(2)=-0.5 * DYpm ; S(3)=0.5d0 *DYpm ; S(4)= 0.5d0 * DYpm 
        T(1)=-0.5 * DZpm ; T(2)= 0.5 * DZpm ; T(3)=0.5d0 *DZpm ; T(4)=-0.5d0 * DZpm 
        D(1)=DZpm; D(2)=DYpm;D(3)= DZpm; D(4)=DYpm

            nv=1
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN_pm + (k - 1)   * DZpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_snew(nv,levmax,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
                    SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k)  + SOURCE(neqs:neqf)

                enddo
            enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN_pm + (iconst - 1) * DZpm
        Z = ZMIN_pm + (iconst - 1) * DZpm
        si = sign(1,iplane)
        DIAG = sqrt(DXpm**2+DYpm**2)
        E1    = 0.d0; E1(1) =  1.d0  
        E2    = 0.d0; E2(2) = -1.d0 * si
        E3    = 0.d0; E3(3) = -1.d0 * si
        COSB  = 0.d0; COSB(2)=1.d0 ; COSB(4)=-1.d0 
        SINB  = 0.d0; SINB(1)=1.d0 ; SINB(3)=-1.d0 
        AREA  = DXpm * DYpm
        S(1)=-0.5 * DXpm ; S(2)=-0.5 * DXpm ; S(3)=0.5d0 *DXpm ; S(4)= 0.5d0 * DXpm 
        T(1)=-0.5 * DYpm ; T(2)= 0.5 * DYpm ; T(3)=0.5d0 *DYpm ; T(4)=-0.5d0 * DYpm 
        D(1)=DYpm; D(2)=DXpm;D(3)= DYpm; D(4)=DXpm

            nv=1
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN_pm + (i - 1)   * DXpm
                    Y  = YMIN_pm + (j - 1) * DYpm

                    XO(1) = X; XO(2) = Y; XO(3) = Z

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1
                    call tree_calc_3d_snew(nv,levmax,leafstart,branch,XO,SOURCE,neqs,neqf,&
                                        DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)


                    SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + SOURCE(neqs:neqf)

                enddo
            enddo

    endif


End Subroutine calc_boundinf_3d_lev_snew

Recursive Subroutine tree_calc_3d_snew(nv,nlev,leafstart,branch,XO,SOURCE,neqs,neqf,&
                     DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
   Implicit None
   integer, intent(in) :: nlev,nv,neqs,neqf,branch
   integer, intent(inout) :: leafstart
   double precision, intent(in)  :: DIAG,E1(3),E2(3),E3(3),PI4
   double precision, intent(in)  :: COSB(4),SINB(4),S(4),T(4),D(4)
   double precision, intent(in)  :: EPSS
   integer         , intent(in)  :: ISING,NSIDE
   double precision, intent(in)  :: XO(3)
   double precision, intent(inout) :: SOURCE(neqf)
   
   integer                         :: newlev, nleaf, leaffin,ibranch,leafs,leaff
   double precision                ::r, DS, RG(3),FIS,RATIO

   leaffin   = min(4*nv,nbound_lev(nlev))!in case of ilev_t.eq.-1!leafstart + 2**(nlev) -1
   leafstart = 4*nv -3
   if (nlev.eq.levmax) then 
      leafstart=1
      leaffin = nbound_lev(nlev)
   endif
   newlev    = nlev -1 
   ibranch   = 0

   do nleaf = leafstart,leaffin
      ibranch = ibranch + 1
      RG(1) = xs_lev(nleaf,nlev,2) 
      RG(2) = ys_lev(nleaf,nlev,2) 
      RG(3) = zs_lev(nleaf,nlev,2) 
      r      = sqrt((XO(1)-RG(1))**2 + (XO(2)-RG(2))**2 + (XO(3)-RG(3))**2)
      DS     = ds_lev(nleaf,nlev,2)
      RATIO  = r/sqrt(DS)
      if (RATIO.lt.10.and.nlev.gt.0) then 
          call tree_calc_3d_snew(nleaf,newlev,leafstart,ibranch,XO,SOURCE,neqs,neqf,&
                            DIAG,E1,E2,E3,COSB,SINB,S,T,D,NSIDE,EPSS,ISING,PI4)
      else 
       !write(*,*) newlev,ratio,r,sqrt(DS)
        if (nlev.ne.levmax) then 
          SOURCE = SOURCE  - source_bound_lev(nleaf,neqs:neqf,nlev,2) * DS/(PI4*r)
        else
          call FSOUR_A4 (XO,RG,E1,E2,E3,&
                         S,T,D,SINB,COSB,&
                         DIAG,DS,PI4,NSIDE,EPSS,ISING,FIS)
          SOURCE = SOURCE + source_bound_lev(nv,neqs:neqf,nlev,2) * FIS
        endif
      endif

   enddo

End Subroutine tree_calc_3d_snew

Subroutine calc_boundinf_3d_lev_newper(iplane,iconst,Ns,Nf,Ns2,Nf2,neqs,neqf,XMIN,YMIN,ZMIN)
    Implicit None

    integer, intent(in) :: iplane, iconst,Ns,Nf,Ns2,Nf2,neqs,neqf
    double precision,intent(in) :: XMIN,YMIN,ZMIN
    double precision    :: X, Y, XR, YR,Z, ZR,r, a,b,ra,rb,greenint,racos,rasin,DS,SOURCE(1:neqf)
    integer             :: i, j, k,nv
    integer             :: leafstart,branch
    !-->Y=constant plane
    if (abs(iplane).eq.1) then
        Y = YMIN + (iconst - 1) * DYpm
        !calculate bc's of all sources on the specified plane defined at iconst
      ! do nv = 1, nbound_lev(0)
            nv=1
            do k = Ns2,Nf2
                do i = Ns,Nf
                    Z  = ZMIN + (k - 1)   * DZpm
                    X  = XMIN + (i - 1)   * DXpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, iconst, k) = SOL_pm(neqs:neqf,i, iconst, k) + SOURCE(neqs:neqf)/PI4 
                enddo
            enddo
     !  enddo
        !-->X=constant plane
    else if (abs(iplane).eq.2) then
        X = XMIN + (iconst - 1) * DXpm

     !  do nv = 1, nbound_lev(0)
            nv=1
            do k = Ns2,Nf2
                do j = Ns,Nf !Because corners are calculated twice
                    Z  = ZMIN + (k - 1)   * DZpm
                    Y  = YMIN + (j - 1) * DYpm


                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1

                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,iconst, j , k) = SOL_pm(neqs:neqf,iconst,j, k)  + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
     !  enddo

    else if (abs(iplane).eq.3) then
        Z = ZMIN + (iconst - 1) * DZpm

       !do nv = 1, nbound_lev(0)
            nv=1
            do j = Ns2,Nf2 !Because corners are calculated twice
                do i = Ns,Nf
                    X  = XMIN + (i - 1)   * DXpm
                    Y  = YMIN + (j - 1) * DYpm

                    SOURCE    =0.d0
                    leafstart = 0
                    branch    = 1
                    call tree_calc_3d_new(nv,levmax,leafstart,branch,X,Y,Z,SOURCE,neqs,neqf)
                    SOL_pm(neqs:neqf,i, j , iconst) = SOL_pm(neqs:neqf,i,j, iconst) + SOURCE(neqs:neqf)/PI4 

                enddo
            enddo
       !enddo

    endif


End Subroutine calc_boundinf_3d_lev_newper

!Subroutine PHIELS calculates the potential induced by constant panels'
! XO is the point of calculation X1,Y1 is the first corner of the constant panel
! DS is the area of the face cosb,sinb give the direction assumed.
Subroutine PHIELS(X0,Y0,X1,Y1,DS,COSB,SINB,PI2,PHILS)
     Implicit None 
     double precision,intent(in) :: X0,Y0,X1,Y1,DS,COSB,SINB,PI2
     double precision,intent(out):: PHILS
     double precision            :: AKSIL,HTAL,TA1,TA2
     AKSIL= (X0-X1)*COSB+(Y0-Y1)*SINB !this is the vector X0-X1 in local coor
     HTAL =-(X0-X1)*SINB+(Y0-Y1)*COSB
     TA1=AKSIL ! vector XO-X1 from the first point
     TA2=AKSIL-DS!vector XO-X2 (since in l.c. X2=X1+DS then X0-X2 = XO-X1 - DS
     if  (dabs(HTAL).gt.1.d-08) then
          PHILS= 1./(PI2) * ( TA1*DLOG(DSQRT(TA1**2+HTAL**2))&     
                               -TA2*DLOG(DSQRT(TA2**2+HTAL**2))&
                               +HTAL*(DATAN(TA1/HTAL)-DATAN(TA2/HTAL))-DS)
     else
          if(abs(TA1).lt.1d-08) then 
             PHILS=1.d0/(2.d0*PI2)*(-TA2*DLOG(TA2**2) -2.d0*DS)
          else if(abs(TA2).lt.1d-08) then 
             PHILS=1.d0/(2.d0*PI2)*(TA1*DLOG(TA1**2)   -2d0*DS)
          else 
            PHILS=1.d0/(2.d0*PI2)*(TA1*DLOG(TA1**2)-TA2*DLOG(TA2**2)-2.d0*DS)
          endif
     end if

     return
End Subroutine PHIELS
                      
!-----------------------------------------------------------------------
!
!
!     Chapter 1. FLAT CONSTANT SOURCE ELEMENT
!     ---------------------------------------
!
!-----------------------------------------------------------------------
!
!  Subr      :FSOUR_A4 
   SUBROUTINE FSOUR_A4 (XO,RG,E1,E2,E3,&
                        S,T,D,SINB,COSB,&
                        DIAG,AREA,PI4,NSIDE,EPSS,ISING,FIS)


   implicit none

   integer, intent(in)           :: NSIDE, ISING
   double precision, intent(in)  :: DIAG, AREA, PI4, EPSS
   double precision, intent(in)  :: XO(3),RG(3),E1(3),E2(3),E3(3),S(4),T(4), &
                                    D(4),SINB(4),COSB(4)
   double precision, intent(out) :: FIS
   
   integer             :: K,K1,K2
   double precision    :: X,Y,Z,RO,RATIO,AREA1,RO3,P,Q,FIL,&
                          XK,YK,ZK,A1,A2,AK,ZP,AZ,R(4),E(4),H(4),UL(3),TOO(4),TINY



   TINY=1d-10
   X  = (XO(1)-RG(1))*E1(1) + (XO(2)-RG(2))*E1(2) + (XO(3)-RG(3))*E1(3)
   Y  = (XO(1)-RG(1))*E2(1) + (XO(2)-RG(2))*E2(2) + (XO(3)-RG(3))*E2(3)
   Z  = (XO(1)-RG(1))*E3(1) + (XO(2)-RG(2))*E3(2) + (XO(3)-RG(3))*E3(3)
  
  
   RO    = sqrt (X*X+Y*Y+Z*Z)
   RATIO = RO/DIAG
   FIS   =0
   UL=0
   if (RATIO.gt.10) then
     AREA1 = AREA/PI4
     FIS   = - AREA1/RO
     return
   endif


   UL (1) = 0.d0
   UL (2) = 0.d0
   FIL    = 0.d0

   ZP =   1.d0
   if (Z.lt.0.0d0) ZP = -1.d0
   AZ = abs(Z)
   AZ = dmax1 (EPSS,AZ)
   Z  = AZ*ZP
   ZK    = Z*Z

   do K = 1,4
   XK    = (X-S(K)) * (X-S(K))
   YK    = (Y-T(K)) * (Y-T(K))
   H (K) = (X-S(K)) * (Y-T(K))
   R (K) = sqrt (XK + YK + ZK)
   E (K) = ZK + XK
   enddo

   do K = 1,4
      if (K.ne.NSIDE)  then
        K1    = K
        K2    = K+1
        if (K.eq.4) K2 = 1
        A2 = R(K1)+R(K2)+D(K)
        A1 = R(K1)+R(K2)-D(K)
        if (abs(A1*A2).lt.TINY) cycle
        AK = dlog(A1/A2)
        TOO  (K) = -(X-S(K1))*SINB(K) + (Y-T(K1))*COSB(K)
        UL (1) = UL (1) + (T(K2)-T(K1))*AK/D(K)
        UL (2) = UL (2) - (S(K2)-S(K1))*AK/D(K)
        FIL    = FIL - TOO(K)*AK
      endif
   enddo

   UL(3) = 0.5d0 * PI4
   if (ISING.ne.0) go to 4
   AZ = dmax1 (TINY,AZ)
   Z  = AZ*ZP
   
   UL(3) = 0.d0
   do K = 1,4
   if (K.ne.NSIDE.and.abs(COSB(K)).ge.TINY)  then
       K1 = K
       K2 = K+1
       if (K1.eq.4) K2 = 1
       UL(3) = UL(3)                                                 &
         + atan ( (SINB(K)*E(K1)-COSB(K)*H(K1))/(Z*COSB(K)*R(K1)) )  &
         - atan ( (SINB(K)*E(K2)-COSB(K)*H(K2))/(Z*COSB(K)*R(K2)) )
   endif
   enddo   

 4 FIS    = FIS   +(FIL + Z*UL(3) ) / PI4

   return
   END SUBROUTINE FSOUR_A4

!---- Subroutine :GPANEL ----------------------------------------------


!        Data :  XCORN(3,i)     global co-ordinates of the vertices
!                IYNCP          if =1 then the c.p. is calculated
!       Output:  SK(4),TK(4)    local co-ordinates of the vertices
!                R(3)           g. c. of the center of the local system
!                EX,EY,EZ(3)    the local base
!                DIAG           diameter of the element
!                AREA           area
!                XCP(3)         g. c. of the control point
!                AMOMX,AMOMY    cartesian 1st order moments
!                AMOMXX,AMOMXY  cartesian 2nd order moments
!                       AMOMYY
!                D(4)           lengths of the element's sides
!                CBA(4), SBA(4) cosines and sines of the element's sides
!                               with respect to the local system
!                NCORN          =0 then the element is a quadrilateral
!                                k then the k-th side of the element is
!                                  of zero length
!
!

 Subroutine  GPANEL( XCORN , R     , EX     , EY     , EZ     ,&
                     SK    , TK    , DIAG   , AREA   , XCP    ,&
                     D     , CBA   , SBA    ,&
                     IYNCP , NCORN , IERR                      )


   implicit none
      
   real(8),intent(in)  :: XCORN(3,4)

   integer,intent(in)  :: IYNCP
     
   real(8),intent(out) :: SK(4),TK(4),R(3),EX(3),EY(3),EZ(3),               &
                          DIAG,AREA,XCP(3),&
                          D(4),CBA(4),SBA(4)
     
   integer,intent(out)::NCORN,IERR

   real(8) :: XM(3),RR(4) , RX (4) , RY (4) , DD (4),A(4),XCORN1(3,4),  &
              DEX,DEY,DN,DE1,SO,TO,    &
              SCP,TCP,U,V,Ux,Uy,Vx,Vy,EPS,DET,  &
              DETS,DETT,DSCP,DTCP,TINY,TINYs

   integer :: J,K,NOD,K1,K2,ITER,L
     TINY=1d-014
     TINYs=1d-014

!-- Define the local coordinate base     EX(), EY(), EZ()
!          the element's center          XM()
!      and the maximum diagonal          DIAG
      IERR = 0
      EX (1:3) = XCORN(1:3,3) - XCORN(1:3,1)
      EY (1:3) = XCORN(1:3,4) - XCORN(1:3,2)
      XM (1:3) = 0.25d0*( XCORN(1:3,1)+XCORN(1:3,2)+XCORN(1:3,3)+XCORN(1:3,4) )
      DEX   =  dsqrt  ( EX(1)*EX(1)+EX(2)*EX(2)+EX(3)*EX(3) )
      DEY   =  dsqrt  ( EY(1)*EY(1)+EY(2)*EY(2)+EY(3)*EY(3) )
      DIAG  =  dmax1  ( DEX , DEY )
      EZ(1) = -EX(2)*EY(3) + EX(3)*EY(2)
      EZ(2) = -EX(3)*EY(1) + EX(1)*EY(3)
      EZ(3) = -EX(1)*EY(2) + EX(2)*EY(1)
      DN    =  dsqrt  ( dot_product (EZ,EZ ) )
   if (DN <= TINYs ) then
      write (*,*) 'GPANEL:NCORN=-1, DN=',DN
     write(*,*) XCORN(1,1), XCORN(2,1), XCORN(3,1)
     write(*,*) XCORN(1,2), XCORN(2,2), XCORN(3,2)
     write(*,*) XCORN(1,3), XCORN(2,3), XCORN(3,3)
     write(*,*) XCORN(1,4), XCORN(2,4), XCORN(3,4)
     NCORN = -1
     return
   endif
      EZ(1:3) =  EZ(1:3)/DN
      EY(1)   = -EX(2)*EZ(3) + EX(3)*EZ(2)
      EY(2)   = -EX(3)*EZ(1) + EX(1)*EZ(3)
      EY(3)   = -EX(1)*EZ(2) + EX(2)*EZ(1)
      DEY     = dsqrt ( dot_product (EY,EY) )
      EX(1:3) = EX(1:3) / DEX
      EY(1:3) = EY(1:3) / DEY


!-- Define the plane coordinates    XCORN1(3,.)
      DE1 = dot_product( EZ(1:3), XM(1:3)-XCORN(1:3,1) )
      do  NOD = 1, 4
         XCORN1 (1:3,NOD) = XCORN(1:3,NOD) + EZ(1:3)* ( (-1.d0)**(NOD-1) ) * DE1
      enddo ! NOD


!-- Define the local coordinates    SK() , TK()
   do J = 1, 4
      SK (J) = dot_product( EX(1:3),XCORN1(1:3,J)-XM(1:3))
      TK (J) = dot_product( EY(1:3),XCORN1(1:3,J)-XM(1:3))
   enddo
      SO = ((TK(1)-TK(2))*SK(4)+(TK(4)-TK(1))*SK(2))/(3.d0*(TK(2)-TK(4)))
      TO = -TK(1)/3.d0

   SK(1:4) = SK(1:4)-SO
   TK(1:4) = TK(1:4)-TO

   R( 1:3) = EX(1:3)*SO + EY(1:3)*TO + XM(1:3)

      NCORN = 0
   do K = 1, 4
      L=K+1
      if (K == 4) L=1
      D   (K) = dsqrt ( (SK(L)-SK(K))**2 + (TK(L)-TK(K))**2 )
      if ( D(K).gt.TINY  ) then
        SBA (K) = (TK(L)-TK(K))/D(K)
        CBA (K) = (SK(L)-SK(K))/D(K)
      else
        SBA (K) =0.d0
        CBA (K) =0.d0
        NCORN   = K
      endif
   enddo


!-- Calculate the element's control point
      if (IYNCP == 0) then
         XCP (1:3) = ( XCORN(1:3,1)+XCORN(1:3,2)+XCORN(1:3,3)+XCORN(1:3,4) ) / 4.d0
         goto 16
      endif

      SCP  = 0.d0
      TCP  = 0.d0
      EPS  = 0.0001d0
      ITER = 0
 10   ITER = ITER+1
   do K=1,4
      RR (K) = dsqrt ( (SK(K)-SCP)**2+(TK(K)-TCP)**2)
   enddo
      U=0.d0
      V=0.d0
      UX=0.d0
      UY=0.d0
      VX=0.d0
      VY=0.d0

   do K=1,4
      L=K+1
      if (K == 4) L=1
      if (NCORN == K) cycle
      RX(K)=(SCP-SK(K))/RR(K)+(SCP-SK(L))/RR(L)
      RY(K)=(TCP-TK(K))/RR(K)+(TCP-TK(L))/RR(L)
      DD(K)=((RR(K)+RR(L))*(RR(K)+RR(L))-D(K)**2)/2.d0
      A(K)=(RR(K)+RR(L)-D(K))/(RR(K)+RR(L)+D(K))
      A(K)=dlog(A(K))
      U=U+(TK(L)-TK(K))*A(K)/D(K)
      V=V+(SK(K)-SK(L))*A(K)/D(K)
      UX=UX+(TK(L)-TK(K))*RX(K)/DD(K)
      UY=UY+(TK(L)-TK(K))*RY(K)/DD(K)
      VX=VX+(SK(K)-SK(L))*RX(K)/DD(K)
      VY=VY+(SK(K)-SK(L))*RY(K)/DD(K)
   enddo

      DET=UX*VY-UY*VX
      DETS=UY*V-U*VY
      DETT=U*VX-UX*V
      DSCP=DETS/DET
      DTCP=DETT/DET
      SCP=SCP+DSCP
      TCP=TCP+DTCP
      if (ITER >= 10) goto 14
      if ( (DSCP < EPS).and.(DTCP < EPS) ) goto 13
      goto 10
! the maximum number of iterations has been exceeded
 14   IERR= 1
      SCP=0.d0
      TCP=0.d0
! convergence accomplished
 13   XCP(1:3)=R(1:3)+EX(1:3)*SCP+EY(1:3)*TCP

 16   AREA =0.5d0*(SK(3)-SK(1))*(TK(2)-TK(4))


 END Subroutine GPANEL
!--------------------------------------------------------------------

