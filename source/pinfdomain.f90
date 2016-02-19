!--------------------------------------------------------------------------------
!>@file 
!>@brief Find the boundary conditions for the poisson solver using infinite domain assumption
!------------------------------------------------------------------------!
!->Subroutine infdomain                                                  !
!  This subroutine calculates Boundary conditions using the infinite     !
!  Domain Poisson Poblem.                                                !
! - First grad^2(phi) = 0 is solved(no boundary conditions)              !
! - From this solution th(phi)/thn is constructed at the boundaries      !
!   which gives the boundary conditions                                  !
!------------------------------------------------------------------------!
!----------------------------------------------------------------------------------!
Subroutine infdomain(neqs,neqf)
  use MPI
    Implicit None
    integer,intent(in) :: neqs,neqf
    integer            :: NXs,NXf,NYs,NYf,nn,ndum,neq_siz
    integer            :: my_rank,ierr
    double precision   :: et,st

     call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)

    !--Solve PM with zero boundary conditions
    neq_siz = neqf
    allocate(SOL_0_pm(neq_siz,NXpm,NYpm,NZpm))
    allocate(source_bound(neq_siz,Nblocks*(NXpm*2+NYpm*2)),x_s(2,Nblocks*(NXpm*2+NYpm*2)),y_s(2,Nblocks*(NXpm*2+NYpm*2)))
    allocate(d_s((Nblocks*(NXpm*2+NYpm*2))))

    !these variables are used in case we want a higher order source definition(linear for example)
    d_s   = 0.d0
    x_s   = 0.d0
    y_s   = 0.d0
    SOL_0_pm=0.d0
    !sources will reside 2 dummy cells away from the original domain
    source_bound=0
    !We will solve a poisson problem with zero bc's at the point of calculation of the sources
    NXs = 1     !NXs_bl(1) - ndum
    NXf = NXpm  !NXf_bl(1) + ndum
    NYs = 1     !NYs_bl(1) - ndum
    NYf = NYpm
    !-->Solve for zero boundary conditions
    do neq=neqs,neqf
        call solve_eq_0(NXs,NXf,NYs,NYf,neq)
    enddo
    !-->Calculate normal gradient (which defines the sources)
    nbound=0
    call calc_normalderiv(NXs,NXf,NYs,NYf,neqs,neqf)
    deallocate(SOL_0_pm)
    if(my_rank.eq.0) st = MPI_WTIME()
    if(itree.eq.1) then 
      call build_level_nbound(NXs,NXf,NYs,NYf,neqs,neqf)
    else if (itree.eq.2) then 
      call build_level_nbound_new(NXs,NXf,NYs,NYf,neqs,neqf)
    endif
    if(my_rank.eq.0) then 
        et=MPI_WTIME()
        write(199,*) 'tree struc',int((et-st)/60),'m',mod(et-st,60.d0),'s'
    endif
    !-->SOL_pm_0 not needed for anything else

    !-->Using the sources calculated above define the correct bc's to the extended domain(BS law)
    NXs =1        !NXs_bl(nb)
    NXf =NXpm     !NXf_bl(nb)
    NYs =1        !NYs_bl(nb)
    NYf =NYpm     !NYf_bl(nb)
    if(my_rank.eq.0) st = MPI_WTIME()
    if (itree.eq.1) then 
       call Bounds2d_lev(ibctyp,NXs,NXf,NYs,NYf,neqs,neqf)
       deallocate(source_bound_lev,ds_lev,xs_lev,ys_lev,nbound_lev)
       deallocate(source_bound,d_s,x_s,y_s)
    else if (itree.eq.2) then
       call Bounds2d_lev_new(ibctyp,NXs,NXf,NYs,NYf,neqs,neqf)
       deallocate(source_bound_lev,ds_lev,xs_lev,ys_lev,nbound_lev,ilev_t)
       deallocate(source_bound,d_s,x_s,y_s)
    else 
       call Bounds2d(ibctyp,NXs,NXf,NYs,NYf,neqs,neqf)
       deallocate(source_bound,d_s,x_s,y_s)
    endif
     if(my_rank.eq.0) et=MPI_WTIME()
     if(my_rank.eq.0) write(199,*) 'Bounds',int((et-st)/60),'m',mod(et-st,60.d0),'s'

End Subroutine infdomain

!------------------------------------------------------------------------!
!->Subroutine infdomain                                                  !
!  This subroutine calculates Boundary conditions using the infinite     !
!  Domain Poisson Poblem.                                                !
! - First grad^2(phi) = 0 is solved(no boundary conditions)              !
! - From this solution th(phi)/thn is constructed at the boundaries      !
!   which gives the boundary conditions                                  !
!------------------------------------------------------------------------!
Subroutine infdomain_3d(neqs,neqf)
  use MPI
    Implicit None
    integer,intent(in) :: neqs,neqf
    integer            :: NXs,NXf,NYs,NYf,NZs,NZf,nn,ndum,neq,nworkb,neq_siz

    integer            :: my_rank,ierr,n,ibctypn,ip,jp,kp
    double precision   :: et,st,DDX,DDY,DDZ
    double precision,allocatable :: xs_lev0(:,:,:),ys_lev0(:,:,:),zs_lev0(:,:,:)

    call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
    nworkb = 2 * NXpm *  NYpm + 2 * NXpm * NZpm + 2 * NZPm* Nypm
    neq_siz = neqf
    allocate(SOL_0_pm(neq_siz,NXpm,NYpm,NZpm))
    allocate(source_bound(neq_siz,Nblocks*(nworkb)),x_s(4,Nblocks*(nworkb)),y_s(4,Nblocks*(nworkb)))
    allocate(z_s(4,Nblocks*(nworkb)))

    allocate(d_s((Nblocks*(nworkb))))
    d_s   = 0.d0
    x_s   = 0.d0
    y_s   = 0.d0
    z_s   = 0.d0
    SOL_0_pm=0.d0
    !sources will reside 2 dummy cells away from the original domain
    ndum = 2
    ndum = 2**levmax
    nb   = 1

    !We will solve a poisson problem with zero bc's at the point of calculation of the sources
    NXs = 1         !NXs_bl(nb) - ndum
    NXf = NXpm      !NXf_bl(nb) + ndum
    NYs = 1         !NYs_bl(nb) - ndum
    NYf = NYpm      !NYf_bl(nb) + ndum
    NZs = 1         !NZs_bl(nb) - ndum
    NZf = NZpm      !NZf_bl(nb) + ndum
    !-->Solve for zero boundary conditions
    do neq = neqs,neqf
        call solve_eq_0_3d(NXs,NXf,NYs,NYf,NZs,NZf,neq)
    enddo
    !-->Calculate normal gradient
    nbound=0
    call calc_normalderiv_3d(NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    deallocate(SOL_0_pm)
    if(my_rank.eq.0) st=MPI_WTIME()
    if(itree.eq.1) then
       call build_level_nbound_3d(NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    else if(itree.eq.2) then 
       call build_level_nbound_3d_new(NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
    endif
    if(my_rank.eq.0) then 
        et=MPI_WTIME()
        write(199,*) 'tree struc',int((et-st)/60),'m',mod(et-st,60.d0),'s'
    endif
    !-->SOL_pm_0 not needed for anything else

    !-->Using the sources calculated above define the correct bc's to the extended domain(BS law)
    NXs = 1     !NXs_bl(1)
    NXf = NXpm  !NXf_bl(1)
    NYs = 1     !NYs_bl(1)
    NYf = NYpm  !NYf_bl(1)
    NZs = 1     !NZs_bl(1)
    NZf = NZpm  !NZf_bl(1)
    if(my_rank.eq.0) st=MPI_WTIME()
    if (itree.eq.1) then 
       call Bounds3d_lev(ibctyp,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
       deallocate(source_bound_lev,ds_lev,xs_lev,ys_lev,zs_lev,nbound_lev)
       deallocate(source_bound,d_s,x_s,y_s,z_s)
    else if (itree.eq.2) then 
       call Bounds3d_lev_new(ibctyp,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
       deallocate(source_bound_lev,ds_lev,xs_lev,ys_lev,zs_lev,nbound_lev,ilev_t)
       deallocate(source_bound,d_s,x_s,y_s,z_s)
    else 
       call Bounds3d(ibctyp,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
       deallocate(source_bound,d_s,x_s,y_s,z_s)
    endif
    if(my_rank.eq.0) then 
        et=MPI_WTIME()
        write(199,*) 'Bounds',int((et-st)/60),'m',mod(et-st,60.d0),'s'
    endif
 
End Subroutine infdomain_3d
!-------------------------------------------------------------------------!
!-> Subroutine calc_normalderiv                                           !
!   This Subroutine calculated normal derivate at the boundary of the     !
!   domain.This together with the appropriate Green function calculate    !
!   the boundary conditions at the domain boundaries.                     !
!   For the calculation of the derivative  a fourth order one sided       !
!   difference approximation is used:                                     !
!   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
!   Note that one [X,Y,Z]MIN boundaries its is a forwardd difference (and !
!   thus with the difference scheme is with a minus but because the       !
!   the sources definitions points outward the minus sign is cancelled    !
!   and so the scheme is the same for MIN and MAX  boundaries             !
!-------------------------------------------------------------------------!
Subroutine calc_normalderiv(NXs, NXf, NYs, NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf,neqs,neqf
    integer           :: i, j,j1,i1,neq
    double precision  :: a1, a2, a3, a4, a5,psi1,psi2


    !Sources are defined in NXs,NXf,Nys,NYf using SOL_pm0(zero bc solution)
    a1 =  25.d0/12.d0
    a2 = -4.d0
    a3 =  3.d0
    a4 = -4.d0/3.d0
    a5 =  0.25d0
    i = Nxs
    do j = NYs,NYf-1
        nbound = nbound + 1
        j1 = j + 1


        do neq = neqs,neqf
            psi1 =  1./DXpm * (a1 * SOL_0_pm (neq,i,j,1)  + a2 * SOL_0_pm (neq,i+1,j,1)  + a3 * SOL_0_pm (neq,i+2,j,1)  + &
                a4 * SOL_0_pm (neq,i+3,j,1)  + a5 * SOL_0_pm (neq,i+4,j,1) )
            psi2 =  1./DXpm * (a1 * SOL_0_pm (neq,i,j1,1)  + a2 * SOL_0_pm (neq,i+1,j1,1)  + a3 * SOL_0_pm (neq,i+2,j1,1)  + &
                a4 * SOL_0_pm (neq,i+3,j1,1)  + a5 * SOL_0_pm (neq,i+4,j1,1) )
            source_bound(neq,nbound) =  0.5d0 * (psi1 + psi2)
        enddo

        x_s(1,nbound)      = XMIN_pm + (i - 1) * DXpm
        x_s(2,nbound)      = XMIN_pm + (i - 1) * DXpm
        y_s(1,nbound)      = YMIN_pm + (j - 1) * DYpm!+ 0.5d0 * DYpm
        y_s(2,nbound)      = YMIN_pm + (j)     * DYpm!+ 0.5d0 * DYpm
        d_s(nbound)        = DYpm

    enddo

    !---XMAX BOUNDARY----
    i = NXf
    do j = NYs,NYf-1
        nbound = nbound+ 1
        j1 = j + 1


        do neq = neqs,neqf
            psi1  = 1./DXpm * (a1 * SOL_0_pm (neq,i,j,1)  + a2 * SOL_0_pm (neq,i-1,j,1)  + a3 * SOL_0_pm (neq,i-2,j,1)  + &
                a4 * SOL_0_pm (neq,i-3,j,1)  + a5 * SOL_0_pm (neq,i-4,j,1) )
            psi2  = 1./DXpm * (a1 * SOL_0_pm (neq,i,j1,1)  + a2 * SOL_0_pm (neq,i-1,j1,1)  + a3 * SOL_0_pm (neq,i-2,j1,1)  +&
                a4 * SOL_0_pm (neq,i-3,j1,1)  + a5 * SOL_0_pm (neq,i-4,j1,1) )
            source_bound(neq,nbound) =  0.5d0 * (psi1 + psi2)
        enddo

        x_s(1,nbound)    = XMIN_pm + (i - 1) * DXpm
        x_s(2,nbound)    = XMIN_pm + (i - 1) * DXpm
        y_s(1,nbound)    = YMIN_pm + (j)  * DYpm
        y_s(2,nbound)    = YMIN_pm + (j - 1) * DYpm
        d_s(nbound)      = DYpm

    enddo

    j1=0

    !---YMIN BOUNDARY----(fn points at the y direction)
    j = NYs
    do i = NXs,NXf-1
        nbound = nbound + 1
        i1 = i + 1


        do neq = neqs,neqf
            psi1 =  1./DYpm * (a1 * SOL_0_pm (neq,i,j,1)  + a2 * SOL_0_pm (neq,i,j+1,1)  + a3 * SOL_0_pm (neq,i,j+2,1)  +&
                a4 * SOL_0_pm (neq,i,j+3,1)  + a5 * SOL_0_pm (neq,i,j+4,1) )
            psi2 =  1./DYpm * (a1 * SOL_0_pm (neq,i1,j,1)  + a2 * SOL_0_pm (neq,i1,j+1,1)  + a3 * SOL_0_pm (neq,i1,j+2,1)  + &
                a4 * SOL_0_pm (neq,i1,j+3,1)  + a5 * SOL_0_pm (neq,i1,j+4,1) )
            source_bound(neq,nbound) = 0.5d0 * (psi1 + psi2)
        enddo


        x_s(1,nbound)    = XMIN_pm + (i ) * DXpm!+ 0.5d0 * DXpm
        x_s(2,nbound)    = XMIN_pm + (i - 1) * DXpm!+ 0.5d0 * DXpm
        y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
        y_s(2,nbound)    = YMIN_pm + (j - 1) * DYpm
        d_s(nbound)      = DXpm

    enddo

    !---YMAX BOUNDARY----
    j = NYf
    do i = NXs,NXf-1
        nbound = nbound + 1
        i1 = i + 1


        do neq = neqs,neqf
            psi1 = 1./DYpm * (a1 * SOL_0_pm (neq,i,j,1)  + a2 * SOL_0_pm (neq,i,j-1,1)  + a3 * SOL_0_pm (neq,i,j-2,1)  + &
                a4 * SOL_0_pm (neq,i,j-3,1)  + a5 * SOL_0_pm (neq,i,j-4,1) )
            psi2 = 1./DYpm * (a1 * SOL_0_pm (neq,i1,j,1)  + a2 * SOL_0_pm (neq,i1,j-1,1)  + a3 * SOL_0_pm (neq,i1,j-2,1)  + &
                a4 * SOL_0_pm (neq,i1,j-3,1)  + a5 * SOL_0_pm (neq,i1,j-4,1) )
            source_bound(neq,nbound) = 0.5d0 * (psi1 + psi2)
        enddo

        x_s(1,nbound)    = XMIN_pm + (i - 1) * DXpm! + 0.5d0 * DXpm
        x_s(2,nbound)    = XMIN_pm + (i) * DXpm! + 0.5d0 * DXpm
        y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
        y_s(2,nbound)    = YMIN_pm + (j - 1) * DYpm
        d_s(nbound)      = DXpm

    enddo

End Subroutine calc_normalderiv

!-------------------------------------------------------------------------!
!-> Subroutine calc_normalderiv                                           !
!   This Subroutine calculated normal derivate at the boundary of the     !
!   domain.This together with the appropriate Green function calculate    !
!   the boundary conditions at the domain boundaries.                     !
!   For the calculation of the derivative  a fourth order one sided       !
!   difference approximation is used:                                     !
!   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
!   Note that one [X,Y,Z]MIN boundaries its is a forwardd difference (and !
!   thus with the difference scheme is with a minus but because the       !
!   the sources definitions points outward the minus sign is cancelled    !
!   and so the scheme is the same for MIN and MAX  boundaries             !
!-------------------------------------------------------------------------!
Subroutine calc_normalderiv_3d(NXs, NXf, NYs, NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer           :: i, j, k, j1, i1, k1,neq
    double precision  :: a1, a2, a3, a4, a5,psi1,psi2,psi3,psi4

    a1 =  25.d0/12.d0
    a2 = -4.d0
    a3 =  3.d0
    a4 = -4.d0/3.d0
    a5 =  0.25d0
    !Sources are defined AT CELL CENTERS
    !---XMIN BOUNDARY----(fn points at the x direction)
    i = Nxs
    do k = NZs, NZf - 1
        do j = NYs,NYf - 1
            nbound = nbound + 1
            j1 = j + 1
            k1 = k + 1

            do neq = neqs,neqf
                psi1 =  1./DXpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm (neq,i+1,j,k)   + a3 * SOL_0_pm (neq,i+2,j,k)   + &
                    a4 * SOL_0_pm (neq,i+3,j,k)   + a5 * SOL_0_pm (neq,i+4,j,k))
                psi2 =  1./DXpm * (a1 * SOL_0_pm(neq,i,j1,k)  + a2 * SOL_0_pm(neq,i+1,j1,k)  + a3 * SOL_0_pm (neq,i+2,j1,k)  + &
                    a4 * SOL_0_pm (neq,i+3,j1,k)  + a5 * SOL_0_pm (neq,i+4,j1,k))
                psi3 =  1./DXpm * (a1 * SOL_0_pm(neq,i,j1,k1) + a2 * SOL_0_pm(neq,i+1,j1,k1) + a3 * SOL_0_pm (neq,i+2,j1,k1) + &
                    a4 * SOL_0_pm (neq,i+3,j1,k1) + a5 * SOL_0_pm (neq,i+4,j1,k1))
                psi4 =  1./DXpm * (a1 * SOL_0_pm(neq,i,j,k1)  + a2 * SOL_0_pm(neq,i+1,j,k1)  + a3 * SOL_0_pm (neq,i+2,j,k1)  + &
                    a4 * SOL_0_pm (neq,i+3,j,k1)  + a5 * SOL_0_pm (neq,i+4,j,k1))
                source_bound(neq,nbound) =  0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(1,nbound)      = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)      = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)      = XMIN_pm + (i - 1) * DXpm
            x_s(4,nbound)      = XMIN_pm + (i - 1) * DXpm
                                                         
            y_s(1,nbound)      = YMIN_pm + (j-1)   * DYpm
            y_s(2,nbound)      = YMIN_pm + (j-1)   * DYpm
            y_s(3,nbound)      = YMIN_pm + (j)     * DYpm
            y_s(4,nbound)      = YMIN_pm + (j)     * DYpm
                                                         
            z_s(1,nbound)      = ZMIN_pm + (k-1)   * DZpm
            z_s(2,nbound)      = ZMIN_pm + (k)     * DZpm
            z_s(3,nbound)      = ZMIN_pm + (k)     * DZpm
            z_s(4,nbound)      = ZMIN_pm + (k-1)   * DZpm

            d_s(nbound)        = DYpm * DZpm

        enddo
    enddo
    !---XMAX BOUNDARY----
    i = NXf
    do k = NZs, NZf - 1
        do j = NYs,NYf -1
            nbound = nbound+ 1
            j1 = j + 1
            k1 = k + 1

            do neq=neqs,neqf
                psi1  = 1./DXpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm(neq,i-1,j,k)   + a3 * SOL_0_pm(neq,i-2,j,k)   + &
                    a4 * SOL_0_pm(neq,i-3,j,k)   + a5 * SOL_0_pm(neq,i-4,j,k) )
                psi2  = 1./DXpm * (a1 * SOL_0_pm(neq,i,j1,k)  + a2 * SOL_0_pm(neq,i-1,j1,k)  + a3 * SOL_0_pm(neq,i-2,j1,k)  + &
                    a4 * SOL_0_pm(neq,i-3,j1,k)  + a5 * SOL_0_pm(neq,i-4,j1,k) )
                psi3  = 1./DXpm * (a1 * SOL_0_pm(neq,i,j1,k1) + a2 * SOL_0_pm(neq,i-1,j1,k1) + a3 * SOL_0_pm(neq,i-2,j1,k1) + &
                    a4 * SOL_0_pm(neq,i-3,j1,k1) + a5 * SOL_0_pm(neq,i-4,j1,k1) )
                psi4  = 1./DXpm * (a1 * SOL_0_pm(neq,i,j,k1)  + a2 * SOL_0_pm(neq,i-1,j,k1)  + a3 * SOL_0_pm(neq,i-2,j,k1)  + &
                    a4 * SOL_0_pm(neq,i-3,j,k1)  + a5 * SOL_0_pm(neq,i-4,j,k1) )
                source_bound(neq,nbound) =  0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(4,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(1,nbound)    = XMIN_pm + (i - 1) * DXpm
                                                          
            y_s(4,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(3,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(2,nbound)    = YMIN_pm + (j)     * DYpm
            y_s(1,nbound)    = YMIN_pm + (j)     * DYpm
                                                          
            z_s(4,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(3,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(2,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(1,nbound)    = ZMIN_pm + (k - 1) * DZpm

            d_s(nbound)      = DYpm * DZpm

        enddo
    enddo
    j1=0

    !---YMIN BOUNDARY----(fn points at the y direction)
    j = NYs
    do k = NZs,NZf - 1
        do i = NXs,NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            k1 = k + 1

            do neq=neqs,neqf
                psi1 =  1./DYpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm(neq,i,j+1,k)   + a3 * SOL_0_pm(neq,i,j+2,k)   + &
                    a4 * SOL_0_pm(neq,i,j+3,k)   + a5 * SOL_0_pm(neq,i,j+4,k))
                psi2 =  1./DYpm * (a1 * SOL_0_pm(neq,i1,j,k)  + a2 * SOL_0_pm(neq,i1,j+1,k)  + a3 * SOL_0_pm(neq,i1,j+2,k)  + &
                    a4 * SOL_0_pm(neq,i1,j+3,k)  + a5 * SOL_0_pm(neq,i1,j+4,k))
                psi3 =  1./DYpm * (a1 * SOL_0_pm(neq,i1,j,k1) + a2 * SOL_0_pm(neq,i1,j+1,k1) + a3 * SOL_0_pm(neq,i1,j+2,k1) + &
                    a4 * SOL_0_pm(neq,i1,j+3,k1) + a5 * SOL_0_pm(neq,i1,j+4,k1))
                psi4 =  1./DYpm * (a1 * SOL_0_pm(neq,i,j,k1)  + a2 * SOL_0_pm(neq,i,j+1,k1)  + a3 * SOL_0_pm(neq,i,j+2,k1)  + &
                    a4 * SOL_0_pm(neq,i,j+3,k1)  + a5 * SOL_0_pm(neq,i,j+4,k1))
                source_bound(neq,nbound) = 0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(1,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)    = XMIN_pm + (i)     * DXpm
            x_s(4,nbound)    = XMIN_pm + (i)     * DXpm
                                                                         
            y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(2,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(3,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(4,nbound)    = YMIN_pm + (j - 1) * DYpm
                                                                         
            z_s(1,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(2,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(3,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(4,nbound)    = ZMIN_pm + (k - 1) * DZpm

            d_s(nbound)      = DXpm * DZpm

        enddo
    enddo

    !---YMAX BOUNDARY----
    j = NYf
    do k = NZs,NZf - 1
        do i = NXs,NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            k1 = k + 1

            do neq =neqs,neqf
                psi1 = 1./DYpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm(neq,i,j-1,k)   + a3 * SOL_0_pm(neq,i,j-2,k)   + &
                    a4 * SOL_0_pm (neq,i,j-3,k)   + a5 * SOL_0_pm (neq,i,j-4,k))
                psi2 = 1./DYpm * (a1 * SOL_0_pm(neq,i1,j,k)  + a2 * SOL_0_pm(neq,i1,j-1,k)  + a3 * SOL_0_pm(neq,i1,j-2,k)  + &
                    a4 * SOL_0_pm (neq,i1,j-3,k)  + a5 * SOL_0_pm (neq,i1,j-4,k))
                psi3 = 1./DYpm * (a1 * SOL_0_pm(neq,i1,j,k1) + a2 * SOL_0_pm(neq,i1,j-1,k1) + a3 * SOL_0_pm(neq,i1,j-2,k1) + &
                    a4 * SOL_0_pm (neq,i1,j-3,k1) + a5 * SOL_0_pm (neq,i1,j-4,k1))
                psi4 = 1./DYpm * (a1 * SOL_0_pm(neq,i,j,k1)  + a2 * SOL_0_pm(neq,i,j-1,k1)  + a3 * SOL_0_pm(neq,i,j-2,k1)  + &
                    a4 * SOL_0_pm (neq,i,j-3,k1)  + a5 * SOL_0_pm (neq,i,j-4,k1))
                source_bound(neq,nbound) = 0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(4,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)    = XMIN_pm + (i)     * DXpm
            x_s(1,nbound)    = XMIN_pm + (i)     * DXpm
                                                                         
            y_s(4,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(3,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(2,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
                                                                         
            z_s(4,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(3,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(2,nbound)    = ZMIN_pm + (k)     * DZpm
            z_s(1,nbound)    = ZMIN_pm + (k - 1) * DZpm

            d_s(nbound)      = DXpm * DZpm

        enddo
    enddo
    !---ZMIN BOUNDARY----
    k1=0
    k = NZs
    do j = NYs,NYf - 1
        do i = NXs,NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            j1 = j + 1

            do neq =neqs,neqf
                psi1 = 1./DZpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm(neq,i,j,k+1)   + a3 * SOL_0_pm(neq,i,j,k+2)   + &
                    a4 * SOL_0_pm (neq,i,j,k+3)   + a5 * SOL_0_pm (neq,i,j,k+4))
                psi2 = 1./DZpm * (a1 * SOL_0_pm(neq,i1,j,k)  + a2 * SOL_0_pm(neq,i1,j,k+1)  + a3 * SOL_0_pm(neq,i1,j,k+2)  + &
                    a4 * SOL_0_pm (neq,i1,j,k+3)  + a5 * SOL_0_pm (neq,i1,j,k+4))
                psi3 = 1./DZpm * (a1 * SOL_0_pm(neq,i1,j1,k) + a2 * SOL_0_pm(neq,i1,j1,k+1) + a3 * SOL_0_pm(neq,i1,j1,k+2) + &
                    a4 * SOL_0_pm (neq,i1,j1,k+3) + a5 * SOL_0_pm (neq,i1,j1,k+4))
                psi4 = 1./DZpm * (a1 * SOL_0_pm(neq,i,j1,k)  + a2 * SOL_0_pm(neq,i,j1,k+1)  + a3 * SOL_0_pm(neq,i,j1,k+2)  + &
                    a4 * SOL_0_pm (neq,i,j1,k+3)  + a5 * SOL_0_pm (neq,i,j1,k+4))
                source_bound(neq,nbound) = 0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(1,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)    = XMIN_pm + (i)     * DXpm
            x_s(4,nbound)    = XMIN_pm + (i)     * DXpm
                                                                                     
            y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(2,nbound)    = YMIN_pm + (j)     * DYpm
            y_s(3,nbound)    = YMIN_pm + (j)     * DYpm
            y_s(4,nbound)    = YMIN_pm + (j - 1) * DYpm
                                                                                     
            z_s(1,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(2,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(3,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(4,nbound)    = ZMIN_pm + (k - 1) * DZpm

            d_s(nbound)      = DXpm * DYpm

            x_s(1,nbound)      = XMIN_pm + (i - 1) * DXpm
            y_s(1,nbound)      = YMIN_pm + (j - 1) * DYpm!+ 0.5d0 * DYpm
            z_s(1,nbound)      = ZMIN_pm + (k - 1) * DZpm!+ 0.5d0 * DYpm
        enddo
    enddo
    !---ZMAX BOUNDARY
    k = NZf
    do j = NYs,NYf - 1
        do i = NXs,NXf - 1
            nbound = nbound + 1
            i1 = i + 1
            j1 = j + 1

            do neq =neqs,neqf
                psi1 = 1./DZpm * (a1 * SOL_0_pm(neq,i,j,k)   + a2 * SOL_0_pm(neq,i,j,k-1)   + a3 * SOL_0_pm(neq,i,j,k-2)   + &
                    a4 * SOL_0_pm (neq,i,j,k-3)   + a5 * SOL_0_pm (neq,i,j,k-4))
                psi2 = 1./DZpm * (a1 * SOL_0_pm(neq,i1,j,k)  + a2 * SOL_0_pm(neq,i1,j,k-1)  + a3 * SOL_0_pm(neq,i1,j,k-2)  + &
                    a4 * SOL_0_pm (neq,i1,j,k-3)  + a5 * SOL_0_pm (neq,i1,j,k-4))
                psi3 = 1./DZpm * (a1 * SOL_0_pm(neq,i1,j1,k) + a2 * SOL_0_pm(neq,i1,j1,k-1) + a3 * SOL_0_pm(neq,i1,j1,k-2) + &
                    a4 * SOL_0_pm (neq,i1,j1,k-3) + a5 * SOL_0_pm (neq,i1,j1,k-4))
                psi4 = 1./DZpm * (a1 * SOL_0_pm(neq,i,j1,k)  + a2 * SOL_0_pm(neq,i,j1,k-1)  + a3 * SOL_0_pm(neq,i,j1,k-2)  + &
                    a4 * SOL_0_pm (neq,i,j1,k-3)  + a5 * SOL_0_pm (neq,i,j1,k-4))
                source_bound(neq,nbound) = 0.25d0 * (psi1 + psi2 + psi3 + psi4)
            enddo

            x_s(4,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(3,nbound)    = XMIN_pm + (i - 1) * DXpm
            x_s(2,nbound)    = XMIN_pm + (i)     * DXpm
            x_s(1,nbound)    = XMIN_pm + (i)     * DXpm
                                                                                     
            y_s(4,nbound)    = YMIN_pm + (j - 1) * DYpm
            y_s(3,nbound)    = YMIN_pm + (j)     * DYpm
            y_s(2,nbound)    = YMIN_pm + (j)     * DYpm
            y_s(1,nbound)    = YMIN_pm + (j - 1) * DYpm
                                                                                     
            z_s(4,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(3,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(2,nbound)    = ZMIN_pm + (k - 1) * DZpm
            z_s(1,nbound)    = ZMIN_pm + (k - 1) * DZpm

            d_s(nbound)      = DXpm * DYpm

        enddo
    enddo
End Subroutine calc_normalderiv_3d


!This Subroutine builds the nbounds_lev matrrix using the values calulated at the finer level.
Subroutine build_level_nbound(NXs,NXf,NYs,NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf,neqs,neqf
    integer          :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin, ires, leafacc
    double precision :: x, y, s, source(neqf), xc, yc, sc, sourcec(neqf)

    leafmax=1
    do lev=1,levmax
       leafmax= 2**lev + leafmax
    enddo
    allocate (source_bound_lev(nbound,neqs:neqf,0:levmax,leafmax)); source_bound_lev     = 0.d0
    allocate (xs_lev(nbound,0:levmax,leafmax)); xs_lev     = 0.d0
    allocate (ys_lev(nbound,0:levmax,leafmax)); ys_lev     = 0.d0
    allocate (ds_lev(nbound,0:levmax,leafmax)); ds_lev     = 0.d0  
    allocate (nbound_lev(0:levmax))             ; nbound_lev = 0
     
    nbound_lev(0)   = nbound
    xs_lev(:, levmax, 1) = 0.5d0 * (x_s(1,:) + x_s(2,:))
    ys_lev(:, levmax, 1) = 0.5d0 * (y_s(1,:) + y_s(2,:))
    ds_lev(:, levmax, 1) = d_s(:)
    do neq=neqs,neqf
       source_bound_lev(:,neq,0,1) = source_bound(neq,:)
    enddo

    
    icount   = 0
    leafcount= 0
    leafacc  = 0
    istep    = 2**levmax
    i = Nxs
    if (mod(NXf-NXs,istep).ne.0.or.mod(NYf-NYs,istep).ne.0) then 
       write(*,*) 'can divide to tree level'
       STOP
    endif
    do j = NYs,NYf-istep,istep
        icount = icount + 1
        !nleaf = 1 is the lev element info
        leaffin=leafmax
        do lev = levmax,0,-1
            leafcount=leafacc + (j-NYs+1)!to start from Nys leafcount + 1
            leafstart = leaffin - 2**lev + 1
               do nleaf = leaffin,leafstart,-1
                  ires = 2 ** (levmax -lev)
                  xc=0; yc=0; sc=0; sourcec=0
                  do k = 1,ires
                     x  = xs_lev(leafcount,levmax,1)
                     y  = ys_lev(leafcount,levmax,1)
                     s  = ds_lev(leafcount,levmax,1)
                     source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
                     
                     xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                     leafcount=leafcount +1
                  enddo
                  xs_lev(icount,0,nleaf) = xc/float(ires)
                  ys_lev(icount,0,nleaf) = yc/float(ires)
                  ds_lev(icount,0,nleaf) = sc
                  source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                  leaffin=leafstart-1
               enddo
         enddo
    enddo
    leafacc = (NYf-1)-(NYs) + 1 + leafacc !because in fine grid NYs,NYf-1
    !---XMAX BOUNDARY----
    i = NXf
    do j = NYs,NYf-istep,istep
        icount = icount + 1
        !nleaf = 1 is the lev element info
        leaffin=leafmax
        do lev = levmax,0,-1
            leafcount=leafacc + (j-NYs+1)!to start from Nys leafcount + 1
            leafstart = leaffin - 2**lev + 1
               do nleaf = leaffin,leafstart,-1
                  ires = 2 ** (levmax -lev)
                  xc=0; yc=0; sc=0; sourcec=0
                  do k = 1,ires
                     x  = xs_lev(leafcount,levmax,1)
                     y  = ys_lev(leafcount,levmax,1)
                     s  = ds_lev(leafcount,levmax,1)
                     source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
                     
                     xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                     leafcount=leafcount +1
                  enddo
                  xs_lev(icount,0,nleaf) = xc/float(ires)
                  ys_lev(icount,0,nleaf) = yc/float(ires)
                  ds_lev(icount,0,nleaf) = sc
                  source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                  leaffin=leafstart-1
               enddo
         enddo

    enddo
    leafacc = (NYf-1)-(NYs) + 1 + leafacc !because in fine grid NYs,NYf-1
   
   
    !---YMIN BOUNDARY----(fn points at the y direction)
    j = NYs
    do i = NXs,NXf-istep,istep
        icount = icount + 1
        !nleaf = 1 is the lev element info
        leaffin=leafmax
        do lev = levmax,0,-1
            leafcount=leafacc + (i-NXs+1)!to start from Nys leafcount + 1
            leafstart = leaffin - 2**lev + 1
               do nleaf = leaffin,leafstart,-1
                  ires = 2 ** (levmax -lev)
                  xc=0; yc=0; sc=0; sourcec=0
                  do k = 1,ires
                     x  = xs_lev(leafcount,levmax,1)
                     y  = ys_lev(leafcount,levmax,1)
                     s  = ds_lev(leafcount,levmax,1)
                     source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
                     
                     xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                     leafcount=leafcount +1
                  enddo
                  xs_lev(icount,0,nleaf) = xc/float(ires)
                  ys_lev(icount,0,nleaf) = yc/float(ires)
                  ds_lev(icount,0,nleaf) = sc
                  source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                  leaffin=leafstart-1
               enddo
         enddo
    enddo
   
    leafacc = (NXf-1)-(NXs) + 1 + leafacc !because in fine grid NYs,NYf-1
    !---YMAX BOUNDARY----
    j = NYf
    do i = NXs,NXf-istep,istep
        icount = icount + 1
        !nleaf = 1 is the lev element info
        leaffin=leafmax
        do lev = levmax,0,-1
            leafcount=leafacc + (i-NXs+1)!to start from Nys leafcount + 1
            leafstart = leaffin - 2**lev + 1
               do nleaf = leaffin,leafstart,-1
                  ires = 2 ** (levmax -lev)
                  xc=0; yc=0; sc=0; sourcec=0
                  do k = 1,ires
                     x  = xs_lev(leafcount,levmax,1)
                     y  = ys_lev(leafcount,levmax,1)
                     s  = ds_lev(leafcount,levmax,1)
                     source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
                     
                     xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                     leafcount=leafcount +1
                  enddo
                  xs_lev(icount,0,nleaf) = xc/float(ires)
                  ys_lev(icount,0,nleaf) = yc/float(ires)
                  ds_lev(icount,0,nleaf) = sc
                  source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                  leaffin=leafstart-1
               enddo
         enddo
    enddo

    leafacc = (NXf-1)-(NXs) + 1 + leafacc !because in fine grid NYs,NYf-1
    if(leafacc.ne.nbound) then 
       write(*,*) 'leafacc problem'
       STOP
    endif
    nbound_lev(0)=icount
  ! do i=1,nbound_lev(0)
  !    lev=0
  !    istep =2
  !    leafstart = 0
  !    do lev = 0,levmax
  !        open(14,file='lev'//char(48+lev),access='APPEND')
  !        leafstart = leafstart +2**max(0,(lev-1))
  !        leaffin   = leafstart +2**(lev) -1
  !        do j=leafstart,leaffin
  !           write(14,*) xs_lev(i,0,j),ys_lev(i,0,j)
  !        enddo
  !    enddo
  ! enddo
End Subroutine build_level_nbound



!-------------------------------------------------------------------------!
!-> Subroutine calc_normalderiv                                           !
!   This Subroutine calculated normal derivate at the boundary of the     !
!   domain.This together with the appropriate Green function calculate    !
!   the boundary conditions at the domain boundaries.                     !
!   For the calculation of the derivative  a fourth order one sided       !
!   difference approximation is used:                                     !
!   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
!-------------------------------------------------------------------------!
Subroutine build_level_nbound_3d(NXs, NXf, NYs, NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer          :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin, ires, leafacc
    integer          :: nleaflev, nleafroot
    integer          :: m, l, il, im, mm, nn, iresroot, NNX, NNY, NNZ, NNX0, NNY0, NNZ0
    double precision :: x, y, z, s, source(neqf), xc, yc, zc, sc, sourcec(neqf)
    double precision,allocatable :: xs_tmp(:,:),ds_tmp(:),s_tmp(:,:)


    leafmax=1
    do lev=1,levmax
       leafmax= 4**lev + leafmax
    enddo
    istep    = 2**levmax
    NNX0  = (NXf - istep -NXs)/istep  +1
    NNY0  = (NYf - istep -NYs)/istep  +1
    NNZ0  = (NZf - istep -NZs)/istep  +1
    allocate (nbound_lev(0:levmax))             ; nbound_lev = 0
    nbound_lev(0) = 2 * NNX0*NNY0 + 2 * NNY0*NNZ0+ 2 * NNX0*NNZ0
    allocate (source_bound_lev(nbound_lev(0),neqs:neqf,0:levmax,leafmax)); source_bound_lev     = 0.d0
    allocate (xs_lev(nbound_lev(0),0:1,leafmax)); xs_lev     = 0.d0
    allocate (ys_lev(nbound_lev(0),0:1,leafmax)); ys_lev     = 0.d0
    allocate (zs_lev(nbound_lev(0),0:1,leafmax)); ys_lev     = 0.d0
    allocate (ds_lev(nbound_lev(0),0:1,leafmax)); ds_lev     = 0.d0  
     
    allocate (xs_tmp(nbound,3),ds_tmp(nbound),s_tmp(nbound,1:neqf))
     xs_tmp(1:nbound,1) = 0.25d0 * (x_s(1,:) + x_s(2,:) + x_s(3,:) + x_s(4,:))
     xs_tmp(1:nbound,2) = 0.25d0 * (y_s(1,:) + y_s(2,:) + y_s(3,:) + y_s(4,:))
     xs_tmp(1:nbound,3) = 0.25d0 * (z_s(1,:) + z_s(2,:) + z_s(3,:) + z_s(4,:))
    ds_tmp(1:nbound) = d_s(:)
    do neq=neqs,neqf
       s_tmp(1:nbound,neq)  = source_bound(neq,1:nbound)
    enddo

   
!tree defined at cells
    NNZ = NZf  - NZs 
    NNY = NYf  - NYs 
    NNX = NXf  - NXs 

    icount   = 0
    leafcount= 0
    leafacc  = 0
    if (mod(NNX,istep).ne.0.or.mod(NNY,istep).ne.0.or.mod(NNZ,istep).ne.0) then 
       write(*,*) 'can divide to tree level,',NNX,NNY,NNZ
       STOP
    endif
!   !It might seem a  little bit complex but it isnt.
!   !tree in 2d surfaces means that the tree has 4 leaf's at each level
!   we start the tree building using the 0-th level  nodes(0th=coarser)
!   since we know the information at levmax(finer)  we loop the finer grid
!   using istep. istep is  2 because,2  in 2 in j make the 4 leafs
!   the stucture is
!   0-th node -->nleaf
!   n leaf is the accumulative number of leaf of the 0th node
!      1
!     /\
!    2  3
!   /\  /\
!  4 5  6 7
! 
!
! so we know that nleaf 4 of 0th node 4 refers at the second level(it's easy to find out that)
!
! the structure is build from levmax to  0th level
! we find the first leaf  for each  level and and we add up one for each leaf on the level until we finish.
! because the structure on the finer level is known in a stuctured 2d grid fashion
! in order to assign the correct nodes to the leads we need to loop through the grid.
! Also we know that each level has members both in i and j direction 
! that's we we want to scan using 2d grid patches
! 
! ires is the resolution of each level.so finer 
!      level has resolution of 1(each member has is the finest you can have)
! 
! for level  levmax -1 ires = 4 in 3d or ires 2 in 2d
!  
! in 3d we split the ires in 2 direction again because of the structured grid fine data we have.
! we need htat splitting beacuse the data for each level is reconstructed from the original fine data
! So for level levmax -1 we need to cycle through the 4 elements(in 3d) ,2 in idirection 2 in j direction
! for level levmax -2 we need to cycle through 16 elements.
!
! Why didn't i use the 4 elements of the previous level?I wanted to use the original data.even thought because
! of the linear interpolation of the data it could be the same using the 4 elements of lev -1
!
! source(lev-1) = Sum(source(levmax)*Area(levmax))/Sum(Area(levmax))
!
! Because it's surface is a 2d grid with it's own numbering. we use leafacc to know the accumulated number 
! when changing sides (we map x_s(nod) to a specific leaf)
!
!
!
! Don't be afraid of  this :leafcount = leafacc + (k- NZs -1 + im  + (m-1) ) * NNY + (j -Nys + il +(l-1))
! node = (j-1)*NY + k
! k,j give the global index : think k - NZs -1 as j-1
!                             think j - NYs    as i
! The global index cycles at 0th level steps(coarser) so
!                            im,il give the starting node for each level
!                            (m-1),(l-1) give the node of the elemmenet for each level(-1 beacuse im,il have one)

    !---XMIN BOUNDARY----(fn points at the x direction)
    i  = Nxs
    do k = NZs, NZf - istep,istep
        do j = NYs,NYf - istep,istep
            icount = icount + 1
            !nleaf = 1 is the lev element info
            leaffin=leafmax
            do lev = levmax,0,-1
                leafstart = leaffin - 4**lev + 1 -1
                nleaflev  = 4**lev !number of leaves per  level
                nleafroot = sqrt(float(nleaflev))
                ires = 4 ** (levmax -lev)
                !resolution of level(i.e. finer level has a ires=1 next level ires=4
                !ires is always a perfect square.i need iresroot because i have to distinguish i from j
                iresroot = sqrt(float(ires))
                im=0;
                nleaf = leafstart
                do mm = 1,nleafroot!jdirection
                  im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                  il=0;
                  do nn = 1, nleafroot ! i direction
                     il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                     xc=0; yc=0; zc=0;sc=0; sourcec=0
                     nleaf = nleaf + 1
                     do m = 1,iresroot
                        do l = 1,iresroot
                            !we use l-1 and m-1 because the 1 is already asigned at il,im
                            !sources at points
                            !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                            !sources at cells
                            leafcount = leafacc + (k - NZs -1 +im + (m-1)) * NNY + (j - NYs + il +(l-1))
                            x  = xs_tmp(leafcount,1)
                            y  = xs_tmp(leafcount,2)
                            z  = xs_tmp(leafcount,3)
                            s  = ds_tmp(leafcount)
                            source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                            
                            
                            xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                            sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                         enddo
                        enddo
                      xs_lev(icount,0,nleaf) = xc/float(ires)
                      ys_lev(icount,0,nleaf) = yc/float(ires)
                      zs_lev(icount,0,nleaf) = zc/float(ires)
                      ds_lev(icount,0,nleaf) = sc
                      source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                      leaffin=leafstart
                    enddo
                 enddo
             enddo
        enddo
    enddo
!do i=1,icount
!      lev=0
!      leafstart = 0

!      do lev = 0,levmax
!          open(14,file='lev'//char(48+lev),access='APPEND')
!          leafstart = leafstart +4**max(0,(lev-1))
!          leaffin   = leafstart +4**(lev) -1
!          if (lev.eq.0) leaffin=1
!          do j=leafstart,leaffin
!             write(14,'(3(e28.17,1x))') xs_lev(i,0,j),ys_lev(i,0,j),zs_lev(i,0,j)
!          enddo
!         close(14)
!      enddo
!   enddo
!   do i=1,nbound
!      write(1511,*) xs_tmp(i,1),xs_tmp(i,2),xs_tmp(i,3)
!   enddo

    leafacc = leafacc + NNZ * NNY
    !---XMAX BOUNDARY----
    i = NXf
    do k = NZs, NZf - istep, istep
        do j = NYs,NYf - istep, istep
            icount = icount + 1
            !nleaf = 1 is the lev element info
            leaffin=leafmax
            do lev = levmax,0,-1
                leafstart = leaffin - 4**lev + 1 -1
                nleaflev  = 4**lev !number of leaves per  level
                nleafroot = sqrt(float(nleaflev))
                ires = 4 ** (levmax -lev)
                !resolution of level(i.e. finer level has a ires=1 next level ires=4
                !ires is always a perfect square.i need iresroot because i have to distinguish i from j
                iresroot = sqrt(float(ires))
                im=0;
                nleaf = leafstart
                do mm = 1,nleafroot!jdirection
                  im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                  il=0;
                  do nn = 1, nleafroot ! i direction
                     il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                     xc=0; yc=0; zc=0;sc=0; sourcec=0
                     nleaf = nleaf + 1
                     do m = 1,iresroot
                        do l = 1,iresroot
                            !we use l-1 and m-1 because the 1 is already asigned at il,im
                            !sources at points
                           !leafcount = leafacc + (k - NZs  +im + (m-1)) * NNY + (j - NYs + 1 + il +(l-1))
                            !sources at cells
                            leafcount = leafacc + (k - NZs -1 +im + (m-1)) * NNY + (j - NYs + il +(l-1))
                            x  = xs_tmp(leafcount,1)
                            y  = xs_tmp(leafcount,2)
                            z  = xs_tmp(leafcount,3)
                            s  = ds_tmp(leafcount)
                            source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                            
                            
                            xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                            sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                         enddo
                        enddo
                      xs_lev(icount,0,nleaf) = xc/float(ires)
                      ys_lev(icount,0,nleaf) = yc/float(ires)
                      zs_lev(icount,0,nleaf) = zc/float(ires)
                      ds_lev(icount,0,nleaf) = sc
                      source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                      leaffin=leafstart
                    enddo
                 enddo
             enddo

        enddo
    enddo

    leafacc = leafacc + NNZ * NNY
    !---YMIN BOUNDARY----(fn points at the y direction)
    j = NYs
    do k = NZs,NZf - istep , istep
        do i = NXs,NXf - istep, istep
            icount = icount + 1
            !nleaf = 1 is the lev element info
            leaffin=leafmax
            do lev = levmax,0,-1
                leafstart = leaffin - 4**lev + 1 -1
                nleaflev  = 4**lev !number of leaves per  level
                nleafroot = sqrt(float(nleaflev))
                ires = 4 ** (levmax -lev)
                !resolution of level(i.e. finer level has a ires=1 next level ires=4
                !ires is always a perfect square.i need iresroot because i have to distinguish i from j
                iresroot = sqrt(float(ires))
                im=0;
                nleaf = leafstart
                do mm = 1,nleafroot!jdirection
                  im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                  il=0;
                  do nn = 1, nleafroot ! i direction
                     il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                     xc=0; yc=0; zc=0;sc=0; sourcec=0
                     nleaf = nleaf + 1
                     do m = 1,iresroot
                        do l = 1,iresroot
                            !we use l-1 and m-1 because the 1 is already asigned at il,im
                            !sources at points
                            !leafcount = leafacc + (k - NZs  +im + (m-1)) * NNY + (i - NXs +1+ il +(l-1))
                            !sources at cells
                            leafcount = leafacc + (k - NZs -1 +im + (m-1)) * NNX + (i - NXs + il +(l-1))
                            x  = xs_tmp(leafcount,1)
                            y  = xs_tmp(leafcount,2)
                            z  = xs_tmp(leafcount,3)
                            s  = ds_tmp(leafcount)
                            source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                            
                            
                            xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                            sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                         enddo
                        enddo
                      xs_lev(icount,0,nleaf) = xc/float(ires)
                      ys_lev(icount,0,nleaf) = yc/float(ires)
                      zs_lev(icount,0,nleaf) = zc/float(ires)
                      ds_lev(icount,0,nleaf) = sc
                      source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                      leaffin=leafstart
                    enddo
                 enddo
             enddo
        enddo
    enddo


    leafacc = leafacc + NNZ * NNX
    !---YMAX BOUNDARY----
    j = NYf
    do k = NZs,NZf - istep, istep
        do i = NXs,NXf - istep, istep
           icount = icount + 1
           !nleaf = 1 is the lev element info
           leaffin=leafmax
           do lev = levmax,0,-1
               leafstart = leaffin - 4**lev + 1 -1
               nleaflev  = 4**lev !number of leaves per  level
               nleafroot = sqrt(float(nleaflev))
               ires = 4 ** (levmax -lev)
               !resolution of level(i.e. finer level has a ires=1 next level ires=4
               !ires is always a perfect square.i need iresroot because i have to distinguish i from j
               iresroot = sqrt(float(ires))
               im=0;
               nleaf = leafstart
               do mm = 1,nleafroot!jdirection
                 im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                 il=0;
                 do nn = 1, nleafroot ! i direction
                    il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                    xc=0; yc=0; zc=0;sc=0; sourcec=0
                    nleaf = nleaf + 1
                    do m = 1,iresroot
                       do l = 1,iresroot
                           !we use l-1 and m-1 because the 1 is already asigned at il,im
                           !sources at points
                           !leafcount = leafacc + (k - NZs  +im + (m-1)) * NNY + (i - NXs +1+ il +(l-1))
                           !sources at cells
                           leafcount = leafacc + (k - NZs -1 +im + (m-1)) * NNX + (i - NXs + il +(l-1))
                           x  = xs_tmp(leafcount,1)
                           y  = xs_tmp(leafcount,2)
                           z  = xs_tmp(leafcount,3)
                           s  = ds_tmp(leafcount)
                           source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                           
                           
                           xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                           sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                        enddo
                       enddo
                     xs_lev(icount,0,nleaf) = xc/float(ires)
                     ys_lev(icount,0,nleaf) = yc/float(ires)
                     zs_lev(icount,0,nleaf) = zc/float(ires)
                     ds_lev(icount,0,nleaf) = sc
                     source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                     leaffin=leafstart
                   enddo
                enddo
            enddo
        enddo
    enddo

    leafacc = leafacc + NNZ * NNX
    
    !---ZMIN BOUNDARY----
    k = NZs
    do j = NYs,NYf - istep, istep
        do i = NXs,NXf - istep, istep
           icount = icount + 1
           !nleaf = 1 is the lev element info
           leaffin=leafmax
           do lev = levmax,0,-1
               leafstart = leaffin - 4**lev + 1 -1
               nleaflev  = 4**lev !number of leaves per  level
               nleafroot = sqrt(float(nleaflev))
               ires = 4 ** (levmax -lev)
               !resolution of level(i.e. finer level has a ires=1 next level ires=4
               !ires is always a perfect square.i need iresroot because i have to distinguish i from j
               iresroot = sqrt(float(ires))
               im=0;
               nleaf = leafstart
               do mm = 1,nleafroot!jdirection
                 im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                 il=0;
                 do nn = 1, nleafroot ! i direction
                    il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                    xc=0; yc=0; zc=0;sc=0; sourcec=0
                    nleaf = nleaf + 1
                    do m = 1,iresroot
                       do l = 1,iresroot
                           !we use l-1 and m-1 because the 1 is already asigned at il,im
                           !sources at points
                           !leafcount = leafacc + (j - NYs +im + (m-1)) * NNX + (i - NXs +1+ il +(l-1))
                           !sources at cells
                           leafcount = leafacc + (j - NYs -1 +im + (m-1)) * NNX + (i - NXs + il +(l-1))
                           x  = xs_tmp(leafcount,1)
                           y  = xs_tmp(leafcount,2)
                           z  = xs_tmp(leafcount,3)
                           s  = ds_tmp(leafcount)
                           source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                           
                           
                           xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                           sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                        enddo
                       enddo
                     xs_lev(icount,0,nleaf) = xc/float(ires)
                     ys_lev(icount,0,nleaf) = yc/float(ires)
                     zs_lev(icount,0,nleaf) = zc/float(ires)
                     ds_lev(icount,0,nleaf) = sc
                     source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                     leaffin=leafstart
                   enddo
                enddo
            enddo
        enddo
    enddo
    leafacc = leafacc + NNY * NNX
    !---ZMAX BOUNDARY
    k = NZf
    do j = NYs,NYf - istep, istep
        do i = NXs,NXf - istep, istep
           icount = icount + 1
           !nleaf = 1 is the lev element info
           leaffin=leafmax
           do lev = levmax,0,-1
               leafstart = leaffin - 4**lev + 1 -1
               nleaflev  = 4**lev !number of leaves per  level
               nleafroot = sqrt(float(nleaflev))
               ires = 4 ** (levmax -lev)
               !resolution of level(i.e. finer level has a ires=1 next level ires=4
               !ires is always a perfect square.i need iresroot because i have to distinguish i from j
               iresroot = sqrt(float(ires))
               im=0;
               nleaf = leafstart
               do mm = 1,nleafroot!jdirection
                 im= (mm-1)*2**(levmax - lev) +1 !jdirection starting node
                 il=0;
                 do nn = 1, nleafroot ! i direction
                    il= (nn-1)*2**(levmax - lev)  +1!idireciton starting node
                    xc=0; yc=0; zc=0;sc=0; sourcec=0
                    nleaf = nleaf + 1
                    do m = 1,iresroot
                       do l = 1,iresroot
                           !we use l-1 and m-1 because the 1 is already asigned at il,im
                           !sources at points
                           !leafcount = leafacc + (j - NYs +im + (m-1)) * NNX + (i - NXs +1+ il +(l-1))
                           !sources at cells
                           leafcount = leafacc + (j - NYs -1 +im + (m-1)) * NNX + (i - NXs + il +(l-1))
                           x  = xs_tmp(leafcount,1)
                           y  = xs_tmp(leafcount,2)
                           z  = xs_tmp(leafcount,3)
                           s  = ds_tmp(leafcount)
                           source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                           
                           
                           xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                           sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
                        enddo
                       enddo
                     xs_lev(icount,0,nleaf) = xc/float(ires)
                     ys_lev(icount,0,nleaf) = yc/float(ires)
                     zs_lev(icount,0,nleaf) = zc/float(ires)
                     ds_lev(icount,0,nleaf) = sc
                     source_bound_lev(icount,neqs:neqf,0,nleaf) = sourcec(neqs:neqf)/sc
                     leaffin=leafstart
                   enddo
                enddo
            enddo

        enddo
    enddo
    
    leafacc = leafacc + NNY * NNX
    if (nbound_lev(0).ne.icount) then 
        write(*,*) 'nbound_lev problem',nbound_lev(0),icount
        stop
    endif
     deallocate(xs_tmp,ds_tmp,s_tmp)
End Subroutine build_level_nbound_3d


Subroutine build_level_nbound_new(NXs,NXf,NYs,NYf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf,neqs,neqf
    integer          :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin, ires, leafacc
    integer          :: ncountlev(0:levmax)
    double precision :: x, y, s, source(neqf), xc, yc, sc, sourcec(neqf)
    character *16    :: filout

    leafmax=1
    leafmax=2
    allocate (source_bound_lev(nbound,neqs:neqf,0:levmax,leafmax)); source_bound_lev     = 0.d0
    allocate (xs_lev(nbound,0:levmax,leafmax)); xs_lev     = 0.d0
    allocate (ys_lev(nbound,0:levmax,leafmax)); ys_lev     = 0.d0
    allocate (ds_lev(nbound,0:levmax,leafmax)); ds_lev     = 0.d0  
    allocate (ilev_t(nbound,0:levmax))         ; ilev_t = 0
    allocate (nbound_lev(0:levmax))             ; nbound_lev = 0
     
    nbound_lev(0)   = nbound
    xs_lev(:, levmax, 1) = 0.5d0 * (x_s(1,:) + x_s(2,:))
    ys_lev(:, levmax, 1) = 0.5d0 * (y_s(1,:) + y_s(2,:))
    ds_lev(:, levmax, 1) = d_s(:)
    do neq=neqs,neqf
       source_bound_lev(:,neq,levmax,1) = source_bound(neq,:)
    enddo

    
    icount   = 0
    leafcount= 0
    leafacc  = 0
    istep    = 2**levmax
    nbound_lev=0;ncountlev=0
    i = Nxs

    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do j=NYs,NYf-istep,istep
          !for all leaves in 
          xc=0; yc=0;sc=0; sourcec=0
          do k=1,istep
             leafcount=leafcount+1
             x  = xs_lev(leafcount,levmax,1)
             y  = ys_lev(leafcount,levmax,1)
             s  = ds_lev(leafcount,levmax,1)
             source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
             xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
          enddo
          icount = icount+1
          xs_lev(icount,lev,2) = xc/float(istep)
          ys_lev(icount,lev,2) = yc/float(istep)
          ds_lev(icount,lev,2) = sc
          source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
          ilev_t(icount,lev)=1
       enddo
       if (mod(NYf-Nys,istep).ne.0.or.NYf-NYs.le.istep) then 
          if (lev.eq.0) then
             write(*,*) 'something is very wrong with levels.lev0 should not be in here'
             STOP
          endif
          icount = icount +1 
          ilev_t(icount,lev)=-1
       endif
       ncountlev(lev)=icount
    enddo
    nbound_lev=ncountlev
   !do lev=0,levmax
   !    write(*,*) ds_lev(1,lev,2)
   !enddo

    i = Nxf
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do j=NYs,NYf-istep,istep
          !for all leaves in 
          xc=0; yc=0; sc=0; sourcec=0
          do k=1,istep
             leafcount=leafcount+1
             x  = xs_lev(leafcount,levmax,1)
             y  = ys_lev(leafcount,levmax,1)
             s  = ds_lev(leafcount,levmax,1)
             source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
             xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
          enddo
          icount = icount+1
          xs_lev(icount,lev,2) = xc/float(istep)
          ys_lev(icount,lev,2) = yc/float(istep)
          ds_lev(icount,lev,2) = sc
          source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
          ilev_t(icount,lev)=1
       enddo
       if (mod(NYf-Nys,istep).ne.0.or.NYf-NYs.le.istep) then 
          if (lev.eq.0) then
             write(*,*) 'something is very wrong with levels.lev0 should not be in here'
             STOP
          endif
          icount = icount +1 
          ilev_t(icount,lev)=-1
       endif
       ncountlev(lev)=icount
    enddo
    nbound_lev=ncountlev

    j = NYs
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do i=NXs,NXf-istep,istep
          !for all leaves in 
          xc=0; yc=0; sc=0; sourcec=0
          do k=1,istep
             leafcount=leafcount+1
             x  = xs_lev(leafcount,levmax,1)
             y  = ys_lev(leafcount,levmax,1)
             s  = ds_lev(leafcount,levmax,1)
             source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
             xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
          enddo
          icount = icount+1
          xs_lev(icount,lev,2) = xc/float(istep)
          ys_lev(icount,lev,2) = yc/float(istep)
          ds_lev(icount,lev,2) = sc
          source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
          ilev_t(icount,lev)=1
       enddo
       if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
          if (lev.eq.0) then
             write(*,*) 'something is very wrong with levels.lev0 should not be in here'
             STOP
          endif
          icount = icount +1 
          ilev_t(icount,lev)=-1
       endif
       ncountlev(lev)=icount
    enddo
    nbound_lev=ncountlev

   

    j = NYf
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do i=NXs,NXf-istep,istep
          !for all leaves in 
          xc=0; yc=0; sc=0; sourcec=0
          do k=1,istep
             leafcount=leafcount+1
             x  = xs_lev(leafcount,levmax,1)
             y  = ys_lev(leafcount,levmax,1)
             s  = ds_lev(leafcount,levmax,1)
             source(neqs:neqf) = source_bound_lev(leafcount,neqs:neqf,levmax,1)
             xc = x + xc;yc= y+ yc; sc= s+ sc;sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
          enddo
          icount = icount+1
          xs_lev(icount,lev,2) = xc/float(istep)
          ys_lev(icount,lev,2) = yc/float(istep)
          ds_lev(icount,lev,2) = sc
          source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
          ilev_t(icount,lev)=1
       enddo
       if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
          if (lev.eq.0) then
             write(*,*) 'something is very wrong with levels.lev0 should not be in here'
             STOP
          endif
          icount = icount +1 
          ilev_t(icount,lev)=-1
       endif
       ncountlev(lev)=icount
    enddo
    nbound_lev=ncountlev
   
    if (nbound_lev(0).ne.nbound) then 
       write(*,*) 'something is very wrong with levels.lev0 bounds should be nbounds'
       STOP

    endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!   do lev=0,levmax
!      write(filout,'(a,i2.2)') 'lev',lev
!      open(27,file=filout)
!      do k=1,nbound_lev(lev)
!         if(ilev_t(k,lev).eq.1) write(27,*)xs_lev(k,lev,2),ys_lev(k,lev,2)
!      enddo
!   enddo
!stop
End Subroutine build_level_nbound_new

!-------------------------------------------------------------------------!
!-> Subroutine calc_normalderiv                                           !
!   This Subroutine calculated normal derivate at the boundary of the     !
!   domain.This together with the appropriate Green function calculate    !
!   the boundary conditions at the domain boundaries.                     !
!   For the calculation of the derivative  a fourth order one sided       !
!   difference approximation is used:                                     !
!   1/h((25/12)f(i,j)-4(i-1,j) + 3f(i-2,j) - 4/3 f(i-3),j) + 0.25f(i-4,j) !
!-------------------------------------------------------------------------!
Subroutine build_level_nbound_3d_new(NXs, NXf, NYs, NYf,NZs,NZf,neqs,neqf)
    Implicit None
    integer,intent(in):: NXs, NXf, NYs, NYf, NZs, NZf,neqs,neqf
    integer          :: icount, istep, lev, nleaf, leafcount, leafmax, leafstart, leaffin, ires, leafacc
    integer          :: nleaflev, nleafroot,ncountlev(0:levmax)
    integer          :: m, l, il, im, mm, nn, iresroot, NNX, NNY, NNZ, NNX0, NNY0, NNZ0
    double precision :: x, y, z, s, source(neqf), xc, yc, zc, sc, sourcec(neqf)
    double precision,allocatable :: xs_tmp(:,:),ds_tmp(:),s_tmp(:,:)
    character *16    :: filout


    leafmax=2
    istep    = 2**levmax
    
    allocate (nbound_lev(0:levmax))             ; nbound_lev = 0
    nbound_lev(0) = nbound
    allocate (source_bound_lev(nbound_lev(0),neqs:neqf,0:levmax,leafmax)); source_bound_lev     = 0.d0
    allocate (xs_lev(nbound_lev(0),0:levmax,leafmax)); xs_lev     = 0.d0
    allocate (ys_lev(nbound_lev(0),0:levmax,leafmax)); ys_lev     = 0.d0
    allocate (zs_lev(nbound_lev(0),0:levmax,leafmax)); ys_lev     = 0.d0
    allocate (ds_lev(nbound_lev(0),0:levmax,leafmax)); ds_lev     = 0.d0  
    allocate (ilev_t(nbound_lev(0),0:levmax)); ilev_t     = 0
     
    allocate (xs_tmp(nbound,3),ds_tmp(nbound),s_tmp(nbound,1:neqf))
     xs_tmp(1:nbound,1) = 0.25d0 * (x_s(1,:) + x_s(2,:) + x_s(3,:) + x_s(4,:))
     xs_tmp(1:nbound,2) = 0.25d0 * (y_s(1,:) + y_s(2,:) + y_s(3,:) + y_s(4,:))
     xs_tmp(1:nbound,3) = 0.25d0 * (z_s(1,:) + z_s(2,:) + z_s(3,:) + z_s(4,:))
    ds_tmp(1:nbound) = d_s(:)
    do neq=neqs,neqf
       s_tmp(1:nbound,neq)  = source_bound(neq,1:nbound)
    enddo

   
!   !It might seem a  little bit complex but it isnt.
!   !tree in 2d surfaces means that the tree has 4 leaf's at each level
!   we start the tree building using the 0-th level  nodes(0th=coarser)
!   since we know the information at levmax(finer)  we loop the finer grid
!   using istep. istep is  2 because,2  in 2 in j make the 4 leafs
!   the stucture is
!   0-th node -->nleaf
!   n leaf is the accumulative number of leaf of the 0th node
!      1
!     /\
!    2  3
!   /\  /\
!  4 5  6 7
! 
!
! so we know that nleaf 4 of 0th node 4 refers at the second level(it's easy to find out that)
!
! the structure is build from levmax to  0th level
! we find the first leaf  for each  level and and we add up one for each leaf on the level until we finish.
! because the structure on the finer level is known in a stuctured 2d grid fashion
! in order to assign the correct nodes to the leads we need to loop through the grid.
! Also we know that each level has members both in i and j direction 
! that's we we want to scan using 2d grid patches
! 
! ires is the resolution of each level.so finer 
!      level has resolution of 1(each member has is the finest you can have)
! 
! for level  levmax -1 ires = 4 in 3d or ires 2 in 2d
!  
! in 3d we split the ires in 2 direction again because of the structured grid fine data we have.
! we need htat splitting beacuse the data for each level is reconstructed from the original fine data
! So for level levmax -1 we need to cycle through the 4 elements(in 3d) ,2 in idirection 2 in j direction
! for level levmax -2 we need to cycle through 16 elements.
!
! Why didn't i use the 4 elements of the previous level?I wanted to use the original data.even thought because
! of the linear interpolation of the data it could be the same using the 4 elements of lev -1
!
! source(lev-1) = Sum(source(levmax)*Area(levmax))/Sum(Area(levmax))
!
! Because it's surface is a 2d grid with it's own numbering. we use leafacc to know the accumulated number 
! when changing sides (we map x_s(nod) to a specific leaf)
!
!
!
! Don't be afraid of  this :leafcount = leafacc + (k- NZs -1 + im  + (m-1) ) * NNY + (j -Nys + il +(l-1))
! node = (j-1)*NY + k
! k,j give the global index : think k - NZs -1 as j-1
!                             think j - NYs    as i
! The global index cycles at 0th level steps(coarser) so
!                            im,il give the starting node for each level
!                            (m-1),(l-1) give the node of the elemmenet for each level(-1 beacuse im,il have one)


!tree defined at cells
    NNZ = NZf  - NZs 
    NNY = NYf  - NYs 
    NNX = NXf  - NXs 

    icount   = 0
    leafcount= 0
    leafacc  = 0
    ncountlev=0
!----------
    nbound_lev=0
    i  = Nxs
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do k = NZs, NZf - istep,istep
          do j = NYs,NYf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NYf-NYs,istep).ne.0.or.NYf-NYs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here 1'
                  write(*,*) istep,Nyf-NYs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NZf-NZs,istep).ne.0.or.NZf-NZs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here 2'
                  write(*,*) istep,NZf-NZs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev

!   do lev=0,levmax
!      write(filout,'(a,i2.2)') 'lev',lev
!      open(27,file=filout)
!      do k=1,nbound_lev(lev)
!         if(ilev_t(k,lev).eq.1) write(27,*)xs_lev(k,lev,2),ys_lev(k,lev,2),zs_lev(k,lev,2)
!      enddo
!   enddo
!stop
!do i=1,icount
!      lev=0
!      leafstart = 0

!      do lev = 0,levmax
!          open(14,file='lev'//char(48+lev),access='APPEND')
!          leafstart = leafstart +4**max(0,(lev-1))
!          leaffin   = leafstart +4**(lev) -1
!          if (lev.eq.0) leaffin=1
!          do j=leafstart,leaffin
!             write(14,'(3(e28.17,1x))') xs_lev(i,0,j),ys_lev(i,0,j),zs_lev(i,0,j)
!          enddo
!         close(14)
!      enddo
!   enddo
!   do i=1,nbound
!      write(1511,*) xs_tmp(i,1),xs_tmp(i,2),xs_tmp(i,3)
!   enddo

    !---XMAX BOUNDARY----
    i  = Nxs
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do k = NZs, NZf - istep,istep
          do j = NYs,NYf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NYf-NYs,istep).ne.0.or.NYf-NYs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here 3'
                  write(*,*) istep,NYf-NYs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NZf-NZs,istep).ne.0.or.NZf-NZs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here 4'
                  write(*,*) istep,NZf-NZs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev

    !---YMIN BOUNDARY----
    j  = NYs
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do k = NZs, NZf - istep,istep
          do i = NXs,NXf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here 5'
                  write(*,*) istep,NXf-NXs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NZf-NZs,istep).ne.0.or.NZf-NZs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here 6'
                  write(*,*) istep,NZf-NZs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev


    !---YMAX BOUNDARY----
    j  = NYf
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do k = NZs, NZf - istep,istep
          do i = NXs,NXf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here 7'
                  write(*,*) istep,NXf-NXs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NZf-NZs,istep).ne.0.or.NZf-NZs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here 8'
                  write(*,*) istep,NZf-NZs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev

    !---ZMIN BOUNDARY----
    k  = NZs
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do j= NYs, NYf - istep,istep
          do i = NXs,NXf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here'
                  write(*,*) istep,NXf-NXs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NYf-NYs,istep).ne.0.or.NYf-NYs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here'
                  write(*,*) istep,NYf-NYs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev

    !---ZMAX BOUNDARY----
    k  = NZf
    do lev=0,levmax
       istep    = 2**lev
       icount=nbound_lev(lev)
       leafcount=nbound_lev(0)
       do j = NYs, NYf - istep,istep
          do i = NXs,NXf - istep,istep
                        !we use l-1 and m-1 because the 1 is already asigned at il,im
                        !sources at points
                        !leafcount = leafacc + (k - NZs +im + (m-1)) * NNY + (j - NYs+1 + il +(l-1))
                        !sources at cells
             xc=0;yc=0;zc=0;sc=0;sourcec=0
             do im =1,istep
                leafcount = leafcount+1
                x  = xs_tmp(leafcount,1)
                y  = xs_tmp(leafcount,2)
                z  = xs_tmp(leafcount,3)
                s  = ds_tmp(leafcount)
                source(neqs:neqf) = s_tmp(leafcount,neqs:neqf)
                
                
                xc = x + xc;yc = y + yc;zc = z + zc; sc= s + sc
                sourcec(neqs:neqf)=source(neqs:neqf)*s + sourcec(neqs:neqf)
             enddo
             icount = icount + 1
             xs_lev(icount,lev,2) = xc/float(istep)
             ys_lev(icount,lev,2) = yc/float(istep)
             zs_lev(icount,lev,2) = zc/float(istep)
             ds_lev(icount,lev,2) = sc
             source_bound_lev(icount,neqs:neqf,lev,2) = sourcec(neqs:neqf)/sc
             ilev_t(icount,lev)=1
           enddo!j
            if (mod(NXf-NXs,istep).ne.0.or.NXf-NXs.le.istep) then 
               if (lev.eq.0) then
                  write(*,*) 'something is very wrong with levels.lev0 should not be in here'
                  write(*,*) istep,NXf-NXs
                  STOP
               endif
               icount = icount +1 
               ilev_t(icount,lev)=-1
            endif
        enddo!k
        if (mod(NYf-NYs,istep).ne.0.or.NYf-NYs.le.istep) then 
           if (lev.eq.0) then
              write(*,*) 'something is very wrong with levels.lev0 should not be in here'
                  write(*,*) istep,NYf-NYs
              STOP
           endif
           icount = icount +1 
           ilev_t(icount,lev)=-1
        endif
        ncountlev(lev)=icount
    enddo!lev
    nbound_lev=ncountlev

    if (nbound_lev(0).ne.nbound) then 
        write(*,*) 'nbound_lev problem',nbound_lev(0),icount
        stop
    endif



     deallocate(xs_tmp,ds_tmp,s_tmp)
End Subroutine build_level_nbound_3d_new
