!--------------------------------------------------------------------------------
!>@file 
!>@brief Poisson solver library.This source code file(pmlib.f90) contains the library 
!!for solving poisson problem on a structured grid with constant DX,DY,DZ
!--------------------------------------------------------------------------------

!>@brief This module defines the variables that will be used internally in the library
!>       All variables are private
module pmlib
    use mkl_dfti
    double precision , save              :: PI, PI2, PI4
    double precision, save               :: XMIN_pm ,XMAX_pm, YMIN_pm, YMAX_pm, ZMIN_pm, ZMAX_pm
    double precision, save               :: DXpm, DYpm, DZpm,DXpm2,DYpm2,DZpm2

    integer ,save                        :: NVR, NXpm, NYpm, NZpm,ND
    integer ,save                        :: NXs_bl(10),NYs_bl(10),NXf_bl(10),NYf_bl(10),NZs_bl(10),NZf_bl(10),NBlocks
 
    integer, save                        :: nbound,levmax
!Here pointers are defined which will be assigned in the external data to save up space
    double precision,pointer             :: SOL_pm(:,:,:,:), RHS_pm(:,:,:,:),QP(:,:),XP(:,:)
    double precision,allocatable         :: SOL_0_pm(:,:,:,:), source_bound(:,:),x_s(:,:),y_s(:,:),z_s(:,:),d_s(:),cos_s(:),sin_s(:)
    double precision,allocatable,save    :: source_bound_lev(:,:,:,:),xs_lev(:,:,:),ys_lev(:,:,:),zs_lev(:,:,:),ds_lev(:,:,:)
    integer,allocatable,save             :: nbound_lev(:),ilev_t(:,:)

    private ::PI,PI2,PI4,XMIN_pm,XMAX_pm,YMIN_pm,YMAX_pm,ZMIN_pm,ZMAX_pm,DXpm,DYpm,DZpm,NVR,NXpm,NYpm,NZPm,ND
    private ::NXs_bl,NYs_bl,NXf_bl,NYf_bl,NZs_bl,NZf_bl,NBlocks,DXpm2,DYpm2,DZpm2
    private ::SOL_pm,RHS_pm,SOL_0_pm,QP,XP
    private ::source_bound,x_s,y_s,z_s,d_s,cos_s,sin_s
    private ::source_bound_lev,xs_lev,ys_lev,zs_lev,ds_lev
    private ::nbound,ilev_t

contains
!--------------------------------------------------------------------------------
!>@function
! SUBROUTINE  pmesh    
!>
!>@author Papis
!>
!>@brief       
!>This is the main subroutine of the Poisson solver library 
!>The input of the subroutine is DSOL_pm,DRHS_pm,DQP,DXP and d velocities.These variables are assigned the 
!>specified pointers.What is also needed is
!> Dpm,NN,NN_bl,Nblocks,Xbound which define the grid
!> Dpm(3)   is DX,DY,DZ
!> NN(3)    is the size of extend domain which coincides with the size of DSOL_pm,DRHS_pm
!> NN_bl(6) is the size of the original grid (smaller domain) in which the poisson problem will be solved
!>          (NN_bl(1:3),the starting nodes of the grid(X,Y,Z) with respect to the extended domain)
!>          (NN_bl(4:6),the last nodes of the grid(X,Y,Z)  with respect to the extended domain)
!>Nblocks is deprecated
!>Xbound(6) Xmin,Ymin,Zmin,Xmax,Ymax,Zmax of the extended domain
!>ibctyp (1 for bio savart law) (2 for infinite domain boundary conditions) 
!>neqs,neqf is an option if we want to solve more than one equation.neqs,neqf shoud coresspond to 
!>DSOL_pm(:,:,:,neqs,neqf)...
!>iynbc is in case we want to add externally bc's for our equation.0 means normal solve (which means the 
!>solution at boundaries is from ibctyp.1 means keeps the solution at the boundaries at what is defined
!> externally)
!>IMPORTANT NOTE: The following library assumes that there are two domains.
!>                -->the original domain and 
!>                -->an extended domain(that's why NN differ from NN_bl)
!>                The solution is found on the extended domain.
!>                IF you don't want an extended domain then set NN_bl(1)=1,NN_bl(4)=NN(1)..NN_bl(2)..
!>                EXTERNALLY
!REVISION HISTORY
!> TODO_dd
!>
!>@param [in] 
!>@param [out] 
!--------------------------------------------------------------------------------
    Subroutine pmesh(DSOL_pm,DRHS_pm,DQP,DXP,Xbound,Dpm,NN,NN_bl,ND,Nblocks,ibctyp,neqs,neqf,iynbc,NVR,itree,levmax)
        Implicit None
        integer,intent(in)                      :: ibctyp,neqs,neqf,iynbc,NVR,itree,levmax
        double precision,intent(in)             :: Xbound(6),Dpm(3)
        integer, intent(in)                     :: NN_bl(6),NN(3),ND,Nblocks
        double precision,intent(inout),target   :: DSOL_pm(:,:,:,:), DRHS_pm(:,:,:,:),DQP(:,:),DXP(:,:)
        double precision                        :: XPM,YPM ,velx, vely
        integer                                 :: i , j, k ,N11,N12,N13,N14,N15 , nb, NXs, NYs, NXf, NYf, NZs,NZf,neq


        double precision :: xi,yi,ksi1,ksi2,th1,th2,w1,w2

        double precision :: R, DX, DY, GreenF, nv
        !--> Define PI
        PI   = 4.d0 * atan(1.d0)
        PI2  = 2.d0 * PI
        PI4  = 4.d0 * PI

        !-->Pass the external data to the namespace of the library
        XMIN_pm = Xbound(1); YMIN_pm = Xbound(2); ZMIN_pm = Xbound(3)
        XMAX_pm = Xbound(4); YMAX_pm = Xbound(5); ZMAX_pm = Xbound(6)
        DXpm    = Dpm(1)   ; DYpm    = Dpm(2)   ; DZpm    = Dpm(3)
        DXpm2   = 2*DXpm   ; DYpm2   = 2*DYpm   ; DZpm2   = 2*DZpm
        !NXs_bl,NXf_bl refer to the starting ending node of the domain we want to solve
        NXpm      = NN(1)   ; NYpm      = NN(2)   ; NZpm      = NN(3)
        NXs_bl(1) = NN_bl(1); NYs_bl(1) = NN_bl(2); NZs_bl(1) = NN_bl(3)
        NXf_bl(1) = NN_bl(4); NYf_bl(1) = NN_bl(5); NZf_bl(1) = NN_bl(6)
        !Assign the pointers to the external data
        SOL_pm     => DSOL_pm    ; RHS_pm     => DRHS_pm    ; QP => DQP; XP => DXP
        !--> 

        !iynbc 1 normal poisson solve.(Calculation of bc's is done here)
        if(iynbc.eq.1) then
            !NXs,NXf is the points we want to calculate boundary conditions
            if (ND.eq.2) then
                if (ibctyp.eq.1) then
                    !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
                    NXs = 1              ! NXs_bl(1)
                    NXf = NXpm           !  NXf_bl(1)
                    NYs = 1              !NYs_bl(1)
                    NYf = NYpm           !NYf_bl(1)
                    call Bounds2D(ibctyp,NXs,NXf,NYs,NYf,neqs,neqf)
                else
                    !Infinite domain boundary conditions(asume zero bc' at the boundary
                    call infdomain(neqs,neqf)
                endif
            else if (ND.eq.3) then
                if (ibctyp.eq.1) then
                    !Biot Savart boundary conditions(set on the ext. domain because that's what will be solved)
                    NXs = 1           !NXs_bl(1)
                    NXf = NXpm        !NXf_bl(1)
                    NYs = 1           !NYs_bl(1)
                    NYf = NYpm        !NYf_bl(1)
                    NZs = 1           !NZs_bl(1)   
                    NZf = NZpm        !NZf_bl(1)
                    call Bounds3D(ibctyp,NXs,NXf,NYs,NYf,NZs,NZf,neqs,neqf)
                else
                    !Infinite domain boundary conditions(asume zero bc' at the boundary
                    call infdomain_3D(neqs,neqf)
                endif
            endif
            
            !Boundary calculations have finished now we have the boundary conditions at NXs,NXf...
            if (ND.eq.2) then
                nb  = 1
                NXs = 1           !NXs_bl(nb)
                NXf = NXpm        !NXf_bl(nb)
                NYs = 1           !NYs_bl(nb)
                NYf = NYpm        !NYf_bl(nb)
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq=neqs,neqf
                    call solve_eq(NXs,NXf,NYs,NYf,neq)
                enddo
            else
                nb =1
                NXs = 1         !NXs_bl(1)
                NXf = NXpm      !NXf_bl(1)
                NYs = 1         !NYs_bl(1)
                NYf = NYpm      !NYf_bl(1)
                NZs = 1         !NZs_bl(1)
                NZf = NZpm
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq=neqs,neqf
                    call solve_eq_3D(NXs,NXf,NYs,NYf,NZs,NZf,neq)
                enddo
            endif
        else
        !In case iynbc = 0 then the boundary calculations are asummed ok from the external data and thus
        !going fo the final poisson solve.
        !IMPORTANT:: because the original idea was to have a solution at the smaller domain (extension of 
        !the domain used for solving the small).The poisson solver is now solved for the small domain,since
        !boundary conditions are assumed ok.
            if (ND.eq.2) then
                nb=1
                NXs = NXs_bl(nb)
                NXf = NXf_bl(nb)
                NYs = NYs_bl(nb)
                NYf = NYf_bl(nb)
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq=neqs,neqf
                    call solve_eq(NXs,NXf,NYs,NYf,neq)
                enddo
            else
                nb =1
                NXs = NXs_bl(1)
                NXf = NXf_bl(1)
                NYs = NYs_bl(1)
                NYf = NYf_bl(1)
                NZs = NZs_bl(1)
                NZf = NZf_bl(1)
                !Solve poisson problem with the boundary conditions set at Nxs,NXf.....
                do neq=neqs,neqf
                    call solve_eq_3D(NXs,NXf,NYs,NYf,NZs,NZf,neq)
                enddo
            endif
        endif

        nullify(SOL_pm,RHS_pm,QP,XP)
    contains
        include 'pmbound.f90'
        include 'pinfdomain.f90'
        include 'pmsolve.f90'

    End Subroutine pmesh

!--------------------------------------------------------------------------------
!>@function
! SUBROUTINE  definepm   
!>
!>@author Papis
!>
!>@brief       
!>definepm  defines the characteristics of the poisson solver grid
!> the Input  is Dpm(3) (DX,DY,DZ)
!>               ndum   (the number of Dpm's which the domain will be extended see above)
!>               Xbound(6) Xmin,Ymin,Zmin,Xmax,Ymax,Zmax of the original domain 
!> the Output is Xbound(6) but of the extended domain
!>               NN(3) (NX,NY,NZ) the nodes of the extended domain (1-->NXpm)
!>               NN_bl(6) the starting and ending nodes of the original domain(defined by xbound at input)
!>               (NXs,NYs,NZs,NXf,NYf,NZf) 
!>
!REVISION HISTORY
!> TODO_dd
!>
!>@param [in] 
!>@param [out] 
!--------------------------------------------------------------------------------
  Subroutine definepm(itype,Xbound,Dpm,ND,ndum,nsize,NN,NN_bl)
        Implicit None
        integer,intent(in)             :: itype,ND,nsize(3)
        integer,intent(in)             :: ndum
        double precision,intent(inout) :: Xbound(6)
        double precision,intent(inout) :: Dpm(3)
        integer,intent(out)            :: NN(3),NN_bl(6)
        integer :: Nblocks,ndum_new(3),nn1,nn2
        double precision :: Xbound_old(6)
        !-> Define Pmesh X,Y,Z min/max boundaries
        if (ND.eq.2) then
            Xbound(6)=0
            Xbound(3)=0
        endif

        ndum_new=ndum
!Itype 1 extends the domain by ndum_new cells
        if (itype.eq.1) then 
            !extend domain by ndum cells
            Xbound(1) = Xbound(1) - ((ndum_new(1)) * Dpm(1))
            Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
            
            Xbound(2) = Xbound(2) - ((ndum_new(2)) * Dpm(2))
            Xbound(5) = Xbound(5) + ((ndum_new(2)) * Dpm(2))
            if(ND.eq.3) then 
               Xbound(3) = Xbound(3) - ((ndum_new(3)) * Dpm(3))
               Xbound(6) = Xbound(6) + ((ndum_new(3)) * Dpm(3))
            endif
            
            !Find number of nodes with the Dpm given from input
            NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1
            NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/ (Dpm(2)))) + 1
            if(ND.eq.3) NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/ (Dpm(3)))) + 1
        else if (itype.eq.4) then 
            Xbound(1) = Xbound(1) - ((ndum_new(1)) * Dpm(1))
            Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
            
            NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1 
            ndum_new(1)  = nsize(1) - mod(NN(1)-1,nsize(1))
           !if (mod(ndum_new(1),2).ne.0) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(1),NN(1)
         ! !   stop
           !endif
           !if (mod(ndum_new(2),2).ne.0) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(2),NN(2)
         ! !   stop
           !endif
           !if (mod(ndum_new(3),2).ne.0.and.ND.eq.3) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(3),NN(3)
         ! !   stop
           !endif
               ndum_new(1)= ndum_new(1)
               Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
               ndum_new(1)= ndum_new(1)+ndum
               NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1
               NN_bl(1) = ndum + 1
               NN_bl(4) = NN(1) - ndum
               return

            

        else if (itype.eq.2) then 
!Itype 2 extends the domain by ndum_new cells and changes Dpm so that the number
!of cells are divided exactly by nsize
            Xbound(1) = Xbound(1) - ((ndum_new(1)) * Dpm(1))
            Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
            
            Xbound(2) = Xbound(2) - ((ndum_new(2)) * Dpm(2))
            Xbound(5) = Xbound(5) + ((ndum_new(2)) * Dpm(2))

            if(ND.eq.3) then 
               Xbound(3) = Xbound(3) - ((ndum_new(3)) * Dpm(3))
               Xbound(6) = Xbound(6) + ((ndum_new(3)) * Dpm(3))
            endif
            
            !Find number of nodes with the Dpm given from input
            NN(1)  = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1
            NN(2)  = int(nint(abs(Xbound(5) - Xbound(2))/ (Dpm(2)))) + 1
            if(ND.eq.3)NN(3)  = int(nint(abs(Xbound(6) - Xbound(3))/ (Dpm(3)))) + 1
            NN(1)  = NN(1)+nsize(1)-mod(NN(1)-1,nsize(1)) 
            NN(2)  = NN(2)+nsize(2)-mod(NN(2)-1,nsize(2)) 
            if(ND.eq.3)NN(3)  = NN(3)+nsize(3)-mod(NN(3)-1,nsize(3)) 
            Dpm(1) = (abs(Xbound(4) - Xbound(1))/ (NN(1)-1))
            Dpm(2) = (abs(Xbound(5) - Xbound(2))/ (NN(2)-1))
            if(ND.eq.3)Dpm(3) = (abs(Xbound(6) - Xbound(3))/ (NN(3)-1))
          ! write(*,*) 'New Dpm(1),Dpm(2),Dpm(3)'
          ! write(*,*) Dpm(1),Dpm(2),Dpm(3)
          ! write(*,*) NN


        else if (itype.eq.3) then 
!Itype 3 extends the domain by ndum_new cells and adds dummy cells at both directions
!so that the total cells are divided by nsize
            Xbound(1) = Xbound(1) - ((ndum_new(1)) * Dpm(1))
            Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
            
            Xbound(2) = Xbound(2) - ((ndum_new(2)) * Dpm(2))
            Xbound(5) = Xbound(5) + ((ndum_new(2)) * Dpm(2))
            
            if(ND.eq.3) then
               Xbound(3) = Xbound(3) - ((ndum_new(3)) * Dpm(3))
               Xbound(6) = Xbound(6) + ((ndum_new(3)) * Dpm(3))
            endif
            NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1 
            NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/ (Dpm(2)))) + 1 
            if(ND.eq.3) NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/ (Dpm(3)))) + 1 
            ndum_new(1)  = nsize(1) - mod(NN(1)-1,nsize(1))
            ndum_new(2)  = nsize(2) - mod(NN(2)-1,nsize(2))
            ndum_new(3)  = nsize(3) - mod(NN(3)-1,nsize(3))
           !if (mod(ndum_new(1),2).ne.0) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(1),NN(1)
         ! !   stop
           !endif
           !if (mod(ndum_new(2),2).ne.0) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(2),NN(2)
         ! !   stop
           !endif
           !if (mod(ndum_new(3),2).ne.0.and.ND.eq.3) then 
           !   write(*,*) 'error sizes',ndum_new,nsize(3),NN(3)
         ! !   stop
           !endif
            if (mod(ndum_new(1),2).eq.0) then
               ndum_new(1)= ndum_new(1)/2
               Xbound(1) = Xbound(1) - ((ndum_new(1)) * Dpm(1))
               Xbound(4) = Xbound(4) + ((ndum_new(1)) * Dpm(1))
               ndum_new(1)= ndum_new(1)+ndum
               NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1
               NN_bl(1) = ndum + 1
               NN_bl(4) = NN(1) - ndum
            else
             !write(*,*) 'aaaaa1'
               nn1 = mod(ndum_new(1),2)
               nn2 = int(ndum_new(1))/int(2)
               Xbound(1) = Xbound(1) - ((nn1+nn2) * Dpm(1))
               Xbound(4) = Xbound(4) + ((nn2) * Dpm(1))
               NN(1) = int(nint(abs(Xbound(4) - Xbound(1))/ (Dpm(1)))) + 1
              !NN_bl(1) = ndum_new(1)+ ndum + 1
              !NN_bl(1) = ndum_new(1)+ ndum + 1
               NN_bl(1) = ndum + 1
               NN_bl(4) = NN(1) - ndum
            endif
            if (mod(ndum_new(2),2).eq.0) then 
               ndum_new(2)= ndum_new(2)/2
               Xbound(2) = Xbound(2) - ((ndum_new(2)) * Dpm(2))
               Xbound(5) = Xbound(5) + ((ndum_new(2)) * Dpm(2))
               ndum_new(2)=ndum_new(2) +ndum !add the initial dummy cells
               NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/ (Dpm(2)))) + 1
              !NN_bl(2) = ndum_new(2) + 1
              !NN_bl(5) = NN(2) - ndum_new(2)
               NN_bl(2) = ndum+1
               NN_bl(5) = NN(2) - ndum
            else
            ! write(*,*) 'aaaaa2'
               nn1 = mod(ndum_new(2),2)
               nn2 = int(ndum_new(2))/int(2)
               Xbound(2) = Xbound(2) - ((nn1+nn2) * Dpm(2))
               Xbound(5) = Xbound(5) + ((nn2) * Dpm(2))
               NN(2) = int(nint(abs(Xbound(5) - Xbound(2))/ (Dpm(2)))) + 1
              !NN_bl(2) = ndum_new(2)+ ndum + 1
              !NN_bl(5) = NN(2) - ndum

               NN_bl(2) = ndum + 1
               NN_bl(5) = NN(2) - ndum
            endif
            if(ND.eq.3) then
            if (mod(ndum_new(3),2).eq.0) then 
               ndum_new(3)= ndum_new(3)/2
               Xbound(3) = Xbound(3) - ((ndum_new(3)) * Dpm(3))
               Xbound(6) = Xbound(6) + ((ndum_new(3)) * Dpm(3))
               ndum_new(3)=ndum_new(3) +ndum !add the initial dummy cells
               NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/ (Dpm(3)))) + 1
              !NN_bl(3) = ndum_new(3) + 1
              !NN_bl(6) = NN(3) - ndum_new(3)
               NN_bl(3) =  ndum+1
               NN_bl(6) = NN(3) - ndum
            else 
             !write(*,*) 'aaaaa3'
               nn1 = mod(ndum_new(3),2)
               nn2 = int(ndum_new(3))/int(2)
               Xbound(3) = Xbound(3) - ((nn1+nn2) * Dpm(3))
               Xbound(6) = Xbound(6) + ((nn2) * Dpm(3))
               NN(3) = int(nint(abs(Xbound(6) - Xbound(3))/ (Dpm(3)))) + 1
              !NN_bl(3) = ndum_new(3) + ndum + 1
              !NN_bl(6) = NN(3) - ndum
               NN_bl(3) = ndum + 1
               NN_bl(6) = NN(3) - ndum
            endif
            endif
            
            if(ND.eq.2) then 
               NN(3)=1
               NN_bl(3) = 1
               NN_bl(6) = 1
            endif
            return

            
        endif
        !2d do not have Z nodes
        if(ND.eq.2)NN(3)=1

        !Define the nodes which the original domain starts (corresponding to Xbound_in)
        NN_bl(1) = ndum_new(1) + 1
        NN_bl(4) = NN(1) - ndum_new(1)

        NN_bl(2) = ndum_new(2) + 1
        NN_bl(5) = NN(2) - ndum_new(2)

        NN_bl(3) = ndum_new(3) + 1
        NN_bl(6) = NN(3) - ndum_new(3)

        !2d do not have Z nodes
        if(ND.eq.2) then
            NN_bl(3) = 1
            NN_bl(6) = 1
        endif
    End Subroutine definepm


end module pmlib
