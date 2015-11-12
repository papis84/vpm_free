!-------Test problems
subroutine fUi_HillsVortex_1(CP,a,us,z0,Uind,Grad,Defm,Vort)
        ! Argumens
        integer, parameter :: MK=8
        real(MK),dimension(3), intent(in) :: CP !< Control Point
        real(MK) :: a  !< hill vortex radius
        real(MK) :: us !< hills self induced velocity us=-5/2 U0  (ie convects along -ez)
        real(MK) :: z0 !< Position on the z axis
        real(MK),dimension(3),intent(out) :: Uind !< Induced velocity
        real(MK),dimension(9),intent(out) :: Grad !< Gradient
        real(MK),dimension(3),intent(out) :: Defm !< Deformation
        real(MK),dimension(3),intent(out) :: Vort !< Vorticity
        ! Variables
        real(MK)              :: r3d,rho
        real(MK)              :: x,y,z
        real(MK)              :: urho
        real(MK)              :: uz
        real(MK),dimension(3) :: e_phi  ! Tangential coordinate
        real(MK),dimension(3) :: e_rho  ! Cylindrical radial
        real(MK)              :: om_phi
        real(MK)              :: defm_phi
        !
        x=CP(1)
        y=CP(2)
        z=CP(3)-z0  ! Note the allowed offset

        r3d=sqrt(x**2 + y**2 + z**2)
        rho=sqrt(x**2 + y**2)
        ! Tangential coordinate (tangent to (y,z))
        if (rho/a<1.e-12) then
            e_phi(1:3)=0.0_MK
            e_rho(1:3)=0.0_MK
        else
            e_phi(1)=-y/rho
            e_phi(2)= x/rho
            e_phi(3)=0
            e_rho(1)= x/rho
            e_rho(2)= y/rho
            e_rho(3)=0
        endif
        if(r3d<a) then
            ! --- Inside the sphere
            ! Velocity
            uz     =3._MK/5._MK*us*(1._MK-(2._MK*rho**2+z**2)/(a**2))+us*2._MK/5._MK ;
            urho   =3._MK/5._MK*us*(rho*z)/(a**2);
            ! Vorticity (along e_phi)
            om_phi= 3._MK*us*rho/a**2  != -15/2  * (u0 rho)/a^2 = 3 us rho/a**2
            ! Deformation
            defm_phi= 9._MK/5._MK*us**2/a**2*rho*z ! =45/4 * uo^2/a^2 * rho z
            ! Gradient
            Grad(1:9)=0.0_MK !TODO
        else
            ! --- Outside of the sphere
            ! Velocity
            uz    =2._MK/5._MK*us* (((a**2)/(z**2+rho**2))**(5._MK/2._MK))*(2._MK*z**2-rho**2)/(2._MK*a**2);
            urho  =3._MK/5._MK*us*rho*z/(a**2)*(((a**2)/(z**2+rho**2))**(5._MK/2._MK));
            ! Vorticity
            om_phi   = 0.0_MK
            ! Deformation
            defm_phi = 0.0_MK
            ! Gradient
            Grad(1:9)=0.0_MK !TODO
        endif
        !
        Uind(1) = urho * e_rho(1)
        Uind(2) = urho * e_rho(2)
        Uind(3) = uz
        !
        Defm(1) =  defm_phi * e_phi(1)
        Defm(2) =  defm_phi * e_phi(2)
        Defm(3) = 0.0_MK
        !
        Vort(1) = om_phi * e_phi(1)
        Vort(2) = om_phi * e_phi(2)
        Vort(3) = 0.0_MK

    end subroutine



Subroutine hill_assign(NN,NN_bl,Xbound,Dpm,RHS_pm_bl,neqpm)
    use test_app
     integer,intent(in) :: NN(3),NN_bl(6),neqpm
     double precision,intent(in) :: Xbound(6),Dpm(3)
     double precision,intent(inout) :: RHS_pm_bl(neqpm,NN(1),NN(2),NN(3))
     double precision :: CP(3),Uind(3),Grad(9),Defm(3),Vort(3)
     integer          :: i,j,k
     logical          :: I_EXIST
      RHS_pm_bl(1:3,:,:,:)=0.0d0
      INQUIRE (FILE='hillref.plt', EXIST=I_EXIST)
      if(I_EXIST.eqv..false.) then 
         open(1,file='hillref.dat')
         WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "U" "V" "W" "VORTX" "VORTY" "VORTZ"'
         WRITE(1,*)'ZONE I=',NN_bl(4)-NN_bl(1)+1,' J=',NN_bl(5)-NN_bl(2)+1,&
            ' K=',NN_bl(6)-NN_bl(3) + 1 ,' F=POINT'
      endif
      allocate(analytic_sol(6,NN(1),NN(2),NN(3)))
        do k=NN_bl(3), NN_bl(6)
            do j=NN_bl(2), NN_bl(5)
                do i=NN_bl(1),NN_bl(4)
                    CP(1)=Xbound(1)+(I-1)*Dpm(1)
                    CP(2)=Xbound(2)+(J-1)*Dpm(2)
                    CP(3)=Xbound(3)+(K-1)*Dpm(3)
                    Uind=0.0d0; Grad=0.0d0; Defm=0.d0; Vort=0.0d0
                    call fUi_HillsVortex_1(CP,1.0d0,-1.0d0,0.0d0,Uind,Grad,Defm,Vort)
                    RHS_pm_bl(1:3,i,j,k)=-Vort(1:3)

                    analytic_sol(1:3,i,j,k)= Uind(1:3)
                    analytic_sol(4:6,i,j,k)= Defm(1:3)

                    write(1,'(9(e28.17,1x))') CP(1),CP(2),CP(3),Uind(1),Uind(2),Uind(3),&
                                               Vort(1),Vort(2),Vort(3)
                enddo
            enddo
        enddo
        close(1)
    !   ---FOR PLOTTING PURPOSES ONLY
    if (I_EXIST.eqv..false.) then 
    call system('~/bin/preplot hillref.dat >/dev/null')
    call system('rm hillref.dat')
    endif

End Subroutine hill_assign

Subroutine hill_error(NN,NN_bl,Xbound,Dpm,SOL_pm,velvrx_pm,velvry_pm,velvrz_pm)
    use test_app
    integer,intent(in) :: NN(3),NN_bl(6)
    double precision,intent(in)  :: Xbound(6),Dpm(3)
    double precision,intent(in)  :: SOL_pm(7,NN(1),NN(2),NN(3))
    double precision,intent(in)  :: velvrx_pm(NN(1),NN(2),NN(3))
    double precision,intent(in)  :: velvry_pm(NN(1),NN(2),NN(3))
    double precision,intent(in)  :: velvrz_pm(NN(1),NN(2),NN(3))
    double precision,allocatable :: error(:,:,:,:)
    double precision :: CP(3),mean_err(7),max_err(7)
    integer          :: i,j,k,l,NNX,NNY,NNZ

    

    allocate(error(7,NN(1),NN(2),NN(3)))
    open(1,file='hillerror.dat')
    WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "Err-u" "Err-v "Err-w" "Err-defmx" &
                      "Err-defmy" "Err-defmz" "Err-Vmag"'
    WRITE(1,*)'ZONE I=',NN_bl(4)-NN_bl(1)-1,' J=',NN_bl(5)-NN_bl(2)-1,&
       ' K=',NN_bl(6)-NN_bl(3) - 1 ,' F=POINT'
    
    max_err=0; mean_err=0
    do k=NN_bl(3) + 1, NN_bl(6) - 1
        do j=NN_bl(2) + 1, NN_bl(5) - 1
            do i=NN_bl(1) + 1,NN_bl(4) -1
               CP(1)=Xbound(1)+(I-1)*Dpm(1)
               CP(2)=Xbound(2)+(J-1)*Dpm(2)
               CP(3)=Xbound(3)+(K-1)*Dpm(3)
               
               error(1,i,j,k)= abs(velvrx_pm(i,j,k) - analytic_sol(1,i,j,k))
               error(2,i,j,k)= abs(velvry_pm(i,j,k) - analytic_sol(2,i,j,k))
               error(3,i,j,k)= abs(velvrz_pm(i,j,k) - analytic_sol(3,i,j,k))

               error(4,i,j,k)= abs(SOL_pm(1,i,j,k) - analytic_sol(4,i,j,k))
               error(5,i,j,k)= abs(SOL_pm(2,i,j,k) - analytic_sol(5,i,j,k))
               error(6,i,j,k)= abs(SOL_pm(3,i,j,k) - analytic_sol(6,i,j,k))
               error(7,i,j,k)= abs(sqrt(velvrx_pm(i,j,k)**2 + velvry_pm(i,j,k)**2 +  &
                                        velvrz_pm(i,j,k)**2) -sqrt(analytic_sol(1,i,j,k)**2 +&
                                        analytic_sol(2,i,j,k)**2 +analytic_sol(3,i,j,k)**2))
               write(1,'(10(e28.17,1x))') CP(1),CP(2),CP(3),error(1,i,j,k),error(2,i,j,k), &
                                          error(3,i,j,k),error(4,i,j,k),error(5,i,j,k),    &
                                          error(6,i,j,k),error(7,i,j,k)

              do l=1,7
                 max_err(l) = max(max_err(l),error(l,i,j,k))
              enddo
            enddo
        enddo
    enddo
    
    write(*,*)'----Maximum Velocity Error-----:',max_err(7)*100,'%'
    call system('~/bin/preplot hillerror.dat >/dev/null')
    call system('rm hillerror.dat')
    return
End Subroutine hill_error

Subroutine definevort(RHS_pm,Xbound,Dpm,NN,NN_bl)
  use test_app
  Implicit None
 
  double precision, intent(in)   :: Xbound(6),Dpm(3)
  integer,intent(in)             :: NN(3),NN_bl(6)
  double precision, intent(inout):: RHS_pm(7,NN(1),NN(2),NN(3))
  double precision               :: X(4), XC(4), Y(4), YC(4), ksi1,ksi2,th1,th2,xi,yi,w1,w2,Vol,ANG,dens1,dens2
  integer                        :: i ,j ,k, NXpm1, NYpm1, NZpm1, NCELLSpm, npar,ndumc
  double precision               :: xc1,xc2,yc1,yc2,rmax,r,th,dth,dr,ths,rs,PI
  logical          :: I_EXIST
    xc1=0;xc2=0;
    yc1=-1.5;yc2=1.5
    PI = 4.*datan(1.d0)
    allocate(analytic_sol(1,NN(1),NN(2),NN(3)));analytic_sol=0.d0
    INQUIRE (FILE='vortref.plt', EXIST=I_EXIST)
    if(I_EXIST.eqv..false.) then 
       open(1,file='vortref.dat')
       WRITE(1,'(a190)')'VARIABLES = "X" "Y" "VORT" "PSI"'
       WRITE(1,*)'ZONE I=',NN_bl(4)-NN_bl(1)+1,' J=',NN_bl(5)-NN_bl(2)+1,&
          ' K=',1,' F=POINT'
    endif
    do j = NN_bl(2),NN_bl(5)
        do i = NN_bl(1),NN_bl(4)
            xi=Xbound(1) + (i-1)*Dpm(1)
            yi=Xbound(2) + (j-1)*Dpm(2)
            ksi1=sqrt((xi-xc1)**2+(yi-yc1)**2)
            ksi2=sqrt((xi-xc2)**2+(yi-yc2)**2)
            th1=atan2((yi-yc1),xi-xc1)
            th2=atan2((yi-yc2),xi-xc2)
            if (th1.lt.0.d0) th1=th1+2.d0*PI
            if (th2.lt.0.d0) th2=th2+2.d0*PI
            w1   = (2.d0 - ksi1**2) * exp (0.5d0*(1.d0-ksi1**2))
            w2   =-(2.d0 - ksi2**2 )* exp (0.5d0*(1.d0-ksi2**2))
            RHS_pm(1,i,j,1)=-(w1+w2)
            analytic_sol(1,i,j,1) = exp(0.5d0*(1.d0-ksi1**2))-exp (0.5d0*(1.d0-ksi2**2))
            write(1,'(4(e28.17,1x))') xi,yi,-RHS_pm(1,i,j,1),analytic_sol(1,i,j,1)
         enddo
    enddo

    if (I_EXIST.eqv..false.) then 
    call system('~/bin/preplot vortref.dat >/dev/null')
    call system('rm vortref.dat')
    endif
 End Subroutine definevort

Subroutine vort_error(NN,NN_bl,Xbound,Dpm,SOL_pm)
    use test_app
    integer,intent(in) :: NN(3),NN_bl(6)
    double precision,intent(in)  :: Xbound(6),Dpm(3)
    double precision,intent(in)  :: SOL_pm(7,NN(1),NN(2),NN(3))
    double precision,allocatable :: error(:,:,:,:)
    double precision :: CP(3),mean_err(7),max_err(7),analytic_max
    integer          :: i,j,k,l,NNX,NNY,NNZ

    

    allocate(error(1,NN(1),NN(2),NN(3)))
    open(1,file='vorterror.dat')
    WRITE(1,'(a190)')'VARIABLES = "X" "Y" "Z" "Err-psi" '
    WRITE(1,*)'ZONE I=',NN_bl(4)-NN_bl(1)+1,' J=',NN_bl(5)-NN_bl(2)+1,&
       ' K=', 1 ,' F=POINT'
    
      analytic_max = maxval(abs(analytic_sol(1,:,:,:)))
      max_err=0; mean_err=0
        do j=NN_bl(2) , NN_bl(5) 
            do i=NN_bl(1) ,NN_bl(4)
               CP(1)=Xbound(1)+(I-1)*Dpm(1)
               CP(2)=Xbound(2)+(J-1)*Dpm(2)
               

               error(1,i,j,1)= abs(SOL_pm(1,i,j,1) - analytic_sol(1,i,j,1))
               write(1,'(10(e28.17,1x))') CP(1),CP(2),CP(3),error(1,i,j,1) 
               max_err(1) = max(max_err(1),error(1,i,j,1))
            enddo
        enddo
    
    write(*,*)'----Maximum Phi Error-----:',max_err(1)/analytic_max*100,'%'
    call system('~/bin/preplot vorterror.dat >/dev/null')
    call system('rm vorterror.dat')
    return
End Subroutine vort_error

