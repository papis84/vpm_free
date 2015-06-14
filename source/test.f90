Program test_pm
use vpm_lib
use MPI
integer ::NVR
double precision,allocatable:: XPR(:,:),QPR(:,:),UPR(:,:),GPR(:,:)
double precision ::Vref
integer :: my_rank,np,ierr,i

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)


if (my_rank.eq.0) then 
    open(1,file='part10M.dat',form='unformatted')
    read(1) NVR
    read(1) Vref
    allocate(XPR(3,NVR),QPR(4,NVR))
    allocate(UPR(3,NVR),GPR(4,NVR))
    write(*,*) 'NVR=',NVR,Vref
    do i=1,NVR
       read(1) XPR(1,i),XPR(2,i),XPR(3,i),QPR(1,i),QPR(2,i),QPR(3,i)
    enddo

QPR(1:3,:) = -QPR(1:3,:) * Vref
QPR(4,:) =Vref
endif
call vpm(XPR,QPR,UPR,GPR,NVR,3)

call MPI_FINALIZE(ierr)
End Program test_pm
