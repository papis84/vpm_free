Program test_pm
use vpm_lib
use vpm_size
use test_mod
use MPI
double precision ::Vref
integer :: my_rank,np,ierr,i

call MPI_INIT(ierr)
call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)


if (my_rank.eq.0) then 
    open(1,file='part10M.dat',form='unformatted')
    read(1) NVR_ext
    read(1) Vref
    allocate(XPR(3,NVR_ext),QPR(4,NVR_ext))
    allocate(UPR(3,NVR_ext),GPR(3,NVR_ext))
    write(*,*) 'NVR=',NVR_ext,Vref
    do i=1,NVR_ext
       read(1) XPR(1,i),XPR(2,i),XPR(3,i),QPR(1,i),QPR(2,i),QPR(3,i)
    enddo

QPR(1:3,:) = -QPR(1:3,:) * Vref
QPR(4,:) =Vref
endif
call vpm_init(XPR,QPR,UPR,GPR,NVR_ext,3)
call define_sizes(1)
call remesh_particles_3d(1)
do i=1,1000
  call vpm(XPR,QPR,UPR,GPR,NVR_ext,3,i)
  if (mod(i,1).eq.0) call remesh_particles_3d(1)
enddo

call MPI_FINALIZE(ierr)
End Program test_pm
