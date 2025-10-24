program gera_num_arestas
!Program that creates a file with the Element-Edge relationship,
!with information about Elements, Nodes and Edges from Tetgen files.
!Autor: Anderson Almeida
implicit none

integer :: i,j,k,l,m,n,nel,aux,cnt,ntedj,nglob
integer :: ni,nj,outl
integer :: nedj=6, nnoele=4
integer,allocatable :: mel(:,:), edj(:,:), cedj(:,:), tedj(:,:)
integer,allocatable :: qtd(:), gedgs(:,:), conts(:)

integer,allocatable :: ind_edg_glob(:,:),bd(:)
integer :: n1, n2, flbd=-1
integer:: N_bd, lixo
real(8), parameter:: small=1.0d-5
real(4) :: auxr
real(8):: xi, xf, yi, yf, zi, zf
real(8), allocatable:: coords(:,:)
character(len=30) :: files

!----
call GET_COMMAND_ARGUMENT(1, files) 

!----


!Reading the coordinates of the mesh ends, to generate the edges on the border.
!
open(99, file=trim(files)//'_vfe.info', status='old', action='read')
read(99,*) xi, xf
read(99,*) yi, yf
read(99,*) zi, zf
close(99)
print *, xi, xf
print *, yi, yf
print *, zi, zf


open(10,file=trim(files)//'.1.ele',status='old',action='read')
read(10,*)nel
allocate( mel(nel,nnoele), edj(nel,nedj), cedj(nel,12))
do i = 1,nel
   read(10,*)aux,mel(i,1:nnoele)
   cedj(i,1)  = mel(i,1); cedj(i,2)  = mel(i,2)
   cedj(i,3)  = mel(i,1); cedj(i,4)  = mel(i,3)
   cedj(i,5)  = mel(i,1); cedj(i,6)  = mel(i,4)
   cedj(i,7)  = mel(i,2); cedj(i,8)  = mel(i,3)
   cedj(i,9)  = mel(i,4); cedj(i,10) = mel(i,2)
   cedj(i,11) = mel(i,3); cedj(i,12) = mel(i,4)
end do
close(10)

open(20,file=trim(files)//'.1.edge',status='old',action='read')
read(20,*)ntedj
allocate(tedj(ntedj,2)) ; 
allocate( bd(ntedj) ) ; 
bd = 0;
do i = 1,ntedj
   read(20,*)aux,tedj(i,:),bd(i)
end do
close(20)

open(30,file=trim(files)//'.1.node',status='old',action='read')
read(30,*)nglob
allocate( coords(Nglob, 3) )
do i = 1, Nglob
   read(30,*) lixo, coords(i,1:3)
end do
close(30)

allocate(qtd(nglob)); 
qtd = 0;
N_bd = 0

do i = 1,ntedj
   n1 = tedj(i,1) ; 
   n2 = tedj(i,2)
   if(n1 < n2)then
      qtd(n1) = qtd(n1) + 1
   else
      qtd(n2) = qtd(n2) + 1
   end if

!  Selects edges of the boundary.
!
   if ( (abs(coords(N1,3)-zi) < small .and. abs(coords(N2,3)-zi) < small ) .or. &
         (abs(coords(N1,3)-zf) < small .and. abs(coords(N2,3)-zf) < small ) )  then
         N_bd = N_bd + 1
         bd(N_bd) = i
   else if ( (abs(coords(N1,2)-yi) < small .and. abs(coords(N2,2)-yi) < small) .or. &
               (abs(coords(N1,2)-yf) < small .and. abs(coords(N2,2)-yf) < small) )  then
         N_bd = N_bd + 1
         bd(N_bd) = i
   else if ( (abs(coords(N1,1)-xi) < small .and. abs(coords(N2,1)-xi) < small) .or. &
               (abs(coords(N1,1)-xf) < small .and. abs(coords(N2,1)-xf) < small) )  then
         N_bd = N_bd + 1
         bd(N_bd) = i
   end if

end do

open(77, file=trim(files)//'_vfe.front', status='replace', action='write' )

print *, 'N_bd', N_bd
write(77,'(I7)') N_bd
do i = 1, N_bd
   write(77,*) bd(i)
end do
close(77)

!print*,'Highest Value:',maxval(qtd)

allocate(gedgs(nglob,maxval(qtd)), conts(nglob))
allocate(ind_edg_glob(nglob,maxval(qtd)))

gedgs=0; conts=0;

do i = 1,ntedj
   n1 = tedj(i,1)
   n2 = tedj(i,2)
   if(n1 < n2)then
      conts(n1) = conts(n1) + 1
      gedgs(n1,conts(n1)) = i
      ind_edg_glob(n1,conts(n1)) = n2
   else
      conts(n2) = conts(n2) + 1
      gedgs(n2,conts(n2)) = i
      ind_edg_glob(n2,conts(n2)) = n1
   end if
end do

do i = 1,nel
   m = 0
   do j = 1,nedj
      ni = cedj(i,m+1); nj = cedj(i,m+2)
      l = 0
      do k = 1,conts(ni) 
           l = l + 1
           if( nj == ind_edg_glob(ni,k) )then        
              edj(i,j)=gedgs(ni,k)
              exit
           end if
      end do
      if( l == conts(ni) )then 
         aux = ni
         ni = nj
         nj = aux
         do k = 1,conts(ni) 
            if( nj == ind_edg_glob(ni,k) )then        
               edj(i,j)=gedgs(ni,k)
               exit
            end if
         end do
      end if
      m = m + 2
   end do
end do

!--Writing the file that stores the Element-Edge relationship

open(30,file=trim(files)//'.1.nedge',status='replace',action='write')

write(30,*)nel
do i = 1,nel
   write(30,'(6i12," ")') edj(i,:)
end do
!print*,maxval(edj)
close(30)



end program

