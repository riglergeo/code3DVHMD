!This code creates the base nodes for tetgen to generate the finite element mesh.
!The nodes for the total region, the anomalous body, and the layers (3D block) will be created, 
!as well as the nodes at the measurement points.
!Author: Cícero Regis
!date of creation:03/10/2021
!
!
module subrotinas_do_gera_poly
implicit none
integer,parameter:: re = kind(1.d0) 
real(re),parameter:: pi = 3.1415926535897932384626433832795_re
save

contains

!subroutine for reading information from the block.in file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine leitura( arq_in, Nl, z_layers, N_obs, xp, yp, z_obs, N_block, x0_bl, y0_bl, z0_bl, lx_bl, ly_bl, lz_bl, Tam_cubo )
implicit none
integer:: j
integer,intent(out):: Nl
integer, intent(out):: N_obs
real(re), intent(out):: Tam_cubo
real(re), allocatable, dimension(:), intent(out):: z_layers, xp, yp
real(re), intent(out):: z_obs
real(re):: dx, dy
integer, intent(out):: N_block
real(re), allocatable, dimension(:), intent(out):: x0_bl, y0_bl, z0_bl, lx_bl, ly_bl, lz_bl
real(re):: Cam_ar
real(re):: h      
real(re):: x_i, y_i, x_f, y_f
real(re):: lixo
character(50),intent(out):: arq_in
character(5):: header

call GET_COMMAND_ARGUMENT(1, arq_in)

open(33, file=arq_in, status='old', action='read')
do
	read(33,*) header
	if (header=='Backg') exit
end do
read(33,*) lixo, Cam_ar
read(33,*) Nl
allocate( z_layers(-1:Nl) )
z_layers(-1) = - Cam_ar
z_layers(0) = 0.d0
do j = 1, Nl
	read(33,*) lixo, lixo, h
	z_layers(j) = z_layers(j-1) + h
end do
Tam_cubo = h


rewind(33)
do
    read(33,*) header
    if ( header == 'Block' ) exit
end do
read(33,*) N_block
allocate( x0_bl(N_block), y0_bl(N_block), z0_bl(N_block) )
allocate( lx_bl(N_block), ly_bl(N_block), lz_bl(N_block) )
do j = 1, N_block
    read(33,*) x0_bl(j), y0_bl(j), z0_bl(j)
    read(33,*) lx_bl(j), ly_bl(j), lz_bl(j)
    read(33,*)
end do

rewind(33)
do
	read(33,*) header
	if (header=='Perfi') exit
end do
read(33,*) N_obs
read(33,*) x_i, y_i
read(33,*) x_f, y_f
read(33,*) z_obs
!read(33,*) offset
close(33)

open(44,file='pontos_obs.dat',status='replace')
write(44,*) N_obs
allocate( xp(N_obs), yp(N_obs) )
dx = (x_f - x_i) / (N_obs-1)
dy = (y_f - y_i) / (N_obs-1)
do j = 1, N_obs
    xp(j) = x_i + (j-1)*dx
    yp(j) = y_i + (j-1)*dy
    write(44,*) xp(j), yp(j)
end do
write(44,*) z_obs
close(44)

end subroutine leitura
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module subrotinas_do_gera_poly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gera_poly_em34
!
!
use subrotinas_do_gera_poly
implicit none
integer:: Nl              
integer:: N_block
integer:: N_obs
real(re), allocatable:: z_layers(:)
real(re):: Tam_cubo, Xi_B, xf_B, yi_B, yf_B, zi_B, zf_B  
real(re):: x0, y0
real(re):: x_obs, y_obs, z_obs, offset, cos_tt, sin_tt
real(re):: ar
real(re), allocatable, dimension(:):: q, xp, yp
real(re), allocatable, dimension(:):: x0_bl, y0_bl, z0_bl, lx_bl, ly_bl, lz_bl
integer:: NNos, Nnos_tot, Nfac, Nf_tot     
integer:: No_bl_1, No1_linha
integer:: i, j, k, No_ex(8)
character(50):: arq_in, arq_poly, arq_info
integer:: Posicao_ponto

call leitura( arq_in, Nl, z_layers, N_obs, xp, yp, z_obs, N_block, x0_bl, y0_bl, z0_bl, lx_bl, ly_bl, lz_bl, Tam_cubo )
allocate( q(0:Nl+N_block) )
write(*,*) Tam_cubo

!Creates the .poly file with the same name as the input file.
posicao_ponto = index(arq_in,'.',.true.) - 1					
if (posicao_ponto == -1) posicao_ponto = len_trim(arq_in)
arq_poly = arq_in(1:posicao_ponto)//'.poly'				
arq_info = arq_in(1:posicao_ponto)//'_vfe.info'

open(22, file=arq_poly, status='replace', action='write')


x0 = 0.5d0 * (xp(N_obs) + xp(1))
y0 = 0.5d0 * (yp(N_obs) + yp(1))

! Border coordinates.
!
xi_B = x0 - Tam_cubo
xf_B = x0 + Tam_cubo
yi_B = y0 - Tam_cubo
yf_B = y0 + Tam_cubo
zi_B = z_layers(-1)
zf_B = z_layers(Nl)

Nnos_tot =  4*(Nl+2) + 8*(N_block) + 4*(N_obs) !+ 8

!Write the border nodes, numbering the corners of the grid by column.
!
write(22,*) '# Nodes'
Nnos = 0
write(22,*) Nnos_tot, 3, 0, 0
do j = -1, Nl
	Nnos = Nnos + 1
	write(22,*) Nnos, xi_b, yi_b, z_layers(j)
end do
do j = -1, Nl
	Nnos = Nnos + 1
	write(22,*) Nnos, xf_b, yi_b, z_layers(j)
end do
do j = -1, Nl
	Nnos = Nnos + 1
	write(22,*) Nnos, xf_b, yf_b, z_layers(j)
end do
do j = -1, Nl
	Nnos = Nnos + 1
	write(22,*) Nnos, xi_b, yf_b, z_layers(j)
end do

!creates the block nodes.
!
No_bl_1 = Nnos + 1
do j = 1, N_block
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)- 0.5d0*lx_bl(j)), (y0_bl(j) - 0.5d0*ly_bl(j)), (z0_bl(j) - 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)+ 0.5d0*lx_bl(j)), (y0_bl(j) - 0.5d0*ly_bl(j)), (z0_bl(j) - 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)+ 0.5d0*lx_bl(j)), (y0_bl(j) + 0.5d0*ly_bl(j)), (z0_bl(j) - 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)- 0.5d0*lx_bl(j)), (y0_bl(j) + 0.5d0*ly_bl(j)), (z0_bl(j) - 0.5d0*lz_bl(j))

    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)- 0.5d0*lx_bl(j)), (y0_bl(j) - 0.5d0*ly_bl(j)), (z0_bl(j) + 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)+ 0.5d0*lx_bl(j)), (y0_bl(j) - 0.5d0*ly_bl(j)), (z0_bl(j) + 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)+ 0.5d0*lx_bl(j)), (y0_bl(j) + 0.5d0*ly_bl(j)), (z0_bl(j) + 0.5d0*lz_bl(j))
    Nnos = Nnos + 1
    write(22,*) Nnos, (x0_bl(j)- 0.5d0*lx_bl(j)), (y0_bl(j) + 0.5d0*ly_bl(j)), (z0_bl(j) + 0.5d0*lz_bl(j))
end do

ar = 5.0d-2
cos_tt = (xp(N_obs)-xp(1)) / sqrt( (xp(N_obs)-xp(1))**2 + (yp(N_obs)-yp(1))**2 )
sin_tt = (yp(N_obs)-yp(1)) / sqrt( (xp(N_obs)-xp(1))**2 + (yp(N_obs)-yp(1))**2 )

!Observation line nodes
!A regular tetrahedron centered at each observation point.
do i = 1, N_obs
    x_obs = xp(i) 
    y_obs = yp(i) 

    Nnos = Nnos + 1
    write(22,*) Nnos, (x_obs - ar/2.d0), (y_obs - sqrt(3.d0)*ar/6.d0), (z_obs + sqrt(6.d0)*ar/12.d0)
    Nnos = Nnos + 1
    write(22,*) Nnos, (x_obs + ar/2.d0), (y_obs - sqrt(3.d0)*ar/6.d0), (z_obs + sqrt(6.d0)*ar/12.d0)
    Nnos = Nnos + 1
    write(22,*) Nnos, (x_obs), (y_obs + sqrt(3.d0)*ar/3.d0), (z_obs + sqrt(6.d0)*ar/12.d0)
    Nnos = Nnos + 1
    write(22,*) Nnos, (x_obs), (y_obs), (z_obs - sqrt(6.d0)*ar/4.d0)

 
end do

!	Facets
write(22,*) '# Facets'
Nf_tot = 4*(Nl+1) + Nl + 2 + 6*(N_block) !+ 6
Nfac = 0
write(22,*) Nf_tot, 0

!creates the vertical faces
do i = 1, 3
	Do j = 1, Nl+1
		Nfac = Nfac + 1
		write(22,*) 1, 0
		k = (i-1)*(Nl+2) + j
		write(22,*) 4, k, k+1, k+Nl+3, k+Nl+2
	end do
end do
i = 4
Do j = 1, Nl+1
	Nfac = Nfac + 1
	write(22,*) 1, 0
	k = (i-1)*(Nl+2) + j
	write(22,*) 4, k, k+1, j+1, j
end do

!creates the horizontal faces
!
do j = 1, Nl+2
	Nfac = Nfac+1
	write(22,*) 1, 0
	write(22,*) 4, j, j+Nl+2, j+2*(Nl+2), j+3*(Nl+2)
end do

!creates the faces of the blocks
!
do j = 1, N_block
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1, No_bl_1+1, No_bl_1+2, No_bl_1+3
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1+4, No_bl_1+5, No_bl_1+6, No_bl_1+7
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1, No_bl_1+4, No_bl_1+5, No_bl_1+1
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1+1, No_bl_1+5, No_bl_1+6, No_bl_1+2
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1+2, No_bl_1+6, No_bl_1+7, No_bl_1+3
    Nfac = Nfac + 1
    write(22,*) 1, 0
    write(22,*) 4, No_bl_1+3, No_bl_1+7, No_bl_1+4, No_bl_1
    No_bl_1 = No_bl_1 + 8
end do



! Holes
!
write(22,*) '# Holes'
write(22,*) 0		! without holes.

! Regions
!
q =-1
q(Nl+1:Nl+N_block)=2.d0
write(22,*) '# Regions'
write(22,*) Nl + 1 + N_block !+ 1
do j = 0, Nl
	write(22,*) j, xi_B+1.d0 , yi_B+1, (z_layers(j-1)+z_layers(j))*0.5d0, j, q(j)
end do
do j = 1, N_block
    write(22,*) Nl+j, x0_bl(j), y0_bl(j), z0_bl(j), Nl+j, q(Nl+j)
end do

print *, 'Nnos, Nnos_tot', Nnos, Nnos_tot
print *, 'Nfac, Nf_tot', Nfac, Nf_tot

!Output of border coordinates.
!
open(44, file=arq_info, status='replace', action='write')
write(44,*) xi_B, xf_B
write(44,*) yi_B, yf_B
write(44,*) zi_B, zf_B
close(44)

end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
