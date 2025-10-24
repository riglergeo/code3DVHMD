!Main code with its subroutine module.
module VFE3D_Hfield
implicit none

integer, parameter	       :: dbl = kind(1.d0)
real(dbl), parameter	  :: pi = 3.1415926535897932384626433832795_dbl 
real(dbl), parameter	  :: eps0 = 8.85418781762039d-12
real(dbl), parameter	  :: mu0 = 4.d-7 * pi
real(dbl), parameter      :: small = 1.d-80
real(dbl), parameter      :: big = 1.0d80
complex(dbl), parameter    :: zero = (0.d0,0.d0)
complex(dbl), parameter    :: ic = (0.d0, 1.d0)
!========variables properties of the environment========================== 
integer 		  :: N_layer, N_block, N_ele, N_edge, N_node, N_front, N_obs
real(dbl)		  ::  moment, z_obs, rt(3)
real(dbl), allocatable	:: x_obs(:), y_obs(:)
real(dbl), allocatable  ::  sig(:,:), coord(:,:), sig0(:)
integer, allocatable	:: nos(:,:), region(:), arestas(:,:), fronteira(:)
real(8), allocatable	:: esp(:), pc(:)
real(8)                 :: w, freq, offset
real(dbl)		      	:: z_rc, h0
real(8), dimension(3)   :: r0
complex(8),parameter    :: ci=(0.d0,1.d0)
complex(dbl)              :: sigma, zeta
real(dbl),Parameter       :: mi0=1.2566370614359173d-6
integer, parameter                 :: N_kk=201
!==============================================================
complex(8)		::  uz, uesp
real(8)                 ::  Fih(2), Fiv(2), senrad, r, kr, dx, dy, dz,xr,zr
real(8)                 ::  dxr, dyr, dx2,dy2, dxy, dxysin, dxycos,ddxr, ddyr 
real(8)			:: cosrad			
integer		        :: cam
integer   		:: navalHx, navalHy, navalHz
!================paramiter=========================
integer, parameter :: sp = kind(1.e0)
integer, parameter :: dp = kind(1.d0)
integer, parameter :: qp = selected_real_kind(30)
real(dp), parameter :: mu = 4.d-7 * pi
real(dp), parameter :: epslon = 8.85d-12
!==================modulo util_nr================================
	integer(4),parameter::npar_arth=16,npar2_arth=8	
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=100
subroutine Ep_dmh(moment, sig0, w, Fi, r0, r1, r2, Ep)
implicit none
real(dbl), intent(in)		:: moment, sig0(0:n_layer), w, Fi(2)
real(dbl), intent(in)		:: r0(3), r1(3), r2(3)
complex(dbl), intent(out)	:: Ep
complex(dbl)			:: s_Ex, s_Ey, s_Ez
integer			:: navalEx, navalEy, navalEz
integer			:: i, j, c
real(dbl)                  	:: norma
real(dbl)			:: x,y,z
real(dbl)                         :: ab(N_kk), pej0(N_kk), pej1(N_kk)
complex(dbl)	      		   :: Ex, Ey, Ez  
complex(dbl)            		   :: Kernel_aux1, Kernel_aux2, kernel_aux3, aux4
complex(dbl)            		   :: kernel_Ez, aux_rte
real(dbl)               		   :: pej0kr
complex(dbl)            		   :: expuzp, expuzn, uz, Asu, rtmexp, uesp
complex(dbl)                   	   :: tanhe
complex(dbl),dimension(0:n_layer)    :: u, Y_int, Z_int, F, A, k_no, E
complex(dbl),dimension(n_layer)      :: Y_ap, Z_ap
complex(dbl),dimension(0:n_layer-1)  :: Rte, Rtm,uh
real(dbl)                      	    :: senrad, r, kr, dx, dy, dz
real(dbl)                      	  ::  dxr, dyr, dx2,dy2, dxy, dxysin, dxycos 
real(dbl)			    :: cosrad, ddxr,ddyr			
integer		    	 :: cam      
call J0J1Wer( ab, pej0, pej1 )
! Midpoint:
   x = (r1(1) + r2(1))*0.5d0
   y = (r1(2) + r2(2))*0.5d0 
   z_rc = (r1(3) + r2(3))*0.5d0
 
   !Distance between source and nodes 1 e 2:
    x = x - r0(1)
    y = y - r0(2)
   h0=r0(3)	

   r=sqrt(x**2+y**2)
   if(r<10.d-2) then
    
      r = 10.d-2
   end if
   cosrad = Fi(1)
   senrad = Fi(2) 
	
    dxr  = (x/r)**2
    dyr  = (y/r)**2
    ddxr  = (x/r)
    ddyr  = (y/r) 
    dx2  = (1.d0 - 2*dxr)/ r
    dy2  = (1.d0 - 2*dyr)/ r
    dxy  = (x*y)/(r**2)
    dxysin = (-2*dxy/r)*senrad
    dxycos = (2*dxy/r)*cosrad
   
   sigma= 1.d-12 + ci*w*8.85*1.d-12
   
   zeta=ci*w*mu0 
  s_Ex = (0.d0, 0.d0)
  s_Ey = (0.d0, 0.d0)
  s_Ez = (0.d0, 0.d0)
Do j=1,N_kk
	kr = (ab(j))/r
	zeta=ci*w*mu0 
	k_no(0) = sqrt(-zeta*sigma)
	k_no(1:n_layer) = sqrt(-zeta*sig0(1:n_layer))


  	u(0:n_layer)=sqrt(kr**2-k_no(0:n_layer)**2) 	 
	Z_int(0)=u(0)/sigma
	Z_int(1:n_layer)=u(1:n_layer)/sig0(1:n_layer)
  	Y_int(0:n_layer)=u(0:n_layer)/zeta
    	Y_ap(n_layer)=Y_int(n_layer) 	
	Z_ap(n_layer)=Z_int(n_layer)
	DO i=n_layer-1,1,-1 
 		tanhe = tanh(u(i)*esp(i)) 
		Y_ap(i)=Y_int(i)*((Y_ap(i+1)+Y_int(i)*tanhe)/(Y_int(i)+Y_ap(i+1)*tanhe))		
		Z_ap(i)=Z_int(i)*((Z_ap(i+1)+Z_int(i)*tanhe)/(Z_int(i)+Z_ap(i+1)*tanhe))
	END DO
 	
	!Reflection coefficient:
	DO i=0,n_layer-1
		Rte(i)=(Y_int(i)-Y_ap(i+1))/(Y_int(i)+Y_ap(i+1))
		Rtm(i)=(Z_int(i)-Z_ap(i+1))/(Z_int(i)+Z_ap(i+1))
	END DO

  	uh(0:n_layer-1)=u(0:n_layer-1)*esp(0:n_layer-1)
	F(0) = exp(u(0)*h0)/(4*pi)
	A(0) = F(0)*k_no(0)**2
	F(1) = F(0)*((1+Rte(0))/(1+Rte(1)*exp(-2*uh(1))))
	A(1) = A(0)*((1+Rtm(0))/(1+Rtm(1)*exp(-2*uh(1))))
	DO i=2,n_layer-1	
		F(i)=F(i-1)*exp(uh(i))*((1+Rte(i-1))*exp(-uh(i-1))/(1+Rte(i)*exp(-2*uh(i))))
		A(i)=A(i-1)*exp(uh(i))*((1+Rtm(i-1))*exp(-uh(i-1))/(1+Rtm(i)*exp(-2*uh(i))))
	END DO
	
	A(n_layer)=A(n_layer-1)*(1+Rtm(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))
	F(n_layer)=F(n_layer-1)*(1+Rte(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))
   IF(z_rc<0.) THEN 
	cam=0			
   Else
	DO c=1,n_layer-1			
		IF(z_rc<pc(c))THEN
			cam=c
			exit
		END IF
	END DO
  		IF(z_rc>=pc(n_layer-1))  cam=n_layer
   End if
     
    IF(cam<n_layer) THEN 
		expuzp = exp(u(cam)*(z_rc-(pc(cam)+ esp(cam))))
		expuzn = exp(-u(cam)*(z_rc-(pc(cam)- esp(cam))))
		Asu    = A(cam)/(sig0(cam)*u(0))
		rtmexp = Rtm(cam)*expuzp			
		Kernel_aux1 =  F(cam)*(expuzn + Rte(cam)*expuzp)*zeta
		Kernel_aux3 =  Asu*u(cam)*(expuzn - rtmexp)
		kernel_Ez   = - Asu*(kr**2)*(expuzn + rtmexp)		
	ELSE 
		uz = u(n_layer)*(z_rc-pc(n_layer-1)) 
		expuzn = exp(-uz)
		Asu    = A(n_layer)/(sig0(n_layer)*u(0))
		kernel_aux1 = F(n_layer)*expuzn*zeta
		kernel_aux3 = Asu*u(n_layer)*expuzn      
		kernel_Ez   = -Asu*((kr)**2)*expuzn						   		
	END IF   
	pej0kr = pej0(j)*kr              	
     Ey =  (kernel_aux1*(pej1(j)*(dx2*cosrad + dxysin) + pej0kr*(dxr*cosrad + &
              senrad*dxy)) + kernel_aux3*(pej1(j)*(-dy2*cosrad + dxysin)+ pej0kr*(-dyr*cosrad &
              + dxy*senrad)))
               
			
	Ex = ((kernel_aux1*(pej1(j)*(dxycos - dy2*senrad) + pej0kr*(-dxy*cosrad -dyr*senrad)) + &
	    kernel_aux3*(pej1(j)*(dxycos + dx2*senrad) + pej0kr*(-dxy*cosrad + dxr*senrad))))			 			
		      
	Ez = ((kernel_Ez*pej1(j)*(cosrad*(-sqrt(dyr)) + senrad*ddxr)))

	S_Ey = S_Ey + Ey 
	S_Ex = S_Ex + Ex
	S_Ez = S_Ez + Ez
 end do
 S_Ey = S_Ey/r
 S_Ex = S_Ex/r
 S_Ez = S_Ez/r

   dx= r2(1)-r1(1)
   dy= r2(2)-r1(2)
   dz= r2(3)-r1(3)
   norma = sqrt( dx**2 + dy**2 + dz**2)
! Projection of the electric field onto the edge: 
  Ep = (s_Ex*dx + s_Ey*dy + s_Ez*dz )/norma

End Subroutine
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=100
subroutine campo_direto_esth(moment,sig0, w, r, Hp )   !campo primario do meio estratificado
implicit none 
real(dbl), intent(in)	 :: moment, w, r
real(dbl), intent(in)      :: sig0(0:N_layer)
integer   		 :: navalHx, navalHy, navalHz
complex(dbl)		 ::s_Hx, s_Hy, s_Hz
complex(dbl), intent(out)  :: Hp(2)
integer		  :: i, j, c
real(dbl)                   :: norma, kr
real(dbl)                   :: y, z, x
real(dbl)                         :: ab(N_kk), pej0(N_kk), pej1(N_kk)
complex(dbl)	      		   :: Hy  
complex(dbl)            		   :: Kernel_aux4, Kernel_aux5
real(dbl)               		   :: pej0kr
complex(dbl)            		   :: expuzp, expuzn, uz, Asu, rtmexp, uesp
complex(dbl)                   	   :: tanhe
complex(dbl),dimension(0:n_layer)    :: u, Y_int, Z_int, F, A, k_no, E
complex(dbl),dimension(n_layer)      :: Y_ap, Z_ap
complex(dbl),dimension(0:n_layer-1)  :: Rte, Rtm,uh

     
call J0J1Wer( ab, pej0, pej1 )
   
   sigma= 1.d-12 + ci*w*8.85*1.d-12
   zeta=ci*w*mi0 

 s_Hy = (0.d0, 0.d0)
Do j=1,N_kk
	kr = (ab(j))/r

	zeta=ci*w*mu0 
	k_no(0) = sqrt(-zeta*sigma)
	k_no(1:n_layer) = sqrt(-zeta*sig0(1:n_layer))

  	u(0:n_layer)=sqrt(kr**2-k_no(0:n_layer)**2) 	 
	Z_int(0)=u(0)/sigma
	Z_int(1:n_layer)=u(1:n_layer)/sig0(1:n_layer)
  	Y_int(0:n_layer)=u(0:n_layer)/zeta
    	Y_ap(n_layer)=Y_int(n_layer) 	
	Z_ap(n_layer)=Z_int(n_layer)
	DO i=n_layer-1,1,-1 
 		tanhe = tanh(u(i)*esp(i)) 
		Y_ap(i)=Y_int(i)*((Y_ap(i+1)+Y_int(i)*tanhe)/(Y_int(i)+Y_ap(i+1)*tanhe))		
		Z_ap(i)=Z_int(i)*((Z_ap(i+1)+Z_int(i)*tanhe)/(Z_int(i)+Z_ap(i+1)*tanhe))
	END DO
 	
	! Reflection coefficient:
	DO i=0,n_layer-1
		Rte(i)=(Y_int(i)-Y_ap(i+1))/(Y_int(i)+Y_ap(i+1))
		Rtm(i)=(Z_int(i)-Z_ap(i+1))/(Z_int(i)+Z_ap(i+1))
	END DO

  	uh(0:n_layer-1)=u(0:n_layer-1)*esp(0:n_layer-1)
	F(0) = exp(u(0)*h0)/(4*pi)
	A(0) = F(0)*k_no(0)**2
	F(1) = F(0)*((1+Rte(0))/(1+Rte(1)*exp(-2*uh(1))))
	A(1) = A(0)*((1+Rtm(0))/(1+Rtm(1)*exp(-2*uh(1))))
	DO i=2,n_layer-1	
		F(i)=F(i-1)*exp(uh(i))*((1+Rte(i-1))*exp(-uh(i-1))/(1+Rte(i)*exp(-2*uh(i))))
		A(i)=A(i-1)*exp(uh(i))*((1+Rtm(i-1))*exp(-uh(i-1))/(1+Rtm(i)*exp(-2*uh(i))))
	END DO
	
	A(n_layer)=A(n_layer-1)*(1+Rtm(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))
	F(n_layer)=F(n_layer-1)*(1+Rte(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))
                    	
 	expuzn = exp(-u(0)*(-0.5d0))
	Kernel_aux4 =  F(0)*(expuzn - Rte(0)*exp(u(0)*(-0.5d0)))*u(0)		
	Kernel_aux5 =  A(0)*(expuzn + Rtm(0)*exp(u(0)*(-0.5d0)))/u(0) 
	Hy = (kernel_aux4*(-(1.d0/r))*pej1(j) + kernel_aux5*(((-1.d0)/r)*pej1(j) + (pej0(j)*kr)))                     	
	      
    S_Hy = S_Hy + Hy          
 end do
 
 S_Hy = S_Hy/r

   Hp(1) = S_Hx
   Hp(2) = S_Hy	
End Subroutine
!======================================================================================
subroutine leituras()
implicit none
character(50):: arq_edge, arq_node, arq_ele
character(50):: arq_in, arq_arestas, arq_Nnn, arq_front
character(50):: modelo
character(5):: Header
integer:: i, j, k, temp, Nregions
real:: lixo, sigma0
call GET_COMMAND_ARGUMENT(1, modelo)

k = len_trim(modelo)

arq_in =   modelo(1:k)//'.in'
arq_edge = modelo(1:k)//'.1.edge'   !file with number of the edges and yours nodes
arq_node = modelo(1:k)//'.1.node'
arq_ele = modelo(1:k)//'.1.ele'
arq_arestas = modelo(1:k)//'.1.nedge' ! file with number of elements and yours edges
arq_front = modelo(1:k)//'_vfe.front'

! reading the model
!!
open(11, file=arq_in, status='old', action='read')


do
	read(11,*) Header
	if (Header=='Block') exit
end do
read(11,*) N_block
rewind(11)

do
	read(11,*) Header
	if (header=='Backg') then
       		read(11,*) sigma0, lixo
		read(11,*) N_layer
        	Nregions = N_layer + N_block
        	allocate( sig0(0:n_layer), esp(0:N_layer), sig(Nregions, 2))		
		sig0(0) = 1.d0/sigma0		
		esp(0) = 0.d0
        	do i = 1, N_layer
			read(11,*) sig(i,1:2), esp(i)
			sig0(i) = sig(i,1)
			sig0(i) = 1./sig0(i)		
		end do
		exit
	end if
end do
	
allocate(pc(0:n_layer-1))

   ! layer interface:
    pc(0) = 0.d0
    pc(1) = esp(1)
    DO i = 2, N_layer-1
	pc(i) = pc(i-1) + esp(i)
   END DO
!Reading the blocks
!
rewind(11)
do
	read(11,*) Header
	if (Header=='Block') exit
end do
read(11,*) N_block
do i = 1, N_block
    read(11,*)
    read(11,*) 
    read(11,*) sig(N_layer+i,1:2) 
end do
sig = 1.d0 / sig

rewind(11)
do
	read(11,*) Header
	if (Header=='Perfi') exit
end do
read(11,*) N_obs
read(11,*)
read(11,*)
read(11,*) z_obs
read(11,*) rt
read(11,*) moment
read(11,*) offset, freq
close(11)

open(22,file='Malha/pontos_obs.dat', status='old')
read(22,*) N_obs
allocate( x_obs(N_obs), y_obs(N_obs) )
do i = 1, N_obs
    read(22,*) x_obs(i), y_obs(i)
end do
close(22)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reading nodes by element.
!
open(22, file=arq_ele, status='old', action='read')
read(22,*) N_ele
allocate( Nos(4,N_ele), region(N_ele) ) !each element is one region in the layer or in the block
do i = 1, N_ele
	read(22,*) temp, Nos(:,i), region(i) !fill of the vector nos(nodes of the element, number of the element)
										 !the vector region (attribute 0,1,2) for each iézimo element
end do
close(22)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reading node coordinates.
!
open(33, file=arq_node, status='old', action='read')
read(33,*) N_node
allocate( coord(3,N_node) )
do j = 1, N_node
	read(33,*) k, (coord(i,j),i=1,3)
end do
close(33)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Reading edges by element
!
open(44, file=arq_arestas, status='old', action='read')
read(44,*) N_ele
allocate( arestas(6,N_ele) )
do i = 1, N_ele
	read(44,*) arestas(:,i)
end do
close(44)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(55, file=arq_edge, status='old', action='read')
read(55,*) N_edge
close(55)

open(66, file=arq_front, status='old', action='read')
read(66,*) N_front
allocate( fronteira(N_front) )
do i = 1, N_front
	read(66,*) fronteira(i)
end do
close(66)


end subroutine leituras

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Cond_de_front_fonte( N_edge, N_obs, N_front, fronteira, V_fonte )
implicit none
integer, intent(in):: N_edge, N_obs, N_front
integer, intent(in):: fronteira(N_front)
complex(dbl), intent(inout):: V_fonte(N_edge, N_obs)
integer:: i

do i = 1, N_front
	V_fonte( fronteira(i), : ) = zero
end do

end subroutine Cond_de_front_fonte
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Cond_de_front( Nnn, N_edge, N_front, ia, fronteira, M_global )
implicit none
integer, intent(in):: Nnn, N_edge, N_front
integer, intent(in):: ia(N_edge+1), fronteira(N_front)
complex(dbl), intent(inout):: M_global(Nnn)
real(dbl), parameter:: big = 1.0d80
integer:: i

do i = 1, N_front
	M_global( ia(fronteira(i)) ) = big
end do

end subroutine Cond_de_front

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Montagem_final( N_edge, N_ele, Nn_estimate, Nnn, ja_temp, global_temp, ia, ja, M_global )
implicit none
integer, intent(in):: N_edge, N_ele, Nn_estimate, Nnn
integer, intent(in):: ja_temp(Nn_estimate, N_edge)
complex(dbl), intent(in):: global_temp(Nn_estimate, N_edge)
integer, intent(inout):: ia(N_edge+1), ja(Nnn)
complex(dbl), intent(inout):: M_global(Nnn)
!
integer:: i, Temp_1, Temp_2

Temp_1 = 1
Temp_2 = 0
do i = 1, N_edge
	Temp_2 = Temp_2 + ia(i)
	ja(Temp_1 : Temp_2) = ja_temp(1:ia(i),i)
	M_global(Temp_1 : Temp_2) = global_temp(1:ia(i),i)
	Temp_1 = Temp_2 + 1
end do

Temp_1 = ia(1)
ia(1) = 1
do i = 2, N_edge+1
	Temp_2 = ia(i)
	ia(i) = Temp_1 + ia(i-1)
	Temp_1 = Temp_2
end do

end subroutine Montagem_final
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Remove_zeros( N_edge, N_ele, Nn_estimate, ia, ja_temp, global_temp )
implicit none
integer, intent(in):: N_edge, N_ele, Nn_estimate
integer, intent(inout):: ia(N_edge+1), ja_temp(Nn_estimate, N_edge)
complex(dbl), intent(inout):: global_temp(Nn_estimate, N_edge)
real(dbl), parameter:: small = 1.d-80
integer:: i, j, k

do i = 1, N_edge
	j = 0
	do
		j = j + 1
		if ( abs(global_temp(j,i)) < small ) then
			if ( j<ia(i) ) then
				do k = j, ia(i)
					global_temp(k,i) = global_temp(k+1,i)
					ja_temp(k,i) = ja_temp(k+1,i)
				end do
				ia(i) = ia(i) - 1
			else
				ia(i) = ia(i) - 1
				exit
			end if
		end if
		if (j == ia(i)) exit
	end do
end do

end subroutine Remove_zeros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Monta_colunas( N_edge, Nn_estimate, Mel, ar, ia, ja_temp, Global_temp )
!mounting of the tempory global matriz 
! 
implicit none

integer, intent(in):: N_edge, Nn_estimate
complex(dbl), intent(in):: Mel(6,6)
integer, intent(in):: ar(6)
integer, intent(inout):: ia(N_edge+1), ja_temp(Nn_estimate, N_edge)
complex(dbl), intent(inout):: global_temp(Nn_estimate, N_edge)
!
integer:: i, j, k, pos, Temp_1, Temp_2
complex(dbl):: Temp_gl_1, Temp_gl_2
integer:: lin, col

do i = 1, 6
	do j = i, 6
		k = 1
		if (ar(i)<=ar(j)) then
			lin = ar(i)
			col = ar(j)
		else
			lin = ar(j)
			col = ar(i)
		end if
			do
				if (ja_temp(k,lin) == 0) then
					ja_temp(k,lin) = col
					ia(lin) = ia(lin) + 1
					Global_temp(k,lin) = Mel(i,j)
					exit
				else if ( col == ja_temp(k,lin) ) then
					Global_temp(k,lin) = Global_temp(k,lin) + Mel(i,j)
					exit
				else if ( col < ja_temp(k,lin) ) then
					pos = k
					Temp_2 = ja_temp(pos,lin)
					ja_temp(pos,lin) = col
					ia(lin) = ia(lin) + 1
					Temp_gl_2 = Global_temp(pos,lin)
					Global_temp(pos,lin) = Mel(i,j)
					do
						pos = pos + 1
						if ( ja_temp(pos,lin) == 0 ) then
							ja_temp(pos,lin) = Temp_2
							Global_temp(pos,lin) = Temp_gl_2
							exit
						else
							Temp_1 = ja_temp(pos,lin)
							ja_temp(pos,lin) = Temp_2
							Temp_2 = Temp_1

							Temp_gl_1 = Global_temp(pos,lin)
							Global_temp(pos,lin) = Temp_gl_2
							Temp_gl_2 = Temp_gl_1
						end if
					end do
					exit
				else 
					k = k + 1
				end if
			end do
	end do
end do
	
end subroutine Monta_colunas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Matriz_el( node, x, y, z, reg, zet0, eta, D_eta, Mel, G )
implicit none
integer, intent(in):: reg, node(4)
real(dbl), intent(in):: x(4), y(4), z(4)
complex(dbl), intent(in):: zet0, eta(2), D_eta(2)
complex(dbl), intent(out):: Mel(6,6), G(6,6)
integer:: i, j, k
real(dbl):: z12, z13, z14, z23, z24, z34
real(dbl):: y12, y13, y14, y23, y24, y34
real(dbl):: TV_E, l_V720(6), V360, TV_G(6)
real(dbl):: fm(4,4)
real(dbl):: a(4), b(4), c(4), d(4), V, l(6)
real(dbl):: cd(6), db(6), bc(6), bb(6), cc(6), dd(6)
real(dbl):: F(6,6)
complex(dbl):: E(6,6) 


z12 = z(1) - z(2)
z13 = z(1) - z(3)
z14 = z(1) - z(4)
z23 = z(2) - z(3)
z24 = z(2) - z(4)
z34 = z(3) - z(4)

y12 = y(1) - y(2)
y13 = y(1) - y(3)
y14 = y(1) - y(4)
y23 = y(2) - y(3)
y24 = y(2) - y(4)
y34 = y(3) - y(4)

l(1) = dsqrt( (x(2)-x(1))**2 + y12**2 + z12**2 )
l(2) = dsqrt( (x(3)-x(1))**2 + y13**2 + z13**2 )
l(3) = dsqrt( (x(4)-x(1))**2 + y14**2 + z14**2 )
l(4) = dsqrt( (x(3)-x(2))**2 + y23**2 + z23**2 )
l(5) = dsqrt( (x(4)-x(2))**2 + y24**2 + z24**2 )
l(6) = dsqrt( (x(4)-x(3))**2 + y34**2 + z34**2 )

k=0
do i=1,3
	do j=i+1,4
		k=k+1
		if (node(i)>node(j)) then
			l(k) = -l(k)
		end if
	end do
end do

!!!!!!!!!!!!!!!!!!!!!! V -> Element volume.
V = abs( x(1)*( -y(2)*z34 + y(3)*z24 - y(4)*z23 ) + &
		  x(2)*(  y(1)*z34 - y(3)*z14 + y(4)*z13 ) + &
		  x(3)*( -y(1)*z24 + y(2)*z14 - y(4)*z12 ) + &
		  x(4)*(  y(1)*z23 - y(2)*z13 + y(3)*z12 ) ) / 6.d0

b(1) = -y(2)*z34 + y(3)*z24 - y(4)*z23
b(2) =  y(1)*z34 - y(3)*z14 + y(4)*z13
b(3) = -y(1)*z24 + y(2)*z14 - y(4)*z12
b(4) =  y(1)*z23 - y(2)*z13 + y(3)*z12

 c(1) =  x(2)*z34 - x(3)*z24 + x(4)*z23
 c(2) = -x(1)*z34 + x(3)*z14 - x(4)*z13
 c(3) =  x(1)*z24 - x(2)*z14 + x(4)*z12
 c(4) = -x(1)*z23 + x(2)*z13 - x(3)*z12

d(1) = -x(2)*(y34) + x(3)*(y24) - x(4)*(y23)
d(2) =  x(1)*(y34) - x(3)*(y14) + x(4)*(y13)
d(3) = -x(1)*(y24) + x(2)*(y14) - x(4)*(y12)
d(4) =  x(1)*(y23) - x(2)*(y13) + x(3)*(y12)

!!!!!!!!!!!!!!!!!!!!!! Eij.
 cd(1) = c(1)*d(2) - d(1)*c(2)
 cd(2) = c(1)*d(3) - d(1)*c(3)
 cd(3) = c(1)*d(4) - d(1)*c(4)
 cd(4) = c(2)*d(3) - d(2)*c(3)
 cd(5) = c(2)*d(4) - d(2)*c(4)
 cd(6) = c(3)*d(4) - d(3)*c(4)

db(1) = d(1)*b(2) - b(1)*d(2)
db(2) = d(1)*b(3) - b(1)*d(3)
db(3) = d(1)*b(4) - b(1)*d(4)
db(4) = d(2)*b(3) - b(2)*d(3)
db(5) = d(2)*b(4) - b(2)*d(4)
db(6) = d(3)*b(4) - b(3)*d(4)

bc(1) = b(1)*c(2) - c(1)*b(2)
bc(2) = b(1)*c(3) - c(1)*b(3)
bc(3) = b(1)*c(4) - c(1)*b(4)
bc(4) = b(2)*c(3) - c(2)*b(3)
bc(5) = b(2)*c(4) - c(2)*b(4)
bc(6) = b(3)*c(4) - c(3)*b(4)

TV_E = 324.d0*(V**3)
do i = 1, 6
	do j = i, 6
		E(i,j) = l(i)*l(j) * ( (cd(i)*cd(j) + db(i)*db(j))/eta(1) + bc(i)*bc(j)/eta(2) ) / TV_E
!		E(j,i) = E(i,j)
	end do
end do

!!!!!!!!!!!!!!!!!!!!!! Fij.
do i = 1, 4
	do j = i, 4
		fm(i,j) = b(i)*b(j) + c(i)*c(j) + d(i)*d(j)
	end do
end do
v360 = 360.d0*V
l_V720 = l / (720.d0*V)

F(1,1) =  l(1)**2 / V360 * ( fm(2,2) - fm(1,2) + fm(1,1) ) 
F(1,2) =  l_v720(1)*l(2) * ( 2*fm(2,3) - fm(1,2) - fm(1,3) + fm(1,1) )
F(1,3) =  l_v720(1)*l(3) * ( 2*fm(2,4) - fm(1,2) - fm(1,4) + fm(1,1) )
F(1,4) =  l_v720(1)*l(4) * ( fm(2,3) - fm(2,2) - 2*fm(1,3) + fm(1,2) )
F(1,5) = -l_v720(1)*l(5) * ( fm(2,2) - fm(2,4) - fm(1,2) + 2*fm(1,4) )
F(1,6) =  l_v720(1)*l(6) * ( fm(2,4) - fm(2,3) - fm(1,4) + fm(1,3)   )
F(2,2) =  l(2)**2 / V360 * ( fm(3,3) - fm(1,3) + fm(1,1) )
F(2,3) =  l_v720(2)*l(3) * ( 2*fm(3,4) - fm(1,3) - fm(1,4) + fm(1,1) )
F(2,4) =  l_v720(2)*l(4) * ( fm(3,3) - fm(2,3) - fm(1,3) + 2*fm(1,2) )
F(2,5) = -l_v720(2)*l(5) * ( fm(2,3) - fm(3,4) - fm(1,2) + fm(1,4)   )
F(2,6) =  l_v720(2)*l(6) * ( fm(1,3) - fm(3,3) - 2*fm(1,4) + fm(3,4) )
F(3,3) =  l(3)**2 / V360 * ( fm(4,4) - fm(1,4) + fm(1,1) )
F(3,4) =  l_v720(3)*l(4) * ( fm(3,4) - fm(2,4) - fm(1,3) + fm(1,2)   )
F(3,5) = -l_v720(3)*l(5) * ( fm(2,4) - fm(4,4) - 2*fm(1,2) + fm(1,4) )
F(3,6) =  l_v720(3)*l(6) * ( fm(4,4) - fm(3,4) - fm(1,4) + 2*fm(1,3) )
F(4,4) =  l(4)**2 / V360 * ( fm(3,3) - fm(2,3) + fm(2,2) )
F(4,5) = -l_v720(4)*l(5) * ( fm(2,3) - 2*fm(3,4) - fm(2,2) + fm(2,4) )
F(4,6) =  l_v720(4)*l(6) * ( fm(3,4) - fm(3,3) - 2*fm(2,4) + fm(2,3) )
F(5,5) =  l(5)**2 / V360 * ( fm(2,2) - fm(2,4) + fm(4,4) )
F(5,6) = -l_v720(5)*l(6) * ( fm(2,4) - 2*fm(2,3) - fm(4,4) + fm(3,4) )
F(6,6) =  l(6)**2 / V360 * ( fm(4,4) - fm(3,4) + fm(3,3) )

Mel = E + zet0*F

!!!!!!!!!!!!!!!!!!!!! source vector
if (reg > N_layer) then
	bb(1) = b(2) - b(1)
	bb(2) = b(3) - b(1)
	bb(3) = b(4) - b(1)
	bb(4) = b(3) - b(2)
	bb(5) = b(4) - b(2)
	bb(6) = b(4) - b(3)

	cc(1) = c(2) - c(1)
	cc(2) = c(3) - c(1)
	cc(3) = c(4) - c(1)
	cc(4) = c(3) - c(2)
	cc(5) = c(4) - c(2)
	cc(6) = c(4) - c(3)

	dd(1) = d(2) - d(1)
	dd(2) = d(3) - d(1)
	dd(3) = d(4) - d(1)
	dd(4) = d(3) - d(2)
	dd(5) = d(4) - d(2)
	dd(6) = d(4) - d(3)

	TV_G = l / (432.d0 * V**2)
	do i = 1, 6
		do j = 1, 6
			G(i,j) = TV_G(j)*l(i) * ( (bb(j)*cd(i) + cc(j)*db(i))*D_eta(1)/eta(1) + dd(j)*bc(i)*D_eta(2)/eta(2) )
		end do
	end do

end if

end subroutine Matriz_el

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine campo_primarioh( N_edge, node, x, y, z, rt, m, sig0, w, Fi, ar, E_prim )
implicit none
integer, intent(in):: N_edge, node(4), ar(6)
real(dbl), intent(in):: x(4), y(4), z(4), rt(3), m, w, Fi(2)
real(dbl), intent(in):: sig0(0:N_layer)
real(dbl):: r1(3), r2(3)
complex(dbl), intent(inout):: E_prim(N_edge)
integer:: i, j, k

k = 0
do i = 1, 3
	do j = i+1, 4
		if (node(i)<node(j)) then
			r1(1) = x(i)
			r1(2) = y(i)
			r1(3) = z(i)
			r2(1) = x(j)
			r2(2) = y(j)
			r2(3) = z(j)
		else
			r1(1) = x(j)
			r1(2) = y(j)
			r1(3) = z(j)
			r2(1) = x(i)
			r2(2) = y(i)
			r2(3) = z(i)
		end if
	
		k = k + 1
        if (E_prim( ar(k) ) == (0.d0,0.d0) ) then
	    	
	     call Ep_dmh( moment, sig0, w, Fi, rt, r1, r2, E_prim(ar(k)))
		
        end if
	end do
end do

end subroutine campo_primarioh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Ep_dmxy( m, s0, w, Fi, r0, r1, r2, Ep )
!
! Projection of the electric field from a magnetic dipole in the (x,y) plane onto an 
! arbitrary edge defined by its node coordinates.
!
! m  	 -> Input: Dipole moment.
! s0     -> Primary medium conductivity.
! w      -> 2*pi*freq
! Fi 	 -> Input: Fi(1)=cos(Dip); Fi(2)=sin(Dip). Dip: Angle of the moment vector relative to the x axis.
! r0 	 -> Input: vector with the coordinates of the source.
! r1, r2 -> Input: vector with the oordinates of the nodes that define the edge.
! Ep 	 -> Output: Projection of the electric field at the center of the edge.
! Global constants:
!   pi, mu0, eps0
!
implicit none
real(dbl), intent(in):: m, s0, w, Fi(2), r0(3), r1(3), r2(3)
complex(dbl), intent(out):: Ep
real(dbl):: x1, x2, y1, y2, z1, z2, l, R
complex(dbl):: ikr, e0

x1 = r1(1) - r0(1)
x2 = r2(1) - r0(1)

y1 = r1(2) - r0(2)
y2 = r2(2) - r0(2)

z1 = r1(3) - r0(3)
z2 = r2(3) - r0(3)

l = dsqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )            ! Length of the edge.
R = dsqrt( (x1+x2)**2 + (y1+y2)**2 + (z1+z2)**2 ) / 2.d0     ! Radial distance from the source to the center of the edge.

ikr = cdsqrt( cmplx( -w**2 * mu0 * eps0, w * mu0 * s0, dbl ) ) * R

E0 = m * (cmplx(0.d0, w*mu0, dbl)) * (1.d0 + ikr) * cdexp(-ikr) / (4.d0*pi*l*R**3)

Ep = E0 * ( (z1*y2 - z2*y1)*Fi(1) + (x1*z2 - x2*z1)*Fi(2) )

end subroutine Ep_dmxy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcula_obs( node, x,y,z, H, xr, yr, zr, fb_v )
implicit none
integer, intent(in):: node(4)
real(dbl), intent(in):: x(4), y(4), z(4)
real(dbl), intent(in):: xr, yr, zr
complex(dbl), intent(in):: H(6)
!complex(dbl), intent(out):: obs
integer:: i, j, k
real(dbl):: z12, z13, z14, z23, z24, z34
real(dbl):: y12, y13, y14, y23, y24, y34
real(dbl):: a(4), b(4), c(4), d(4), V, l(6)
real(dbl):: fb_n(4) ! FBN: Nodal base function.
complex(dbl):: fb_v(3) ! FBV: Vector basis function.

z12 = z(1) - z(2)
z13 = z(1) - z(3)
z14 = z(1) - z(4)
z23 = z(2) - z(3)
z24 = z(2) - z(4)
z34 = z(3) - z(4)

y12 = y(1) - y(2)
y13 = y(1) - y(3)
y14 = y(1) - y(4)
y23 = y(2) - y(3)
y24 = y(2) - y(4)
y34 = y(3) - y(4)

l(1) = dsqrt( (x(2)-x(1))**2 + y12**2 + z12**2 )
l(2) = dsqrt( (x(3)-x(1))**2 + y13**2 + z13**2 )
l(3) = dsqrt( (x(4)-x(1))**2 + y14**2 + z14**2 )
l(4) = dsqrt( (x(3)-x(2))**2 + y23**2 + z23**2 )
l(5) = dsqrt( (x(4)-x(2))**2 + y24**2 + z24**2 )
l(6) = dsqrt( (x(4)-x(3))**2 + y34**2 + z34**2 )

!!!!!!!!!!!!!!!!!!!!!! V -> Element volume.
V = abs( x(1)*( -y(2)*z34 + y(3)*z24 - y(4)*z23 ) + &
			x(2)*(  y(1)*z34 - y(3)*z14 + y(4)*z13 ) + &
			x(3)*( -y(1)*z24 + y(2)*z14 - y(4)*z12 ) + &
			x(4)*(  y(1)*z23 - y(2)*z13 + y(3)*z12 ) ) / 6.d0

!!!!!!!!!!!!!!!!!!!!!! coefficients a, b, c, d.

a(1) = x(2)*( y(3)*z(4)-y(4)*z(3) ) + x(3)*( y(4)*z(2)-y(2)*z(4) ) + x(4)*( y(2)*z(3)-y(3)*z(2) )
a(2) = x(3)*( y(1)*z(4)-y(4)*z(1) ) + x(4)*( y(3)*z(1)-y(1)*z(3) ) + x(1)*( y(4)*z(3)-y(3)*z(4) )
a(3) = x(4)*( y(1)*z(2)-y(2)*z(1) ) + x(1)*( y(2)*z(4)-y(4)*z(2) ) + x(2)*( y(4)*z(1)-y(1)*z(4) )
a(4) = x(1)*( y(3)*z(2)-y(2)*z(3) ) + x(2)*( y(1)*z(3)-y(3)*z(1) ) + x(3)*( y(2)*z(1)-y(1)*z(2) )

b(1) = -y(2)*z34 + y(3)*z24 - y(4)*z23
b(2) =  y(1)*z34 - y(3)*z14 + y(4)*z13
b(3) = -y(1)*z24 + y(2)*z14 - y(4)*z12
b(4) =  y(1)*z23 - y(2)*z13 + y(3)*z12

 c(1) =  x(2)*z34 - x(3)*z24 + x(4)*z23
 c(2) = -x(1)*z34 + x(3)*z14 - x(4)*z13
 c(3) =  x(1)*z24 - x(2)*z14 + x(4)*z12
 c(4) = -x(1)*z23 + x(2)*z13 - x(3)*z12

d(1) = -x(2)*(y34) + x(3)*(y24) - x(4)*(y23)
d(2) =  x(1)*(y34) - x(3)*(y14) + x(4)*(y13)
d(3) = -x(1)*(y24) + x(2)*(y14) - x(4)*(y12)
d(4) =  x(1)*(y23) - x(2)*(y13) + x(3)*(y12)

do i = 1, 4
	fb_n(i) = ( a(i) + b(i)*xr + c(i)*yr + d(i)*zr ) / (36 * V**2)
	if (fb_n(i) < 0.d0) print *, 'ALERTA!!',fb_n(i)
end do

k = 0
fb_v = (0.d0,0.d0)
do i = 1, 3
	do j = i+1, 4
		k = k+1
		if (node(i)>node(j)) l(k) = -l(k)
		fb_v(1) = fb_v(1) + ( fb_n(i)*b(j) - fb_n(j)*b(i) ) * l(k) * H(k)
		fb_v(2) = fb_v(2) + ( fb_n(i)*c(j) - fb_n(j)*c(i) ) * l(k) * H(k)
		fb_v(3) = fb_v(3) + ( fb_n(i)*d(j) - fb_n(j)*d(i) ) * l(k) * H(k)
	end do
end do

end subroutine calcula_obs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function localiza_obs(A,B,C,D,P) result(isornot)
implicit none
real(dbl), dimension(3), intent(in) :: A, B, C, D, P
integer :: isornot
real(dbl) :: u(3), v(3), w(3), n(3), AP(3), sinal
u = B - A
v = C - A
w = D - A
AP= P - A
n = (/ u(2)*v(3)-u(3)*v(2), u(3)*v(1)-u(1)*v(3), u(1)*v(2)-u(2)*v(1) /) !curl of u x v
sinal = dot_product(w,n) * dot_product(AP,n) !inner product of vectors.
if (sinal < 0.d0) then
	isornot = 0
else
	isornot = 1
end if
end function localiza_obs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Campo_direto( m, s0, w, r, Hp )
!
! primary field: DMV, Hz.
!
implicit none
real(dbl), intent(in):: m, s0(0:N_layer), w, r
complex(dbl), intent(out):: Hp(2)
complex(dbl):: ikr

ikr = cdsqrt( cmplx( -w**2 * mu0 * eps0, w * mu0 * s0(0), dbl ) ) * r
Hp(1) =  m * cdexp(-ikr) * (1.d0+ikr) / (2.d0*pi*r**3)     ! Coaxial
Hp(2) = -m * cdexp(-ikr) * (ikr**2 + ikr + 1.d0) / (4.d0*pi*r**3)   ! Coplanar

end subroutine Campo_direto
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine J0J1Wer( absc, wJ0, wJ1 )
    !	Be careful! The values ​​identified as abscissas are actually the j values of exp(s1+(i-1)*del) ou 10**(s1+(i-1)*del)
    !	(depending on the technique used by the filter creator)
        real(dbl), intent(out) :: absc( 201 ), wJ0( 201 ), wJ1( 201 )
    
       absc = (/8.653980893285999343D-04, 9.170399868578506730D-04, 9.717635708540675208D-04, 1.029752738345341345D-03, &
                     1.091202360259058363D-03, 1.156318936280037154D-03, 1.225321288786770891D-03, 1.298441298197734088D-03, &
                     1.375924682198844014D-03, 1.458031821470676028D-03, 1.545038634690233410D-03, 1.637237505747725754D-03, &
                     1.734938266294200121D-03, 1.838469236921894991D-03, 1.948178330476121357D-03, 2.064434221206373505D-03, &
                     2.187627583685520637D-03, 2.318172405660454734D-03, 2.456507379246002792D-03, 2.603097375137131825D-03, &
                     2.758435004793547123D-03, 2.923042275846299189D-03, 3.097472346289414230D-03, 3.282311383351382318D-03, &
                     3.478180533293271509D-03, 3.685738008752825322D-03, 3.905681300649099623D-03, 4.138749522080598965D-03, &
                     4.385725892093550113D-03, 4.647440367666977670D-03, 4.924772432759196537D-03, 5.218654053788364042D-03, &
                     5.530072811478751842D-03, 5.860075219597365645D-03, 6.209770241733286907D-03, 6.580333017937914017D-03, &
                     6.973008813749223544D-03, 7.389117204870794542D-03, 7.830056511567909730D-03, 8.297308497682591086D-03, &
                     8.792443350058264107D-03, 9.317124955107483966D-03, 9.873116490254248145D-03, 1.046228634904101083D-02, &
                     1.108661441981135046D-02, 1.174819873906770597D-02, 1.244926254186268406D-02, 1.319216173291635694D-02, &
                     1.397939280356631266D-02, 1.481360122115480543D-02, 1.569759031904567614D-02, 1.663433071714528685D-02, &
                     1.762697030458526201D-02, 1.867884481811322647D-02, 1.979348905174003331D-02, 2.097464873531334692D-02, &
                     2.222629312193481060D-02, 2.355262832652090660D-02, 2.495811146033087569D-02, 2.644746560896089893D-02, &
                     2.802569570413705746D-02, 2.969810534264444649D-02, 3.147031460891126092D-02, 3.334827896114080786D-02, &
                     3.533830924445728605D-02, 3.744709289831891358D-02, 3.968171642946591998D-02, 4.204968922592208086D-02, &
                     4.455896879207719985D-02, 4.721798748965116282D-02, 5.003568087440296575D-02, 5.302151772380819111D-02, &
                     5.618553185661335353D-02, 5.953835585119460899D-02, 6.309125677603129312D-02, 6.685617405236506106D-02, &
                     7.084575957628110043D-02, 7.507342023504076645D-02, 7.955336296053971967D-02, 8.430064247129397115D-02, &
                     8.933121186338736919D-02, 9.466197622039211612D-02, 1.003108494224146802D-01, 1.062968143451746283D-01, &
                     1.126399866514113390D-01, 1.193616823889895873D-01, 1.264844896228644322D-01, 1.340323443416230609D-01, &
                     1.420306108936851552D-01, 1.505061672234652703D-01, 1.594874951939309615D-01, 1.690047762990831426D-01, &
                     1.790899930879975566D-01, 1.897770366412597776D-01, 2.011018204609656135D-01, 2.131024011570106513D-01, &
                     2.258191063352318062D-01, 2.392946701171655144D-01, 2.535743767468323084D-01, 2.687062127671349665D-01, &
                     2.847410282772538936D-01, 3.017327078129410922D-01, 3.197383514239506286D-01, 3.388184665571116749D-01, &
                     3.590371713898611872D-01, 3.804624102975323052D-01, 4.031661821784713884D-01, 4.272247824042622599D-01, &
                     4.527190592081252185D-01, 4.797346853730769523D-01, 5.083624461328500876D-01, 5.386985442530570767D-01, &
                     5.708449233178135573D-01, 6.049096103082166609D-01, 6.410070786239048246D-01, 6.792586327676195523D-01, &
                     7.197928159854947161D-01, 7.627458422329327359D-01, 8.082620539176789132D-01, 8.564944069583262376D-01, &
                     9.076049847882753374D-01, 9.617655430324436594D-01, 1.019158086687097953D+00, 1.079975481742401877D+00, &
                     1.144422103303022853D+00, 1.212714522384784388D+00, 1.285082233695329812D+00, 1.361768426844476521D+00, &
                     1.443030803575896748D+00, 1.529142443766406734D+00, 1.620392723103028176D+00, 1.717088285521651159D+00, &
                     1.819554073675149652D+00, 1.928134420893808709D+00, 2.043194208307565596D+00, 2.165120091018537973D+00, &
                     2.294321797444361266D+00, 2.431233506198735572D+00, 2.576315305136151590D+00, 2.730054737463883274D+00, &
                     2.892968440116886253D+00, 3.065603879901345419D+00, 3.248541193241110570D+00, 3.442395135709440002D+00, &
                     3.647817147897402634D+00, 3.865497544561224075D+00, 4.096167834405143537D+00, 4.340603178295351583D+00, &
                     4.599624994165753655D+00, 4.874103717369282940D+00, 5.164961725750846000D+00, 5.473176439271492555D+00, &
                     5.799783604600056819D+00, 6.145880775710010013D+00, 6.512631002177978523D+00, 6.901266737578352739D+00, &
                     7.313093981108034214D+00, 7.749496666359123154D+00, 8.211941311987910552D+00, 8.701981949908562441D+00, &
                     9.221265347572652260D+00, 9.771536541883744320D+00, 1.035464470334362019D+01, 1.097254935013656407D+01, &
                     1.162732693303372145D+01, 1.232117781324618733D+01, 1.305643365667540934D+01, 1.383556526940939513D+01, &
                     1.466119090079532228D+01, 1.553608504199116247D+01, 1.646318774956322173D+01, 1.744561452546165370D+01, &
                     1.848666678657500739D+01, 1.958984295904654616D+01, 2.075885023463468926D+01, 2.199761702862400270D+01, &
                     2.331030618115176267D+01, 2.470132894631230158D+01, 2.617535981604930129D+01, 2.773735222865149197D+01, &
                     2.939255521463919862D+01, 3.114653103598043415D+01, 3.300517387791185087D+01, 3.497472965617872376D+01, &
                     3.706181700625468523D+01, 3.927344952507589682D+01, 4.161705934003131091D+01, 4.410052208441296528D+01, &
                     4.673218336325475519D+01, 4.952088679849771324D+01, 5.247600374772711262D+01, 5.560746479634944706D+01, &
                     5.892579312903905020D+01, 6.244213989259677078D+01, 6.616832166905808776D+01, 7.011686018497628936D+01, &
                     7.430102439032442874D+01, 7.873487504841911289D+01, 8.343331198671111792D+01, 8.841212416722630962D+01, &
                     9.368804274491759543D+01/)
    
        wJ0 = (/2.940900904253498815D+00,-1.601154970027019786D+01, 3.574488144287594338D+01,-3.775710631592443178D+01, &
                       9.347313619702582344D+00, 1.903333998229986435D+01,-1.266160545445113073D+01,-1.274822633210828471D+01, &
                       1.499040669570122830D+01, 1.128114232815630835D+01,-3.141814347069891511D+01, 2.185520874118305557D+01, &
                       4.204688908817206361D+00,-2.039396651211594147D+01, 1.703214350560853418D+01,-3.858035665251962953D+00, &
                      -5.664929771474184861D+00, 6.236013601339658763D+00,-8.349084962847823643D-01,-5.076648434675173682D+00, &
                       8.069437229646323928D+00,-7.693517618528387558D+00, 5.264540454351635645D+00,-2.378232295562936471D+00, &
                       9.488986688091295696D-02, 1.211665432819921451D+00,-1.641788465968915700D+00, 1.505783472315168181D+00, &
                      -1.120376089133631847D+00, 7.202627618735454318D-01,-4.293308092527763353D-01, 2.880519253847271255D-01, &
                      -2.765909110132210857D-01, 3.540771325316399709D-01,-4.702415607316920432D-01, 5.886079250132733032D-01, &
                      -6.804407177384068639D-01, 7.365579045711867501D-01,-7.516232294769148448D-01, 7.332518369461620278D-01, &
                      -6.840940481215002089D-01, 6.145375558739786248D-01,-5.263752170486778459D-01, 4.305607160730071659D-01, &
                      -3.307504282674695872D-01, 2.429511351824290011D-01,-1.749591466871897039D-01, 1.456366850874822316D-01, &
                      -1.595894233786172844D-01, 2.280395161025021433D-01,-3.413711397379035062D-01, 4.950015142033588056D-01, &
                      -6.617826176259098414D-01, 8.243458750063681340D-01,-9.468449720632085009D-01, 1.012979648696065826D+00, &
                      -9.946084380621768029D-01, 8.932243609342348512D-01,-7.019325270091921753D-01, 4.480078578272302381D-01, &
                      -1.456381084707468188D-01,-1.604285443826692914D-01, 4.506359421881006022D-01,-6.820827325803067165D-01, &
                       8.491728004462286705D-01,-9.252486609867097700D-01, 9.260533367474167443D-01,-8.402031897556576645D-01, &
                       6.982570225067336045D-01,-4.946654104012083719D-01, 2.669962789856030194D-01,-1.042528939201309117D-02, &
                      -2.319933928270803414D-01, 4.654858099944349514D-01,-6.436180173401860882D-01, 7.780252822887177011D-01, &
                      -8.272840865922600484D-01, 8.197748606729908794D-01,-7.271116816505861502D-01, 5.992292042110317629D-01, &
                      -4.191980233493379782D-01, 2.522483458861179417D-01,-8.121597742330896597D-02,-2.652177515646544220D-02, &
                       1.028351781199190462D-01,-8.958878867947789315D-02, 4.230973855067972356D-02, 8.503909842995018009D-02, &
                      -2.142619566112575480D-01, 3.865561094837716705D-01,-5.104727910435372662D-01, 6.353554205789004872D-01, &
                      -6.678424785380696616D-01, 6.769052092946485910D-01,-5.740789051340543514D-01, 4.533865101372042683D-01, &
                      -2.311110970314779189D-01, 2.495118155831516429D-02, 2.508326503987150513D-01,-4.617568849058834579D-01, &
                       7.064398492531979157D-01,-8.419497320020807862D-01, 9.894625411708843910D-01,-1.001995584015555441D+00, &
                       1.027606495406535592D+00,-9.137586391106879979D-01, 8.347545540345707726D-01,-6.273557157264747497D-01, &
                       4.895538883108097039D-01,-2.415916082478662963D-01, 1.027651673270286170D-01, 1.284590129111473078D-01, &
                      -2.133464520148750654D-01, 3.796869867196227544D-01,-3.706104772018420923D-01, 4.433984947033918766D-01, &
                      -3.235750823461782666D-01, 2.994380356487793549D-01,-7.888021235143724552D-02,-1.910472056257344828D-02, &
                       3.043509756878240990D-01,-4.298489130677841108D-01, 7.233655581626995401D-01,-8.161936207950235556D-01, &
                       1.054311731055406431D+00,-1.057859070622381603D+00, 1.189505276805579825D+00,-1.071819355511486993D+00, &
                       1.077289848919564141D+00,-8.458743245226547636D-01, 7.443908003078970603D-01,-4.461507052923350813D-01, &
                       2.868208498373531756D-01, 8.293917061440093941D-03,-1.688393814622326516D-01, 3.930003549460732160D-01, &
                      -5.148213607425790039D-01, 6.232961968619017412D-01,-6.962255871378753014D-01, 6.749680427277989780D-01, &
                      -7.157027353444623818D-01, 5.769543174180707945D-01,-6.151602424632486299D-01, 3.892025319337197864D-01, &
                      -4.506933666706669506D-01, 1.804687704243498336D-01,-2.727244303574414830D-01, 1.174346932607069939D-02, &
                      -1.137579568068424057D-01,-7.811421895993872488D-02, 1.779393134675293781D-02,-8.675295089134033022D-02, &
                       1.335093001151206882D-01,-5.514694464669991220D-02, 2.441624548660136784D-01,-5.518933795483266236D-02, &
                       3.328770195966135881D-01,-1.532145946078544430D-01, 3.579812932040986606D-01,-3.542418291888540516D-01, &
                       3.138078445399646865D-01,-5.530414472417091165D-01, 2.996207005996078809D-01,-5.735406549188309944D-01, &
                       4.429014477665748628D-01,-3.892397979507834505D-01, 6.319796721730875921D-01,-3.224362675105876264D-01, &
                       5.280557016061210307D-01,-5.925152169592526885D-01, 3.097933933610434454D-01,-6.161398121500974989D-01, &
                       5.824585099080735739D-01,-2.863603061506506120D-01, 5.579069742591284964D-01,-5.790177379949956737D-01, &
                       1.778203303660013113D-01,-2.248159329174458654D-01, 4.570279405020428731D-01,-2.082694991220853942D-01, &
                      -1.639463813834285966D-01, 1.219471051626438846D-01, 1.315907246465166380D-01,-1.561922951483753763D-01, &
                      -9.000173966500184253D-02, 3.555651953426591794D-01,-4.594671107387718334D-01, 4.051183609553802301D-01, &
                      -2.826594312508529105D-01, 1.664628231342878961D-01,-8.573386548433219179D-02, 3.940897144775467459D-02, &
                      -1.632396319789397934D-02, 6.095399479649578345D-03,-2.034197210251703809D-03, 5.959585007671103435D-04, &
                      -1.489396688112877424D-04, 3.040429693056684047D-05,-4.735729642331119566D-06, 4.982902654694162416D-07, &
                      -2.645962918790746080D-08/)
    
        wJ1 = (/-2.594301879688918743D-03, 2.339490069567079153D-02,-9.478827695764932559D-02, 2.269100107268227362D-01, &
                       -3.508900104631907935D-01, 3.515648778060092017D-01,-1.994483426338835574D-01, 9.293484158449249674D-03, &
                        8.048790331434139966D-02,-5.691508909482047296D-02, 2.075411236619714023D-02,-5.275228907131672418D-02, &
                        1.348060780914397960D-01,-1.939017682982676627D-01, 1.854125006516780250D-01,-1.238162114378456441D-01, &
                        5.435433909585837137D-02,-1.270327896980103649D-02, 7.560629247467166685D-03,-2.713648547163816441D-02, &
                        5.383176198347366936D-02,-7.463014212375274070D-02, 8.413254388281748986D-02,-8.280914671407160754D-02, &
                        7.393816085156304507D-02,-6.119942568117115594D-02, 4.746907894866293082D-02,-3.450474036399159977D-02, &
                        2.313104032185092293D-02,-1.354138672902785445D-02, 5.604638952874172256D-03, 9.486899997092870969D-04, &
                       -6.382314415042419746D-03, 1.091601370650388016D-02,-1.467479730415811347D-02, 1.771829574440562591D-02, &
                       -2.001990396213982823D-02, 2.152720536808908763D-02,-2.214634467574278301D-02, 2.181839780171816040D-02, &
                       -2.049006667874660528D-02, 1.819226618052163791D-02,-1.498058384418661862D-02, 1.101097898088004151D-02, &
                       -6.444381466006598655D-03, 1.527769500667217105D-03, 3.533045890696502513D-03,-8.470577184056853060D-03, &
                        1.310891443998619434D-02,-1.722186749311952272D-02, 2.071070403809925978D-02,-2.340849145682623311D-02, &
                        2.528895672008608583D-02,-2.620974028011485712D-02, 2.617463680378217042D-02,-2.501974545161460978D-02, &
                        2.276381800709064221D-02,-1.923696679546721411D-02, 1.454384146700743972D-02,-8.605818651281786635D-03, &
                        1.738844503682727928D-03, 5.947517515245405624D-03,-1.385885800081448037D-02, 2.171703574368676753D-02, &
                       -2.872814698074024550D-02, 3.461838232135941440D-02,-3.858941536214401807D-02, 4.060565074712448042D-02, &
                       -4.004001014668581715D-02, 3.725276070158690944D-02,-3.184536604257427045D-02, 2.462644816656397312D-02, &
                       -1.539955242663786951D-02, 5.437594699869644464D-03, 5.315602443775078317D-03,-1.514798704965643165D-02, &
                        2.414409295479183482D-02,-3.032245442310742278D-02, 3.416403682006976389D-02,-3.375211066545635158D-02, &
                        3.045425146826126819D-02,-2.272683240323554107D-02, 1.317397219169818418D-02,-6.119848802195663844D-04, &
                       -1.115685958739305421D-02, 2.344848055710842955D-02,-3.169397030448774938D-02, 3.817520279763731567D-02, &
                       -3.809651169502087376D-02, 3.552394037620545952D-02,-2.568176960339346032D-02, 1.487169880055576841D-02, &
                        1.997965754252607837D-03,-1.644546122863936588D-02, 3.480415781283187349D-02,-4.677537576184215284D-02, &
                        6.111820051573001178D-02,-6.590953690231306228D-02, 7.332112878498593667D-02,-6.941405787722150500D-02, &
                        7.037058951810715168D-02,-5.921683716330301134D-02, 5.655416878814193554D-02,-4.097804417398832194D-02, &
                        3.810136479256000241D-02,-2.058552251762190213D-02, 2.012009960638341116D-02,-1.904932224648239053D-03, &
                        5.324133880430357950D-03, 1.357085636140308721D-02,-5.592398994521355186D-03, 2.589415408390925710D-02, &
                       -1.292224858788840886D-02, 3.561835458678683924D-02,-1.700336052804993919D-02, 4.302272543873187499D-02, &
                       -1.764501119714675936D-02, 4.776824187802861110D-02,-1.404967795338326469D-02, 4.909105531804618811D-02, &
                       -5.268363386320836297D-03, 4.641444115511134810D-02, 9.124170345995139680D-03, 3.991380510079723526D-02, &
                        2.859680534219202416D-02, 3.073032173944831302D-02, 5.151734601415858261D-02, 2.097520864952298614D-02, &
                        7.511832720419686638D-02, 1.354148686126329000D-02, 9.576070647813632319D-02, 1.140272504589405489D-02, &
                        1.097360846564945647D-01, 1.641156037670641471D-02, 1.142607074214286172D-01, 2.822722328175052836D-02, &
                        1.079741026261866604D-01, 4.397838549427537935D-02, 9.072532008032478668D-02, 5.865331390992817306D-02, &
                        6.310580311949447185D-02, 6.592890500839264367D-02, 2.624768245556352575D-02, 5.954092122676180737D-02, &
                       -1.789442088462365327D-02, 3.555984804371448149D-02,-6.551973841104229146D-02,-4.489823587807304991D-03, &
                       -1.089659607197719371D-01,-5.052377390460756346D-02,-1.345089168583019634D-01,-8.286532978746213862D-02, &
                       -1.229744228392625760D-01,-7.728334735138943368D-02,-5.831402820224085293D-02,-2.088719549713474385D-02, &
                        5.205792146829352207D-02, 6.258826423151815643D-02, 1.537593434905694390D-01, 9.979165672991972824D-02, &
                        1.545936043232943868D-01, 1.442374184540378551D-02, 1.426608350698001064D-02,-1.462978689171700042D-01, &
                       -1.285617304470880462D-01,-1.519560295563492092D-01,-4.700754607069132507D-02, 1.106161196544740571D-01, &
                        1.168039951766200457D-01, 2.366820009670577984D-01,-1.283956636561735531D-01, 1.406024618127853579D-02, &
                       -4.011626826391421763D-01, 2.656020926376266300D-01,-1.291178852844470371D-01, 5.190017993719785450D-01, &
                       -5.307617026829943851D-01, 2.388177316428383157D-01,-4.213777919843004205D-01, 7.425023205388556757D-01, &
                       -6.039491377476802203D-01, 2.966551251977549986D-01,-3.157763414193298090D-01, 5.842678377857475347D-01, &
                       -7.375386163194860289D-01, 6.321787190293964853D-01,-3.871768116388242809D-01, 1.635153404860406612D-01, &
                       -3.134030444480916111D-02,-2.007964312111037639D-02, 2.747060309825222202D-02,-1.969027573928272892D-02, &
                        1.080160200622757964D-02,-4.917741573140644099D-03, 1.900872502605262596D-03,-6.228037791451519409D-04, &
                        1.698280169898911708D-04,-3.718088395106743311D-05, 6.139140918052498045D-06,-6.796809441105533029D-07, &
                        3.780807126974975489D-08/)
    
    end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Ep_dmxz( m, s0, w, Fi, r0, r1, r2, Ep )
! !
! ! Projection of the electric field from a magnetic dipole in the (x,z) plane onto an 
! ! arbitrary edge defined by its node coordinates.
! !
! ! m  	 -> Input: Dipole moment.
! ! s0     -> Primary medium conductivity.
! ! w      -> 2*pi*freq
! ! Fi 	 -> Input: Fi(1)=cos(Dip); Fi(2)=sin(Dip). Dip: Angle of the moment vector relative to the vertical (z).
! ! r0 	 -> Input: vector with the coordinates of the source.
! ! r1, r2 -> Input: vector with the oordinates of the nodes that define the edge.
! ! Ep 	 -> Output: Projection of the electric field at the center of the edge.
! ! Global constants:
! !   pi, mu0, eps0
! !
 implicit none
 real(dbl), intent(in):: m, s0, w, Fi(2), r0(3), r1(3), r2(3)
 complex(dbl), intent(out):: Ep
 real(dbl):: x1, x2, y1, y2, z1, z2, l, R
 complex(dbl):: ikr, e0

 x1 = r1(1) - r0(1)
 x2 = r2(1) - r0(1)

 y1 = r1(2) - r0(2)
 y2 = r2(2) - r0(2)

 z1 = r1(3) - r0(3)
 z2 = r2(3) - r0(3)

 l = dsqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )            ! Length of the edge,.
 R = dsqrt( (x1+x2)**2 + (y1+y2)**2 + (z1+z2)**2 ) / 2.d0     ! Radial distance from the source to the center of the edge.

 ikr = cdsqrt( cmplx( -w**2 * mu0 * eps0, w * mu0 * s0, dbl ) ) * R

 E0 = -m * (cmplx(0.d0, w*mu0, dbl)) * (1.d0 + ikr) * cdexp(-ikr) / (4.d0*pi*l*R**3)

 Ep = E0 * ( (x1*y2 - x2*y1)*Fi(1) + (y1*z2 - y2*z1)*Fi(2) )

 end subroutine Ep_dmxz
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Ep_dmxyy( m, s0, w, Fi, r0, r1, r2, Ep )
!
! Projection of the electric field from a magnetic dipole in the (x,y) plane onto an 
! arbitrary edge defined by its node coordinates.
!
! m  	 -> Input: Dipole moment.
! s0     -> Primary medium conductivity.
! w      -> 2*pi*freq
! Fi 	 -> Input: Fi(1)=cos(Dip); Fi(2)=sin(Dip). Dip: Angle of the moment vector relative to the x axis.
! r0 	 -> Input: vector with the coordinates of the source.
! r1, r2 -> Input: vector with the oordinates of the nodes that define the edge.
! Ep 	 -> Output: Projection of the electric field at the center of the edge.
! Global constants:
!   pi, mu0, eps0
!
implicit none
real(dbl), intent(in):: m, s0, w, Fi(2), r0(3), r1(3), r2(3)
complex(dbl), intent(out):: Ep
real(dbl):: x1, x2, y1, y2, z1, z2, l, R
complex(dbl):: ikr, e0

x1 = r1(1) - r0(1)
x2 = r2(1) - r0(1)

y1 = r1(2) - r0(2)
y2 = r2(2) - r0(2)

z1 = r1(3) - r0(3)
z2 = r2(3) - r0(3)

l = dsqrt( (x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2 )            ! Length of the edge.
R = dsqrt( (x1+x2)**2 + (y1+y2)**2 + (z1+z2)**2 ) / 2.d0     ! Radial distance from the source to the center of the edge.

ikr = cdsqrt( cmplx( -w**2 * mu0 * eps0, w * mu0 * s0, dbl ) ) * R

E0 = m * (cmplx(0.d0, w*mu0, dbl)) * (1.d0 + ikr) * cdexp(-ikr) / (4.d0*pi*l*R**3)

Ep = E0 * ( (z1*y2 - z2*y1)*Fi(1) + (x1*z2 - x2*z1)*Fi(2) )

end subroutine Ep_dmxyy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Ep_dmv(moment, sig0, w, Fi, r0, r1, r2, Ep)
implicit none
real(dbl), intent(in)		:: moment, sig0(0:n_layer), w, Fi(2)
real(dbl), intent(in)		:: r0(3), r1(3), r2(3)
complex(dbl), intent(out)	:: Ep
complex(dbl)			:: s_Ex, s_Ey, s_Ez
integer			:: navalEx, navalEy, navalEz
integer			:: i, j, c
real(dbl)                  	:: norma
real(dbl)			:: x,y,z
real(dbl)                         :: ab(N_kk), pej0(N_kk), pej1(N_kk)
complex(dbl)	      		   :: Ex, Ey  
complex(dbl)            		   :: Kernel_aux1, Kernel_aux2, kernel_aux3, aux4
complex(dbl)            		   :: kernel_Ez, aux_rte
real(dbl)               		   :: pej0kr
complex(dbl)            		   :: expuzp, expuzn, uz, Asu, rtmexp, uesp
complex(dbl)                   	   :: tanhe
complex(dbl),dimension(0:n_layer)    :: u, Y_int, Z_int, F, A, k_no, E
complex(dbl),dimension(n_layer)      :: Y_ap, Z_ap
complex(dbl),dimension(0:n_layer-1)  :: Rte, Rtm
real(dbl)                      	    :: senrad, r, kr, dx, dy, dz
real(dbl)                      	  ::  dxr, dyr, dx2,dy2, dxy, dxysin, dxycos 
real(dbl)			    :: cosrad			
integer		    	 :: cam      
call J0J1Wer( ab, pej0, pej1 )
! midpoint:
   x = (r1(1) + r2(1))*0.5d0
   y = (r1(2) + r2(2))*0.5d0 
   z_rc = (r1(3) + r2(3))*0.5d0
 
   !Distance between source and nodes 1 e 2:
    x = x - r0(1)
    y = y - r0(2)    
   h0=r0(3)	
   ! ray:   
   r=sqrt(x**2+y**2)
   if(r<10.d-2) then
      r = 10.d-2
   end if	
   dxr  = x/r
   dyr  = y/r
   
   sigma= 1.d-12 + ci*w*8.85*1.d-12

   zeta=ci*w*mu0 
 s_Ex = (0.d0, 0.d0)
 s_Ey = (0.d0, 0.d0)
 s_Ez = (0.d0, 0.d0)
Do j=1,N_kk
	kr = (ab(j))/r
	zeta=ci*w*mu0 
	k_no(0) = sqrt(-zeta*sigma)
	k_no(1:n_layer) = sqrt(-zeta*sig0(1:n_layer))

  	Z_int(0)=u(0)/sigma
  	u(0:n_layer)=sqrt(kr**2-k_no(0:n_layer)**2) 	 
  	Y_int(0:n_layer)=u(0:n_layer)/zeta
    	Y_ap(n_layer)=Y_int(n_layer) 	

	DO i=n_layer-1,1,-1 
 		tanhe = tanh(u(i)*esp(i)) 
		Y_ap(i)=Y_int(i)*((Y_ap(i+1)+Y_int(i)*tanhe)/(Y_int(i)+Y_ap(i+1)*tanhe))		
	END DO
 	
	! Reflection coefficient:
	DO i=0,n_layer-1
		Rte(i)=(Y_int(i)-Y_ap(i+1))/(Y_int(i)+Y_ap(i+1))
	END DO

	E(0)=(exp(u(0)*h0)*(zeta))/(4*pi*u(0))
	 E(1)=E(0)*((1+Rte(0))/(1+Rte(1)*exp(-2*u(1)*esp(1))))
 	 DO i=2,n_layer-1	
		E(i)=E(i-1)*exp(u(i)*esp(i))*((1+Rte(i-1))*exp(-u(i-1)*esp(i-1))/(1+Rte(i)*exp(-2*u(i)*esp(i))))
	 END DO
  	E(n_layer)=E(n_layer-1)*(1+Rte(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))
   IF(z_rc<0.) THEN 
	cam=0			
   Else
	DO c=1,n_layer-1			
		IF(z_rc<pc(c))THEN
			cam=c
			exit
		END IF
	END DO
  		IF(z_rc>=pc(n_layer-1))  cam=n_layer
   End if
       IF(cam<n_layer) THEN 
		expuzp = exp(u(cam)*(z_rc-(pc(cam)+ esp(cam))))
		expuzn = exp(-u(cam)*(z_rc-(pc(cam)- esp(cam))))			
		Kernel_aux1 =  E(cam)*(expuzn + Rte(cam)*expuzp)		
	ELSE 
		uz = u(n_layer)*(z_rc-pc(n_layer-1)) 
		expuzn = exp(-uz)
		kernel_aux1 = E(n_layer)*expuzn						   		
	END IF   
                    	
    Ey = (kernel_aux1*(-dxr)*(kr)**2) 
    Ex =(kernel_aux1*(dyr)*(kr)**2)
	       
    S_Ey = S_Ey + Ey*pej1(j) 
    S_Ex = S_Ex + Ex*pej1(j)         
 end do
 S_Ey = S_Ey/r
 S_Ex = S_Ex/r

   dx= r2(1)-r1(1)
   dy= r2(2)-r1(2)
   dz= r2(3)-r1(3)
   norma = sqrt( dx**2 + dy**2 + dz**2)
! Projection of the electric field onto the edge: 
  Ep = (s_Ex*dx + s_Ey*dy)/norma
End Subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine campo_primariov( N_edge, node, x, y, z, rt, m, sig0, w, Fi, ar, E_prim )
implicit none
integer, intent(in):: N_edge, node(4), ar(6)
real(dbl), intent(in):: x(4), y(4), z(4), rt(3), m, sig0(0:N_layer), w, Fi(2)
real(dbl):: r1(3), r2(3)
complex(dbl), intent(inout):: E_prim(N_edge)
integer:: i, j, k

k = 0
do i = 1, 3
	do j = i+1, 4
		if (node(i)<node(j)) then
			r1(1) = x(i)
			r1(2) = y(i)
			r1(3) = z(i)
			r2(1) = x(j)
			r2(2) = y(j)
			r2(3) = z(j)
		else
			r1(1) = x(j)
			r1(2) = y(j)
			r1(3) = z(j)
			r2(1) = x(i)
			r2(2) = y(i)
			r2(3) = z(i)
		end if

		k = k + 1
        if (E_prim( ar(k) ) == (0.d0,0.d0) ) then
	call Ep_dmv(moment, sig0, w, Fi, rt, r1, r2, E_prim(ar(k)))	
        end if
	end do
end do

end subroutine campo_primariov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine campo_direto_estv(moment,sig0, w, r, Hp )
implicit none 
real(dbl), intent(in)	 :: moment, w, r
real(dbl), intent(in)      :: sig0(0:N_layer)
integer   		 :: navalHx, navalHy, navalHz
complex(dbl)		 ::s_Hx, s_Hy, s_Hz
complex(dbl), intent(out)  :: Hp(2)
integer		  :: i, j, c
real(dbl)                   :: norma, kr
real(dbl)                   :: y, z, x
real(dbl)                         :: ab(N_kk), pej0(N_kk), pej1(N_kk)
complex(dbl)	      		   :: Hz  
complex(dbl)            		   :: Kernel_aux1, Kernel_aux2, kernel_aux3, aux4
complex(dbl)            		   :: kernel_Ez, aux_rte
real(dbl)               		   :: pej0kr
complex(dbl)            		   :: expuzp, expuzn, uz, Asu, rtmexp, uesp
complex(dbl)                   	   :: tanhe
complex(dbl),dimension(0:n_layer)    :: u, Y_int, Z_int, k_no, E
complex(dbl),dimension(n_layer)      :: Y_ap, Z_ap
complex(dbl),dimension(0:n_layer-1)  :: Rte, Rtm 

     
call J0J1Wer( ab, pej0, pej1 )
   
   sigma= 1.d-12 + ci*w*8.85*1.d-12
   zeta=ci*w*mi0 

 s_Hz = (0.d0, 0.d0)
Do j=1,N_kk
	kr = (ab(j))/r
	zeta=ci*w*mu0 
	k_no(0) = sqrt(-zeta*sigma)
	k_no(1:n_layer) = sqrt(-zeta*sig0(1:n_layer))

  	Z_int(0)=u(0)/sigma
  	u(0:n_layer)=sqrt(kr**2-k_no(0:n_layer)**2) 	 
  	Y_int(0:n_layer)=u(0:n_layer)/zeta
    	Y_ap(n_layer)=Y_int(n_layer) 	

	DO i=n_layer-1,1,-1 
 		tanhe = tanh(u(i)*esp(i)) 
		Y_ap(i)=Y_int(i)*((Y_ap(i+1)+Y_int(i)*tanhe)/(Y_int(i)+Y_ap(i+1)*tanhe))		
	END DO
 	
	!Reflection coefficient:
	DO i=0,n_layer-1
		Rte(i)=(Y_int(i)-Y_ap(i+1))/(Y_int(i)+Y_ap(i+1))
	END DO

	E(0)=(exp(u(0)*h0)*(zeta))/(4*pi*u(0))
	 E(1)=E(0)*((1+Rte(0))/(1+Rte(1)*exp(-2*u(1)*esp(1))))
 	 DO i=2,n_layer-1	
		E(i)=E(i-1)*exp(u(i)*esp(i))*((1+Rte(i-1))*exp(-u(i-1)*esp(i-1))/(1+Rte(i)*exp(-2*u(i)*esp(i))))
	 END DO
  	E(n_layer)=E(n_layer-1)*(1+Rte(n_layer-1))*exp(-u(n_layer-1)*esp(n_layer-1))

	expuzp = exp(u(0)*(-0.5d0))
	expuzn = exp(-u(0)*(-0.5d0))  
	Hz = (E(0)*(expuzn + Rte(0)*expuzp)*(1/zeta)*(kr)**3)
	       
    S_Hz = S_Hz + Hz*pej0(j)          
 end do
 
   S_Hz = S_Hz/r
   Hp(1) = S_Hx  
   Hp(2) = S_Hz	
End Subroutine

subroutine para_range(n1, n2, nprocs, irank,ista, iend)
implicit none
integer,intent(in) :: n1, n2, nprocs, irank
integer,intent(out) :: ista, iend
integer :: iwork1, iwork2
iwork1 = (n2 - n1 + 1) / nprocs
iwork2 = MOD(n2 - n1 + 1, nprocs)
ista = irank * iwork1 + n1 + MIN(irank, iwork2)
iend = ista + iwork1 - 1
if (iwork2 > irank) iend = iend + 1
end subroutine


end module VFE3D_Hfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
INCLUDE 'mkl_pardiso.f90'
program perfilagem3D
use VFE3D_Hfield
implicit none
include 'mpif.h'
integer, parameter:: Nn_estimate = 100
integer::  Nnn, N_el_sec
integer:: i, j, k, Conta, ind_loop, maior, elem, porc
integer, allocatable::  ia(:), ja(:)
integer, allocatable:: ja_temp(:,:), ind_ele_obs(:)
integer:: node(4), edge(6), ar(6)
integer:: sinal
integer:: c
real(dbl):: x(4), y(4), z(4), z_med
real(dbl)::  yr, xc, yc, zc, r_c
real(dbl):: ang_dipoleh, ang_dipolev
real(dbl):: Na(3), Nb(3), Nc(3), Nd(3), P(3)
complex(dbl), allocatable:: global_temp(:,:), G_ij(:,:,:)
complex(dbl), allocatable:: admit(:,:), M_global(:),Hs(:,:), Hss(:), Hsh(:,:), Hsv(:,:)
complex(dbl):: Mel(6,6), eta(2), D_eta(2), G_el(6,6), Eph(6), Epv(6)
complex(dbl):: zet0, eta0
complex(dbl):: obs(3), Hph(3), Hpv(3)
complex(dbl), allocatable:: E_primh(:), E_primv(:), Hcx(:), Hcp(:), H_pri(:), H_sec(:)
real(dbl)		:: cont
real(dbl) :: start, finish
!
real:: t1, t2, t3
!
! PARDISO VARIABLES:
integer, parameter:: MTYPE = 6
integer, parameter:: MAXFCT = 1,  MNUM = 1,   MSGLVL = 1
integer:: PHASE, ERROR
integer:: PT(64)
!TYPE(MKL_PARDISO_HANDLE):: PT(64)
integer:: IPARM(64)
integer, allocatable::  PERM(:)
!%%%%%  MPI variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer myrank, nprocs, ierr, irank, ireq, errorcode
integer, allocatable, dimension(:):: jjsta, jjlen, iireq
integer ista, iend, jsta, jend, vec_base(6)
integer istatus (MPI_STATUS_SIZE)
integer, allocatable, dimension(:):: vector_int
real(dbl), allocatable, dimension(:):: vector_real
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call cpu_time(t1)
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

call leituras()

allocate( admit(0:N_layer+N_block,2), vector_real(6+N_layer) )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%PARALLELIZING PIECE IN THE PROCESS 0
if (myrank == 0) then !%%%%%%%%%%%%%%%%%
w = 2.d0*pi*freq
zet0 = cmplx(0._dbl, w*mu0, dbl)
eta0 = cmplx(sig0(0), w*eps0, dbl)
admit(0,:) = eta0
do i = 1, N_layer+N_block
	admit(i,:) = cmplx(sig(i,1:2), w*eps0, dbl)
end do

!N_el_sec = count(region/=0)
N_el_sec = count(region>N_layer)
allocate( ia(N_edge+1), ja_temp(Nn_estimate,N_edge) )
allocate( global_temp(Nn_estimate,N_edge), G_ij(6,6,N_el_sec))
ia = 0
ja_temp = 0
global_temp = zero
k = 0
do i = 1, N_ele
	node = Nos(:,i)
	x = coord( 1, node(1:4) ) !the matrix nos(4, n_ele) is read and each element transfer your node to
	y = coord( 2, node(1:4) ) !vector node after the vector node(4) transfer to vectors x(4), .. 
	z = coord( 3, node(1:4) ) !y(4), z(4) the coordinates of each node that belong to element i
	z_med=sum(z)/4            !gruping for your coordinates xyz each is a vector.

 	
   	 IF(z_med<0.) THEN 
		j=0			
   	Else
		DO c=1,N_layer-1 
			IF(z_med<pc(c))THEN  !test of the interface of the layer
				j=c
				exit
			END IF
		END DO
  		IF(z_med>=pc(N_layer-1))  j=N_layer 
   		End if

	eta = admit( region(i), : ) ! here conductivity
	D_eta = admit( region(i), : ) - admit( j, : ) ! difference between the camada and block

	call Matriz_el( node, x, y, z, region(i), zet0, eta, D_eta, Mel, G_el )

	! Local matrices for calculating the source vector.
	if ( region(i)> N_layer ) then
		k = k + 1
		G_ij(:,:,k) = G_el   ! here is identify the elementa matrix  with element belong to block.
	end if
	
	ar = arestas(:,i) ! this vector contain the edges of the element_i of the loop current
	
	call Monta_colunas( N_edge, Nn_estimate, Mel, ar, ia, ja_temp, Global_temp ) !mouting of tempory global matriz
end do
!!!***************** read all elements, finish

maior = -1
do ind_loop = 1, N_edge
	conta = count(ja_temp(:,ind_loop)/=0) !count for column how many element not null
	if (conta>maior) then
		maior = conta
	end if
end do

print *, 'Links:', maior

!!!*****************

call Remove_zeros( N_edge, N_ele, Nn_estimate, ia, ja_temp, global_temp )

NNN = sum(ia)

allocate( ja(Nnn), M_global(Nnn) )

call Montagem_final( N_edge, N_ele, Nn_estimate, Nnn, ja_temp, global_temp, ia, ja, M_global )

deallocate( ja_temp, global_temp )

call Cond_de_front( Nnn, N_edge, N_front, ia, fronteira, M_global )

open(666, file='hy.dat', status='replace')

open(222, file='hz.dat', status='replace')


!!!!!!LOCATE THE ELEMENTS OF THE OBSERVATION POINTS

allocate( ind_ele_obs(N_obs) )

do i = 1, N_obs
	xr = x_obs(i) 
	yr = y_obs(i)
	zr = z_obs

	do j = 1, N_ele
		if ( region(j) == 0 ) then
			node = Nos(:,j)
			x = coord( 1, node(1:4) )
			y = coord( 2, node(1:4) )
			z = coord( 3, node(1:4) )

			!SELECT THE OBSERVATION ELEMENT
			if ( abs((z(1)+z(2)+z(3)+z(4))/4.d0 - zr) < 0.5d0 ) then
				Na = (/ x(1), y(1), z(1) /)
				Nb = (/ x(2), y(2), z(2) /)
				Nc = (/ x(3), y(3), z(3) /)
				Nd = (/ x(4), y(4), z(4) /)
				P = (/ xr, yr, zr /)
				if ( localiza_obs( Na, Nb, Nc, Nd, P )==1 .and. localiza_obs( Na, Nc, Nd, Nb, P )==1 .and. &
						localiza_obs( Na, Nb, Nd, Nc, P )==1 .and. localiza_obs( Nb, Nc, Nd, Na, P )==1 ) then
						ind_ele_obs(i) = j
						exit
				end if
			end if

		end if
	end do
end do

allocate( Hsh(N_edge,N_obs), Hsv(N_edge,N_obs), Hss(N_edge*N_obs) )
allocate( E_primh(N_edge), E_primv(N_edge) )
allocate(Hcx(N_obs), Hcp(N_obs), H_pri(N_obs), H_sec(N_obs))

call pardisoinit( pt, mtype, iparm )
iparm(6) = 1    ! The solver stores the solution on the right-hand side b.

allocate( PERM(N_edge) )

!phase = 11    ! 11: Analysis
!phase = 22    ! 22: Numerical factorization
phase = 12     ! 12: Analysis, numerical factorization
call pardiso( pt, maxfct, mnum, mtype, phase, N_edge, M_global, ia, ja, perm, N_obs, iparm, msglvl, Hsh, Hss, error )

!!!!!! RESOLVE VMD
ang_dipolev=0.d0
ang_dipoleh = pi/2.d0
Fih(1) = cos(ang_dipoleh)
Fih(2) = sin(ang_dipoleh)

Fiv(1) = cos(ang_dipolev)
Fiv(2) = sin(ang_dipolev)

Hsh = zero 
Hss = zero 
Hsv = zero 


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!variables integer and real, they'll be sending to other process

vector_real(1)=w
vector_real(2:3)=Fih
vector_real(4:5)=Fiv
vector_real(6:6+n_layer)=sig0

endif !!end of the process 0
!000000000000000000000000000000000000000000000000000000000000
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

call MPI_BARRIER(MPI_COMM_WORLD, ierr)

if(myrank /= 0)then
	N_el_sec = count(region>N_layer)
    allocate(G_ij(6,6,N_el_sec))
    allocate( Hsh(N_edge,N_obs))
    allocate( Hsv(N_edge,N_obs))
    allocate( E_primh(N_edge))
    allocate( E_primv(N_edge))
	Hsh=zero
    Hsv=zero
endif

CALL MPI_BCAST(vector_real, 6+N_layer, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

w= vector_real(1)
Fih= vector_real(2:3)
Fiv= vector_real(4:5)
sig0= vector_real(6:6+n_layer)

CALL MPI_BCAST(G_ij, 6*6*N_el_sec, MPI_COMPLEX16, 0, MPI_COMM_WORLD, ierr)

ALLOCATE (jjsta(0:nprocs-1))
ALLOCATE (jjlen(0:nprocs-1))
ALLOCATE (iireq(0:nprocs-1))

ireq=0
iireq=0

DO irank = 0, nprocs - 1
	CALL para_range(1, n_obs, nprocs, irank, jsta, jend)
	jjsta(irank) = jsta
	jjlen(irank) = n_edge * (jend - jsta + 1)
ENDDO

call para_range(1, n_obs, nprocs, myrank, ista, iend)

do i = ista, iend
	rt(1) = x_obs(i) - offset
	xr = x_obs(i)
    yr = y_obs(i)
 	zr = z_obs        
	k = 0
   	E_primh = (0.d0,0.d0)
    E_primv = (0.d0,0.d0)
	do j = 1, N_ele

		if ( region(j) > N_layer ) then
		
			node = Nos(:,j)
			x = coord( 1, node(1:4) )
			y = coord( 2, node(1:4) )
			z = coord( 3, node(1:4) )
		
			call campo_primarioh( N_edge, node, x, y, z, rt, moment, sig0, w, Fih, arestas(:,j), E_primh )
            call campo_primariov( N_edge, node, x, y, z, rt, moment, sig0, w, Fiv, arestas(:,j), E_primv )
			k = k + 1
            Eph = E_primh(arestas(:,j))
			Eph = matmul(G_ij(:,:,k), Eph)
			Hsh( arestas(:,j), i ) = Hsh( arestas(:,j), i ) + Eph

            Epv = E_primv(arestas(:,j))
			Epv = matmul(G_ij(:,:,k), Epv)
			Hsv( arestas(:,j), i ) = Hsv( arestas(:,j), i ) + Epv
		end if
	end do
end do
! Data together 
IF (myrank == 0) THEN
	DO irank = 1, nprocs - 1
		CALL MPI_IRECV(Hsh(1,jjsta(irank)), jjlen(irank),MPI_COMPLEX16,&
							irank, 1, MPI_COMM_WORLD, iireq(irank),ierr)
		CALL MPI_IRECV(Hsv(1,jjsta(irank)), jjlen(irank),MPI_COMPLEX16,&
							irank, 1, MPI_COMM_WORLD, iireq(irank),ierr)
	ENDDO
	DO irank = 1, nprocs - 1
		CALL MPI_WAIT(iireq(irank), istatus, ierr)
	ENDDO

ELSEif (myrank /= 0)then

	CALL MPI_ISEND(Hsh(1,ista), jjlen(myrank), MPI_COMPLEX16,&
					0, 1, MPI_COMM_WORLD, ireq, ierr)
	CALL MPI_ISEND(Hsv(1,ista), jjlen(myrank), MPI_COMPLEX16,&
					0, 1, MPI_COMM_WORLD, ireq, ierr)
	CALL MPI_WAIT(ireq, istatus, ierr)
ENDIF

DEALLOCATE (jjsta, jjlen, iireq)

if (myrank == 0)then

PHASE = 33  ! 33: Solve, iterative refinement

call pardiso( pt, maxfct, mnum, mtype, phase, N_edge, M_global, ia, ja, perm, N_obs, iparm, msglvl, Hsh, Hss, error )

call pardiso( pt, maxfct, mnum, mtype, phase, N_edge, M_global, ia, ja, perm, N_obs, iparm, msglvl, Hsv, Hss, error )


!print *, 'Resolveu DM.   Erro:', error
!write(*,*) 'Offset:', offset

call Campo_direto_esth( moment, sig0, w, offset, Hph )
call Campo_direto_estv(moment, sig0, w, offset, Hpv )

do i = 1, N_obs
	elem = ind_ele_obs(i)
	node = Nos(:,elem)
	x = coord( 1, node(1:4) )
	y = coord( 2, node(1:4) )
	z = coord( 3, node(1:4) )
	xr = x_obs(i)
	yr = y_obs(i)
    zr = z_obs

    call calcula_obs( node, x,y,z, Hsh(arestas(:,elem),i), xr, yr, zr, obs )
    Hcx(i) = obs(2) + Hph(2)
    H_sec(i) = obs(2)
    H_pri(i) = Hph(2)
    write(666,'(3(ES20.8,x))')x_obs(i)-0.5d0*offset, real(Hcx(i)),aimag(Hcx(i))

	call calcula_obs( node, x,y,z, Hsv(arestas(:,elem),i), xr, yr, zr, obs )
    Hcx(i) = obs(3) + Hpv(2)
    H_sec(i) = obs(2)
    H_pri(i) = Hpv(2)
    write(222,'(3(ES20.8,x))')x_obs(i)-0.5d0*offset, real(Hcx(i)),aimag(Hcx(i))

end do

print *, 'Fim'
endif
CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
CALL MPI_FINALIZE(ierr)

end program

