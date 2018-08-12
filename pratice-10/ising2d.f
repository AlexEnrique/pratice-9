!
!      Simulacao do modelo de Ising bidimensional usando metodo de Metropolis
!
       program ising
       implicit real*8 (a-h,o-z)
       implicit integer (i-n)
       dimension s(1000,1000)
       dimension sM(10000000)
       dimension ener(10000000)
       dimension c(500)

!      entrada (arquivo input)
       open (1, FILE='entrada', STATUS='UNKNOWN')
       read (1,*) Nt        !passos MC termalizacao
       read (1,*) Ne        !passos MC equilibrio
       read (1,*) N         !tamanho da rede
       read (1,*) istart    !condicao inicial
       read (1,*) beta      !beta (=1/kT)
       read (1,*) cJ        !J (termo de troca)
       read (1,*) H         ! campo magnetico externo

!      arquivos de saida
       open (8, FILE='saida', STATUS='UNKNOWN')
       open (2, FILE='resultadosE_M', STATUS='UNKNOWN') ! resultados Energia/spin e Magnetizacao/spin
       open (3, FILE='correlacaoM', STATUS='UNKNOWN')  ! correlacao magnetizacao
       open (4, FILE='correlacaoE', STATUS='UNKNOWN')  ! correlacao energia


       ! inicializacao
       if (istart.eq.1) then
	 do j = 1, N
         do i = 1, N
	   s(i,j) = 1.d0
	 enddo
	 enddo
       else
	 do j = 1, N
         do i = 1, N
           r = rand()
	   if (r.le.0.5d0) then
	     s(i,j) = 1.d0
	   else
	     s(i,j) = -1.d0
	   endif
	 enddo
	 enddo
       endif


! Termalizacao
       do k = 1, Nt
         do j = 1, N
         do i = 1, N
           ! calcula vizinhos impondo cond. period.
           i1 = N - mod(N-i+1,N)
           i2 = 1 + mod(i,N)
           j1 = N - mod(N-j+1,N)
           j2 = 1 + mod(j,N)
           sj = s( i1, j ) + s( i2, j ) +
     .          s( i, j1 ) + s( i, j2 )

           ! calcula probabilidade de flipar s(i)
           p = exp( -2.d0 * beta * ( cJ * s(i,j) * sj + H * s(i,j) ) )
           if (p.ge.1.d0) then
             s(i,j) = -s(i,j)
           else
             r = rand()
             if (r.le.p) then
               s(i,j) = -s(i,j)
             endif
           endif
         end do
         end do
       end do
! fim termalizacao


       ! loop principal
       aM = 0.d0
       aM1 = 0.d0
       aM2 = 0.d0
       aE = 0.d0
       na = 0
       do k = 1, Ne
         sM(k) = 0.d0
         ener(k) = 0.d0
	 do j = 1, N
         do i = 1, N
	   ! calcula vizinhos impondo condicao periodica de contorno
	   i1 = N - mod(N-i+1,N)
	   i2 = 1 + mod(i,N)
	   j1 = N - mod(N-j+1,N)
	   j2 = 1 + mod(j,N)
	   sj = s( i1, j ) + s( i2, j ) +
     .          s( i, j1 ) + s( i, j2 )

	   ! calcula probabilidade de flipar s(i)
	   p = exp( -2.d0 * beta * ( cJ * s(i,j) * sj + H * s(i,j) ) )
	   if (p.ge.1.d0) then
	     s(i,j) = -s(i,j)
	     na = na + 1
	   else
	     r = rand()
	     if (r.le.p) then
	       s(i,j) = -s(i,j)
	       na = na + 1
	     endif
	   endif
           ! adiciona contribuicoes `a magnetizacao e energia
	   sM(k) = sM(k) + s(i,j)
           if (i.ge.2.and.j.ge.2) then
             ener(k) = ener(k) - cJ * s(i,j) * ( s(i-1,j) + s(i,j-1) )
     .	                                                - H * s(i,j)
           endif
           if (i.eq.N) then
             ener(k) = ener(k) - cJ * s(1,j) * s(N,j) - H * s(1,j)
           endif
           if (j.eq.N) then
             ener(k) = ener(k) - cJ * s(i,1) * s(i,N)
           endif
         enddo
	 enddo
	 sM(k) = sM(k) / dfloat(N**2)
	 ener(k) = ener(k) / dfloat(N**2)
	 aM = aM + sM(k)
	 aM1 = aM1 + sqrt(sM(k)*sM(k))
	 aM2 = aM2 + sM(k)*sM(k)
	 aE = aE + ener(k)
	 write (2,*) ener(k), sM(k)  !, sqrt( sM(k)*sM(k) )
       enddo

       accept = dfloat(na)/dfloat(Ne*N**2)
       aM = aM/dfloat(Ne)
       aM1 = aM1/dfloat(Ne)
       aM2 = aM2/dfloat(Ne)
       aE = aE/dfloat(Ne)

! calcula sigma
       sigmaM = 0.d0
       sigmaM1 = 0.d0
       sigmaM2 = 0.d0
       sigmaE = 0.d0
       do i = 1, Ne
         sigmaM = sigmaM + (sM(i) - aM)**2
         sigmaM1 = sigmaM1 + (sqrt(sM(i)**2) - aM1)**2
         sigmaM2 = sigmaM2 + (sM(i)**2 - aM2)**2
         sigmaE = sigmaE + (ener(i) - aE)**2
       enddo

       sigmaM = sqrt(sigmaM / dfloat(Ne))
       sigmaM1 = sqrt(sigmaM1 / dfloat(Ne))
       sigmaM2 = sqrt(sigmaM2 / dfloat(Ne))
       sigmaE = sqrt(sigmaE / dfloat(Ne))

! solucao exata a campo zero
       exactM0 = ( 1.d0 - (sinh(2.d0*beta*cJ))**(-4.d0) )**(0.125d0)

! saida
       write (8,*) " "
       write (8,*) "ising2d METROPOLIS"
       write (8,*) " ",Ne," iteracoes   N  =",N
       write (8,*) "  J  =",cJ,"   beta  =",beta,"   H  =",H
       write (8,*) "  aceitacao  =",accept*100,"%"
       write (8,*) " "

       write (8,*) "MAGNETIZACAO: ", aM, sigmaM / sqrt(dfloat(Ne))
       write (8,*) "Mabs: ", aM1, sigmaM1 / sqrt(dfloat(Ne))
       write (8,*) "Msqrt: ", sqrt(aM2), sigmaM2 / sqrt(dfloat(Ne))
       if (H.eq.0.d0) then
         write (8,*) "exato (H=0)  =", exactM0, 
     .                      "    diferenca  =", exactM0-aM1
       endif


       write (8,*) "correlacoes:"

       do k = 0, 20
         c(k) = 0.d0
         do i = 1, Ne-k
	   c(k) = c(k) + ( sM(i) - aM ) * ( sM(i+k) - aM )
	 enddo
	 c(k) = c(k)/dfloat(Ne-k)
	 c(k) = c(k)/sigmaM**2
	 write (8,*) k, c(k)
	 write (3,*) k, abs(c(k))
       enddo

       write (8,*) " "
       write (8,*) "ENERGIA: ", aE, sigmaE / sqrt(dfloat(Ne))

       write (8,*) "correlacoes:"

       do k = 0, 20
	 c(k) = 0.d0
	 do i = 1, Ne-k
	   c(k) = c(k) + ( ener(i) - aE ) * ( ener(i+k) - aE )
	 enddo
	 c(k) = c(k)/dfloat(Ne-k)
	 c(k) = c(k)/sigmaE**2
	 write (8,*) k, c(k)
	 write (4,*) k, abs(c(k))
       enddo

       
       end program
