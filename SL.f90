      program SL
      USE DISLIN
      implicit none
      ! parametros
      real*8, parameter :: DFREC = 0.05 ! [Hertz]
      real*8, parameter :: TW = 8.0 ! [s] !también el tiempo máximo
      real*8, parameter :: t0 = 0.0 ! [s]
      integer, parameter :: dwn_N = 256
      integer, parameter :: frec_N = 256
      integer, parameter :: time_N = 1024
      real*8, parameter :: Qq = 1000
      integer, parameter :: Npixeles_X = 80
      integer, parameter :: Npixeles_Z = 80
      real*8, parameter :: cotaX = 5.0 ![m]
      real*8, parameter :: cotaZ = 5.0 ![m]
      real, parameter :: escala = 8.0
      
      logical, parameter :: impEspectros = .false.
      logical, parameter :: impSismogramas = .false.
      logical, parameter :: impPelicula = .true.
      ! ———————————————————————————
      real*8 :: DK,Dt,OME,OMEI,FREC
      complex*16 :: cOME
      real*8, dimension(2*dwn_N) :: k_
      complex*16, dimension(time_N) :: t_,Uo,S
      
      real*8, parameter :: PI = real(4.0d0*ATAN(1.0d0),8)
      complex*16, parameter :: UI = cmplx(0.0d0,1.0d0,8), &
                               UR = cmplx(1.0d0,0.0d0,8), &
                               Z0 = cmplx(0.0d0,0.0d0,8)
      ! tipos
      type t_campo
        complex*16 :: ux,uz,p
      end type t_campo
      
      type Punto2d
        real*8 :: x,z
      end type Punto2d
      
      type t_receptor
        type(Punto2d) :: center
        integer :: layer !0 liquido, 1 solido
        type(t_campo), dimension (2*dwn_N) :: campo
        type(t_campo), dimension (frec_N+1, Npixeles_X) :: pixeles
      end type t_receptor
      type(t_receptor), dimension(Npixeles_Z),target :: receptor
      
      real*8, dimension(Npixeles_X),target :: pixel_coordsX
      real*8, pointer :: z,x
      complex*16, dimension(:), pointer :: ptCampo

      !geometria
      ! del estrato liquido
        real*8 :: h, zF, zR
        complex*16 :: RHO,ETA,alfaF
      ! del semiespacio
        complex*16 :: ALFA,BETA,RHOs,AMU,LAMBDA
        complex*16 :: GAMMA,NU
      
      ! contadores
      integer :: i,j,ik,iOME
      
      ! auxiliares
      complex*16, dimension(4,4) :: M
      complex*16, dimension(4) :: A,B,work
      complex*16 :: eietah, b1y3
      integer, dimension(4) :: ipiv
      integer :: info
      integer :: lwork,status,n_maxtime
      real*8 :: signo
      real :: factor
      complex*16, dimension (2*dwn_N) :: auxKp,auxKx,auxKz
      type(t_campo), dimension(Npixeles_X,Npixeles_Z,time_N),target :: fotogramas
      character(LEN=100) :: titleN,yAx, CTIT
      
      real, dimension(:,:,:), allocatable :: xvmat,yvmat
      real, dimension(Npixeles_X) :: xpray
      real, dimension(Npixeles_Z) :: ypray
      
      real :: minX,maxX,minY,maxY,xstep,zstep
      ! DATOS
      
      ! geometría estrato liquido
      zF = -1.0 ![m]
      zR = 1.0  ![m]
      h = zR - zF !espesor del estrato liquido
      RHO = 0.5 * UR ![T/m3]
      alfaF = 1.0 * UR ![m/s]
      
      ! semiespacio
      ALFA = 2.0 * UR ![m/s]
      BETA = 1.0 * UR ![m/s]
      RHOs = 1.0 * UR ![T/m3]
      AMU = RHOs * BETA**2.0
      LAMBDA = RHOs * ALFA**2.0 - 2.0_8 * AMU
      
      DK = 2.0_8*pi*(DFREC * frec_N)/(beta * dwn_N)* 1.1
      Dt = (1.0) / (real(time_N) * DFREC)
      ! --------------------
      
       write(6,'(/,a,F9.7)') "DK = ",DK
       write(6,'(a,F9.7)') "delta X = ", real(pi / (dwn_N * DK),4)
       write(6,'(a,EN12.2E2,a)') "L = ",2*pi/DK, "m"
       write(6,'(a,I0,/)') 'N. wavenumbers: ', dwn_N
       write(6,'(a,F9.7,a,F15.5)') "dt = ",Dt," Tmax=",Dt* time_N
       write(6,'(a,I0,a,F8.4,a,F12.4,a,/,a,F8.1,/)') & 
           'N. frequencies: ', frec_N,'  @',DFREC,'Hertz :. Fmax = ', & 
           frec_N*DFREC,'Hertz','Atenuation Q = ',Qq
           
      ! --------------------
      ! receptores
      do i=1, Npixeles_Z
        receptor(i)%center%z = zF + (i-1) * cotaZ/(Npixeles_Z-1)
        receptor(i)%center%x = 0.0_8
        receptor(i)%layer = 0
        if (receptor(i)%center%z .ge. zR) receptor(i)%layer = 1
!       print*,"z",i," ",receptor(i)%center%z,"  L=",receptor(i)%layer 
        receptor(i)%campo(1:2*dwn_N)%ux = Z0
        receptor(i)%campo(1:2*dwn_N)%uz = Z0
        receptor(i)%campo(1:2*dwn_N)%p = Z0
        receptor(i)%pixeles(1: frec_N, 1: Npixeles_X)%ux = Z0
        receptor(i)%pixeles(1: frec_N, 1: Npixeles_X)%uz = Z0
        receptor(i)%pixeles(1: frec_N, 1: Npixeles_X)%p = Z0
      end do!
      do i= 0,Npixeles_X-1
        pixel_coordsX(i+1) = -cotaX / 2.0_8 + i * cotaX/(Npixeles_X-1)
!       print*,"x",i+1," ",pixel_coordsX(i+1)
      end do
!     stop
      ! vector k_
      ! cero
        k_(1) = real(dk * 0.01,8)
      ! positivos
      do ik = 2, dwn_N+1
        k_(ik) = real(ik-1,8) * dk
      end do
      ! negativos
      do ik = dwn_N+2,2*dwn_N
        k_(ik) = (ik - 2* dwn_N - 1) * dk
      end do
!     do ik=1,2*dwn_N
!      print*,ik,k_(ik)
!     end do;stop
      ! vector t_
      t_(1: time_N) = z0
      t_(1) = exp(cmplx(0.0,-0.01* dfrec * t0*(2*pi),8))
      do i = 2,frec_N+1
        t_(i) = exp(cmplx(0.0,-(i-1)* dfrec * t0*(2*pi),8))
      end do!
      do i = time_N-frec_N+2,time_N
        t_(i) = exp(cmplx(0.0,-(i-time_N-1)* dfrec * t0*(2*pi),8))
      end do
      
      factor = sqrt(real(2*dwn_N))
      
      call ricker(Uo,Dt,pi,time_N)
      ! frec loop 
      do iOME = 1, frec_N+1
        write(6,'(A)', ADVANCE = "NO") "X"
        FREC = DFREC*real(iOME-1) !Hz
        if (iOME .eq. 1)  FREC = 0.1_8 * DFREC
        OME=2.0*PI*FREC !rad/s
        OMEI = - 2.0*PI/TW
        cOME = CMPLX(OME, OMEI,8)
        cOME = cOME * cmplx(1.0, -1.0/2.0/Qq,8) !histeretic damping
        

        ! k loop 
        do ik = 1,2*dwn_N
        gamma = sqrt(cOME**2.0_8/ALFA**2.0_8 - k_(ik)**2.0_8)
           nu = sqrt(cOME**2.0_8/BETA**2.0_8 - k_(ik)**2.0_8)
          eta = sqrt(cOME**2.0_8/alfaF**2.0_8 - k_(ik)**2.0_8)
          if(aimag(gamma).gt.0.0)gamma = conjg(gamma)
          if(aimag(nu).gt.0.0)nu= conjg(nu)
          if(aimag(eta).gt.0.0)eta= conjg(eta)
          
          eietah = exp(-UI*eta*h)
        !   matriz

        ! p = 0
        M(1,1) = UR
        M(1,2) = eietah
        M(1,3) = Z0
        M(1,4) = Z0
        
        ! uz^f = uz^s
        M(2,1) = (-UI * eta) / (rho * come**2) * eietah
        M(2,2) = (UI * eta) / (rho * come**2)
        M(2,3) = UI * gamma
        M(2,4) = UI * k_(ik)
        
        ! sigma_zz = -p
        M(3,1) = eietah
        M(3,2) = UR
        M(3,3) = -lambda * k_(ik)**2 - gamma**2 * (lambda + 2.0 * amu)
        M(3,4) = - 2.0 * amu * k_(ik) * nu
        
        ! simga_zx = 0
        M(4,1) = Z0
        M(4,2) = Z0
        M(4,3) = -2.0 * amu * k_(ik) * gamma
        M(4,4) = amu * (nu**2 - k_(ik)**2)
        
        ! terminos de fuente
        b1y3 = - (rho * come**2 * dk) / (4.0 * pi * UI * eta) 
        
        B(1) = b1y3 * exp(-UI * eta * abs(zF))
        B(2) = dk / (4.0 * pi) * exp(-UI * eta * abs(zR)) !* sig(zr)
        B(3) = b1y3 * exp(-UI * eta * abs(zR))
        B(4) = Z0
        
        !   coeficientes
        lwork = 4*4
        call zgetrf(4,4,M,4,ipiv,info)
        if(info .ne. 0) stop "Problem at LU factorization of matrix "
        call zgetri(4,M,4,ipiv,work,lwork,info)
        if(info .ne. 0) stop "Problem at inverse of matrix "
        A = matmul(M,B)
        !       call showMNmatrixZ(4,1,A,"  A  ",6)
        !       stop "ok"
        !   campo en receptores 
        do i=1, Npixeles_Z
          if (receptor(i)%layer .eq. 0) then ! liquido
            z => receptor(i)%center%z
            ! ec (7)
            receptor(i)%campo(ik)%p = & 
            (rho * come**2.0 * dk)/(4.0*pi*UI*eta) * exp(-UI*eta*abs(z)) + &
            A(1) * exp(-UI*eta*(z-zF)) + A(2) * exp(UI*eta*(z-zR))
            
            ! ec (8)
            receptor(i)%campo(ik)%uz = &
            (- signo(z) * dk)/(4.0 * pi) * exp(-UI * eta * abs(z)) + & 
            (UI * eta)/(rho * come**2.0) * & 
            (A(2)*exp(UI*eta*(z-zR)) - A(1)*exp(-UI*eta*(z-zF)))
            
            receptor(i)%campo(ik)%ux = -(&
            (dk)/(4.0 * pi * eta) * exp(-UI * eta * abs(z)) + & 
            (UI)/(rho * come**2.0) * & 
            (A(2)*exp(UI*eta*(z-zR)) + A(1)*exp(-UI*eta*(z-zF)))&
            ) * k_(ik)
          elseif (receptor(i)%layer .eq. 1) then ! sólido
            z => receptor(i)%center%z
            ! ec (10)
            receptor(i)%campo(ik)%uz = &
            -UI * gamma * A(3) * exp(-UI * gamma *(z-zR)) &
            -UI * k_(ik) * A(4) * exp(-UI * nu *(z-zR))
            
            receptor(i)%campo(ik)%ux = &
            -UI * k_(ik) * A(3) * exp(-UI * gamma *(z-zR)) &
            +UI * nu * A(4) * exp(-UI * nu *(z-zR))
          end if
        end do !i campo en pixeles coord Z
        end do ! k loop
        
        ! El campo desplazamientos es par e impar en x
      ! -e^-ikx
          do i=1, Npixeles_Z
          if (receptor(i)%layer .eq. 0) then ! liquido.
            do j=1,Npixeles_X
              x => pixel_coordsX(j)
              
              ! recordar campo sin fase horizontal
              auxkp(1:2*dwn_N) = receptor(i)%campo(1:2*dwn_N)%p 
              auxkx(1:2*dwn_N) = receptor(i)%campo(1:2*dwn_N)%ux
              auxkz(1:2*dwn_N) = receptor(i)%campo(1:2*dwn_N)%uz
              !   fase horizontal
              do ik = 1, 2*dwn_N
                auxkp(ik) = auxkp(ik) * exp(-UI * k_(ik) * x)
                auxkx(ik) = auxkx(ik) * exp(-UI * k_(ik) * x)
                auxkz(ik) = auxkz(ik) * exp(-UI * k_(ik) * x)
              end do ! ik
              
             ! integrar k
             auxkp(1:2*dwn_N) = auxkp(1:2*dwn_N) * factor
             auxkx(1:2*dwn_N) = auxkx(1:2*dwn_N) * factor
             auxkz(1:2*dwn_N) = auxkz(1:2*dwn_N) * factor
             call fork(2*dwn_N,auxkp(1:2*dwn_N),+1)
             call fork(2*dwn_N,auxkx(1:2*dwn_N),+1)
             call fork(2*dwn_N,auxkz(1:2*dwn_N),+1)
             auxkp(1:2*dwn_N) = auxkp(1:2*dwn_N) / factor
             auxkx(1:2*dwn_N) = auxkx(1:2*dwn_N) / factor
             auxkz(1:2*dwn_N) = auxkz(1:2*dwn_N) / factor
             
             ! guardar resultado en frecuencia
             receptor(i)%pixeles(iome,j)%p = auxkp(1)
             receptor(i)%pixeles(iome,j)%ux = auxkx(1)
             receptor(i)%pixeles(iome,j)%uz = auxkz(1)
            end do ! j
          elseif (receptor(i)%layer .eq. 1) then ! sólido.
             do j=1,Npixeles_X
              x => pixel_coordsX(j)
              
              !   fase horizontal
              auxkx(1:2*dwn_N) = receptor(i)%campo(1:2*dwn_N)%ux
              auxkz(1:2*dwn_N) = receptor(i)%campo(1:2*dwn_N)%uz
              do ik = 1,2*dwn_N
                auxkx(ik) = auxkx(ik) * exp(-UI * k_(ik) * x)
                auxkz(ik) = auxkz(ik) * exp(-UI * k_(ik) * x)
              end do ! ik
              
             ! integrar k
             auxkx(1:2*dwn_N) = auxkx(1:2*dwn_N) * factor
             auxkz(1:2*dwn_N) = auxkz(1:2*dwn_N) * factor
             call fork(2*dwn_N,auxkx(1:2*dwn_N),+1)
             call fork(2*dwn_N,auxkz(1:2*dwn_N),+1)
             auxkx(1:2*dwn_N) = auxkx(1:2*dwn_N) / factor
             auxkz(1:2*dwn_N) = auxkz(1:2*dwn_N) / factor
             
             ! guardar resultado en frecuencia
             receptor(i)%pixeles(iome,j)%ux = auxkx(1)
             receptor(i)%pixeles(iome,j)%uz = auxkz(1)
            end do ! j
          end if
          end do ! i coord Z
      end do !frec loop
      
      
      ! espectros 
      call system("mkdir outs")
      CALL chdir("outs")
      
      if (impEspectros) print*,"imprimiendo espectros"
      call system("mkdir espectros")
      call chdir("espectros",status)
      if (status .eq. 0) call system("rm *.*")
        
      fotogramas(1:Npixeles_X,1:Npixeles_Z,1:time_N)%p = z0
      fotogramas(1:Npixeles_X,1:Npixeles_Z,1:time_N)%ux = z0
      fotogramas(1:Npixeles_X,1:Npixeles_Z,1:time_N)%uz = z0
      
      do i=1, Npixeles_Z
        if (receptor(i)%layer .eq. 0) then ! liquido.
          z => receptor(i)%center%z
          do j=1,Npixeles_X
          x => pixel_coordsX(j)
!         print*,"(",x,z,")"
          !p
          ptCampo => receptor(i)%pixeles(1:frec_N+1,j)%p
          ! crepa
          fotogramas(i,j,1:frec_N+1)%p = ptCampo(1:frec_N+1)
          fotogramas(i,j,time_N-frec_N+2:time_N)%p = conjg(ptCampo(frec_N:2:-1))
          fotogramas(i,j,1:time_N)%p = fotogramas(i,j,1:time_N)%p * t_
          fotogramas(i,j,1:time_N)%p = fotogramas(i,j,1:time_N)%p * Uo
          ptCampo => fotogramas(i,j,1:time_N)%p
          if (impEspectros) then
            write(yAx,'(a)') '$p__$ [T/m2]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               'f_p__[',j,',',i,'].pdf' 
            call plotXYcomp(ptCampo(1: frec_N+1),real(DFREC,4), frec_N+1,titleN, & 
            'frec[hz] ',yAx, CTIT ,1200,800,0.0)
          end if
          
          !ux
          ptCampo => receptor(i)%pixeles(1:frec_N+1,j)%ux
          ! crepa
          fotogramas(i,j,1:frec_N+1)%ux = ptCampo(1:frec_N+1)
          fotogramas(i,j,time_N-frec_N+2:time_N)%ux = conjg(ptCampo(frec_N:2:-1))
          fotogramas(i,j,1:time_N)%ux = fotogramas(i,j,1:time_N)%ux * t_
          fotogramas(i,j,1:time_N)%ux = fotogramas(i,j,1:time_N)%ux * Uo
          ptCampo => fotogramas(i,j,1:time_N)%ux
           if (impEspectros) then
            write(yAx,'(a)') '$ux_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               'f_ux_[',j,',',i,'].pdf' 
            call plotXYcomp(ptCampo(1: frec_N+1),real(DFREC,4), frec_N+1,titleN, & 
            'frec[hz] ',yAx, CTIT ,1200,800,0.0)
           end if
            
          !uz
          ptCampo => receptor(i)%pixeles(1:frec_N+1,j)%uz
          ! crepa
          fotogramas(i,j,1:frec_N+1)%uz = ptCampo(1:frec_N+1)
          fotogramas(i,j,time_N-frec_N+2:time_N)%uz = conjg(ptCampo(frec_N:2:-1))
          fotogramas(i,j,1:time_N)%uz = fotogramas(i,j,1:time_N)%uz * t_
          fotogramas(i,j,1:time_N)%uz = fotogramas(i,j,1:time_N)%uz * Uo
          ptCampo => fotogramas(i,j,1:time_N)%uz
           if (impEspectros) then
            write(yAx,'(a)') '$uz_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               'f_uz_[',j,',',i,'].pdf' 
            call plotXYcomp(ptCampo(1: frec_N+1),real(DFREC,4), frec_N+1,titleN, & 
            'frec[hz] ',yAx, CTIT ,1200,800,0.0)
            end if
            end do ! j
        elseif (receptor(i)%layer .eq. 1) then ! sólido.
          z => receptor(i)%center%z
          do j=1,Npixeles_X
          x => pixel_coordsX(j)
!         print*,"(",x,z,")"
          
          !ux
          ptCampo => receptor(i)%pixeles(1:frec_N+1,j)%ux
          ! crepa
          fotogramas(i,j,1:frec_N+1)%ux = ptCampo(1:frec_N+1)
          fotogramas(i,j,time_N-frec_N+2:time_N)%ux = conjg(ptCampo(frec_N:2:-1))
          fotogramas(i,j,1:time_N)%ux = fotogramas(i,j,1:time_N)%ux * t_
          fotogramas(i,j,1:time_N)%ux = fotogramas(i,j,1:time_N)%ux * Uo
          ptCampo => fotogramas(i,j,1:time_N)%ux
           if (impEspectros) then
            write(yAx,'(a)') '$ux_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               'f_ux_[',j,',',i,'].pdf' 
            call plotXYcomp(ptCampo(1: frec_N+1),real(DFREC,4), frec_N+1,titleN, & 
            'frec[hz] ',yAx, CTIT ,1200,800,0.0)
            end if
            
          !uz
          ptCampo => receptor(i)%pixeles(1:frec_N+1,j)%uz
          ! crepa
          fotogramas(i,j,1:frec_N+1)%uz = ptCampo(1:frec_N+1)
          fotogramas(i,j,time_N-frec_N+2:time_N)%uz = conjg(ptCampo(frec_N:2:-1))
          fotogramas(i,j,1:time_N)%uz = fotogramas(i,j,1:time_N)%uz * t_
          fotogramas(i,j,1:time_N)%uz = fotogramas(i,j,1:time_N)%uz * Uo
          ptCampo => fotogramas(i,j,1:time_N)%uz
           if (impEspectros) then
            write(yAx,'(a)') '$uz_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               'f_uz_[',j,',',i,'].pdf' 
            call plotXYcomp(ptCampo(1: frec_N+1),real(DFREC,4), frec_N+1,titleN, & 
            'frec[hz] ',yAx, CTIT ,1200,800,0.0)
            end if
            end do ! j
          end if
      end do ! i coord Z
      call chdir("..")
      

      ! sismogramas 
      factor = sqrt(1.0* time_N)
        if (impSismogramas) print*,"imprimiendo sismogramas"
        call system("mkdir sismograma")
        call chdir("sismograma",status)
        if (status .eq. 0) call system("rm *.*")
      
      !tiempo maximo para graficar
         n_maxtime = int(TW/dt)
         if(TW .lt. dt) n_maxtime = 2*frec_N
         if(TW .gt. time_N * real(dt,4)) n_maxtime = time_N
      
      do i=1, Npixeles_Z
        if (receptor(i)%layer .eq. 0) then ! liquido.
          z => receptor(i)%center%z
          do j=1,Npixeles_X
          x => pixel_coordsX(j)
!         print*,"(",x,z,")"
          !p
          S = fotogramas(i,j,1:time_N)%p
          S = S * factor
          call fork(time_N,S,+1)
          S = S / factor
!         S = S * dfrec * sqrt(1.0* time_N)  
          S = S * exp(- OMEI * Dt*((/(i,i=0, time_N-1)/)))
          fotogramas(i,j,1:time_N)%p = S
          if (impSismogramas) then
            write(yAx,'(a)') '$p__$ [T/m2]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               't_p__[',j,',',i,'].pdf' 
            call plotXYcomp(S(1:n_maxtime),real(Dt,4), n_maxtime,titleN, & 
            'time[sec]',yAx, CTIT ,1200,800,0.0)
          end if
          
          !ux
          S = fotogramas(i,j,1:time_N)%ux
          S = S * factor
          call fork(time_N,S,+1)
          S = S / factor
          S = S * dfrec * sqrt(1.0* time_N) 
          S = S * exp(- OMEI * Dt*((/(i,i=0, time_N-1)/)))
          fotogramas(i,j,1:time_N)%ux = S
          if (impSismogramas) then
            write(yAx,'(a)') '$ux_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               't_ux_[',j,',',i,'].pdf' 
            call plotXYcomp(S(1:n_maxtime),real(Dt,4), n_maxtime,titleN, & 
            'time[sec]',yAx, CTIT ,1200,800,0.0)
          end if
            
          !uz
          S = fotogramas(i,j,1:time_N)%uz
          S = S * factor
          call fork(time_N,S,+1)
          S = S / factor
          S = S * dfrec * sqrt(1.0* time_N) 
          S = S * exp(- OMEI * Dt*((/(i,i=0, time_N-1)/)))
          fotogramas(i,j,1:time_N)%uz = S
          if (impSismogramas) then
            write(yAx,'(a)') '$uz_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               't_uz_[',j,',',i,'].pdf' 
            call plotXYcomp(S(1:n_maxtime),real(Dt,4), n_maxtime,titleN, & 
            'time[sec]',yAx, CTIT ,1200,800,0.0)
          end if
            end do ! j
        elseif (receptor(i)%layer .eq. 1) then ! sólido.
          z => receptor(i)%center%z
          do j=1,Npixeles_X
          x => pixel_coordsX(j)
!         print*,"(",x,z,")"
          
          !ux
          S = fotogramas(i,j,1:time_N)%ux
          S = S * factor
          call fork(time_N,S,+1)
          S = S / factor
          S = S * dfrec * sqrt(1.0* time_N) 
          S = S * exp(- OMEI * Dt*((/(i,i=0, time_N-1)/)))
          fotogramas(i,j,1:time_N)%ux = S
          if (impSismogramas) then
            write(yAx,'(a)') '$ux_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               't_ux_[',j,',',i,'].pdf' 
            call plotXYcomp(S(1:n_maxtime),real(Dt,4), n_maxtime,titleN, & 
            'time[sec]',yAx, CTIT ,1200,800,0.0)
          end if
            
          !uz
          S = fotogramas(i,j,1:time_N)%uz
          S = S * factor
          call fork(time_N,S,+1)
          S = S / factor
          S = S * dfrec * sqrt(1.0* time_N) 
          S = S * exp(- OMEI * Dt*((/(i,i=0, time_N-1)/)))
          fotogramas(i,j,1:time_N)%uz = S
          if (impSismogramas) then
            write(yAx,'(a)') '$uz_$ [m]'
            write(CTIT,'(a,F7.2,a,F7.2,a)')'(', x,' , ',z,')'
            write(titleN,'(a,I0,a,I0,a)') & 
               't_uz_[',j,',',i,'].pdf' 
            call plotXYcomp(S(1:n_maxtime),real(Dt,4), n_maxtime,titleN, & 
            'time[sec]',yAx, CTIT ,1200,800,0.0)
          end if
            end do ! j
          end if
      end do ! i coord Z
      call chdir("..")

      ! pelicula con los desplazamientos
      if (impPelicula) then
        print*,'imprimiendo fotogramas'
        call system("mkdir fotogramas")
        call chdir("fotogramas",status)
        if (status .eq. 0) call system("rm *.*")
        
      allocate(xvmat(Npixeles_X, Npixeles_Z,n_maxtime))
      allocate(yvmat(Npixeles_X, Npixeles_Z,n_maxtime))
      do j=1,Npixeles_X
        x => pixel_coordsX(j)
        xpray(j) = real(x,4)
      end do!
      
      do i=1, Npixeles_Z
        z => receptor(i)%center%z
        ypray(i) = real(z,4)
        do j=1,Npixeles_X
          xvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime)%ux * escala,4)
          yvmat(j,i,1:n_maxtime) = real(fotogramas(i,j,1:n_maxtime)%uz * escala,4)
        end do
      end do
      
      minx = minval(xpray)
      maxx = maxval(xpray)
      miny = minval(ypray)
      maxy = maxval(ypray)
      xstep = real(abs(xpray(1))/3.0,4)
      zstep = real(max(ypray(2)-ypray(1),real(int( (maxY-minY) / 10.))),4)
      
      CALL METAFL('PNG')
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL PAGE (int(3100,4),int(2400,4))
      call imgfmt('RGB')
      call winsiz(int(1100,4),int(800,4)) !1200,800
      CALL SCRMOD('REVERS') !fondo blanco
      
      do i=1,n_maxtime
      write(titleN,'(a,I0,a)') 'foto_',i,'.png'
      CALL SETFIL(trim(titleN))
      CALL DISINI()
      CALL BMPFNT ('SIMPLEX')
           !the position of an axis system.
      CALL axspos (int(300,4) ,int(2200,4)) ! Lower left corner
      call axslen (int(2000,4), int(2000,4)) !size of the axis system.
      call name('X [m]','X')
      call name('Z [m]','Y')
      
      call graf(real(minX,4),real(maxX,4),real(minX,4),real(xstep,4), & 
                 real(maxY,4),real(minY,4),real(maxY,4),real(-zstep,4))
      
      call height(50) ! de los caracteres
      call vecclr(-2) ! color de las puntas de flecha activado
      
      CALL VECOPT(real(1.0,4),'SCALE')
      
      call vecmat(xvmat(:,:,i),yvmat(:,:,i),Npixeles_X,Npixeles_Z,xpray,ypray,int(4221))
      
      call color ('FORE')
      call rline(real(minx,4),real(zR,4),real(maxx,4),real(zR,4))
      
      call errmod ("all", "off")
      call disfin
      end do ! i=1,n_maxtime
      
      write(titleN,'(a)')'ffmpeg -i foto_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p video.mp4'
      call system(trim(titleN))
      call chdir("..")
      call system('cp fotogramas/video.mp4 video.mp4')
      call system('rm fotogramas/video.mp4')
      
      end if ! impPelicula
      write(6,*)"hello"
      end program
      

      function signo(z)
      real*8 :: signo,errT
      real*8, intent(in) :: z
      errT = 0.0001_8
      if (abs(z) > errT ) then
            signo = real(z / ABS(z))
          else
            signo = 0.0_8
          end if
      end function
      
      subroutine showMNmatrixZ(m,n,MAT,name,outpf)
      integer, intent(in) :: m,n,outpf
      complex*16, dimension(m,n), intent(in) :: MAT
      integer :: i,j 
      character(LEN=5), intent(in) :: name
      
      write(outpf,'(A)') trim(name)
      do i = 1,m
        do j = 1,n
          write(outpf,'(A,E15.5,A)',advance='no') "(",REAL(MAT(i,j)),","
          write(outpf,'(E15.5,A)',advance='no') AIMAG(MAT(i,j)),"i) "
        end do
        write(outpf,'(A)',advance='yes')''
!       write(outpf,'(A)',advance='yes')''
      end do
      end subroutine!
      
      SUBROUTINE FORK(LX,CX,SIGNI)
      implicit none
      integer, intent(in) :: LX,SIGNI
      COMPLEX*16 :: CARG,CW,CTEMP 
      complex*16,intent(inout) :: CX(LX)
      real*8, parameter :: pi = 4.*ATAN(1.)
      real*8 :: SC
      integer :: i,j,m,istep,l
      
      J=1
      SC=DSQRT(real(1.0,8)/real(LX,8))
      DO 30 I=1,LX
      IF(I > J)GO TO 10
      CTEMP=CX(J)*cmplx(SC,0.0,8)
      CX(J)=CX(I)*cmplx(SC,0.0,8)
      CX(I)=CTEMP
   10 M=LX/2
   20 IF(J <= M)GO TO 30
      J=J-M
      M=M/2
      IF(M >= 1)GO TO 20
   30 J=J+M
      L=1
   40 ISTEP=2*L
      DO 50 M=1,L
      CARG=cmplx(0.0,(pi*real(SIGNI*(M-1)))/real(L),8)  
      CW=EXP(CARG)
      DO 50 I=M,LX,ISTEP
      CTEMP=CW*CX(I+L)
      CX(I+L)=CX(I)-CTEMP
   50 CX(I)=CX(I)+CTEMP
      L=ISTEP
      IF(L < LX)GO TO 40
      RETURN
      END subroutine fork
      
      subroutine plotXYcomp(y_in,Dt,n,titleN,xAx,yAx,CTIT,W,H,ma)
      ! (Uo,Dt,size(Uo),'FIGURE_NAME.pdf','time[sec]','amplitude',1200,800) 
      USE DISLIN
!     use glovars, only : pi
      implicit none
      real, intent(in)                              :: Dt,ma
      integer, intent(in)                           :: n,H,W
      character(LEN=9), intent(in)                  :: xAx
      character(LEN=100), intent(in)                :: yAx
      character(LEN=100)                            :: titleN
      COMPLEX*16, DIMENSION(n), intent(in) :: y_in
      complex,    dimension(n)             :: y
      
      real, dimension(n) :: x
      real maxY,minY,xstep,ystep,maxYc,minYc,val
      integer :: i,nPow10x,nPow10y,signo
!     integer :: Modo,nx,ny
      character(LEN=100) :: dumb
      CHARACTER(LEN=30) :: CBUF
      character(LEN=100) :: CTIT
!     integer*4 :: lentitle
!     character(LEN=100),parameter :: f1='(F50.16,2x,F50.16)'
      
      
!     allocate(x(n))
!     print*,size(y)
      
      ! para que se vea el texto en los ejes si es muy pequeño en formato F3.1
      nPow10x = 0
      if (Dt *(n-1) < 0.6) then
        do i = 1,10
          if (Dt *(n-1)*(10.0**i) > 1.0) then
            exit
          end if 
        end do 
        nPow10x = i
        
      elseif (Dt * (n-1) > 6000.) then  
        do i = 1,10
          if (Dt *(n-1)*(10.0**(-i)) < 1000.0) then
            exit
          end if
        end do
        nPow10x = -i
      end if
      DO i = 1,n
        x(i) = Dt*(i-1)*(10.0**(nPow10x))
        y(i) = cmplx(real(y_in(i)),aimag(y_in(i)),4) 
      END DO
      
      if (ma .eq. 0) then  
      minY=MINVAL(real(y(:)),1)
!     write(6,*)"MinReal Val= ",minY
      maxY=MAXVAL(real(y(:)),1)
!     write(6,*)"MaxReal Val= ",maxY
      minYc=MINVAL(aimag(y(:)),1)
!     write(6,*)"MinComplex Val= ",minYc
      maxYc=MAXVAL(aimag(y(:)),1)
!     write(6,*)"MaxComplex Val= ",maxYc
      minY =MIN(minY,minYc)
      maxY =MAX(maxY,maxYc)
      else
      maxy = ma
      miny = -ma
      end if
      
      val = max(abs(miny),abs(maxy))
      signo = 1
      if (miny*maxy < 0) then
        if (val - maxy > 0.0001) signo = -1
      else
        if (maxy < 0.0)  signo = -1
      end if
      nPow10y = 0
      do i = 1,10
!       if (
      end do
      
!     print*,"plotting"
! Dislin plotting routines ...
      CALL METAFL('PDF') !define intended display XWIN,PS,EPS,PDF
!     write(titleN,'(a,a)') trim(titleN),'.eps' 
!     titleN = trim(titleN)
!     titleN = trim(titleN)
      CALL SETFIL(trim(adjustl(titleN)))
!     print*,"file: ",trim(adjustl(titleN))
      call filmod('DELETE') ! para sobreescribir el archivo
      CALL SETPAG('DA4P')
      CALL PAGE(int(W+1200,4),int(H+350,4))
      CALL PAGMOD('NONE')
      CALL DISINI() ! dislin initialize
      call errmod ("all", "off") 
!     CALL PAGERA() ! plots a page border
      CALL COMPLX ! sets a complex font
      CALL TEXMOD ('ON') ! latex!!
      CALL HWFONT()
      CALL axspos (int(380,4) ,int(H+100,4)) !the position of an axis system. Lower left corner
      call axslen (int(W+650,4) ,int(H,4)) !size of the axis system.
      if (nPow10x .ne. 0) then
       write(dumb,'(a,a,I0)') trim(xAx),'x10^',(nPow10x *(-1))
       call name(trim(dumb),'X')
      else
       call name(trim(xAx),'X') 
      end if
      call name(trim(yAx),'Y') 
      call labels('EXP','Y')
      call labdig(int(2,4),'X') !number of decimal places for labels
      call labdig(int(1,4),'Y')
      call ticks (int(10,4) ,'X') 
      call ticks (int(5,4) ,'Y') 
!     call titlin ( titleN , 1 )
      xstep = x(n)/6.0 ! incremen for labels
      ystep = (maxY-minY)/6.0
      
      if (maxY * minY < 0.0) then
       maxY = max(abs(minY),abs(maxY))
       minY = -1.0 * maxY
       ystep = (maxY-minY)/2.0 !solo 3 etiquetas
       maxY = maxY - mod(maxY,ystep)
       minY = -1.0 * maxY
       ystep = (maxY-minY)/6.0
       maxY = maxY + ystep
       minY = -1.0 * maxY
      end if
      
      
      
      
!     print*,"maxmin= ",maxy,miny
      call graf(real(x(1),4), & 
                real(x(n)+x(2),4), & 
                real(x(1),4), & 
                real(xstep,4), & 
                real(minY,4), & 
                real(maxY,4), & 
                real(minY,4), & 
                real(ystep,4)) 
      
!     call title 
      call color ('RED') 
      call curve(real(x,4) ,real(y,4) ,int(n,4))
      call color('BLUE')
      call curve(real(x,4), real(aimag(y),4), int(n,4))
      call color ('FORE') 
      call dash() !sets a dashed line style
      call xaxgit() !plots the line Y=0
      
      call legini(CBUF,int(2,4),int(20,4))
 !     nx = nxposn(x(n)*n + x(n)*n / 20.)
 !     ny = nyposn(minY + (maxY-minY)*0.7)
 !     print*,nx
 !     print*,ny
      if (maxval(abs(aimag(y))) .gt. 0.05*maxval(abs(real(y)))) then
      call legpos(int(1840,4),int(720,4))
      write(dumb,'(a)') 'Re(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(1,4))
      write(dumb,'(a)') 'Im(z)'
 !     print*,dumb
      call leglin(CBUF,dumb,int(2,4))
      call legtit('') ! or '' for nothing
      call legend(CBUF,int(2,4))
      end if
      
!     write(CTIT,'(a,ES11.4E2,a)') 'dt=',Dt,' seg'
!     lentitle = NLMESS(CTIT)
!     CALL MESSAG(CTIT,int((1700),4),int(200,4))
      CALL MESSAG (CTIT, int(900,4), int(110,4))
!     call errmod ("protocol", "off") !suppress dislin info
      call disfin()      
      
!     print*,'plotted ',trim(titleN)
      end subroutine plotXYcomp
      
      subroutine ricker(Uo,Dt,pi,time_N)
      implicit none
      
      integer, intent(in) :: time_N
      complex*16, dimension(time_N) :: Uo
      real*8, intent(in) :: Dt,pi
      
      real, parameter :: Ts = 0.5
      real, parameter :: Tp = 0.4
      
      integer :: i
      real*16 :: A
      
      Uo(1: time_N) = 0
      
      A = pi*(-Ts) / Tp
      Uo(1) = cmplx((A*A-0.5)* exp(- A * A),0.,8) 
      do i = 2, time_N/2+1
        A = pi*(Dt*(i-1)-Ts) / Tp
        A = A * A
        A = (A - 0.5_16) * exp(- A) 
        if ( abs(A) .gt. 0.0001_16) then
           Uo(i) = cmplx(real(A,8),0.0_8,8)
        end if
      end do
      
      call fork(size(Uo),Uo,-1)
      
      end subroutine ricker

