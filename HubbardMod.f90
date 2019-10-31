Module HubbardMod
!use omp_lib

CONTAINS

SUBROUTINE HubbardGreen(Ue,nAvgUPfinal,nAVGDWfinal,  atomsCam,nimpurezas, Uarray,Vbordarr, Eoarray, &
                            Ey,Ez, fileDatos, fileDatos2)
!1 sola camada !! encontra nAvgUP, nAvgDW
! saida:  nAvgUPfinal,nAVGDWfinal

Use Routines1
Use Routines2
Use Routines3
Use Routines6
Use Routines7
Use RedeMod
Use transporteMod

!REAL(8), PARAMETER :: Vy = 0.2d0, Vz=0.0d0, hzeeman =0.0d0 !hzeeman: contiene campoB fraco
    ! campoE esta en unidades talque a un voltaje=t(voltio~energia), una carga |qe|=1 en
    !entre dos placas separadas acc, el electron ganar? energia cinetica t (=1).
    ! acc en t?rminos de t=gamma0=1, vf=3ta/2
    ! la separacion entre camadas es 3,4Amstrong, la separacion entre atomos de carbono en una misma
    ! camada es 1,42Amstrong

IMPLICIT NONE

!CHARACTER (LEN=*), INTENT(IN):: endereco,fileEnergFerm, fileDatos, fileLDOSup, fileLDOSdown, filedosUP, filedosDOWN, &
!                             fileDatSitios, fileQcamadas, fileLDOSup2, fileLDOSdown2, fileTUP, fileTDW,fileDatos2

INTEGER, INTENT(IN) :: atomsCam, nimpurezas
REAL(8), dimension(atomsCam), INTENT(IN) :: Uarray,Vbordarr,Eoarray
!*********************************************************************************************************************
    REAL(8), dimension(atomsCam), INTENT(OUT) :: nAvgUPfinal,nAVGDWfinal,Ue
!                    REAL(8), dimension(atomsCam), INTENT(INOUT) :: nAvgUPfinal,nAVGDWfinal,Ue
!*********************************************************************************************************************
CHARACTER (LEN=*), INTENT(IN)::  fileDatos, fileDatos2 !,endereco

!REAL(8), INTENT(IN)          :: U, Vy
!REAL(8), INTENT(IN)          :: Ey, Ez, Flux
REAL(8), INTENT(IN)          :: Ey, Ez
REAL(8)            :: angleDeg, E !vuelven a ser calculados
!REAL(8)            :: Bfield

!                     U = 1.0d0,           &

!REAL(8), PARAMETER::  gamma0 = 1.0d0,&!-2.92d0,          & ! = VpzpzPi (NOTACION DE HERNAN, PARA h-BN)
!                      gamma1 = gamma0/10.0d0,   &
!                      !gamma1 = 0.0d0,   &
!                      !gamma3 = gamma1/3.0d0,    &
!                      gamma3 = 0.00d0,&
!                      !g2Vznh = gamma0/10.0d0,   &! segundos vizinhos
!                      g2Vznh = 0.0d0,           &! segundos vizinhos

REAL(8), PARAMETER::  acc    = 1.0d0,           &
                      az     = 3.4d0/1.42d0,    &

                      !Vy     = 0.0d0,           &
                      !Vz     = 0.0d0,           &
                      hzeeman= 0.0d0

!INTEGER, PARAMETER :: atomsCam = 16 ,             &!Build again!!
INTEGER, PARAMETER :: Camadas  = 1
INTEGER ::   atomosT !atomsCam: átomos por camada.


!!==------------------------------------------------------
!    !CHARACTER (LEN=*), PARAMETER :: conf= "ANTIFERRO"
!    !CHARACTER (LEN=*), PARAMETER :: conf= "FERRO"
!    CHARACTER (LEN=*), PARAMETER :: conf= "NCB"
!!==------------------------------------------------------


REAL(8), PARAMETER :: imagEn     =  1.0D-4,       &!!era 1.0D-4
                        Eni      = -8.5d0,        & !-8.0d0
                        Enf      = -Eni

INTEGER, PARAMETER :: Nintervalos   = 80,         & !80 !200 !400, 70 ,, 0.2
                        izero       = 5*Nintervalos,    &
                        Ndizimac    = -1!30!22!-1

!REAL(8) :: campoEy, campoEz

REAL(8),    PARAMETER :: Pi = 2*dACOS(0.0d0)
COMPLEX(8), PARAMETER :: I  = dcmplx(0,1) !=sqrt(-1)

COMPLEX(8), dimension(atomsCam, atomsCam) :: Gup, Gdown, GRup, GRdown
complex(8),    dimension(atomsCam, atomsCam) :: T00,T,TD
REAL(8),    dimension(atomsCam, atomsCam) :: temp,temp2,temp3,temp4
complex(8),    dimension(atomsCam, atomsCam) :: Ident
REAL(8),    dimension(atomsCam)          :: nAvgUP,nAVGDOWN,nup0,ndown0,Y,Z,V, UeOriginal
REAL(8),  allocatable, dimension(:)    :: Vhex
REAL(8),    dimension(2*izero)          :: Energia
REAL(8),    dimension(2*izero, atomsCam) :: LDOSup,LDOSdown

REAL(8)    :: En, cargaTup, cargaTdown, cargaT, &
              cTupL1, cTdownL1, Efermi0, Efermi, Qlayer1, Qlayer2, QUPl1, QDWl1, QUPl2, QDWl2
real(8)    :: errorMax,normaU,alpha,Vo,Ncharges
INTEGER    :: k,j,s,hu,ja, cont,natom,ni

COMPLEX(8), dimension(2*izero)  :: TransmisionUP, TransmisionDW !!CAMBIAR COMPLEX POR REAL CASO EL PROGRAMA ESTE OK
COMPLEX(8)                       :: a,TRANSMIup, TRANSMIdw !!CAMBIAR COMPLEX POR REAL CASO EL PROGRAMA ESTE OK

integer, parameter:: lunDat=20, lunEn=21, lunfilesOut=22,lunT00=100,lunT=101,lunTD=102,lunEnTotal=103

    real(8) :: t1,t2,t3,t4,t5,t6
    character(1) :: atm

    logical :: isFermiLevelMoving
    namelist /fermiLevel/ isFermiLevelMoving




V=0.0d0

atomosT  = atomsCam*1 !atomsCam: átomos por camada. Considerando 1 camad
natom=atomsCam
!_________________________________________________________________________
!-------------------------------------------------------------------------
open (unit=lunDat, file=fileDatos)

    !Flux = 0.0d-3,             & !Flux es adimensional, = FluxHexagono/FluxoMin, FluxMin=h/e
    !Flux = 33.0d0/1000.0d0,             &
!Bfield = Flux/(1.267d-5)     !& !! Bfield en TESLA!!
    !!!Dfluxo = Bfield*0.422d-5

!__________________ Inicializando Potencial Electrostatico _______________
!-------------------------------------------------------------------------
Ue=0.0d0
!____________________Campo Ey ____________________________________________
!-------------------------------------------------------------------------
!call CentrarEmYarmachair(Y, Ue, atomosT, Camadas, atomsCam, acc, Ey, lunDat)

call Ue_Ydirection(Y, Ue, natom, Ey, lunDat)

!call CentrarEmYZStrain(Y,Z,Ue, atomsCam, Ey,Ez, lunDat) !Sinal de U corregido
!! YA NO ES NECESARIO centrar en Z !!
!! *******   EL CAMPO "Ey,Ez" DEBE ESTAR EN ENERGY[eV]/ANGST ROM *******  !!!

!call CentrarEmY(Y, Ue, atomosT, Camadas, atomsCam, acc, Ey, lunDat)
!salida Y, Ue (de tamanho atomosT)
!Y: posiciones de los sitios en el eje Y. Sistema centrado
!Ue: potencial electrostatico en cada sitio
!******OUTPUT : Y entra como [alat], output como [Angstrom]
!******OUTPUT : Ue sai como [eV]



!___________________ Campo Ez _____________________________________________

!call CentrarEmZ(Z, Ue, atomosT, Camadas, atomsCam, az, Ez, lunDat)
!salida Z, Ue (de tama?o atomosT)
!Z: posiciones de los sitios en el eje Z. Sistema centrado
!Ue: potencial electrost?tico en cada sitio

!_________________________  Nanofita ZigZag  _____________________________
!-------------------------------------------------------------------------

  !!call ZigZag(T00, T, TD, Ident, gamma0, gamma1, gamma3, g2Vznh, atomosT, atomsCam, Camadas,Y,Bfield,acc)
  !call ZigZag(T00, T, TD, Ident, gamma0, gamma1, gamma3, g2Vznh, atomosT, atomsCam, Camadas)
  !!salida : T00, T, TD, Ident
     !**************************************************************************************************************
      call RedeLeonor(T00, T, TD, Ident, natom) ! da T00, T, TD en unidades eV !!! segun el file "hoppingsMod.f90"
      !!salida : T00, T, TD, Ident
     !**************************************************************************************************************



if(Ey.ne.0.d0) then
   angleDeg = (aTan(Ez/Ey))*(180/Pi) !esta en degrees, grados sexagesimales
else
   angleDeg = 90.d0 !esta en degrees, grados sexagesimales
endif

E        = sqrt(Ey**2  + Ez**2)

!call writeVal(Eoarray,Uarray, gamma0, gamma1, gamma3, g2Vznh, acc, az, angleDeg,E, Ey,Ez,hzeeman, &
!              atomsCam, Camadas, atomosT, imagEn, T00, T, TD, lunDat)
call writeValLeo(Eoarray,Uarray,  angleDeg,E, Ey,Ez, hzeeman, &
                    atomsCam, Camadas, atomosT, imagEn, lunDat)


!_________________________  Green ________________________________________
!-------------------------------------------------------------------------
Gup   = 0.0d0
Gdown = 0.0d0

!______________inicializando nAvgUP, nAvgDOWN_____________________________
!-------------------------------------------------------------------------
!call nAvgYaCalculado(nAvgUP,nAvgDOWN,  atomosT,fileDatSitios) !caso que los nAvg ya estén calculados
!call nAvgInicial(nAvgUP, nAvgDOWN, atomosT, Camadas, atomsCam, conf, lunDat)

!*********************************************************************************************************************
            call nAvgInicial(nAvgUP, nAvgDOWN, atomosT,nimpurezas, Camadas, atomsCam, lunDat,Y)
!*********************************************************************************************************************
!nAvgUP = nAvgUPfinal
!nAVGDOWN = nAVGDWfinal


!salida: nAvgUP, nAvgDOWN
!_______________________Energia __________________________________________
!-------------------------------------------------------------------------
Energia =  Energy(Eni, Enf, Nintervalos, izero, lunEn)
!obtiene los puntos Energia para los cuales se calcularán las propiedades
!se imprime en la salida.

!-------------------------------------------------------------------------
!_____________________ AUTOCONSISTENCIA __________________________________
!-------------------------------------------------------------------------
Efermi   = 0.0d0;
cont     = 1
hu       = 1

normaU = norm2(Uarray)
if (normaU < 1.0d0) then
    alpha = normaU/2.0d0
else
    alpha = 0.5d0
endif



!**********************************************************************************************
!open(unit=27, file='X_nAvg')
!    do j=1,natom
!        read(27,*) nAvgUP(j), nAvgDOWN(j)
!    enddo
!close(27)
!open(27, file='X_Efermi')
!    read(27,*) Efermi
!close(27)



!    do j=1,natom
!        if ( (j<11).or.(j>18)  ) Ue(j)=0.0d0
!        print*, Y(j),Vbordarr(j), Ue(j), Eoarray(j)+Vbordarr(j)+ Ue(j)
!    enddo
!stop

!**********************************************************************************************


!**********************************************************************************************
Ncharges=0.d0
open(unit=145, file='energiasSK')
    do  j=1,natom
        read(145,*) atm,t1,t2,t3,t4,t5,t6
        !    print*, atm,t1,t2,t3,t4,t5,t6
        !    stop
        if (atm=='B') Ncharges = Ncharges +0.0d0
        if (atm=='C') Ncharges = Ncharges +1.0d0
        if (atm=='N') Ncharges = Ncharges +2.0d0
!        if (atm=='B') Ncharges = Ncharges +1.0d0
!        if (atm=='C') Ncharges = Ncharges +1.0d0
!        if (atm=='N') Ncharges = Ncharges +1.0d0

    enddo
close(145)
print*, Ncharges
!stop
!**********************************************************************************************


UeOriginal = Ue
DO WHILE(cont==1)

    if(atomosT > 30) then
        open(unit=135,file='navg_temporal.dat')
            do j=1,atomosT
                write(135,*) nAvgUP(j), nAVGDOWN(j)
            enddo
        close(135)
    endif

!do hu=1,10
    nup0   = nAvgUP
    ndown0 = nAvgDOWN
    Efermi0= Efermi

    !_____________________________________________________________________
    Ue = UeOriginal

!    call PotencialLineas(V, nAvgUP,nAVGDOWN,atomsCam,nimpurezas,Y)
!    call VPuntualBlindado3(V, natom,nAvgUP,nAVGDOWN)
!    Ue = Ue + V

!     Ue = Ue + Ee_e(atomosT, nAvgUP, nAvgDOWN) ! ESTA DEFINIDO EN Routines7.f90

!**    call PoissonNanotubo(V,Vhex, nAvgUP,nAVGDOWN,atomsCam, Camadas, atomosT)
!    call Poisson2(V,Vhex, nAvgUP,nAVGDOWN,atomsCam, Camadas, atomosT)
    !salida: V, potencial em cada sitio (consideramos q_electron=+1, corregir, --> V tiene terminos de energia)
    !salida: Vhex
!**    Ue = Ue + V
    !_____________________________________________________________________

    !_____________________________________________________________________
    !------------------ calculo de <n> integrando no plano complexo ------


open(unit=141, file='myinputRibbon.in')
    read(141, fermiLevel) !se lee isFermiLevelMoving
close(141)

IF (isFermiLevelMoving) then
!*******************************************************************************
    call ocupacao(nAvgUP,nAvgDOWN, Efermi, T00,T,TD,Ident,Eoarray,Uarray,Vbordarr,atomosT,Ndizimac, &
                Ue, hzeeman,imagEn,Ncharges)
    ! "ocupacao" shifted Efermi and calculated new values of nAvgUP, nAvgDOWN
ELSE
!*******************************************************************************
!    call ocupacao(nAvgUP,nAvgDOWN, Efermi, T00,T,TD,Ident,Eoarray,Uarray,Vbordarr,atomosT,Ndizimac,&
!               Ue, hzeeman,imagEn,lunDat)
    print*, 'integrando en el plano complejo...'
    call IntegraPlanoCmplx(nAvgUP, nAvgDOWN, atomosT,hzeeman, Eoarray,Uarray,Vbordarr, Ue, Efermi, T00,T,TD,Ident, Ndizimac)
    print*, 'hemos simplificado para el caso en que Efermi no se mueve!'
!*******************************************************************************
ENDIF


    write(*,*),     'calculo de <n> integrando no plano complexo'
    write(lunDat,*),'calculo de <n> integrando no plano complexo'
    cargaT   = 0.0d0
    cargaTup = 0.0d0 ;    cargaTdown = 0.0d0
    cTdownL1 = 0.0d0 ;    cTupL1     = 0.0d0

    do j=0, Camadas-1
        write(*,*),      'Camada No', j+1
        write(lunDat,*), 'Camada No', j+1
        ja = j*atomsCam
            write(*,199) ' ----spinUP---- | ---spinDOWN--- | ---sSuma------   atomo,    iteracao : '
            write(lunDat,199) ' ----spinUP---- | ---spinDOWN--- | ---sSuma------   atomo,    iteracao : '
        do s=1,atomsCam
            write(*,200)       nAvgUP(s+ja),'|' ,nAvgDOWN(s+ja),'|',nAvgUP(s+ja) + nAvgDOWN(s+ja), s,hu
            write(lunDat,200)  nAvgUP(s+ja),'|' ,nAvgDOWN(s+ja),'|',nAvgUP(s+ja) + nAvgDOWN(s+ja), s,hu
            !write(lunDat,201)  nAvgUP(s+ja),nAvgDOWN(s+ja),nAvgUP(s+ja) + nAvgDOWN(s+ja), s*0.1,hu*0.1
            cargaTup   = cargaTup +   nAvgUP(s+ja)
            cargaTdown = cargaTdown + nAvgDOWN(s+ja)
        enddo
        if ((Camadas==2).and.(j==0)) then
            cTupL1   = cargaTup
            cTdownL1 = cargaTdown
        endif
    enddo
199 format(a73)
200 format(f14.8, a3, f14.8, a3, f14.8, i8, i8)
201 format(f14.8, a3, f14.8, a3, f14.8)
!201   format (80(f12.8,',')) !! >=atomosCam*Camadas!!multiplicar por el n?

    cargaT = cargaTup + cargaTdown




    Qlayer1=cTupL1+cTdownL1
    Qlayer2=(cargaTup-cTupL1)+(cargaTdown-cTdownL1)
    QUPl1  = cTupL1
    QDWl1  = cTdownL1
    QUPl2  = cargaTup-cTupL1
    QDWl2  = cargaTdown-cTdownL1

    write(*,199)      '---QTotal UP--- |---QTotal DW--- | ----QTOTAL----                        '
    write(*,201)       cargaTup, '|' ,cargaTdown,'|' , cargaT
    write(lunDat,199) '---QTotal UP--- |---QTotal DW--- | ----QTOTAL----                        '
    write(lunDat,201)  cargaTup, '|' ,cargaTdown,'|' , cargaT

   if (Camadas==2) then
        write(*,*),          'Carga da Camada Superior (N=2): Q-up, Q-down:         ', QUPl2,   QDWl2
        write(*,*),          'Carga da Camada Inferior (N=1): Q-up, Q-down:         ', QUPl2,   QDWl2
        write(*,*),          'total:                                      :         ', Qlayer1, Qlayer2

        write(lunDat,*),     'Carga da Camada Superior (N=2): Q-up, Q-down:         ', QUPl2,   QDWl2
        write(lunDat,*),     'Carga da Camada Inferior (N=1): Q-up, Q-down:         ', QUPl2,   QDWl2
        write(lunDat,*),     'total:                                      :         ', Qlayer1, Qlayer2
    endif

    write(*,*),      ' '
!    write(*,*),      'dQ ==', (cargaT -atomosT)!/atomosT
    write(*,*),      'dQ ==', (cargaT -Ncharges)!/atomosT
    write(*,*),      ' '
!    write(lunDat,*), 'dQ ==', (cargaT -atomosT)!/atomosT
    write(lunDat,*), 'dQ ==', (cargaT -Ncharges)!/atomosT
    write(lunDat,*), ' '

    !________________ criterio de convergencia _______________________________

    cont = continuarH(nAvgUP, nAvgDOWN, nup0, ndown0, atomosT, lunDat)
    !! introducir aqui  LA IDEA DE SERGIO ULLOA, DISMINUIR alpha (incrementar el peso de la ulima solucion) cuando comienza a converger !!!

    if (cont.ne.0) then
        nAvgUP      = (1.0d0-alpha)*nAvgUP     + alpha*nup0
        nAVGDOWN    = (1.0d0-alpha)*nAVGDOWN   + alpha*ndown0
        Efermi      = (1.0d0-alpha)*Efermi     + alpha*Efermi0
    endif

!    if (cont.ne.0) then
!            nAvgUP     = (nAvgUP   + nup0)/2.0d0
!            nAvgDOWN   = (nAvgDOWN + ndown0)/2.0d0
!            Efermi     = (Efermi   + Efermi0)/2.0d0
!    endif
    !-------------------------------------------------------------------------


    hu=hu+1
    IF (hu>100)  then
        errorMax = ErrorMaximo(nAvgUP, nAvgDOWN, nup0, ndown0, atomosT)
        if (errorMax .lt. 2.0d0/100) then
            write(*,*), 'PARECE QUE NO LOGRA DISMINUIR EL ERROR. Error =',errorMax*100,'%'
            write(lunDat,*), 'PARECE QUE NO LOGRA DISMINUIR EL ERROR. Error =',errorMax*100,'%'
            cont=0 !que no continue con las iteraciones, y termine de rodar.
        else
            WRITE(*,*),      'PARECE QUE NO CONVERGE ...'
            WRITE(lunDat,*), 'PARECE QUE NO CONVERGE ...'
!            GO TO 110
        endif
    END IF
    !if (hu==3) stop
    !_______________ caso que converja... ____________________________________

!STOP



    IF (cont==0) then
       call writeOccupation(nAvgUP, nAvgDOWN,atomosT) !Esta en Subrutinas3.f90
    ENDIF


!!! inicio de calculos que no es necesario para encontrar bandas.
    IF (cont==0) THEN
        call write1(nAvgUP, nAvgDOWN, atomosT, Eoarray,Uarray,Vbordarr, Ue,V, hzeeman, lunDat, &
                    cargaTup, cargaTdown, cargaT, Qlayer1, Qlayer2, QUPl1, QDWl1, QUPl2, QDWl2) !just write final data
!!        call write1(nAvgUP, nAvgDOWN, atomosT, Eoarray,Uarray,Vbordarr, Ue,V,Vhex, hzeeman, lunDat, &
!!                    cargaTup, cargaTdown, cargaT, Qlayer1, Qlayer2, QUPl1, QDWl1, QUPl2, QDWl2) !just write final data
        !---------------------------------------------------------------------
        !____ El proceso termino. Solo resta tener puntoas para graficar______
        !---------------------------------------------------------------------
!****** PODRIAS COMENTAR DESDE AQUI:!

         DO j=1, 2*izero     ! size(Energia)=2*izero
            En = Energia(j)

                !!GreenInicial(Gup,Gdown,atomosT,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, Energy, imagEn)
            call GreenInicial(Gup,Gdown,atomosT,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, En, imagEn)
            !Gup   = GreenUP(atomosT,   nAvgDOWN, hzeeman, U, Ue, En, imagEn)
            !Gdown = GreenDOWN(atomosT, nAvgUP,   hzeeman, U, Ue, En, imagEn)!


            GRup   = cGreenRnrmlzd(Gup,  T00,T,TD, Ident, atomosT,Ndizimac)
            GRdown = cGreenRnrmlzd(Gdown,T00,T,TD, Ident, atomosT,Ndizimac)

!        !=====================================================================================
            TRANSMIup = 0.0d0
            TRANSMIdw = 0.0d0
!            call cTransporte(TRANSMIup,TRANSMIdw,Gup,Gdown,GRup,GRdown, T00,T,TD,Ident,atomosT,Ndizimac)
!        !TRANSMISIONup y TRANSMISIONdw adquieren nuevos valores dependiente de la energia En!
!		!la cual fue trasnmitida a traves de los inicializadores Gup e Gdown!
!
!            TransmisionUP(j) = TRANSMIup
!            TransmisionDW(j) = TRANSMIdw
!
!            !-------------------------------------------------------
            do s=1,atomosT
                    LDOSup(j,s)   = -Dimag(GRup(s,s))/Pi
                    LDOSdown(j,s) = -Dimag(GRdown(s,s))/Pi
            enddo


          ENDDO
!
!!          call IntegraPlanoCmplx(nAvgUP, nAvgDOWN, atomosT,hzeeman, Eoarray,Uarray,Vbordarr, Ue, Efermi, T00,T,TD,Ident, Ndizimac)
!!do ni=1, atomosT
!!   print*, nAvgUP(ni), nAvgDOWN(ni), Uarray(ni)
!!enddo
!!stop
!
!!******** HASTA AQUI PODRIAS DESCOMENTAR/COMENTAR

          !=====================================================================================
          !calculando energia segundo Lloyd's formula
            !DeltaE(natoms, ni,Vo, Efermi, imagEn,hzeeman, Eoarray,Uarray,Vbordarr, Ue,nAvgUP,nAVGDOWN)

!            ni = 1
            !Vo = Eoarray(ni)-2.5  ! delta = Vo - Eoarray(ni)
            !WRITE(*,*) 'DeltaE = '
            !write(*,*)  DeltaE(atomsCam, ni,Vo, Efermi, imagEn,hzeeman, Eoarray,Uarray,Vbordarr, Ue,nAvgUP,nAVGDOWN)
!            write(*,*) ''
!            write(*,*) 'Energia total='
!        open(unit=lunEnTotal, file='Energiatotal_convergencia_energia')
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 1000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 2000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 3000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 4000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 5000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 6000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 7000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 8000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 9000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 10000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 11000,lunEnTotal)
!            call SystemEnergy(atomsCam,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
!                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, 12000,lunEnTotal)

!        close(lunEnTotal)




            write(*,*) '-=-=-=-=-=-=-=-='
         !=====================================================================================

            !---------------- solo para verificar:  ------------------------------
            !---------------- integrando en energia ------------------------------
          !!!$OMP PARALLEL DO
          !do s=1,atomosT
          !      nAvgUP(s)   = Gauss5Quad(  LDOSup(1:izero, s)  ,            Eni, 0.0d0, Nintervalos)
          !      nAvgDOWN(s) = Gauss5Quad(  LDOSdown(1:izero, s),            Eni, 0.0d0, Nintervalos)
          !      sUP(s)      = Gauss5Quad(  LDOSup(  izero+1:2*izero, s),  0.0d0, Enf,   Nintervalos)
          !      sDOWN(s)    = Gauss5Quad(  LDOSdown(izero+1:2*izero, s),  0.0d0, Enf,   Nintervalos)
          !end do
          !!!$OMP END PARALLEL DO!

          !do j=0, Camadas-1
          !      write(*,*), 'ocupacao <n> da camada No', j+1
          !      write(lunDat,*), 'ocupacao <n> da camada No', j+1
          !      ja=j*atomsCam
          !      do s=1,atomsCam
          !          write(*,*),      'spin UP:   ', nAvgUP(s+ja),   sUP(s+ja),  nAvgUP(s+ja)   + sUP(s+ja),  s,hu
          !          write(lunDat,*), 'spin UP:   ', nAvgUP(s+ja),   sUP(s+ja),  nAvgUP(s+ja)   + sUP(s+ja),  s,hu
         !          write(*,*),      'spin DOWN: ', nAvgDOWN(s+ja), sDOWN(s+ja),nAvgDOWN(s+ja) + sDOWN(s+ja),  s,hu
          !          write(lunDat,*), 'spin DOWN: ', nAvgDOWN(s+ja), sDOWN(s+ja),nAvgDOWN(s+ja) + sDOWN(s+ja),  s,hu
          !      enddo
          !enddo
          write(*,*),      ' '
          write(lunDat,*), ' '!

!!!*** si vas a descomentar lo siguiente, entonces descomenta tambien el loop anterior, donde se define LDOSup y LDOSdown!!!
 call SalidaDatos(Energia,Efermi,izero,atomsCam, Camadas, LDOSup, LDOSdown, &
                        TransmisionUP, TransmisionDW,fileDatos2)!

!!
!
!!
!
!!!-----------AGREoGADO-------------------------------------------------------
!!       !Efermi = Efermi ! se toma el Efermi ultimo que se calculó
!!       call IntegraPlanoCmplx(nAvgUP, nAvgDOWN, atomosT,hzeeman, U, Ue, Efermi, T00,T,TD,Ident, Ndizimac)
!!       nup0   = nAvgUP
!!       ndown0 = nAvgDOWN
!!            En = Efermi
!!            Gup   = GreenUP(atomosT,   nAvgDOWN, hzeeman, U, Ue, En, imagEn)
!!            Gdown = GreenDOWN(atomosT, nAvgUP,   hzeeman, U, Ue, En, imagEn)
!!            GRup   = cGreenRnrmlzd(Gup,  T00,T,TD, Ident, atomosT,Ndizimac)
!!            GRdown = cGreenRnrmlzd(Gdown,T00,T,TD, Ident, atomosT,Ndizimac)
!!            do s=1,atomosT
!!                    DensEfermiUP(s) = -Dimag(GRup(s,s))/Pi
!!                    DensEfermiDW(s) = -Dimag(GRdown(s,s))/Pi
!!            enddo
!!         DO j=1, 2*izero     ! size(Energia)=2*izero
!!            Efermi = Energia(j)
!!            call IntegraPlanoCmplx(nAvgUP, nAvgDOWN, atomosT,hzeeman, U, Ue, Efermi, T00,T,TD,Ident, Ndizimac)
!!            LDOSup(j,:)   = -nup0   + nAvgUP   + DensEfermiUP !para cada uno de los sitios!
!!            LDOSdown(j,:) = -ndown0 + nAvgDOWN + DensEfermiDW
!!         ENDDO
!!
!!        call SalidaDatos2(izero, fileLDOSup2, fileLDOSdown2, LDOSup, LDOSdown)
!!
!!!----------- FIN DE AGREGADO---------------------------------------------------------------!!
!
!
!    !_________________________________________________________________________
    ENDIF
!!fin para bandas

ENDDO




110 CONTINUE
write(*,*), 'termine HubbardMod'
write(lunDat,*), 'termine HubbardMod'

close(lunDat)
nAvgUPfinal=nAvgUP
nAVGDWfinal= nAVGDOWN

!-------------------------------------------------------------------------
!_________________________ FIN  __________________________________________
!-------------------------------------------------------------------------



CONTAINS

!===============================================================
!===============================================================
SUBROUTINE SystemEnergy(natoms,Efermi, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
                                      Eoarray,Uarray,Vbordarr,Ndizimac, Ue, hzeeman, imagEn, pasos,lunEnTotal)
Use Routines2
    implicit none
    INTEGER,                              INTENT(IN)     :: natoms,Ndizimac,pasos,lunEnTotal
    REAL(8), DIMENSION(natoms),         INTENT(IN)     :: nAvgUP,nAvgDOWN
    REAL(8),                              INTENT(IN)     :: hzeeman, imagEn
    REAL(8), intent(in), dimension(natoms) :: Eoarray, Uarray,Vbordarr
    complex(8), DIMENSION(natoms, natoms),INTENT(IN)     :: T00,T,TD
    complex(8), DIMENSION(natoms, natoms),INTENT(IN)     :: Ident
    REAL(8), DIMENSION(natoms),         INTENT(IN)     :: Ue
    real(8), intent(in) :: Efermi
    real(8) :: integ,Emax,Emin,dE,rho
!    integer, parameter :: pasos =1000
    real(8),dimension(5*pasos) :: Energypoints,Y
    integer :: j

        integ = 0.d0
        Emin = Efermi - 9.0d0
        Emax = Efermi
!        pasos= 3300
!        dE   = (Emax-Emin)/pasos

        Energypoints = Gauss5X(Emin, Emax, pasos)

        DO j=1, 5*pasos
!            En  = Emin + j*dE -Efermi
            En = Energypoints(j)
            rho = DOS(En, nAvgUP,nAvgDOWN,T00,T,TD,Ident,&
                                      Eoarray,Uarray,Vbordarr,natoms,Ndizimac, Ue, hzeeman, imagEn)
            if ( isnan(rho) ) then
                print*, '"DOS" is a NaN at En, j = ',En, j
                rho = 0.d0
            endif
            Y(j) = En*rho
            !integ = integ + dE*En*rho
        ENDDO
        integ =   Gauss5Quad(Y,Emin,Emax,pasos)
            open(unit=321,file='ddddd')
                do j=1,5*pasos
                    write(321,*) Energypoints(j),Y(j)
                enddo
            close(321)

        write(*,*) 'pasos, System Energy = ', pasos, integ
        write(lunEnTotal,*) 'pasos, System Energy = ', pasos, integ
        open(unit=321,file='totalEnergy')
            write(321,*) integ
        close(321)

RETURN
END SUBROUTINE SystemEnergy

!===============================================================
!===============================================================
REAL(8) FUNCTION DeltaE(natoms, ni,Vo, Efermi, imagEn,hzeeman, Eoarray,Uarray,Vbordarr, Ue,nAvgUP,nAVGDOWN)
implicit none
integer, intent(in)                   :: ni,natoms
real(8), intent(in),dimension(natoms) :: Eoarray,Uarray,Vbordarr, Ue,nAvgUP,nAVGDOWN
real(8), intent(in)                   :: Vo,Efermi,imagEn,hzeeman
complex(8), dimension(natoms,natoms)  :: Ident
complex(8), dimension(natoms,natoms)  :: Gup,Gdw,GRup,GRdw
COMPLEX(8), PARAMETER :: I  = dcmplx(0,1) !=sqrt(-1)
REAL(8),    PARAMETER :: Pi = 2*dACOS(0.0d0)
real(8) :: Emax,Emin,dE,En,deltaN, deltaUP_i,deltaDW_i
complex(8) :: suma,factUP,factDW
integer :: j,pasos

    suma         = 0.d0
    Ident        = 0.d0
    do j=1,natoms
        Ident(j,j) = 1.d0
    enddo

    deltaUP_i = Vo - Eoarray(ni)
    deltaDW_i = Vo - Eoarray(ni)

    Emin = Efermi - 10.0d0
    Emax = Efermi
    pasos= 300
    dE   = (Emax-Emin)/pasos
    DO j=1, pasos
        En = Emin + j*dE
        call GreenInicial(Gup,Gdw,natoms,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, En, imagEn)
        GRup = cGreenRnrmlzd(Gup, T00,T,TD, Ident, atomosT,Ndizimac)
        GRdw = cGreenRnrmlzd(Gdw, T00,T,TD, Ident, atomosT,Ndizimac)


        factUP = 1.0d0 - GRup(ni,ni)*deltaUP_i
        factDW = 1.0d0 - GRdw(ni,ni)*deltaDW_i
        suma   = suma + dE*Log(factUP*factDW)
!write(*,*) GRup(ni,ni),Log(factUP*factDW)
!stop
    ENDDO
    DeltaE = dimag(suma)/Pi
    write(*,*) 'deltaE = ', DeltaE

RETURN
END FUNCTION DeltaE
!===============================================================
!===============================================================



!_________________________________________________________________________
!_________________________________________________________________________

SUBROUTINE GreenInicial(Gup,Gdw,atomosT,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, Energy, imagEn)
    IMPLICIT NONE
    INTEGER,                                 INTENT(IN)    :: atomosT
    COMPLEX(8), DIMENSION(atomosT, atomosT), INTENT(INOUT) :: Gup,Gdw
    REAL(8),    dimension(atomosT),          INTENT(IN)    :: nAvgUP, nAvgDOWN, Ue
    REAL(8),                                 INTENT(IN)    :: hzeeman, Energy, imagEn
    REAL(8), dimension(atomosT), INTENT(IN)                :: Eoarray, Uarray,Vbordarr
    COMPLEX(8), PARAMETER                                  :: I  = dcmplx(0,1) !=sqrt(-1)
    integer :: s
    real(8):: nn
    Gup = 0.0d0
    Gdw = 0.0d0
do s=1,atomosT
    !nn = nAvgDOWN(s)*nAvgUP(s)
    Gup(s,s)= 1.0d0/( Energy +I*imagEn  -( Eoarray(s)+Vbordarr(s)+Ue(s)  +Uarray(s)*(nAvgDOWN(s) -0.5d0 ) -hzeeman    ) )
    Gdw(s,s)= 1.0d0/( Energy +I*imagEn  -( Eoarray(s)+Vbordarr(s)+Ue(s)  +Uarray(s)*(nAvgUP(s)   -0.5d0 ) +hzeeman    ) )





!    ! Se retiró Eoarray(s) porque el nivel de Fermi se duplicaba.
!    ! dica : si consideras una energia "Ei" aqui, ya no la consideres en la matriz T00, (tampoco Vbordarr!!) o
!    ! o viceversa (al menos para la autoconsistencia, en el programa de Leonor no influye).

    ! Repara que se colocó 0.0d0 en vez de 0.5 (venia de Wakabayashi) para dar con el nivel de Fermi correcto!!!
    ! o aproximadamente correcto ^^'
enddo
END SUBROUTINE GreenInicial


!_________________________________________________________________________
!_________________________________________________________________________
!FUNCTION GreenUP(atomosT, nAvgDOWN, hzeeman, U, Ue, Efermi, imagEn)
!    IMPLICIT NONE
!    INTEGER,                                INTENT(IN)  :: atomosT
!    REAL(8),    dimension(atomosT),         INTENT(IN)  :: nAvgDOWN, Ue
!    REAL(8),                                INTENT(IN)  :: hzeeman, U, Efermi, imagEn
!    COMPLEX(8), DIMENSION(atomosT, atomosT)             :: GreenUP, Gup
!    COMPLEX(8), PARAMETER                               :: I  = dcmplx(0,1) !=sqrt(-1)
!    integer :: s
!        Gup = 0.0d0
!        !!!$OMP PARALLEL DO
!        do s=1,atomosT
!            Gup(s,s)   = 1.0d0/( Efermi +I*imagEn  -( U*(nAvgDOWN(s) -0.5d0 ) -hzeeman + Ue(s)    ) )
!        enddo
!        !!!$OMP END PARALLEL DO
!        GreenUP = Gup
!        RETURN
!END FUNCTION GreenUP
!-----------------------------

!FUNCTION GreenDOWN(atomosT, nAvgUP, hzeeman, U, Ue, Efermi, imagEn)
!    IMPLICIT NONE
!    INTEGER,                                INTENT(IN)  :: atomosT
!    REAL(8),    dimension(atomosT),         INTENT(IN)  :: nAvgUP, Ue
!    REAL(8),                                INTENT(IN)  :: hzeeman, U, Efermi, imagEn
!    COMPLEX(8), DIMENSION(atomosT, atomosT)             :: GreenDOWN, Gdw
!    COMPLEX(8), PARAMETER                               :: I  = dcmplx(0,1) !=sqrt(-1)
!    integer :: s
!        Gdw = 0.0d0
!        !!!$OMP PARALLEL DO
!        do s=1,atomosT
!            Gdw(s,s)   = 1.0d0/( Efermi  +I*imagEn  -( U*(nAvgUP(s)  - 0.5d0 )    +hzeeman + Ue(s)    ) )
!        enddo
!        !!!$OMP END PARALLEL DO
!       GreenDOWN = Gdw
!        RETURN
!END FUNCTION GreenDOWN
!_________________________________________________________________________
!_________________________________________________________________________







SUBROUTINE IntegraPlanoCmplx(nAvgUP, nAvgDOWN, atomosT,hzeeman, Eoarray,Uarray,Vbordarr, Ue, Efermi, T00,T,TD,Ident, Ndizimac)
USE Routines1
! Integra desde -Inf hasta E=Efermi, pero usa el plano complejo para llevar a cabo la integraci?n.
! Devuelve nAvgUp y nAvgDown recalculados, teniendo como referencia ese Efermi.
! Para integrar, en vez de variar en la variable del eje x, usa el eje Y, que en este caso corresponde
! a la variable imagEn, esto es, imagEn es variable, y E=Efermi se mantiene constante.
!use omp_lib
    IMPLICIT NONE
    INTEGER,    INTENT(IN):: atomosT, Ndizimac

    REAL(8),    Dimension(atomosT),         INTENT(INOUT)   :: nAvgUP, nAvgDOWN
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)      :: T00,T,TD
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)      :: Ident
    REAL(8),    Dimension(atomosT),         INTENT(IN)      :: Ue
    REAL(8),                                INTENT(IN)      :: hzeeman, Efermi
    REAL(8),    dimension(atomosT),         intent(in)      :: Eoarray, Uarray,Vbordarr

    REAL(8),    PARAMETER :: Pi = 2*dACOS(0.0d0)
    REAL(8),    PARAMETER :: yi =0.0d0,  yf =1.0d0
!    INTEGER,    PARAMETER :: Ninterv=20, Ly =5*Ninterv !Ly=size(Y)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!!    INTEGER,    PARAMETER :: Ninterv=20, Ly =5*Ninterv !Ly=size(Y)
    INTEGER,    PARAMETER :: Ninterv=5, Ly =5*Ninterv !Ly=size(Y)  ! <<<< se tuvo que reducir aqui, para disminuir el tiempo :(
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    REAL(8),    DIMENSION(Ly, atomosT)      :: fup,fdown
    COMPLEX(8), DIMENSION(atomosT, atomosT) :: Gup, Gdown, GRup, GRdown
    REAL(8),    DIMENSION(Ly)               :: Y
    REAL(8)                                 :: yy, imagEn2
    INTEGER                                 :: k, s

!inicializar
!Gup=.0d0; Gdown=.0d0
!GRup=.0d0; GRdown=.0d0; Y=.0d0

        !calcular y
            Y  = Gauss5X(yi, yf, Ninterv)
        !inicializando

!            fup     = 0.0d0
!            fdown   = 0.0d0
!            Gup     = 0.0d0
!            Gdown   = 0.0d0
        !funciones de green
            do k=1, Ly
                yy      = Y(k)
                imagEn2 = yy/(1.0d0 - yy)


call GreenInicial(Gup,Gdown,atomosT,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, Efermi, imagEn2)

                !Gup   = GreenUP(atomosT, nAvgDOWN, hzeeman, U, Ue, Efermi, imagEn2)
                !Gdown = GreenDOWN(atomosT, nAvgUP, hzeeman, U, Ue, Efermi, imagEn2)

!!! cada "cGreenRnrmlzd" usa 27 cores, entonces paralelizar lo siguiente seria dividir
!!! los cores. El resultado es que demora mucho mas tiempo paralelizar lo siguiente.
!!! Mejor no paralelizar lo siguiente.
!!!***   !$OMP PARALLEL SECTIONS
!!!***   !$OMP SECTION
                GRup   = cGreenRnrmlzd(Gup,  T00,T,TD, Ident, atomosT, Ndizimac)
!!!***    !$OMP SECTION
                GRdown = cGreenRnrmlzd(Gdown,T00,T,TD, Ident, atomosT, Ndizimac)

!!!***   !$OMP END PARALLEL SECTIONS

                !definiendo el integrando
                    do s=1,atomosT
                        fup(k,s)   = Real(GRup(s,s))/(  (1.0d0 - yy)**2.0d0  )
                        fdown(k,s) = Real(GRdown(s,s))/(  (1.0d0 - yy)**2.0d0  )
                    enddo
            enddo

        !sumar o integraR
            do s = 1, atomosT
                nAvgUP(s)   = 0.5d0 + Gauss5Quad(fup(:,s), yi,yf, Ninterv)/Pi ! integra en cada columna
                nAvgDOWN(s) = 0.5d0 + Gauss5Quad(fdown(:,s), yi,yf, Ninterv)/Pi ! integra en cada columna
            enddo

END SUBROUTINE IntegraPlanoCmplx
!_________________________________________________________________________
!_________________________________________________________________________

SUBROUTINE ocupacao(nAvgUP,nAvgDOWN,Efermi, T00,T,TD,Ident,Eoarray,Uarray,Vbordarr,atomosT,Ndizimac, &
            Ue, hzeeman,imagEn,Ncharges)
    IMPLICIT NONE
    INTEGER,                                INTENT(IN)     :: atomosT,Ndizimac
    REAL(8),    DIMENSION(atomosT),         INTENT(INOUT)  :: nAvgUP,nAvgDOWN,Ue
    REAL(8),                                INTENT(INOUT)  :: Efermi
    REAL(8),                                INTENT(IN)     :: hzeeman,imagEn,Ncharges
    REAL(8),    DIMENSION(atomosT),         INTENT(IN)     :: Eoarray, Uarray,Vbordarr
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)     :: T00,T,TD
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)     :: Ident

    !REAL(8),    DIMENSION(atomosT)      :: UV,tempU ! se usaba con poisson subrotuine
    REAL(8),  allocatable, dimension(:) :: Vhex
    REAL(8),    PARAMETER               :: dNeminimo  =1.0d-3! dNeminimo  =1.0d-5 !!simplificar a: 1.0d-3
    REAL(8),    DIMENSION(atomosT)      :: nAvgUPi, nAvgDOWNi
    REAL(8)                             :: dNe, dNeTest, DOSf, maxdEfermi, ddE, ddq, Efermitemp,temp,valor
    INTEGER                             :: j,k

!    UV=0.0d0  ! se usaba con poisson subrotuine
    j   = 1
    k   = 1
    nAvgUPi   =  nAvgUP
    nAvgDOWNi =  nAvgDOWN
!print*, Ncharges
    !integralPlanoComplexo:
    dNe = dQ_exceso(atomosT,Ndizimac,nAvgUPi,nAvgDOWNi,Efermi,hzeeman,Ue,Eoarray,&
                    Uarray,Vbordarr,T00,T,TD,Ident,Ncharges)
    write(*,*)    'exceso del Numero de electrones en el Ribbon: dNe     =', dNe

    DO WHILE (  ((abs(dNe) > dNeminimo)  .or. (j==1)) .and.(j<3) )
        nAvgUPi   =  nAvgUP
        nAvgDOWNi =  nAvgDOWN
        !tempU = Ue-UV
        DOSf = DOS(Efermi,nAvgUPi,nAvgDOWNi,T00,T,TD,Ident,Eoarray,Uarray,Vbordarr,atomosT,Ndizimac,Ue,hzeeman,imagEn)

                ! REPARA que aqui se ponde Ue-UV para retornar al valor antiguo de Ue, antes de Poisson (o VPuntualBlindado)
        maxdEfermi = -(  dNe/atomosT  )/DOSf

        ddE     = maxdEfermi
        if ( Abs(maxdEfermi)>1.0d0 ) ddE=SIGN(1.0d0,maxdEfermi) ! esto para evitar que se dispare cuando DOSf sea cero

        dNeTest = 0.0d0
        if (dNe==0.d0) then
           valor = 1.d0
        else
           valor = abs((dNeTest-dNe)/dNe)
        endif
!        print*, valor
        DO WHILE ( ( valor.gt.1.0d-3  )  .and.  (abs(Efermi)<20.d0)  .and.  (abs(dNe)>1.d-5) .and. (k < 10) ) !***
        k = k + 1
!        DO WHILE ( ( abs((dNeTest-dNe)/dNe).gt.1.0d-3  )  .and.  (abs(Efermi)<20.d0)   ) !***
!        DO WHILE ( ( abs((dNeTest-dNe)/dNe).gt.1.0d-3  )  .and.  (abs(Efermi)<20.d0)   ) !***
                ! se integra alrededor del nivel de Fermi, en un intervalo de ancho ddE, y se calcula la carga ddq:
                nAvgUPi   =  nAvgUP
                nAvgDOWNi =  nAvgDOWN
                !tempU = Ue-UV
                ddq = dQ_Efermi2(atomosT,Ndizimac,  ddE,Efermi,hzeeman,imagEn, &
                        nAvgUPi,nAvgDOWNi,Ue,Eoarray, Uarray,Vbordarr,   T00,T,TD,Ident ) !Ue0=Ue-UV

                ! suma acumulada hasta llegar a dNe:
                dNeTest = dNeTest + ddq
                ! caso que lo acumulado sobrepase dNe:
                if (abs(dNeTest) > abs(dNe)) then
                        dNeTest = dNeTest - ddq !regresa al acumulo anterior
                        ddE     = ddE/2.0d0     !regresa al punto anterior, pero para integrar en un intervalo ddE/2
                else
                    temp = abs(dNeTest/dNe)
                    if (  (0<temp) .and. (temp<1.0d-4) )      ddE = 1.0d3*ddE
                    if (  (1.0d-4<temp) .and. (temp<1.0d-3) ) ddE = 1.0d2*ddE
                    if (  (1.0d-3<temp) .and. (temp<1.0d-2) ) ddE = 1.0d1*ddE
                    !if (  (1.0d-2<temp)                     ) ddE = ddE
                    if ( Abs(ddE)>1.0d0 ) ddE=SIGN(1.0d0,ddE) ! esto para evitar que se dispare
                    Efermi  = Efermi  + ddE ! caso contrario, el nivel de Fermi se desplaza
                endif
                write(*,900)      'dNeTest, dNe, Efermi, ddE:     ', dNeTest, dNe, Efermi, ddE
                if (ddq==0.0d0) dNeTest=dNe
        ENDDO
            !----------------- fin del desplaziento del NIVEL DE FERMI ---------------------
            write(*,*),      'final:dNeTest,dNe,Efermi,ddE:  ', dNeTest, dNe, Efermi, ddE
            write(*,*),      ''
            write(lunDat,*), 'final:dNeTest,dNe,Efermi,ddE:  ', dNeTest, dNe, Efermi, ddE
            write(lunDat,*), ''


        !-------------------------------------------------------------------------------------------
        ! cuenta el numero de veces que se ha tenido que cambiar el nivel de Fermi para hacer que dQ <= Qminimo.
        j=j+1

!==================================================================
!===== OSCILACIONES ===============================================
!
!*        if (j>50) go to 150   ! para retirarse en caso que se este oscilando ... :
!
!          if (  (hu<50).and.(j>1)  ) then
!                goto 150
!          elseif (  (hu>=50).and.(j>2)  )  then
!                goto 150
!          endif
!
!==================================================================
        nAvgUPi   =  nAvgUP
        nAvgDOWNi =  nAvgDOWN
        dNe = dQ_exceso(atomosT,Ndizimac,nAvgUPi,nAvgDOWNi,Efermi,hzeeman,Ue,Eoarray,&
                        Uarray,Vbordarr,T00,T,TD,Ident,Ncharges)
        write(*,*)    'exceso del Numero de electrones en el Ribbon: dNe     =', dNe
   ENDDO

150 CONTINUE
    !salida:
    nAvgUP   = nAvgUPi
    nAvgDOWN = nAvgDOWNi


900 format(a31, 4f12.6)
END SUBROUTINE ocupacao


!==================================================================
REAL FUNCTION dQ_exceso(atomosT,Ndizimac,nAvgUPi,nAvgDOWNi,Efermi,hzeeman,Ue,Eoarray, &
                        Uarray,Vbordarr,T00,T,TD,Ident,Ncharges)
    implicit none
    INTEGER,                                INTENT(IN)     :: atomosT,Ndizimac
    REAL(8),    DIMENSION(atomosT),         INTENT(INOUT)  :: nAvgUPi,nAvgDOWNi
    REAL(8),                                INTENT(IN)     :: Efermi,hzeeman,Ncharges
    REAL(8),    DIMENSION(atomosT),         INTENT(in)     :: Ue,Eoarray, Uarray,Vbordarr
    complex(8),    DIMENSION(atomosT, atomosT),INTENT(IN)  :: T00,T,TD,Ident

            !_______________ integrando en el plano complejo __________________________________________________________
        write(*,*)      '**nuevo nivel de Fermi:', Efermi, '**'

        call IntegraPlanoCmplx(nAvgUPi,nAvgDOWNi,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi,T00,T,TD,Ident,Ndizimac)

        ! se integra desde -Infinito hasta Efermi y devuelve nAvgUPi, nAvgDOWNi RECALCULADOS!!!.
        ! Para el calculo se usa funciones de Green, en donde entran las energias debidos al pot.electrico incluidas en Ue
        !
        ! Para estos nuevos valores de nAvgUpi y nAvgDOWNi hay un nuevo valor para dNe:

!print*, Ncharges
!dQ_exceso =  Sum(  nAvgUPi + nAvgDOWNi   )  -   ( atomosT +Bcharge )! dNe es el ***exceso del NUMERO DE ELECTRONES*** en el ribbon todo.
dQ_exceso =  Sum(  nAvgUPi + nAvgDOWNi   )  -   Ncharges ! dNe es el ***exceso del NUMERO DE ELECTRONES*** en el ribbon todo.
!********************************************************************


return
END FUNCTION  dQ_exceso


!==================================================================
REAL(8) FUNCTION dQ_Efermi(atomosT,Ndizimac,  ddE,Efermi,hzeeman,imagEn, &
                        nAvgUPtemp,nAvgDOWNtemp,Ue0,Eoarray, Uarray,Vbordarr,   T00,T,TD,Ident ) !Ue0=Ue-UV
    implicit none
    INTEGER,                                INTENT(IN)     :: atomosT,Ndizimac
    REAL(8),                                INTENT(IN)     :: ddE, Efermi,hzeeman,imagEn
    REAL(8),    DIMENSION(atomosT),         INTENT(IN)     :: nAvgUPtemp, nAvgDOWNtemp,Ue0,Eoarray, Uarray,Vbordarr
    complex(8),    DIMENSION(atomosT, atomosT),INTENT(IN)     :: T00,T,TD,Ident

    INTEGER,    PARAMETER               :: Ninterv2 = 5,    Lx2 = 5*Ninterv2
    REAL(8),    DIMENSION(Lx2)          :: puntosE, puntosDOS
    integer :: h

                puntosE=.0d0; puntosDOS=.0d0;
                if (ddE>0) then
                        puntosE = Gauss5X(    Efermi  , Efermi +ddE, Ninterv2 )
                else
                        puntosE = Gauss5X( Efermi +ddE,    Efermi  , Ninterv2 )
                endif

                !puntos en el eje Y:
                do h=1,Lx2
                    puntosDOS(h) = DOS(puntosE(h), nAvgUPtemp,nAvgDOWNtemp,T00,T,TD,Ident, &
                                       Eoarray,Uarray,Vbordarr,atomosT,Ndizimac, Ue0, hzeeman, imagEn)
                    ! REPARA que aqui se ponde Ue-UV para retornar al valor antiguo de Ue, antes de Poisson (o VPuntualBlindado)
 !                   write(*,*) puntosE(h), puntosDOS(h)
                enddo

                puntosDOS = puntosDOS*atomosT! correccion, debido a que la integral de DOS en todo el nanoribbon es 1
                ! es 1 segun como lo definiste en tu funcion DOS ^^'

                ! INTEGRA, USANDO LOS PUNTOS X e Y:
                if (ddE>0) then
                        dQ_Efermi    = -Gauss5Quad(puntosDOS, Efermi, Efermi +ddE, Ninterv2)
                else
                        dQ_Efermi    =  Gauss5Quad(puntosDOS, Efermi +ddE, Efermi, Ninterv2)
                endif

return
END FUNCTION dQ_Efermi

REAL(8) FUNCTION dQ_Efermi2(atomosT,Ndizimac,  ddE,Efermi,hzeeman,imagEn, &
                        nAvgUP,nAvgDOWN,Ue0,Eoarray, Uarray,Vbordarr,   T00,T,TD,Ident ) !Ue0=Ue-UV
    implicit none
    INTEGER,                                INTENT(IN)     :: atomosT,Ndizimac
    REAL(8),                                INTENT(IN)     :: ddE, Efermi,hzeeman,imagEn
    REAL(8),    DIMENSION(atomosT),         INTENT(IN)     :: nAvgUP, nAvgDOWN,Ue0,Eoarray, Uarray,Vbordarr
    complex(8),    DIMENSION(atomosT, atomosT),INTENT(IN)     :: T00,T,TD,Ident
REAL(8),    DIMENSION(atomosT) ::  nAvgUPtemp,nAvgDOWNtemp
    INTEGER,    PARAMETER               :: Ninterv2 = 5,    Lx2 = 5*Ninterv2
    REAL(8),    DIMENSION(Lx2)          :: puntosE, puntosDOS
    integer :: h
    REAL(8) :: Q1,Q2

        nAvgUPtemp   = nAvgUP
        nAvgDOWNtemp = nAvgDOWN
if (ddE>0) then
    call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi,T00,T,TD,Ident,Ndizimac)
    Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  )
        nAvgUPtemp   = nAvgUP
        nAvgDOWNtemp = nAvgDOWN
    call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+ddE,T00,T,TD,Ident,Ndizimac)
    Q2 = Sum(  nAvgUPtemp + nAvgDOWNtemp  )
    dQ_Efermi2 = -(Q2 - Q1)

    if (dQ_Efermi2>0.0d0) then
        nAvgUPtemp   = nAvgUP
        nAvgDOWNtemp = nAvgDOWN
        dQ_Efermi2 = dQ_Efermi(atomosT,Ndizimac,ddE,Efermi,hzeeman,imagEn,nAvgUPtemp,nAvgDOWNtemp,&
                                Ue0,Eoarray,Uarray,Vbordarr,T00,T,TD,Ident) !Ue0=Ue-UV
        print*, dQ_Efermi2,' = dQ_Efermi2'
        if (dQ_Efermi2>0.0d0) then
            print*,'dQ_Efermi2 es tambien positivo!'
            stop
        endif
    endif


else
    call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+ddE,T00,T,TD,Ident,Ndizimac)
    Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  )
        nAvgUPtemp   = nAvgUP
        nAvgDOWNtemp = nAvgDOWN
    call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi,T00,T,TD,Ident,Ndizimac)
    Q2 = Sum(  nAvgUPtemp + nAvgDOWNtemp  )
    dQ_Efermi2 = (Q2 - Q1)

    if (dQ_Efermi2<0.0d0) then
        nAvgUPtemp   = nAvgUP
        nAvgDOWNtemp = nAvgDOWN
        dQ_Efermi2 = dQ_Efermi(atomosT,Ndizimac,ddE,Efermi,hzeeman,imagEn,nAvgUPtemp,nAvgDOWNtemp,&
                                Ue0,Eoarray,Uarray,Vbordarr,T00,T,TD,Ident) !Ue0=Ue-UV
        print*, dQ_Efermi2,' = dQ_Efermi2'
        if (dQ_Efermi2<0.0d0) then
            print*,'dQ_Efermi2 es tambien Negativo!'
            if (abs(dQ_Efermi2)<1.0d-3) then
                dQ_Efermi2 =0.0d0
            else
                stop
            endif
        endif
    endif
    if (dQ_Efermi2<0.0d0) then
        print*, 'error dQ_Efermi2 negativo, incrementar N en IntegralCpmplx?'



call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+6*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'6*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+5*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'5*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+4*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'4*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+3*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'3*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+2*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'2*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi+1*ddE,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'1*ddE'
    nAvgUPtemp   = nAvgUP
    nAvgDOWNtemp = nAvgDOWN
call IntegraPlanoCmplx(nAvgUPtemp,nAvgDOWNtemp,atomosT,hzeeman,Eoarray,Uarray,Vbordarr,Ue,Efermi,T00,T,TD,Ident,Ndizimac)
Q1 = Sum(  nAvgUPtemp + nAvgDOWNtemp  ) ; print*, Q1,'0*ddE'

Q1=dQ_Efermi(atomosT,Ndizimac,ddE,Efermi,hzeeman,imagEn,nAvgUPtemp,nAvgDOWNtemp,Ue0,Eoarray,Uarray,Vbordarr,T00,T,TD,Ident) !Ue0=Ue-UV
dQ_Efermi2 = Q1
print*, Q1,' = dQ_Efermi2'
if (dQ_Efermi2<0.0d0) then
    print*,'dQ_Efermi2 es tambien negativo'
    stop
endif

    endif
endif
return
END FUNCTION dQ_Efermi2





!_________________________________________________________________________
!_________________________________________________________________________

REAL(8) FUNCTION DOS(En, nAvgUP,nAvgDOWN,T00,T,TD,Ident,Eoarray,Uarray,Vbordarr,atomosT,Ndizimac, Ue, hzeeman, imagEn)
    IMPLICIT NONE
    INTEGER,                              INTENT(IN)     :: atomosT,Ndizimac
    REAL(8), DIMENSION(atomosT),         INTENT(IN)     :: nAvgUP,nAvgDOWN
    REAL(8),                              INTENT(IN)     :: hzeeman,  En, imagEn
    REAL(8), intent(in), dimension(atomosT) :: Eoarray, Uarray,Vbordarr
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)     :: T00,T,TD
    complex(8), DIMENSION(atomosT, atomosT),INTENT(IN)     :: Ident
    REAL(8), DIMENSION(atomosT),         INTENT(IN)     :: Ue

    REAL(8),    PARAMETER :: Pi = 2*dACOS(0.0d0)
    COMPLEX(8), DIMENSION(atomosT, atomosT) :: Gup, Gdown, GRup, GRdown
    !INTEGER,    PARAMETER                   :: Ninterv=15, Lx=5*Ninterv !Lx=size(X)
    REAL(8)                                 :: dens, NewEfermi, dup, ddw
    INTEGER                                 :: s

    !inicializar
    GRup=.0d0; GRdown=.0d0;

            NewEfermi = En


call GreenInicial(Gup,Gdown,atomosT,nAvgUP,nAvgDOWN, hzeeman, Eoarray,Uarray,Vbordarr, Ue, NewEfermi, imagEn)
            !Gup   = GreenUP(atomosT, nAvgDOWN, hzeeman, U, Ue, NewEfermi, imagEn)
            !Gdown = GreenDOWN(atomosT, nAvgUP, hzeeman, U, Ue, NewEfermi, imagEn)

            GRup   = cGreenRnrmlzd(Gup,T00,T,TD, Ident, atomosT,Ndizimac)
            GRdown = cGreenRnrmlzd(Gdown,T00,T,TD, Ident, atomosT,Ndizimac)

            dens = 0.0d0
            do s=1,atomosT
                dup = -Dimag(    GRup(s,s)  )
                ddw = -Dimag(    GRdown(s,s))
                if ( (dup<0.0d0).or.(ddw<0.0d0) ) then
                    WRITE(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
                    WRITE(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

                    WRITE(*,*),'---deu_DENSIDADE_NEGATIVA!---En,dup,dw = ', En, dup, ddw
                    dup=0.0d0
                    ddw=0.0d0
                    WRITE(*,*) '**********************************************'
                    WRITE(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
                    WRITE(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

                    !stop
                endif
                dens = dens + (dup + ddw)/(atomosT*Pi)
                !dens = dens -Dimag(    GRup(s,s)  + GRdown(s,s)  )/(atomosT*Pi)
            enddo
            DOS = dens
            return
END FUNCTION DOS





!_________________________________________________________________________
!_________________________________________________________________________
! FUNCTION cGreenRnrmlzd(G,T00,T,TD, Ident, atomosT,Ndizimac)
! 'c' stands for cluster
FUNCTION cGreenRnrmlzd(G,T00,T,TD, Ident, n, Nd) ! 'c' stands for cluster
    use  Routines1
    use  Routines2
    implicit none
    integer,                    INTENT(IN)  :: n, Nd
    complex(8), DIMENSION(n,n), INTENT(IN)  :: T00,T,TD
    complex(8), DIMENSION(n,n), INTENT(IN)  :: Ident
    COMPLEX(8), DIMENSION(n,n), INTENT(IN)  :: G

    COMPLEX(8), DIMENSION(n,n)  :: cGreenRnrmlzd,GR,GR0, TR, TRD, Z, ZzD, &
                                   TRs, TRDs,InvGR, &
                                   tempM, tempMa, tempMb, tempMc, tempMd
    integer :: j, continuaGR,k,ii

!INTEGER, SAVE :: CONTA = 0
!CONTA = CONTA +1
!PRINT*, CONTA

        continuaGR = 1

        !tempM = cProd(n,G,T00)
                tempM = cProd(n,G,T00)

        !GR  = cRenormQ(Ident, tempM, G, n) !8888888888
                GR = G
                call cRenormQroutine(Ident, tempM, GR, n) ! << G es actualizado (overwritten)

        TR  = T
        TRD = TD
        GR0 = 0.0d0


    IF (Nd==-1) THEN !dizimacoes livres
        j=1
        do while (continuaGR==1)  !j=1,Ndizimac

!!!!!!****  !$OMP PARALLEL SECTIONS
!!!!!!****  !$OMP SECTION
            Z      = cProd(n,GR,TR)    ! Z(N-1)   = GR(N-1)*TR(N-1)
!!!!!!****  !$OMP SECTION
            ZzD    = cProd(n,GR,TRD) ! ZzD(N-1) = GR(N-1)*TRD(N-1)
!!!!!!****  !$OMP SECTION
            tempMa = cProd(n,TR,GR)
            !TR    = cProd(n,tempMa, TR) !   TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
!!!!!!****  !$OMP SECTION
            tempMb = cProd(n,TRD,GR)
            !TRD   = cProd(n,tempMb, TRD) ! TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
!!!!!!****  !$OMP END PARALLEL SECTIONS


!!!!!!****  !$OMP PARALLEL SECTIONS
!!!!!!****  !$OMP SECTION
            TR     = cProd(n,tempMa, TR) !   TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
!!!!!!****  !$OMP SECTION
            TRD    = cProd(n,tempMb, TRD) ! TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
!!!!!!****  !$OMP SECTION
            tempMc = cProd(n,Z,ZzD)
!!!!!!****  !$OMP SECTION
            tempMd = cProd(n,ZzD,Z)
!!!!!!****  !$OMP END PARALLEL SECTIONS


            !tempM = cProd(n,Z,ZzD) + cProd(n,ZzD,Z)
            tempM = tempMc + tempMd

            !GR = cRenormQ(Ident, tempM , GR, n)
                  call cRenormQroutine(Ident, tempM, GR, n) ! << GR es actualizado (overwritten)


            !GR = cRenormQ(Ident, cProd(n,Z,ZzD) + cProd(n,ZzD,Z), GR, n)
            !!! GR=inv(  Ident - Z*ZzD - ZzD*Z  )*GR;  %GR(N)= inv(...)*GR(N-1)
            if (j>1) continuaGR = continuaDizimando(GR, GR0, n)
            GR0 = GR
            j=j+1
        enddo

    ELSE
        !!$OMP PARALLEL DO
        do j=1,Nd
            Z  = matmul(GR,TR)    ! Z(N-1)   = GR(N-1)*TR(N-1)
            ZzD= matmul(GR,TRD) ! ZzD(N-1) = GR(N-1)*TRD(N-1)

            tempM = matmul(TR,GR)
            TR = matmul(tempM, TR) !   TR(N) = TR(N-1)*GR(N-1)*TR(N-1)

            tempM = matmul(TRD,GR)
            TRD= matmul(tempM, TRD) ! TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)

            tempM = matmul(Z,ZzD) + matmul(ZzD,Z)
            GR = RenormQ(Ident, tempM , GR, n)
            !!! GR=inv(  Ident - Z*ZzD - ZzD*Z  )*GR;  %GR(N)= inv(...)*GR(N-1)
        enddo
       !!$OMP END PARALLEL DO
    ENDIF
    cGreenRnrmlzd = GR!

    RETURN

END FUNCTION cGreenRnrmlzd

!_________________________________________________________________________
!_________________________________________________________________________





!_________________________________________________________________________
!_________________________________________________________________________
!_________________________ dizimacao _____________________________________
!-------------------------------------------------------------------------
FUNCTION GreenRnrmlzd(G,T00,T,TD, Ident, atomosT,Ndizimac)
    implicit none
    integer,                                 INTENT(IN)  :: atomosT,Ndizimac
    complex(8), DIMENSION(atomosT, atomosT), INTENT(IN)  :: T00,T,TD
    complex(8),    DIMENSION(atomosT, atomosT), INTENT(IN)  :: Ident
    COMPLEX(8), DIMENSION(atomosT, atomosT), INTENT(IN)  :: G

    COMPLEX(8), DIMENSION(atomosT, atomosT)  :: GreenRnrmlzd,GR,GR0, TR, TRD, Z, ZzD,tempMa,tempMb,tempMc
    integer :: j, continuaGR

        continuaGR = 1
        tempMa = matmul(G,T00)
        GR  = RenormQ(Ident, tempMa, G, atomosT)
        TR  = T
        TRD = TD
        GR0 = 0.0d0

    IF (Ndizimac==-1) THEN !dizimacoes livres
        j=1
        do while (continuaGR==1)  !j=1,Ndizimac
            Z  = matmul(GR,TR)    ! Z(N-1)   = GR(N-1)*TR(N-1)
            ZzD= matmul(GR,TRD) ! ZzD(N-1) = GR(N-1)*TRD(N-1)
            TR = matmul(matmul(TR,GR), TR) !   TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
            TRD= matmul(matmul(TRD,GR), TRD) ! TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
            tempMa  =  matmul(Z,ZzD)
            tempMb  =  matmul(ZzD,Z)
            tempMc  = tempMa + tempMb
            GR = RenormQ(Ident, tempMc, GR, atomosT)
            !!! GR=inv(  Ident - Z*ZzD - ZzD*Z  )*GR;  %GR(N)= inv(...)*GR(N-1)
            if (j>1) continuaGR = continuaDizimando(GR, GR0, atomosT)
            GR0 = GR
            j=j+1
        enddo
    ELSE
        !!$OMP PARALLEL DO
        do j=1,Ndizimac
            Z  = matmul(GR,TR)    ! Z(N-1)   = GR(N-1)*TR(N-1)
            ZzD= matmul(GR,TRD) ! ZzD(N-1) = GR(N-1)*TRD(N-1)
            TR = matmul(matmul(TR,GR), TR) !   TR(N) = TR(N-1)*GR(N-1)*TR(N-1)
            TRD= matmul(matmul(TRD,GR), TRD) ! TRD(N) = TRD(N-1)*GR(N-1)*TRD(N-1)
            GR = RenormQ(Ident, matmul(Z,ZzD) + matmul(ZzD,Z), GR, atomosT)
            !!! GR=inv(  Ident - Z*ZzD - ZzD*Z  )*GR;  %GR(N)= inv(...)*GR(N-1)
        enddo
        !!$OMP END PARALLEL DO
    ENDIF
    GreenRnrmlzd = GR
    RETURN

END FUNCTION GreenRnrmlzd

!_________________________________________________________________________
!_________________________________________________________________________


!_________________________________________________________________________
!_________________________________________________________________________
END SUBROUTINE HubbardGreen
!_________________________________________________________________________
END Module HubbardMod
