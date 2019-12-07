MODULE neighbors
contains

FUNCTION getDistance2D( xi, yi, xj, yj ) RESULT (distance2D)
   IMPLICIT NONE
   REAL(8), INTENT(IN) :: xi, yi, xj, yj
   REAL(8) :: distance2D
   distance2D = DSQRT( ((xi - xj)**2) + ((yi - yj)**2) )
END FUNCTION getDistance2D


FUNCTION isFirstNeighbor(distance2D) RESULT (ans)
   IMPLICIT NONE
   REAL(8), INTENT(IN) :: distance2D
   LOGICAL :: ans
   REAL(8), parameter:: paramRed = 2.4595121d0              , &      !unidades Amstrong=10^-10 metro
                        dCC      = paramRed/sqrt(3.0d0)

   if ( ( ( 0.9 * dCC ) < distance2D ) .AND. ( distance2D < (1.1 * dCC) ) ) then
      ans = .TRUE.
   else 
      ans = .FALSE.
   endif
END FUNCTION isFirstNeighbor


SUBROUTINE writeNeighbors2D(nAtom, X, Y, Xs, Ys, fileOut)
   INTEGER, INTENT(IN) :: nAtom
   CHARACTER(len=*), INTENT(IN) :: fileOut
   REAL(8), DIMENSION(nAtom), INTENT(IN) :: X, Y, Xs, Ys
   INTEGER :: i, j
   INTEGER, PARAMETER :: lunOut = 100
   REAL(8) :: distance2D
   
   open(lunOut, file=fileOut)
   do i = 1, nAtom
      do j = 1, nAtom
         if (i .ne. j) then
            distance2D = getDistance2D( X(i), Y(i), Xs(j), Ys(j) )
            if ( isFirstNeighbor( distance2D ) ) then
               WRITE(lunOut, *) i, j, X(i), Y(i), Xs(j), Ys(j)
            endif
         endif
      enddo
   enddo
   close(lunOut)

END SUBROUTINE writeNeighbors2D


SUBROUTINE getUnrepeated(nAtom, fileName)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nAtom
   CHARACTER(len=*), INTENT(IN) :: fileName
   INTEGER :: i, j
   INTEGER, PARAMETER :: lun = 100, lunOutUnrepeat = 101
   REAL(8) :: xi, yi, xsj, ysj
   INTEGER, DIMENSION(nAtom, nAtom) :: indxUnrepeated
   
   indxUnrepeated = -1
   open(lun, file=fileName)
   open(lunOutUnrepeat, file="unrepeated.txt")
   do while (.TRUE.)
      read(lun, *, end=999)  i, j, xi, yi, xsj, ysj
      indxUnrepeated(i, j) = 1
      if ( indxUnrepeated(j, i) .NE. 1 ) then ! if it was not defined...
         WRITE(lunOutUnrepeat, *) i, j, xi, yi, xsj, ysj
      endif
   enddo
   999 CONTINUE
   close(lun)
   close(lunOutUnrepeat)
END SUBROUTINE getUnrepeated


SUBROUTINE shiftPositions(nAtom, X, Y, Xs, Ys, distX, distY)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nAtom
   REAL(8), INTENT(IN) :: distX, distY
   REAL(8), DIMENSION(nAtom), INTENT(INOUT) :: X, Y, Xs, Ys
   integer :: i
   do i = 1, nAtom
      Xs(i) = X(i) + distX
      Ys(i) = Y(i) + distY
   enddo
END SUBROUTINE shiftPositions


SUBROUTINE getNeighbors(nAtom, fileAtoms, normCh, normT)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nAtom
   CHARACTER(len=*), INTENT(IN) :: fileAtoms
   REAL(8), INTENT(IN) :: normCh, normT
   REAL(8) :: lenCellX, lenCellY
   INTEGER, PARAMETER :: lun = 101
   INTEGER :: i
   REAL(8), DIMENSION(nAtom) :: X, Y, Xs, Ys
   REAL(8), DIMENSION(2) :: R
   REAL(8) :: z, distX, distY
   REAL(8), parameter:: paramRed = 2.4595121d0 
   
   open(lun, file=fileAtoms)
      do i = 1, nAtom
         read(lun, *) X(i), Y(i), z
      enddo
   close(lun)
   
   lenCellX = paramRed * normCh
   lenCellY = paramRed * normT

   Xs = X
   Ys = Y
   CALL writeNeighbors2D(nAtom, X, Y, Xs, Ys, "neighbors00.txt")
   CALL getUnrepeated(nAtom, "neighbors00.txt")

   R(1) = 0.d0
   R(2) = 1.d0
   distX = lenCellX * R(1)
   distY = lenCellY * R(2)
   CALL shiftPositions(nAtom, X, Y, Xs, Ys, distX, distY)
   CALL writeNeighbors2D(nAtom, X, Y, Xs, Ys, "neighbors01.txt")

   R(1) = 0.d0
   R(2) = -1.d0
   distX = lenCellX * R(1)
   distY = lenCellY * R(2)
   CALL shiftPositions(nAtom, X, Y, Xs, Ys, distX, distY)
   CALL writeNeighbors2D(nAtom, X, Y, Xs, Ys, "neighbors10.txt")


   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   ! for Daiara's test:
   !
   R(1) = 1.d0
   R(2) = 0.d0
   distX = lenCellX * R(1)
   distY = lenCellY * R(2)
   CALL shiftPositions(nAtom, X, Y, Xs, Ys, distX, distY)
   CALL writeNeighbors2D(nAtom, X, Y, Xs, Ys, "neighbors01x.txt")

   R(1) =-1.d0
   R(2) = 0.d0
   distX = lenCellX * R(1)
   distY = lenCellY * R(2)
   CALL shiftPositions(nAtom, X, Y, Xs, Ys, distX, distY)
   CALL writeNeighbors2D(nAtom, X, Y, Xs, Ys, "neighbors10x.txt")
   !
   !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

END SUBROUTINE getNeighbors


SUBROUTINE bondings2D()
   IMPLICIT NONE
   integer :: i, j
   REAL(8) :: xA, yA, xB, yB
   open(50, file="TeX2-atoms-bondings")

   open(100, file="unrepeated.txt")
   do while (.TRUE.)
      read(100, *, end=9997)  i, j, xA, yA, xB, yB
      write(50, 10000) '\draw [thin] (axis cs: ', xA,', ', yA,') -- (axis cs: ',xB,', ', yB,') node[]{};'
   enddo
   9997 CONTINUE
   close(100)

   open(100, file="neighbors01.txt")
   do while (.TRUE.)
      read(100, *, end=9998)  i, j, xA, yA, xB, yB
      write(50, 10000) '\draw [thin] (axis cs: ', xA,', ', yA,') -- (axis cs: ',xB,', ', yB,') node[]{};'
   enddo
   9998 CONTINUE
   close(100)

   open(100, file="neighbors10.txt")
   do while (.TRUE.)
      read(100, *, end=9999)  i, j, xA, yA, xB, yB
      write(50, 10000) '\draw [thin] (axis cs: ', xA,', ', yA,') -- (axis cs: ',xB,', ', yB,') node[]{};'
   enddo
   9999 CONTINUE
   close(100)


   10000 format(a, f12.3, a, f12.3, a, f12.3, a, f12.3, a)
   close(100)
   close(50)
END SUBROUTINE bondings2D



END MODULE neighbors
