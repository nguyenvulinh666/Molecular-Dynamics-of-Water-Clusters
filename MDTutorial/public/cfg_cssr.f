      program cfg_cssr
c
c*********************************************************************
c
c     dl_poly routine for converting CONFIG files to CSSR files
c     suitable for input to CERIUS 2 by MSI.
c
c     origin of coordinates shifts from cell center to corner
c     to agree with CERIUS 2 format. But no shift for system with
c     imcon of 0, 4 and 5.
c
c     rotate the system (cell and atoms) to make the final cell
c     matrix have the form
c
c                       | A11 A12 A13 |
c                       |  0  A22 A23 |
c                       |  0   0  A33 |
c
c     which is used by CERIUS 2 of MSI to orientate unit cell.
c
c     this version developed from the original con_cssr.f
c     by d. collins at daresbury laboratory feb. 1996
c
c     author - xianglong yuan   alfred university
c     january 1998
c
c*********************************************************************
c

      implicit real*8 (a-h,o-z)

      parameter (mxatms=10000)

      logical find

      character*80 cfgname
      character*8 atmnam(mxatms)

      dimension cell(9)
      dimension celorg(9),cvt(9)
      dimension xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension xorg(mxatms),yorg(mxatms),zorg(mxatms)

c
c open CONFIG and config.cssr file

      nconf=15    ! channel 15 for CONFIG      file
      ncssr=16    ! channel 16 for config.cssr file

      inquire(file='CONFIG',exist=find)
      if(.not.find) then
        write(*,'(a)') 'CONFIG file not found'
        stop
      endif
      open(nconf,file='CONFIG',status='old')
      open(ncssr,file='config.cssr',status='unknown')

c
c read cfgname, levcfg, imcon

      read(nconf,'(a)') cfgname
      write(*,'(a)') cfgname
      read(nconf,'(2i10)') levcfg, imcon

c
c check levcfg

      if(levcfg.eq.0) then
        write(*,'(/,a)')'File contains coordinates only'
      elseif(levcfg.eq.1) then
        write(*,'(/,a)')'File contains coordinates & velocities'
      elseif(levcfg.eq.2) then
        write(*,'(/,a)')'File contains coordinates, velocities & forces'
      else
        write(*,"(/,'Error - invalid levcfg',i10)") levcfg
        stop
      endif

c
c check imcon and calculate cell matrix

      if(imcon.eq.0) then

        write(*,'(a)') 'No Periodic Boundary Condition'
        write(ncssr,'(/)')

      elseif(imcon.ge.1.and.imcon.le.6) then

        do i=1,3
          read(nconf,*) (celorg(i*3-3+j),j=1,3) 
        end do                  ! matrix element: cell(i,j) = cell(i*3-3+j)

        call cellrotate(celorg,cvt)
                            ! calculate conversion matrix to rotate system
        do i=1,3
        do j=1,3
          cell(i*3-3+j)=0.0d0
          do k=1,3
            cell(i*3-3+j)=cell(i*3-3+j)+celorg(i*3-3+k)*cvt(k*3-3+j)
          enddo
        enddo
        enddo               ! calculate new matrix of rotated cell

        call cellrite2cssr(cell,ncssr)
                            ! calc. and write new cell vectors to CSSR file
      else

        write(*,'(/,a)') 'Error - invalid imcon'
        stop

      endif

c
c read atom coordinates

      write(*,'(/,a)') 'Reading coordinates ...'

      do i = 1,10000
         read(nconf,'(a8)',end=100) atmnam(i)
         read(nconf,*) xorg(i),yorg(i),zorg(i) 
         if(levcfg.gt.0) read(nconf,*)
         if(levcfg.gt.1) read(nconf,*)
      end do

      write(*,'(a)') 'Max No. of atoms exceeded - program terminating'
      stop

100   write(*,'(a)') 'Coordinates read correctly'
      natms = i -1

c
c calculate new coordinates of rotated atoms

      if(imcon.eq.1.or.imcon.eq.2.or.imcon.eq.3.or.imcon.eq.6) then

        do i=1,natms
          xxx(i)=xorg(i)*cvt(1)+yorg(i)*cvt(4)+zorg(i)*cvt(7)
          yyy(i)=xorg(i)*cvt(2)+yorg(i)*cvt(5)+zorg(i)*cvt(8)
          zzz(i)=xorg(i)*cvt(3)+yorg(i)*cvt(6)+zorg(i)*cvt(9)
        enddo
                           ! calculate new coordinates of rotated atoms

        xx0=0.5*cell(1)+0.5*cell(4)+0.5*cell(7)
        yy0=0.5*cell(2)+0.5*cell(5)+0.5*cell(8)
        zz0=0.5*cell(3)+0.5*cell(6)+0.5*cell(9)
                           ! calculate vectors to shift origin from
                           ! corner (-0.5,-0.5,-0.5) to center (0,0,0)
      elseif(imcon.eq.0.or.imcon.eq.4.or.imcon.eq.5) then

        xx0 = 0.0d0
        yy0 = 0.0d0
        zz0 = 0.0d0
                    ! no shift for no periodic boundary, truncated
                    ! octahedron boundary and rhombic dodecahedron
      endif         ! boundary conditions
 
c
c write remaining part of CSSR file

      ncood = 1      ! coordinate system: Cartesian
      iz = 0         ! connectivity
      chge = 0.0     ! charge

      write(*,"('No. of atoms',25x,i10)") natms
      write(ncssr,"(2i4,1x,'Created from CONFIG file')")natms,ncood
      write(ncssr,'(a80)') cfgname

      do i=1,natms
        write(ncssr,'(i4,1x,a2,4x,3(f9.5,1x),8i4,1x,f7.3)')
     +    i,atmnam(i)(1:2),xxx(i)+xx0,yyy(i)+yy0,zzz(i)+zz0,
     +    iz,iz,iz,iz,iz,iz,iz,iz,chge
      enddo

      write(*,'(a)') 'good end'

      stop
      end


      subroutine cellrite2cssr(cell,ncssr)
c
c*********************************************************************
c
c     calculate the cell parameters from the cell matrix and write
c     to CSSR file of CERIUS 2
c
c     author - xianglong yuan   alfred university
c     january 1998
c
c*********************************************************************
c

      implicit real*8 (a-h,o-z)

      parameter (pi=3.14159265358979d0)

      dimension cell(9)

c
c calculate cell parameters

      aa=cell(1)*cell(1)+cell(2)*cell(2)+cell(3)*cell(3)
      bb=cell(4)*cell(4)+cell(5)*cell(5)+cell(6)*cell(6)
      cc=cell(7)*cell(7)+cell(8)*cell(8)+cell(9)*cell(9)
      ab=cell(1)*cell(4)+cell(2)*cell(5)+cell(3)*cell(6)
      ac=cell(1)*cell(7)+cell(2)*cell(8)+cell(3)*cell(9)
      bc=cell(4)*cell(7)+cell(5)*cell(8)+cell(6)*cell(9)

      a=sqrt(aa)
      b=sqrt(bb)
      c=sqrt(cc)

c equation to calculate the angle 'gamma'
c
c    AB*AB = OA*OA + OB*OB - 2*OA*OB*cos(gamma)
c    AB*AB = OA*OA - 2*ab + OB*OB
c       ==> ab = a*b*cos(gamma)

      alpha=acos(bc/b/c)*180.0d0/pi
      beta =acos(ac/a/c)*180.0d0/pi
      gamma=acos(ab/a/b)*180.0d0/pi

      write(ncssr,'(38x,3f8.3)') a,b,c
      write(ncssr,"(21x,3f8.3,4x,'SPGR=  1 P 1         OPT = 1')")
     +  alpha,beta,gamma 

      end


      subroutine cellrotate(celorg,cvt)
c
c*********************************************************************
c
c     calculate the conversion matrix to rotate the cell and make
c     the new cell matrix have the form
c
c                       | A11 A12 A13 |
c                       |  0  A22 A23 |
c                       |  0   0  A33 |
c
c     which is used by CERIUS 2 of MSI to orientate unit cell
c
c     author - xianglong yuan   alfred university
c     january 1998
c
c*********************************************************************
c
      implicit real*8 (a-h,o-z)

      dimension cell(9),celorg(9),cel_inverse(9),cvt(9)

c
c calculate parameters of original cell

      aa=celorg(1)*celorg(1)+celorg(2)*celorg(2)+celorg(3)*celorg(3)
      bb=celorg(4)*celorg(4)+celorg(5)*celorg(5)+celorg(6)*celorg(6)
      cc=celorg(7)*celorg(7)+celorg(8)*celorg(8)+celorg(9)*celorg(9)
      ab=celorg(1)*celorg(4)+celorg(2)*celorg(5)+celorg(3)*celorg(6)
      ac=celorg(1)*celorg(7)+celorg(2)*celorg(8)+celorg(3)*celorg(9)
      bc=celorg(4)*celorg(7)+celorg(5)*celorg(8)+celorg(6)*celorg(9)

c
c calculate elements of new rotated cell
c
c for both original and new cell, length of OA, OB, OC, AB, AC or BC
c is not changed, i.e., aa, bb, cc, ab, ab, or bc is the same

      cell(9)=sqrt(cc)                                   ! OC = OC'
      cell(8)=0.0d0
      cell(7)=0.0d0
      cell(6)=bc/cell(9)                                 ! BC = BC'
      cell(5)=sqrt(bb-cell(6)*cell(6))                   ! OB = OB'
      cell(4)=0.0d0
      cell(3)=ac/cell(9)                                 ! AC = AC'
      cell(2)=(ab-cell(3)*cell(6))/cell(5)               ! AB = AB'
      cell(1)=sqrt(aa-cell(2)*cell(2)-cell(3)*cell(3))   ! OA = OA'

c
c calculate the inverse matrix of original cell

      call inverse(celorg,cel_inverse,3)

c
c calculate the conversion matrix
c
c   [celorg]*[cvt] = [celnew]
c  ==> [cvt] = [celorg inverse]*[celnew]

      do i=1,3
      do j=1,3
        cvt(i*3-3+j)=0.0d0     ! matrix element: cvt(i,j) = cvt(i*3-3+j)
        do k=1,3
          cvt(i*3-3+j)=cvt(i*3-3+j)+cel_inverse(i*3-3+k)*cell(k*3-3+j)
        enddo
      enddo
      enddo

      end

      subroutine inverse(rrr,ppp,n)
c
c*********************************************************************
c
c     invert matrix R to its inverse P, i.e., R*P = I
c
c     use Faddeev-Leverrier matrix inversion method
c
c     author - xianglong yuan   alfred university
c     january 1998
c
c*********************************************************************
c
      implicit real*8 (a-h,o-z)

      parameter (nmax=10)

      dimension rrr(*),ppp(*),qqq(nmax*nmax)

      if(n.gt.nmax) then
        write(*,'(a,a)') 'inversion matrix exceeds max.',
     x    'in inverse subroutine'
        stop
      endif

      nn=n*n

c
c     set Q = 0

      do i=1,nn
        qqq(i)=0.0d0
      enddo

c
c     set leading polynomial coefficient to unity

      C = 1.0d0

c
c     find the n P matrices

      do l=1,n

        do i=1,nn               ! form P = Q + C*I
          ppp(i)=qqq(i)
        enddo
        do j=1,n
          ppp(j*n-n+j)=ppp(j*n-n+j) + C
        enddo           ! matrix element: ppp(i,j)=ppp(i*n-n+j)

        do i=1,n                ! Q = R*P
        do j=1,n
          qqq(i*n-n+j)=0.0d0
          do k=1,n
            qqq(i*n-n+j)=qqq(i*n-n+j)+rrr(i*n-n+k)*ppp(k*n-n+j)
          enddo
        enddo
        enddo

        C = 0.0d0               ! form C = (-1/L)*trace(Q)
        do k=1,n
          C = C - qqq(k*n-n+k)
        enddo
        C = C/float(l)

      enddo

c     stop if R is singular

      if(abs(C).lt.1.0d-5) then
        write(*,'(a)') 'matrix is singular!'
        stop
      endif

c     form P = R inverse

      do i=1,nn
        ppp(i) = -ppp(i)/C
      enddo

      end

#########################################################################

README:

Examples

1). For CONFIG file (a-quartz.cfg1) with cell matrix form

            | A11 A12  0  |
            |  0  A22  0  |,
            |  0   0  A33 |

both converted CSSR files, 'a-quartz1-con.cssr' from 'con_cssr.f' and
'a-quartz1-cfg.cssr' from 'cfg_cssr.f' work properly. This is because the
cell matrix has the same format used in CERIUS 2.

2). For CONFIG file (a-quartz.cfg2) with cell matrix form

            | A11  0   0  |
            | A21 A22  0  |,
            |  0   0  A33 |

the CSSR file 'a-quartz2-con.cssr' made from 'con_cssr.f' does not work
properly. The translation by CERIUS 2 distroys the quartz structure.
However, the rotated CSSR file 'a-quartz2-cfg.cssr' made from 'cfg_cssr.f'
works perfectly and is exactly the same as the previous
'a-quartz1-cfg.cssr'.

3). Actually, both CONFIG file (a-quartz.cfg1 and a-quartz.cfg2) are the
same, except for the different cell orientation. The CSSR input format for 
cell vectors cannot tell the cell orientation. The cell orientation of
CERIUS 2 is fixed.


Xianglong Yuan

#########################################################################
file: a-quartz.cfg1

CONFIG: alpha-quartz SiO2 with cell matrix-1
         0         3
        4.2547830000       -2.4565000000        0.0000000000
        0.0000000000        4.9130000000        0.0000000000
        0.0000000000        0.0000000000        5.4050000000
O2-              1
   -0.3676132512       -0.9295396000        -1.545289500    
O2-              2
    0.9888115692        -2.310092600        0.2563591500    
O2-              3
     1.506193182       -0.4451178000         2.058061850    
O2-              4
   -0.9888115692        0.1464074000         1.545289500    
O2-              5
    -1.506193182         2.011382200       -0.2563591500    
O2-              6
    0.3676132512         1.526960400        -2.058061850    
Si4+             7
   -0.1272180117        -2.383050650       -0.9008513500    
Si4+             8
    -2.127391500         1.081351300        0.9008513500    
Si4+             9
    0.1272180117        0.7344935000E-01    -2.702500000    

#########################################################################
file: a-quartz.cfg2

CONFIG: alpha-quartz SiO2 with cell matrix-2
         0         3
        4.9130000000        0.0000000000        0.0000000000
       -2.4565000000        4.2547830000        0.0000000000
        0.0000000000        0.0000000000        5.4050000000
O2-              1
    0.1464074000       -0.9888115692        -1.545289500    
O2-              2
     2.011382200        -1.506193182        0.2563591500    
O2-              3
     1.526960400        0.3676132512         2.058061850    
O2-              4
   -0.9295396000       -0.3676132512         1.545289500    
O2-              5
    -2.310092600        0.9888115692       -0.2563591500    
O2-              6
   -0.4451178000         1.506193182        -2.058061850    
Si4+             7
     1.081351300        -2.127391500       -0.9008513500    
Si4+             8
    -2.383050650       -0.1272180117        0.9008513500    
Si4+             9
    0.7344935000E-01    0.1272180117        -2.702500000    

#################################################################################
file: a-quartz1-con.cssr

                                         4.913   4.913   5.405
                       90.000  90.000 120.000    SPGR=  1 P 1         OPT = 1
   9   1 Created from CONFIG file
     2nd Title Line :
   1 O2     -0.36761  -0.92954  -1.54529    0   0   0   0   0   0   0   0   0.000
   2 O2      0.98881  -2.31009   0.25636    0   0   0   0   0   0   0   0   0.000
   3 O2      1.50619  -0.44512   2.05806    0   0   0   0   0   0   0   0   0.000
   4 O2     -0.98881   0.14641   1.54529    0   0   0   0   0   0   0   0   0.000
   5 O2     -1.50619   2.01138  -0.25636    0   0   0   0   0   0   0   0   0.000
   6 O2      0.36761   1.52696  -2.05806    0   0   0   0   0   0   0   0   0.000
   7 Si     -0.12722  -2.38305  -0.90085    0   0   0   0   0   0   0   0   0.000
   8 Si     -2.12739   1.08135   0.90085    0   0   0   0   0   0   0   0   0.000
   9 Si      0.12722   0.07345  -2.70250    0   0   0   0   0   0   0   0   0.000

#################################################################################
file: a-quartz1-cfg.cssr

                                         4.913   4.913   5.405
                       90.000  90.000 120.000    SPGR=  1 P 1         OPT = 1
   9   1 Created from CONFIG file
CONFIG: alpha-quartz SiO2 with cell matrix-1                                    
   1 O2      1.75978   0.29871   1.15721    0   0   0   0   0   0   0   0   0.000
   2 O2      3.11620  -1.08184   2.95886    0   0   0   0   0   0   0   0   0.000
   3 O2      3.63358   0.78313   4.76056    0   0   0   0   0   0   0   0   0.000
   4 O2      1.13858   1.37466   4.24779    0   0   0   0   0   0   0   0   0.000
   5 O2      0.62120   3.23963   2.44614    0   0   0   0   0   0   0   0   0.000
   6 O2      2.49500   2.75521   0.64444    0   0   0   0   0   0   0   0   0.000
   7 Si      2.00017  -1.15480   1.80165    0   0   0   0   0   0   0   0   0.000
   8 Si      0.00000   2.30960   3.60335    0   0   0   0   0   0   0   0   0.000
   9 Si      2.25461   1.30170   0.00000    0   0   0   0   0   0   0   0   0.000

#################################################################################
file: a-quartz2-con.cssr

                                         4.913   4.913   5.405
                       90.000  90.000 120.000    SPGR=  1 P 1         OPT = 1
   9   1 Created from CONFIG file
     2nd Title Line :
   1 O2      0.14641  -0.98881  -1.54529    0   0   0   0   0   0   0   0   0.000
   2 O2      2.01138  -1.50619   0.25636    0   0   0   0   0   0   0   0   0.000
   3 O2      1.52696   0.36761   2.05806    0   0   0   0   0   0   0   0   0.000
   4 O2     -0.92954  -0.36761   1.54529    0   0   0   0   0   0   0   0   0.000
   5 O2     -2.31009   0.98881  -0.25636    0   0   0   0   0   0   0   0   0.000
   6 O2     -0.44512   1.50619  -2.05806    0   0   0   0   0   0   0   0   0.000
   7 Si      1.08135  -2.12739  -0.90085    0   0   0   0   0   0   0   0   0.000
   8 Si     -2.38305  -0.12722   0.90085    0   0   0   0   0   0   0   0   0.000
   9 Si      0.07345   0.12722  -2.70250    0   0   0   0   0   0   0   0   0.000

#################################################################################
file: a-quartz2-cfg.cssr

                                         4.913   4.913   5.405
                       90.000  90.000 120.000    SPGR=  1 P 1         OPT = 1
   9   1 Created from CONFIG file
CONFIG: alpha-quartz SiO2 with cell matrix-2                                    
   1 O2      1.75978   0.29871   1.15721    0   0   0   0   0   0   0   0   0.000
   2 O2      3.11620  -1.08184   2.95886    0   0   0   0   0   0   0   0   0.000
   3 O2      3.63358   0.78313   4.76056    0   0   0   0   0   0   0   0   0.000
   4 O2      1.13858   1.37466   4.24779    0   0   0   0   0   0   0   0   0.000
   5 O2      0.62120   3.23963   2.44614    0   0   0   0   0   0   0   0   0.000
   6 O2      2.49500   2.75521   0.64444    0   0   0   0   0   0   0   0   0.000
   7 Si      2.00017  -1.15480   1.80165    0   0   0   0   0   0   0   0   0.000
   8 Si      0.00000   2.30960   3.60335    0   0   0   0   0   0   0   0   0.000
   9 Si      2.25461   1.30170   0.00000    0   0   0   0   0   0   0   0   0.000

#################################################################################


