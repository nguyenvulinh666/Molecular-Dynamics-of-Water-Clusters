      program con_cssr

c*********************************************************************
c
c	dl_poly routine for converting CONFIG files to CSSR files
c	suitable for input to CERIUS 2 by MSI
c
c	copyright daresbury laboratory feb 1996
c	author d. collins feb. 1996
c
c*********************************************************************

      implicit real*8 (a-h,o-z)

      logical find, lperiodic, firstcall

      character*80 string, title
      character*11 outfile
      character*8 alabel

      dimension alabel(10000)
      dimension xlat(3,3),xbas(10000,3),rorigin(3)

      data lfsearch/.false./, firstcall/.true./

c
c STEP 1  open CONFIG      file   channel 7
c         open config.cssr file   channel 8

c         read title
c         read keytrj, imcon
c
      inquire(file='CONFIG',exist=find)
      if (.not.find) then
        write(*,'(''CONFIG file not found'')')
        stop
      end if
      open(7,file='CONFIG',status='old')
      open(8,file='config.cssr',status='unknown')

      read(7,'(a)',err=501, end=502) string
      write(*,'(a)') string
      title=string
      read(7,'(2i10)',err=503, end=504) keytrj, imcon

c perform checks on keytrj & imcon and set iformat

c 6 format possibilities
c set iformat as follows
c   1    non-periodic - coords only
c   2    non-periodic - coords & velocities
c   3    non-periodic - coords, velocities & forces
c   4    periodic - coords only
c   5    periodic - coords & velocities
c   6    periodic - coords, velocities & forces

c check keytrj

      if (keytrj.eq.0) then
         write(*,'(/,''File contains coordinates only'')')
         iformat=1
      else if (keytrj.eq.1) then
         write(*,'(/,''File contains coordinates & velocities'')')
         iformat=2
      else if (keytrj.eq.2) then
         write(*,'(/,''File contains coordinates, velocities'',
     +             '' & forces'')')
         iformat=3
      else
         write(*,'(/,''Error - invalid keytrj'')')
         stop
      end if


      if (imcon.eq.0) then
         write(*,'(''No Periodic Boundary Conditions'')')
         write(8,'(/)')
         icerkey = 1

      else if (imcon.ge.1.and.imcon.le.6) then

         icerkey = 1
         call cell

      else
         write(*,'(/,''Error - invalid imcon'')')
         stop
      end if


c
c STEP 2 read & store atomic coordinates
c

      write(*,'(/,''Reading coordinates ...'')')
      print*,'iformat =', iformat
      goto(1001,1002,1003) iformat

 1001 print*,'coords only'
      do i = 1,10000
         read(7,'(a)',err=505,end=506) alabel(i)
         read(7,*,err=507,end=508) (xbas(i,j), j=1,3)
c         write(*,'(a)') alabel(i)
c         print*, i
      end do
      goto 2000

 1002 print*,'coords & velocities'
      do i = 1,10000
         read(7,'(a)',err=505,end=506) alabel(i)
         read(7,*,err=507,end=508) (xbas(i,j), j=1,3)
         read(7,*)
      end do
      goto 2000

 1003 print*,'coords, velocities & forces'
      do i = 1,10001
         read(7,'(a)',err=505,end=506) alabel(i)
         read(7,*,err=507,end=508) (xbas(i,j), j=1,3)
         read(7,'(/)')
      end do

 2000 write(*,'(''Max no. of atoms exceeded - program terminating'')')
      stop

 506  print*,'Coords read correctly'
      natom = i -1

c determine no atoms & check for corruption to file
      write(*,'(''No. of atoms'',25x,i10)') natom


c
c STEP 3 write remaining part of cssr file
c        natom
c        coordinates
c
      iz = 0
      rz=0.0d0

      write(8,'(i4,3x,i1,1x''Created from CONFIG file'')')
     +       natom, icerkey
      write(8,'(5x,''2nd Title Line :'')')

      do i = 1, natom

         write(8,'(i4,1x,a2,4x,3(f9.5,1x),8(i4),1x,f7.3)')
     +        i, alabel(i)(1:2), (xbas(i,j), j=1,3),
     +   iz,iz,iz,iz,iz,iz,iz,iz,rz

      end do

      print*,'good end'
      stop


c print error messages

 501  write(*,'(''Error in reading title'')')
      stop

 502  write(*,'(''Empty file'')')
      stop

 503  write(*,'(''Error in reading Record(2) - levcfg & imcon'')')
      stop

 504  write(*,'(''EOF detected while reading Record(2) '',
     +           '' - levcfg & imcon'')')
      stop

 505  write(*,'(''Error in reading atom labels'')')
      stop


 507  write(*,'(''Error while reading coordinates'')')
      stop

 508  write(*,'(''EOF detected while reading coordinates'')')
      stop

      end

      subroutine cell

      implicit real*8 (a-h,o-z)


      dimension xlat(3,3)
      dimension ra(3),rb(3),r1(3),r2(3)

      data  pi/3.14159265358979d0/,zero/0.0d0/,ang/180.0d0/

      adtb = zero
      adtc = zero
      bdtc = zero

c      open(3,file='HISTORY',status='old')
c      read(3,*)
c read lattice vectors and convert to unit cell dimensions

      do na=1,3
        read(7,*,err=501,end=502) (xlat(i,na),i=1,3)
      end do

      do na = 1,3

        ra(na) = zero

        do nb=1,3
          ra(na) = xlat(nb,na)*xlat(nb,na) + ra(na)
        end do

        ra(na) = sqrt(ra(na))
        adtb = xlat(na,1) * xlat(na,2) + adtb
        adtc = xlat(na,1) * xlat(na,3) + adtc
        bdtc = xlat(na,2) * xlat(na,3) + bdtc

      end do

      rb(3) = acos(adtb/(ra(1)*ra(2)))
      rb(1) = acos(bdtc/(ra(2)*ra(3)))
      rb(2) = acos(adtc/(ra(1)*ra(3)))

      rlat = 1.0e0
      do na=1,3
        r2(na) = rb(na) * ang/pi
        r1(na) = ra(na) * rlat
      end do

         write(8,'(38x,3f8.3)') r1(1), r1(2), r1(3)
         write(8,'(21x,3f8.3,4x,''SPGR=  1 P 1         OPT = 1'')')
     +             r2(1), r2(2), r2(3)

c       write(12,102) r1(1), r1(2), r1(3)
c       write(12,1021) r2(1), r2(2), r2(3)
c 102   format(38x,3(f8.3))
c 1021  format(21x,3(f8.3))

      return
 501  write(*,'(''Error while reading periodic data'')')
      stop

 502  write(*,'(''EOF while reading periodic data'')')
      stop

      end


