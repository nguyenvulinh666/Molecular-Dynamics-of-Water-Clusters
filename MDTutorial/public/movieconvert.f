      subroutine movieconvert
c
c     subroutine to write a CERIUS 2 trajectory file
c
c     Thijs J.H. Vlugt
c     Department of Chemical Engineering
c     University of Amsterdam
c     Nieuwe Achtergracht 166
c     1018 WV Amsterdam
c     The Netherlands
c     email  : tampert@molsim.chem.uva.nl
c
      implicit none

      character*4       hdr,hdr2
      character*80      eex
      character*5       l
      character*7       ll
      real*8            rr(25000),yy(25000),zz(25000)
      integer*4         a,b(20),i,kkk

      open(1,file='movie.trj',form='unformatted',status='unknown')

c    movie.bin is a sample file with coordinates in binary format; see next

      open(2,file='movie.bin',form='unformatted',status='unknown')

      hdr='mdtr'
      do i=1,20
         b(i)=0
      enddo
      b(1)=2010

      write(1) hdr,(b(i),i=1,20)

      a=1
      eex='created by cerius2                               '
     &    //'                               '
      write(1) a,eex

      a=0
      write(1) a

      do i=1,8
         b(i)=0
      enddo
      write(1) (b(i),i=1,8)

c    natom is the number of atoms

      b(1)=1
      b(2)=natom
      b(3)=natom
      hdr='mode'
      hdr2='   '
      write(1) b(1),b(2),b(3),hdr,hdr2
      write(1) b(2),(i,i=1,natom)

      b(1)=7
      b(2)=5
      ll='notitle'
      l='nopar'
      write(1) b(1),ll
      write(1) b(2),l

      do kkk=1,nframe

         do i=1,58
            rr(i)=0.0d0
         enddo

         do i=1,6
            b(i)=0
         enddo

         write(1) dble(kkk*1.0d-2),kkk,rr(2),(rr(i),i=3,58),(b(i),i=1,6)
         write(1) (rr(i),i=1,12)

         read(2)  (rr(i),i=1,natom)
         read(2)  (yy(i),i=1,natom)
         read(2)  (zz(i),i=1,natom)

         write(1) (rr(i),i=1,natom)
         write(1) (yy(i),i=1,natom)
         write(1) (zz(i),i=1,natom)

      enddo

      close(1)
      close(2)

      return
      end
