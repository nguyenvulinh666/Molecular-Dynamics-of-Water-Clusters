      subroutine hread
     x  (new,history,cfgname,name,iflg,imcon,keytrj,natms,
     x   nstep,tstep,cell,chge,weight,xyz,vel,frc,status)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading the formatted history file 
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c     single processor version
c
c     wl
c     2004/08/31 14:53:45
c     1.2
c     Exp
c     
c***********************************************************************
c     
      
      implicit none

      logical new

      character*80 cfgname
      character*40 history
      character*8 name(*),step
      
      integer iflg,imcon,keytrj,natms,nstep,status,nhist,ktrj,matms,i,j

      real*8 tstep,vx,vy,vz,fx,fy,fz

      real*8 cell(9)
      real*8 chge(*),weight(*)
      real*8 xyz(3,*),vel(3,*),frc(3,*)

      save ktrj
      data nhist/77/

      status=0

c     open history file if new job

      if(new)then
      
        open(nhist,file=history,status='old',err=100)
        
        read(nhist,'(a80)',err=200) cfgname
        write(*,'(a,a)')'History file header: ',cfgname
        read(nhist,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)then

            write(*,'(a)')'Error - no velocities in file'
            status=-1
            close (nhist)
            new=.true.
            return

          endif
          if(keytrj.gt.1)then

            write(*,'(a)')'Error - no forces in file'
            status=-2
            close (nhist)
            new=.true.
            return

          endif
        endif

        new=.false.

      endif
        
      read(nhist,'(a8,4i10,f12.6)',end=200)
     x     step,nstep,matms,ktrj,imcon,tstep

      if(matms.gt.natms)then

        write(*,'(a)')'Error - too many atoms in MD cell'
        write(*,'(a,i6,a)')'File contains',matms,' atoms'
        status=-3
        close (nhist)
        new=.true.
        return

      endif

      natms=matms
      
      if(imcon.gt.0) then

         read(nhist,*,end=200) cell(1),cell(2),cell(3)
         read(nhist,*,end=200) cell(4),cell(5),cell(6)
         read(nhist,*,end=200) cell(7),cell(8),cell(9)

      endif

      do i = 1,natms

        read(nhist,'(a8,i10,2f12.6)',end=200)
     x    name(i),j,weight(i),chge(i)
        read(nhist,*,end=200) xyz(1,i),xyz(2,i),xyz(3,i)
        if(keytrj.ge.1)then
          read(nhist,*,end=200) vel(1,i),vel(2,i),vel(3,i)
        else if(ktrj.ge.1)then
          read(nhist,*,end=200) vx,vy,vz
        endif
        if(keytrj.ge.2)then
          read(nhist,*,end=200) frc(1,i),frc(2,i),frc(3,i)
        else if(ktrj.ge.2)then
          read(nhist,*,end=200) fx,fy,fz
        endif

      enddo

      if(iflg.lt.0)then

        close (nhist)
        new=.true.

      endif

      return

  100 continue

      write(*,'(a)')'Error - History file not found'
      status=-4
      return

  200 continue
      write(*,'(a)')'Warning - end of History file encountered'
      close (nhist)
      iflg=-1
      new=.true.

      return
      end
