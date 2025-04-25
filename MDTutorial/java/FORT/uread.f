      subroutine uread
     x  (new,history,cfgname,name,iflg,imcon,keytrj,natms,
     x   nstep,tstep,cell,chge,weight,xyz,vel,frc,status)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading unformatted history files
c     
c     double precision, single processor version
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c     
c     wl
c     2004/08/31 14:53:45
c     1.3
c     Exp
c     
c***********************************************************************
c     
      
      implicit none
      
      logical new
      
      character*80 cfgname
      character*40 history
      character*8 name(*)
      integer iflg,imcon,keytrj,natms,nstep,nhist,status,matms,i,ktrj
      real*8 tstep,datms,dstep,trjkey,dimcon
      
      real*8 cell(9)
      real*8 chge(*),weight(*)
      real*8 xyz(3,*),vel(3,*),frc(3,*)
      data nhist/77/

      status=0

c     open the history file if new job
      
      if(new)then
        
        open(nhist,file=history,form='unformatted',
     x       status='old',err=100)
        
        read(nhist,err=200) cfgname
        write(*,'(a,a)')'History file header: ',cfgname
        read(nhist,end=200) datms
        matms=nint(datms)
        if(matms.gt.natms)then
          
          write(*,'(a)')'Error - too many atoms in MD cell'
          write(*,'(a,i6,a)')'File contains',matms,' atoms'
          status=-1
          close (nhist)
          return
          
        endif
        natms=matms
        read(nhist,end=200) (name(i),i=1,natms)
        read(nhist,end=200) (weight(i),i=1,natms)
        read(nhist,end=200) (chge(i),i=1,natms)
        
        new=.false.
        
      endif
      
      read(nhist,end=200)dstep,datms,trjkey,dimcon,tstep
      nstep=nint(dstep)
      ktrj=nint(trjkey)
      imcon=nint(dimcon)
      if(keytrj.gt.ktrj)then
        
        if(ktrj.eq.0)then

          write(*,'(a)')'Error - no velocities in file'
          status=-2
          close (nhist)
          new=.true.
          return

        endif
        if(keytrj.gt.1)then

          write(*,'(a)')'Error - no forces in file'
          status=-3
          close (nhist)
          new=.true.
          return

        endif
      endif
      
      if(imcon.gt.0) read(nhist,end=200) cell
      
      read(nhist,end=200) (xyz(1,i),i = 1,natms)
      read(nhist,end=200) (xyz(2,i),i = 1,natms)
      read(nhist,end=200) (xyz(3,i),i = 1,natms)
      
      if(keytrj.ge.1)then
        read(nhist,end=200) (vel(1,i),i = 1,natms)
        read(nhist,end=200) (vel(2,i),i = 1,natms)
        read(nhist,end=200) (vel(3,i),i = 1,natms)
      else if(ktrj.ge.1)then
        read(nhist,end=200)
        read(nhist,end=200)
        read(nhist,end=200)
      endif
      if(keytrj.ge.2)then
        read(nhist,end=200) (frc(1,i),i = 1,natms)
        read(nhist,end=200) (frc(2,i),i = 1,natms)
        read(nhist,end=200) (frc(3,i),i = 1,natms)
      else if(ktrj.ge.2)then
        read(nhist,end=200)
        read(nhist,end=200)
        read(nhist,end=200)
      endif

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
