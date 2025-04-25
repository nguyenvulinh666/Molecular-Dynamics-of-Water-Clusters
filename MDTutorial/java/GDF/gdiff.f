      subroutine gdiff
     x  (dl_fmt,atnam1,atnam2,fname,title,name,newnam,mxtim,ncrs,
     x  mxstr,mxrad,mxatms,natms,nconf,isampl,iocrs,numfil,tstep,rcut,
     x 	imd,msm,avcell,gcross,gzero,xyz,vel,frc,weight,chge,cell,bcell,
     x  xy0,xy1,xy2,status)

c*********************************************************************
c     
c     dl_poly program to calculate cross term of Van Hove density 
c     correlation function for selected atoms 
c     
c     copyright daresbury laboratory 1995
c     author  w.smith march 1995
c     
c     wl
c     2006/11/01 09:49:18
c     1.4
c     Exp
c
c*********************************************************************
      
      implicit none
      
      real*8 pi
      parameter (pi=3.1415926536d0)

      logical dl_fmt,all,new,lflag
      character*2 fnum
      character*40 fname
      character*80 title
      integer mxrad,ncrs,mxtim,mxatms,natms,nconf,isampl,iocrs,ii,jj
      integer status,numfil,nogcrs,i,j,k,m,lsr,msr,ipass,iflg
      integer iconf,keytrj,imcon,nstep,matms,nsgcrs,npts,mxstr,nta
      integer ntb,nstp0
      integer imd(mxtim),msm(mxtim)
      character*8 atnam1,atnam2,name(mxatms),newnam(mxstr)
      real*8 tstep,rcut,rcut2,delr,timstp,uuu,vvv,www,det,rmsx,rmsy
      real*8 rmsz,rsq,rrr,time,volm,rnorm,f1,f2,tcut
      real*8 gcross(mxtim,mxrad),gzero(mxtim,mxstr,3)
      real*8 xyz(3,mxatms),vel(3,mxatms),frc(3,mxatms),rcell(9)
      real*8 weight(mxatms),chge(mxatms),avcell(9),cell(9),bcell(10)
      real*8 xy0(3,mxstr),xy1(3,mxstr),xy2(3,mxstr)

      status=0
      lflag=.true.
      call strcut(atnam1,8)
      call strcut(atnam2,8)
      call strcut(fname,40)

      all=(atnam1.eq."ALL".or.atnam2.eq."ALL".or.
     x     atnam1.eq."all".or.atnam2.eq."all")

      if(all)then

         atnam1="ALL"
         atnam2="ALL"

      endif

c     open report file

      open(9,file='GDFPRG.rep')

      write(9,'(a)')
     x  'van Hove Distinct Correlation Function Program'
      

c     check on specified control variables
      
      if(ncrs.gt.mxtim)then
        
        ncrs=mxtim
        write(9,'(a,i6)')
     x    'Warning - Correlation array length set to ',ncrs
        
      endif

      if(mod(ncrs,iocrs).ne.0)then
        
        ncrs=max(iocrs*(ncrs/iocrs),1)
        write(9,'(a,i6)')
     x    'Warning - 1st dimension of Gcross array reset to ',ncrs
        
      endif
      
      rcut2=rcut**2
      delr=rcut/dble(mxrad)
      nogcrs=max(ncrs/iocrs,1)
      
      write(9,'(a,a40)')'Label  of target HISTORY file : ',fname
      write(9,'(a,a8)')'Label  of first atom          : ',atnam1
      write(9,'(a,a8)')'Label  of second atom         : ',atnam2
      write(9,'(a,i6)')'Length of correlation arrays  : ',ncrs
      write(9,'(a,i6)')'Number of requested configs.  : ',nconf
      write(9,'(a,i6)')'Sampling interval             : ',isampl
      write(9,'(a,i6)')'Interval between origins      : ',iocrs
      write(9,'(a,f6.3)')'Selected correlation radius   : ',rcut
      write(9,'(a,f6.3)')'Bin width in correlation      : ',delr

c     initialise gcross arrays

      do j=1,mxtim
        msm(j)=0
        do i=1,mxrad
          
          gcross(j,i)=0.d0
          
        enddo
      enddo

c     initialize average cell vectors

      do i=1,9

         avcell(i)=0.d0

      enddo

c     multiple pass over HISTORY file

      do ipass=0,isampl-1

        lsr=0
        msr=0
        nsgcrs=0
        new=.true.
        
c     set default cell vectors
        
        do i=1,9
          cell(i)=0.d0
        enddo
        cell(1)=1.d0
        cell(5)=1.d0
        cell(9)=1.d0

        iflg=0
        keytrj=0
        
        do iconf=0,nconf-1

          if(dl_fmt)then

            call hread
     x        (new,fname,title,name,iflg,imcon,keytrj,natms,
     x        nstep,timstp,cell,chge,weight,xyz,vel,frc,status)

          else

            call uread
     x        (new,fname,title,name,iflg,imcon,keytrj,natms,
     x        nstep,timstp,cell,chge,weight,xyz,vel,frc,status)

          endif

c     check cutoff radius

          if(iconf.eq.0)then

            call dcell(cell,bcell)
            
            tcut=0.5d0*min(bcell(7),bcell(8),bcell(9))
            if(rcut.gt.tcut)then
              
              write(9,'(a,f10.6)')'Warning - RDF cutoff reset to',tcut
              rcut=tcut
              rcut2=rcut**2
              delr=rcut/dble(mxrad)
              
            endif

          endif
          if(iconf.eq.0)nstp0=nstep
          if(iconf.eq.1)then

            if(timstp.lt.1.d-8)timstp=1.d0
            tstep=timstp*dble(nstep-nstp0)
            if(ipass.eq.0)write(9,'(a,f6.3)')
     x        'Time interval between configs : ',tstep

          endif

          if(iflg.lt.0)go to 100

          if(imcon.ge.1.and.imcon.le.3)then

            call invert(cell,rcell,det)
            
          else
            
            write(9,'(a)')
     x        'Error - incorrect boundary condition'
            status=-2
            return
            
          endif

          j=0
          do i=1,natms

            if(all.or.name(i).eq.atnam1.or.name(i).eq.atnam2)then

              j=j+1
              if(j.gt.mxstr)then

                if(lflag.and.iconf.eq.0)write(9,'(a,2i5)')
     x            'Warning - too many atoms of selected type',j,mxstr
                lflag=.false.

              else

                newnam(j)=name(i)
                xy2(1,j)=xyz(1,i)*rcell(1)+xyz(2,i)*rcell(4)+
     x            xyz(3,i)*rcell(7)
                xy2(2,j)=xyz(1,i)*rcell(2)+xyz(2,i)*rcell(5)+
     x            xyz(3,i)*rcell(8)
                xy2(3,j)=xyz(1,i)*rcell(3)+xyz(2,i)*rcell(6)+
     x            xyz(3,i)*rcell(9)
                
              endif

            endif

          enddo

          matms=j

c     running average of cell vectors

          f1=(dble(iconf)/dble(iconf+1))
          f2=1.d0/dble(iconf+1)
          do i=1,9
             
             avcell(i)= f1*avcell(i)+f2*cell(i)
             
          enddo

          if(iconf.eq.ipass)then

c     set initial positions of atoms

            do i=1,matms

               xy0(1,i)=xy2(1,i)
               xy0(2,i)=xy2(2,i)
               xy0(3,i)=xy2(3,i)

            enddo

          else if(iconf.gt.ipass)then
            
c     accumulate incremental distances
            
            do i=1,matms
               
               uuu=xy2(1,i)-xy1(1,i)
               vvv=xy2(2,i)-xy1(2,i)
               www=xy2(3,i)-xy1(3,i)
            
               uuu=uuu-nint(uuu)
               vvv=vvv-nint(vvv)
               www=www-nint(www)
            
               xy0(1,i)=xy0(1,i)+uuu
               xy0(2,i)=xy0(2,i)+vvv
               xy0(3,i)=xy0(3,i)+www
              
            enddo
            
          endif
          
          do i=1,matms
            
             xy1(1,i)=xy2(1,i)
             xy1(2,i)=xy2(2,i)
             xy1(3,i)=xy2(3,i)
            
          enddo
          
c     calculate cross correlation function
          
          if(mod(iconf,isampl).eq.ipass)then
            
            if(mod(nsgcrs,iocrs).eq.0)then
              
              lsr=min(lsr+1,nogcrs)
              msr=mod(msr,nogcrs)+1
              imd(msr)=1
              do i=1,matms
                
                gzero(msr,i,1)=xy0(1,i)
                gzero(msr,i,2)=xy0(2,i)
                gzero(msr,i,3)=xy0(3,i)
                
              enddo
              
            endif
            
            nsgcrs=nsgcrs+1
            
            do j=1,lsr
              
              m=imd(j)
              imd(j)=m+1
              msm(m)=msm(m)+1
              do ii=1,matms
                if(all.or.atnam2.eq.newnam(ii))then
                  do jj=1,matms
                    
                    if((ii.ne.jj).and.(all.or.atnam1.eq.newnam(jj)))then
                      
                      uuu=xy0(1,ii)-gzero(j,jj,1)
                      vvv=xy0(2,ii)-gzero(j,jj,2)
                      www=xy0(3,ii)-gzero(j,jj,3)
                      
                      uuu=uuu-nint(uuu)
                      vvv=vvv-nint(vvv)
                      www=www-nint(www)
                      
                      rmsx=uuu*avcell(1)+vvv*avcell(4)+www*avcell(7)
                      rmsy=uuu*avcell(2)+vvv*avcell(5)+www*avcell(8)
                      rmsz=uuu*avcell(3)+vvv*avcell(6)+www*avcell(9)
                      
                      rsq=rmsx**2+rmsy**2+rmsz**2
                      
                      if(rsq.lt.rcut2)then
                        
                        k=int(sqrt(rsq)/delr)+1
                        gcross(m,k)=gcross(m,k)+1.d0
                        
                      endif
                      
                    endif
                    
                  enddo
                endif
              enddo
              
            enddo
            
          endif

        enddo

  100   continue

      enddo

      if(iflg.lt.0)iconf=iconf-1
      npts=min(nsgcrs,ncrs)
      write(9,'(a,i6)')'Number of configurations sampled: ',iconf
      
c     normalise and print correlation function

      volm=abs(avcell(1)*(avcell(5)*avcell(9)-avcell(6)*avcell(8))
     x  +avcell(4)*(avcell(3)*avcell(8)-avcell(2)*avcell(9))
     x  +avcell(7)*(avcell(2)*avcell(6)-avcell(3)*avcell(5)))

      if(all)then

         nta=matms
         ntb=matms

      else

         nta=0
         ntb=0
         do i=1,matms

            if(newnam(i).eq.atnam1)nta=nta+1
            if(newnam(i).eq.atnam2)ntb=ntb+1
            
         enddo

      endif

      tstep=tstep*dble(isampl)

c     open the correlation data files

      numfil=numfil+1
      if(numfil.gt.99)then
         fnum="XX"
         close (88)
         open(88,file='HOVGDF.'//fnum)
      else
         write(fnum,'(i2.2)')numfil
         close(88)
         open(88,file='HOVGDF.'//fnum)
      endif
      write(88,'(a80)')title
      if(atnam1.eq.atnam2) then

        rnorm=volm/(4.d0*pi*delr**3*dble(nta)*dble(nta))
        write(88,'(2a8,1pe14.6)')atnam1,atnam1,rcut
        
      else 
        
        rnorm=volm/(4.d0*pi*delr**3*dble(nta)*dble(ntb))
        write(88,'(2a8,1pe14.6)')atnam1,atnam2,rcut
        
      endif

      write(88,'(2i10)')npts,mxrad
        
      do j=1,npts
          
        do i=1,mxrad
          
          gcross(j,i)=(rnorm*gcross(j,i)/dble(msm(j)))
     x      /((dble(i)-0.5d0)**2+1.d0/12.d0)
          
        enddo
        
      enddo
      
      do j=1,npts
          
        time=tstep*dble(j)
        write(88,'(a8,1p,e14.6)')'time    ',time
        
        do i=1,mxrad
          
          rrr=delr*(dble(i)-0.5d0)
          write(88,'(1p,2e14.6)')rrr,gcross(j,i)
          
        enddo
        
      enddo
      
      close(88)
      write(9,'(a)')'G diff file HOVGDF.'//fnum//' created'
      write(9,'(a,i5)')'Number of functions per file :',npts
      close (9)

      return
      end
      subroutine strcut(a,n)

      character*1 a(n)
      logical y
      
      y=.false.
      do i=1,n

        if(ichar(a(i)).eq.0)y=.true.
        if(y)a(i)=" "

      enddo
      return
      end
