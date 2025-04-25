      subroutine start
     x  (idnode,mxnode,mcore,msp,natm,ncon,kmax,ntime,ngap,lden,lcor,
     x   ltim,lchg,kill,tstep)

c***********************************************************************
c
c     write control variables and check parameters
c
c     author w smith 1994
c     copyright daresbury laboratory 1994
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      logical lden,ltim,lcor,kill,lchg

      if(idnode.eq.0)then

         if(lchg)then
            write(3,'(a)')'Charge density option selected'
         else
            write(3,'(a)')'Particle density option selected'
         endif

      endif

c     check on control parameters

      natm=msp
      iman=(msp+mxnode-1)/mxnode
      klim=((2*kmax+1)**3-1)/2
      nkt=(klim+mxnode-1)/mxnode
      m1=0
      m2=0
      m3=0
      if(lden)m1=6*iman*(kmax+2)+6*klim
      if(lcor)m2=4*nkt*ntime+2*klim+2*ntime
      if(ltim)m3=ntime*(2*nkt+19)
      kcore=max0(m1,m2,m3)
      if(mcore.lt.kcore)then

         if(idnode.eq.0)write(3,'(a,i8,a)')
     x        'Error - insufficient core allocated. ',kcore,
     x        ' words required'
         kill=.true.
         
      endif

      if(idnode.eq.0)then

         write(3,'(a,i10)')
     x     'Maximum size of dynamic core area         = ',mcore
         write(3,'(a,i10)')
     x     'Maximum number of configurations          = ',ncon
         write(3,'(a,i10)')
     x     'Largest k vector component                = ',kmax
         write(3,'(a,i10)')
     x     'Length of correlation arrays              = ',ntime
         write(3,'(a,i10)')
     x     'Interval between time origins             = ',ngap

      endif

      m=ngap
      do i=1,20

         m=m/2
         if(m.eq.0)go to 100

      enddo
 100  m=2**(i-1)

      if(ngap.ne.m)then

         if(idnode.eq.0)write(3,'(a,i6)')
     x        'Warning - origin interval reset to',m
         ngap=m

      endif

      if(ngap.gt.1)then

         m=ntime
         do i=1,30
            
            m=m/ngap
            if(m.eq.0)go to 200
            
         enddo
 200     m=ngap**(i-1)
         
         if(ntime.ne.m)then
            
            if(idnode.eq.0)write(3,'(a,i6)')
     x           'Warning - length of correlation arrays reset to',m
            ntime=m
            
         endif

      endif

      return
      end
