      subroutine fft(kill,idnode,ind,isw,ndiv,key,aaa,wfft,bbb)
c***********************************************************************
c     
c     fast fourier transform routine
c     
c     copyright daresbury laboratory 1994
c
c     author w smith
c
c***********************************************************************
      
      implicit real*8(a-h,o-z)
      parameter (irite=6)
      logical check,kill
      complex*16 aaa(ndiv),bbb(ndiv),wfft(ndiv),ttt
      dimension key(ndiv)
      data tpi/6.2831853072d0/
      
c     
c     check that array is of suitable length
      nt=1
      check=.true.
      do i=1,20
         nt=2*nt
         if(nt.eq.ndiv)then
            check=.false.
            nu=i
         endif
      enddo
      if(check)then
         if(idnode.eq.0)
     x    write(3,'(a)')'Error - number of points not a power of two'
         kill=.true.
         return
      endif
      
      if(ind.gt.0)then
c     
c     set reverse bit address array
         
         do kkk=1,ndiv
            iii=0
            jjj=kkk-1
            do j=1,nu
               jj2=jjj/2
               iii=2*(iii-jj2)+jjj
               jjj=jj2
            enddo
            key(kkk)=iii+1
         enddo
c     
c     initialise complex exponential factors
         
         tpn=tpi/dble(ndiv)
         arg=0.d0
         np1=ndiv+1
         np2=ndiv/2
         wfft(1)=(1.d0,0.d0)
         do i=1,np2
            arg=tpn*dble(i)
            wfft(i+1)=cmplx(cos(arg),sin(arg))
            wfft(np1-i)=conjg(wfft(i+1))
         enddo
         
         return
      endif
      
c     
c     take conjugate of exponentials if required
      
      if(isw.lt.0)then
         
         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo
         
      endif
      
c     
c     take copy input array
      
      do i=1,ndiv
         bbb(i)=aaa(i)
      enddo
c     
c     perform fourier transform
      
      kkk=0
      nu1=nu-1
      np2=ndiv/2
      do l=1,nu

  100    do i=1,np2
            iii=key(kkk/2**nu1+1)
            kk1=kkk+1
            k12=kk1+np2
            ttt=bbb(k12)*wfft(iii)
            bbb(k12)=bbb(kk1)-ttt
            bbb(kk1)=bbb(kk1)+ttt
            kkk=kkk+1
         enddo
         kkk=kkk+np2
         if(kkk.lt.ndiv)go to 100
         kkk=0
         nu1=nu1-1
         np2=np2/2

      enddo
c     
c     unscramble the fft using bit address array
      
      do kkk=1,ndiv
         iii=key(kkk)
         if(iii.gt.kkk)then
            ttt=bbb(kkk)
            bbb(kkk)=bbb(iii)
            bbb(iii)=ttt
         endif
      enddo
c     
c     restore exponentials to unconjugated values if necessary
      
      if(isw.lt.0)then
         
         do i=1,ndiv
            wfft(i)=conjg(wfft(i))
         enddo
         
      endif
      
      return
      end
