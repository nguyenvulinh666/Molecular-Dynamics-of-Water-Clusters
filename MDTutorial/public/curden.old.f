c***********************************************************************
c
c     daresbury laboratory ccp5 program for the calculation
c     of correlation functions. this version calculates the
c     fourier transform of the current density in space and
c     time
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,lcor,ltim
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common space(mcore)
      call second(timrun)
c
c     read in control variables
      call start
      call second(timrun)
      write(irite,10)timrun
   10 format('time after start  = ',f15.8)
c
c     calculate fourier transform of current density
      if(lden)call fourdc
      call second(timrun)
      write(irite,20)timrun
   20 format('time after fourdc = ',f15.8)
c
c     calculate current density correlation function
      if(lcor)call correl
      call second(timrun)
      write(irite,30)timrun
   30 format('time after correl = ',f15.8)
c
c     calculate correlation function fourier transform
      if(ltim)call curfor
      call second(timrun)
      write(irite,40)timrun
   40 format('time after curfor = ',f15.8)
c
c     print out results
      call output
      stop
      end
      blockdata
c***********************************************************************
c
c     set program parameters
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      data iread,irite,iconf,isave0,isave1,isave2/5,6,7,8,9,10/
      end
      subroutine start
c***********************************************************************
c
c     read in control variables and check parameters
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      logical lden,ltim,lcor,kill
      dimension title(10)
      data twopi/6.2831853072d0/
      open(iread,file='cur_in')
      kill=.false.
      read(iread,10)title
   10 format(10a8)
      write(irite,20)title
   20 format(
     x       ' ccp5 program to  ',/,                  
     x       ' calculate space  ',/,                  
     x       ' and time fourier ',/,
     x       ' transforms of    ',/,                  
     x       ' current density  ',/,                  
     x       ' w.smith march 82 ',/,
     x       10a8,/,)

      read(iread,*)nsteps
      read(iread,*)lden
      read(iread,*)lcor
      read(iread,*)ltim
      read(iread,*)time
c
c     adjust number of timesteps
      nor=nsteps/ngap
      nsteps=ngap*nor
c
c     check on control parameters
      write(irite,30)mcore,nsp,nsteps,kmax,klim,nkt,ntime,time
   30 format('maximum size of dynamic core area         = ',i12,/,  
     x       'number of species in simulation box       = ',i12,/,  
     x       'number of timesteps in simulation data    = ',i12,/,  
     x       'largest k vector component                = ',i12,/,  
     x       'number of (positive) k vectors            = ',i12,/,  
     x       'number of time fourier transforms in core = ',i12,/,  
     x       'number of time steps in fourier transform = ',i12,/,  
     x       'magnitude of time step in results  (ps)   = ',e12.4)  
c
c     calculate omega interval
      if(time.gt.0.d0)go to 35
      write(irite,32)
   32 format('error - time step not assigned nonzero value')        
      kill=.true.
      go to 38
   35 omeg=twopi/(time*ntime*4)
      write(irite,36)omeg
   36 format('magnitude  omega  step in results  (1/ps) = ',e12.4)
   38 m1=3*nsp*(kmax+2)+4*klim
      m2=ntime*(8*nkt+1)+4*klim
      m3=ntime*(4*nkt+15)
      m4=4*ntime*nkt+4*klim
      kcore=2*max0(m1,m2,m3,m4)
      if(mcore.ge.kcore)go to 50
      write(irite,40)kcore
   40 format('error - isufficient core allocated. ',i8,' words',    
     x           ' required')
      kill=.true.
   50 jlim=((2*kmax+1)**3+1)/2
      if(klim.eq.jlim)go to 70
      write(irite,60)
   60 format('error - parameter klim .ne. ((2*kmax+1)**3+1)/2')     
      kill=.true.
   70 if(mod(klim,nkt).eq.0)go to 90
      write(irite,80)
   80 format('error - param. klim not divisible by param. nkt')     
      kill=.true.
   90 if(nbits(ntime).eq.1)go to 110
      write(irite,100)
  100 format('error - parameter ntime not of the form 2**n')        
      kill=.true.
  110 if(nsteps.gt.ntime)go to 130
      write(irite,120)
  120 format('error - parameter ntime greater than nsteps')         
      kill=.true.
  130 if(nbits(ngap).eq.1)go to 150
      write(irite,140)
  140 format('error - parameter ngap not of the form 2**m')         
      kill=.true.
  150 if(kill)stop
      return
      end
      subroutine fourdc
c***********************************************************************
c
c     calculate spatial fourier transform of current density
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      character*8 word
      common/cmcrds/x(nsp),y(nsp),z(nsp)
      common/cmvelo/vx(nsp),vy(nsp),vz(nsp)
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common ekx(nsp,0:kmax),eky(nsp,0:kmax),ekz(nsp,0:kmax),elc(nsp),
     x     emc(nsp),enc(nsp),vxyz(4,klim)
      logical lden,ltim,lcor
      complex*16 ekx,eky,ekz,elc,emc,enc,vxyz
      complex*16 argc,ckx,cky,ckz
      data cl,pi,twopi/2.d0,3.1415926536d0,6.2831853072d0/
      rcl=twopi/cl
      nnn=nsp
      rnsp=sqrt(1.d0/dble(nsp))
c
c     read in configuration data
      open(iconf,file='HISTORY',status='OLD')
      open(isave0,file='CUR_SPC',form='unformatted')

      read(iconf,'(a)')word
      read(iconf,'(a)')word
      do 170 kstep=1,nsteps
      read(iconf,'(a)')word
      read(iconf,*)boxsiz
      boxsiz=0.5d0*boxsiz
      read(iconf,'(a)')word
      read(iconf,'(a)')word
c
c     read in configuration data
      do i=1,nsp
      read(iconf,'(a)')word
      read(iconf,*)x(i),y(i),z(i)
      x(i)=x(i)/boxsiz
      y(i)=y(i)/boxsiz
      z(i)=z(i)/boxsiz
      read(iconf,*)vx(i),vy(i),vz(i)
      enddo
c
c     calculate fourier exponential terms
      do 10 i=1,nsp
      ekx(i,0)=(1.d0,0.d0)
      eky(i,0)=(1.d0,0.d0)
      ekz(i,0)=(1.d0,0.d0)
      elc(i)=cos(rcl*x(i))+(0.d0,1.d0)*sin(rcl*x(i))
      emc(i)=cos(rcl*y(i))+(0.d0,1.d0)*sin(rcl*y(i))
      enc(i)=cos(rcl*z(i))+(0.d0,1.d0)*sin(rcl*z(i))
   10 continue
      do 20 l=1,kmax
      do 20 i=1,nsp
      ekx(i,l)=ekx(i,l-1)*elc(i)
      eky(i,l)=eky(i,l-1)*emc(i)
      ekz(i,l)=ekz(i,l-1)*enc(i)
   20 continue
c
c     start loop over k vectors
      kkk=1
      mmin=0
      lmin=1
      vxyz(1,kkk)=(0.d0,0.d0)
      vxyz(2,kkk)=(0.d0,0.d0)
      vxyz(3,kkk)=(0.d0,0.d0)
      vxyz(4,kkk)=(0.d0,0.d0)
      do 160 n=0,kmax
      rn=dble(n)
      do 30 i=1,nsp
      enc(i)= ekz(i,n)
   30 continue
      do 150 m=mmin,kmax
      rm=dble(m)
      if(m)60,40,40
   40 do 50 i=1,nsp
      emc(i)= eky(i,m)*enc(i)
   50 continue
      go to 80
   60 do 70 i=1,nsp
      emc(i)= conjg(eky(i,-m))*enc(i)
   70 continue
   80 do 140 l=lmin,kmax
      rl=dble(l)
      if(l)110,90,90
   90 do 100 i=1,nsp
      elc(i)= ekx(i,l)*emc(i)
  100 continue
      go to 130
  110 do 120 i=1,nsp
      elc(i)= conjg(ekx(i,-l))*emc(i)
  120 continue
  130 kkk=kkk+1
      cnm=rnsp/sqrt(rl*rl+rm*rm+rn*rn)
      call crdc(nnn,vx,1,elc,1,argc)
      ckx= cnm*conjg(argc)
      call crdc(nnn,vy,1,elc,1,argc)
      cky= cnm*conjg(argc)
      call crdc(nnn,vz,1,elc,1,argc)
      ckz= cnm*conjg(argc)
      vxyz(1,kkk)=rl*ckx+rm*cky+rn*ckz
      vxyz(2,kkk)=rm*ckz-rn*cky
      vxyz(3,kkk)=rn*ckx-rl*ckz
      vxyz(4,kkk)=rl*cky-rm*ckx
  140 continue
      lmin=-kmax
  150 continue
      mmin=-kmax
  160 continue
c
c     store fourier transform of current density
      write(isave0)vxyz
  170 continue

      close(isave0)
      close(iconf)

      return
      end
      subroutine correl
c***********************************************************************
c
c     calculate current density correlation function
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      complex*16 ckr,ckr0,cfkt
      logical lden,ltim,lcor
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common cfkt(ntime,4*nkt),ckr0(ntime,4*nkt),ckr(4*klim),
     x       ind(ntime),num(ntime)
c
c     set control parameters
      nkt4=4*nkt
      norg=ntime/ngap

      open(isave1,file='CUR_COR',form='unformatted')

c
c     loop over current density components
      ibgn=0
      do 110 kkk=1,klim,nkt

      open(isave0,file='CUR_SPC',form='unformatted')

      lor=0
      mor=0
c
c     initialise cfkt array
      do 10 k=1,nkt4
      do 10 l=1,ntime
      cfkt(l,k)=(0.d0,0.d0)
   10 continue
      do 20 l=1,ntime
   20 num(l)=0
c
c     start of loop over time steps
      do 80 n=0,nsteps-1
      read(isave0)ckr
      if(mod(n,ngap).gt.0)go to 50
      lor=min0(lor+1,norg)
      mor=mod(mor,norg)+1
      ind(mor)=1
      do 40 k=1,nkt4
   40 ckr0(mor,k)=conjg(ckr(k+ibgn))
   50 do 70 l=1,lor
      m=ind(l)
      ind(l)=m+1
      num(m)=num(m)+1
      do 60 k=1,nkt4
   60 cfkt(m,k)=cfkt(m,k)+ckr(k+ibgn)*ckr0(l,k)
   70 continue
   80 continue
c
c     normalise correlation functions
      do 90 k=1,nkt4,4
      rnorm1=real(cfkt(1,k))
      if(rnorm1.gt.0.d0)rnorm1=dble(num(1))/rnorm1
      rnorm2=real(cfkt(1,k+1)+cfkt(1,k+2)+cfkt(1,k+3))
      if(rnorm2.gt.0.d0)rnorm2=dble(num(1))/rnorm2
      do 90 l=1,ntime
      cfkt(l,k)=rnorm1*cfkt(l,k)/dble(num(l))
      cfkt(l,k+1)=rnorm2*cfkt(l,k+1)/dble(num(l))
      cfkt(l,k+2)=rnorm2*cfkt(l,k+2)/dble(num(l))
      cfkt(l,k+3)=rnorm2*cfkt(l,k+3)/dble(num(l))
   90 continue
c
c     store correlation functions in disc file
      write(isave1)cfkt
      ibgn=ibgn+nkt4

      close(isave0)

  110 continue

      close(isave1)

      return
      end
      subroutine curfor
c***********************************************************************
c
c    calculate fourier transform of current density correlation function
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,ltim,lcor
      complex*16 cor,work,fft
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common cor(ntime,4*nkt),work(10*ntime),fft(4*ntime),www(ntime)
      dimension ckw(2*ntime,4*nkt)
      equivalence (cor(1),ckw(1))
      data tpi/6.2831853072d0/,a0,a1,a2/0.42d0,0.50d0,0.08d0/
      iset=klim/nkt
      nkt4=4*nkt
      ntime2=2*ntime
      ntime4=4*ntime

      open(isave2,file='CUR_FFT',form='unformatted')

c
c     initialise cray complex fast fourier transform routine
      ind=1
      isgn=-1
      call cfft2(ind,isgn,ntime4,fft,work,fft)
      ind=0
c
c     set up window function (blackman function)
      arg=tpi/dble(ntime2)
      do 10 i=1,ntime
      ccc=cos(arg*dble(i+ntime-1))
   10 www(i)=a0-a1*ccc+a2*(2.d0*ccc**2-1.d0)
c
c     loop over correlation functions (longitudinal and transverse)
      do 50 kkk=1,klim,nkt

      open(isave1,file='CUR_COR',form='unformatted')

      read(isave1)cor
      do 40 i=1,nkt4
c
c     apply window function
      do 20 j=1,ntime
   20 fft(j)=www(j)*cor(j,i)
      do 25 j=ntime+1,ntime4
   25 fft(j)=(0.d0,0.d0)
      fft(1)=fft(1)/2.d0
c
c     apply complex fourier transform
      call cfft2(ind,isgn,ntime4,fft,work,fft)
c
c     store fourier coefficients
      do 30 j=1,ntime2
   30 ckw(j,i)=real(fft(j))
   40 continue
c
c     save fourier coefficients
      write(isave2)ckw

      close (isave1)

   50 continue

      close (isave2)

      return
      end
      subroutine output
c***********************************************************************
c
c     print out results
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=116000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,ltim,lcor
      complex*16 buff,bufd,bufc,sxyz,slng,strn
      common/cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common/cmcntl/nsteps,lden,lcor,ltim,time,omeg
      common sxyz(4,klim),buff(4*ntime*nkt)
c     note. dimension buff must be 4*max0(ntime*nkt,klim)
      dimension bufd(4,klim),bufc(4*ntime*nkt),bufs(8*ntime*nkt)
      equivalence (buff(1),bufd(1)),(buff(1),bufc(1)),(buff(1),bufs(1))
c
c     set printing control constants
      nkt4=4*nkt
      klim4=4*klim
      ksqm=2*kmax+1
      ntime2=2*ntime
      ntime4=4*ntime
c
c     print out fourier transforms of current density
      if(.not.lden)go to 110
c     calculate mean density

      open(isave0,file='CUR_SPC',form='unformatted')

      do k=1,klim
      sxyz(1,k)=(0.d0,0.d0)
      sxyz(2,k)=(0.d0,0.d0)
      sxyz(3,k)=(0.d0,0.d0)
      sxyz(4,k)=(0.d0,0.d0)
      enddo
      do i=1,nsteps
      read(isave0)bufd
      do k=1,klim
      sxyz(1,k)=sxyz(1,k)+bufd(1,k)
      sxyz(2,k)=sxyz(2,k)+bufd(2,k)
      sxyz(3,k)=sxyz(3,k)+bufd(3,k)
      sxyz(4,k)=sxyz(4,k)+bufd(4,k)
      enddo
      do k=1,klim
      sxyz(1,k)=sxyz(1,k)/dble(nsteps)
      sxyz(2,k)=sxyz(2,k)/dble(nsteps)
      sxyz(3,k)=sxyz(3,k)/dble(nsteps)
      sxyz(4,k)=sxyz(4,k)/dble(nsteps)
      enddo
      enddo

c     write out the density summary

      write(irite,'(a,a,a)')'    L    M    N','  Re(rho(k))  Im(rho(k))'
      k=1
      mmin=0
      lmin=1
      do n=0,kmax
      do m=mmin,kmax
      do l=lmin,kmax
      k=k+1
      write(irite,'(3i5,1p,2e12.4)')l,m,n,sxyz(1,k)
      write(irite,'(3i5,1p,2e12.4)')l,m,n,sxyz(2,k)
      write(irite,'(3i5,1p,2e12.4)')l,m,n,sxyz(3,k)
      write(irite,'(3i5,1p,2e12.4)')l,m,n,sxyz(4,k)
      enddo
      lmin=-kmax
      enddo
      mmin=-kmax
      enddo
      close (isave0)
c
c     print out correlation functions
  110 if(.not.lcor)go to 160
      k0=1
      kkk=1
      do 150 kk=1,klim,nkt

      open(isave1,file='CUR_COR',form='unformatted')

      read(isave1)bufc
      kt=kk+nkt
  120 if(kkk.ge.kt)go to 146
      if(kkk.gt.1)then
      k=kkk+klim-2
      n=k/ksqm**2
      m=k/ksqm-n*ksqm
      l=k-ksqm*(n*ksqm+m)
      n=n-kmax
      m=m-kmax
      l=l-kmax
      write(irite,130)l,m,n
  130 format('current correlation functions for k vector',
     x   3i5,/,8x,'Time',' Re(Cl(k,t))',' Im(Cl(k,t))',
     x  ' Re(Ct(k,t))',' Im(Ct(k,t))')
      ibase=4*(kkk-k0)
      do 140 i=1,ntime
      ttt=dble(i-1)*time
      slng=bufc(ibase*ntime+i)
      strn=bufc((ibase+1)*ntime+i)+bufc((ibase+2)*ntime+i)+
     x  bufc((ibase+3)*ntime+i)
  140 write(irite,145)ttt,slng,strn
  145 format(1p,5e12.4)
      endif
      kkk=kkk+1
      go to 120
  146 k0=k0+nkt
  150 continue
c
c     print out Fourier transforms of correlation functions
  160 if(.not.ltim)go to 210
      k0=1
      kkk=1
      do 200 kk=1,klim,nkt

      open(isave2,file='CUR_FFT',form='unformatted')

      read(isave2)bufs
      kt=kk+nkt
  170 if(kkk.ge.kt)go to 195
      if(kkk.gt.1)then
      k=kkk+klim-2
      n=k/ksqm**2
      m=(k-n*ksqm**2)/ksqm
      l=k-ksqm*(n*ksqm+m)
      n=n-kmax
      m=m-kmax
      l=l-kmax
      write(irite,180)l,m,n
  180 format('Fourier transforms of current correlations for k vector',
     x   3i5,/,7x,'Omega',' FT(Cl(k,t))',' FT(Ct(k,t))')
      ibase=4*(kkk-k0)
      do 185 i=1,ntime2,2
      omg=dble(i-1)*omeg
      ftcl=bufs(ibase*ntime2+i)
      ftct=bufs((ibase+1)*ntime2+i)+bufs((ibase+2)*ntime2+i)+
     x  bufs((ibase+3)*ntime2+i)
  185 write(irite,190)omg,ftcl,ftct
  190 format(1p,3e12.4)
      endif
      kkk=kkk+1
      go to 170
  195 k0=k0+nkt
  200 continue
  210 continue

      close(isave0)
      close(isave1)
      close(isave2)

      return
      end
      subroutine crdc(nnn,aaa,iii,bbb,jjj,arg)
c***********************************************************************
c
c     utility routine for dot product of a real array 'a' and a complex
c     array 'b'.                                 w. smith march 1982.
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      complex*16 bbb(*),arg
      dimension aaa(*)
      i=1
      j=1
      arg=(0.d0,0.d0)
      do k=1,nnn
        arg=arg+aaa(i)*bbb(j)
        i=i+iii
        j=j+jjj
      enddo
      return
      end
      subroutine cfft2(ind,isg,npt,aaa,wrk,bbb)
c***********************************************************************
c
c     fast fourier transform routine
c
c***********************************************************************
      implicit real*8(a-h,o-z)                                          
      complex*16 aaa(npt),bbb(npt),wrk(npt),ttt,ci,cj
      data tpi/6.2831853072d0/
      nu=0
      nt=1
10    nu=nu+1
      nt=2*nt
      if(nt.eq.npt)go to 30
      if(nt.lt.npt)go to 10
      write(6,20)
20    format('error - number of points not a power of two')         
      stop
c
c     initialise complex exponential factors
30    if(ind.eq.0)go to 80
      ci=(0.d0,1.d0)
      cj=(0.d0,-1.d0)
      tpn=tpi/dble(npt)
      arg=0.d0
      np1=npt+1
      np2=npt/2
      wrk(1)=(1.d0,0.d0)
      if(isg)60,40,40
40    do 50 i=1,np2
      arg=arg+tpn
      wrk(i+1)=cdexp(arg*ci)
      wrk(np1-i)=dconjg(wrk(i+1))
50    continue
      return
60    do 70 i=1,np2
      arg=arg+tpn
      wrk(i+1)=cdexp(arg*cj)
      wrk(np1-i)=dconjg(wrk(i+1))
70    continue
      return
c
c     perform fourier transform
80    kkk=0
      nu1=nu-1
      np2=npt/2
      do 85 i=1,npt
85    bbb(i)=aaa(i)
      do 110 l=1,nu
90    do 100 i=1,np2
      iii=0
      jjj=kkk/2**nu1
      do 95 j=1,nu
      jj2=jjj/2
      iii=2*(iii-jj2)+jjj
95    jjj=jj2
      iii=iii+1
      kk1=kkk+1
      k12=kk1+np2
      ttt=bbb(k12)*wrk(iii)
      bbb(k12)=bbb(kk1)-ttt
      bbb(kk1)=bbb(kk1)+ttt
      kkk=kkk+1
100   continue
      kkk=kkk+np2
      if(kkk.lt.npt)go to 90
      kkk=0
      nu1=nu1-1
      np2=np2/2
110   continue
      do 120 kkk=1,npt
      iii=0
      jjj=kkk-1
      do 115 j=1,nu
      jj2=jjj/2
      iii=2*(iii-jj2)+jjj
115   jjj=jj2
      iii=iii+1
      if(iii.le.kkk)go to 120
      ttt=bbb(kkk)
      bbb(kkk)=bbb(iii)
      bbb(iii)=ttt
120   continue
      return
      end
      function  ssum(nnn,ccc,iii)
      implicit real*8(a-h,o-z)
      dimension ccc(*)
      j=1
      ssum=0.d0
      do i=1,nnn
        ssum=ssum+ccc(j)
        j=j+iii
      enddo
      return
      end
      function nbits(n)
c
c**********************************************************************
c
c     dl_poly utility for counting number of set bits in an integer N
c
c     author w smith
c
c     copyright daresbury laboratory 1994
c
c**********************************************************************
c
      m=n
      nbits=0

      do i=1,128

         if(m.eq.0)return

         nbits=nbits+m-2*(m/2)
         m=m/2

      enddo

      return
      end
      function second(time)
      implicit real*8(a-h,o-z)
      save init
      data init/0/
      if(init.eq.0)init=mclock()
      itime=mclock()-init
      time=dble(itime)/1000.d0
      return
      end
