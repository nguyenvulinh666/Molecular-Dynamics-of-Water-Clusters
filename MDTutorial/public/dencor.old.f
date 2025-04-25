c***********************************************************************
c
c     daresbury laboratory ccp5 program for the calculation
c     of correlation functions. this version calculates the
c     fourier transform of the particle density in space and
c     time
c
c***********************************************************************
c
      implicit real*8(a-h,o-z)
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      logical lden,lcor,ldsf
      common space(mcore)
      call second(timrun)
c
c     read in control variables
      call input
      call second(timrun)
      write(irite,10)timrun
   10 format('time after input  = ',f15.8)
c
c     calculate fourier transform of particle density
      if(lden)call forden
      call second(timrun)
      write(irite,20)timrun
   20 format('time after forden = ',f15.8)
c
c     calculate the correlation function
      if(lcor)call corden
      call second(timrun)
      write(irite,30)timrun
   30 format('time after corden = ',f15.8)
c
c     calculate the dynamic structure factor
      if(ldsf)call dynstr
      call second(timrun)
      write(irite,40)timrun
   40 format('time after dynstr = ',f15.8)
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
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      data iread,irite,iconf,isave0,isave1,isave2/5,6,7,8,9,10/
      end
      subroutine input
c***********************************************************************
c
c     read in control variables
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      logical lden,ldsf,lcor,kill
      dimension title(10)
      data twopi/6.2831853072d0/
      kill=.false.
      open(iread,file='den_in')
      read(iread,10)title
   10 format(10a8)
      write(irite,20)title
   20 format(' ccp5 program to  ',/,
     x       ' calculate space  ',/,
     x       ' and time fourier ',/,
     x       ' transforms of    ',/,
     x       ' particle density ',/,
     x       ' w.smith march 82 ',/,
     x       10a8,/)
      read(iread,*)nsteps
      read(iread,*)lden
      read(iread,*)ldsf
      read(iread,*)lcor
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
     x       'number of correlation functions in core   = ',i12,/,
     x       'number of time steps in correlation fn.   = ',i12,/,
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
   38 m1=3*nsp*(kmax+2)+klim
      m2=ntime*(2*nkt+1)+klim
      m3=ntime*(nkt+15)
      m4=ntime*nkt+klim
      kcore=2*max0(m1,m2,m3,m4)
      if(mcore.ge.kcore)go to 50
      write(irite,40)kcore
   40 format('error - insufficient core allocated. ',i8,' words',
     x       ' required')
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
      subroutine forden
c***********************************************************************
c
c     calculate spatial fourier transform of particle density
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      character*8 word
      common /cmcrds/x(nsp),y(nsp),z(nsp)
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common ekx(nsp,0:kmax),eky(nsp,0:kmax),ekz(nsp,0:kmax),
     x       elc(nsp),emc(nsp),enc(nsp),rekr(klim)
      logical lden,ldsf,lcor
      complex*16 ekx,eky,ekz,elc,emc,enc,rekr,argc
      data cl,pi,twopi/2.d0,3.1415926536d0,6.2831853072d0/
      rcl=twopi/cl
      rnsp=sqrt(1.d0/dble(nsp))
      open(iconf,file='HISTORY',status='OLD')
      open(isave0,file='DEN_SPC',form='unformatted')
      read(iconf,'(a)')word
      read(iconf,'(a)')word
      do 180 kstep=1,nsteps
      read(iconf,'(a)')word
      read(iconf,*)boxsiz
      read(iconf,'(a)')word
      read(iconf,'(a)')word
      boxsiz=boxsiz*0.5d0
c
c     read in configuration data
      do i=1,nsp
      read(iconf,'(a)')word
      read(iconf,*)x(i),y(i),z(i)
      x(i)=x(i)/boxsiz
      y(i)=y(i)/boxsiz
      z(i)=z(i)/boxsiz
      read(iconf,*)vx,vy,vz
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
      do 25 l=1,kmax
      do 20 i=1,nsp
      ekx(i,l)=ekx(i,l-1)*elc(i)
      eky(i,l)=eky(i,l-1)*emc(i)
      ekz(i,l)=ekz(i,l-1)*enc(i)
   20 continue
   25 continue
c
c     start loop over k vectors
      kkk=1
      mmin=0
      lmin=1
      rekr(kkk)=(1.d0,0.d0)/rnsp
      do 160 n=0,kmax
      do 30 i=1,nsp
      enc(i)= ekz(i,n)
   30 continue
      do 150 m=mmin,kmax
      if(m)60,40,40
   40 do 50 i=1,nsp
      emc(i)= eky(i,m)*enc(i)
   50 continue
      go to 80
   60 do 70 i=1,nsp
      emc(i)= conjg(eky(i,-m))*enc(i)
   70 continue
   80 do 140 l=lmin,kmax
      if(l)110,90,90
   90 do 100 i=1,nsp
      elc(i)= ekx(i,l)*emc(i)
  100 continue
      go to 130
  110 do 120 i=1,nsp
      elc(i)= conjg(ekx(i,-l))*emc(i)
  120 continue
  130 kkk=kkk+1
      call csum(nsp,elc,1,argc)
      rekr(kkk)=rnsp*conjg(argc)
  140 continue
      lmin=-kmax
  150 continue
      mmin=-kmax
  160 continue
c
c     store fourier transform of density
      write(isave0)rekr
  180 continue
      close (iconf)
      close (isave0)
      return
      end
      subroutine corden
c***********************************************************************
c
c     calculate density correlation function
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,ldsf,lcor
      complex*16 ckr,ckr0,cfkt
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common cfkt(ntime,nkt),ckr0(ntime,nkt),ckr(klim),
     x       ind(ntime),num(ntime)
      open(isave1,file='DEN_COR',form='unformatted')
c
c     calculate working parameters
      norg=ntime/ngap
c
c     loop over components of transform of density
      ibgn=0
      do 70 kkk=1,klim,nkt

      open(isave0,file='DEN_SPC',form='unformatted')

      lor=0
      mor=0
c
c     initialise cfkt array
      do 10 k=1,nkt
      do 10 l=1,ntime
      cfkt(l,k)=(0.d0,0.d0)
   10 continue
      do 15 l=1,ntime
   15 num(l)=0
c
c     loop over time steps
      do 50 n=0,nsteps-1
      read(isave0)ckr
      if(mod(n,ngap).gt.0)go to 30
      lor=min0(lor+1,norg)
      mor=mod(mor,norg)+1
      ind(mor)=1
      do 20 k=1,nkt
   20 ckr0(mor,k)=conjg(ckr(k+ibgn))
   30 do 45 l=1,lor
      m=ind(l)
      ind(l)=m+1
      num(m)=num(m)+1
      do 40 k=1,nkt
   40 cfkt(m,k)=cfkt(m,k)+ckr(k+ibgn)*ckr0(l,k)
   45 continue
   50 continue
c
c     normalise correlation function
      do 60 k=1,nkt
      rnorm=dble(num(1))/real(cfkt(1,k))
      do 60 l=1,ntime
      cfkt(l,k)=rnorm*cfkt(l,k)/dble(num(l))
   60 continue
c
c     write out correlation functions to disc
      write(isave1)cfkt
      ibgn=ibgn+nkt

      close(isave0)

   70 continue

      close (isave1)

      return
      end
      subroutine dynstr
c***********************************************************************
c
c     calculate dynamic structure factor
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,ldsf,lcor
      complex*16 cor,work,fft
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common cor(ntime,nkt),work(10*ntime),fft(4*ntime),www(ntime)
      dimension skw(2*ntime,nkt)
      equivalence (cor(1),skw(1))
      data tpi/6.2831853072d0/,a0,a1,a2/0.42d0,0.50d0,0.08d0/
      ntime2=2*ntime
      ntime4=4*ntime

      open(isave1,file='DEN_COR',form='unformatted')
      open(isave2,file='DEN_FFT',form='unformatted')

c
c     initialise cray fast fourier transform routine
      indx=1
      isgn=-1
      call cfft2(indx,isgn,ntime4,fft,work,fft)
      indx=0
c
c     set up the window function (blackman function)
      arg=tpi/dble(ntime2)
      do 10 i=1,ntime
      ccc=cos(arg*dble(i+ntime-1))
   10 www(i)=a0-a1*ccc+a2*(2.d0*ccc**2-1.d0)

c
c     loop over correlation functions for all rho(k) terms
      do 50 kkk=1,klim,nkt
      read(isave1)cor
c
c     calculate dynamic structure factors
      do 40 i=1,nkt
c
c     apply window function for fourier transform
      do 20 j=1,ntime
   20 fft(j)=www(j)*cor(j,i)
      do 25 j=ntime+1,ntime4
   25 fft(j)=(0.d0,0.d0)
      fft(1)=fft(1)/2.d0
c
c     apply complex fourier transform
      call cfft2(indx,isgn,ntime4,fft,work,fft)
c
c     store dynamic structure factor
      do 30 j=1,ntime2
   30 skw(j,i)=real(fft(j))
c
c     calculate area under skw function
      skw(ntime2,i)=omeg*(2.d0*ssum(ntime/2,skw(1,i),2)-skw(1,i)
     x                   +4.d0*ssum(ntime/2,skw(2,i),2)-skw(ntime-1,i)
     x                   -4.d0*skw(ntime,i))/3.d0
   40 continue
c
c     save dynamic structure factors
      write(isave2)skw
   50 continue

      close (isave1)
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
      parameter (mcore=30000,nsp=216,ntime=128)
      parameter (kmax=1,klim=14,nkt=14,ngap=1)
      logical lden,ldsf,lcor
      complex*16 buff,bufd,bufc,sxyz
      common /cmchnl/iread,irite,iconf,isave0,isave1,isave2
      common /cmcntl/nsteps,lden,ldsf,lcor,time,omeg
      common sxyz(klim),buff(ntime*nkt)
c     note. dimension buff must be max0(ntime*nkt,klim)
      dimension bufd(klim),bufc(ntime*nkt),bufs(ntime*nkt*2)
      equivalence (buff(1),bufd(1)),(buff(1),bufc(1)),(buff(1),bufs(1))

c
c     set printing control constants
      ksqm=2*kmax+1
      ntime2=2*ntime

c
c     print out fourier transforms of particle density
      if(.not.lden)go to 110
c     calculate mean density

      open(isave0,file='DEN_SPC',form='unformatted')

      do k=1,klim
      sxyz(k)=(0.d0,0.d0)
      enddo
      do i=1,nsteps
      read(isave0)bufd
      do k=1,klim
      sxyz(k)=sxyz(k)+bufd(k)
      enddo
      enddo
      do k=1,klim
      sxyz(k)=sxyz(k)/dble(nsteps)
      enddo

c     write out the density summary

      write(irite,'(a,a,a)')'    L    M    N','        S(k)',
     x     '  Re(rho(k))  Im(rho(k))'
      k=1
      mmin=0
      lmin=1
      do n=0,kmax
      do m=mmin,kmax
      do l=lmin,kmax
      k=k+1
      sok=sxyz(k)*conjg(sxyz(k))
      write(irite,'(3i5,1p,3e12.4)')l,m,n,sok,sxyz(k)
      enddo
      lmin=-kmax
      enddo
      mmin=-kmax
      enddo
      close (isave0)
c
c     print out intermediate scattering function
  110 if(.not.lcor)go to 160
      k0=1
      kkk=1
      do 150 kk=1,klim,nkt

      open(isave1,file='DEN_COR',form='unformatted')

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
  130 format('intermediate scattering function for k vector',
     x       3i5)
      do 140 i=1,ntime
  140 write(irite,145)dble(i-1)*time,bufc((kkk-k0)*ntime+i)
  145 format(1p,3e12.4)
      endif
      kkk=kkk+1
      go to 120
  146 k0=k0+nkt
  150 continue
c
c     print out dynamic structure factors
  160 if(.not.ldsf)go to 210
      k0=1
      kkk=1
      do 200 kk=1,klim,nkt

      open(isave2,file='DEN_FFT',form='unformatted')

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
      write(irite,180)l,m,n,bufs((kkk-k0+1)*ntime2)
  180 format('dynamic structure factor for k vector ',
     x       3i5,/,'(area under curve =', 1pe12.4,')')
      do 185 i=1,ntime2,2
  185 write(irite,190)dble(i-1)*omeg,bufs((kkk-k0)*ntime2+i)
  190 format(1p,2e12.4)
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
      subroutine  csum(nnn,ccc,iii,arg)
      implicit real*8(a-h,o-z)
      complex*16 ccc(*),arg
      j=1
      arg=(0.d0,0.d0)
      do i=1,nnn
        arg=arg+ccc(j)
        j=j+iii
      enddo
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
