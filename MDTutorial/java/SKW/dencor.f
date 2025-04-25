      subroutine dencor
     x     (lden,lcor,ltim,lchg,fname,mcore,msp,natm,ncon,kmax,
     x     ntime,ngap,numfil,xxx,yyy,zzz,chge,space,status)

c***********************************************************************
c
c     daresbury laboratory dl_poly program for the calculation
c     of the fourier transform of the density in space and time
c
c     original dencor program written by w.smith march 82
c     adapted for dl_poly by w.smith november 1994
c     parallel version 
c
c     author w smith 1994
c     copyright daresbury laboratory 1994
c     adapted by w smith for SDK/CERIUS interface february 1998
c
c***********************************************************************
      implicit real*8(a-h,o-z)

      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)

      integer status,numfil
      character*2 fnum
      character*80 title
      logical lden,lcor,ltim,lchg,kill
      dimension xxx(msp),yyy(msp),zzz(msp),chge(msp)
      dimension man(10),space(mcore)
      data man/10*0/

      status=0
      kill=.false.
      if(numfil.gt.99)then
         fnum="XX"
      else
         write(fnum,'(i2.2)')numfil
      endif

c     determine parallel machine characteristics

      idnode=0
      mxnode=1

c     open report file

      open(3,file="SKWPRG.rep")

      if(idnode.eq.0)then

         write(3,'(a)')
     x   'DL_POLY Density Correlation Program'

      endif

      call start
     x  (idnode,mxnode,mcore,msp,natm,ncon,kmax,ntime,ngap,lden,lcor,
     x   ltim,lchg,kill,tstep)
      if(kill)then

         if(idnode.eq.0)write(3,'(a)')
     x        'Error - terminating job. See above errors'
         status=-1
         close(3)
         return

      endif

      klim=((2*kmax+1)**3-1)/2

c     calculate fourier transform of density
      if(lden)then
        iman=(natm+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+2*iman
        man(3)=man(2)+2*iman
        man(4)=man(3)+2*iman
        man(5)=man(4)+2*klim
        man(6)=man(5)+2*klim
        man(7)=man(6)+2*klim
        man(8)=man(7)+2*iman*(kmax+1)
        man(9)=man(8)+2*iman*(kmax+1)
        man(10)=man(9)+2*iman*(kmax+1)
        call forden
     x     (kill,fname,title,fnum,idnode,mxnode,msp,
     x      natm,ncon,kmax,lchg,tstep,chge,
     x      xxx,yyy,zzz,space(man(1)),space(man(2)),
     x      space(man(3)),space(man(4)),space(man(5)),
     x      space(man(6)),space(man(7)),space(man(8)),
     x      space(man(9)))
        if(kill)then
           
           if(idnode.eq.0)write(3,'(a)')
     x          'Error - rho(k,t) calculation failure'
           status=-2
           close(3)
           return
           
        endif

      endif

c     calculate density correlation function

      if(lcor)then
        nkt=(klim+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+ntime
        man(3)=man(2)+ntime
        man(4)=man(3)+2*nkt*ntime
        man(5)=man(4)+2*nkt*ntime
        man(6)=man(5)+2*klim
        call correl
     x      (fnum,kill,title,idnode,
     x      mxnode,ncon,kmax,ntime,ngap,tstep,space(man(1)),
     x      space(man(2)),space(man(3)),space(man(4)),space(man(5)))
        if(kill)then
           
           if(idnode.eq.0)write(3,'(a)')
     x          'Error - F(k,t) calculation failure'
           status=-3
           close(3)
           return
           
        endif
      endif

c     calculate correlation function fourier transform

      if(ltim)then
        nkt=(klim+mxnode-1)/mxnode
        man(1)=1
        man(2)=man(1)+2*ntime
        man(3)=man(2)+ntime
        man(4)=man(3)+2*nkt*ntime
        man(5)=man(4)+8*ntime
        man(6)=man(5)+8*ntime
        call denfft
     x      (fnum,kill,title,idnode,
     x      mxnode,kmax,ntime,ngap,tstep,space(man(1)),
     x      space(man(2)),space(man(3)),space(man(4)),space(man(5)))
        if(kill)then
           
           if(idnode.eq.0)write(3,'(a)')
     x          'Error - S(k,w) calculation failure'
           status=-4
           close(3)
           return

        endif
      endif

      close(3)
      return
      end
