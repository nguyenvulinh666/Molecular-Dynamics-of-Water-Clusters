      subroutine forden
     x     (kill,fname,title,fnum,idnode,mxnode,msp,natm,ncon,
     x      kmax,lchg,tstep,chge,xxx,yyy,zzz,elc,emc,enc,sxyz,
     x      dxyz,buff,ekx,eky,ekz)
 
c***********************************************************************
c     
c     calculate spatial fourier transform of density
c     
c     author w smith 1994
c     copyright daresbury laboratory 1994
c
c***********************************************************************
      implicit real*8(a-h,o-z)
      parameter (iconf=77,iwork=88)
      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)
      character*2 fnum
      character*8 atname,step
      character*40 fname
      character*80 title
      logical lchg,kill
      complex*16 ekx,eky,ekz,elc,emc,enc,sxyz,dxyz,buff
      dimension xxx(*),yyy(*),zzz(*)
      dimension ekx(*),eky(*),ekz(*),buff(*)
      dimension elc(*),emc(*),enc(*),dxyz(*),sxyz(*)
      dimension cell(9),rcell(9),chge(*)
      data cl/1.d0/
      rcl=tpi/cl
      klim=((2*kmax+1)**3-1)/2

c     zero the particle density accumulator

      do k=1,klim
         sxyz(k)=(0.d0,0.d0)
      enddo

c     open history and density files
      
      open(iconf,file=fname,status='old',err=200)
      open(iwork,file='SPCDEN.'//fnum,form='unformatted')
      
      read(iconf,'(a80)')title
      if(idnode.eq.0)then
         write(3,'(a)')'HISTORY file header:'
         write(3,'(a80)')title
      endif
      read(iconf,'(2i10)')levtrj,imcon

c     check file particulars
      
      if(imcon.eq.0)then
        if(idnode.eq.0)
     x       write(3,'(a)')
     x       'Error - no periodic boundary in HISTORY file'
        close(iconf)
        kill=.true.
        return
      endif
      
      do nstp=1,ncon

        read(iconf,'(a8,4i10,f12.6)',end=100)
     x        step,nstep,natm,ktrj,imcon,timstp

        if(nstp.eq.1)then

           time=0.d0
           nstp0=nstep

        else 
           
           time=timstp*dble(nstep-nstp0)

        endif

        if(natm.gt.msp)then

           if(idnode.eq.0)
     x          write(3,'(a)')
     x          'Error - too many atoms in MD cell'
           close(iconf)
           kill=.true.
           return

        endif

        rnatm=sqrt(1.d0/dble(natm))

c     read in configuration data
        
        read(iconf,*,end=100)cell(1),cell(2),cell(3)
        read(iconf,*,end=100)cell(4),cell(5),cell(6)
        read(iconf,*,end=100)cell(7),cell(8),cell(9)

c     calculate reciprocal lattice basis

        call invert(cell,rcell,det)
        
c     read atomic data

        do i=1,natm
          
          if(lchg)then
             read(iconf,'(a8,i10,2f12.6)')atname,idx,mass,chge(i)
          else
             read(iconf,'(a8)')atname
          endif
          read(iconf,*)xx,yy,zz
          if(levtrj.gt.0)read(iconf,*)vxx,vyy,vzz
          if(levtrj.gt.1)read(iconf,*)fxx,fyy,fzz
          xxx(i)=xx*rcell(1)+yy*rcell(4)+zz*rcell(7)
          yyy(i)=xx*rcell(2)+yy*rcell(5)+zz*rcell(8)
          zzz(i)=xx*rcell(3)+yy*rcell(6)+zz*rcell(9)
          
        enddo

c     calculate fourier exponential terms
        
        m=0
        do i=idnode+1,natm,mxnode
          m=m+1
          ekx(m)=(1.d0,0.d0)
          eky(m)=(1.d0,0.d0)
          ekz(m)=(1.d0,0.d0)
          elc(m)=cmplx(cos(rcl*xxx(i)),sin(rcl*xxx(i)))
          emc(m)=cmplx(cos(rcl*yyy(i)),sin(rcl*yyy(i)))
          enc(m)=cmplx(cos(rcl*zzz(i)),sin(rcl*zzz(i)))
        enddo
        matm=m
        do l=1,kmax
          do i=1,matm
            ekx(i+l*matm)=ekx(i+matm*(l-1))*elc(i)
            eky(i+l*matm)=eky(i+matm*(l-1))*emc(i)
            ekz(i+l*matm)=ekz(i+matm*(l-1))*enc(i)
          enddo
        enddo

c     start loop over k vectors

        kkk=0
        mmin=0
        lmin=1
        do n=0,kmax
          rn=dble(n)
          do i=1,matm
            enc(i)= ekz(i+n*matm)
          enddo
          if(lchg)then
             j=0
             do i=idnode+1,natm,mxnode
                j=j+1
                enc(j)=enc(j)*chge(i)
             enddo
          endif
          do m=mmin,kmax
            rm=dble(m)
            if(m.ge.0)then
              do i=1,matm
                emc(i)= eky(i+m*matm)*enc(i)
              enddo
            else
              do i=1,matm
                emc(i)= conjg(eky(i-m*matm))*enc(i)
              enddo
            endif
            do l=lmin,kmax
              rl=dble(l)
              if(l.ge.0)then
                do i=1,matm
                  elc(i)= ekx(i+l*matm)*emc(i)
                enddo
              else
                do i=1,matm
                  elc(i)= conjg(ekx(i-l*matm))*emc(i)
                enddo
              endif
              kkk=kkk+1
              dxyz(kkk)=(0.d0,0.d0)
              do i=1,matm
                 dxyz(kkk)=dxyz(kkk)+elc(i)
              enddo
            enddo
            lmin=-kmax
          enddo
          mmin=-kmax
        enddo

c     global sum of fourier transform of density
        
        do k=1,klim
           dxyz(k)=rnatm*conjg(dxyz(k))
           sxyz(k)=sxyz(k)+dxyz(k)
        enddo

c     store fourier transform of density

        if(idnode.eq.0)write(iwork)time,(dxyz(j),j=1,klim)

      enddo
     
  100 continue

      tstep=time/dble(nstp-2)
      nstp=nstp-1

c     close working files

      close (iwork)
      close (iconf)

      write(3,'(a)')'Correlation file SPCDEN.'//fnum//' created'

      return
  200 continue

      if(idnode.eq.0)
     x   write(3,'(a)')'Error - HISTORY file not found'
      kill=.true.

      return

      end
