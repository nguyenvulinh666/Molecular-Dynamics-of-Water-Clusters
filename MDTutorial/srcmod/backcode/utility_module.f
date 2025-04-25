      module utility_module

c***********************************************************************
c     
c     dl_poly module for utility subroutines and functions
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************
      
      implicit none

      contains
      
      subroutine global_sum_forces(natms,mxnode,gxx,gyy,gzz)
      
c***********************************************************************
c     
c     dl_poly subroutine to perform global sum of atomic forces as
c     requred by replicated data strategy
c     
c     copyright - daresbury laboratory 
c     author    - w.smith december 2005
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************
      
      use config_module
      
      implicit none
      
      integer natms,mxnode,i,j
      real(8) gxx(*),gyy(*),gzz(*)

      if(mxnode.gt.1) then
CVAM
CVAM              call VTBEGIN(28, ierr)
CVAM
        j=0
        do i=1,natms
          
          buffer(j+1)=gxx(i)
          buffer(j+2)=gyy(i)
          buffer(j+3)=gzz(i)
          j=j+3
          
        enddo
        call gdsum(buffer(1),3*natms,buffer(3*natms+1))
        j=0
        do i=1,natms
          
          gxx(i)=buffer(j+1)
          gyy(i)=buffer(j+2)
          gzz(i)=buffer(j+3)
          j=j+3
          
        enddo
CVAM
CVAM              call VTEND(28, ierr)
CVAM
      endif
      
      return
      end subroutine global_sum_forces
      
      subroutine images
     x  (imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating the minimum image
c     of atom pairs within a specified MD cell
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith march 1992.
c     T3D optimised version. t.forester july 1994
c     
c     for
c     imcon=0 no boundary conditions apply
c     imcon=1 standard cubic boundaries apply
c     imcon=2 orthorhombic boundaries apply
c     imcon=3 parallelepiped boundaries apply
c     imcon=4 truncated octahedron boundaries apply
c     imcon=5 rhombic dodecahedron boundaries apply
c     imcon=6 x-y parallelogram boundary conditions : no periodicity in z
c     imcon=7 hexagonal prism boundaries apply
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************
      
      use error_module
      
      implicit none

      integer imcon,idnode,mxnode,natms,iatm1,iatm2,i
      real(8) cell,xxx,yyy,zzz,aaa,bbb,ccc,det,rt2,rt3,ssx
      real(8) ssy,ssz,ddd,xss,yss,zss,rcell

      dimension xxx(*),yyy(*),zzz(*)
      dimension cell(9),rcell(9)

      data rt2/1.41421356623d0/,rt3/1.7320508075d0/

CVAM
CVAM      call VTBEGIN(81, ierr)
CVAM

      if(imcon.gt.0) then

c     block indices

        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif
      
      if(imcon.eq.1)then

c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)

        do i=iatm1,iatm2
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
        enddo
        
      else if(imcon.eq.2)then

c     rectangular (slab) boundary conditions
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(5)
        ccc=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(5)*nint(bbb*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ccc*zzz(i))
          
        enddo
        
      else if(imcon.eq.3)then

c     parallelepiped boundary conditions
        
        call invert(cell,rcell,det)
        
        do i=iatm1,iatm2
          
          ssx=(rcell(1)*xxx(i)+rcell(4)*yyy(i)+rcell(7)*zzz(i))
          ssy=(rcell(2)*xxx(i)+rcell(5)*yyy(i)+rcell(8)*zzz(i))
          ssz=(rcell(3)*xxx(i)+rcell(6)*yyy(i)+rcell(9)*zzz(i))
          
          xss=ssx-nint(ssx)
          yss=ssy-nint(ssy)
          zss=ssz-nint(ssz)
          
          xxx(i)=(cell(1)*xss+cell(4)*yss+cell(7)*zss)
          yyy(i)=(cell(2)*xss+cell(5)*yss+cell(8)*zss)
          zzz(i)=(cell(3)*xss+cell(6)*yss+cell(9)*zss)
          
        enddo
        
      else if(imcon.eq.4)then

c     truncated octahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(5)-cell(9)).lt.1.d-6)) call error(idnode,130)
        
        aaa=1.d0/cell(1)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(zzz(i))).ge.
     x      (0.75d0*cell(1)))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(1),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.5)then

c     rhombic dodecahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6)) 
     x    call error(idnode,140)
        
        aaa=1.d0/cell(1)
        bbb=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(bbb*zzz(i))
          
          if((abs(xxx(i))+abs(yyy(i))+abs(rt2*zzz(i))).ge.
     x      cell(1))then
            
            xxx(i)=xxx(i)-0.5d0*sign(cell(1),xxx(i))
            yyy(i)=yyy(i)-0.5d0*sign(cell(1),yyy(i))
            zzz(i)=zzz(i)-0.5d0*sign(cell(9),zzz(i))
            
          endif
          
        enddo
        
      else if(imcon.eq.6) then

c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

        if(abs(det).lt.1.d-6)call error(idnode,120)
        
        det = 1.d0/det

        rcell(1) =  det*cell(5)
        rcell(2) = -det*cell(2)
        rcell(4) = -det*cell(4)
        rcell(5) =  det*cell(1)
        
        do i=iatm1,iatm2

          ssx = rcell(1)*xxx(i) + rcell(4)*yyy(i)
          ssy = rcell(2)*xxx(i) + rcell(5)*yyy(i)

          xss = ssx - nint(ssx)
          yss = ssy - nint(ssy)

          xxx(i)=cell(1)*xss + cell(4)*yss
          yyy(i)=cell(2)*xss + cell(5)*yss

        enddo

      else if(imcon.eq.7) then

c     hexagonal prism boundary conditions
        
        if(abs(cell(1)-rt3*cell(5)).ge.1.d-6)
     x    call error(idnode,135)
        
        aaa=cell(1)/(rt3*2.d0)
        bbb=cell(1)/rt3
        ccc=rt3/cell(1)
        ddd=1.d0/cell(9)
        
        do i=iatm1,iatm2
          
          yyy(i)=yyy(i)-bbb*nint(ccc*yyy(i))
          zzz(i)=zzz(i)-cell(9)*nint(ddd*zzz(i))
          
          if((abs(yyy(i))+abs(rt3*xxx(i))).ge.bbb)then
            
            xxx(i)=xxx(i)-rt3*sign(aaa,xxx(i))
            yyy(i)=yyy(i)-sign(aaa,yyy(i))
            
          endif
          
        enddo
        
      endif
CVAM
CVAM      call VTEND(81, ierr)
CVAM
      return
      end subroutine images

      subroutine config_write(fname,levcfg,imcon,natms,engcfg)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing CONFIG files
c     
c     copyright - daresbury laboratory 
c     author    - w. smith aug 2007
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************
      
      use config_module
      use setup_module
      
      implicit none
      
      character*6 fname
      
      integer i,natms,levcfg,imcon,nstep
      real(8) engcfg

      open(nconf,file=fname,form='formatted')
      
      write(nconf,'(80a1)') cfgname
      write(nconf,'(3i10,1p,g20.12)') levcfg,imcon,natms,engcfg
      if(imcon.gt.0) write(nconf,'(3f20.12)') cell
      
      do i=1,natms
        
        write(nconf,'(a8,i10)') atmnam(i),i
        write(nconf,'(3g20.10)') xxx(i),yyy(i),zzz(i)
        if(levcfg.gt.0)write(nconf,'(3g20.12)')
     x    vxx(i),vyy(i),vzz(i)
        if(levcfg.gt.1)write(nconf,'(3g20.12)') 
     x    fxx(i),fyy(i),fzz(i)
        
      enddo
      
      close (nconf)
      
      return
      end subroutine config_write
      
      subroutine bomb(idnode,nyr,nmn,ndy)

c***********************************************************************
c     
c     dl_poly subroutine to set an expiry date in a compiled program
c     
c     copyright - daresbury laboratory 
c     author    - w. smith    oct 2002
c
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      use setup_module

      implicit none

      logical safe
      integer info(8)
      character*12 dat,tim,zon
      integer idnode,nyr,nmn,ndy

      safe=.true.

      call date_and_time(dat,tim,zon,info)
      
      if(info(1).gt.nyr)then

        safe=.false.

      else if(info(1).eq.nyr)then

        if(info(2).gt.nmn)then

          safe=.false.

        else if(info(2).eq.nmn)then

          if(info(3).ge.ndy)safe=.false.

        endif

      endif

      if(.not.safe)then

        if(idnode.eq.0)write(nrite,'(a,/,a)')
     x    'THE EXPIRY DATE OF THIS EXECUTABLE HAS PASSED.',
     X    'PLEASE CONTACT W.SMITH@DL.AC.UK FOR A NEW LICENCE'

        call exitcomms()

      endif

      return
      end subroutine bomb

      subroutine cpy_rtc(nnn,aaa,bbb)

c**********************************************************************
c
c     dl_poly subroutine for copying a real array into a complex array
c     of the same dimension
c
c     copyright daresbury laboratory 1998
c     author w.smith oct 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) aaa(*)
      complex(8) bbb(*)

      do i=1,nnn

        bbb(i)=cmplx(aaa(i),0.d0,kind=8)

      enddo

      return
      end subroutine cpy_rtc

      function duni()

c*********************************************************************
c     
c     dl_poly random number generator based on the universal
c     random number generator of marsaglia, zaman and tsang
c     (stats and prob. lett. 8 (1990) 35-39.) it must be
c     called once to initialise parameters u,c,cd,cm
c     
c     copyright daresbury laboratory 1992
c     author -  w.smith         july 1992
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c*********************************************************************

      implicit none

      logical new
      integer ir,jr,i,j,k,l,m,ii,jj
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

CVAM
CVAM      call VTBEGIN(134, ierr)
CVAM

      if(new)then

c     initial values of i,j,k must be in range 1 to 178 (not all 1)
c     initial value of l must be in range 0 to 168.

        i=12
        j=34
        k=56
        l=78
c     
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

c     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
CVAM
CVAM      call VTEND(134, ierr)
CVAM
      return
      end function duni

      subroutine ele_prd(nnn,aaa,bbb,ccc)

c**********************************************************************
c
c     dl_poly subroutine for element by element product of
c     a real array (bbb) and a complex array (ccc)
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) bbb(*)
      complex(8) aaa(*),ccc(*)

      do i=1,nnn

        aaa(i)=bbb(i)*ccc(i)

      enddo

      return
      end subroutine ele_prd

      subroutine gauss(natms,vxx,vyy,vzz,idnode,mxnode)

c*********************************************************************
c     
c     dl_poly subroutine for constructing velocity arrays
c     with a gaussian distribution of unit variance.
c     
c     based on the Box-Muller method
c     
c     note - this version uses a universal random number 
c     generator, which generates pseudo-random numbers between
c     0 and 1. it is based on the algorithm of marsaglia, zaman
c     and tsang in: stats and prob. lett. 8 (1990) 35-39.
c     
c     copyright daresbury laboratory 2007
c     author - w. smith         nov  2007
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c*********************************************************************
      use random_module

      use setup_module
      
      implicit none

      integer natms,ierr,i, idnode, mxnode
      real*8 vxx,vyy,vzz,a1,a3,a5,a7,a9,rr1,rr2
      integer iatm1, iatm2
      
      dimension vxx(natms),vyy(natms),vzz(natms)
      
      data a1,a3,a5/3.949846138d0,0.252408784d0,0.076542912d0/
      data a7,a9/0.008355968d0,0.029899776d0/
      
CVAM
CVAM      call VTBEGIN(144, ierr)
CVAM

      iatm1=(idnode*natms)/mxnode+1
      iatm2=((idnode+1)*natms)/mxnode

!!p      do i=1,natms
      do i=iatm1,iatm2
        
        rr1=grnd()
        rr2=grnd()
        vxx(i) = SQRT(-2*LOG(1-rr1))*COS(2.*pi*rr2)

        rr1=grnd()
        rr2=grnd()
        vyy(i) = SQRT(-2*LOG(1-rr1))*COS(2.*pi*rr2)
 
        rr1=grnd()
        rr2=grnd()
        vzz(i) = SQRT(-2*LOG(1-rr1))*COS(2.*pi*rr2)
        
      enddo
CVAM
CVAM      call VTEND(144, ierr)
CVAM

      return

c$$$      integer natms,i
c$$$      real(8) vxx,vyy,vzz,rrr,rr1,rr2
c$$$      
c$$$      dimension vxx(natms),vyy(natms),vzz(natms)
c$$$CVAM
c$$$CVAM      call VTBEGIN(144, ierr)
c$$$CVAM
c$$$      
c$$$c     initialise random number generator
c$$$      
c$$$      rrr=duni()
c$$$      
c$$$c     calculate gaussian random numbers
c$$$      
c$$$      do i=1,2*(natms/2),2
c$$$        
c$$$        rr1=sqrt(-2.d0*log(duni()))
c$$$        rr2=2.d0*pi*duni()
c$$$        vxx(i)=rr1*cos(rr2)
c$$$        vxx(i+1)=rr1*sin(rr2)
c$$$
c$$$        rr1=sqrt(-2.d0*log(duni()))
c$$$        rr2=2.d0*pi*duni()
c$$$        vyy(i)=rr1*cos(rr2)
c$$$        vyy(i+1)=rr1*sin(rr2)
c$$$
c$$$        rr1=sqrt(-2.d0*log(duni()))
c$$$        rr2=2.d0*pi*duni()
c$$$        vzz(i)=rr1*cos(rr2)
c$$$        vzz(i+1)=rr1*sin(rr2)
c$$$        
c$$$      enddo
c$$$      if(mod(natms,2).ne.0)then
c$$$        
c$$$        rr1=sqrt(-2.d0*log(duni()))
c$$$        rr2=2.d0*pi*duni()
c$$$        vxx(natms)=rr1*cos(rr2)
c$$$        vyy(natms)=rr1*sin(rr2)
c$$$        rr1=sqrt(-2.d0*log(duni()))
c$$$        rr2=2.d0*pi*duni()
c$$$        vzz(natms)=rr1*cos(rr2)
c$$$        
c$$$      endif
c$$$CVAM
c$$$CVAM      call VTEND(144, ierr)
c$$$CVAM
      end subroutine gauss

      
      double precision function gasdev()
c implemented as described in numerical recipes
      use random_module

      implicit none
      integer, save :: iset = 0
      real*8, save :: gset
      real*8 fac,rsq,v1,v2
      if(iset==0) then
1       v1=2.*grnd()-1.0d0
        v2=2.*grnd()-1.0d0
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      end if
      end function

      double precision function gamdev(ia)
c implemented as described in numerical recipes
      use random_module

      implicit none
      integer, intent(in) :: ia
      integer j
      real*8 am,e,s,v1,v2,x,y
      if(ia.lt.1)pause 'bad argument in gamdev'
      if(ia.lt.6)then
        x=1.
        do 11 j=1,ia
          x=x*grnd()
11      continue
        x=-log(x)
      else
1         v1=2.*grnd()-1.
          v2=2.*grnd()-1.
        if(v1**2+v2**2.gt.1.)goto 1
          y=v2/v1
          am=ia-1
          s=sqrt(2.*am+1.)
          x=s*y+am
        if(x.le.0.)goto 1
          e=(1.+y**2)*exp(am*log(x/am)-s*y)
        if(grnd().gt.e)goto 1
      endif
      gamdev=x
      end function

      double precision function betadev(ip,iq)
        implicit none
        integer, intent(in) :: ip,iq
c returns a number with distribution probability
c P(x) propto x**(ip-1) * (1.0-x)**(iq-1)
c useful for Tsallis distribution
c
c implemented as described in Schmeiser Babu, Operations Research, 28, 917 (1980)
c (in fact it does not use their method, but a method that they call RG.)
        real*8 :: y,z
        y = gamdev(ip)
        z = gamdev(iq)
        betadev=y/(y+z)
      end function betadev

      double precision function sumnoises(nn)
        implicit none
        integer, intent(in) :: nn
c returns the sum of n independent gaussian noises squared
c (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
        if(modulo(nn,2)==0) then
          sumnoises=2.0*gamdev(nn/2)
        else
          sumnoises=2.0*gamdev((nn-1)/2) + gasdev()**2
        end if
      end function sumnoises

      subroutine invert(a,b,d)

c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      real(8) a,b,d,r

      dimension a(9),b(9)

c     calculate adjoint matrix
      b(1)=a(5)*a(9)-a(6)*a(8)
      b(2)=a(3)*a(8)-a(2)*a(9)
      b(3)=a(2)*a(6)-a(3)*a(5)
      b(4)=a(6)*a(7)-a(4)*a(9)
      b(5)=a(1)*a(9)-a(3)*a(7)
      b(6)=a(3)*a(4)-a(1)*a(6)
      b(7)=a(4)*a(8)-a(5)*a(7)
      b(8)=a(2)*a(7)-a(1)*a(8)
      b(9)=a(1)*a(5)-a(2)*a(4)

c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d

c     complete inverse matrix
      b(1)=r*b(1)
      b(2)=r*b(2)
      b(3)=r*b(3)
      b(4)=r*b(4)
      b(5)=r*b(5)
      b(6)=r*b(6)
      b(7)=r*b(7)
      b(8)=r*b(8)
      b(9)=r*b(9)

      return
      end subroutine invert

      subroutine jacobi(a,v,n)

c***********************************************************************
c     
c     diagonalisation of real symmetric matices by jacobi method
c     
c     input parameters:
c     
c     a(n,n) is the matrix to be diagonalised
c     v(n,n) is the eigenvector matrix
c     n   is the dimension of the matrices
c     
c     jacobi processes lower triangle only (upper triangle unchanged)
c     
c     variable rho sets absolute tolerance on convergence
c     variable tes is a moving tolerance that diminishes
c     on each pass until at true convergence tes<rho
c     
c     author w.smith 1993
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      logical pass
      integer n,i,j,k
      real(8) a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        enddo
        v(i,i)=1.0d0
      enddo

c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl)scl=abs(a(i,i))
      enddo
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        enddo
      enddo

c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho
        
c     jacobi diagonalisation
        
        pass=.true.
        
c     recycle until moving tolerance satisfied
        
        do while(pass)
          
          pass=.false.
          
          do i=2,n
            
            do j=1,i-1
              
              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif
              
            enddo
            
          enddo
          
        enddo

      enddo

c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      subroutine scl_csum(nnn,tot,aaa)

c**********************************************************************
c
c     dl_poly subroutine to calculate the scalar sum of the elements
c     of a complex array
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      complex(8) aaa(*),tot

      tot=(0.d0,0.d0)

      do i=1,nnn

        tot=tot+aaa(i)

      enddo

      return
      end subroutine scl_csum

      subroutine set_block(nnn,ccc,aaa)

c**********************************************************************
c
c     dl_poly subroutine to initialise an array to a single value
c
c     copyright daresbury laboratory 1998
c     author w.smith july 1998
c
c**********************************************************************

      implicit none

      integer i,nnn
      real(8) ccc,aaa(nnn)

      do i=1,nnn,2

        aaa(i)=ccc
        aaa(i+1)=ccc

      enddo
      
      return
      end subroutine set_block

      subroutine shellsort(n,list)

c***********************************************************************
c     
c     dlpoly shell sort routine. 
c     Sorts an array of integers into ascending order
c     
c     copyright daresbury laboratory 1993
c     author - t.forester   november 1993
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      implicit none

      integer n,list,nn,nl,i,j,ix,imax

      dimension list(*)

c     set up sort

      if(n.gt.1) then

c     number of lists

        nl = n/2

c     iterate shell sort

        do while(nl.gt.0)

          do nn = 1,nl
            
c     begin insertion sort on nnth list
            
            do i = nn+nl,n,nl
              
              imax = list(i)
              ix = i
              
c     find location for insertion
              
              j = i
              do while(j.ge.nl+1)
                
                j = j-nl
                if (list(j).gt.imax) then
                  ix = j
                else
                  j =1
                endif
                
              enddo
              
c     insert in index array
              
              do j = i,ix+nl,-nl
                list(j) = list(j-nl)
              enddo
              
              list(ix) = imax
              
            enddo
            
          enddo
        
          nl = nl/2

        enddo
        
      endif

      return
      end subroutine shellsort

      subroutine fcap(lfcap,natms,fmax,temp)
      
c*********************************************************************
c     
c     DLPOLY routinue for limiting the absolute magnitude of
c     forces. Used in equilibration period only
c     
c     copyright daresbury laboratory 1993
c     
c     author -     t. forester march 1993
c     amended-     t. forester  sept 1994
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c*********************************************************************

      use config_module
      
      implicit none

      logical lfcap
      integer natms,i
      real(8) fmax,temp,fmax1,fmax2,fxc,fyc,fzc,fmod,fscale
      
CVAM
CVAM      call VTBEGIN(29, ierr)
CVAM

      if(lfcap) then

c     maximum force permitted
        
        fmax1 = boltz*fmax*temp
        fmax2 = fmax1*fmax1

c     cap forces and conserve linear momentum
        
        fxc = 0.d0
        fyc = 0.d0
        fzc = 0.d0
        
        do i = 1,natms
          
          fmod = fxx(i)**2 + fyy(i)**2 + fzz(i)**2
          
          if(fmod.gt.fmax2) then
            
            fscale = sqrt(fmax2/fmod)
            
            fxx(i) = fxx(i)*fscale
            fyy(i) = fyy(i)*fscale
            fzz(i) = fzz(i)*fscale
            
          endif

c     accummulate forces - to check on momentum conservation
          
          fxc = fxc + fxx(i)
          fyc = fyc + fyy(i)
          fzc = fzc + fzz(i)
          
        enddo

c     ensure net forces sum to zero
        
        fxc = -fxc/dble(natms)
        fyc = -fyc/dble(natms)
        fzc = -fzc/dble(natms)

c     conserve momentum
        
        do i = 1,natms
          
          fxx(i) = fxx(i) + fxc
          fyy(i) = fyy(i) + fyc
          fzz(i) = fzz(i) + fzc
          
        enddo
        
      endif
CVAM
CVAM      call VTEND(29, ierr)
CVAM
      return
      end subroutine fcap

      subroutine freeze(natms)

c***********************************************************************
c     
c     dlpoly routine to quench forces and velocities on 'frozen' atoms
c     replicated data version - blocked data
c     
c     copyright daresbury laboratory 1994
c     author t.forester nov 1994
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      use config_module

      implicit none

      integer natms,i
CVAM
CVAM      call VTBEGIN(31, ierr)
CVAM

      do i = 1,natms
        
        if(lstfrz(i).ne.0) then
          
          vxx(i) = 0.d0
          vyy(i) = 0.d0
          vzz(i) = 0.d0
          fxx(i) = 0.d0
          fyy(i) = 0.d0
          fzz(i) = 0.d0
          
        endif
        
      enddo
CVAM
CVAM      call VTEND(31, ierr)
CVAM
      return
      end subroutine freeze

      subroutine matmul(aaa,bbb,ccc)

c***********************************************************************
c     
c     dlpoly utility to multiply 3x3 matrices
c
c     copyright daresbury laboratory
c     author      w.smith  oct  2005
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c**********************************************************************

      implicit none

      integer i
      real(8) aaa(9),bbb(9),ccc(9),tmp(9)

      tmp(1)=aaa(1)*bbb(1)+aaa(4)*bbb(2)+aaa(7)*bbb(3)
      tmp(2)=aaa(2)*bbb(1)+aaa(5)*bbb(2)+aaa(8)*bbb(3)
      tmp(3)=aaa(3)*bbb(1)+aaa(6)*bbb(2)+aaa(9)*bbb(3)

      tmp(4)=aaa(1)*bbb(4)+aaa(4)*bbb(5)+aaa(7)*bbb(6)
      tmp(5)=aaa(2)*bbb(4)+aaa(5)*bbb(5)+aaa(8)*bbb(6)
      tmp(6)=aaa(3)*bbb(4)+aaa(6)*bbb(5)+aaa(9)*bbb(6)

      tmp(7)=aaa(1)*bbb(7)+aaa(4)*bbb(8)+aaa(7)*bbb(9)
      tmp(8)=aaa(2)*bbb(7)+aaa(5)*bbb(8)+aaa(8)*bbb(9)
      tmp(9)=aaa(3)*bbb(7)+aaa(6)*bbb(8)+aaa(9)*bbb(9)
      
      do i=1,9
        ccc(i)=tmp(i)
      enddo
      
      return
      end subroutine matmul

      subroutine getrotmat(q0,q1,q2,q3,rot)
      
c***********************************************************************
c     
c     dlpoly utility to  construct rotation matrix
c     from quaternions using x convention for euler angles
c
c     copyright daresbury laboratory
c     author      w.smith   mar 2005
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c**********************************************************************

      implicit none
      
      real(8) q0,q1,q2,q3,rot(9)
      
      rot(1)=q0**2+q1**2-q2**2-q3**2
      rot(2)=2.d0*(q1*q2-q0*q3)
      rot(3)=2.d0*(q1*q3+q0*q2)
      rot(4)=2.d0*(q1*q2+q0*q3)
      rot(5)=q0**2-q1**2+q2**2-q3**2
      rot(6)=2.d0*(q2*q3-q0*q1)
      rot(7)=2.d0*(q1*q3-q0*q2)
      rot(8)=2.d0*(q2*q3+q0*q1)
      rot(9)=q0**2-q1**2-q2**2+q3**2
      
      return
      end subroutine getrotmat

      function sdot0(n,aaa,bbb)

c***********************************************************************
c     
c     dlpoly utility to calculate scalar product of two arrays
c
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c**********************************************************************

      implicit none

      integer n,i
      real(8) sdot0,aaa,bbb

      dimension aaa(*),bbb(*)

      sdot0=0.d0

      do i=1,n
        sdot0=sdot0+aaa(i)*bbb(i)
      enddo

      return
      end function sdot0

      function sdot1(natms,idnode,mxnode,aaa,bbb)

c***********************************************************************
c     
c     dlpoly utility to calculate scalar product of two arrays
c     distributed version
c     
c     copyright daresbury laboratory
c     author      w.smith  july 2005
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c**********************************************************************
      
      use config_module
      
      implicit none
      
      integer natms,idnode,mxnode,i,iatm0,iatm1
      real(8) sdot1,aaa,bbb
      
      dimension aaa(*),bbb(*)
      
c     assign block of atoms to processor

      iatm0=(idnode*natms)/mxnode + 1
      iatm1=((idnode+1)*natms)/mxnode
      
      sdot1=0.d0
      
      do i=iatm0,iatm1
        sdot1=sdot1+aaa(i)*bbb(i)
      enddo
      
      if(mxnode.gt.1)then
        buffer(1)=sdot1
        call gdsum(buffer(1),1,buffer(2))
        sdot1=buffer(1)
      endif
      
      return
      end function sdot1

      integer function loc2(i,j)

c*********************************************************************
c
c     calculates double index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************
      
      integer i,j
      
      loc2=(max(i,j)*(max(i,j)-1))/2+min(i,j)
      
      return
      end function loc2

      integer function loc3(i,j,k)

c*********************************************************************
c
c     calculates triple index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith september 2008
c
c*********************************************************************
      
      integer i,j,k,a,b,c,u,v,w
      
      a=max(i,j)
      b=min(a,k)
      c=min(i,j)
      u=max(a,k)
      v=max(b,c)
      w=min(b,c)
      loc3=(u*(u*u-1))/6+(v*(v-1))/2+w
      
      return
      end function loc3

      integer function loc4(i,j,k,l)

c*********************************************************************
c
c     calculates quaduple index array minimum reference
c
c     copyright daresbury laboratory
c     author w.smith september 2008
c
c*********************************************************************
      
      integer i,j,k,l,a,b,c,d,e,f,t,u,v,w
      
      a=max(i,j)
      b=max(k,l)
      c=min(i,j)
      d=min(k,l)
      e=max(c,d)
      f=min(a,b)
      t=max(a,b)
      u=max(e,f)
      v=min(e,f)
      w=min(c,d)
      loc4=((((t+2)*t-1)*t-2)*t)/24+(u*(u*u-1))/6+(v*(v-1))/2+w
      
      return
      end function loc4

      character*3 function intstr3(nnn)

c*********************************************************************
c
c     converts a 3 digit integer to a string "001" etc.
c
c     copyright daresbury laboratory
c     author w.smith november 2005
c
c*********************************************************************

      implicit none

      integer nnn

      write(intstr3,'(i3.3)')nnn

      return
      end function intstr3
      
      subroutine traject
     x  (ltraj,idnode,imcon,istraj,keytrj,natms,nstraj,nstep,tstep)

c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      use setup_module
      use config_module

      implicit none
      
      logical newjob,ltraj
      integer idnode,imcon,istraj,keytrj,natms,nstraj,nstep,i
      real(8) tstep

      save newjob
      data newjob/.true./
      
CVAM
CVAM      call VTBEGIN(70, ierr)
CVAM

      if(ltraj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if(newjob)then
          
          newjob = .false.

          open(nhist,file='HISTORY',position='append')

        endif
        
        if(nstep.eq.nstraj.or.nstep.eq.istraj)then
          
          write(nhist,'(80a1)') cfgname
          write(nhist,'(3i10)') keytrj,imcon,natms
          
        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then
          
          write(nhist,'(a8,4i10,f12.6)') 'timestep',
     x         nstep,natms,keytrj,imcon,tstep

          if(imcon.gt.0) write(nhist,'(3g12.4)') cell

          do i = 1,natms

            write(nhist,'(a8,i10,2f12.6)')
     x        atmnam(i),i,weight(i),chge(i)
            write(nhist,'(1p,3e12.4)') xxx(i),yyy(i),zzz(i)
            if(keytrj.ge.1)then
              write(nhist,'(1p,3e12.4)') vxx(i),vyy(i),vzz(i)
            endif
            if(keytrj.ge.2)then
              write(nhist,'(1p,3e12.4)') fxx(i),fyy(i),fzz(i)
            endif

          enddo

        endif

c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close (nhist)
          newjob = .true.
          
        endif
        
      endif
CVAM
CVAM      call VTEND(70, ierr)
CVAM
      return
      end subroutine traject
      
      subroutine traject_u
     x     (ltraj,idnode,imcon,istraj,keytrj,natms,nstraj,nstep,tstep)
      
c***********************************************************************
c     
c     dl_poly subroutine for writing history file at selected
c     intervals in simulation
c     
c     Unformatted, double precision version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith dec 1992.
c     
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************
      
      use setup_module
      use config_module
      
      implicit none
      
      logical newjob,ltraj
      integer idnode,imcon,istraj,keytrj,natms,nstraj,nstep,i
      real(8) tstep
      
      save newjob
      data newjob/.true./
CVAM
CVAM      call VTBEGIN(163, ierr)
CVAM
      if(ltraj.and.idnode.eq.0)then
        
c     open the history file if new job or file closed
        
        if(newjob)  then
          
          newjob = .false.
          open(nhist,file='HISTORY',form='unformatted',
     x      position='append')
          
        endif
        if(nstep.eq.nstraj.or.nstep.eq.istraj)then
          
          write(nhist) cfgname
          write(nhist) dble(natms)
          write(nhist) (atmnam(i),i=1,natms)
          write(nhist) (weight(i),i=1,natms)
          write(nhist) (chge(i),i=1,natms)
          
        endif
        
        if(mod(nstep-nstraj,istraj).eq.0)then
          
          write(nhist)dble(nstep),dble(natms),dble(keytrj),
     x      dble(imcon),tstep
          
          if(imcon.gt.0) write(nhist) cell
          
          write(nhist) (xxx(i),i = 1,natms)
          write(nhist) (yyy(i),i = 1,natms)
          write(nhist) (zzz(i),i = 1,natms)
          
          if(keytrj.ge.1)then
            write(nhist) (vxx(i),i = 1,natms)
            write(nhist) (vyy(i),i = 1,natms)
            write(nhist) (vzz(i),i = 1,natms)
          endif
          if(keytrj.ge.2)then
            write(nhist) (fxx(i),i = 1,natms)
            write(nhist) (fyy(i),i = 1,natms)
            write(nhist) (fzz(i),i = 1,natms)
          endif
          
        endif
        
c     close history file at regular intervals
        
        if(.not.newjob.and.mod(nstep,ndump).eq.0)then
          
          close (nhist)
          newjob=.true.
          
        endif
        
      endif
CVAM
CVAM      call VTEND(163, ierr)
CVAM
      return
      end subroutine traject_u
      
      subroutine timchk(ktim,time)
      
c***********************************************************************
c     
c     dlpoly timing routine for time elapsed in seconds
c     copyright daresbury laboratory
c     author w.smith nov 2003
c
c     wl
c     2008/12/23 10:29:12
c     1.9
c     Exp
c     
c***********************************************************************

      use setup_module

      implicit none

      logical init
      character*12 dat,tim,zon
      integer idnode,mynode,ktim,day
      real(8) time,told,tsum,tnow
      integer info(8)

      save init,idnode,told,tsum,day

      data init/.true./

   10 format(/,' time elapsed since job start = ',f15.3,' seconds',/)

      call date_and_time(dat,tim,zon,info)
      
      if(init)then

         tsum=0.d0
         time=0.d0
         day=info(3)
         idnode=mynode()
         told=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         init=.false.

      else

         tnow=3600.d0*dble(info(5))+60.d0*dble(info(6))+
     x         dble(info(7))+0.001d0*dble(info(8))
         if(day.ne.info(3))then
           told=told-86400.d0
           day=info(3)
         endif
         tsum=tsum+tnow-told
         told=tnow
         time=tsum

      endif

      if(ktim.gt.0.and.idnode.eq.0) write(nrite,10)time

      return
      end subroutine timchk
      
      end module utility_module
