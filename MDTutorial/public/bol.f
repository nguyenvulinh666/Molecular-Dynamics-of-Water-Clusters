      program bol

c***********************************************************************
c
c     This programme goes through all OW-OW lengths and calculates
c     the third-order spherical harmonic to give the local bond-order
c     parameter for the liquid.
c
c     P.-L. Chau   February 1995
c
c***********************************************************************

c MAKE SURE PARAMETERS ARE LARGE ENOUGH!!
      parameter (mxatms=800,maxfiles=10,maxvoisin=12,maxbin=100)
      parameter (maxnts=3000)
      implicit real*8(a-h,o-y)
      external iabend
      complex*16 zh,zh1
      real*8 hl(-3:3)
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension nindex(mxatms)
      dimension r(mxatms,maxvoisin)
      dimension theta(mxatms,maxvoisin)
      dimension phi(mxatms,maxvoisin)
      dimension nb(mxatms),nbhist(-1:1,maxvoisin),nv(-1:1)
      dimension irecl(maxfiles),nrec(maxfiles),keytrj(maxfiles)
      dimension imcon(maxfiles),nstrun(maxfiles),nstraj(maxfiles)
      dimension istraj(maxfiles),tstep(maxfiles),natms(maxfiles)
      dimension cell(maxfiles,9),rcut(maxfiles),delr(maxfiles)
      dimension cell1(9),hbp(maxnts),hss(maxnts),nvoisin(maxnts,-1:1)
      dimension numcor(maxfiles)
      
      dimension ueberbin(0:maxbin),izentrum(-1:1,maxbin)
      real*8 zentrum(maxbin)

      character*8 atmnam(mxatms)
      character*80 dumpfile(maxfiles)
      integer idnode,mxnode,nr

      open(33,file='boout')
c define constants
      pi = 3.141592653589793d0
c evaluate constant
      fac = 4.0d0*pi/7.0d0
c sample data every INTERVAL timestep(s)
      interval = 10
c define parameters
      cutoff = 3.5d0
      cutoff2 = cutoff*cutoff
      idnode=0
      mxnode=1
      hl2max=-1.0d90
      hl2min=1.0d90
c define parameters that have to be carried over time steps
      nts=0
c zero all variables
      nv(-1) = 0
      nv(1) = 0
      do 20 i = 1,maxnts
         nvoisin(i,1) = 0
         nvoisin(i,-1) = 0
         hbp(i) = 0.0d0
         hss(i) = 0.0d0
 20   continue
      do 25 i = 1,maxvoisin
         nbhist(1,i) = 0
         nbhist(-1,i) = 0
 25   continue
c define lower and upper limits for the bond-order parameter
      do 30 j = -1,1,2
         ueberbin(0)=0.0d0
         do 40 i = 1,maxbin
            ueberbin(i) = (1.0d0/dble(maxbin))*dble(i)
            zentrum(i) = (ueberbin(i-1)+ueberbin(i))/2.0d0
            izentrum(j,i) = 0
 40      continue
 30   continue
c provisions for processing 10 files simultaneously
      nfiles=6
      dumpfile(1)='HIST2.19950517'
      dumpfile(2)='HIST2.19950522b'
      dumpfile(3)='HIST2.19950527a'
      dumpfile(4)='HIST2.19950601'
      dumpfile(5)='HIST2.19950611'
      dumpfile(6)='HIST2.19951013'

      write(0,*)'Files used:'
      write(33,*)'Files used:'
      do 150 i = 1,nfiles
         write(0,*)dumpfile(i)
         write(33,*)dumpfile(i)
 150  continue
      write(33,*)'"bond" cutoff distance =',cutoff

c open input and output files
      do 60 nf=1,nfiles
         open(nf,file=dumpfile(nf),access='direct',recl=4)
         read(nf,rec=1)irecl(nf)
         close(nf)
         open(nf,file=dumpfile(nf),access='direct',recl=irecl(nf))
         read(nf,rec=1)irecl(nf),nrec(nf),keytrj(nf),imcon(nf)
     +        ,nstrun(nf),nstraj(nf),istraj(nf),tstep(nf),natms(nf)
     +        ,(cell(nf,i),i=1,9),rcut(nf),delr(nf)
c calculate number of timesteps
         numcor(nf) = (nstrun(nf)-nstraj(nf))/istraj(nf) + 1
c for integrating v.a.c.f. to calculate D, one truncates at about 10ps.
         write(0,*)'file no.',nf
         write(0,*)' '
         write(0,*)'record length =',irecl(nf)
         write(0,*)'number of data records =',nrec(nf)
         write(0,*)'number of runs =',nstrun(nf)
         write(0,*)'trajectory dump starts at step',nstraj(nf)
         write(0,*)'trajectory dump frequency =',istraj(nf)
         write(0,*)'number of time-steps available =',numcor(nf)
         write(0,*)'time-step size (ps) =',tstep(nf)
         write(0,*)'number of atoms =',natms(nf)
         write(0,*)'keytrj =',keytrj(nf)
         write(0,*)' '
         write(0,*)'boundary conditions =',imcon(nf)
         write(0,10)cell(nf,1),cell(nf,2),cell(nf,3)
         write(0,10)cell(nf,4),cell(nf,5),cell(nf,6)
         write(0,10)cell(nf,7),cell(nf,8),cell(nf,9)
 10      format(3f12.4)

c trap to make sure that files are compatible
         if (nf.gt.1) then
            if ( (imcon(nf).ne.imcon(nf-1)) .or. (istraj(nf).ne.
     +           istraj(nf-1)) .or. (natms(nf).ne.natms(nf-1))
     +           .or. ((tstep(nf)-tstep(nf-1)).gt.1.0d-6) ) then
               write(0,*)'Incompatible files; execution terminated.'
               write(33,*)'Incompatible files; execution terminated.'
            endif
         endif
              
 60   continue

c trap to make sure arrays are large enough
      do 70 i=1,nfiles
         if (natms(i).gt.mxatms) then
            write(0,*)'MXATMS too small; execution terminated.'
            write(33,*)'MXATMS too small; execution terminated.'
         endif
 70   continue
c numf is the cardinal number of the file it is reading
      numf=1
c nr is the cardinal number of the record number in the direct access file. 
c nr-1 is the cardinal number of the timestep.
      nr=1
 99   continue
c set up local variables
      do 80 i=1,9
         cell1(i) = cell(numf,i)
 80   continue
c read in data
      nr=nr+interval
      if (keytrj(numf).eq.0) then
         read(numf,rec=nr) (nindex(j),atmnam(j),xxx(j),yyy(j)
     +        ,zzz(j),j=1,natms(numf))
      else if (keytrj(numf).eq.1) then
         read(numf,rec=nr) (nindex(j),atmnam(j),xxx(j),yyy(j)
     +        ,zzz(j),vxx,vyy,vzz,j=1,natms(numf))
      else if (keytrj(numf).eq.2) then
         read(numf,rec=nr) (nindex(j),atmnam(j),xxx(j),yyy(j)
     +        ,zzz(j),vxx,vyy,vzz,fxx,fyy,fzz
     +        ,j=1,natms(numf))
      else
         continue
      endif
      close(nf)

c special processing for pure water simulations
      do 120 i = 1,natms(numf)
         if (nindex(i).eq.0) nindex(i)=1
 120  continue

c step number
      write(0,*)'step number =',nr-1
c tally timesteps
      nts = nts + 1
c set trap
      if (nts.gt.maxnts) then
         write(0,*)'too many timesteps: adjust MAXNTS and re-run.'
         write(33,*)'too many timesteps: adjust MAXNTS and re-run.'
         stop
      endif

c calculate all other atoms with respect to a central atom
      do 140 j=1,natms(numf)
         if (atmnam(j).ne."OW        ") goto 140
         nc = 1
         do 160 k=1,natms(numf)
c set trap so that it will not calculate self-distance
            if (j.eq.k) goto 160
            if (atmnam(k).ne."OW        ") goto 160
            dxx=xxx(k)-xxx(j)
            dyy=yyy(k)-yyy(j)
            dzz=zzz(k)-zzz(j)
            call images2(j,imcon(numf),idnode,mxnode,1,cell1
     +           ,dxx,dyy,dzz)
            if ((dxx.gt.cutoff).or.(dyy.gt.cutoff).or.(dzz.gt.cutoff))
     +           then
c reject if too far; reset dxx, dyy, dzz and get out of loop
               dxx=0.0d0
               dyy=0.0d0
               dzz=0.0d0
               goto 160
            else
c calculate real distance
               tdist=dxx*dxx+dyy*dyy+dzz*dzz
               if (tdist.gt.cutoff2) then
c also reject if too far; reset dxx, dyy, dzz and get out of loop
                  dxx=0.0d0
                  dyy=0.0d0
                  dzz=0.0d0
                  goto 160
               else
c this is near enough; calculate r, theta and phi. set trap.
                  r(j,nc) = sqrt(tdist)
                  theta(j,nc) = acos(dzz/r(j,nc))
                  if (dxx.eq.0.0d0) then
                     if (dyy.ge.0.0d0) then
                        phi(j,nc) = pi/2.0d0
                     else if (dyy.lt.0.0d0) then
                        phi(j,nc) = 3.0d0*pi/2.0d0
                     endif
                  else
                     if (dxx.gt.0.0d0) then
                        if (dyy.gt.0.0d0) then
                           phi(j,nc) = atan(dyy/dxx)
                        else if (dyy.lt.0.0d0) then
                           phi(j,nc) = 2.0d0*pi + atan(dyy/dxx)
                        else
                           phi(j,nc) = 0.0d0
                        endif
                     else if (dxx.lt.0.0d0) then
                        if (dyy.gt.0.0d0) then
                           phi(j,nc) = pi + atan(dyy/dxx)
                        else if (dyy.lt.0.0d0) then
                           phi(j,nc) = pi + atan(dyy/dxx)
                        else
                           phi(j,nc) = -pi
                        endif
                     endif
                  endif
                  nc=nc+1
               endif
            endif
 160     continue
         nb(j) = nc-1
 140  continue

c calculate the average |Y(theta,phi)|^2 for each atom
      l = 3
      do 200 j = 1,natms(numf)
         if (atmnam(j).ne."OW        ") goto 200
         if (nb(j).eq.0) goto 200
         hl2 = 0.0d0
         do 240 m = -3,3,1
            zh1 = (0.0d0,0.0d0)
            hl(m) = 0.0d0
            do 260 k1 = 1,nb(j)
               call y(l,m,theta(j,k1),phi(j,k1),zh)
               zh1 = zh1 + zh
 260        continue
            hl(m) = zh1 * conjg(zh1)
            hl2 = hl2 + hl(m)
 240     continue
         hl2 = hl2 * fac / dble(nb(j)) / dble(nb(j))
c take the square root of hl2 to get the Q_3 value
         hl2 = sqrt(hl2)
c find max and min hl2
         if (hl2.gt.hl2max) then
            hl2max = hl2
            nbmax = nb(j)
         else if (hl2.lt.hl2min) then
            hl2min = hl2
            nbmin = nb(j)
         endif
c find max and min number of neighbours
         if (nb(j).gt.nbmax) then
            nbmax = nb(j)
         else if (nb(j).lt.nbmin) then
            nbmin = nb(j)
         endif
c separate into solvation sphere and bulk phase water molecules
         if (nindex(j).eq.+1) then
            hbp(nts) = hbp(nts) + hl2
            nvoisin(nts,1) = nvoisin(nts,1) + 1
         else if (nindex(j).eq.-1) then
            hss(nts) = hss(nts) + hl2
            nvoisin(nts,-1) = nvoisin(nts,-1) + 1
         else
            write(33,*)'wrong nindex =',nindex(j),'    abort'
            stop
         endif
c enter count for histogram
         do 300 i1 = 1,maxbin
            if ((hl2.ge.ueberbin(i1-1)).and.(hl2.lt.ueberbin(i1))) then
               izentrum(nindex(j),i1) = izentrum(nindex(j),i1) + 1
            endif
 300     continue
c count the number of neighbours to make a histogram
         nbhist(nindex(j),nb(j)) = nbhist(nindex(j),nb(j)) + 1
 200  continue
c calculate the average value for step nts
      hbp(nts) = hbp(nts) / dble(nvoisin(nts,1))
      hss(nts) = hss(nts) / dble(nvoisin(nts,-1))

c if still in same file, then continue
 90   if (nr.lt.nrec(numf)+1) then
         goto 99
c if one can start another file, do so
      else if (nr.eq.nrec(numf)+1) then
         if (numf.lt.nfiles) then
            numf=numf+1
            nr=1
            goto 99
         else
c if end of last file, then go and print results
            goto 100
         endif
      else
c trap for funny numbers.
         write(0,*)'Wrong nr number. Stop execution.'
         write(33,*)'Wrong nr number. Stop execution.'
         stop
      endif

 100  continue
c calculate the mean of spherical harmonic values.
      hbpav = 0.0d0
      hssav = 0.0d0
      do 580 i = 1,nts
         hbpav = hbp(i) + hbpav
         hssav = hss(i) + hssav
         nv(-1) = nvoisin(i,-1) + nv(-1)
         nv(1) = nvoisin(i,1) + nv(1)
 580  continue
      hbpav = hbpav / dble(nts)
      hssav = hssav / dble(nts)
c calculate standard deviation of the parameters
      ep1 = 0.0d0
      ep2 = 0.0d0
      var1 = 0.0d0
      var2 = 0.0d0
      do 600 i = 1,nts
         s1 = hbp(i) - hbpav
         s2 = hss(i) - hssav
         ep1 = ep1 + s1
         ep2 = ep2 + s2
         p1 = s1*s1
         p2 = s2*s2
         var1 = var1 + p1
         var2 = var2 + p2
 600  continue
      var1 = (var1 - (ep1*ep1/dble(nts))) / (dble(nts-1))
      sdev1 = sqrt(var1)
      var2 = (var2 - (ep2*ep2/dble(nts))) / (dble(nts-1))
      sdev2 = sqrt(var2)

c write out results
      write(33,*)'Number of steps used =',nts
      write(33,*)'Q_3 for bulk-phase water ='
      write(33,*)'       ',hbpav,' +/- ',sdev1
      write(33,*)'Q_3 for solvation-shell water ='
      write(33,*)'       ',hssav,' +/- ',sdev2
      write(33,*)' '
      write(33,*)'max and min :'
      write(33,*)'max bond-order parameter =',hl2max,nbmax
      write(33,*)'min bond-order parameter =',hl2min,nbmin
      write(33,*)'max no. of neighbours =',nbmax
      write(33,*)'min no. of neighbours =',nbmin
 419  format(f10.4,i8)
      write(33,*)' '
      if (nbmax.gt.maxvoisin) then
         write(33,*)' '
         write(33,*)'nbmax > maxvoisin. Execution aborted.'
         stop
      endif
      if (hl2max.gt.1.0d0) then
         write(0,*)'histogram data unreliable because hl2max > 0'
         write(33,*)'histogram data unreliable because hl2max > 0'
      endif
c output histogram data: first the b-o parameter distribution
      write(33,*)'normalized Q_3 distribution, bp + ss:'
      ibp = 0
      iss = 0
      do 440 j = 1,maxbin
         ibp = ibp + izentrum(+1,j)
         iss = iss + izentrum(-1,j)
 440  continue
      do 460 j = 1,maxbin
         bp = dble(izentrum(+1,j))/dble(ibp)
         ss = dble(izentrum(-1,j))/dble(iss)
         write(33,*)zentrum(j),bp,ss
 460  continue

c output no. of neighbours histogram
      write(33,*)'water molecules: no. of neighbours'
      do 420 j=-1,1,2
         write(33,*)'water code =',j,'   normalized frequency:'
         do 400 i=1,maxvoisin
            write(33,399)i,dble(nbhist(j,i))/dble(nv(j))
 399        format(i4,f10.4)
 400     continue
 420  continue

      end

      subroutine y(l,m,theta,phi,zh)
      implicit real*8(a-h,o-y)
      implicit complex*16(z)

c define constants
      pi = 3.141592653589793d0
      z1=(0.0d0,1.0d0)
      zp1=z1*dcmplx(phi)
      zm1=z1*dcmplx(-phi)
      z2=(0.0d0,2.0d0)
      zp2=z2*dcmplx(phi)
      zm2=z2*dcmplx(-phi)
      z3=(0.0d0,3.0d0)
      zp3=z3*dcmplx(phi)
      zm3=z3*dcmplx(-phi)

      if (l.ne.3) then
c abort execution
         write(33,*)'cannot evaluate spherical harmonic'
         stop
      else
         if (m.eq.0) then
            zh = dcmplx( (sqrt(7.0d0/pi))*(5.0d0*(cos(theta))**3
     +           - 3.0d0*cos(theta)) / 4.0d0)
         else if (m.eq.+1) then
            zh = exp(zp1) * dcmplx( (sqrt(2.1d1/pi)) * sin(theta) /
     +           8.0d0 * (1.0d0 - 5.0d0*(cos(theta)**2)) )
         else if (m.eq.-1) then
            zh = -exp(zm1) * dcmplx( (sqrt(2.1d1/pi)) * sin(theta) 
     +           / 8.0d0 * (1.0d0 - 5.0d0*(cos(theta)**2)) )
         else if (m.eq.+2) then
            zh = (1.5d1,0.0d0) * exp(zp2) * dcmplx( (sqrt(7.0d0/
     +           (3.0d1*pi))) * cos(theta) / 4.0d0 * sin(theta)**2 )
         else if (m.eq.-2) then
            zh = (1.5d1,0.0d0) * exp(zm2) * dcmplx( (sqrt(7.0d0/
     +           (3.0d1*pi))) * cos(theta) / 4.0d0 * sin(theta)**2 )
         else if (m.eq.+3) then
            zh = (-5.0d0,0.0d0) * exp(zp3) * dcmplx( (sqrt(7.0d0/
     +           (5.0d0*pi))) * sin(theta)**3 / 8.0d0 )
         else if (m.eq.-3) then
            zh = (5.0d0,0.0d0)* exp(zm3) * dcmplx( (sqrt(7.0d0/
     +           (5.0d0*pi))) * sin(theta)**3 / 8.0d0 )
         else
            write(33,*)'m > l, cannot evaluate spherical harmonic'
            stop
         endif
      endif

      return
      end
