      program nn

c***********************************************************************
c
c     Programme for analysis of water structure
c
c     This programme goes through all OW-OW lengths and calculates
c     the number of nearest neighbours per OW
c
c     P.-L. Chau   January 1996
c
c***********************************************************************

c reduce parameters for this run
c      parameter (mxatms=2000,maxfiles=10,maxvoisin=12)
c      parameter (maxnts=30000)
      parameter (mxatms=2000,maxfiles=10,maxvoisin=12)
      parameter (maxnts=1000)
      implicit real*8(a-h,o-y)
      external iabend
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension nindex(mxatms),now(-1:1)
      dimension r(mxatms,maxvoisin)
      dimension nb(mxatms),nbhist(-1:1,maxvoisin),nv(-1:1)
      dimension irecl(maxfiles),nrec(maxfiles),keytrj(maxfiles)
      dimension imcon(maxfiles),nstrun(maxfiles),nstraj(maxfiles)
      dimension istraj(maxfiles),tstep(maxfiles),natms(maxfiles)
      dimension cell(maxfiles,9),rcut(maxfiles),delr(maxfiles)
      dimension cell1(9),nvoisin(maxnts,-1:1),voisin(maxnts,-1:1)
      dimension snav(-1:1)
      dimension numcor(maxfiles)

      character*8 atmnam(mxatms)
      character*80 dumpfile(maxfiles)
      integer idnode,mxnode,nr

      open(33,file='nnout')
c set trap for underflow etc
c      i = ieee_handler("set", "invalid", iabend)
c      if (i.ne.0) write(33,*)'IEEE trapping not supported.'
c define constants
      pi = 3.141592653589793d0
c define parameters
      cutoff = 3.5d0
      cutoff2 = cutoff*cutoff
      idnode=0
      mxnode=1
c define parameters that have to be carried over time steps
      nts=0
c sample data every INTERVAL timestep(s)
      interval = 10
c zero all variables
      nv(1) = 0
      nv(-1) = 0
      snav(1) = 0.0d0
      snav(-1) = 0.0d0
      nbmax = -10000
      nbmin = 10000
      do 20 i = 1,maxnts
         nvoisin(i,1) = 0
         nvoisin(i,-1) = 0
         voisin(i,1) = 0.0d0
         voisin(i,-1) = 0.0d0
 20   continue
      do 25 i = 1,maxvoisin
         nbhist(1,i) = 0
         nbhist(-1,i) = 0
 25   continue
c provisions for processing 10 files simultaneously
      nfiles=2
      dumpfile(1)='HIST2.19950727'
      dumpfile(2)='HIST2.19950910'

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
c special processing to locate number of OW's.
      now(1) = 0
      now(-1) = 0
      do 130 j = 1,natms(numf)
         if (atmnam(j).eq."OW        ") then
            if (nindex(j).eq.1) then
               now(1) = now(1) + 1
            else if (nindex(j).eq.-1) then
               now(-1) = now(-1) + 1
            else
               write(0,*)'wrong nindex =',nindex(j),'    abort'
               write(33,*)'wrong nindex =',nindex(j),'    abort'
               stop
            endif
         endif
 130  continue
c write out number of OW's each step
      write(0,*)now(1),now(-1)

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
                  nc=nc+1
               endif
            endif
 160     continue
         nb(j) = nc-1
 140  continue

c calculate the average number of near neighbours for this timestep
      do 200 j = 1,natms(numf)
         if (atmnam(j).ne."OW        ") goto 200
         if (nb(j).eq.0) goto 200
c find max and min number of neighbours
         if (nb(j).gt.nbmax) then
            nbmax = nb(j)
         else if (nb(j).lt.nbmin) then
            nbmin = nb(j)
         endif
c separate into solvation sphere and bulk phase water molecules
         if (nindex(j).eq.+1) then
            nvoisin(nts,1) = nvoisin(nts,1) + nb(j)
            nv(1) = nv(1) + 1
         else if (nindex(j).eq.-1) then
            nvoisin(nts,-1) = nvoisin(nts,-1) + nb(j)
            nv(-1) = nv(-1) + 1
         else
            write(0,*)'wrong nindex =',nindex(j),'    abort'
            write(33,*)'wrong nindex =',nindex(j),'    abort'
            stop
         endif
c count the number of neighbours to make a histogram
         nbhist(nindex(j),nb(j)) = nbhist(nindex(j),nb(j)) + 1
 200  continue
c calculate the average value for step nts
      do 220 j = -1,1,2
         voisin(nts,j) = dble(nvoisin(nts,j)) / dble(now(j))
 220  continue

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

c      if (nr.le.50) then
c         goto 99
c      else
c         goto 100
c      endif

 100  continue

c calculate the mean of nearest neighbours
      do 540 i = 1,nts
         snav(1) = voisin(i,1) + snav(1)
         snav(-1) = voisin(i,-1) + snav(-1)
 540  continue
      do 560 i = -1,1,2
         snav(i) = snav(i) / dble(nts)
 560  continue
c calculate standard deviation of the parameters
      ep1 = 0.0d0
      ep2 = 0.0d0
      var1 = 0.0d0
      var2 = 0.0d0
      do 600 i = 1,nts
         s1 = voisin(i,1) - snav(1)
         s2 = voisin(i,-1) - snav(-1)
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
      write(33,*)'no. of near neighbours for bulk-phase water ='
      write(33,*)'       ',snav(1),' +/- ',sdev1
      write(33,*)'no. of near neighbours for solvation-shell water ='
      write(33,*)'       ',snav(-1),' +/- ',sdev2
      write(33,*)' '
      write(33,*)'max and min :'
      write(33,*)'max no. of neighbours =',nbmax
      write(33,*)'min no. of neighbours =',nbmin
 419  format(f10.4,i8)
      write(33,*)' '
      if (nbmax.gt.maxvoisin) then
         write(33,*)' '
         write(33,*)'nbmax > maxvoisin. Execution aborted.'
         stop
      endif

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

      integer function iabend(sig,code,sigcontext)
      integer sig,code,sigcontext(5)
      write(0,*)'problems here: sig =',sig,'   code =',code
      write(33,*)'problems here: sig =',sig,'   code =',code
      call exit(0)
      end
