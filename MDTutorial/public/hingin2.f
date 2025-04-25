      program hingin2

c***********************************************************************
c
c     This programme goes through all water molecules and calculates
c     the average number of hydrogen bonds and the broken bond
c     fraction. It also calculates the solute coordination number.
c
c     P.-L. Chau   October 1995
c
c***********************************************************************

c just for some quick runs on pure water
      parameter (mxatms=1200,maxfiles=10,maxnts=5000)
c      parameter (mxatms=2000,maxfiles=10,maxnts=300000)
      parameter (maxnh=8,mxtheta=50,maxr=150)
      implicit real*8(a-h,o-y)
      external iabend
      real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
      dimension nindex(mxatms)
c maximum no of bonds set to twice the number of maximum no of atoms.
c perfect tetrahedral networks is twice, so liquids are even fewer.
      dimension nbd(mxatms*2),nba(mxatms*2),nindex2(mxatms*2)
      dimension r(mxatms*2),theta(mxatms*2)
      dimension bl(0:maxr),zbl(maxr)
      dimension ba(0:mxtheta),zba(mxtheta)
      dimension hb(3),ho(3)
      dimension numhb(-1:1,maxnh),ibin(-1:1,maxr,mxtheta)
      dimension irecl(maxfiles),nrec(maxfiles),keytrj(maxfiles)
      dimension imcon(maxfiles),nstrun(maxfiles),nstraj(maxfiles)
      dimension istraj(maxfiles),tstep(maxfiles),natms(maxfiles)
      dimension cell(maxfiles,9),rcut(maxfiles),delr(maxfiles)
      dimension cell1(9)
      dimension numcor(maxfiles)
      real*4 anbp(maxnts),anss(maxnts),solco(maxnts),t,prob
      
      character*8 atmnam(mxatms)
      character*80 dumpfile(maxfiles)
      integer idnode,mxnode,nr

      open(33,file='hgout3.v4')
c set trap for underflow etc
c      i = ieee_handler("set", "invalid", iabend)
c      if (i.ne.0) write(33,*)'IEEE trapping not supported.'
c define constants
      pi = 3.141592653589793d0
c define constants
c define parameters that have to be carried over time-steps
      nsolco1 = 0
      anbp1 = 0.0d0
      anss1 = 0.0d0
      nts = 0
c upper and lower limits of theta
      thetau = pi
c defintion 1: 130 - 180 degrees
c      thetal = pi*1.3d0/1.8d0
c defintion 2: 150 - 180 degrees
c      thetal = pi*1.5d0/1.8d0
c defintion 3: 2.4 radians to pi radians
      thetal = 2.4d0
c maximum OW-OW and OW-HW distance, upper and lower limits of r, hydrogen 
c bond length
      oomax = 3.5d0
      oomax2 = oomax*oomax
c definitions 1 -- 3: bond length = 2.5 AA
c      ohmax = 2.5d0
c definitions 4: bond length = 2.15 AA
      ohmax = 2.15d0
      ohmax2 = ohmax*ohmax
      ru = ohmax
      rl = 1.5d0
c initialize the statistical stacks, first the number of hydrogen bonds
c for each molecule
      do 40 i = -1,1,2
         do 60 j = 1,maxnh
            numhb(i,j) = 0
 60      continue
c the bond length and bond angle are entered into a two-dimensional
c array
         do 80 j = 1,maxr
            do 110 k = 1,mxtheta
               ibin(i,j,k) = 0
 110        continue
 80      continue
 40   continue
c define the upper and lower bins for bond lengths and class-mark
      bl(0)=rl
      rint = (ru-rl)/dble(maxr)
      do 120 k = 1,maxr
         bl(k) = rint*dble(k) + rl
         zbl(k) = (bl(k-1)+bl(k))/2.0d0
 120  continue
c define the upper and lower bins for bond angles and class-mark
      ba(0)=thetal
      tint = (thetau-thetal)/dble(mxtheta)
      do 140 k = 1,mxtheta
         ba(k) = tint*dble(k) + thetal
         zba(k) = (ba(k-1)+ba(k))/2.0d0
 140  continue

c define parameters
      idnode=0
      mxnode=1

c provisions for processing 10 files simultaneously
      nfiles=1
      dumpfile(1)='HISTORY.19950815'
      
      write(0,*)'Files used:'
      write(33,*)'Files used:'
      do 150 i = 1,nfiles
         write(0,*)dumpfile(i)
         write(33,*)dumpfile(i)
 150  continue
      write(0,*)'Changing the limits of r and theta to see...'
      write(33,*)'Changing the limits of r and theta to see...'
      write(0,*)'lower and upper limits of r =',rl,ru
      write(0,*)'upper and lower limits of theta =',thetal,thetau
      write(33,*)'lower and upper limits of r =',rl,ru
      write(33,*)'upper and lower limits of theta =',thetal,thetau

c open input and output files
      do 160 nf=1,nfiles
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
              
 160  continue

      interval = 1
      write(0,*)'sample every ',interval,' steps.'
      write(33,*)'sample every ',interval,' steps.'
c numf is the cardinal number of the file it is reading
      numf=1
c nr is the cardinal number of the record number in the direct access file. 
c nr-1 is the cardinal number of the timestep.
      nr=1
 99   continue
c set up local variables
      do 180 i=1,9
         cell1(i) = cell(numf,i)
 180  continue
c read in data every INTERVAL steps
      nr=nr+interval
      write(0,*)'step no. =',nr-1
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
      do 190 i = 1,natms(numf)
         if (nindex(i).eq.0) nindex(i)=1
 190  continue

c tally time-steps
      nts = nts + 1

c first count number of OXYGEN atoms in the first solvation shell.
c this is the solute coordination number
      nsolco = 0
      nbulk = 0
      do 200 i = 1,natms(numf)
         if ((atmnam(i).eq."OW        ").and.(nindex(i).eq.-1)) then
            nsolco = nsolco + 1
         else if ((atmnam(i).eq."OW        ").and
     +           .(nindex(i).eq.1)) then
            nbulk = nbulk + 1
         endif
 200  continue

c hydrogen bond count set to zero initially
      icount = 0
c go through each molecule.
      do 220 i = 1,natms(numf)
c if oxygen, locate hydrogen atoms within hydrogen-bonding distance
         if (atmnam(i).eq."OW        ") then
c locate all nearest neighbours.
            do 240 j = 1,natms(numf)
c set trap for self-distance
               if (i.eq.j) goto 240
c determine distances between oxygen atoms
               if (atmnam(j).eq."OW        ") then
c evaluate OW-OW distance. maximum is defined by oomax
                  dxx = xxx(i) - xxx(j)
                  dyy = yyy(i) - yyy(j)
                  dzz = zzz(i) - zzz(j)
                  call images2(j,imcon(numf),idnode,mxnode,1,cell1
     +                 ,dxx,dyy,dzz)
                  if ((abs(dxx).gt.oomax).or.(abs(dyy).gt.oomax).or
     +                 .(abs(dzz).gt.oomax)) then
c reject if too far; get out of loop
                     goto 240
                  else
c calculate real distance
                     tdist=dxx*dxx+dyy*dyy+dzz*dzz
                     if (tdist.gt.oomax2) then
c also reject if too far; get out of loop
                        goto 240
                     else
c now calculate OW-HW distance. maximum is defined by ohmax
                        do 260 k = 1,2
                           dxx = xxx(i) - xxx(j+k)
                           dyy = yyy(i) - yyy(j+k)
                           dzz = zzz(i) - zzz(j+k)
                           call images2(j,imcon(numf),idnode,mxnode
     +                          ,1,cell1,dxx,dyy,dzz)
                           if ((abs(dxx).gt.ohmax).or.(abs(dyy).gt
     +                          .ohmax).or.(abs(dzz).gt.ohmax)) then
c beyond maximum OW-HW distance; get out of loop
                              goto 260
                           else
c calculate real distance
                              tdist=dxx*dxx+dyy*dyy+dzz*dzz
                              if (tdist.gt.ohmax2) then
c beyond maximum OW-HW distance; get out of loop
                                 goto 260
                              else
c accept but now calculate angle theta2
                                 hb(1) = dxx
                                 hb(2) = dyy
                                 hb(3) = dzz
                                 ho(1) = xxx(j) - xxx(j+k)
                                 ho(2) = yyy(j) - yyy(j+k)
                                 ho(3) = zzz(j) - zzz(j+k)
                                 call images2(j,imcon(numf),idnode
     +                                ,mxnode,1,cell1,ho(1),ho(2)
     +                                ,ho(3))
c calculate scalar product and magnitude of hb and oh vectors, and thus
c evaluate the angle theta2.
                                 scaprod = 0.0d0
                                 hbm = 0.0d0
                                 hom = 0.0d0
                                 do 280 l = 1,3
                                    scaprod = scaprod + hb(l)*ho(l)
                                    hbm = hbm + hb(l)*hb(l)
                                    hom = hom + ho(l)*ho(l)
 280                             continue
                                 hbm = sqrt(hbm)
                                 hom = sqrt(hom)
                                 theta2 = acos(scaprod/hbm/hom)
c is theta within the defined bounds? if so, accept and put in stack
                                 if ((theta2.ge.thetal).and.(theta2.le
     +                                .thetau)) then
                                    icount = icount + 1
                                    r(icount) = sqrt(tdist)
                                    theta(icount) = theta2
                                    nba(icount) = i
                                    nbd(icount) = j+k
c define three kinds of hydrogen bonds: pure solvation-shell/solvation-shell
c hydrogen bond:
                                    if ((nindex(i).eq.-1).and.
c     +                                   (nindex(j+k).eq.-1)) then
     +                                   (nindex(j).eq.-1)) then
                                       nindex2(icount) = -1
c or pure bulk-phase/bulk-phase hydrogen bond
                                    else if ((nindex(i).eq.1).and.
c     +                                   (nindex(j+k).eq.1)) then
     +                                   (nindex(j).eq.1)) then
                                       nindex2(icount) = +1
c or mixed bulk-phase/solvation-shell hydrogen bond
                                    else if (((nindex(i).eq.1).and.
c     +                                   (nindex(j+k).eq.-1)).or.
c     +                                   ((nindex(i).eq.-1).and.
c     +                                   (nindex(j+k).eq.1))) then
     +                                   (nindex(j).eq.-1)).or.
     +                                   ((nindex(i).eq.-1).and.
     +                                   (nindex(j).eq.1))) then
                                       nindex2(icount) = -3
                                    endif
                                 endif
                              endif
                           endif
 260                    continue
                     endif
                  endif
               endif
 240        continue
         endif
 220  continue

c go through all hydrogen bonds, count them, and divide by the number of water
c molecules in the bulk phase or the solvation shell
      nss = 0
      nbp = 0
      do 400 i=1,icount
         if (nindex2(i).eq.-1) then
            nss = nss + 2
         else if (nindex2(i).eq.1) then
            nbp = nbp + 2
         else if (nindex2(i).eq.-3) then
            nbp = nbp + 1
            nss = nss + 1
         endif
 400  continue
c calculate the average number of hydrogen bonds per bulk-phase water molecule.
      anbp(nts) = float(nbp) / float(nbulk)
c calculate the average number of hydrogen bonds per solvation-shell water 
c molecule.
      anss(nts) = float(nss) / float(nsolco)

c statistics that have to be carried across time-steps
c nsolco1 is the total solute coordination number, which will be divided by 
c the number of time-steps to yield the average solute coordination number.
      nsolco1 = nsolco1 + nsolco
      solco(nts) = nsolco
c anbp1 is the total average number of hydrogen bonds per bulk-phase water
c molecule. anss is that per solvation-shell water molecule. They both have
c to be divided by the number of time-steps to yield the average number of
c hydrogen bonds per water molecule in the bulk phase and the solvation shell.
      anbp1 = anbp1 + anbp(nts)
      anss1 = anss1 + anss(nts)

c if still in same file, then continue
      if (nr.lt.nrec(numf)+1) then
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

c calculate results
      solco1 = float(nsolco1) / float(nts)
      anbp1 = anbp1 / float(nts)
      anss1 = anss1 / float(nts)

c calculate standard deviation of the parameters
      do 600 i = 1,nts
         s1 = solco(i) - solco1
         s2 = anbp(i) - anbp1
         s3 = anss(i) - anss1
         ep1 = ep1 + s1
         ep2 = ep2 + s2
         ep3 = ep3 + s3
         p1 = s1*s1
         p2 = s2*s2
         p3 = s3*s3
         var1 = var1 + p1
         var2 = var2 + p2
         var3 = var3 + p3
 600  continue
      var1 = (var1 - (ep1*ep1/float(nts))) / (float(nts-1))
      sdev1 = sqrt(var1)
      var2 = (var2 - (ep2*ep2/float(nts))) / (float(nts-1))
      sdev2 = sqrt(var2)
      var3 = (var3 - (ep3*ep3/float(nts))) / (float(nts-1))
      sdev3 = sqrt(var3)

c Student t-test
c      call tutest(anbp,nts,anss,nts,t,prob)

c write out results
      write(0,*)'number of time-steps =',nts
      write(33,*)'number of time-steps =',nts
      write(0,*)' '
      write(33,*)' '

      write(0,601)solco1,sdev1
      write(0,602)anbp1,sdev2
      write(0,603)anss1,sdev3
      write(0,*)'probability of difference due to chance =',prob
 601  format('average solute coordination number =',f8.4,' +/- ',f9.6)
 602  format('av. no. of H-bonds per bulk-phase molecule =',f9.6
     +     ,' +/- ',f9.6)
 603  format('av. no. of H-bonds per solvation-shell molecule =',f9.6
     +     ,' +/- ',f9.6)

      write(33,*)'average solute coordination number =',solco1
     +     ,' +/- ',sdev1
      write(33,*)'average no. of H-bonds per bulk-phase molecule ='
     +     ,anbp1,' +/- ',sdev2
      write(33,*)'average no. of H-bonds per solvation-shell molecule ='
     +     ,anss1,' +/- ',sdev3
      write(33,*)'probability of difference due to chance =',prob

c close FORTRAN units
      close(33)
      
      end


