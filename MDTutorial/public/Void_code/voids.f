      program voids


c-------------------------------------------------------------------------------
c     Copiright-computer modelling lab., University of Modena and Reggio Emilia.
c     Author- Dr. Pedone Alfonso Jenuary 2005
c
c     The voids program is available free of charge to academic institutions. 
c     Copies should be obtained from the author or with the DL_POLY package.
c
c     No claim is made that the program is free from errors, the users are responsible 
c     for checking the validity of their results.
c
c     Details of the algorithm and its application can be found into the paper:
c     "Void size distribution in MD-modelled silica glass structure"
c     G.Malavasi, M.C.Menziani,A.Pedone,U.Segre, J. of Non-Cryst. Solids, 352 (2006) 285-296
c 
c     Program for the construction of Delaunay and Voronoi tessellation
c     of a set of atoms with different radii contained in a box simulation with cubic or 
c     parallelepiped PBC applyed.
c     The program is used to analyze the void size distribution in amorphous structures.
c     The algorithm creates three output files; 
c     interstices.out contains interstitial size distributions with (Ri) and without (Riw) 
c     overlaps and bottleneck distribution (Rb)
c     tdf.out contains total distribution function of voronoi sites that lye in the empty
c     space
c     voids.trj contains delaunay S-simplexes, the voronoi sites  and interstitial 
c     sphere associated to every delaunay S-simplex.    
c     The connectivity matrix and bottlenecks for every configuration sampled.
c
c     The program is for the DL_POLY HISTORY file, the suboutine uread is used for 
c     unformatted file, the subroutine hread for formatted HISTORY file 
c     The input file 'voidstat.inp' is attached to the program. 
c
c-------------------------------------------------------------------------------- 	 

c      declarations

      implicit real*8 (a-h,o-z)
	integer*2 mxatms,mxspec,mxcan,mxdt,itdfdiv
	parameter (itdfdiv=500)
	parameter (mxatms=1600)
	parameter (mxspec=11)
	parameter (mxcan=230)
	parameter (mxdt=10000)

	

	real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
	real*8 chge(mxatms),weight(mxatms)
      real*8 vxx(mxatms),vyy(mxatms),vzz(mxatms)
      real*8 fxx(mxatms),fyy(mxatms),fzz(mxatms)
	real*8 tstep

	real*8 ri(mxdt),riw(mxdt),vx(mxdt),vy(mxdt),vz(mxdt),volumv(mxdt)
	real*8 px(mxcan),py(mxcan),pz(mxcan),ps(mxcan),cell(9)
	real*8 BOTTLE(4,mxdt)
	real*8 AA(4,4),BB(4),ratom(mxspec),ratom2(mxspec),radii(mxspec)
	real*8 alpha,beta,gamma
	real*8 radius,vt,rxij,ryij,rzij,rijsq,cutoff,distance
	real*8 tdfcut,dis2bin
	real*8 vbox,det,distv,smr,sum,dr,cost,denmean,trfactor,facpi
	real*8 nvoid(-20:90),nvoidw(-20:90),nbottle(-20:90)
	real*8 tdfden(itdfdiv),vol(itdfdiv)

	integer*2 delaunay(4,mxdt),indx(4),temp(4),NABLSTVOR(4,mxdt)
	integer*2 tag(mxcan),tagbas(mxcan)
	integer*2 atomindex(mxatms)
	integer*2 i,j,k,l,m,ii,jj,kk,ll,mm,n,nn,nd,ncount,nc,s
	integer*2 ncan,nver,numspec,iatom,natms,nvn,non,jk
	integer*2 ncnt,ncnt1,nspec
	integer*2 mmbas(mxspec),ltpbas(mxatms)

	character*8 atmnam(mxatms),speclab(mxspec)
	character*8 elename(mxspec),label(mxspec)
	character*40 filein
	character title*80,headline*40,formatt*40
	character cfgname*80,yesno*4

	logical new

	data elename/'O2-     ','Si4+    ','P5+     ','Na+     ',
     &             'Ca2+    ','Li+     ','Zn2+    ','K+      ',
     &             'Al3+    ','La3+    ','Mg2+    '/
      data ratom2/1.6,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1/
	data speclab/mxspec*' '/
	data mmbas/mxspec*0/
	data ltpbas/mxatms*0/
	data cell/9*0.0/
	
	

*******************************************************************************************
C     READ CONFIGURATION FILE
*******************************************************************************************

      write(6,"('************************************************')")
      write(6,"('*  This program is to calculate void size      *')")
      write(6,"('*  distribution in glasses modelled by means   *')")
	write(6,"('*  of MD simulations.                          *')")
	write(6,"('*                                              *')")
      write(6,"('*  PhD, Alfonso Pedone                         *')")
	write(6,"('*  University of Modena and Reggio Emilia      *')")
      write(6,"('************************************************')")


c*************************************************************************
c     initialization, open history file and read 1_st set of data
c*************************************************************************

      keytrj=0      ! only atom coordinates is needed
	
	nipt=55
      open(nipt,file='voidstat.inp')
      read(nipt,'(a)')headline
      read(nipt,'(a)')filein
	read(nipt,'(a)')formatt
      read(nipt,*)nstart,nend,interval
      read(nipt,*)rc
      read(nipt,'(i3)')nspec
      do i=1,nspec
        read(nipt,'(a8,f6.3)')label(i),radii(i)
      enddo
	read(nipt,'(a)')yesno !write output file with voronoi trajectories
      close(nipt)

	cutoff=rc**2

	
c-----------------------------------------------------------------------------
c     Decreasing the radius of all the particles to avoid atomic overlapping.
c     The atomic radius of the smaller particle is the constant used.
c     The right atomic radii of the particles are stored into the ratm2 vector
c     while the fictitious radii are stored into ratom array
c-----------------------------------------------------------------------------
      do i=1,mxspec
	do j=1,nspec

	if(elename(i).eq.label(j))then
	ratom2(i)=radii(j)
	endif

      enddo
	enddo

	cost=MINVAL(ratom2)
	
	do i=1,mxspec   !1
	ratom(i)=ratom2(i)-cost
	enddo           !1
c-----------------------------------------------------------------------------

C get the information of number of atoms and total number of steps

      if(formatt.eq.'formatted')then
      call hread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

	 else

	 call uread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
	endif

      
	
	write(6,"(60('-'))")
	write(6,"('Input history file name : ',a40)") filein
	write(6,"('Title of the file :',/,a80)")title
	write(6,"('Lattice parameter is given as :')")
	write(6,'(3f20.10)')(cell(i),i=1,3)
	write(6,'(3f20.10)')(cell(i),i=4,6)
	write(6,'(3f20.10)')(cell(i),i=7,9)
	write(6,"(60('-'))")
c       write(6,*)(atmnam(i),i=1,natms)
	write(6,"('Actual time for each time step :',f10.5)")tstep
	write(6,"(60('-'))")

****************************************************************
c get statistices of input atoms
c buildinf of typeindex(natms) and atomindex(natms) arrays
c typeindex contains 0,1 or 2 meaning that i-th atom is an
c oxygen, network former or modifier
c atomidex contains indexies that represent the species of atoms   
c atmnam is loaded by HISTORY file and contains the label of all 
c atoms in the box, these are ordered in species's block
c oxygens first, then silicon etc...
****************************************************************

      do i=1,natms
                if(atmnam(i).eq.elename(1)) then
                        atomindex(i)=1
                else if(atmnam(i).eq.elename(2)) then
                        atomindex(i)=2
                else if(atmnam(i).eq.elename(3)) then
                        atomindex(i)=3
                else if(atmnam(i).eq.elename(4)) then
                        atomindex(i)=4
                else if(atmnam(i).eq.elename(5)) then
                        atomindex(i)=5
                else if(atmnam(i).eq.elename(6)) then
                        atomindex(i)=6
                else if(atmnam(i).eq.elename(7)) then
                        atomindex(i)=7
                else if(atmnam(i).eq.elename(8)) then
                        atomindex(i)=8
                else if(atmnam(i).eq.elename(9)) then
                        atomindex(i)=9
                else if(atmnam(i).eq.elename(10)) then
                        atomindex(i)=10
                else if(atmnam(i).eq.elename(11)) then
                        atomindex(i)=11
                else
        write(6,"('Element name in history is not found in program')")
	  write(6,"('You have to update the elename vector in the   ')")
	write(6,"(' source code and compile it again               ')")
                stop
                end if
        end do   ! i

***********************************************************************
c     First, you save the species label, # of
C     atoms of those types, and flag(=ltpbas)
c speclab(numspec) contains species label present into simulation box
c specindex(numspec) contains species type index 0,1,2
c mmbas(numspec) contains into i-th position atoms number of i-th specie
***********************************************************************

      speclab(1)=atmnam(1)
      mmbas(1)=1
      ltpbas(1)=1
      numspec=1
      new=.true.

      do i=2,natms
      do j=1,numspec

      if(atmnam(i).eq.speclab(j)) then
      mmbas(j)=mmbas(j)+1
      ltpbas(i)=j
      new=.false.
      end if

      end do


      if(new) then
      numspec=numspec+1
      mmbas(numspec)=mmbas(numspec)+1
      ltpbas(i)=numspec
      speclab(numspec)=atmnam(i)
      end if
      new=.true.

      end do

	write(6,"('Total number of Atoms :',i10)")natms
      write(6,"('index    species_label      number')")
      do i=1,numspec
      write(6,'(i3,10x,a8, 7x,i5)') i,speclab(i),mmbas(i)
      end do 
	write(6,"(60('-'))")

c --------------------------------------------------------------

	nkount=1
	nfirst=nstep
      do while(iflg.ne.-1)
  
      if(formatt.eq.'formatted')then
      call hread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

	 else

	 call uread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
	endif

      nkount=nkount+1
	if(nkount.eq.2) interval=nstep-nfirst
      nfinal=nstep
	enddo    ! HISTORY while statement

	write(6,"(60('-'))")
	write(6,"('Starting time step :',i10)")nfirst
	write(6,"('Final time step    :',i10)")nfinal
	write(6,"('Interval between saved timestep :',i10)")interval
	ntotal=(nfinal-nfirst)/interval+1
	write(6,"('Total number of time step saved :',i10)")ntotal 
	write(6,"(60('-'))")
	ntotal=(nend-nstart)/interval+1
      write(6,"('There are ',i4,' configurations.')")ntotal
      
      nacc=ntotal

	noutput=ntotal/nacc
	if(mod(ntotal,nacc).gt.0) noutput=noutput+1
	write(6,"(60('-'))")
	write(6,"('You will have ',i5,' output files')")noutput 
	write(6,"(60('-'))")

	

C-------End of preliminary ---------------------------------------
	
	vbox=cell(1)*cell(5)*cell(9)

      free=0.0
	do k=-20,60
	nvoid(k)=0.
	nvoidw(k)=0.
	nbottle(k)=0.
	enddo


	do k=1,itdfdiv
      tdfden(k)=0.0
	enddo

      nvn=0
	non=0
      nvert=0
      if(yesno.eq.'yes')open(36,file='voids.trj')

C-------Start of main program-----------------------------------



C      print*,'iflg before setting ',iflg
	iflg=0

       if(formatt.eq.'formatted')then
      call hread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

	 else

	 call uread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
	endif
 
      
C --------------------------------------------------------------------
C READ from history file till iflg=-1 and calculate several structural
C properties
C---------------------------------------------------------------------


      write(6,"(' You may need to wait for a while...')")

      ncnt=0
      ncnt1=0
      do while(iflg.ne.-1)
	

      ncnt=ncnt+1

      if(nstep.ge.nstart.and.nstep.le.nend) then  ! 12-30-1997
	
      ncnt1=ncnt1+1
      
c     initialization
*************************************************************main

***********************************************************************
C     Loops on all atoms, we construct VP for each atom
c     main part of the program 
***********************************************************************

      sum=0.0     !  counter of the space covered during tessellation

	s=0         !  number of delaunay tetrahedron found
	do jk=1,MXDT
	volumv(jk)=0.
	RI(jk)=0.
	RIW(jk)=0.
	do kk=1,4
	DELAUNAY(kk,jk)=0
	enddo
	enddo


      do i=1,natms   !3
	iatom=atomindex(i)
c------------------------------------------------------------------
c     candidates selection, only atoms inside a cutoff distance
c     from the i-th atom are used, this improves the efficiency of 
c     program.
c     Keep in mind!! if the cutoff is too short the tessellation is 
c     wrong
c-------------------------------------------------------------------

      ncan=0
	do j=1,natms   !4

	if(j.ne.i)then   !5
	rxij=xxx(j)-xxx(i)
	ryij=yyy(j)-yyy(i)
	rzij=zzz(j)-zzz(i)

	if(abs(rxij).GT.cell(1)/2)rxij=rxij-sign(cell(1),rxij)
	if(abs(ryij).GT.cell(5)/2)ryij=ryij-sign(cell(5),ryij)
	if(abs(rzij).GT.cell(9)/2)rzij=rzij-sign(cell(9),rzij)

      rijsq=rxij**2+ryij**2+rzij**2

c------------------------------------------------------------------------
c     The reference system has been shifted on the i-th atom
c     all the coordinates of candidates are referred to the i atom
c     
c     If the j-th atom is inside the cutoff distance, it is a candidate
c     for the constryction of the voronoi polyhedron of the i atom
c     or similarly for the construction of the set of delaunay 
c     tetrahedra sharing the i atom.
c------------------------------------------------------------------------
     
      if(rijsq.LT.cutoff)then !6
	ncan=ncan+1
	
	  if(ncan.gt.mxcan)then  !7
	  write(*,*)'too many candidates'
	  stop
	  endif    !7

      ps(ncan)=rijsq !square distances of the candidates from i
	px(ncan)=rxij  !coordinates of candidates
	py(ncan)=ryij
	pz(ncan)=rzij
      tag(ncan)=j    ! indexes of atoms
	tagbas(ncan)=atomindex(j) !indexes of species
      
	endif !6

      endif !5
	enddo !4

c------------------------------------------------------------------------
c     loops on candidates, we explore the quadruplets of all the DT's 
c     sharing the i atom
c------------------------------------------------------------------------

      do l=1,ncan-2 !8
	if(tag(l).gt.i)then
	do m=l+1,ncan-1 !9
	if(tag(m).gt.i)then
	do k=m+1,ncan  !10
	if(tag(k).gt.i)then
	

c----------------------------------------------------------
c     construction of the linear system AX=B
c     necessary to find the equation of the sphere touching
c     the surfaces of the i,l,m,k quadruplets of atoms
c----------------------------------------------------------

	AA(1,1)=0.
      AA(2,1)=px(l)
      AA(3,1)=px(m)
      AA(4,1)=px(k)

	AA(1,2)=0.
      AA(2,2)=py(l)
      AA(3,2)=py(m)
      AA(4,2)=py(k)

	AA(1,3)=0.
      AA(2,3)=pz(l)
      AA(3,3)=pz(m)
      AA(4,3)=pz(k)

	AA(1,4)=1.0
	AA(2,4)=1.
      AA(3,4)=1.		
	AA(4,4)=1.

	BB(1)=ratom(iatom)**2
	BB(2)=-(PS(l)-ratom(TAGBAS(l))**2)
	BB(3)=-(PS(m)-ratom(TAGBAS(m))**2)
	BB(4)=-(PS(k)-ratom(TAGBAS(k))**2)

C-----------------------------------------------------------------
c     controllare nel compilatore intel per linux
c     le librerie che risolvono i sistemi lineari
c     e il determinante di una matrice e vedere come si fanno
c     le chiamate
C-----------------------------------------------------------------

c----------------------------------------
C     LU decomposition
c----------------------------------------
	call ludcmp(AA,4,4,indx,det)
c----------------------------------------
C     determinat of the A matrix 
c----------------------------------------
	
	do kk=1,4 !11
	det=det*AA(kk,kk)
	enddo !11

c----------------------------------------------------------------
c     if DET(A)=0 the matrix is singular and there are not 
c     solutions, this often happens with crystalline structures
c     because four points could lie in the same plane
c     but in  this cases we can't find a DT
c----------------------------------------------------------------

      if(abs(det).GT.0.0001)then !12

c     resolution of the linear system
      
	call lubksb(AA,4,4,indx,BB)

c     calculation of the radius of the empty sphere touching
c     the quadruplet i,tag(l),tag(m),tag(k)

	radius=BB(1)**2+BB(2)**2+BB(3)**2-4*BB(4)

	     if(radius.LT.0)then !13
	     write(*,*)'warning - radius^2 <0 for the quadruplet :'
	     write(*,*)i,tag(l),tag(m),tag(k)
	     radius=0
	     endif !13
      
	radius=sqrt(radius)/2

c     calculation of the center of the sphere

      alpha=-BB(1)/2
	beta=-BB(2)/2
	gamma=-BB(3)/2

c----------------------------------------------------------------
c     we check that no other ii surface particles intersect the 
c     sphere just calculated
c----------------------------------------------------------------

      nd=1 ! flag of the sphere, if nd=1 the sphere is empty
	     ! elseif nd=-1 the sphere touching the i-l-m-k quadruplet
	     ! is not empty end the quadruplet is not a DT

      do ii=1,ncan !14
	if(ii.ne.k.and.ii.ne.m.and.ii.ne.l)then !15

c     distance between the center of the sphere and the guess particle

      distv=(alpha-px(ii))**2+(beta-py(ii))**2+(gamma-pz(ii))**2

c     sum between the radius of the guess particle and the radius
c     of the sphere that touch the quadruplet i-l-m-k

      smr=ratom(tagbas(ii))**2+radius**2

	       if(smr.gt.distv)then !16
             nd=-1
	       goto 11 
	       endif !16


	endif !15
	enddo !14
   11 continue
      
	if(nd.eq.1)then !17

c     the i-l-m-k quadruplet form a DT
      
c	write(30,*)i,tag(l),tag(m),tag(k)
	temp(1)=i
	temp(2)=tag(l)
	temp(3)=tag(m)
	temp(4)=tag(k)

c     the temp array is sorted in increasing order, this is useful
c     in the follow because we want store only DT not previously foun
      
c	call SORT(temp,4)


      
	s=s+1

	if(s.gt.mxdt)then
	write(*,'(a30,i7)')'too many voronoi vertices found',s
	write(*,*)'you need to increase the mxdt parameter'
	  stop
	endif

	  do nn=1,4 !22
	  delaunay(nn,s)=temp(nn)
	  enddo !22

c     the coordinates of the voronoi vertexes are stored
c     the reference system is shifted on the corner of the box

      vx(s)=xxx(i)+alpha
	vy(s)=yyy(i)+beta
	vz(s)=zzz(i)+gamma
c     radius of the interstitial voids
      ri(s)=radius-cost
	riw(s)=radius-cost
c     volume of the DT (vt) and volume of the space covered (sum)

      vt=abs(det)/6

      call volume(temp,VT,ratom2,VV,xxx,yyy,zzz,mxatms,cell,atomindex)
      
	volumv(s)=VV

	sum=sum+vt
      if(sum.ge.vbox)goto 13 ! the entire space is covered

	
	endif !17
	endif !12


      endif
	enddo !10
	endif
	enddo !9
	endif
	enddo !8
	enddo !3
   13 continue

      nver =s  ! number of voronoi vertexes
	print*,'step =',ncnt1,'n° ver =',nver

      nvn=nvn+nver
c-------------------------------------------------------
c     sorting voronoi radius into increasing order
c--------------------------------------------------------

      call SORT2(ri,riw,mxdt,nver,vx,vy,vz,delaunay)

c     eliminate the overlaid	

      call OVERLAID(riw,cell,mxdt,nver,vx,vy,vz)

	
	call TDF(riw,cell,mxdt,nver,vx,vy,vz,tdfden,nv)

	nvert=nvert+nv
	
        call CONNETTIVITY(NVER,delaunay,NABLSTVOR,BOTTLE,nbottle,
     &            atomindex,xxx,yyy,zzz,ratom,cell,mxatms,cost,mxdt)

	if(yesno.eq.'yes')then

	write(36,"('timestep = ',i6,'  nver = ',i6)")nstep,nver
	write(36,'(a)')
     &'  NV      DT verteces       rv        vx        vy        vz'   
      do i=1,nver
	write(36,'(5i5,2x,4f7.3)')i,(delaunay(j,i),j=1,4),ri(i),vx(i),
     &      vy(i),vz(i)
	enddo
	write(36,'(a)')'CONNECTIVITY MATRIX'
	write(36,"('vertex    1°     2°     3°     4°  neighbours')")
	do i=1,nver
	write(36,'(5i6)')i,(NABLSTVOR(j,i),j=1,4)

	enddo

      write(36,'(a)')'BOTTLENECK MATRIX'
	write(36,"('vertex    1°     2°     3°     4° bottleneck radii')")
	do i=1,nver
	write(36,'(i5,2x,4f7.3)')i,(BOTTLE(j,i),j=1,4)
	enddo

	endif


      dr=0.05
	do mm=1,nver !24

	free=free+volumv(mm)

	kk=anint(ri(mm)/dr)
	nvoid(kk)=nvoid(kk)+1.0
	if(riw(mm).gt.-10.)then
	nvoidw(kk)=nvoidw(kk)+1.0
	non=non+1
	endif
	enddo !24
			   

**************************************************************************
c     end main part of the algorithm
*************************************************************end main	
	end if   				 ! 12-30-1997


       if(formatt.eq.'formatted')then
      call hread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

	 else

	 call uread
     x  (filein,title,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
	endif

      enddo    ! HISTORY

************************************************************
C     PRINT RESULTS
************************************************************
 
 
      free=free/dble(ncnt1)

      write(6,"('-----------------------------------------------')")
	write(6,"('*                                             *')")
      write(*,"('* AVERAGE FREE VOLUME =',f10.3)")free
	write(*,"('* SIMULATION CELL VOLUME =',f10.3)")vbox

	free=free/vbox

	write(*,"('* FRACTIONAL FREE VOLUME (FFV) =',f7.3)")free
	write(6,"('*                                             *')")
	write(6,"('-----------------------------------------------')")
 
 	open(34,file='interstices.out')
	write(34,'(a)')
     &	' dr        ri       riw         rb'    
	do k=-20,90
        a=nvoid(k)/dble(ncnt1)
        b=nvoidw(k)/dble(ncnt1)
        c=nbottle(k)/dble(ncnt1)
	write(34,'(4f9.3)')k*dr,a,b,c
	enddo
	close(34)
 
****************************************************************************
c     TDF VOIDS NORMALIZATION 
****************************************************************************
     
	tdfcut=8.0d0
	dis2bin=dble(itdfdiv)/tdfcut !inverso di dr

      dr=1.0d0/dis2bin
	dr2=0.5d0*dr
      facpi=(4.0d0/3.0d0)*3.141592654

	do j=1,itdfdiv
	r=dr*j
	r1=r+dr2
	r2=r-dr2
	vol(j)=facpi*(r1**3-r2**3)
	enddo
      
	denmean=nvert/(vbox*dble(ncnt1))

      do i=1,itdfdiv-1
	distance=i*dr
      trfactor=4.0*distance*3.141592654d0*denmean
c     trfactor=1.0
      tdfden(i)=2.0*tdfden(i)*trfactor/(vol(i)*denmean*nvert)        
      end do

	open(35,file='tdf.out')
	write(35,'(a)')'  r       tdf  '
	do k=1,itdfdiv-1
	write(35,'(f7.3,3x,f9.3)')k*0.016,tdfden(k)
	enddo
	close(35)

      write(6,"('END PROGRAM')") 
      stop
      end

***********************************************************************
C     STARTING SUBROUTINE
***********************************************************************

      subroutine CONNETTIVITY(NVER,DELA,NABLSTVOR,BOTTLE,nbottle,
     &            atomindex,xxx,yyy,zzz,ratom,cell,mxatms,cost,mxdt)

**********************************************************************
*     this subroutine builds two new arrays NABLSTVOR(NVER,8) 
*     and BOTTLE(NVER,8) in which there are connettivity 
*     between voronoi vertices and the bottleneck joining them.
**********************************************************************
      INTEGER*2 NVER,mxdt,kk,i,j,k,m,nn,mxatms,nj
      INTEGER*2 NABLSTVOR(4,mxdt),templab(3),cord(mxdt) 
	REAL*8 BOTTLE(4,mxdt),BX,BY,BZ,cell(9),cost,dr,RB
	real*8 nbottle(-20:60)
	real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms),ratom(11)
	INTEGER*2 DELA(4,mxdt),atomindex(mxatms)

	do i=1,NVER
	cord(i)=0
	do j=1,4
	NABLSTVOR(j,i)=0
	BOTTLE(j,i)=0.0
	enddo
	enddo

      nj=0
	do i=1,NVER-1
	do j=i+1,NVER
	nn=0
	do k=1,4
	do m=1,4
	if(DELA(k,i).eq.DELA(m,j))then
	nn=nn+1
      templab(nn)=DELA(k,i)
	endif
	enddo !m
	enddo !k

	if(nn.eq.3)then     ! shared face between i-th and j-th vertices

	call bottleneck(templab,BX,BY,BZ,RB,atomindex,xxx,yyy,zzz,
     &                ratom,cell,mxatms,cost)
	

	cord(i)=cord(i)+1
	cord(j)=cord(j)+1

	NABLSTVOR(cord(i),i)=j
	NABLSTVOR(cord(j),j)=i

	BOTTLE(cord(i),i)=RB
	BOTTLE(cord(j),j)=RB
	
	dr=0.05
      kk=anint(RB/dr)
	nbottle(kk)=nbottle(kk)+1.0

	endif
	enddo !j
	enddo !i

	return
	end

 

      Subroutine bottleneck(lab,BX,BY,BZ,RB,atomindex,xxx,yyy,zzz,ratom,
     &            cell,mxatms,cost)

***************************************************************************
*     this subroutine calculates the radii of spheres that are able to 
*     pass through the face defined by the three vertices contain into 
*     lab arrays
***************************************************************************

      integer*2 lab(3),INDX(4),mxatms
	real*8 RB,RADIUS
	real*8 X1(3),X2(3),X3(3)
	
	real*8 AA(4,4),BB(4),PS2,PS3
	real*8 A,B,C,D
	integer*2 atomindex(mxatms)
	real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms),ratom(11),cost,cell(9)
	
	
	

	X1(1)=0.
	X1(2)=0.
	X1(3)=0.

	X2(1)=xxx(lab(2))-xxx(lab(1))
	X2(2)=yyy(lab(2))-yyy(lab(1))
	X2(3)=zzz(lab(2))-zzz(lab(1))
	
	IF(abs(X2(1)).GT.cell(1)/2)X2(1)=X2(1)-SIGN(cell(1),X2(1))
	IF(abs(X2(2)).GT.cell(5)/2)X2(2)=X2(2)-SIGN(cell(5),X2(2))
	IF(abs(X2(3)).GT.cell(9)/2)X2(3)=X2(3)-SIGN(cell(9),X2(3))

	PS2=X2(1)**2+X2(2)**2+X2(3)**2

      X3(1)=xxx(lab(3))-xxx(lab(1))
	X3(2)=yyy(lab(3))-yyy(lab(1))
	X3(3)=zzz(lab(3))-zzz(lab(1))

	IF(abs(X3(1)).GT.cell(1)/2)X3(1)=X3(1)-SIGN(cell(1),X3(1))
	IF(abs(X3(2)).GT.cell(5)/2)X3(2)=X3(2)-SIGN(cell(5),X3(2))
	IF(abs(X3(3)).GT.cell(9)/2)X3(3)=X3(3)-SIGN(cell(9),X3(3))

	PS3=X3(1)**2+X3(2)**2+X3(3)**2
c     components of the vector normal to the plane between three points

      A=X2(2)*X3(3)-X2(3)*X3(2)
	B=X3(1)*X2(3)-X2(1)*X3(3)
	C=X2(1)*X3(2)-X2(2)*X3(1)	
      D=0.
c     searching of the voronoi vertex between three particles 
c     in the plane (triangular face of the delaunay tetrahedron
 
	AA(1,1)=-2*X1(1)
	AA(2,1)=-2*X2(1)
	AA(3,1)=-2*X3(1)
	AA(4,1)=A

	AA(1,2)=-2*X1(2)
	AA(2,2)=-2*X2(2)
	AA(3,2)=-2*X3(2)
      AA(4,2)=B

	AA(1,3)=-2*X1(3)
	AA(2,3)=-2*X2(3)
      AA(3,3)=-2*X3(3)
      AA(4,3)=C

	AA(1,4)=1.
	AA(2,4)=1.
	AA(3,4)=1.
	AA(4,4)=0.

	BB(1)=ratom(atomindex(lab(1)))**2
	BB(2)=-(PS2-ratom(atomindex(lab(2)))**2)
	BB(3)=-(PS3-ratom(atomindex(lab(3)))**2)
	BB(4)=-D

	call ludcmp(AA,4,4,INDX,D)
	call lubksb(AA,4,4,INDX,BB)

      RADIUS=BB(1)**2+BB(2)**2+BB(3)**2 -BB(4)
	if(RADIUS.lt.0)RADIUS=0.0
	RADIUS=SQRT(RADIUS)
	RB=RADIUS-cost
      
      return
	end

      subroutine hread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)

c
c***********************************************************************
c
c     dl_poly subroutine for reading the formatted history file
c
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c
c     single processor version
c
c     wl
c     1996/02/15 14:33:26
c     1.1.1.1
c     Exp
c
c***********************************************************************
c

      implicit real*8(a-h,o-z)

      logical new

      character*80 cfgname
      character*40 history
      character*8 atmnam(*),step

      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)

      save new

      data new/.true./,nhist_f/77/

C      iflg=0

c     open history file if new job

      if(iflg.eq.0)then

        open(nhist_f,file=history,status='old',err=100)

        read(nhist_f,'(a80)',err=200) cfgname
        write(*,'(a,a)')'# History file header: ',cfgname
        read(nhist_f,'(2i10)',end=200) ktrj,imcon
        if(keytrj.gt.ktrj)then

          if(ktrj.eq.0)write(*,'(a)')'# error - no velocities in file'
          if(keytrj.gt.1)write(*,'(a)')'# error - no forces in file'
          stop

        endif

        new=.false.

      endif



      read(nhist_f,'(a8,4i10,f12.6)',end=200)
     x     step,nstep,natms,ktrj,imcon,tstep

C      if(natms.ne.matms)then

C        write(*,'(a)')'# error - incorrect number of atoms in file'
C        write(*,'(a,i6,a)')'# file contains',matms,' atoms'
C        stop

C      endif

      if(imcon.gt.0) read(nhist_f,'(3g12.4)',end=200) cell

      do i = 1,natms
        read(nhist_f,'(a8,i10,2f12.6)',end=200)
     x    atmnam(i),j,weight(i),chge(i)
        read(nhist_f,'(1p,3e12.4)',end=200) xxx(i),yyy(i),zzz(i)
        if(keytrj.ge.1)then
          read(nhist_f,'(1p,3e12.4)',end=200) vxx(i),vyy(i),vzz(i)
        else if(ktrj.ge.1)then
          read(nhist_f,'(1p,3e12.4)',end=200) vx,vy,vz
        endif
        if(keytrj.ge.2)then
          read(nhist_f,'(1p,3e12.4)',end=200) fxx(i),fyy(i),fzz(i)
        else if(ktrj.ge.2)then
          read(nhist_f,'(1p,3e12.4)',end=200) fx,fy,fz
        endif
      enddo

      iflg=1

      return

  100 continue

      write(*,'(a)')'# error - History file not found'
      stop

  200 continue
c      write(*,'(a)')'# warning - end of History file encountered'
      close (nhist_f)
      iflg=-1

      return
      end

	subroutine uread
     x  (history,cfgname,atmnam,iflg,imcon,keytrj,natms,nstep,tstep,
     x   cell,chge,weight,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz)
      
c     
c***********************************************************************
c     
c     dl_poly subroutine for reading unformatted history files
c     
c     double precision, single processor version
c     
c     copyright - daresbury laboratory 1996
c     author    - w. smith jan 1996.
c     
c***********************************************************************
c     
      
      real*8 datms,dimcon,trjkey,dstep,tstep
      integer*2 ktrj,keytrj,nhist_u,natms,i,j
	integer*2 imcon,nstep,iflg
      character*80 cfgname
      character*40 history
      character*8 atmnam(*)
      
      dimension cell(9)
      dimension chge(*),weight(*)
      dimension xxx(*),yyy(*),zzz(*)
      dimension vxx(*),vyy(*),vzz(*)
      dimension fxx(*),fyy(*),fzz(*)
      
      data nhist_u/77/

c     open the history file if new job
      
      if(iflg.eq.0)then
        
        open(nhist_u,file=history,form='unformatted',
     x       status='old',err=100)
	
        read(nhist_u,err=300) cfgname
        write(6,'(a,a)')'# History file header: ',cfgname
        read(nhist_u,end=300) datms
	natms=nint(datms)
        read(nhist_u,end=300) (atmnam(i),i=1,natms)
        read(nhist_u,end=300) (weight(i),i=1,natms)
        read(nhist_u,end=300) (chge(i),i=1,natms)

        
      endif
      
      read(nhist_u,end=300)dstep,datms,trjkey,dimcon,tstep
      nstep=nint(dstep)
      ktrj=nint(trjkey)
      imcon=nint(dimcon)
      if(keytrj.gt.ktrj)then
        
        if(ktrj.eq.0)write(6,'(a)')'# error - no velocities in file'
        if(keytrj.gt.1)write(6,'(a)')'# error - no forces in file'
        stop

      endif
      
      if(imcon.gt.0) read(nhist_u,end=300) cell

	natms=nint(datms)
      
      read(nhist_u,end=300) (xxx(i),i = 1,natms)
      read(nhist_u,end=300) (yyy(i),i = 1,natms)
      read(nhist_u,end=300) (zzz(i),i = 1,natms)
      
      if(keytrj.ge.1)then
        read(nhist_u,end=300) (vxx(i),i = 1,natms)
        read(nhist_u,end=300) (vyy(i),i = 1,natms)
        read(nhist_u,end=300) (vzz(i),i = 1,natms)
      else if(ktrj.ge.1)then
        read(nhist_u,end=300)
        read(nhist_u,end=300)
        read(nhist_u,end=300)
      endif
      if(keytrj.ge.2)then
        read(nhist_u,end=300) (fxx(i),i = 1,natms)
        read(nhist_u,end=300) (fyy(i),i = 1,natms)
        read(nhist_u,end=300) (fzz(i),i = 1,natms)
      else if(ktrj.ge.2)then
        read(nhist_u,end=300)
        read(nhist_u,end=300)
        read(nhist_u,end=300)
      endif
      
      iflg=1

      return

  100 continue

      write(6,'(a)')'# error - History file not found'
      stop

  300 continue
c      write(6,'(a)')'# warning - end of History file encountered'
      close (nhist_u)
      iflg=-1

      return
      end


	
 	SUBROUTINE ludcmp(a, n, np, indx, d)

C ***************************************************************************
C **  'NUMERICAL RECIPES IN FORTRAN', Second Edition                       **
C **  by W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery    **
C **  Chapter 2, p. 38                                                     **
C ***************************************************************************

	INTEGER*2 n, np, indx(n), NMAX
	REAL*8 d, a(np,np), TINY
	PARAMETER (NMAX = 500, TINY = 1.0e-20)
	INTEGER*2 i, imax, j, k
	REAL*8 aamax, dum, sum, vv(NMAX)

	d = 1.
	do i = 1,n
	aamax = 0
	do j = 1,n
	if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
	enddo ! j
	if (aamax.eq.0.) pause 'singular matrix in ludcmp'
	vv(i) = 1./aamax
	enddo ! i
	do j = 1,n
	do i = 1,j-1
	sum = a(i,j)
	do k = 1,i-1
	sum = sum-a(i,k)*a(k,j)
	enddo ! k
	a(i,j)=sum
	enddo ! i
	aamax = 0.
	do i = j,n
	sum = a(i,j)
	do k=1,j-1
	sum=sum-a(i,k)*a(k,j)
	enddo ! k
	a(i,j)=sum
	dum=vv(i)*abs(sum)
	if(dum.ge.aamax)then
	imax=i
	aamax=dum
	endif
	enddo ! i
	if(j.ne.imax)then
	do k = 1,n
	dum=a(imax,k)
	a(imax,k)=a(j,k)
	a(j,k)=dum
	enddo ! k
	d=-d
	vv(imax)=vv(j)
	endif
	indx(j)=imax
	if(a(j,j).eq.0.)a(j,j)=TINY
	if(j.ne.n)then
	dum=1./a(j,j)
	do i=j+1,n
	a(i,j)=a(i,j)*dum
	enddo ! i
	endif
	enddo ! j
	return
	END
	



	SUBROUTINE lubksb(a, n, np, indx, b)

C ***************************************************************************
C **  'NUMERICAL RECIPES IN FORTRAN', Second Edition                       **
C **  by W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery    **
C **  Chapter 2, p. 39                                                     **
C ***************************************************************************

	INTEGER*2 n, np, indx(n)
	REAL*8 a(np,np),b(n)
	INTEGER*2 i, ii, j, ll
	REAL*8 sum

	ii=0
	do i=1,n
	ll=indx(i)
	sum=b(ll)
	b(ll)=b(i)
	if(ii.ne.0)then
	do j=ii,i-1
	sum=sum-a(i,j)*b(j)
	enddo
	else if(sum.ne.0.)then
	ii=i
	endif
	b(i)=sum
	enddo
	do i=n,1,-1
	sum=b(i)
	do j=i+1,n
	sum=sum-a(i,j)*b(j)
	enddo
	b(i)=sum/a(i,i)
	enddo
	return
	END



      SUBROUTINE SORT ( RS,NN )

C    *******************************************************************
C    ** ROUTINE TO SORT NEIGHBOURS INTO INCREASING ORDER OF DISTANCE  **
C    **                                                               **
C    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
C    *******************************************************************

        INTEGER*2  NN
        INTEGER*2    RS(NN),RSI
        

        LOGICAL CHANGE
        INTEGER*2 I, ITOP, I1
        

C    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NN - 1

1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP

              I1 = I + 1

              IF ( RS(I) .GT. RS(I1) ) THEN

                 RSI = RS(I)
                 
                 RS(I) = RS(I1)
                 
                 RS(I1) = RSI
                 
                 CHANGE = .TRUE.

              ENDIF

100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

        RETURN
        END

         SUBROUTINE SORT2( RI,RIW,mxdt,nver,vx,vy,vz,delaunay)

C    *******************************************************************
C    ** ROUTINE TO SORT NEIGHBOURS INTO DECREASING ORDER OF DISTANCE  **
C    **                                                               **
C    ** FOR SIMPLICITY WE USE A BUBBLE SORT - OK FOR MAXCAN SMALL.    **
C    *******************************************************************

        INTEGER*2 MXDT, NVER
        REAL*8    VX(MXDT), VY(MXDT), VZ(MXDT), RI(MXDT),RIW(MXDT)
        INTEGER*2 delaunay(4,MXDT)

        LOGICAL CHANGE
        INTEGER*2 I, ITOP, I1, DTI(4),k1
        REAL*8    VXI, VYI, VZI, RII,RIIW
	

C    *******************************************************************

        CHANGE = .TRUE.
        ITOP = NVER - 1

1000    IF ( CHANGE .AND. ( ITOP .GE. 1 ) ) THEN

           CHANGE = .FALSE.

           DO 100 I = 1, ITOP

              I1 = I + 1

              IF ( RI(I) .LT. RI(I1) ) THEN

                 VXI = VX(I)
                 VYI = VY(I)
                 VZI = VZ(I)
                 RII = RI(I)
	           RIIW=RIW(I)
                 
	           do k1=1,4
	           DTI(k1)=delaunay(k1,I)
	           enddo

                 VX(I) = VX(I1)
                 VY(I) = VY(I1)
                 VZ(I) = VZ(I1)
                 RI(I) = RI(I1)
                 RIW(I)=RIW(I1)

                 do k1=1,4
	           delaunay(k1,I)=delaunay(k1,I1)
	           enddo

                 VX(I1) = VXI
                 VY(I1) = VYI
                 VZ(I1) = VZI
                 RI(I1) = RII

                 RIW(I1)=RIIW
                 do k1=1,4
	           delaunay(k1,I1)=DTI(k1)
	           enddo

	
	           CHANGE = .TRUE.

              ENDIF

100        CONTINUE

           ITOP = ITOP - 1
           GOTO 1000

        ENDIF

        RETURN
        END

      SUBROUTINE OVERLAID(RIW,cell,mxdt,nver,vx,vy,vz)

*************************************************************************
C     ruotine for eliminating the overlaids
*************************************************************************

      INTEGER*2 NVER,MXDT,i,j
	REAL*8 RIW(MXDT),cell(9),VX(MXDT),VY(MXDT),VZ(MXDT)
      REAL*8 DI,RR,VXI,VYI,VZI


	do i=1,NVER-1

      if(RIW(i).GT.0)then

	do j=i+1,NVER
	if(RIW(j).GE.0)then

	VXI=VX(j)-VX(i)
	VYI=VY(j)-VY(i)
	VZI=VZ(j)-VZ(i)

	IF(abs(VXI).GT.cell(1)/2)VXI=VXI-SIGN(cell(1),VXI)
	IF(abs(VYI).GT.cell(5)/2)VYI=VYI-SIGN(cell(5),VYI)
	IF(abs(VZI).GT.cell(9)/2)VZI=VZI-SIGN(cell(9),VZI)

	DI=VXI**2+VYI**2+VZI**2
	DI=SQRT(DI)

	RR=RIW(j)+RIW(i)

	IF(DI.LT.RR)RIW(j)=-10.

	endif
	enddo
	endif
	enddo

	return
	end

      SUBROUTINE TDF(RIW,cell,mxdt,nver,vx,vy,vz,tdfden,nv)

************************************************************************************
C     routine to calculate the total distribution function of voronoi vertex
************************************************************************************

      INTEGER*2 itdfdiv,i,j,k,NVER,nba,nez,MXDT,nv
	PARAMETER (itdfdiv=500)
	REAL*8 VXI,VYI,VZI,VX(MXDT),VY(MXDT),VZ(MXDT)
	REAL*8 cell(9),facpi,r,tdfden(itdfdiv),vc
	REAL*8 tdfcut,dis2bin,dr,vol(itdfdiv),dr2,r1,r2,distance
      REAL*8 trfactor,RIW(mxdt),DI


	tdfcut=8.0
	dis2bin=real(itdfdiv)/tdfcut
	dr=1/dis2bin
	dr2=0.5*dr
      facpi=4*3.1415926/3
      


	do j=1,itdfdiv
	
	tdfden(j)=0.0
	enddo

      nv=0
	do i=1,NVER-1 !3
	if(RIW(i).gt.0)then !4
      nv=nv+1
	do j=i+1,NVER !1
	if(RIW(j).gt.0)then !2
      
	VXI=VX(j)-VX(i)
	VYI=VY(j)-VY(i)
	VZI=VZ(j)-VZ(i)

	IF(abs(VXI).GT.cell(1)/2)VXI=VXI-SIGN(cell(1),VXI)
	IF(abs(VYI).GT.cell(5)/2)VYI=VYI-SIGN(cell(5),VYI)
	IF(abs(VZI).GT.cell(9)/2)VZI=VZI-SIGN(cell(9),VZI)

	DI=VXI**2+VYI**2+VZI**2
	nba=idnint(sqrt(DI)*dis2bin)
	nez=MIN(itdfdiv,nba)
	if(nez.eq.0)nez=1.
	tdfden(nez)=tdfden(nez)+1.0
	endif !2
	enddo !1
	endif !4
	enddo !3

	


	return
	end



      Subroutine volume(DT,VT,ratom2,VV,RX,RY,RZ,mxatms,cell,atomindex)

************************************************************
C     this program calculates the void volume inside a 
C     tetrahedron in which vertices are spheres with different radii
C     24/12/2004
***************************************************************

      integer*2 mxatms
      REAL*8 RXIJ,RYIJ,RZIJ
	REAL*8 PS(4),PX(4),PY(4),PZ(4),RA(4),PX2(3),PY2(3),PZ2(3)
      REAL*8 COS23,COS24,COS34,SEN23,SEN24,SEN34
	REAL*8 VO(4,4),VS(4),VI(6),M(3),CS(6),ACS(6)
	REAL*8 AREA(4),PI
	PARAMETER (PI=3.1411592654)
	REAL*8 R1,R2,DET,VT,VV,VSP,VIP
	REAL*8 dd,RD,T1,T2,cost,VL
      INTEGER*2 DT(4)
	INTEGER*2 i,j,k,kk,nc,n1,n
	REAL*8 RX(mxatms),RY(mxatms),RZ(mxatms),ratom2(11),cell(9)
	INTEGER*2 atomindex(mxatms)


	do j=1,4
      
	  RXIJ=RX(DT(j))-RX(DT(1))
	  RYIJ=RY(DT(j))-RY(DT(1))
	  RZIJ=RZ(DT(j))-RZ(DT(1))

	IF(abs(RXIJ).GT.cell(1)/2)RXIJ=RXIJ-SIGN(cell(1),RXIJ)
	IF(abs(RYIJ).GT.cell(5)/2)RYIJ=RYIJ-SIGN(cell(5),RYIJ)
	IF(abs(RZIJ).GT.cell(9)/2)RZIJ=RZIJ-SIGN(cell(9),RZIJ)
	  
	  PS(j)=RXIJ**2+RYIJ**2+RZIJ**2
	  PX(j)=RXIJ
	  PY(j)=RYIJ
	  PZ(j)=RZIJ
c     print*,PX(j),PY(j),PZ(j),PS(j)
      enddo

C     calcolus volume of the tetrahedron

C     calcolus dihedral angles

      nc=0

	do i=1,3
	n=0
	   do j=1,4
	   if(i.ne.j)then
 	   n=n+1
	   PX2(n)=PX(j)-PX(i)
         PY2(n)=PY(j)-PY(i)
	   PZ2(n)=PZ(j)-PZ(i)
	   M(n)=SQRT(PX2(n)**2+PY2(n)**2+PZ2(n)**2)
	   endif
         enddo

         COS23=(PX2(1)*PX2(2)+PY2(1)*PY2(2)+PZ2(1)*PZ2(2))/(M(1)*M(2))
	   COS24=(PX2(1)*PX2(3)+PY2(1)*PY2(3)+PZ2(1)*PZ2(3))/(M(1)*M(3))
	   COS34=(PX2(2)*PX2(3)+PY2(2)*PY2(3)+PZ2(2)*PZ2(3))/(M(2)*M(3))

	   SEN23=SQRT(1-COS23**2)
	   SEN24=SQRT(1-COS24**2)
	   SEN34=SQRT(1-COS34**2)

      IF(i.eq.1)THEN

	nc=nc+1
	CS(nc)=(COS34-COS24*COS23)/(SEN24*SEN23)  ! dihedral angle of the 1-2 edge
      
	nc=nc+1
	CS(nc)=(COS24-COS23*COS34)/(SEN34*SEN23)  ! 1-3

	nc=nc+1
	CS(nc)=(COS23-COS24*COS34)/(SEN24*SEN34)  ! 1-4

	ELSEIF(i.eq.2)THEN

      nc=nc+1
	CS(nc)=(COS24-COS23*COS34)/(SEN34*SEN23)  ! 2-3

	nc=nc+1
	CS(nc)=(COS23-COS24*COS34)/(SEN24*SEN34)  ! 2-4

	ELSEIF(i.eq.3)THEN

	nc=nc+1
	CS(nc)=(COS23-COS24*COS34)/(SEN24*SEN34)  ! 3-4

	ENDIF
	ENDDO !i



 	do j=1,6

	  ACS(j)=ACOS(CS(j))
c	print*,j,ACS(j),CS(j)

	enddo
	
*****************************************************************
C     the six dihedral angles are referred at the edges
C       1     =     1-2
C       2     =     1-3
C       3     =     1-4
C       4     =     2-3
C       5     =     2-4
C       6     =     3-4
*****************************************************************

*****************************************************************
C     now we have to calculate the solid angle at each vertex
C     tetrahedra vertex               dihedral angles necessary
C             1                                1,2,3
C             2                                1,4,5
C             3                                2,4,6
C             4                                3,5,6
*****************************************************************

C     Solid angles by GIRARD formula

      AREA(1)=ACS(1)+ACS(2)+ACS(3)-PI
	AREA(2)=ACS(1)+ACS(4)+ACS(5)-PI
	AREA(3)=ACS(2)+ACS(4)+ACS(6)-PI
	AREA(4)=ACS(3)+ACS(5)+ACS(6)-PI

C     calculus of the portion of spheres inside the tetrahedron

      do i=1,4

	  VS(i)=(AREA(i)*ratom2(atomindex(DT(i)))**3)/3

c	print*, i,AREA(i)

      enddo

C     intersection volumes between the particles (spheres)

      n1=0

	do i=1,3
	  do j=i+1,4

	  RXIJ=PX(j)-PX(i)
	  RYIJ=PY(j)-PY(i)
        RZIJ=PZ(j)-PZ(i)

	  dd= RXIJ**2+ RYIJ**2+ RZIJ**2

	  dd=SQRT(dd)  ! distance between atoms i-j

	  n1=n1+1
        
	  cost=PI/(12*dd)

	  R1=ratom2(atomindex(DT(i)))
	  R2=ratom2(atomindex(DT(j)))
	  RD=R1+R2

	IF(dd.LT.RD)THEN

	  T1=(R1+R2-dd)**2
	  T2=(dd**2 +2*dd*R2 -3*R2**2 +2*dd*R1 +6*R2*R1 -3*R1**2)
        
C     volume of the lens creates by the spheres intersection

        VL=cost*T1*T2

C     volume of the portion of lens inside the tetrahedron

        VI(n1)=VL*ACS(n1)/(2*PI)

c	print*,n1,VL,ACS(n1),VI(n1)

      ELSE

	  VI(n1)=0.

	ENDIF

	enddo !i
	enddo !j

C     total volume of the spheres inside the tetrahedron

      VSP=0.
	
	do k=1,4
	VSP=VSP+VS(k)
	enddo

C     total volume of the portions of lens inside the tetrahedron

	VIP=0.

	do j=1,6
	VIP=VIP+VI(j)
	enddo
	     
C     therefore, the void volume inside the tetrahedron is:

      VV=VT-VSP+VIP

	IF(VV.GT.VT)VV=0.
      if(VV.LT.0)VV=0.
c	write(34,*) VT,VSP,VIP,VV

	return
	end


	      
	      
