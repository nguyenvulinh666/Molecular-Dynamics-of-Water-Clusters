      program cluster

**********************************************************************************
*    
*     Copiright-computer modelling lab., University of Modena and Reggio Emilia.
*     Author- Dr. Pedone Alfonso Jenuary 2005
*
*     The cluster program is available free of charge to academic institutions. 
*     Copies should be obtained from the author or with the DL_POLY package.
*
*     No claim is made that the program is free from errors, the user is responsible 
*     for checking the validity of their results.
*
*     Details of the algorithm and its application can be found into the paper: 
*     "A tool for the prediction of crystalline phases obtained from controlled 
*     crystallization of glasses".
*     G.Lusvardi,G.Malavasi,L.Menabue,M.C.Menziani,A.Pedone and U.Segre
*     J.Phys.Chem., web Release Date: 22-oct-2005 DOI: 10.1021/jp0546857
* 
*     Aim of the program:
*
*     This is a FORTRAN interactive program able to search cluster of ions 
*     (crystallization nucleus) with stoichiometry resembling a the crystal phase that
*     you are looking for inside the structure of glasses modelled by MD simulation.
*
*     The source data is assumed to be formatted and compatible with the CONFIG and 
*     REVCON file written by the subroutine REVIVE.f of the DL_POLY package. 
*
*     An explanation of the CONFIG file format is available on line to the following 
*     link: http://www.cse.clrc.ac.uk/msi/software/DL_POLY/MANUALS/USRMAN3/node113.html
*
*     However, a input file for the cluster program is attached with the program.
*
*     After you run the executable the program ask you the following questions:
*
*     1.Name of the REVCON or CONFIG input file you want to analyze.
*     2.Name of the output file 
*     3.The stoichiometry of the crystal phase you are looking for.
*
*     For example, if you want to search a crystallization nucleus of the phase NaCaPO4
*     inside a multi-component glass with O,Si,Na,Ca and P ions, the program will ask you:
*
*     >>Insert number of O ions
*     >>4.0           (this should be your answer)
*     >>Insert number of Si ions
*     >>0.0           (the crystal phase does not contain Si atoms)
*     >>Insert number of Na ions
*     >>1.0
*     and so on.
*
*
*     OUTPUT file:
*
*     Record 1: Number of the atom that lies at the centre of the spherical cluster.
*     Record 2-3: Cartesian coordinates of the central atom
*     Record 4: Radius of the cluster
*     Record 5-n:It is listed the composition and minimum stoichiometry of the cluster 
*     with the smaller displacement found.
*     Record n+1:displacement of the cluster
*
*     The subsequent records contain the same information for bigger clusters.  
*
*     15/09/2005
***********************************************************************************

      integer mxatms,mxspec,numspec,natms
	integer imcon,boundary,nstep,nary
	integer mxatnump,i,j,k,n,iatom2
	parameter (mxatms=5000,mxspec=11)
	parameter (mxatnump=mxatms*(mxatms+1)/2)
	character*8 speclab(mxspec)
      integer mmbas(mxspec),ltpbas(mxatms)
      integer ii,jj,nhis
	character*80 title,fname,fileout
	character*8 atmnam(mxatms)
	real*8 xxx(mxatms),yyy(mxatms),zzz(mxatms)
	real*8 xdf(mxatnump), ydf(mxatnump), zdf(mxatnump)
      real*8 rr2(mxatms*(mxatms+1)/2)
	real*8 vxx,vyy,vzz,fxx,fyy,fzz
	real*8 cell(9),phstoy(mxspec),clstoy(mxspec),clstoymin(mxspec)
	real*8 natmin(mxspec),natom(mxspec)
	real*8 tstep,dmin,rcut,rcut2,rmin,sum,rmin2
      
	logical new


c-----------------------------------------------------------------
c     read revcon input file
c-----------------------------------------------------------------
	
	write(*,'(a)')'Enter Input file name:'
	read(5,'(a)')fname
	write(*,'(a)')'Enter output file name:'
      read(5,'(a)')fileout
	
	
	


      open(33,file=fname)
      read(33,'(a80)')title
      read(33,*)imcon,boundary,nstep,tstep
      read(33,'(3(6x,f14.10))')cell(1),cell(2),cell(3) 
      read(33,'(3(6x,f14.10))')cell(4),cell(5),cell(6)      
      read(33,'(3(6x,f14.10))')cell(7),cell(8),cell(9)

	if(imcon.eq.0)then
	do i=1,mxatms
      read(33,'(a4,i18)',end=112)atmnam(i),n
      read(33,'(3(f20.8))')xxx(i),yyy(i),zzz(i)
	enddo
      elseif(imcon.eq.1)then
      do i=1,mxatms
      read(33,'(a4,i18)',end=112)atmnam(i),n
      read(33,'(3(f20.8))')xxx(i),yyy(i),zzz(i)
      read(33,'(3(f20.8))')vxx,vyy,vzz
      enddo
	elseif(imcon.eq.2)then      
      do i=1,mxatms
      read(33,'(a4,i18)',end=112)atmnam(i),n
      read(33,'(3(f20.8))')xxx(i),yyy(i),zzz(i)
      read(33,'(3(f20.8))')vxx,vyy,vzz
      read(33,'(3(f20.8))')fxx,fyy,fzz
      enddo
	endif
	close(33)
     
  112 natms=i-1

C------------------------------------------------------------------
C     get statistic of input atoms
C------------------------------------------------------------------

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

	write(*,*)'ENTER STOICHIOMETRY OF THE CRYSTAL PHASE :'
	do j=1,numspec
	write(*,"('INSERT NUMBER OF ',a4,'IONS')")speclab(j)
	read(5,*)phstoy(j)
	enddo

c---------------------------------------------------------------------
C     initialization
c---------------------------------------------------------------------

c     distance between each pair of atoms

      do i=1,natms
      do j=i,natms
      k=j*(j-1)/2+i
      xdf(k)=xxx(i)-xxx(j)
      ydf(k)=yyy(i)-yyy(j)
      zdf(k)=zzz(i)-zzz(j)

      enddo
      enddo
      nary=natms*(natms+1)/2
      
c     periodic boundary conditions

      call images(imcon,0,1,nary,cell,xdf,ydf,zdf)

c     square-distances between atoms

      do i=1,nary
      rr2(i)=xdf(i)*xdf(i)+ydf(i)*ydf(i)+zdf(i)*zdf(i)
      enddo

c-----------------------------------------------------------------------

      rcut=3.0
      n=0
      open(34,file=fileout)
      do while(rcut.le.8.0)
      n=n+1

	rcut2=rcut**2

      rmin2=1000.0
      do i=1,natms
	if(atmnam(i).eq.'O2-      ')then
	
	do ii=1,numspec
      clstoy(ii)=0.0
	enddo
	
	do j=1,natms
	k=idx(i,j)
	if(rr2(k).lt.rcut2)then
	do ii=1,numspec
	if(atmnam(j).eq.speclab(ii))then
	clstoy(ii)=clstoy(ii)+1.
	endif
	enddo !ii
	endif
	enddo !j

c******************************************************************************
c     you can choose to normalize the minimum stoichiometry of
c     the cluster respect the oxygen.
c	do ii=1,numspec
c	if(speclab(ii).eq.'O2-')then
c	dmin=clstoy(ii)/phstoy(ii)
c	endif
c	enddo
c******************************************************************************

c******************************************************************************
c     cluster minimum stoichiometry normalized respect the number of the ions species
c     present in smaller amount


      dmin=100.
	do ii=1,numspec
	if(phstoy(ii).ne.0)then
	if(clstoy(ii).lt.dmin.and.clstoy(ii).ne.0.0)dmin=clstoy(ii)
	endif
	enddo
*******************************************************************************

	do ii=1,numspec
	clstoymin(ii)=clstoy(ii)/dmin
	enddo

	sum=0.0
	

	do jj=1,numspec
	sum=sum+abs((clstoymin(jj)-phstoy(jj)))
	enddo
  

	rmin=sum

c     we store only the composition of the best cluster	
	if(rmin.lt.rmin2)then
	rmin2=rmin
	iatom2=i
	do k=1,numspec
	natom(k)=clstoy(k)
	natmin(k)=clstoymin(k)
	enddo !k
	endif

	endif
	enddo !i
     
C     print composition of the cluster with minimum displacement
C     at the radius r     

	write(34,*)'Number of the central atom :',iatom2
	write(34,*)'coordinates ','    x    ','    y    ','    z    '
	write(34,'(10x,3(2x,f8.4))')xxx(iatom2),yyy(iatom2),zzz(iatom2)
      write(34,"('Cutoff : ',f6.3)")rcut
	
	write(34,'(a)')' Atoms  Composition  Min.Stoychiometry :'
      do k=1,numspec
        write(34,'(a8,3x,f8.3,f8.3)')speclab(k),natom(k),natmin(k)
      end do !k
	write(34,*)'displacement=',rmin2
	write(34,"(40('-'))")

      rcut=rcut+0.5
	enddo !while


	stop
	end

******************************************************************************
C     starting subroutine
******************************************************************************

      function idx(i,j)
      
	integer i,j 

      if(i.gt.j) then
        idx=i*(i-1)/2+j
      else
        idx=j*(j-1)/2+i
      endif

      return

      end

      subroutine images
     x  (imcon,idnode,mxnode,natms,cell,xxx,yyy,zzz)
      
c     
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
c     
c     note: in all cases the centre of the cell is at (0,0,0)
c     warning - replicated data version: does not re-merge 
c     coordinate arrays
c     
c     wl
c     1996/02/15 14:32:57
c     1.1.1.1
c     Exp
c     
c***********************************************************************
c     
      
      real*8 xxx(*),yyy(*),zzz(*)
      real*8 cell(9),rcell(9)
	integer incom,idnode,iatm1,iatm2
	integer mxnode,natms,i,imcon
	real*8 aaa,bbb,ccc,det,rt2
	real*8 xss,yss,zss,ssx,ssy,ssz

      if(imcon.gt.0) then

c     
c     block indices

        iatm1 = (idnode*natms)/mxnode+1
        iatm2 = ((idnode+1)*natms)/mxnode

      endif
      
      if(imcon.eq.1)then
c     
c     standard cubic boundary conditions
        
        
        aaa=1.d0/cell(1)


        do i=iatm1,iatm2
          xxx(i)=xxx(i)-cell(1)*nint(aaa*xxx(i))
          yyy(i)=yyy(i)-cell(1)*nint(aaa*yyy(i))
          zzz(i)=zzz(i)-cell(1)*nint(aaa*zzz(i))
        enddo
        
      else if(imcon.eq.2)then
c     
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
c     
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
c     
c     truncated octahedral boundary conditions
        
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(5)-cell(9)).lt.1.d-6)) then
          print *,'error-130'
          stop
        endif
        
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
c     
c     rhombic dodecahedral boundary conditions
        
        rt2=sqrt(2.d0)
        if(.not.(abs(cell(1)-cell(5)).lt.1.d-6.and.
     x    abs(cell(9)-cell(1)*rt2).lt.1.d-6))then 
          print *,'error-140'
          stop
        endif
        
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
c     
c     x-y boundary conditions 

        det = cell(1)*cell(5) - cell(2)*cell(4)

        if(abs(det).lt.1.d-6) then
          print *,'error-120'
          stop
        endif
        
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

      endif

      return
      end


      subroutine invert(a,b,d)
c     
c***********************************************************************
c     
c     dl_poly subroutine to invert a 3 * 3 matrix using cofactors
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith       april 1992
c     
c     wl
c     1996/02/15 14:32:58
c     1.1.1.1
c     Exp
c     
c***********************************************************************
c     
      
      real*8 a(9),b(9)
	real*8 d,r
c     
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
c     
c     calculate determinant
      d=a(1)*b(1)+a(4)*b(2)+a(7)*b(3)
      r=0.d0
      if(abs(d).gt.0.d0)r=1.d0/d
c     
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
      end

      subroutine getname(string,length,lout,mark)

c***********************************************************************
c
c     DL_POLY subroutine to get characters without blanks from a string
c     return total number of characters by variable lout
c
c     september 1996   - x. yuan
c
c***********************************************************************

      character*(*) string
      character*1 mark
      integer lout,length

      call strip(string,length)

      lout=0

      do while(string(lout+1:lout+1).ne.mark)

        lout=lout+1

      enddo

      return

      end

      subroutine strip(string,length)

c***********************************************************************
c     
c     DL_POLY routine to strip blanks from start of a string
c     maximum length is 255 characters
c     
c     copyright daresbury laboratory 1993
c     author   t.forester       july 1993
c     
c     wl
c     1996/02/15 14:33:22
c     1.1.1.1
c     Exp
c     
c***********************************************************************

      character*(*) string
	integer imax,i,j,length
      
      imax = min(length,255)
      do i = 1,imax

        if(string(1:1).eq.' ') then

          do j = 1,imax-1

            string(j:j) = string(j+1:j+1)

          enddo

          string(imax:imax) = ' '

        endif

      enddo

      return
      end
	


      