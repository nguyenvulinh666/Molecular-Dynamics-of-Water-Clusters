      module driver_module

c***********************************************************************
c     
c     dl_poly module for defining simulation driver routines
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.10
c     Exp
c     
c***********************************************************************
      
      use config_module
      use forces_module
      use nlist_builders_module
      use optimiser_module
      use temp_scalers_module
      
      implicit none
      
      integer, parameter :: mxpass=250
      
      contains

      subroutine molecular_dynamics
     x  (lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero,newlst,
     x  stropt,cycle,ltad,lsolva,lfree,lghost,idnode,imcon,
     x  keyfce,keyfld,keyshl,keystr,keytol,kmax1,kmax2,kmax3,multt,
     x  mxnode,natms,ngrp,nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,
     x  nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,ntptbp,
     x  ntpter,ntpvdw,ntshl,ntteth,ntcons,numrdf,nsolva,isolva,
     x  alpha,delr,dlrpot,drewd,elrc,engang,engbnd,engcpe,engdih,
     x  engfbp,engfld,enginv,engshl,engsrp,engtbp,engter,engtet,
     x  epsq,fmax,opttol,rctter,rcut,rcutfb,rcuttb,rprim,rvdw,shlke,
     x  engcfg,temp,tstep,virang,virbnd,vircpe,virdih,virfbp,virfld,
     x  virinv,virlrc,virmet,virshl,virsrp,virtbp,virter,virtet,volm,
     x  engmet,virtot)
      
c***********************************************************************
c     
c     dl_poly subroutine for controlling subroutine calls in a standard
c     molecular dynamics simulation
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c     wl
c     2008/12/23 10:29:12
c     1.10
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical lfcap,lgofr,lneut,lnsq,loglnk,loptim,lzeql,lzero
      logical newlst,stropt,cycle,ltad,lsolva,lfree,lghost
      
      integer idnode,imcon,keyfce,keyfld,keyshl,keytol,ntcons
      integer keystr,kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp
      integer nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,nsteql
      integer ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet
      integer ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,nsolva
      integer isolva
      
      real(8) alpha,delr,dlrpot,drewd,elrc,engang,engbnd
      real(8) engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp
      real(8) engtbp,engter,engtet,epsq,fmax,opttol,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp
      real(8) tstep,virang,virbnd,vircpe,virdih,virfbp
      real(8) virfld,virinv,virlrc,virmet,virshl,virsrp
      real(8) virtbp,virter,virtet,volm,engmet,virtot
      
c     construct verlet neighbour list
      
      call nlist_driver
     x  (newlst,lneut,lnsq,loglnk,ltad,natms,idnode,mxnode,imcon,
     x  nneut,keyfce,rcut,delr,tstep)
      
c     calculate atomic forces
      
      call force_manager
     x  (newlst,lneut,lnsq,lgofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x  lghost,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x  numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,
     x  ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,
     x  keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,
     x  isolva,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x  vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x  rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x  virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x  virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x  engcfg,fmax,temp)
      
c     frozen atoms option
      
      call freeze(natms)
      
c     structure optimisation
      
      if(loptim.or.lzero)then
        
        call optimisation_selector
     x    (loptim,stropt,lzero,idnode,mxnode,natms,imcon,ntcons,
     x    nscons,ngrp,ntfree,keystr,keytol,engcfg,tstep,opttol)
        
        if(stropt.and.idnode.eq.0)
     x    write(nrite,"(/,/,1x,'structure optimisation converged ',
     x    'at step ',i6,/,/,/,1x,120('-'))") nstep
        
        cycle=(cycle.and.(.not.stropt))
        
      endif
      
c     total virial (excluding constraint virial and c.o.m virial)
c     for npt routines     note: virsrp already includes virlrc
      
      virtot=vircpe+virsrp+virbnd+virtbp+virter+virfld+
     x  virang+virshl+virtet+virmet
      
      return
      end subroutine molecular_dynamics
      
      subroutine shell_relaxation
     x  (lfcap,lgofr,lneut,lnsq,loglnk,lzeql,newlst,ltad,lsolva,
     x  lfree,lghost,idnode,imcon,keyfce,keyfld,keyshl,
     x  kmax1,kmax2,kmax3,multt,mxnode,natms,nhko,nlatt,nneut,
     x  nospl,nstbgr,nstep,nsteql,ntangl,ntbond,ntdihd,ntinv,
     x  ntpfbp,ntpmet,ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,
     x  ntpmls,nsolva,isolva,alpha,delr,dlrpot,drewd,elrc,engang,
     x  engbnd,engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp,
     x  engtbp,engter,engtet,epsq,fmax,rctter,rcut,rcutfb,rcuttb,
     x  rprim,rvdw,shlke,engcfg,temp,tstep,virang,virbnd,vircpe,
     x  virdih,virfbp,virfld,virinv,virlrc,virmet,virshl,virsrp,
     x  virtbp,virter,virtet,volm,engmet,virtot,rlxtol,pass0,
     x  pass1,pass2)
      
c***********************************************************************
c     
c     dl_poly subroutine for controlling subroutine calls in a
c     relaxed shell molecular dynamics simulation
c
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2006
c     
c     wl
c     2008/12/23 10:29:12
c     1.10
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lfcap,lgofr,lneut,lnsq,loglnk,lzeql,ltad
      logical newlst,relaxed,shgofr,lsolva,lfree,lghost

      integer idnode,imcon,keyfce,keyfld,keyshl
      integer kmax1,kmax2,kmax3,multt,mxnode,natms
      integer nhko,nlatt,nneut,nospl,nstbgr,nstep,nsteql
      integer ntangl,ntbond,ntdihd,ntinv,ntpfbp,ntpmet
      integer ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf
      integer keyrlx,ntpmls,pass,nsolva,isolva

      real(8) alpha,delr,dlrpot,drewd,elrc,engang,engbnd
      real(8) engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp
      real(8) engtbp,engter,engtet,epsq,fmax,rctter
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp
      real(8) tstep,virang,virbnd,vircpe,virdih,virfbp
      real(8) virfld,virinv,virlrc,virmet,virshl,virsrp
      real(8) virtbp,virter,virtet,volm,engmet,virtot,rlxtol
      real(8) pass0,pass1,pass2
      
      pass=0
      keyrlx=0
      shgofr=lgofr
      relaxed=.false.
      
      do while(.not.relaxed.and.pass.le.mxpass)
        
c     construct verlet neighbour list
        
        call nlist_driver
     x    (newlst,lneut,lnsq,loglnk,ltad,natms,idnode,mxnode,imcon,
     x    nneut,keyfce,rcut,delr,tstep)
        
c     calculate atomic forces
        
        call force_manager
     x    (newlst,lneut,lnsq,shgofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x    lghost,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x    numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,
     x    ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,
     x    keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,
     x    isolva,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x    vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x    rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x    virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x    virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x    engcfg,fmax,temp)
        
c     frozen atoms option
        
        call freeze(natms)
        
c     total virial (excluding constraint virial and c.o.m virial)
c     for npt routines     note: virsrp already includes virlrc
        
        virtot=vircpe+virsrp+virbnd+virtbp+virter+virfld+
     x    virang+virshl+virtet+virmet
        
c     relaxed shell option
        
        call relax_shells
     x    (relaxed,keyrlx,idnode,mxnode,natms,ntpmls,tstep,
     x    rlxtol)
        
        if(relaxed)then
          
          pass1=pass0*pass1
          pass0=pass0+1.d0
          pass1=pass1/pass0+pass/pass0
          pass2=max(dble(pass),pass2)
          
        endif
        
        pass=pass+1
        if(pass.gt.mxpass)call error(idnode,1950)
        shgofr=.false.
        
c     end of shell relaxation
        
      enddo
      
      return
      end subroutine shell_relaxation

      subroutine minimiser
     x  (lfcap,lneut,lnsq,loglnk,lzeql,newlst,idnode,imcon,keyfce,
     x  keyfld,keyshl,keytol,kmax1,kmax2,kmax3,multt,mxnode,natms,
     x  ngrp,nhko,nlatt,nneut,nospl,nscons,ntcons,nstbgr,nstep,
     x  nsteql,ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet,
     x  ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,alpha,delr,dlrpot,
     x  drewd,elrc,engang,engbnd,engcpe,engdih,engfbp,engfld,enginv,
     x  engshl,engsrp,engtbp,engter,engtet,epsq,fmax,opttol,rctter,
     x  rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp,tstep,
     x  virang,virbnd,vircpe,virdih,virfbp,virfld,virinv,virlrc,
     x  virmet,virshl,virsrp,virtbp,virter,virtet,volm,engmet,
     x  virtot,sigma,tolnce,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine for controlling subroutine calls in a
c     minimisation simulation
c
c     copyright - daresbury laboratory
c     author    - w. smith    may 2007
c     
c     wl
c     2008/12/23 10:29:12
c     1.10
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lfcap,lneut,lnsq,loglnk,lzeql,newlst,stropt,shgofr
      logical conopt,newjob,ltad,lsolva,lfree,lghost
      
      integer idnode,imcon,keyfce,keyfld,keyshl,keystr,pass,i
      integer kmax1,kmax2,kmax3,multt,mxnode,natms,ngrp,keytol
      integer nhko,nlatt,nneut,nospl,nscons,nstbgr,nstep,nsteql
      integer ntangl,ntbond,ntdihd,ntfree,ntinv,ntpfbp,ntpmet
      integer ntptbp,ntpter,ntpvdw,ntshl,ntteth,numrdf,ntcons
      integer fail,nsolva,isolva

      real(8) alpha,delr,dlrpot,drewd,elrc,engang,engbnd
      real(8) engcpe,engdih,engfbp,engfld,enginv,engshl,engsrp
      real(8) engtbp,engter,engtet,epsq,fmax,opttol,rctter,sigma
      real(8) rcut,rcutfb,rcuttb,rprim,rvdw,shlke,engcfg,temp
      real(8) tstep,virang,virbnd,vircpe,virdih,virfbp
      real(8) virfld,virinv,virlrc,virmet,virshl,virsrp,tolnce
      real(8) virtbp,virter,virtet,volm,engmet,virtot,engcon
      real(8) cfgmin,engunit,hnorm,grad0,grad1,ff1,sgn
      
      real(8), allocatable :: sxx(:),syy(:),szz(:)
      
      save cfgmin,newjob
      
      data newjob/.true./
      
      pass=0
      keystr=0
      shgofr=.false.
      stropt=.false.
      
c     dummy variables
      
      ltad=.false.
      lsolva=.false.
      lfree=.false.
      lghost=.false.
      nsolva=0
      isolva=1
      
c$$$c     diagnostic printing (not usually active)
c$$$
c$$$      if(idnode.eq.0)then
c$$$        
c$$$        write(nrite,"(1x,120('-'),
c$$$     x    /,/,1x,'    pass',5x,'eng_cfg',5x,'eng_vdw',5x,'eng_cou',
c$$$     x    5x,'eng_bnd',5x,'eng_ang',5x,'eng_dih',5x,'eng_tet',
c$$$     x    5x,'eng_met',/,1x,120('-'))")
c$$$        
c$$$      endif
      
      do while(.not.stropt.and.pass.lt.mxpass)
        
        pass=pass+1
        
c     construct verlet neighbour list
        
        call nlist_driver
     x    (newlst,lneut,lnsq,loglnk,ltad,natms,idnode,mxnode,imcon,
     x    nneut,keyfce,rcut,delr,tstep)
        
c     calculate atomic forces
        
        call force_manager
     x    (newlst,lneut,lnsq,shgofr,lzeql,loglnk,lfcap,lsolva,lfree,
     x    lghost,idnode,mxnode,natms,imcon,nstep,nstbgr,nsteql,
     x    numrdf,keyfce,kmax1,kmax2,kmax3,nhko,nlatt,ntpvdw,
     x    ntpmet,nospl,multt,nneut,ntptbp,ntpfbp,ntpter,keyshl,
     x    keyfld,ntbond,ntangl,ntdihd,ntinv,ntteth,ntshl,nsolva,
     x    isolva,delr,dlrpot,engcpe,engsrp,epsq,rcut,rprim,rvdw,
     x    vircpe,virsrp,alpha,drewd,volm,engmet,virmet,elrc,virlrc,
     x    rcuttb,engtbp,virtbp,rcutfb,engfbp,virfbp,rctter,engter,
     x    virter,engbnd,virbnd,engang,virang,engdih,virdih,enginv,
     x    virinv,engtet,virtet,engshl,shlke,virshl,engfld,virfld,
     x    engcfg,fmax,temp)
        
c     frozen atoms option
        
        call freeze(natms)
        
c     total virial (excluding constraint virial and c.o.m virial)
c     for npt routines     note: virsrp already includes virlrc
        
        virtot=vircpe+virsrp+virbnd+virtbp+virter+virfld+
     x    virang+virshl+virtet+virmet
        
c     conjugate gradient structure optimisation
        
        call strucopt
     x    (stropt,keystr,keytol,idnode,mxnode,natms,ntcons,nscons,
     x    imcon,ngrp,ntfree,tstep,opttol,engcfg,hnorm,grad0,grad1,
     x    ff1,sgn)
        
c$$$c     diagnostic printing (not usually active)
c$$$
c$$$        if(idnode.eq.0)then
c$$$          
c$$$          write(nrite,"(1x,i8,1p,8e12.4)")
c$$$     x      pass,engcfg/engunit,engsrp/engunit,engcpe/engunit,
c$$$     x      engbnd/engunit,engang/engunit,engdih/engunit,engtet/
c$$$     x      engunit,engmet/engunit
c$$$          write(nrite,"(1x,120('-'))")
c$$$          
c$$$        endif

c     end of structure minimisation
        
      enddo
      
c     ensure constraints are satisfied
      
      if(stropt.and.ntcons.gt.0)then
        
        allocate(sxx(mxatms),syy(mxatms),szz(mxatms),stat=fail)
        if(fail.ne.0)call error(idnode,9999)
        
c     store current forces
        
        do i=1,natms
          
          sxx(i)=fxx(i)
          syy(i)=fyy(i)
          szz(i)=fzz(i)
          
        enddo
        
        keystr=0
        conopt=.false.
        
        do while(.not.conopt.and.pass.lt.mxpass)

          pass=pass+1
          engcon=0.d0
          
          do i=1,natms
            
            fxx(i)=0.d0
            fyy(i)=0.d0
            fzz(i)=0.d0
            
          enddo
          
c     conjugate gradient structure optimisation of constraint bonds
          
          call strucopt
     x      (conopt,keystr,keytol,idnode,mxnode,natms,ntcons,nscons,
     x      imcon,ngrp,ntfree,tstep,opttol,engcon,hnorm,grad0,grad1,
     x      ff1,sgn)
          
        enddo
        
c     restore current forces
        
        do i=1,natms
          
          fxx(i)=sxx(i)
          fyy(i)=syy(i)
          fzz(i)=szz(i)
          
        enddo
        
        deallocate(sxx,syy,szz,stat=fail)
        
      endif
      
c     write data summary
      
      if(idnode.eq.0)then
        
        if(stropt)then
          
          write(nrite,'(1x,"minimisation converged after ",i6," cycles"
     x      ," energy minimum: ",1pe12.4)')pass,engcfg/engunit
          
        else
          
          write(nrite,'(1x,"minimisation NOT converged after ",i6,
     x      " cycles")')pass
          
        endif
        
        write(nrite,"(1x,120('-'))")
        
      endif
      
c     reset velocities after structure optimisation
      
      call regauss(idnode,imcon,mxnode,natms,ngrp,nscons,ntcons,
     x  ntshl,keyshl,sigma,temp,tolnce)
      
c     write out minimised structure if lowest obtained so far
      
      if(newjob.or.cfgmin.gt.engcfg)then
        
        if(idnode.eq.0)call config_write('CFGMIN',0,imcon,natms,engcfg)
        cfgmin=engcfg
        newjob=.false.
        
      endif
      
      return
      end subroutine minimiser
      
      end module driver_module
