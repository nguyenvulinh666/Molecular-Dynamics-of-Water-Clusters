      ! SIMPLIFIED STOCHASTIC-THERMOSTATS ROUTINE, 
      ! only global rescaling, langevin and GLE

      subroutine stochastic(keystochastic,idnode,mxnode,imcon,natms,
     x     ngrp,ntshl,temp,tstep,taut,taus,taup,taub,chip,keyens,
     x     degfre,lgdof,ascale,engke,engrot,engstoc)

      use setup_module
      use config_module
      use rigid_body_module
      use core_shell_module
      use ensemble_tools_module

      implicit none

      integer, intent(in)    :: keystochastic                  ! switch for various stochastic thermostats
      integer, intent(in)    :: idnode,mxnode,imcon,natms,ngrp ! atoms and nodes indexes
      real*8,  intent(in)    :: temp                           ! target 
      real*8,  intent(in)    :: degfre                         ! number of degrees of freedom, keeping in account frozen
                                                               ! atoms and bonds
      integer, intent(in)    :: lgdof,keyens                   ! number of stochastic DOF for generalized colored langevin
      real*8,  intent(in)    :: ascale                         ! GLE scaling factor for Ap
      real*8,  intent(in)    :: taut                           ! thermostat time-scale
      real*8,  intent(in)    :: taus                           ! shell-scaling time
      real*8,  intent(in)    :: taup                           ! barostat timescale 
      real*8,  intent(in)    :: taub                           ! barostat thermostat time
      real*8,  intent(inout) :: chip                           ! barostat velocity
      real*8,  intent(in)    :: tstep   	                     ! time step 
      real*8,  intent(inout) :: engke, engrot,engstoc          ! time step 
      

      integer, intent(in)    :: ntshl
      integer iatm1,iatm2,i,j,mshl, ios

      real*8 sigma, kt
      real*8 scale,scale1,scale2 
      real*8 resamplekin
!      real*8 :: getkin,gasdev,resamplekin,getkinr
      real*8 pmass,ekinbar
      real*8 sumke0,sumke,twork
      real*8 vx,vy,vz
      real(8), save, allocatable :: fixcell(:)

c     GLE propagator matrix & temporaries
      real*8, allocatable, save :: gT(:,:), gS(:,:)
      
      logical :: first = .true.

      real*8 :: vcx,vcy,vcz,wwx,wwy,wwz,wcm, wi,wj
      integer :: idof

CVAM
CVAM      call VTBEGIN(64, ierr)
CVAM
c     block indices
      iatm1=(idnode*natms)/mxnode+1
      iatm2=((idnode+1)*natms)/mxnode
      idof=int(degfre)

      if (.not. allocated(fixcell)) then
        allocate(fixcell(9))
        open(666,file='FIXBOX',status='OLD',iostat=ios)
        fixcell=1
        if (ios.eq.0) then
          read(666,*) fixcell(:)         
        end if
        close(666)
      endif

      kt=temp*boltz
      sigma=kt*degfre*0.5d0
      
c     kinetic energy of barostat  
      ekinbar=0.0    
      pmass=2.d0*sigma*taup**2
      if (keyens.eq.5) then
        ekinbar=ekinbar+chip**2
      else if (keyens.eq.7) then
        do i=1,9
          ekinbar=ekinbar+eta(i)**2
        enddo
      endif
      ekinbar=ekinbar*0.5d0*pmass
c     total kinetic energy, before thermostat
c     we recompute, as with new Trotter decomposition I am not sure of what
c     we are getting when the routine is called.
      sumke0=getkin(natms,idnode,mxnode)+ekinbar
      twork=sumke0

      
      select case(keystochastic)
      case(0)
c	**********************************************
c       *        stochastic velocity rescale         *
c	**********************************************
        if(sumke0.gt.1.d-6) then
          sumke=0.0d0
          if(idnode.eq.0) then
            sumke=resamplekin(sumke0,sigma,idof,taut/tstep)
          end if
          if(mxnode.gt.1) call gdsum(sumke,1,buffer)
          scale=sqrt(sumke/sumke0)
        else
          scale=1.0
          sumke=sumke0
        end if

c       rescale the velocities
!!p        do i=1,natms
        do i=iatm1,iatm2
          vxx(i)=scale*vxx(i)
          vyy(i)=scale*vyy(i)
          vzz(i)=scale*vzz(i)
        enddo

      case(1)
c	**********************************************
c       *        standard langevin thermostat        *
c	**********************************************

        if(taut/tstep>0.1) then
          scale1 = exp(-tstep/taut)
          scale2 = sqrt(kT*(1.0-exp(-2.0*tstep/taut)))
        else
          scale1 = 0.0
          scale2 = sqrt(kT)
        end if

        do i=iatm1,iatm2
            vxx(i)=scale1*vxx(i)+scale2/sqrt(weight(i))*gasdev()
            vyy(i)=scale1*vyy(i)+scale2/sqrt(weight(i))*gasdev()
            vzz(i)=scale1*vzz(i)+scale2/sqrt(weight(i))*gasdev()
        end do

      case(2)
c	**********************************************
c       *        colored  langevin thermostat        *
c	**********************************************
c     first, let's compute the propagator
        if (.not. allocated(gT)) then
          allocate(gS(lgdof+1,lgdof+1))
          allocate(gT(lgdof+1,lgdof+1))

          call general_matrices(lgdof,tstep,ascale,kT,gT,gS)
c          write(0,*) gT
c          write(0,*) gS
c     if one wants to print out the work of the thermostat
!          open(242,file='THERMOWORK')
        end if
c     then we scale velocities to get >>mass-scaled momenta<<
         do i=iatm1,iatm2
            if (weight(i).gt.0) then ! frakkin' massless charges!
              wi=sqrt(weight(i))
              vxx(i)=vxx(i)*wi
              vyy(i)=vyy(i)*wi
              vzz(i)=vzz(i)*wi
            end if 
        end do

C shell-model removed for a speedup!

        !TENTATIVE PARALLEL IMPLEMENTATION
!        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,grxx(iatm1,1),vxx(iatm1))
!        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,gryy(iatm1,1),vyy(iatm1))
!        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,grzz(iatm1,1),vzz(iatm1))

        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,grxx(1,iatm1),vxx(iatm1))
        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,gryy(1,iatm1),vyy(iatm1))
        call gs_evolve(lgdof,iatm2-iatm1,gS,gT,grzz(1,iatm1),vzz(iatm1))
         
!        do i=iatm1,iatm2
!           call gs_evolve(gS,gT,lgdof,grxx(i,:),vxx(i))
!           call gs_evolve(gS,gT,lgdof,gryy(i,:),vyy(i))
!           call gs_evolve(gS,gT,lgdof,grzz(i,:),vzz(i))
!        end do
            
        ! now we project out the change in aux. velocities for rigid bodies
        if(ngrp.gt.0) then
        do i=iatm1,iatm2 !from ms momenta to velocities, taking care of massless charges
          if (weight(i).eq.0) then
            wi=0
          else 
            wi=1./sqrt(weight(i))
          endif
          grxx(:,i)=grxx(:,i)*wi
          gryy(:,i)=gryy(:,i)*wi
          grzz(:,i)=grzz(:,i)*wi  
        enddo
        do i=1,lgdof
        !this merge seems to be really necessary, as rb atoms need not be stored on same node
          if(mxnode.gt.1) call merge(idnode,mxnode,natms,mxbuff,
     x          grxx(i,:),gryy(i,:),grzz(i,:),buffer)
          scale=1.
          call vsimplescaleg(idnode,mxnode,imcon,natms,
     x                ngrp,scale,grxx(i,:),gryy(i,:),grzz(i,:))
        enddo
        do i=iatm1,iatm2 !and back to ms momenta
          wi=sqrt(weight(i))
          grxx(:,i)=grxx(:,i)*wi
          gryy(:,i)=gryy(:,i)*wi
          grzz(:,i)=grzz(:,i)*wi  
        enddo
        endif
            
c     then we scale back to velocities
         do i=iatm1,iatm2
            wi=1./sqrt(weight(i))
            vxx(i)=vxx(i)*wi
            vyy(i)=vyy(i)*wi
            vzz(i)=vzz(i)*wi
        end do
      end select
      if(mxnode.gt.1) then
          call merge(idnode,mxnode,natms,mxbuff,vxx,vyy,vzz,buffer)
      end if    
      
c     this is to update quaternions and remove com motion and angular momentum
      scale=1.0
      call vsimplescaleg(idnode,mxnode,imcon,natms,ngrp,scale, 
     x                   vxx, vyy, vzz)

c     if requested, apply a langevin piston thermostat to the cell's DOF
      if(taub>0.0) then
        if(taub/tstep>0.1) then
          scale1 = exp(-tstep/taub)
          scale2 = sqrt(kT*(1.0-exp(-2.0*tstep/taub)))
        else
          scale1 = 0.0
          scale2 = sqrt(kT)
        end if
        if (keyens.eq.5) then
          if(idnode.eq.0) then
            chip=scale1*chip+scale2/sqrt(pmass)*gasdev()
          else
            chip=0.0
          end if
          call gdsum(chip,1,buffer)
        else if (keyens.eq.7) then
          do i=1,9
          if(idnode.eq.0) then
            eta(i)=scale1*eta(i)+fixcell(i)*scale2/sqrt(pmass)*gasdev()
          else
            eta(i)=0.0
          end if
          enddo
          call gdsum(eta,9,buffer)
        end if
      end if


c     total kinetic energy, after thermostat
      ekinbar=0.0    
      if (keyens.eq.5) then
        ekinbar=ekinbar+chip**2
      else if (keyens.eq.7) then
        do i=1,9
          ekinbar=ekinbar+eta(i)**2
        enddo
      endif
      ekinbar=ekinbar*0.5d0*pmass
      sumke=getkin(natms,idnode,mxnode)+ekinbar

      twork=twork-sumke
!      write(242,*) twork

      engstoc=engstoc+twork

c     rotational energy
      engrot=getkinr(ngrp,idnode,mxnode)

c     translational energy
      engke=sumke-engrot-ekinbar

CVAM
CVAM      call VTEND(64, ierr)
CVAM
      return
      end

      function resamplekin(kk,sigma,ndeg,taut)
        use utility_module
        implicit none
        real*8 :: resamplekin
        real*8, intent(in) :: kk,sigma,taut
        integer, intent(in) :: ndeg

        real*8 :: fact,rr

        if(taut>0.1) then
          fact=exp(-1.0/taut)
        else
          fact=0.0
        end if
        rr = gasdev()

        resamplekin = kk + (1.0-fact)*
     x                  (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk)
     x              + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-fact)*fact)

      end function resamplekin

      subroutine general_matrices(lgdof, tstep, ascale, kT, gT, gS)
      use setup_module
      implicit none
      integer, intent(in) :: lgdof
      real*8, intent(in)  :: tstep, kT, ascale
      real*8, intent(out) :: gT(lgdof+1,lgdof+1),gS(lgdof+1,lgdof+1)
      
      integer :: i,j,ios
      real*8 :: gA(lgdof+1,lgdof+1)

      !read in A matrix
      open(121,file='STOCH-A',status='OLD',iostat=ios)
      if (ios.ne.0) write(0,*) "Unable to read STOCH-A file for GLE"
      do i=1,lgdof+1
         read(121,*) gA(i,:)
      enddo
      close(121)
      gA=gA*ascale
      !now we have the deterministic propagator
      call matrix_exp(-tstep*gA, lgdof+1,15,15,gT)
            
      !if there is no STOCH-C file, assume a diagonal covariance matrix
      !we expect to have the covariance matrix in Kelvin!

      open(121,file='STOCH-C',status='OLD',iostat=ios)
      if (ios.ne.0) then
        gS=-matmul(gT,transpose(gT))
        do i=1,lgdof+1
           gS(i,i)=gS(i,i)+1
        enddo
        gS=gS*kT
      else
        do i=1,lgdof+1
          read(121,*) gS(i,:)
        enddo
        gS=gS*boltz ! now the covariance is in energy units...
        if (abs(gS(1,1)-kT).gt.1e-5) then
          write(nrite,*)  " WARNING! STOCH-C classical limit ",
     x    " mismatches with temperature in CONTROL."
          write(nrite,*) " Rescaling",gS(1,1)/boltz, " to ",kT/boltz
          gA=gA*kT/gS(1,1)
          gS=gS*kT/gS(1,1)
          !we have to recompute the deterministic propagator...
          call matrix_exp(-tstep*gA, lgdof+1,15,15,gT)
        end if
        gS=gS-matmul(gT,matmul(gS,transpose(gT)))
      end if
      close(121)

      !now we must cholesky-decompose: gA contains SS^t
      gA=gS
      call cholesky(gA, gS, lgdof+1)
      end subroutine
      
      subroutine matrix_exp(M, n, j, k, EM)
      integer, intent(in)  :: n, j, k
      real*8, intent(in)   :: M(n,n)
      real*8, intent(out)   :: EM(n,n)
      
      real *8 :: tc(j+1), tmp(n,n), SM(n,n)
      integer p, i
      tc(1)=1
      do i=1,j
         tc(i+1)=tc(i)/dble(i)
      enddo
      
      SM=M*(1./2.**k)
      EM=0.
      do i=1,n
         EM(i,i)=tc(j+1)
      enddo
      
      !taylor exp of scaled matrix
      do p=j,1,-1
         EM=matmul(SM,EM);
         do i=1,n
            EM(i,i)=EM(i,i)+tc(p)
         enddo
      enddo

      do p=1,k
         EM=matmul(EM,EM)
      enddo
      end subroutine

  ! brute-forcei stabilized cholesky decomposition
      subroutine cholesky(SST, S, n)
        integer, intent(in)  :: n
        real*8, intent(in)   :: SST(n,n)
        real*8, intent(out)   :: S(n,n)
        real*8 :: D(n), L(n,n)
        integer i,j,k
        D=0. 
        L=0.
        do i=1,n
           L(i,i)=1.
           D(i)=SST(i,i)
           do j=1,i-1
              L(i,j)=SST(i,j);
              do k=1,j-1
                 L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k)
              enddo
              if (D(j).ne.0.) L(i,j)=L(i,j)/D(j)
           enddo
           do k=1,i-1
              D(i)=D(i)-L(i,k)*L(i,k)*D(k)
           end do
        enddo
        S=0.
        do i=1,n
        do j=1,i
           if (D(j)<=0) write(nrite,*) "NEG VALUE IN CHOLESKY ",D(j)
           if (D(j)>0.) S(i,j)=S(i,j)+L(i,j)*sqrt(D(j))
        enddo
        enddo
      end subroutine cholesky

      subroutine gs_evolve(lgdof, nat, gS, gT, gr, lp)
      use utility_module 
      implicit none
      integer, intent(in)    :: lgdof, nat
      real*8, intent(in)     :: gS(lgdof+1,lgdof+1),gT(lgdof+1,lgdof+1)
      real*8, intent(inout)     :: gr(lgdof,nat), lp(nat)
      real*8, save, allocatable ::  geta(:,:), grp(:,:)    ! we save, so these are not mallocated at each call
      integer k, i

      if (.not. allocated(geta) ) then
        allocate(geta(lgdof+1,nat))
        allocate(grp(lgdof+1, nat))
      end if

      do i=1,nat
      do k=1,lgdof+1
         grp(k,i)=gasdev()
      end do
      enddo
      
!      geta=matmul(gS,grp)
      call dgemm('n','n',lgdof+1,nat,lgdof+1,1.0d0,gS,lgdof+1,
     x            grp,lgdof+1,0.0d0,geta,lgdof+1)

      do i=1,nat
      grp(1,i)=lp(i)
      do k=1,lgdof
         grp(k+1,i)=gr(k,i)
      end do
      enddo
      
!      geta=matmul(gT,grp)+geta
      call dgemm('n','n',lgdof+1,nat,lgdof+1,1.0d0,gT,lgdof+1,
     x           grp,lgdof+1,1.0d0,geta,lgdof+1)

      do i=1,nat
      lp(i)=geta(1,i)
      do k=1,lgdof
         gr(k,i)=geta(k+1,i)
      end do
      enddo
      end subroutine

      subroutine vsimplescaleg(idnode,mxnode,imcon,natms,ngrp,scale,
     x                        hxx,hyy,hzz)
c*********************************************************************
c     it is vscaleg, modified so as to simply rescale the velocities with a given factor
c     
c     dl_poly subroutine for scaling the velocity arrays to the
c     desired temperature
c     
c     zeroes angular momentum in non-periodic system.
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1992.
c     author - w.smith july 1992
c     amended - t.forester oct 1993
c     amended - t.forester dec 1994 : block data
c     
c     wl
c     2005/12/02 14:22:19
c     1.3
c     Exp
c     
c*********************************************************************
      
      use setup_module
      use config_module
      use rigid_body_module

      implicit none
      real *8, intent(inout)::hxx(natms), hyy(natms), hzz(natms)
      integer idnode,mxnode,imcon,natms,ngrp,ierr,iatm1,iatm2,i
      real*8 roti,rotinv,cmx,cmy,cmz,cmvx,cmvy,cmvz,sysmas
      real*8 amx,amy,amz,det,scale,rsq,wxx,wyy,wzz

      dimension roti(9),rotinv(9)

CVAM
CVAM      call VTBEGIN(64, ierr)
CVAM

c     block indices

      iatm1=(idnode*natms)/mxnode+1
      iatm2=((idnode+1)*natms)/mxnode

c     calculate centre of mass position and motion from the system
      
      cmx=0.d0
      cmy=0.d0
      cmz=0.d0
      cmvx=0.d0
      cmvy=0.d0
      cmvz=0.d0
      sysmas=0.d0
      
      do i=iatm1,iatm2
        
        if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6) then

          cmx=cmx+weight(i)*xxx(i)
          cmy=cmy+weight(i)*yyy(i)
          cmz=cmz+weight(i)*zzz(i)
          sysmas=sysmas+weight(i)
          cmvx=cmvx+hxx(i)*weight(i)
          cmvy=cmvy+hyy(i)*weight(i)
          cmvz=cmvz+hzz(i)*weight(i)
          
        endif

      enddo
      
      if(mxnode.gt.1) then
        buffer(8)=sysmas
        buffer(9)=cmx
        buffer(10)=cmy
        buffer(11)=cmz
        buffer(12)=cmvx
        buffer(13)=cmvy
        buffer(14)=cmvz
        call gdsum(buffer(8),7,buffer(1))
        sysmas= buffer(8) 
        cmx=buffer(9) 
        cmy=buffer(10) 
        cmz=buffer(11) 
        cmvx=buffer(12) 
        cmvy=buffer(13) 
        cmvz=buffer(14) 
      endif

      cmx=cmx/sysmas
      cmy=cmy/sysmas
      cmz=cmz/sysmas
      
      cmvx=cmvx/sysmas
      cmvy=cmvy/sysmas
      cmvz=cmvz/sysmas
      
c     remove centre of mass motion  
      
      do i=1,natms
        
        if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6) then

          hxx(i)=hxx(i)-cmvx
          hyy(i)=hyy(i)-cmvy
          hzz(i)=hzz(i)-cmvz
          
        else

          hxx(i)=0.d0
          hyy(i)=0.d0
          hzz(i)=0.d0

        endif

      enddo
      
c     zero angular momentum about centre of mass - non-periodic system
      
      if(imcon.eq.0) then

c     move to centre of mass origin
        
!!p        do i=1,natms
!!p          
!!p          xxx(i)=xxx(i)-cmx
!!p          yyy(i)=yyy(i)-cmy
!!p          zzz(i)=zzz(i)-cmz
!!p          
!!p        enddo
        
c     angular momentum accumulators
        
        amx=0.d0
        amy=0.d0
        amz=0.d0

c     rotational inertia accumulators
        
        do i=1,9
          roti(i)=0.d0
        enddo
        
        do i=iatm1,iatm2
          
          amx=amx+weight(i)*(yyy(i)*hzz(i)-zzz(i)*hyy(i))
          amy=amy+weight(i)*(zzz(i)*hxx(i)-xxx(i)*hzz(i))
          amz=amz+weight(i)*(xxx(i)*hyy(i)-yyy(i)*hxx(i))
          
          rsq=xxx(i)**2+yyy(i)**2+zzz(i)**2
          roti(1)=roti(1)+weight(i)*(xxx(i)*xxx(i)-rsq)
          roti(2)=roti(2)+weight(i)* xxx(i)*yyy(i)
          roti(3)=roti(3)+weight(i)* xxx(i)*zzz(i)
          roti(5)=roti(5)+weight(i)*(yyy(i)*yyy(i)-rsq)
          roti(6)=roti(6)+weight(i)* yyy(i)*zzz(i)
          roti(9)=roti(9)+weight(i)*(zzz(i)*zzz(i)-rsq)
          
        enddo

c     complete rotational inertia matrix
        
        roti(4)=roti(2)
        roti(7)=roti(3)
        roti(8)=roti(6)

c     globally sum

        if(mxnode.gt.1) then
          buffer(13)=amx
          buffer(14)=amy
          buffer(15)=amz
          do i=1,9
            buffer(15+i)=roti(i)
          enddo
          call gdsum(buffer(13),12,buffer(1))
          amx= buffer(13) 
          amy= buffer(14) 
          amz= buffer(15) 
          do i=1,9
            roti(i)=buffer(15+i)
          enddo
        endif

c     invert rotational inertia matrix
        
        call invert (roti,rotinv,det)

c     correction to angular velocity
        
        wxx=rotinv(1)*amx+rotinv(2)*amy+rotinv(3)*amz
        wyy=rotinv(4)*amx+rotinv(5)*amy+rotinv(6)*amz
        wzz=rotinv(7)*amx+rotinv(8)*amy+rotinv(9)*amz

c     correction to linear velocity
        
        do i=1,natms

          if(lstfrz(i).eq.0.and.weight(i).gt.1.d-6) then

            hxx(i)=hxx(i)+(wyy*zzz(i)-wzz*yyy(i))
            hyy(i)=hyy(i)+(wzz*xxx(i)-wxx*zzz(i))
            hzz(i)=hzz(i)+(wxx*yyy(i)-wyy*xxx(i))
            
          endif

        enddo

c     reset positions to original reference frame
        
!!p        do i=1,natms
!!p          
!!p          xxx(i)=xxx(i)+cmx
!!p          yyy(i)=yyy(i)+cmy
!!p          zzz(i)=zzz(i)+cmz
!!p          
!!p        enddo
        
      endif

c     apply temperature scaling
      
      do i=1,natms
        
        hxx(i)=scale*hxx(i)
        hyy(i)=scale*hyy(i)
        hzz(i)=scale*hzz(i)
        
      enddo
      if(ngrp.gt.0) then
        call quatqnch(idnode,imcon,mxnode,natms,ngrp,
     x                hxx,hyy,hzz)

      endif
CVAM
CVAM      call VTEND(64, ierr)
CVAM
      return
      end

      subroutine shlqnch(idnode,mxnode,ntshl,temp,tstep,taus)
      
c*********************************************************************
c     
c     dl_poly subroutine for quenching the internal bond energies
c     in ions defined by shell model
c     
c     copyright - daresbury laboratory 1994
c     author w.smith july  1994
c     
c     wl
c     2005/12/02 14:22:19
c     1.3
c     Exp
c
c*********************************************************************
      
      use setup_module
      use config_module
      use core_shell_module
      
      implicit none

      integer idnode,mxnode,ntshl,ishl1,ishl2,i,j,k,m,ierr
      real*8 temp,pke,rmu,dvx,dvy,dvz,tmx,tmy,tmz,scl
      real*8 taus,tstep


CVAM
CVAM      call VTBEGIN(62, ierr)
CVAM

c     permitted core-shell internal kinetic energy 
      
      pke=boltz*temp*1.d-4

c     block indices

      ishl1 = (idnode*ntshl)/mxnode+1
      ishl2 = ((idnode+1)*ntshl)/mxnode

c     calculate core and shell velocities from total momentum
      
      m=0
      do k=ishl1,ishl2
        
        m=m+1
        
        i=listshl(m,2)
        j=listshl(m,3)

        rmu=(weight(i)*weight(j))/(weight(i)+weight(j))
        
        dvx=vxx(j)-vxx(i)
        dvy=vyy(j)-vyy(i)
        dvz=vzz(j)-vzz(i)

!        scl=sqrt(pke/(rmu*(dvx*dvx+dvy*dvy+dvz*dvz)))
! rescale the relative core-shell velocity
! with Langevin thermostat with temp=0 K
! and time scale taut
        scl=exp(-tstep/taus)


        tmx=weight(i)*vxx(i)+weight(j)*vxx(j)
        tmy=weight(i)*vyy(i)+weight(j)*vyy(j)
        tmz=weight(i)*vzz(i)+weight(j)*vzz(j)
        
        vxx(i)=tmx/(weight(i)+weight(j))-scl*rmu*dvx/weight(i)
        vxx(j)=tmx/(weight(i)+weight(j))+scl*rmu*dvx/weight(j)
        vyy(i)=tmy/(weight(i)+weight(j))-scl*rmu*dvy/weight(i)
        vyy(j)=tmy/(weight(i)+weight(j))+scl*rmu*dvy/weight(j)
        vzz(i)=tmz/(weight(i)+weight(j))-scl*rmu*dvz/weight(i)
        vzz(j)=tmz/(weight(i)+weight(j))+scl*rmu*dvz/weight(j)

      enddo

      if(mxnode.gt.1) call shlmerge(idnode,mxnode,ntshl)
CVAM
CVAM      call VTEND(62, ierr)
CVAM
      return
      end

      subroutine quatqnch(idnode,imcon,mxnode,natms,ngrp,
     x                    hxx,hyy,hzz)

c***********************************************************************
c     
c     dlpoly subroutine to convert atomic velocities to rigid body 
c     c.o.m. and angular velocity
c     
c     parallel replicated data version : block data
c     
c     copyright daresbury laboratory 1993.
c     author   - t.forester nov 1993.
c     author   - t.forester dec 1994 : block data.
c     
c     wl
c     2005/12/02 14:22:19
c     1.3
c     Exp
c     
c***********************************************************************

      use setup_module
      use config_module
      use rigid_body_module

      implicit none
      real *8, intent(inout)::hxx(natms), hyy(natms), hzz(natms)
      integer idnode,imcon,mxnode,natms,ngrp,fail,ierr,ig,jr,id
      integer igrp1,igrp2,i,j
      real*8 rot,wxx,wyy,wzz

      dimension rot(9)
      real*8, allocatable :: xxt(:),yyt(:),zzt(:)

      data fail/0/
CVAM
CVAM      call VTBEGIN(59, ierr)
CVAM
c     allocate work arrays 

      allocate (xxt(mxatms),yyt(mxatms),zzt(mxatms),stat=fail)
      if(fail.ne.0)call error(idnode,1780)

c     block indices for groups

      igrp1=(idnode*ngrp)/mxnode+1
      igrp2=((idnode+1)*ngrp)/mxnode

c     translate atomic velocites to com velocity & angular velocity

      jr=0
      do ig=igrp1,igrp2

        gvxx(ig)=0.d0
        gvyy(ig)=0.d0
        gvzz(ig)=0.d0
        omx(ig)=0.d0
        omy(ig)=0.d0
        omz(ig)=0.d0

        id=lstgtp(ig)

        do j=1,numgsit(id)

          jr =jr+1
          i =lstrgd(jr)

c     centre of mass momentum

          gvxx(ig)=gvxx(ig)+weight(i)*hxx(i)
          gvyy(ig)=gvyy(ig)+weight(i)*hyy(i)
          gvzz(ig)=gvzz(ig)+weight(i)*hzz(i)

c     distance to c.o.m of molecule

          xxt(jr)=xxx(i)-gcmx(ig)
          yyt(jr)=yyy(i)-gcmy(ig)
          zzt(jr)=zzz(i)-gcmz(ig)

        enddo

c     centre of mass velocity

        gvxx(ig)=gvxx(ig)/gmass(id)
        gvyy(ig)=gvyy(ig)/gmass(id)
        gvzz(ig)=gvzz(ig)/gmass(id)

      enddo

      call images(imcon,0,1,jr,cell,xxt,yyt,zzt)

      jr=0
      do ig=igrp1,igrp2

c     rotational matrix

        rot(1)=q0(ig)**2+q1(ig)**2-q2(ig)**2-q3(ig)**2
        rot(2)=2.d0*(q1(ig)*q2(ig)-q0(ig)*q3(ig))
        rot(3)=2.d0*(q1(ig)*q3(ig)+q0(ig)*q2(ig))
        rot(4)=2.d0*(q1(ig)*q2(ig)+q0(ig)*q3(ig))
        rot(5)=q0(ig)**2-q1(ig)**2+q2(ig)**2-q3(ig)**2
        rot(6)=2.d0*(q2(ig)*q3(ig)-q0(ig)*q1(ig))
        rot(7)=2.d0*(q1(ig)*q3(ig)-q0(ig)*q2(ig))
        rot(8)=2.d0*(q2(ig)*q3(ig)+q0(ig)*q1(ig))
        rot(9)=q0(ig)**2-q1(ig)**2-q2(ig)**2+q3(ig)**2
        
c     angular momentum accumulators

        wxx=0.d0
        wyy=0.d0
        wzz=0.d0

        id=lstgtp(ig)

        do j=1,numgsit(id)

          jr =jr+1
          i =lstrgd(jr)

          wxx=wxx+weight(i)*(yyt(jr)*hzz(i)-zzt(jr)*hyy(i))
          wyy=wyy+weight(i)*(zzt(jr)*hxx(i)-xxt(jr)*hzz(i))
          wzz=wzz+weight(i)*(xxt(jr)*hyy(i)-yyt(jr)*hxx(i))

        enddo

c     angular velocity in body fixed frame

        omx(ig)=(rot(1)*wxx+rot(4)*wyy+rot(7)*wzz)*rotinx(id,2)
        omy(ig)=(rot(2)*wxx+rot(5)*wyy+rot(8)*wzz)*rotiny(id,2)
        omz(ig)=(rot(3)*wxx+rot(6)*wyy+rot(9)*wzz)*rotinz(id,2)
        
        jr=jr-numgsit(id)
        do j=1,numgsit(id)
          
          jr=jr +1
          i=lstrgd(jr)
          
c     site velocity in body frame 

          wxx=omy(ig)*gzz(id,j)-omz(ig)*gyy(id,j)
          wyy=omz(ig)*gxx(id,j)-omx(ig)*gzz(id,j)
          wzz=omx(ig)*gyy(id,j)-omy(ig)*gxx(id,j)

c     new atomic velocites in lab frame

          hxx(i)=rot(1)*wxx+rot(2)*wyy+rot(3)*wzz+gvxx(ig)
          hyy(i)=rot(4)*wxx+rot(5)*wyy+rot(6)*wzz+gvyy(ig)
          hzz(i)=rot(7)*wxx+rot(8)*wyy+rot(9)*wzz+gvzz(ig)

        enddo

      enddo

      if(mxnode.gt.1) then
        call merge(idnode,mxnode,ngrp,mxbuff,gvxx,gvyy,gvzz,buffer)
        call merge(idnode,mxnode,ngrp,mxbuff,omx,omy,omz,buffer)
        call merge1(idnode,mxnode,natms,lstme,hxx,hyy,hzz,buffer)

      endif
c     deallocate work arrays

      deallocate (xxt,yyt,zzt,stat=fail)
CVAM
CVAM      call VTEND(59, ierr)
CVAM
      return
      end

