
      module external_field_module
      
c***********************************************************************
c     
c     dl_poly module for defining external field potential arrays
c     copyright - daresbury laboratory
c     author    - w. smith    oct 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************
      
      use config_module
      use error_module
      use parse_module
      use setup_module
      
      implicit none
      
      real(8), allocatable :: prmfld(:)
      
      save prmfld
      
      contains
      
      subroutine alloc_fld_arrays(idnode)
      
      implicit none
      
      integer fail,idnode
      
      data fail/0/
      
      allocate (prmfld(mxfld),stat=fail)
      if(fail.ne.0)call error(idnode,1200)
      
      end subroutine alloc_fld_arrays
      
      subroutine define_external_field
     x  (safe,lunits,idnode,keyfld,engunit)
      
c***********************************************************************
c     
c     dl_poly subroutine to define external fields
c     
c     copyright - daresbury laboratory
c     author    - w. smith    aug 2003
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical safe,lunits
      character*8 keyword
      character*1 message(80)
      integer idnode,keyfld,nfld,i,k,idum
      real(8) engunit
      
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      
      call strip(record,lenrec)
      call lowcase(record,lenrec)
      call copystring(record,message,80)
      call getword(keyword,record,4,lenrec)
      
      if(keyword(1:4).eq.'elec') then
        keyfld =1 
      elseif(keyword(1:4).eq.'oshr') then
        keyfld=2
      elseif(keyword(1:4).eq.'shrx') then
        keyfld=3
      elseif(keyword(1:4).eq.'grav') then
        keyfld=4
      elseif(keyword(1:4).eq.'magn') then
        keyfld=5
      elseif(keyword(1:4).eq.'sphr') then
        keyfld=6
      elseif(keyword(1:4).eq.'zbnd') then
        keyfld=7
      else
        if(idnode.eq.0) write(nrite,*) message
        call error(idnode,454)
      endif
      
      do i = 1,mxfld
        prmfld(i)=0.d0
      enddo
      
      nfld=intstr(record,lenrec,idum)
      if(nfld.eq.0)nfld=5
      call getrec(safe,idnode,nfield)
      if(.not.safe)return
      do k=1,nfld
        
        prmfld(k)=dblstr(record,lenrec,idum)
        if(idum.gt.lenrec.and.k.lt.nfld)then
          
          call getrec(safe,idnode,nfield)
          if(.not.safe)return
          
        endif
        
      enddo
      
      if(idnode.eq.0) then
        
        write(nrite,"(/,/,1x,'external field key ',13x,a4,
     x    /,/,30x,'external field parameters')") keyword(1:4)
        write(nrite,"(2(/,1x,1p,5e15.5))") prmfld
        
      endif      
      
c     convert to internal units
      
      if(keyfld.eq.1.or.keyfld.eq.4.or.keyfld.eq.5) then
        
        if(.not.lunits)call error(idnode,6)
        
        do i = 1,3
          prmfld(i) = prmfld(i)*engunit
        enddo
        
      elseif(keyfld.eq.2.or.keyfld.eq.6.or.keyfld.eq.7) then
        
        prmfld(1) = prmfld(1)*engunit
        
      endif
      
      return
      end subroutine define_external_field
      
      subroutine extnfld
     x  (idnode,imcon,keyfld,mxnode,natms,engfld,virfld)
      
c***********************************************************************
c     
c     dl_poly routine for application of an external field
c     
c     replicated data version / block data
c     
c     copyright daresbury laboratory 1993
c     author -    t.forester october 1993
c     amended-    t.forester dec 1994
c     
c     wl
c     2008/12/23 10:29:12
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      integer idnode,imcon,keyfld,mxnode,natms,iatm1,iatm2,i
      real(8) engfld,virfld,rz,rrr,gamma,zdif
      
CVAM
CVAM      call VTBEGIN(27, ierr)
CVAM
      
c     energy and virial accumulators 
      
      engfld = 0.d0
      virfld = 0.d0
      
c     block indices
      
      iatm1 = (idnode*natms)/mxnode + 1
      iatm2 = ((idnode+1)*natms)/mxnode
      
      if(keyfld.eq.1) then
        
c     electric field: prmfld(1-3) are field components
        
        do i = iatm1,iatm2
          
          fxx(i) = fxx(i)+ chge(i)*prmfld(1)
          fyy(i) = fyy(i)+ chge(i)*prmfld(2)
          fzz(i) = fzz(i)+ chge(i)*prmfld(3)
          
        enddo
        
      elseif(keyfld.eq.2) then
        
c     oscillating shear: orthorhombic box:  Fx = a*cos(b.2.pi.z/L)
        
        rz = 2.d0*pi/cell(9)
        
        do i = iatm1,iatm2
          
          fxx(i) = fxx(i) + prmfld(1)*cos(prmfld(2)*zzz(i)*rz)
          
        enddo
        
      elseif(keyfld.eq.3.and.imcon.eq.6) then
        
c     continuous shear of walls : 2D periodic box (imcon=6)
c     shear rate = prmfld(1) angstrom per ps for atoms at
c     abs(z) > prmfld(2)
        
        do i=iatm1,iatm2
          
          if(abs(zzz(i)).gt.prmfld(2)) then
            
            vxx(i) = 0.5d0*sign(prmfld(1),zzz(i))
            
          endif
          
        enddo
        
      elseif(keyfld.eq.4) then
        
c     gravitational field: field components given by prmfld(1-3)
        
        do i =iatm1,iatm2
          
          fxx(i) = fxx(i) + prmfld(1)*weight(i)
          fyy(i) = fyy(i) + prmfld(2)*weight(i)
          fzz(i) = fzz(i) + prmfld(3)*weight(i)
          
        enddo
        
      elseif(keyfld.eq.5) then
        
c     magnetic field: field components given by prmfld(1-3)
        
        do i = iatm1,iatm2
          
          fxx(i)=fxx(i)+(vyy(i)*prmfld(3)-vzz(i)*prmfld(2))
     x      *chge(i)
          fyy(i)=fyy(i)+(vzz(i)*prmfld(1)-vxx(i)*prmfld(3))
     x      *chge(i)
          fzz(i)=fzz(i)+(vxx(i)*prmfld(2)-vyy(i)*prmfld(1))
     x      *chge(i)
          
        enddo
        
      elseif(keyfld.eq.6) then
        
c     containing sphere : r^(-n) potential
        
        do i = iatm1,iatm2
          
          rrr = sqrt(xxx(i)**2+yyy(i)**2+zzz(i)**2)
          if(rrr.gt.prmfld(4)) then
            rrr = prmfld(2) - rrr
            if(rrr.lt.0.d0) rrr = 0.1d0
            
            gamma  = prmfld(1)*rrr**(-prmfld(3))
            engfld = engfld + gamma
            
            gamma = -prmfld(3)*gamma/((prmfld(2)-rrr)*rrr)
            
            fxx(i)=fxx(i)+ gamma*xxx(i)
            fyy(i)=fyy(i)+ gamma*yyy(i)
            fzz(i)=fzz(i)+ gamma*zzz(i)
            
          endif
          
        enddo
        
      elseif(keyfld.eq.7) then
        
c     repulsive wall (harmonic) starting at z0
        
        do i = iatm1,iatm2
          
          if(prmfld(3)*zzz(i).gt.prmfld(3)*prmfld(2)) then
            
            zdif = zzz(i) - prmfld(2)
            gamma = -prmfld(1)*zdif
            
            fzz(i) = fzz(i) + gamma
            engfld = engfld - gamma*zdif/2.
            
          endif
          
        enddo
        
      else
        
c     unidentified field potential error exit
        
        call error(idnode,454)
        
      endif
      
c     global sum of external field potential and virial
      
      if(mxnode.gt.1)then
        
        buffer(1)=engfld
        buffer(2)=virfld
        call gdsum(buffer(1),2,buffer(3))
        engfld=buffer(1)
        virfld=buffer(2)
        
      endif
CVAM
CVAM      call VTEND(27, ierr)
CVAM
      return
      end subroutine extnfld
      
      end module external_field_module
