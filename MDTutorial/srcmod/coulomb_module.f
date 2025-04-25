      module coulomb_module

c***********************************************************************
c     
c     dl_poly module for defining coulomb terms
c     copyright - daresbury laboratory
c     
c     author    - w. smith    sep 2003
c     adapted for solvation, free energy and excitation
c               - p.-a. cazade oct 2007
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************

      use config_module
      use ewald_module
      use pair_module
      use property_module
      use setup_module
      use solvation_module

      implicit none

      contains
      
      subroutine coul0
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     1/r potential, no truncation or damping
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester february 1993
c     stress tensor - t.forester may 1994
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea,rsq
      real(8) chgprd,rrr,coul,fcoul,fi,fx,fy,fz

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(84, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0     
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)/epsq*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik

c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif

          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10)then
            
c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              
c     coulomb potential and force
              
              coul=chgprd/rrr
              fcoul=coul/rsq
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-coul
                
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate forces
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
              
                strs1=strs1+xdf(m)*fx
                strs2=strs2+xdf(m)*fy
                strs3=strs3+xdf(m)*fz
                strs5=strs5+ydf(m)*fy
                strs6=strs6+ydf(m)*fz
                strs9=strs9+zdf(m)*fz
                
              endif
              
            endif
            
          endif
          
        enddo
        
c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)
        
c     complete stress tensor
        
        stress(1)=stress(1)+strs1
        stress(2)=stress(2)+strs2
        stress(3)=stress(3)+strs3
        stress(4)=stress(4)+strs2
        stress(5)=stress(5)+strs5
        stress(6)=stress(6)+strs6
        stress(7)=stress(7)+strs3
        stress(8)=stress(8)+strs6
        stress(9)=stress(9)+strs9

      endif
CVAM
CVAM      call VTEND(84, ierr)
CVAM
      return
      end subroutine coul0

      subroutine coul1
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a standard coulomb potential truncated at rcut
c     and shifted to zero at rcut.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1992
c     author    - w. smith december 1992.
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     stress tensor t.forester may 1994
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea,rsq
      real(8) fi,chgprd,omega,egamma,fx,fy,fz,rrr

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(85, ierr)
CVAM
      lskip=(lfree.or.lghost)

c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      chgea=chge(iatm)*r4pie0/epsq

      if(abs(chgea).gt.1.d-10)then

c     start of primary loop for forces evaluation
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10) then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)

c     calculate potential energy and virial

              omega=chgprd*(rcut-rrr)/(rrr*rcut)
              egamma=chgprd/(rrr*rsq)
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+omega
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-omega
                  omega=lambda1*omega
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+omega
                  omega=lambda2*omega
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate potential energy and virial
              
                engcpe=engcpe+omega
                vircpe=vircpe-egamma*rsq
                
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+omega
              
c     calculate forces
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
              
                strs1=strs1+xdf(m)*fx
                strs2=strs2+xdf(m)*fy
                strs3=strs3+xdf(m)*fz
                strs5=strs5+ydf(m)*fy
                strs6=strs6+ydf(m)*fz
                strs9=strs9+zdf(m)*fz
                
              endif

            endif
            
          endif

        enddo

c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)

c     complete stress tensor
        
        stress(1)=stress(1)+strs1
        stress(2)=stress(2)+strs2
        stress(3)=stress(3)+strs3
        stress(4)=stress(4)+strs2
        stress(5)=stress(5)+strs5
        stress(6)=stress(6)+strs6
        stress(7)=stress(7)+strs3
        stress(8)=stress(8)+strs6
        stress(9)=stress(9)+strs9

      endif
CVAM
CVAM      call VTEND(85, ierr)
CVAM
      return
      end subroutine coul1

      subroutine coul2
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a distance dependant dielectric `constant'.
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1993
c     author    - t. forester    april 1993
c     stress tensor added - t.forester may 1994
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,kkk
      real(8) engcpe,vircpe,rcut,epsq
      real(8) fi,rcsq,strs1,strs2,strs3,strs5,strs6,strs9,chgea
      real(8) chgprd,rsq,rrsq,coul,egamma,fx,fy,fz
      
      dimension fi(3)
      
CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(86, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      
c     initialise stress tensor accumulators
      
      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0
      
c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      chgea=chge(iatm)/epsq*r4pie0
      if(abs(chgea).gt.1.d-10)then
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
c     start of primary loop for forces evaluation
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10)then
            
c     calculate interatomic distance
            
            rsq=rsqdf(m)
            
c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
c     calculate potential energy and Virial
              
              rrsq=1.d0/rsq
              coul=chgprd*rrsq
              egamma=2.d0*coul*rrsq
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  egamma=lambda1*egamma
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  egamma=lambda2*egamma
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate potential energy and Virial
                
                engcpe=engcpe+coul
                vircpe=vircpe-2.d0*coul
                
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
              
c     calculate forces
                
                fx=egamma*xdf(m)
                fy=egamma*ydf(m)
                fz=egamma*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz

                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
              
c     calculate stress tensor
                
                strs1=strs1+xdf(m)*fx
                strs2=strs2+xdf(m)*fy
                strs3=strs3+xdf(m)*fz
                strs5=strs5+ydf(m)*fy
                strs6=strs6+ydf(m)*fz
                strs9=strs9+zdf(m)*fz
                
              endif
              
            endif
            
          endif
          
        enddo
        
c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)
        
c     complete stress tensor
        
        stress(1)=stress(1)+strs1
        stress(2)=stress(2)+strs2
        stress(3)=stress(3)+strs3
        stress(4)=stress(4)+strs2
        stress(5)=stress(5)+strs5
        stress(6)=stress(6)+strs6
        stress(7)=stress(7)+strs3
        stress(8)=stress(8)+strs6
        stress(9)=stress(9)+strs9
        
      endif
CVAM
CVAM      call VTEND(86, ierr)
CVAM
      return
      end subroutine coul2

      subroutine coul3
     x  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic force.
c     reaction field  potential
c     Ref: M Neumann, J Chem Phys, 82, 5633, (1985)
c     adapted for fennell-gezelter coulombic model
c     by w.smith june 2007
c     Ref: CJ Fennell and JD Gezelter, J Chem Phys, 
c     124, 234104, (2006)
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    - t. forester february 1995
c     stress tensor - t.forester   feb 1995
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,l,kkk
      real(8) engcpe,vircpe,rcut,epsq,vcon,fcon,rdr,ppp
      real(8) fi,rcsq,b0,rfld0,rfld1,rfld2,strs1,strs2,strs3
      real(8) strs5,strs6,strs9,chgea,chgprd,rsq,coul,omega
      real(8) fx,fy,fz,fcoul,rrr,vk0,vk1,vk2,gk0,gk1,gk2,t1,t2
      dimension fi(3)
CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(87, ierr)
CVAM
      lskip=(lfree.or.lghost)
        
c     reaction field terms
      
      b0=2.d0*(epsq-1.d0)/(2.d0*epsq+1.d0)
      rfld0=b0/rcut**3
      rfld1=(1.d0+b0*0.5d0)/rcut
      rfld2=rfld0*0.5d0
      
c     screened coulomb terms
        
      vcon=erc(mxegrd-4)+rfld2*rcut**2-rfld1
      fcon=rcut*fer(mxegrd-4)-rfld0*rcut
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2
      rdr=dble(mxegrd-4)/rcut

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     start of primary loop for forces evaluation
      
      chgea=chge(iatm)*r4pie0
      
      if(abs(chgea).gt.1.d-10)then
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik

c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              l=int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              
c     calculate potential energy using 3-point interpolation
              
              vk0=erc(l)
              vk1=erc(l+1)
              vk2=erc(l+2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              omega=t1+(t2-t1)*ppp*0.5d0-vcon+fcon*(rrr-rcut)
              coul=chgprd*(omega+rfld2*rsq-rfld1)

c     calculate forces using 3-point interpolation
              
              gk0=fer(l)
              gk1=fer(l+1)
              gk2=fer(l+2)
              t1=gk0+(gk1-gk0)*ppp
              t2=gk1+(gk2-gk1)*(ppp-1.0d0)
              fcoul=chgprd*((t1+(t2-t1)*ppp*0.5d0)-fcon/rrr-rfld0)
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate coulombic energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-fcoul*rsq
              
c     calculate solvation energy
              
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate coulombic force
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
              
                strs1=strs1+xdf(m)*fx
                strs2=strs2+xdf(m)*fy
                strs3=strs3+xdf(m)*fz
                strs5=strs5+ydf(m)*fy
                strs6=strs6+ydf(m)*fz
                strs9=strs9+zdf(m)*fz
                
              endif
              
            endif

          endif
          
        enddo
        
c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)

c     complete stress tensor
        
        stress(1)=stress(1)+strs1
        stress(2)=stress(2)+strs2
        stress(3)=stress(3)+strs3
        stress(4)=stress(4)+strs2
        stress(5)=stress(5)+strs5
        stress(6)=stress(6)+strs6
        stress(7)=stress(7)+strs3
        stress(8)=stress(8)+strs6
        stress(9)=stress(9)+strs9

      endif
CVAM
CVAM      call VTEND(87, ierr)
CVAM
      return
      end subroutine coul3

      subroutine coul4
     X  (lsolva,lfree,lghost,iatm,ik,engcpe,vircpe,rcut,epsq)

c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic forces 
c     assuming a force shifted coulomb potential.
c     adapted for fennell-gezelter coulombic model
c     by w.smith may 2007
c     Ref: CJ Fennell and JD Gezelter, J Chem Phys, 
c     124, 234104, (2006)
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory 1995
c     author    -  t.forester october  1995
c     adapted   - p.-a. cazade oct 2007: solvation, free energy etc
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none
      
      logical lsolva,lfree,lghost,lselect,lskip
      integer iatm,ik,m,jatm,l,kkk
      real(8) engcpe,vircpe,rcut,epsq,vcon,fcon,rdr,ppp
      real(8) fi,rcsq,strs1,strs2,strs3,strs5,strs6,coul
      real(8) strs9,chgea,chgprd,rsq,rrr,omega,fcoul,fx,fy,fz
      real(8) vk0,vk1,vk2,gk0,gk1,gk2,t1,t2

      dimension fi(3)

CDIR$ CACHE_ALIGN fi
CVAM
CVAM      call VTBEGIN(88, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     screened coulomb terms
        
      vcon=erc(mxegrd-4)
      fcon=rcut*fer(mxegrd-4)
      rdr=dble(mxegrd-4)/rcut
      
c     set cutoff condition for pair forces
      
      rcsq=rcut**2

c     initialise stress tensor accumulators

      strs1=0.d0
      strs2=0.d0
      strs3=0.d0
      strs5=0.d0
      strs6=0.d0
      strs9=0.d0

c     initialise potential energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
      chgea=chge(iatm)*r4pie0/epsq

      if(abs(chgea).gt.1.d-10)then

c     start of primary loop for forces evaluation
        
        fi(1)=fxx(iatm)
        fi(2)=fyy(iatm)
        fi(3)=fzz(iatm)
        
        do m=1,ik
          
c     atomic index and charge product
          
          jatm=ilist(m)
          
          if(lskip)then
            if(atm_fre(iatm)*atm_fre(jatm).eq.2)cycle
          endif
          
          chgprd=chgea*chge(jatm)
          if(abs(chgprd).gt.1.d-10)then

c     calculate interatomic distance
            
            rsq=rsqdf(m)

c     apply truncation of potential
            
            if(rcsq.gt.rsq)then
              
              rrr=sqrt(rsq)
              l=int(rrr*rdr)
              ppp=rrr*rdr-dble(l)
              
c     calculate potential energy using 3-point interpolation
              
              vk0=erc(l)
              vk1=erc(l+1)
              vk2=erc(l+2)
              t1=vk0+(vk1-vk0)*ppp
              t2=vk1+(vk2-vk1)*(ppp-1.0d0)
              omega=t1+(t2-t1)*ppp*0.5d0
              coul=chgprd*(omega-vcon+fcon*(rrr-rcut))
              
c     calculate forces using 3-point interpolation
              
              gk0=fer(l)
              gk1=fer(l+1)
              gk2=fer(l+2)
              t1=gk0+(gk1-gk0)*ppp
              t2=gk1+(gk2-gk1)*(ppp-1.0d0)
              fcoul=chgprd*((t1+(t2-t1)*ppp*0.5d0)-fcon/rrr)
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(iatm),atmolt(jatm))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(iatm).ne.1).and.(atm_fre(jatm).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(iatm)+atm_fre(jatm).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc(kkk)=cou_exc(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(iatm).eq.1).or.(atm_fre(jatm).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(iatm).eq.2).or.(atm_fre(jatm).eq.2))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then

c     calculate the coulombic energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-fcoul*rsq
                
c     calculate solvation energy
                
                if(lsolva)cou_sol(kkk)=cou_sol(kkk)+coul
                
c     calculate coulombic forces
                
                fx=fcoul*xdf(m)
                fy=fcoul*ydf(m)
                fz=fcoul*zdf(m)
                
                fi(1)=fi(1)+fx
                fi(2)=fi(2)+fy
                fi(3)=fi(3)+fz
                
                fxx(jatm)=fxx(jatm)-fx
                fyy(jatm)=fyy(jatm)-fy
                fzz(jatm)=fzz(jatm)-fz
                
c     calculate stress tensor
                
                strs1=strs1+xdf(m)*fx
                strs2=strs2+xdf(m)*fy
                strs3=strs3+xdf(m)*fz
                strs5=strs5+ydf(m)*fy
                strs6=strs6+ydf(m)*fz
                strs9=strs9+zdf(m)*fz
                
              endif
              
            endif
            
          endif

        enddo

c     load temps back to fxx(iatm) etc
        
        fxx(iatm)=fi(1)
        fyy(iatm)=fi(2)
        fzz(iatm)=fi(3)

c     complete stress tensor
        
        stress(1)=stress(1)+strs1
        stress(2)=stress(2)+strs2
        stress(3)=stress(3)+strs3
        stress(4)=stress(4)+strs2
        stress(5)=stress(5)+strs5
        stress(6)=stress(6)+strs6
        stress(7)=stress(7)+strs3
        stress(8)=stress(8)+strs6
        stress(9)=stress(9)+strs9

      endif
CVAM
CVAM      call VTEND(88, ierr)
CVAM
      return
      end subroutine coul4
      
      subroutine coul_nsq
     x  (lsolva,lfree,lghost,idnode,mxnode,natms,imcon,epsq,rcut,
     x  engcpe,vircpe)
      
c***********************************************************************
c     
c     dl_poly subroutine for calculating coulombic potential and forces
c     for the all-pairs algorithm beyond the range of the normal cutoff
c     i.e. the 'tertiary' forces.  frozen atom option included
c     
c     to be used with multiple_nsq
c     
c     parallel replicated data version
c     
c     copyright - daresbury laboratory
c     author    - w.smith  august 2008
c     
c     wl
c     2008/12/23 10:29:11
c     1.4
c     Exp
c     
c***********************************************************************
      
      implicit none

      logical lsolva,lfree,lghost,lselect,lskip
      integer natms,idnode,mxnode,imcon,ibig,i,last,mpm2
      integer npm2,m,ii,j,idum,kkk
      real(8) engcpe,epsq,rcut,vircpe,rsq,rrr,chgprd,fcoul,coul,rct2

CVAM
CVAM      call VTBEGIN(x, ierr)
CVAM
      lskip=(lfree.or.lghost)
      
c     zero energy and virial
      
      engcpe=0.d0
      vircpe=0.d0
      
c     zero force arrays
      
      do i=1,natms
        
        flx(i)=0.d0
        fly(i)=0.d0
        flz(i)=0.d0
        
      enddo
      
c     zero stress tensor
      
      do i=1,9
        stresl(i)=0.d0
      enddo
      
c     zero solvation and  excitation accumulators
      
      if(lsolva)then
        
        cou_sol_lng(:)=0.d0
        
        if(lghost)then
          
          cou_exc_lng(:)=0.d0
          
        endif
        
      endif
      
c     set control variables
      
      last=natms
      mpm2=natms/2
      npm2=(natms-1)/2
      
c     set cutoff radius
      
      rct2=rcut**2
      
c     outer loop over atoms
      
      do m=1,mpm2
        
        if(m.gt.npm2)last=mpm2
        
c     inner loop over atoms
        
        ii=0
        do i=idnode+1,last,mxnode
          
c     calculate atom indices
          
          j=i+m
          if(j.gt.natms)j=j-natms
          
          if(lskip)then
            if(atm_fre(i)*atm_fre(j).eq.2)cycle
          endif
          
          if(lstfrz(i)*lstfrz(j).eq.0)then
            
            ii=ii+1
            
c     calculate interatomic displacements
            
            xdf(ii)=xxx(i)-xxx(j)
            ydf(ii)=yyy(i)-yyy(j)
            zdf(ii)=zzz(i)-zzz(j)
            
          endif
          
        enddo
        
c     apply minimum image convention
        
        call images(imcon,0,1,ii,cell,xdf,ydf,zdf)

c     calculate coulomb terms
        
        ii=0
        
        do i=idnode+1,last,mxnode
          
c     calculate atom indices
          
          j=i+m
          if(j.gt.natms)j=j-natms
          
          if(lskip)then
            if(atm_fre(i)*atm_fre(j).eq.2)cycle
          endif
          
          ii=ii+1
          if(lstfrz(i).eq.0.or.lstfrz(j).eq.0)then
            
c     reject frozen atoms and calculate interatomic distance
            
            rsq=xdf(ii)**2+ydf(ii)**2+zdf(ii)**2
            
c     running check of neighbour list array capacity
            
            if(rsq.ge.rct2)then
              
              chgprd=chge(i)*chge(j)*r4pie0/epsq
              rrr=sqrt(rsq)
              
c     calculate potential energy and force
              
              coul=chgprd/rrr
              fcoul=coul/rsq
              
c     set selection control
              
              lselect=.true.
              
c     set double index
              
              if(lsolva)kkk=loc2(atmolt(i),atmolt(j))
              
              if(lghost)then
                
c     selected excitation option
                
                if((atm_fre(i).ne.1).and.(atm_fre(j).ne.1))then
                  
c     reset selection control
                  
                  lselect=(atm_fre(i)+atm_fre(j).eq.0)
                  
c     calculate solvation energy
                  
                  if(lsolva)cou_exc_lng(kkk)=cou_exc_lng(kkk)+coul
                  
                endif
                
              elseif(lfree)then
                
c     selected free energy option
                
                if((atm_fre(i).eq.1).or.(atm_fre(j).eq.1))then
                  
c     set hamiltonian mixing parameter
                  
                  cou_fre=cou_fre-coul
                  coul=lambda1*coul
                  fcoul=lambda1*fcoul
                  
                elseif((atm_fre(i).eq.2).or.(atm_fre(j).eq.2))then
                  
c     set hamiltonian mixing parameter

                  cou_fre=cou_fre+coul
                  coul=lambda2*coul
                  fcoul=lambda2*fcoul
                  
                endif
                
              endif
              
              if(lselect)then
                
c     calculate potential energy and virial
              
                engcpe=engcpe+coul
                vircpe=vircpe-coul

c     calculate solvation energy
              
                if(lsolva)cou_sol_lng(kkk)=cou_sol_lng(kkk)+coul
                
c     calculate forces
                
                flx(i)=flx(i)+fcoul*xdf(ii)
                fly(i)=fly(i)+fcoul*ydf(ii)
                flz(i)=flz(i)+fcoul*zdf(ii)             
                
                flx(j)=flx(j)-fcoul*xdf(ii)
                fly(j)=fly(j)-fcoul*ydf(ii)
                flz(j)=flz(j)-fcoul*zdf(ii)  
                
c     stress tensor
                
                stresl(1)=stresl(1)+xdf(ii)*fcoul*xdf(ii)
                stresl(2)=stresl(2)+xdf(ii)*fcoul*ydf(ii)
                stresl(3)=stresl(3)+xdf(ii)*fcoul*zdf(ii)
                stresl(5)=stresl(5)+ydf(ii)*fcoul*ydf(ii)
                stresl(6)=stresl(6)+ydf(ii)*fcoul*zdf(ii)
                stresl(9)=stresl(9)+zdf(ii)*fcoul*zdf(ii)
                
              endif
              
            endif
            
          endif
          
        enddo
          
      enddo

c     complete stress tensor

      stresl(4)=stresl(2)
      stresl(7)=stresl(3)
      stresl(8)=stresl(6)
      
CVAM
CVAM      call VTEND(x,ierr)
CVAM
      return
      end subroutine coul_nsq
      
      end module coulomb_module
