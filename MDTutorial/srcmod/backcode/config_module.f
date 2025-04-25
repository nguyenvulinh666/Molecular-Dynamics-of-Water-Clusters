      module config_module

c***********************************************************************
c     
c     dl_poly module for defining simulation configuration data
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2007/07/31 10:02:56
c     1.3
c     Exp
c     
c***********************************************************************

      use setup_module
      use error_module
      
      implicit none

      character*1 cfgname(80)
      character*1 sysname(80)
      real(8) cell(9),rcell(9),celprp(10)
      real(8) eta(9),stress(9),stresl(9),strcns(9),strbod(9)
      
      character*8, allocatable :: atmnam(:)
      real(8), allocatable :: xxx(:),yyy(:),zzz(:)
      real(8), allocatable :: vxx(:),vyy(:),vzz(:)
      real(8), allocatable :: fxx(:),fyy(:),fzz(:)
      real(8), allocatable :: flx(:),fly(:),flz(:)
      real(8), allocatable :: chge(:),weight(:),rmass(:)
      real(8), allocatable :: grxx(:,:),gryy(:,:),grzz(:,:) !! GLE
      integer, allocatable :: ltype(:),lstfrz(:)
      integer, allocatable :: neulst(:),lstneu(:)
      integer, allocatable :: lentry(:),list(:,:)
      integer, allocatable :: lstout(:),link(:)
      integer, allocatable :: lct(:),lst(:)

      real(8), allocatable :: buffer(:)

c PluMeD variables
      logical, save :: lplumed  
      character(len=20) :: plumedfile
      character(len=80) :: parse_file
c end PluMeD variables

      save atmnam,neulst,lstneu,cfgname,sysname
      save cell,xxx,yyy,zzz,vxx,vyy,vzz,fxx,fyy,fzz
      save buffer,weight,chge,ltype,lstfrz,flx,fly,flz
      save lentry,list,lstout,link,lct,lst,celprp,rmass
      save eta,stress,stresl,strcns,rcell
      
      contains
      
      subroutine alloc_config_arrays(idnode)

c***********************************************************************
c     
c     dl_poly subroutine for defining simulation configuration arrays
c     copyright - daresbury laboratory
c     author    - w. smith    sep 2003
c     
c     wl
c     2007/07/31 10:02:56
c     1.3
c     Exp
c     
c***********************************************************************
      
      integer, parameter :: nnn=33

      integer i,fail,idnode
      dimension fail(nnn)

      do i=1,nnn
        fail(i)=0
      enddo

      allocate (xxx(mxatms),stat=fail(1))
      allocate (yyy(mxatms),stat=fail(2))
      allocate (zzz(mxatms),stat=fail(3))
      allocate (vxx(mxatms),stat=fail(4))
      allocate (vyy(mxatms),stat=fail(5))
      allocate (vzz(mxatms),stat=fail(6))
      allocate (fxx(mxatms),stat=fail(7))
      allocate (fyy(mxatms),stat=fail(8))
      allocate (fzz(mxatms),stat=fail(9))
      allocate (weight(mxatms),stat=fail(11))
      allocate (chge(mxatms),stat=fail(12))
      allocate (ltype(mxatms),stat=fail(13))
      allocate (lstfrz(mxatms),stat=fail(14))
      allocate (flx(mxatms),stat=fail(15))
      allocate (fly(mxatms),stat=fail(16))
      allocate (flz(mxatms),stat=fail(17))
      allocate (atmnam(mxatms),stat=fail(18))
      allocate (neulst(mxneut),stat=fail(19))
      allocate (lstneu(mxatms),stat=fail(20))
      allocate (lstout(mxatms),stat=fail(21))
      allocate (lentry(msatms),stat=fail(22))
      allocate (list(msatms,mxlist),stat=fail(23))
      allocate (link(mxatms),stat=fail(24))
      allocate (lct(mxcell),stat=fail(25))
      allocate (lst(mxcell),stat=fail(26))
      allocate (rmass(mxatms),stat=fail(27))
      allocate (buffer(mxbuff),stat=fail(10))
      if (mxlgdof.gt.0) then                    !! GLE
        allocate (grxx(mxlgdof,mxatms), stat=fail(31))
        allocate (gryy(mxlgdof,mxatms), stat=fail(32))
        allocate (grzz(mxlgdof,mxatms), stat=fail(33))      
      endif



      do i=1,nnn
        if(fail(i).gt.0)then
          if(idnode.eq.0)write(nrite,'(10i5)')fail
          call error(idnode,1000)
        endif
      enddo

      end subroutine alloc_config_arrays

      end module config_module
