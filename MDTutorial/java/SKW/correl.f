      subroutine correl
     x (fnum,kill,title,idnode,mxnode,ncon,kmax,ntime,ngap,
     x  tstep,ind,num,cfkt,ckr0,ckr)

c***********************************************************************
c     
c     calculate density correlation function
c     
c     author w smith 1994
c     copyright daresbury laboratory 1994
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      parameter (iwork=88,isave1=99)
      logical kill
      character*2 fnum
      character*80 title
      complex*16 ckr,ckr0,cfkt
      dimension cfkt(*),ckr0(*),ckr(*),ind(*),num(*)

c     open files isave1 and iwork

      open(iwork,file='SPCDEN.'//fnum,form='unformatted',
     x  status='old',err=300)


      if(idnode.eq.0)then

         open(isave1,file='DENFKT.'//fnum)
         write(isave1,'(a80)')title
         write(isave1,'(a)')'F(k,t) Data File'

      endif

c     set control parameters
      norg=ntime/ngap
      klim=((2*kmax+1)**3-1)/2
      nblk=mod(klim,mxnode)
      if(idnode.lt.nblk)then
         iblk=(klim/mxnode+1)
         ibgn=idnode*iblk
      else
         iblk=klim/mxnode
         ibgn=idnode*iblk+nblk
      endif

      lor=0
      mor=0
      
c     initialise arrays

      do l=1,iblk*ntime
        cfkt(l)=(0.d0,0.d0)
      enddo
      do l=1,ntime
        num(l)=0
      enddo

c     start of loop over time steps

      do n=0,ncon

        read(iwork,end=200)time,(ckr(j),j=1,klim)
        if(mod(n,ngap).le.0)then

          lor=min(lor+1,norg)
          mor=mod(mor,norg)+1
          ind(mor)=1

          do k=1,iblk
            ckr0(mor+ntime*(k-1))=conjg(ckr(k+ibgn))
          enddo

        endif

        do l=1,lor

          m=ind(l)
          ind(l)=m+1
          num(m)=num(m)+1
          do k=1,iblk
            cfkt(m+ntime*(k-1))=cfkt(m+ntime*(k-1))+ckr0(l+ntime*(k-1))
     x                          *ckr(k+ibgn)
          enddo

        enddo

      enddo
  200 continue

      last=min(n-1,ntime)
      tstep=time/dble(n-1)
      write(isave1,'(2i10)')iblk,last

c     normalise correlation functions

      do k=1,iblk
        
        rnorm=dble(num(1))/dreal(cfkt(1+(k-1)*ntime))

        do l=1,last
            cfkt(l+(k-1)*ntime)=rnorm*cfkt(l+(k-1)*ntime)/dble(num(l))
        enddo

      enddo

c     store correlation functions in disc file


      if(idnode.eq.0)then
           
         k=0
         mmin=0
         lmin=1
         do n=0,kmax
            do m=mmin,kmax
               do l=lmin,kmax
                  k=k+1

                  write(isave1,'(a8,3i6)')'k-vector',l,m,n
                  do ll=1,last
                     
                     write(isave1,'(1p,3e14.6)')
     x                    tstep*dble(ll-1),cfkt(ll+(k-1)*ntime)
                     
                  enddo
                  
               enddo
               lmin=-kmax
            enddo
            mmin=-kmax
         enddo

      endif

c     close files

      close (isave1)
      close (iwork)

      write(3,'(a)')'Correlation file DENFKT.'//fnum//' created'
      write(3,'(a,i6)')'Number of correlation functions:',iblk

      return
  300 continue

      if(idnode.eq.0)
     x   write(3,'(a)')'Error - SPCDEN.'//fnum//' file not found'
      kill=.true.
      return

      end
