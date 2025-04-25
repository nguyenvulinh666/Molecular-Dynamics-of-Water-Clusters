      subroutine denfft
     x (fnum,kill,title,idnode,mxnode,kmax,ntime,ngap,tstep,key,wind,
     x  cfkt,work,fta)

c***********************************************************************
c
c    calculate fourier transform of density correlation function
c
c     author w smith 1994
c     copyright daresbury laboratory 1994
c
c***********************************************************************

      implicit real*8(a-h,o-z)

      parameter (isave1=99,isave2=55)
      parameter (pi=3.141592653589793d0,tpi=6.283185307179586d0)

      logical kill
      character*2 fnum
      character*80 title
      complex*16 cfkt(*),work(*),fta(*)
      dimension key(*),wind(*)
      data a0,a1,a2/0.42d0,0.50d0,0.08d0/

c     set control parameters

      ntime2=2*ntime
      ntime4=4*ntime
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

c     open files isave1 and isave2

      open(isave1,file='DENFKT.'//fnum,status='old',err=200)
      read(isave1,'(a80)')title
      read(isave1,*)
      read(isave1,'(2i10)')nfns,lfns

      if(idnode.eq.0)then

         open(isave2,file='DENSKW.'//fnum)
         write(isave2,'(a80)')title
         write(isave2,'(a)')'S(k,w) Data File'

      endif

c     initialise complex fast fourier transform routine

      ind=1
      isw=-1
      call fft(kill,idnode,ind,isw,ntime4,key,fta,work,fta)
      ind=0
c     set up window function (blackman function)

      arg=tpi/dble(ntime2)
      do i=1,ntime

         ccc=cos(arg*dble(i+ntime-1))
         wind(i)=a0-a1*ccc+a2*(2.d0*ccc**2-1.d0)

      enddo

c     read the correlation functions from disc

      do i=0,mxnode-1

         if(i.eq.idnode)then

            do k=1,nfns

               read(isave1,*)

               do m=1,lfns
                  
                  read(isave1,'(3e14.6)')
     x                 time,cfkt(m+(k-1)*ntime)
                  
               enddo
               
            enddo
            
         endif
         
      enddo

      tstep=time/(dble(lfns-1))

c     loop over correlation functions

      do k=1,iblk

c     apply window function

         do j=1,ntime4
            
            if(j.le.ntime)then
               
               fta(j)=wind(j)*cfkt(j+ntime*(k-1))

            else
               
               fta(j)=(0.d0,0.d0)
               
            endif
            
         enddo
         
         fta(1)=fta(1)/2.d0
         
c     apply complex fourier transform
         
         call fft(kill,idnode,ind,isw,ntime4,key,fta,work,fta)

c     store fourier coefficients
         
         m=1
         do j=1,ntime2,2

            cfkt(m+ntime*(k-1))=fta(j)
            m=m+1

         enddo

      enddo

c     save fourier coefficients

      write(isave2,'(2i10)')nfns,ntime/4

      omega=tpi/(dble(ntime2)*tstep)

      if(idnode.eq.0)then
           
         k=0
         mmin=0
         lmin=1
         do n=0,kmax
            do m=mmin,kmax
               do l=lmin,kmax
                  k=k+1

                  write(isave2,'(a8,3i6)')'k-vector',l,m,n
                  do ll=1,ntime/4
                     
                     write(isave2,'(1p,3e14.6)')
     x                    omega*dble(ll-1),cfkt(ll+(k-1)*ntime)
                     
                  enddo
                  
               enddo
               lmin=-kmax
            enddo
            mmin=-kmax
         enddo

      endif

c     close files

      close (isave1)
      close (isave2)

      write(3,'(a)')'Correlation file DENSKW.'//fnum//' created'
      write(3,'(a,i6)')'Number of correlation functions:',nfns

      return

  200 continue

      if(idnode.eq.0)
     x   write(3,'(a)')'Error - DENFKT.'//fnum//' file not found'
      kill=.true.
      return

      end
