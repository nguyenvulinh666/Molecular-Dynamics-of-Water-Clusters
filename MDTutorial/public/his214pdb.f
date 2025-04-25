c This program convert the DLPOLY 2.14 HISTORY file to PDB format
c The output is primarily designed to be used in VMD
c
c Author:  Ambarish Kulkarni
c Address: Georgia Institue of Technology, USA
c Email:   gtg156s@mail.gatech.edu       
c
c
	program his2pdb
	implicit none
	integer ierror,natoms,levcfg,imcon,nsteps
	integer i,stepsize,k,writestep,choice,ilen
	real*8 x,y,z,value,q,w,xa,xb,xc,ya,yb,yc,za,zb,zc
      real*8,allocatable ::  r(:,:),v(:,:),f(:,:)
	character (len=100)  in_file,traj,outfile
	character (len=10)  time,str
      character*4 charp
c
1210  format(g2.0)	
900   format (i1,i4,6f12.4)
820   format(a4,2x,i5,x,a2,2x,a1,a3,x,a1,i4,a1,3x,f8.3,f8.3,
     x          f8.3,f6.2,f6.2,6x,a4,a2,a2)
c
c read input file name
      write(*,*) 'Convert DLPOLY HISTORY file to PDB format : '
	write(*,*) 'The input file name:'
	read(*,'(a)') in_file
c
c
c the program can generate a single file with all timesteps or
c a file per timestep
      write(*,*) '(1) Single PDB file/(2) Multiple PDB files :'
      read (*,*) choice
c
c
c read output file name
	write(*,*) 'The output file name (do not add extension) :'
	read(*,'(a)') traj

c  single file name
      if(choice.eq.1) then     
	  outfile = trim(traj)//'.pdb'
	  open (unit=2,file=outfile,status='replace',iostat=ierror)
      endif
c
c  number of simulation steps in the HISTORY file
	write(*,*) 'Enter the number of steps:'
	read(*,*) nsteps
	open (unit=13, file=in_file, status='old')

c
c start reading HISTORY file and writing PDB file
	do k=1,nsteps

c for multiple files
      if (choice.eq.2) then
        call writeint(k,str,ilen)
        outfile = trim(traj)//str(1:ilen)//'.pdb'
        open =
(unit=2,file=outfile,status='replace',iostat=ierror)
      endif
c
c only the first step requires reading file header and other parameters
   	    if (k==1) then
c
c read the file header
	         read (13,*,end=999,err=999)
c
c read trajectory key,periodic boundary key, # atoms
	         read (13,*,end=999,err=999) levcfg,imcon,natoms 
	    endif
c
c read timestep, current timestep, #atoms
c      trajectory key, periodic boundary key, integration timestep      
          read (13,*,end=999,err=999) time,stepsize,natoms,
     x         levcfg,imcon,q
c
c read all components of the cell vector
	    read(13,*) xa,ya,za
	    read(13,*) xb,yb,zb
	    read(13,*) xc,yc,zc
c
c allocate the coordinate array
		if(k == 1) then
	      allocate(r(3,natoms))
	      allocate(v(3,natoms))
	      allocate(f(3,natoms))
	    endif

c
c loop over atoms to read data
		do i = 1,natoms
c
c charp is used to read the atom name
	     read(13,1210) charp
c
c based on trajectory key read in coordinates, velocity and forces
           read(13,*,end=999,err=999) r(1,i),r(2,i),r(3,i)
		 if(levcfg.gt.0)read(13,*,end=999,err=999)v(1,i),v(2,i),v(3,i)
		 if(levcfg.gt.1)read(13,*,end=999,err=999)f(1,i),f(2,i),f(3,i)

c
c  I have used two atom types with names AA+ and BB- in my HISTORY file
c  These can be changed to read in any other atom names. 
	     if(charp.eq.'AA') then
c
c value can be used to set the occupancy which can be used to =
differentiate
c between various atom types in VMD          
		  value=1.0
c
c another visualization method in VMD is to use some preset atom names =
stored
c in its database. Setting charp as any one of those is another method
c of differentiating the atom types
	      charp='N'
	      endif
	      if(charp.eq.'BB') then
              value=2.0
			charp='O' 
              else
              value=3.0
              endif
c
c writestep is used to drop out the time steps which are not required
            writestep=1
c
c begin writing the PDB file
		  if ((k.eq.1).or.(mod(k,writestep).eq.0)) then
	      write(2,820) 'ATOM',i,charp,' ','   ',' ',
     x         1,'A',r(1,i),r(2,i),r(3,i),
     x         value, 2.0,' ','  ','  '
	           endif
c
c end atom loop
	    enddo

c
c print current step on console and write END at each step in PDB file 
		if ((k.eq.1).or.(mod(k,writestep).eq.0)) then
            print*, "Writing step",k," to PDB file"
		  write(2,'(a3)')'END'
          endif
c
c end reading HISTORY and writing PDB file
	enddo
	close(2)
	close(13)

999   continue
      write(*,*) 'File conversion finished!'
	write(*,*) 'Press enter to exit'
	read(*,*)
	end program his2pdb


c write a positive integer into the beginning of a string variable
c return inum = length of string

      subroutine writeint(ivalue,str,inum)
      implicit none
      integer ivalue
      character*(*) str
      integer inum

      if (ivalue.lt.10) then
        write (str,'(i1)') ivalue
          str = '0000'//str(1:1)
        inum = 5
      else if (ivalue.lt.100) then
        write (str,'(i2)') ivalue
          str = '000'//str(1:2)
        inum = 5
      else if (ivalue.lt.1000) then
        write (str,'(i3)') ivalue
          str = '00'//str(1:3)
        inum = 5
      else if (ivalue.lt.10000) then
        write (str,'(i4)') ivalue
          str = '0'//str(1:4)
        inum = 5
      else if (ivalue.lt.100000) then
        write (str,'(i5)') ivalue
        inum = 5
      else if (ivalue.lt.1000000) then
        write (str,'(i6)') ivalue
        inum = 6
      else if (ivalue.lt.10000000) then
        write (str,'(i7)') ivalue
        inum = 7
      else if (ivalue.lt.100000000) then
        write (str,'(i8)') ivalue
        inum = 8
      else if (ivalue.lt.1000000000) then
        write (str,'(i9)') ivalue
        inum = 9
      else
        write (str,'(i10)') ivalue
        inum = 10
      endif

      return
      end 
