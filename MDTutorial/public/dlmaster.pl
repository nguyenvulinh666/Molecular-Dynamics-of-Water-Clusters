#!/usr/bin/perl
# begining of main script
$DATE_OF_REVISION="Tue Feb 18 14:48:58 GMT 2003";

$PREAMBLE[0]=" Program: dlmaster.pl\n";
$PREAMBLE[1]=" Author: Changman Moon <c.moon".'@'."warwick.ac.uk>\n";
$PREAMBLE[2]=" Update: $DATE_OF_REVISION\n";
$PREAMBLE[3]=" Remark: This perl script processes the DL_POLY source files and writes a new makefile, having original Makefile renamed to Makefile.original, for implimentation of multiple DL_POLYs running simultanious with one MPI call. Once this script is run, it will produce neccessary files and make changes. It doesn't alter original source codes\n\n";
$PREAMBLE[4]=" Usage: make <object> \n\n";

print "\n@PREAMBLE\n";

write_makefile(); 
write_dlmaster();
make_groupcomm();
write_ilp();
write_readme();
write_master();

# end of main script

sub write_makefile{
    if ( -e "Makefile" ){
	if ( !-e "Makefile.original" ){
	    unless ( rename("Makefile","Makefile.original")){
		print "Error - fail to rename the original Makefile to Makefile.original ...\n";
		exit;
	    }
	}
	unless ( open (MAKEFILE,"Makefile.original") ){
	    print "Error - fail to open Makefile.original ...\n";exit;
	}else {
	    unless ( open (DLMASTER,"> Makefile")){
		print "Error - fail to make Makefile.dlmaster ...\n";exit;
	    }else{
		$block_found = 0;
		while (<MAKEFILE>) {
		    if (/^# Author:/){
			print DLMASTER "# DLMASTER==============================================================\n";
			print DLMASTER "# Revised: Changman Moon, Feburary 2003          \n";
			print DLMASTER "#          Multiple DL_POLY tasks in one MPI call\n";
			print DLMASTER "#          Dept of Chemistry, Warwick University \n";
			print DLMASTER "# Action:  Original Makefile will be renamed to Makefile.original\n";
			print DLMASTER "# ======================================================================\n";
		    }
		    if (/clean:/i){
			$block_found = 1;
		    }
		    if (/dlpoly.o/){
			$_=~s/dlpoly.o/dlmaster.o dlpoly.o/;
		    }
		    if (/DPP/ && /.dpp/){
			$_="ILP = ./ilp\n$_";
		    }
		    if (/EX/ && /DLPOLY.X/){
			$_=~s/DLPOLY.X/DLPOLY.MASTER/;
		    }
		    if ($block_found && /\*.tmp.f/ && /mpif.h/){
			chomp($_);
			$_="$_ *.dlm.f\n";
		    }
		    if ($block_found && /^#=====/ ){
			print DLMASTER "$_";
			print DLMASTER "# Clean up dlmaster files\n";
			print DLMASTER "clean_dlmaster:\n";
			print DLMASTER "\trm -f groupcomm.h dlmaster.f ilp README MASTER\n\n";
			$block_found=0;
		    }
		    if (/.f/ && /.tmp.f/ && /DPP/ && /CPFLAGS/){
			$_=~s/>/| \$(ILP) >/;
		    }
		    print DLMASTER "$_";
		}
	    }
	}
    }else {
	print "Error - Makefile doesn't exist ...\n";exit;
    }
}

sub make_groupcomm {
    unless (open(GROUPCOMM,">groupcomm.h")){
	print "Error - fail to make groupcomm.h ...\n";exit;
    }else{
	print GROUPCOMM "      INTEGER GROUP_COMM_WORLD\n";
	print GROUPCOMM "      COMMON /GroupMpiComm/ GROUP_COMM_WORLD\n";
    }
}

sub write_ilp {
    unless(open(ILP,">ilp")){
	print "Error - fail to produce ilp - perl in-line-processor ...\n";exit;
    }else{
    print ILP <<END_OF_ILP;
#!/usr/bin/perl -w
#Perl script of in-line-processor to edit DLPOLY source files for DLMASTER

\$i=10;\$file[0]=\$i;
\$dlpoly_flg=0;
\$file_opening=0;
\$file_skip=0;
\$file[1]=\$dlpoly_flg;
\$file[2]=\$file_opening;
\$file[3]=\$file_skip;

while(<>){
    if(/program/i && /dlpoly/){
	\$dlpoly_flg=1;\$file[1]=1;
    }
    if(/'CONFIG'/||/'CONTROL'/||/'FIELD'/||/'TABLE'/||/'REVOLD'/||/'REVCON'/||/'REVIVE'/||/'HISTORY'/||/'STATIS'/||/'RDFDAT'/||/'ZDNDAT'/||/'OUTPUT'/){
	\$file_opening=1;\$file[2]=1;
    }
    if(/program/i && /dlmaster/){
	\$file_skip=1;\$file[3]=1;
    }
    \$file[\$i]=\$_;
    \$i++;
}

write_stdout(\@file);

sub write_stdout {

    my(\$i,\$dlpoly_flg,\$file_opening,\$file_skip);
    \$i=\$_[0];
    \$dlpoly_flg=\$_[1];
    \$file_opening=\$_[2];
    \$file_skip=\$_[3];

    while(defined(\$_[\$i])){

	if(!\$file_skip){
	    \$_[\$i]=~s/MPI_COMM_WORLD/GROUP_COMM_WORLD/;
	}
	if(\$file_opening && \$_[\$i]=~/implicit/){
	    \$_[\$i]="\$_[\$i]\\n";
	    \$_[\$i]="\$_[\$i]      integer MyGroup,GROUP_COMM_WORLD\\n";
	    \$_[\$i]="\$_[\$i]      character(len=40) WorkDir\\n";
	    \$_[\$i]="\$_[\$i]      common /slavegroup/ MyGroup,WorkDir\\n\\n";
	}
	if(\$file_opening && \$_[\$i]=~/include/ && \$_[\$i]=~/dl_params.inc/ && !\$dlpoly_flg){
	    \$_[\$i]="\$_[\$i]\\n";
	    \$_[\$i]="\$_[\$i]      integer MyGroup,GROUP_COMM_WORLD\\n";
	    \$_[\$i]="\$_[\$i]      character(len=40) WorkDir\\n";
	    \$_[\$i]="\$_[\$i]      common /slavegroup/ MyGroup,WorkDir\\n\\n";
	}
	if(\$_[\$i]=~/include/ && \$_[\$i]=~/mpif.h/){
	    \$_[\$i]="\$_[\$i]\\n      include \\\"groupcomm.h\\\"\\n\\n";
	}
	\$_[\$i]=~s/'CONFIG'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'CONFIG'/;
	\$_[\$i]=~s/'CONTROL'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'CONTROL'/;
	\$_[\$i]=~s/'FIELD'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'FIELD'/;
	\$_[\$i]=~s/'TABLE'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'TABLE'/;
	\$_[\$i]=~s/'REVOLD'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'REVOLD'/;
	\$_[\$i]=~s/'REVCON'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'REVCON'/;
	\$_[\$i]=~s/'REVIVE'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'REVIVE'/;
	\$_[\$i]=~s/'HISTORY'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'HISTORY'/;
	\$_[\$i]=~s/'STATIS'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'STATIS'/;
	\$_[\$i]=~s/'RDFDAT'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'RDFDAT'/;
	\$_[\$i]=~s/'ZDNDAT'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'ZDNDAT'/;
	\$_[\$i]=~s/'OUTPUT'/WorkDir\\(1:LEN_TRIM\\(WorkDir\\)\\)\\/\\/'OUTPUT'/;
	if(\$dlpoly_flg && \$_[\$i]=~/program/ && \$_[\$i]=~/dlpoly/){
	    \$_[\$i]=~s/program/subroutine/;
	    \$_[\$i]=~s/dlpoly/dlpoly(GRP_COMM_WORLD,Group,WorkPath)/;
	}
	if(\$dlpoly_flg && \$_[\$i]=~/include/ && \$_[\$i]=~/dl_params.inc/){
	    \$_[\$i]="\$_[\$i]\\n";
	    \$_[\$i]="\$_[\$i]      integer Group,GRP_COMM_WORLD,MyGroup,GROUP_COMM_WORLD\\n";
	    \$_[\$i]="\$_[\$i]      character(len=40) WorkPath,WorkDir\\n";
	    \$_[\$i]="\$_[\$i]      common /GroupMpiComm/ GROUP_COMM_WORLD\\n";
	    \$_[\$i]="\$_[\$i]      common /slavegroup/ MyGroup,WorkDir\\n";
	}
	if(\$dlpoly_flg && \$_[\$i]!~/^c/ && \$_[\$i]!~/if/ && \$_[\$i]!~/do/ && \$_[\$i]=~/      end/){
	    print "      return\\n";
	    print "      end subroutine dlpoly\\n";
	}elsif(\$dlpoly_flg && \$_[\$i]=~/call/ && \$_[\$i]=~/exitcomms/){
	}elsif(\$dlpoly_flg && \$_[\$i]=~/call/ && \$_[\$i]=~/initcomms/){
	}elsif(\$dlpoly_flg && \$_[\$i]=~/data/ && \$_[\$i]=~/safep/ && \$_[\$i]=~/true/){
	    print "\\n";
	    print "      GROUP_COMM_WORLD=GRP_COMM_WORLD\\n";
	    print "      MyGroup=Group\\n";
	    print "      WorkDir=WorkPath\\n";
	    print "\\n";
	}else{
	    if(length(\$_[\$i])>72){
		if(\$_[\$i]=~/open/){
		    \@words=split("//",\$_[\$i]);
		    print "\$words[0]\\n";
		    print "     x //\$words[1]";
		}elsif(\$_[\$i]=~/MPI_ANY_SOURCE/){
		    \@words=split("request",\$_[\$i]);
		    print "\$words[0]\\n";
		    print "     x        request\$words[1]";
		}else{
		    print "\$_[\$i]";
		}
	    }else{
		print "\$_[\$i]";
	    }
	}
	\$i++;
    }
}

END_OF_ILP
}
chmod(0777,"ilp");
}

sub write_dlmaster {
    unless(open(DLMASTER,">dlmaster.f")){
	print "Error - fail to produce dlmaster.f ...\n";exit;
    }else{
    print DLMASTER <<END_OF_DLMASTER;
      PROGRAM dlmaster
c-----------------------------------------------------------------
c     Program: dlmaster.f
c     Author: Changman Moon <c.moon@warwick.ac.uk>
c             Department of Chemistry, University of Warwick
c     Update: $DATE_OF_REVISION
c     Remark: This is a master file that controls multiple DL_POLY 
c             slave tasks running simultaniously.
c-----------------------------------------------------------------

      IMPLICIT NONE

c     use mpi header file
      INCLUDE 'mpif.h'

      INTEGER i,j,k,CharLen
      INTEGER MyId,NumProcs,MyNewId,NumGroupProcs
      INTEGER Group,EachGroup,GRP_COMM_WORLD
      INTEGER NumOfJobs,NumOfRanks,Error
      INTEGER GroupStart,GroupEnd
      INTEGER MyGroupId,NumGroup
      INTEGER, ALLOCATABLE :: GroupRanks(:),AllRanks(:)
      CHARACTER(LEN=40), ALLOCATABLE :: JobDirectory(:)
      CHARACTER GroupDistrib

c     initialise the WORLD MPI communication, my id, size
      CALL MPI_INIT(Error)
      CALL MPI_COMM_RANK( MPI_COMM_WORLD,MyId,Error )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD,NumProcs,Error )


      OPEN (UNIT=7, FILE='MASTER',STATUS='old',ERR=10)
      READ(7,*,END=20,ERR=20) NumOfJobs,GroupDistrib
      IF(NumOfJobs.LE.0) THEN
         IF(MyId.EQ.0) PRINT*,'Error - invalid MASTER records'
         GOTO 20
      END IF
      ALLOCATE (JobDirectory(0:NumOfJobs-1))

      DO i=0,NumOfJobs-1
         READ(7,'(A40)',END=30,ERR=30) 
     &        JobDirectory(i)
         CharLen=LEN_TRIM(JobDirectory(i))
         IF(JobDirectory(i)(CharLen:CharLen).NE.'/') THEN
            CharLen=CharLen+1
            JobDirectory(i)(CharLen:CharLen)='/'
         END IF
      END DO
      
c     check if number of nodes available is valid
      IF(NumProcs.LT.NumOfJobs .OR. MOD(NumProcs,NumOfJobs).NE.0) THEN
         GOTO 40
      END IF

c     allocate groups of ranks arrays
      ALLOCATE (AllRanks(0:NumProcs-1))

      DO i=0,NumProcs-1
         AllRanks(i)=i
      END DO

c     assign groups with appropriate ranks 
      NumOfRanks=NumProcs/NumOfJobs
      ALLOCATE(GroupRanks(NumOfRanks))

c     type of group ranks distribution, which is either distributed or 
c     blockwise. program only checks the first character (d or b).
      IF(GroupDistrib.EQ."d" .OR. GroupDistrib.EQ."D") THEN
         Group=MOD(MyId,NumOfJobs)
         GroupRanks=AllRanks(Group:NumProcs-1:NumOfJobs)
      ELSE IF(GroupDistrib.EQ."b" .OR. GroupDistrib.EQ."B") THEN
         Group=INT(MyId/NumOfRanks)
         GroupStart=Group*NumOfRanks
         GroupEnd=GroupStart+NumOfRanks-1
         GroupRanks=AllRanks(GroupStart:GroupEnd)
      ELSE
         IF(MyId.EQ.0) 
     &   PRINT*,'Error - invalid job distribution type in MASTER'
      END IF

c     PRINT*,'(',Group,')',GroupRanks,'/',AllRanks
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,GRP_COMM_WORLD,Error)
      CALL MPI_GROUP_INCL(GRP_COMM_WORLD,NumOfRanks,GroupRanks,
     &     EachGroup,Error)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,EachGroup,
     &     GRP_COMM_WORLD,Error)
      CALL dlpoly(GRP_COMM_WORLD,Group,JobDirectory(Group))

c     finalise MPI communication
      CALL MPI_FINALIZE(Error)
      
      STOP

 10   CONTINUE
      CALL MPI_FINALIZE(Error)
      IF(MyId.EQ.0) PRINT*,'Error - MASTER file doesn''t exist'
      STOP

 20   CONTINUE
      CLOSE(UNIT=7)
      CALL MPI_FINALIZE(Error)
      STOP

 30   CONTINUE
      CLOSE(UNIT=7)
      CALL MPI_FINALIZE(Error)
      IF(MyId.EQ.0) PRINT*,'Error - invalid number of records in MASTER'
      STOP

 40   CONTINUE
      CALL MPI_FINALIZE(Error)
      IF(MyId.EQ.0) PRINT*, 'Error - invalid number of processes'
      STOP 

      END PROGRAM dlmaster

END_OF_DLMASTER
}
}

sub write_readme{
    unless(open(README,">README")){
	print "Error - fail to produce README ...\n";exit;
    }else{
    print README <<END_OF_README;

README for DLMASTER
===================

Disclaimer of Warranty

This script "dlmaster.pl" should be freely available for academic purposes without communication to the author (C. Moon). No provision is made for comercial usage. The author disclaims all liability for direct or consequential damages resulting from your use of this program and makes no warranties, express or implied, that the codes contained on this program are free of error.

Dr. Changman Moon, Feburary 2003, The department of Chemistry, University of Warwick, Coventry CV4 7AL, U.K.

Copyright(c)

DL_POLY is a molecular dynamics package written by W. Smith, Daresbury Laboratory. 

 
(1) INTRODUCTION

This script assumes you are compiling a MPI version DL_POLY although proper changes to the codes enable other messege passing methods, such as PVM. It DOES NOT ALTER any original DL_POLY source codes in any case except for Makefile. It ONLY processes temporary codes. However, backing up the DL_POLY's original source files is strongly recommended. 

(2) FILES

By executing "dlmaster.pl" under DL_POLY's /source directory, following files will be produced,

    Makefile
    groupcomm.h
    dlmaster.f
    ilp
    README
    MASTER

This script "dlmaster.pl" written in PERL generates a few files and description for the files are ,\n

Makefile - a revised version of Makefile and the original Makefile is saved as Makefile.original.

groupcomm.h - a header file defining MPI group communicator.

dlmaster.f - a master Fortran 90 source file which deals with initiation, distribution, and termination of MPI nodes to the given number of DL_POLY tasks.

ilp - a in_line_processor that edits neccessary DL_POLY source files immeadiately after "dpp" preprocessor taking standard inputs and redirect it's outputs to standard output.

README - this document.

MASTER - a sample task control file.

(3) USAGE
Once you run this script within the DL_POLY /source directory, it will produce a number of files listed above. As usual, compile the DL_POLY sources with "make <object>".

Instead of the usual executable file name, DLPOLY.X, it generates an executable called DLPOLY.MASTER.

You have to write a job file called MASTER that contains a number specifying the number of DL_POLY tasks and that number of directory paths in separate lines. The first record also contains the type of group ranks distribution such as blockwise groupping or distributed groupping Refer to the following example, suppose you have four DL_POLY tasks, running DLPOLY.MASTER on the local directory "./", the task file MASTER must reside on the same local directory and will look like, 

---- sample MASTER file ----
4, blockwise
./task1/
./task2/
./task3/
./task4/
----------------------------

where each directory must contain neccessary DL_POLY input files such as CONFIG, CONTROL, and FIELD at least. The number of nodes available after calling MPI must be no less than the number of task and a multiple of the number of tasks submitted. No blank or comment lines must be present at the begining or between the lines in MASTER. 

(4) RECOVERY OF ORIGINAL STATE (clean, clean_dlmaster)

By executing "make clean", all the temporary and object files will be deleted as usual. Executing "make clean_dlmaster" will delete DLMASTER files excluding Makefile. Therefore, after "make clean_dlmaster" rename Makefile.original to Makefile will restate original source files.

Good Luck!!
-Changman Moon

END_OF_README
}
}

sub write_master{
    unless(open(MASTER,">MASTER")){
	print "Error - fail to produce MASTER ...\n";exit;
    }else{
    print MASTER <<END_OF_MASTER;
4
./task1/
./task2/
./task3/
./task4/
END_OF_MASTER
}
}
