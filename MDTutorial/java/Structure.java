import java.io.*;

public class Structure extends Super
{
	/*
*********************************************************************

dl_poly/java class to define the structure of a configuration
in terms of molecules and other components

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
    Config cfg;
    int natms,nmols,nbnds,nmoltp,nunq,nrept,nshl,ntatm,imcon;
    int[] lbnd,idc,ist,isz,mtp,mst,msz;
    int[][] join,bond;
    Element[] atoms;
    Molecule[] molecules;
    double[][] xyz;
    double[] cell,chge;
    String[] name,unqatm;

    Structure(Config con)
    {
	/*
*********************************************************************

dl_poly/java constructor for ForceField class

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	cfg=con;
	natms=cfg.natms;
	imcon=cfg.imcon;
	cell=cfg.cell;
	xyz=cfg.xyz;
	chge=cfg.chge;
	atoms=cfg.atoms;
	nmols=0;

	// identify unique atoms

	nunq=uniqueAtoms();

	// atomic repeat pattern in CONFIG file

	nrept=numRepeatAtoms();

	// number of core-shell units

	nshl=numCoreShells();
	
	// number of atoms types

	ntatm=numAtomTypes();

	// identify molecules and molecular types in configuration

	if(cfg.nbnds == 0) cfg.getBonds();
	nbnds=cfg.nbnds;
	bond=cfg.bond;
	lbnd=cfg.lbnd;
	join=cfg.join;
	if(nbnds > 0) nmols=molFind();
	if(nmols > 0) nmoltp=molSame();
	if(nmols ==0) nmols=natms/nrept;

    }

    int molFind()
    {
	/*
*********************************************************************

dl_poly/java routine for identifying individual molecules in a 
configuration using cluster analysis

author: w.smith september 2001
copyright: daresbury laboratory

*********************************************************************
*/
	int nmols,icl,jcl,kcl,iatm,jatm,idm,nmax;

	idc=new int[natms];

    // initialise cluster arrays
      
	for (int i=0;i<natms;i++)
	    {
		idc[i]=-1;
	    }
      
	// search for molecules

	nmols=0;
	for(int i=0;i<nbnds;i++)
	    {
		iatm=Math.min(join[0][i],join[1][i]);
		jatm=Math.max(join[0][i],join[1][i]);

		icl=idc[iatm];
		if(icl < 0)
		    {
			icl=iatm;
			idc[iatm]=iatm;
			nmols++;
		    }

		if(idc[jatm] < 0)
		    {
			idc[jatm]=icl;
		    }
		else if (idc[jatm] != icl)
		    {
			jcl=Math.min(icl,idc[jatm]);
			kcl=Math.max(icl,idc[jatm]);
                  
			for(int k=0;k<natms;k++)
			    {
				if(idc[k] == kcl) idc[k]=jcl;
			    }
			icl=jcl;
			nmols--;
		    }

	    }
	
	for (int i=0;i<natms;i++)
	    {
		if(idc[i]<0) nmols++;
	    }

	monitor.println("Number of molecules found: "+BML.fmt(nmols,6));
	if(nmols == 0)return 0;

    // define  molecule locations in CONFIG arrays
      
	ist=new int[nmols];
	isz=new int[nmols];

	for(int i=0;i<nmols;i++)
	    {
		ist[i]=0;
		isz[i]=0;
	    }

	idm=0;
	for(int i=1;i<natms;i++)
	    {
		if(idc[i] != idc[i-1])
		    {	
			isz[idm]=i-ist[idm];
			idm++;
			if(idm == nmols)
			    {
				monitor.println("Error - molecule data not contiguous in CONFIG");
				return -1;
			    }
			ist[idm]=i;
		    }
	    }
	isz[idm]=natms-ist[idm];
	idm++;
    
	nmax=1;
	for(int i=0;i<nmols;i++)
	    {
		nmax=Math.max(nmax,isz[i]);
	    }
	monitor.println("Largest molecule found: "+BML.fmt(nmax,6));

	return nmols;
    }
			
    int molSame()
    {
	/*
*********************************************************************

dl_poly/java routine for identifying sequences of identical molecules
in a configuration (NB this operation valid for dl_poly only)

author: w.smith september 2001
copyright: daresbury laboratory

*********************************************************************
*/
	String molnam;
	boolean same,found;
	int ia,ib,kk,mm,idm,kkk,ja,jb;

	if(nmols == 0) return 0;

	mtp=new int[nmols];
	mst=new int[nmols];
	msz=new int[nmols];

	for(int i=0;i<nmols;i++)
	    {
		mtp[i]=i;
		mst[i]=0;
		msz[i]=0;
	    }
	// now identify equivalent molecules 

	for(int i=1;i<nmols;i++)
	    {
		// check molecule size
		same=true;
	    OUT:
		if(isz[i] == isz[i-1])
		    {
			ia=ist[i];
			ib=ist[i-1];
			// compare corresponding atoms
			for(int j=0;j<isz[i];j++)
			    {
				// check atom types
				if(!(atoms[ia+j].zsym.equals(atoms[ib+j].zsym)))
				    {
					same=false;
					break OUT;
				    }
				// check atom valencies
				if(lbnd[ia+j] != lbnd[ib+j])
				    {
					same=false;
					break OUT;
				    }
			    }
			// check molecular topology

			found=false;
		    OUT1:
			for(int j=0;j<isz[i];j++)
			    {
				kk=ia+j;
				mm=ib+j;
				for(int k=0;k<lbnd[kk];k++)
				    {
					kkk=bond[k][kk]-ist[i];
					for(int m=0;m<lbnd[mm];m++)
					    {
						if(kkk == bond[m][mm]-ist[i-1])
						    {
							found=true;
						    }
					    }
					if(!found) break OUT1;
				    }
			    }
			same=found;
			if(!same) break OUT;
			mtp[i]=mtp[i-1];
		    }
	    }

	// make list of unique molecules

	idm=0;
	mst[0]=0;
	for(int i=1;i<nmols;i++)
	    {
		if(mtp[i] != mtp[i-1])
		    {	
			msz[idm]=i-mst[idm];
			mst[++idm]=i;
		    }
	    }
	msz[idm]=nmols-mst[idm];
	nmoltp=++idm;

	// assign structural details to molecules

	molecules=new Molecule[nmoltp];

	for(int i=0;i<nmoltp;i++)
	    {
		ja=ist[mst[i]];
		jb=isz[mst[i]];
		molnam="Species "+BML.fmt(i,6);
		monitor.println("Analysing structure of "+molnam);
		molecules[i]=new Molecule(ja,jb,imcon,molnam,atoms,cell,chge,xyz);
	    }

	return nmoltp;
    }

    int uniqueAtoms()
    {
	/*
*********************************************************************

dl_poly/java class to determine unique atoms in a configuration

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	int kkk,mxunq=10;
	boolean lnew;

	name=new String[natms];
	unqatm=new String[mxunq];

	for(int i=0;i<natms;i++)
	    {
		name[i]=BML.fmt(atoms[i].zsym,8);
	    }

	// determine unique atom types

	kkk=1;
	unqatm[0]=name[0];

	for(int i=1;i<natms;i++)
	    {
		lnew=true;
		for(int j=0;j<kkk;j++)
		    {
			if(name[i].equals(unqatm[j])) lnew=false;
		    }
		if(lnew)
		    {
			if(kkk == mxunq)
			    {
				String unqnew[]=new String[2*mxunq];
				System.arraycopy(unqatm,0,unqnew,0,mxunq);
				unqatm=unqnew;
				mxunq*=2;
			    }
			unqatm[kkk]=name[i];
			kkk++;
		    }
	    }
	return kkk;
    }

    int numRepeatAtoms()
    {
	/* 
*********************************************************************

dl_poly/java routine to determine the repeat pattern of atoms in a
CONFIG file

copyright daresbury laboratory
author w.smith june 2001

**********************************************************************
*/
	boolean lnum;
	int nnum;

	nnum=1;
	nrept=natms;

    OUT1:
	for(int j=1;j<=natms/2;j++)
	    {
		if(natms%j == 0)
		    {
			lnum=true;
			for(int i=0;i<natms-j;i++)
			    {            
				if(!name[i].equals(name[i+j]))lnum=false;
			    }
			if(lnum)
			    {
				nrept=j;
				nnum=natms/j;
				break OUT1;
			    }
		    }
	    }
	return nrept;
    }

    int numCoreShells()
	{
	    /*
**********************************************************************

dl_poly/java routine to determine number of core-shell units in a
configuration

copyright daresbury laboratory
author w.smith june 2001

**********************************************************************
*/

	    nshl=0;
	    for(int i=0;i<nrept;i++)
		{
		    if(name[i].charAt(4)=='_' && name[i].toLowerCase().charAt(5)=='s')nshl++;
		}
	    return nshl;

	}

    int numAtomTypes()
	{
	    /*
**********************************************************************

dl_poly/java routine to determine number of atom types in a configuration
allowing for core-shell units

copyright daresbury laboratory
author w.smith june 2001

**********************************************************************
*/

	ntatm=0;
	for(int i=0;i<nunq;i++)
	    {
		if(unqatm[i].substring(4).equals("    "))
		    {
			ntatm++;
		    }
		else if(unqatm[i].toLowerCase().substring(4,6).equals("_s"))
		    {
			ntatm++;
		    }
	    }
	return ntatm;

	}

}
