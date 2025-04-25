import java.io.*;
import java.awt.*;

public class Config extends Super
{
	/*
*********************************************************************

dl_poly/java class to define a molecular configuration

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
    GUI home;
    String title;
    int mxatms=100,mxjoin=100,mxcon=20;
    int natms,nmols,nbnds,nmoltp,imcon,levcfg,nvrt,nedge;
    int[] idc,ist,isz,mtp,mst,msz,lbnd;
    int[][] edge,bond,join;
    double[] cell,weight,chge;
    double[][] xyz,vrt;
    Element[] atoms;

    Config()
    {
	/*
*********************************************************************

dl_poly/java constructor to define Config class

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	setValues();
    }

    Config(int newatms)
    {
	/*
*********************************************************************

dl_poly/java constructor to define Config class

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	mxatms=newatms;
	setValues();
    }

    Config(GUI here,String ftype)
    {
	/*
*********************************************************************

dl_poly/java constructor to define Config class

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	home=here;
	natms=0;
	imcon=0;
	levcfg=0;
	title="";
	cell=new double[9];
	xyz=new double[3][mxatms];
	atoms=new Element[mxatms];
	weight=new double[mxatms];
	chge=new double[mxatms];

	if(ftype.equals("CFG"))
	    {
		monitor.println("Select required CONFIG file for input");
		if((fname=selectFileNameBegins(home,"CFG"))!=null)
		    {
			if(!rdCFG(fname)) return;
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if(ftype.equals("XYZ"))
	    {
		monitor.println("Select required XYZ file for input");
		if((fname=selectFileNameEnds(home,"XYZ"))!=null)
		    {
			if(rdXYZ(fname))
			    {
				fname="CFGXYZ."+numxyz;
				if(configWrite(fname)) numxyz++;
			    }
			else
			    return;
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if(ftype.equals("SEQ"))
	    {
		monitor.println("Select required SEQNET file for input");
		if((fname=selectFileNameEnds(home,"SEQ"))!=null)
		    {
			if(rdSEQ(fname))
			    {
				fname="CFGSEQ."+numseq;
				if(configWrite(fname)) numseq++;
			    }
			else
			    return;
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if(ftype.equals("MSI"))
	    {
		monitor.println("Select required CERIUS/MSI file for input");
		if((fname=selectFileNameEnds(home,"MSI"))!=null)
		    {
			if(rdMSI(fname))
			    {
				fname="CFGMSI."+nummsi;
				if(configWrite(fname)) nummsi++;
			    }
			else
			    return;
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if(ftype.equals("MDR"))
	    {
		monitor.println("Select required MDR file for input");
		if((fname=selectFileNameEnds(home,"MDR"))!=null)
		    {
			if(rdMDR(fname))
			    {
				fname="CFGMDR."+nummdr;
				if(configWrite(fname)) nummdr++;
			    }
			else
			    return;
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if(ftype.equals("CONFIG"))
	    {
		if(!rdCFG(ftype)) return;
	    }
	else if(ftype.equals("REVCON"))
	    {
		if(!rdCFG(ftype)) return;
	    }
	setBox();
	
    }

    void setValues()
    {
	/*
*********************************************************************

dl_poly/java routine to define Config variables

copyright - daresbury laboratory 2001
author    - w.smith june 2001

*********************************************************************
*/
	natms=0;
	nbnds=0;
	nvrt=0;
	nedge=0;
	imcon=0;
	levcfg=0;
	nmols=0;
	nmoltp=0;
	title="";
	cell=new double[9];
	xyz=new double[3][mxatms];
	atoms=new Element[mxatms];
	weight=new double[mxatms];
	chge=new double[mxatms];
	lbnd=new int[mxatms];
	bond=new int[mxcon][mxatms];
	join=new int[2][mxjoin];
	for(int i=0;i<9;i++)
	    cell[i]=0;
	for(int i=0;i<mxatms;i++)
	    lbnd[i]=0;
    }

    boolean configWrite(String filename)
    {
	/*
*********************************************************************

dl_poly/java routine to write a DL_POLY CONFIG file

copyright - daresbury laboratory
author    - w.smith october 2000

*********************************************************************
*/
	try
	    {
		DataOutputStream outStream = new DataOutputStream(new FileOutputStream(filename));
		outStream.writeBytes(title+"\n");
		outStream.writeBytes(BML.fmt(0,10)+BML.fmt(imcon,10)+"\n");
		if(imcon>0)
		    {
			outStream.writeBytes(BML.fmt(cell[0],20)+BML.fmt(cell[1],20)+BML.fmt(cell[2],20)+"\n");
			outStream.writeBytes(BML.fmt(cell[3],20)+BML.fmt(cell[4],20)+BML.fmt(cell[5],20)+"\n");
			outStream.writeBytes(BML.fmt(cell[6],20)+BML.fmt(cell[7],20)+BML.fmt(cell[8],20)+"\n");
		    }
		for (int k=0;k<natms;k++)
		    {
			outStream.writeBytes(BML.fmt(atoms[k].zsym,8)+BML.fmt(k+1,10)+BML.fmt(atoms[k].znum,10)+"\n");
			outStream.writeBytes(BML.fmt(xyz[0][k],20)+BML.fmt(xyz[1][k],20)+BML.fmt(xyz[2][k],20)+"\n");
		    }
		outStream.close();
		monitor.println("New CONFIG file created: "+filename);
	    }
	catch(Exception e)
	    {
		monitor.println("Error - writing file: "+filename);
		return false;
	    }
	return true;
    }

    boolean rdCFG(String fname)
    {
	/*
*********************************************************************

dl_poly/java routine to read a DL_POLY CONFIG file

copyright - daresbury laboratory
author    - w.smith june 2001

*********************************************************************
*/
        int i,j,k,m;
        LineNumberReader lnr=null;
        String record="",namstr="";
        double xlo=0,xhi=0,ylo=0,yhi=0,zlo=0,zhi=0,ddd;

        for (i=0;i<9;i++)
            cell[i]=0.0;

	// open the CONFIG file

        try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		title = lnr.readLine();
		monitor.println("File header record: "+title);
		record = lnr.readLine();
		levcfg=BML.giveInteger(record,1);
		imcon =BML.giveInteger(record,2);
		if(imcon > 0)
		    {
			record = lnr.readLine();
			cell[0]=BML.giveDouble(record,1);
			cell[1]=BML.giveDouble(record,2);
			cell[2]=BML.giveDouble(record,3);
			record = lnr.readLine();
			cell[3]=BML.giveDouble(record,1);
			cell[4]=BML.giveDouble(record,2);
			cell[5]=BML.giveDouble(record,3);
			record = lnr.readLine();
			cell[6]=BML.giveDouble(record,1);
			cell[7]=BML.giveDouble(record,2);
			cell[8]=BML.giveDouble(record,3);
		    }
	    
		// read configuration
	    
		i=0;
		j=0;
		k=levcfg+2;
		while((record=lnr.readLine()) != null)
		    {
			i=j/k;
			m=j-k*i;
			if(m == 0)
			    {
				namstr = BML.fmt(BML.giveWord(record,1),8);

				if(i == mxatms)
				    {
					double chg[]=new double[2*mxatms];
					double wgt[]=new double[2*mxatms];
					Element atomz[]=new Element[2*mxatms];
					double uvw[][]=new double[3][2*mxatms];
					for(int n=0;n<mxatms;n++)
					    {
						chg[n]=chge[n];
						wgt[n]=weight[n];
						atomz[n]=new Element(atoms[n].zsym);
						uvw[0][n]=xyz[0][n];
						uvw[1][n]=xyz[1][n];
						uvw[2][n]=xyz[2][n];
					    }
					atoms=atomz;
					xyz=uvw;
					chge=chg;
					weight=wgt;
					mxatms*=2;
				    }
			    
				if(!namstr.equals(""))
				    {
					atoms[i]=new Element(namstr);
					weight[i]=atoms[i].zmas;
					chge[i]=0.0;
				    }
				else
				    {
					monitor.println("Error  - unknown atom type in CONFIG file: "+namstr);
					lnr.close();
					return false;
				    }
			    }
			if(m == 1)
			    {
				xyz[0][i]=BML.giveDouble(record,1);
				xyz[1][i]=BML.giveDouble(record,2);
				xyz[2][i]=BML.giveDouble(record,3);
				if(imcon == 0)
				    {
					xlo=Math.min(xlo,xyz[0][i]);
					ylo=Math.min(ylo,xyz[1][i]);
					zlo=Math.min(zlo,xyz[2][i]);
					xhi=Math.max(xhi,xyz[0][i]);
					yhi=Math.max(yhi,xyz[1][i]);
					zhi=Math.max(zhi,xyz[2][i]);
				    }
			    }
			j++;
		    }
		lnr.close();
	    }
        catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return false;
	    }
        catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return false;
	    }
	natms=++i;
        if(imcon == 0)
	    {
		//imcon=2;
		ddd=Math.pow((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms,(1.0/3.0));
		cell[0]=-(xlo-xhi)+ddd;
		cell[4]=-(ylo-yhi)+ddd;
		cell[8]=-(zlo-zhi)+ddd;
	    }
	else
	    {
		images(natms,imcon,cell,xyz);
	    }
        monitor.println("Selected CONFIG file loaded successfully");
        monitor.println("Number of atoms found: "+natms);
	return true;
    }

    boolean rdXYZ(String fname)
    {
	/*
*********************************************************************
     
dl_poly/java routine to read an XYZ configuration file
     
copyright - daresbury laboratory
author    - w.smith june 2001
     
*********************************************************************
*/      
	LineNumberReader lnr=null;
	String record="",namstr="";
	double xlo,xhi,ylo,yhi,zlo,zhi,ddd;

        for (int i=0;i<9;i++)
            cell[i]=0.0;

	// open XYZ file

	try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		record = lnr.readLine();
		natms = BML.giveInteger(record,1);
		monitor.println("Number of atoms in file : "+natms);
		if(natms > mxatms)
		    {
			mxatms=natms;
			atoms=new Element[mxatms];
			xyz=new double[3][mxatms];
			chge=new double[mxatms];
			weight=new double[mxatms];
		    }
		title = lnr.readLine();
		monitor.println("File header record: "+title);

		// read configuration
		
		xlo=0.0;
		ylo=0.0;
		zlo=0.0;
		xhi=0.0;
		yhi=0.0;
		zhi=0.0;
		
		for(int i=0;i<natms;i++)
		    {
			record = lnr.readLine();
			namstr = BML.fmt(BML.giveWord(record,1),8);
			if(!namstr.equals(""))
			    atoms[i]=new Element(namstr);
			else
			    {
				monitor.println("Error  - unknown atom type in XYZ file: "+namstr);
				lnr.close();
				return false;
			    }
			xyz[0][i]=BML.giveDouble(record,2);
			xyz[1][i]=BML.giveDouble(record,3);
			xyz[2][i]=BML.giveDouble(record,4);
			xlo=Math.min(xlo,xyz[0][i]);
			ylo=Math.min(ylo,xyz[1][i]);
			zlo=Math.min(zlo,xyz[2][i]);
			xhi=Math.max(xhi,xyz[0][i]);
			yhi=Math.max(yhi,xyz[1][i]);
			zhi=Math.max(zhi,xyz[2][i]);
		    }
		ddd=Math.pow(((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms),(1.0/3.0));
		cell[0]=-(xlo-xhi)+ddd;
		cell[4]=-(ylo-yhi)+ddd;
		cell[8]=-(zlo-zhi)+ddd;
		lnr.close();
	    }
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return false;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return false;
	    }
	if(imcon > 0)
	    {
		images(natms,imcon,cell,xyz);
	    }
	monitor.println("Selected XYZ file loaded successfully");
        levcfg=0;
        //imcon=2;

	return true;
    }

    boolean rdSEQ(String fname)
    {
	/*
*********************************************************************
     
dl_poly/java routine to read a SEQNET configuration file
     
copyright - daresbury laboratory
author    - w.smith june 2001
   
*********************************************************************
*/      
	LineNumberReader lnr=null;
	String record="",namstr="";
	String label,header,atom,hetatm,seqres;
	double x,y,z,xlo,xhi,ylo,yhi,zlo,zhi,ddd;

	header="HEADER";
	atom="ATOM";
	hetatm="HETATM";
	seqres="SEQRES";

	natms=0;
	for (int i=0;i<9;i++)
	    cell[i]=0.0;
    
    // open SEQNET file
    
	try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
	    
		// read configuration
	    
		xlo=0.0;
		ylo=0.0;
		zlo=0.0;
		xhi=0.0;
		yhi=0.0;
		zhi=0.0;
	    
		while((record=lnr.readLine()) != null)
		    {
			if(record.indexOf(header)>=0)
			    {
				title=record.substring(10,71);
				monitor.println("File header record: "+title);
			    }
			else if(record.indexOf(atom)>=0 || record.indexOf(hetatm)>=0)
			    {
				if(natms == mxatms)
				    {
					double chg[]=new double[2*mxatms];
					double wgt[]=new double[2*mxatms];
					Element atomz[]=new Element[2*mxatms];
					double uvw[][]=new double[3][2*mxatms];
					for(int n=0;n<mxatms;n++)
					    {
						chg[n]=chge[n];
						wgt[n]=weight[n];
						atomz[n]=new Element(atoms[n].zsym);
						uvw[0][n]=xyz[0][n];
						uvw[1][n]=xyz[1][n];
						uvw[2][n]=xyz[2][n];
					    }
					atoms=atomz;
					xyz=uvw;
					chge=chg;
					weight=wgt;
					mxatms*=2;
				    }
				label=new String(record);
				record=lnr.readLine().substring(12);
				namstr=BML.fmt(BML.giveWord(record,1).trim(),8);
				record=record.substring(18);
				x=BML.giveDouble(record,1);
				y=BML.giveDouble(record,2);
				z=BML.giveDouble(record,3);
				if(label.indexOf(hetatm)>=0 && namstr.equals("O       ")) namstr="OW      ";
				atoms[natms]=new Element(namstr);
				weight[natms]=atoms[natms].zmas;
				chge[natms]=0.0;
				xyz[0][natms]=x;
				xyz[1][natms]=y;
				xyz[2][natms]=z;
				xlo=Math.min(xlo,x);
				ylo=Math.min(xlo,y);
				zlo=Math.min(xlo,z);
				xhi=Math.max(xhi,x);
				yhi=Math.max(yhi,y);
				zhi=Math.max(zhi,z);
				natms++;
			    }
		    }
		lnr.close();
	    }
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return false;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return false;
	    }
	ddd=Math.pow(((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms),(1.0/3.0));
	cell[0]=-(xlo-xhi)+ddd;
	cell[4]=-(ylo-yhi)+ddd;
	cell[8]=-(zlo-zhi)+ddd;
	if(imcon > 0)
	    {
		images(natms,imcon,cell,xyz);
	    }
	monitor.println("Selected SEQNET file loaded successfully");
	levcfg=0;
	//imcon=2;

	return true;
    }

    boolean rdMSI(String fname)
    {
	/*
*********************************************************************

dl_poly/java utility to read a CERIUS 2 configuration file

CERIUS 2 is the copyright of Molecular Simulations Inc

copyright daresbury laboratory
written by w.smith  june 2001

********************************************************************
*/      
        LineNumberReader lnr=null;
        String record="",namstr="";
        int keypbc,i,j,k,m;
        double xlo=0,xhi=0,ylo=0,yhi=0,zlo=0,zhi=0,ddd;

	natms=0;
        for (i=0;i<9;i++)
            cell[i]=0.0;

	// open the CERIUS file
	
        try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		title=lnr.readLine();
		monitor.println("File header record: "+title);
		while((record=lnr.readLine()) != null)
		    {
			if(record.indexOf("Model")>=0)
			    {
				// do nothing in current implementation
			    }
			else if((m=record.indexOf("PeriodicType"))>=0)
			    {          
				record=record.substring(m+12);
				keypbc=BML.giveInteger(record,1);
				if(keypbc==100)imcon=1;
			    }
			else if((m=record.indexOf("A3"))>=0)
			    {
				record=record.substring(m+2);
				cell[0]=BML.giveDouble(record,1);
				cell[1]=BML.giveDouble(record,2);
				cell[2]=BML.giveDouble(record,3);
			    }
			else if((m=record.indexOf("B3"))>=0)
			    {
				record=record.substring(m+2);
				cell[3]=BML.giveDouble(record,1);
				cell[4]=BML.giveDouble(record,2);
				cell[5]=BML.giveDouble(record,3);
			    }
			else if((m=record.indexOf("C3"))>=0)
			    {
				record=record.substring(m+2);
				cell[6]=BML.giveDouble(record,1);
				cell[7]=BML.giveDouble(record,2);
				cell[8]=BML.giveDouble(record,3);
			    }
			else if(record.indexOf("Atom1")>=0)
			    {
				//Ignore!
			    }
			else if(record.indexOf("Atom2")>=0)
			    {
				//Ignore!
			    }
			else if(record.indexOf("Atom")>=0)
			    {
				if(natms == mxatms)
				    {
					double chg[]=new double[2*mxatms];
					double wgt[]=new double[2*mxatms];
					Element atomz[]=new Element[2*mxatms];
					double uvw[][]=new double[3][2*mxatms];
					for(int n=0;n<mxatms;n++)
					    {
						chg[n]=chge[n];
						wgt[n]=weight[n];
						atomz[n]=new Element(atoms[n].zsym);
						uvw[0][n]=xyz[0][n];
						uvw[1][n]=xyz[1][n];
						uvw[2][n]=xyz[2][n];
					    }
					atoms=atomz;
					weight=wgt;
					xyz=uvw;
					chge=chg;
					mxatms*=2;
				    }
				chge[natms]=0.0;
				weight[natms]=0.0;
				xyz[0][natms]=0.0;
				xyz[1][natms]=0.0;
				xyz[2][natms]=0.0;
				namstr="        ";

			    OUT:
				for(k=0;k<100;k++)
				    {
					record=lnr.readLine();
					record=record.trim();
					if(record.charAt(0)==')')
					    {
						natms++;
						break OUT;
					    }
					else if((m=record.indexOf("FFType"))>=0)
					    {
						record=record.substring(m+6);
						namstr=BML.giveWord(record,1);
						if(namstr.equals("H___A")) namstr="H__HB";
						atoms[natms]=new Element(namstr);
					    }
					else if((m=record.indexOf("Charge"))>=0)
					    {
						record=record.substring(m+6);
						chge[natms]=BML.giveDouble(record,1);
					    }
					else if((m=record.indexOf("Mass"))>=0)
					    {
						record=record.substring(m+4);
						weight[natms]=BML.giveDouble(record,1);
					    }
					else if((m=record.indexOf("XYZ"))>=0)
					    {
						record=record.substring(m+3);
						xyz[0][natms]=BML.giveDouble(record,1);
						xyz[1][natms]=BML.giveDouble(record,2);
						xyz[2][natms]=BML.giveDouble(record,3);
						if(imcon == 0)
						    {
							xlo=Math.min(xlo,xyz[0][natms]);
							ylo=Math.min(ylo,xyz[1][natms]);
							zlo=Math.min(zlo,xyz[2][natms]);
							xhi=Math.max(xhi,xyz[0][natms]);
							yhi=Math.max(yhi,xyz[1][natms]);
							zhi=Math.max(zhi,xyz[2][natms]);
						    }
					    }
				    }
			    }
		    }
		lnr.close();
	    }
        catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return false;
	    }
        catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		monitor.println(record);
		return false;
	    }
        if(imcon == 0)
	    {
		//imcon=2;
		ddd=Math.pow((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms,(1.0/3.0));
		cell[0]=-(xlo-xhi)+ddd;
		cell[4]=-(ylo-yhi)+ddd;
		cell[8]=-(zlo-zhi)+ddd;
	    }
	else
	    {
		images(natms,imcon,cell,xyz);
	    }
        monitor.println("CERIUS file loaded successfully");
        monitor.println("Number of atoms found: "+natms);

	return true;
    }

    boolean rdMDR(String fname)
    {
	/*
*********************************************************************
     
dl_poly/java routine to read an MDR configuration file
     
copyright - daresbury laboratory
author    - w.smith june 2001
     
*********************************************************************
*/      
	LineNumberReader lnr=null;
	String record="",namstr="";
	double xlo,xhi,ylo,yhi,zlo,zhi,ddd;

        for (int i=0;i<9;i++)
            cell[i]=0.0;

	// open MDR file

	try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		title = lnr.readLine();
		monitor.println("File header record: "+title);
		record = lnr.readLine();
		natms = BML.giveInteger(record,1);
		monitor.println("Number of atoms in file : "+natms);
		if(natms > mxatms)
		    {
			mxatms=natms;
			atoms=new Element[mxatms];
			xyz=new double[3][mxatms];
			chge=new double[mxatms];
			weight=new double[mxatms];
		    }

		// read configuration
		
		xlo=0.0;
		ylo=0.0;
		zlo=0.0;
		xhi=0.0;
		yhi=0.0;
		zhi=0.0;
		
		for(int i=0;i<natms;i++)
		    {
			record = lnr.readLine();
			namstr = BML.fmt(BML.giveWord(record,1),8);
			if(!namstr.equals(""))
			    atoms[i]=new Element(namstr);
			else
			    {
				monitor.println("Error  - unknown atom type in MDR file: "+namstr);
				lnr.close();
				return false;
			    }
			record = lnr.readLine();
			xyz[0][i]=BML.giveDouble(record,1);
			xyz[1][i]=BML.giveDouble(record,2);
			xyz[2][i]=BML.giveDouble(record,3);
			xlo=Math.min(xlo,xyz[0][i]);
			ylo=Math.min(ylo,xyz[1][i]);
			zlo=Math.min(zlo,xyz[2][i]);
			xhi=Math.max(xhi,xyz[0][i]);
			yhi=Math.max(yhi,xyz[1][i]);
			zhi=Math.max(zhi,xyz[2][i]);
		    }
		ddd=Math.pow(((xhi-xlo)*(yhi-ylo)*(zhi-zlo)/natms),(1.0/3.0));
		cell[0]=-(xlo-xhi)+ddd;
		cell[4]=-(ylo-yhi)+ddd;
		cell[8]=-(zlo-zhi)+ddd;
		lnr.close();
	    }
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return false;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return false;
	    }
	if(imcon > 0)
	    {
		images(natms,imcon,cell,xyz);
	    }
	monitor.println("Selected MDR file loaded successfully");
        levcfg=0;
        //imcon=2;

	return true;
    }

    void switchBonds()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	if(nbnds==0)
	    {
		nbnds=getBonds();
	    }
	else
	    {
		nbnds=0;
	    }
    }

    int getBonds()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	int i,j,m,m1,m2,last;
	double xd,yd,zd,ff,rsq,radi,radj,bndfac;

	double uvw[]=new double[3];
	lbnd=new int[mxatms];
	bond=new int[mxcon][mxatms];

	nbnds=0;
	last=natms;
	m1=natms/2;
	m2=(natms-1)/2;
	bndfac=1.0+bondpc/100.0;
	
	for(i=0;i<natms;i++)
	    lbnd[i]=0;
	
	for(m=1;m<=m1;m++)
	    {
		if(m > m2)last=m1;
		
		for(i=0;i<last;i++)
		    {
			j=i+m;
			if(j>=natms)j=j-natms;
			radi=atoms[i].zrad;
			radj=atoms[j].zrad;
			ff=bndfac*(radi+radj);
			uvw[0]=xyz[0][i]-xyz[0][j];
			uvw[1]=xyz[1][i]-xyz[1][j];
			uvw[2]=xyz[2][i]-xyz[2][j];
			image(imcon,cell,uvw);
			rsq=uvw[0]*uvw[0]+uvw[1]*uvw[1]+uvw[2]*uvw[2];
			if(rsq <= ff*ff)
			    {
				if((lbnd[i]==mxcon)||(lbnd[j]==mxcon))
				    {
					int bnd[][]=new int[2*mxcon][mxatms];
					for (int k=0;k<mxatms;k++)
					    {
						for (int l=0;l<mxcon;l++)
							bnd[l][k]=bond[l][k];
					    }
					mxcon*=2;
					bond=bnd;
				    }
				bond[lbnd[i]][i]=j;
				bond[lbnd[j]][j]=i;
				lbnd[i]+=1;
				lbnd[j]+=1;
				nbnds++;
			    }
		    }

	    }
	defineJoin();

	return nbnds;
    }

    void defineJoin()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	int k,n;

	n=0;
	if(nbnds > 0)
	    {
		mxjoin=nbnds;
		join=new int[2][mxjoin];
		for(int i=0;i<natms;i++)
		    {
			for(int j=0;j<lbnd[i];j++)
			    {
				k=bond[j][i];
				if(k < i)
				    {
					join[0][n]=k;
					join[1][n]=i;
					n++;
				    }
			    }
		    }
	    }
    }

    void boxCELL()
    {
	/*
*********************************************************************

dl_poly/java GUI routine - standard periodic boundaries

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	nvrt=8;
	nedge=12;

	vrt=new double[3][nvrt];
	edge=new int[2][nedge];

	vrt[0][0]=-0.5*(cell[0]+cell[3]+cell[6]);
	vrt[1][0]=-0.5*(cell[1]+cell[4]+cell[7]);
	vrt[2][0]=-0.5*(cell[2]+cell[5]+cell[8]);

	vrt[0][1]= 0.5*(cell[0]-cell[3]-cell[6]);
	vrt[1][1]= 0.5*(cell[1]-cell[4]-cell[7]);
	vrt[2][1]= 0.5*(cell[2]-cell[5]-cell[8]);

	vrt[0][2]= 0.5*(cell[0]+cell[3]-cell[6]);
	vrt[1][2]= 0.5*(cell[1]+cell[4]-cell[7]);
	vrt[2][2]= 0.5*(cell[2]+cell[5]-cell[8]);

	vrt[0][3]=-0.5*(cell[0]-cell[3]+cell[6]);
	vrt[1][3]=-0.5*(cell[1]-cell[4]+cell[7]);
	vrt[2][3]=-0.5*(cell[2]-cell[5]+cell[8]);

	vrt[0][4]=-0.5*(cell[0]+cell[3]-cell[6]);
	vrt[1][4]=-0.5*(cell[1]+cell[4]-cell[7]);
	vrt[2][4]=-0.5*(cell[2]+cell[5]-cell[8]);

	vrt[0][5]= 0.5*(cell[0]-cell[3]+cell[6]);
	vrt[1][5]= 0.5*(cell[1]-cell[4]+cell[7]);
	vrt[2][5]= 0.5*(cell[2]-cell[5]+cell[8]);

	vrt[0][6]= 0.5*(cell[0]+cell[3]+cell[6]);
	vrt[1][6]= 0.5*(cell[1]+cell[4]+cell[7]);
	vrt[2][6]= 0.5*(cell[2]+cell[5]+cell[8]);

	vrt[0][7]=-0.5*(cell[0]-cell[3]-cell[6]);
	vrt[1][7]=-0.5*(cell[1]-cell[4]-cell[7]);
	vrt[2][7]=-0.5*(cell[2]-cell[5]-cell[8]);

	edge[0][0]=0;
	edge[1][0]=1;
	edge[0][1]=1;
	edge[1][1]=2;
	edge[0][2]=2;
	edge[1][2]=3;
	edge[0][3]=3;
	edge[1][3]=0;
	edge[0][4]=4;
	edge[1][4]=5;
	edge[0][5]=5;
	edge[1][5]=6;
	edge[0][6]=6;
	edge[1][6]=7;
	edge[0][7]=7;
	edge[1][7]=4;
	edge[0][8]=0;
	edge[1][8]=4;
	edge[0][9]=1;
	edge[1][9]=5;
	edge[0][10]=2;
	edge[1][10]=6;
	edge[0][11]=3;
	edge[1][11]=7;
    }    

    void boxHEX()
    {
	/*
*********************************************************************

dl_poly/java GUI routine - hexagonal periodic boundaries

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	nvrt=12;
	nedge=18;
	double aaa=cell[0]/3.0;
	double bbb=0.5*aaa*Math.sqrt(3.0);
	double ccc=0.5*cell[8];

	vrt=new double[3][nvrt];
	edge=new int[2][nedge];

	vrt[0][0]= aaa;
	vrt[1][0]= 0.0;
	vrt[2][0]= ccc;

	vrt[0][1]= 0.5*aaa;
	vrt[1][1]= bbb;
	vrt[2][1]= ccc;

	vrt[0][2]=-0.5*aaa;
	vrt[1][2]= bbb;
	vrt[2][2]= ccc;

	vrt[0][3]=-aaa;
	vrt[1][3]= 0.0;
	vrt[2][3]= ccc;

	vrt[0][4]=-0.5*aaa;
	vrt[1][4]=-bbb;
	vrt[2][4]= ccc;

	vrt[0][5]= 0.5*aaa;
	vrt[1][5]=-bbb;
	vrt[2][5]= ccc;

	vrt[0][6]= aaa;
	vrt[1][6]= 0.0;
	vrt[2][6]=-ccc;

	vrt[0][7]= 0.5*aaa;
	vrt[1][7]= bbb;
	vrt[2][7]=-ccc;

	vrt[0][8]=-0.5*aaa;
	vrt[1][8]= bbb;
	vrt[2][8]=-ccc;

	vrt[0][9]=-aaa;
	vrt[1][9]= 0.0;
	vrt[2][9]=-ccc;

	vrt[0][10]=-0.5*aaa;
	vrt[1][10]=-bbb;
	vrt[2][10]=-ccc;

	vrt[0][11]= 0.5*aaa;
	vrt[1][11]=-bbb;
	vrt[2][11]=-ccc;

	edge[0][0]=0;
	edge[1][0]=1;
	edge[0][1]=1;
	edge[1][1]=2;
	edge[0][2]=2;
	edge[1][2]=3;
	edge[0][3]=3;
	edge[1][3]=4;
	edge[0][4]=4;
	edge[1][4]=5;
	edge[0][5]=5;
	edge[1][5]=0;
	edge[0][6]=6;
	edge[1][6]=7;
	edge[0][7]=7;
	edge[1][7]=8;
	edge[0][8]=8;
	edge[1][8]=9;
	edge[0][9]=9;
	edge[1][9]=10;
	edge[0][10]=10;
	edge[1][10]=11;
	edge[0][11]=11;
	edge[1][11]=6;
	edge[0][12]=0;
	edge[1][12]=6;
	edge[0][13]=1;
	edge[1][13]=7;
	edge[0][14]=2;
	edge[1][14]=8;
	edge[0][15]=3;
	edge[1][15]=9;
	edge[0][16]=4;
	edge[1][16]=10;
	edge[0][17]=5;
	edge[1][17]=11;
    }    

    void boxOCT()
    {
	/*
*********************************************************************

dl_poly/java GUI routine - truncated octahedral periodic boundaries

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	nvrt=24;
	nedge=36;

	vrt=new double[3][nvrt];
	edge=new int[2][nedge];

	vrt[0][0]= 0.5*cell[0];
	vrt[1][0]= 0.25*cell[0];
	vrt[2][0]= 0.0;

	vrt[0][1]= 0.5*cell[0];
	vrt[1][1]= 0.0;
	vrt[2][1]= 0.25*cell[0];

	vrt[0][2]= 0.5*cell[0];
	vrt[1][2]=-0.25*cell[0];
	vrt[2][2]= 0.0;

	vrt[0][3]= 0.5*cell[0];
	vrt[1][3]= 0.0;
	vrt[2][3]=-0.25*cell[0];

	vrt[0][4]=-0.5*cell[0];
	vrt[1][4]= 0.25*cell[0];
	vrt[2][4]= 0.0;

	vrt[0][5]=-0.5*cell[0];
	vrt[1][5]= 0.0;
	vrt[2][5]= 0.25*cell[0];

	vrt[0][6]=-0.5*cell[0];
	vrt[1][6]=-0.25*cell[0];
	vrt[2][6]= 0.0;

	vrt[0][7]=-0.5*cell[0];
	vrt[1][7]= 0.0;
	vrt[2][7]=-0.25*cell[0];

	vrt[0][8]= 0.25*cell[0];
	vrt[1][8]= 0.5*cell[0];
	vrt[2][8]= 0.0;

	vrt[0][9]= 0.0;
	vrt[1][9]= 0.5*cell[0];
	vrt[2][9]= 0.25*cell[0];

	vrt[0][10]=-0.25*cell[0];
	vrt[1][10]= 0.5*cell[0];
	vrt[2][10]= 0.0;

	vrt[0][11]= 0.0;
	vrt[1][11]= 0.5*cell[0];
	vrt[2][11]=-0.25*cell[0];

	vrt[0][12]= 0.25*cell[0];
	vrt[1][12]=-0.5*cell[0];
	vrt[2][12]= 0.0;
		   
	vrt[0][13]= 0.0;
	vrt[1][13]=-0.5*cell[0];
	vrt[2][13]= 0.25*cell[0];
		   
	vrt[0][14]=-0.25*cell[0];
	vrt[1][14]=-0.5*cell[0];
	vrt[2][14]= 0.0;
		   
	vrt[0][15]= 0.0;
	vrt[1][15]=-0.5*cell[0];
	vrt[2][15]=-0.25*cell[0];

	vrt[0][16]= 0.25*cell[0];
	vrt[1][16]= 0.0;
	vrt[2][16]= 0.5*cell[0];

	vrt[0][17]= 0.0;
	vrt[1][17]= 0.25*cell[0];
	vrt[2][17]= 0.5*cell[0];
		   
	vrt[0][18]=-0.25*cell[0];
	vrt[1][18]= 0.0;
	vrt[2][18]= 0.5*cell[0];
		   
	vrt[0][19]= 0.0;
	vrt[1][19]=-0.25*cell[0];
	vrt[2][19]= 0.5*cell[0];

	vrt[0][20]= 0.25*cell[0];
	vrt[1][20]= 0.0;
	vrt[2][20]=-0.5*cell[0];

	vrt[0][21]= 0.0;
	vrt[1][21]= 0.25*cell[0];
	vrt[2][21]=-0.5*cell[0];
		   
	vrt[0][22]=-0.25*cell[0];
	vrt[1][22]= 0.0;
	vrt[2][22]=-0.5*cell[0];
		   
	vrt[0][23]= 0.0;
	vrt[1][23]=-0.25*cell[0];
	vrt[2][23]=-0.5*cell[0];
		   
	edge[0][0]=0;
	edge[1][0]=1;
	edge[0][1]=1;
	edge[1][1]=2;
	edge[0][2]=2;
	edge[1][2]=3;
	edge[0][3]=3;
	edge[1][3]=0;
	edge[0][4]=4;
	edge[1][4]=5;
	edge[0][5]=5;
	edge[1][5]=6;
	edge[0][6]=6;
	edge[1][6]=7;
	edge[0][7]=7;
	edge[1][7]=4;
	edge[0][8]=8;
	edge[1][8]=9;
	edge[0][9]=9;
	edge[1][9]=10;
	edge[0][10]=10;
	edge[1][10]=11;
	edge[0][11]=11;
	edge[1][11]=8;
	edge[0][12]=12;
	edge[1][12]=13;
	edge[0][13]=13;
	edge[1][13]=14;
	edge[0][14]=14;
	edge[1][14]=15;
	edge[0][15]=15;
	edge[1][15]=12;
	edge[0][16]=16;
	edge[1][16]=17;
	edge[0][17]=17;
	edge[1][17]=18;
	edge[0][18]=18;
	edge[1][18]=19;
	edge[0][19]=19;
	edge[1][19]=16;
	edge[0][20]=20;
	edge[1][20]=21;
	edge[0][21]=21;
	edge[1][21]=22;
	edge[0][22]=22;
	edge[1][22]=23;
	edge[0][23]=23;
	edge[1][23]=20;
	edge[0][24]=0;
	edge[1][24]=8;
	edge[0][25]=10;
	edge[1][25]=4;
	edge[0][26]=6;
	edge[1][26]=14;
	edge[0][27]=12;
	edge[1][27]=2;
	edge[0][28]=1;
	edge[1][28]=16;
	edge[0][29]=18;
	edge[1][29]=5;
	edge[0][30]=7;
	edge[1][30]=22;
	edge[0][31]=20;
	edge[1][31]=3;
	edge[0][32]=9;
	edge[1][32]=17;
	edge[0][33]=19;
	edge[1][33]=13;
	edge[0][34]=15;
	edge[1][34]=23;
	edge[0][35]=21;
	edge[1][35]=11;

    }    

    void boxDEC()
    {
	/*
*********************************************************************

dl_poly/java GUI routine - rhombic docecahedral periodic boundaries

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	nvrt=14;
	nedge=24;
	double ddd=cell[0]/Math.sqrt(2.0);

	vrt=new double[3][nvrt];
	edge=new int[2][nedge];

	vrt[0][0]= 0.0;
	vrt[1][0]= 0.0;
	vrt[2][0]= ddd;

	vrt[0][1]= 0.5*cell[0];
	vrt[1][1]= 0.0;
	vrt[2][1]= 0.5*ddd;

	vrt[0][2]= 0.0;
	vrt[1][2]= 0.5*cell[0];
	vrt[2][2]= 0.5*ddd;

	vrt[0][3]=-0.5*cell[0];
	vrt[1][3]= 0.0;
	vrt[2][3]= 0.5*ddd;

	vrt[0][4]= 0.0;
	vrt[1][4]=-0.5*cell[0];
	vrt[2][4]= 0.5*ddd;

	vrt[0][5]= 0.5*cell[0];
	vrt[1][5]=-0.5*cell[0];
	vrt[2][5]= 0.0;

	vrt[0][6]= 0.5*cell[0];
	vrt[1][6]= 0.5*cell[0];
	vrt[2][6]= 0.0;

	vrt[0][7]=-0.5*cell[0];
	vrt[1][7]= 0.5*cell[0];
	vrt[2][7]= 0.0;

	vrt[0][8]=-0.5*cell[0];
	vrt[1][8]=-0.5*cell[0];
	vrt[2][8]= 0.0;

	vrt[0][9]= 0.5*cell[0];
	vrt[1][9]= 0.0;
	vrt[2][9]=-0.5*ddd;

	vrt[0][10]= 0.0;
	vrt[1][10]= 0.5*cell[0];
	vrt[2][10]=-0.5*ddd;

	vrt[0][11]=-0.5*cell[0];
	vrt[1][11]= 0.0;
	vrt[2][11]=-0.5*ddd;

	vrt[0][12]= 0.0;
	vrt[1][12]=-0.5*cell[0];
	vrt[2][12]=-0.5*ddd;
		   
	vrt[0][13]= 0.0;
	vrt[1][13]= 0.0;
	vrt[2][13]=-ddd;
		   
	edge[0][0]=0;
	edge[1][0]=1;
	edge[0][1]=0;
	edge[1][1]=2;
	edge[0][2]=0;
	edge[1][2]=3;
	edge[0][3]=0;
	edge[1][3]=4;
	edge[0][4]=1;
	edge[1][4]=5;
	edge[0][5]=1;
	edge[1][5]=6;
	edge[0][6]=2;
	edge[1][6]=6;
	edge[0][7]=2;
	edge[1][7]=7;
	edge[0][8]=3;
	edge[1][8]=7;
	edge[0][9]=3;
	edge[1][9]=8;
	edge[0][10]=4;
	edge[1][10]=5;
	edge[0][11]=4;
	edge[1][11]=8;
	edge[0][12]=9;
	edge[1][12]=5;
	edge[0][13]=9;
	edge[1][13]=6;
	edge[0][14]=10;
	edge[1][14]=6;
	edge[0][15]=10;
	edge[1][15]=7;
	edge[0][16]=11;
	edge[1][16]=7;
	edge[0][17]=11;
	edge[1][17]=8;
	edge[0][18]=12;
	edge[1][18]=5;
	edge[0][19]=12;
	edge[1][19]=8;
	edge[0][20]=13;
	edge[1][20]=9;
	edge[0][21]=13;
	edge[1][21]=10;
	edge[0][22]=13;
	edge[1][22]=11;
	edge[0][23]=13;
	edge[1][23]=12;
    }    

    void setBox()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/

	// select MD cell geometry

	if(imcon == 4)
	    {
		boxOCT();
	    }
	else if(imcon == 5)
	    {
		boxDEC();
	    }
	else if(imcon == 7)
	    {
		boxHEX();
	    }
	else if(imcon > 0)
	    {
		boxCELL();
	    }
    }
    
}
