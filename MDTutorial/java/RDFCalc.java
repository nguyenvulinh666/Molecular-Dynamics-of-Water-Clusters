import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class RDFCalc extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to calculate radial distribution function

copyright - daresbury laboratory
author    - w.smith 2005

*********************************************************************
*/
    static GUI home;
    static RDFCalc job;
    static GraphDraw rdfplt=null;
    static double time0,time1,rcut,delr;
    static String name1,name2;
    static int nconf,lenrdf,mxrad,isampl;
    static JTextField atom1,atom2,history,configs,length,sample,cutoff;
    static JCheckBox format;
    static boolean form;
    static JButton run,close;
    static String[] name;
    static double[] rdf,cell,chge,weight;
    static double[][] xyz,vel,frc;

    // Define the Graphical User Interface

    public RDFCalc()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	setTitle("RDF Calculator");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	// Define the Run button
	
        run = new JButton("Run");
        run.setBackground(butn);
        run.setForeground(butf);
        fix(run,grd,gbc,0,0,1,1);
	
	fix(new JLabel("  "),grd,gbc,1,0,1,1);

	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,2,0,1,1);
	
	// Instruction label 1
	
        JLabel lab1 = new JLabel("Required HISTORY file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
	
	// Name of HISTORY file
	
        history = new JTextField(18);
        history.setBackground(scrn);
        history.setForeground(scrf);
        fix(history,grd,gbc,0,2,3,1);
        
	// History file format

	format=new JCheckBox("File");
	format.setBackground(back);
	format.setForeground(fore);
	fix(format,grd,gbc,0,3,1,1);
        JLabel lab2 = new JLabel("is formatted?",JLabel.LEFT);
        fix(lab2,grd,gbc,1,3,2,1);
	
	// Instruction label 2
	
        JLabel lab3 = new JLabel("Atomic names:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,3,1);
	
	// Name of first atom type
	
        atom1 = new JTextField(8);
        atom1.setBackground(scrn);
        atom1.setForeground(scrf);
        fix(atom1,grd,gbc,0,5,1,1);
        
	// Name of second atom type
	
        atom2 = new JTextField(8);
        atom2.setBackground(scrn);
        atom2.setForeground(scrf);
        fix(atom2,grd,gbc,2,5,1,1);
        
	// Number of configurations
	
        JLabel lab4 = new JLabel("No. configurations:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,6,2,1);
        configs = new JTextField(8);
        configs.setBackground(scrn);
        configs.setForeground(scrf);
        fix(configs,grd,gbc,2,6,1,1);
        
	// MSD array length
	
        JLabel lab5 = new JLabel("RDF array length:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,7,2,1);
        length = new JTextField(8);
        length.setBackground(scrn);
        length.setForeground(scrf);
        fix(length,grd,gbc,2,7,1,1);
        
	// Sampling interval
	
        JLabel lab6 = new JLabel("Sampling interval:",JLabel.LEFT);
        fix(lab6,grd,gbc,0,8,2,1);
        sample = new JTextField(8);
        sample.setBackground(scrn);
        sample.setForeground(scrf);
        fix(sample,grd,gbc,2,8,1,1);
        
	// Cutoff radius
	
        JLabel lab7 = new JLabel("Cutoff radius (A):",JLabel.LEFT);
        fix(lab7,grd,gbc,0,9,2,1);
        cutoff = new JTextField(8);
        cutoff.setBackground(scrn);
        cutoff.setForeground(scrf);
        fix(cutoff,grd,gbc,2,9,1,1);
        
	// Register action buttons
	
	run.addActionListener(this);
	close.addActionListener(this);
        
    }

    public RDFCalc(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for calculating RDFs");
	job=new RDFCalc();
	job.pack();
	job.show();
	name1="ALL";
	name2="ALL";
	fname="HISTORY";
        nconf=1000;
	lenrdf=512;
	isampl=1;
	rcut=7.5;
	form=true;
	atom1.setText(name1);
        atom2.setText(name2);
	history.setText(fname);
	format.setSelected(form);
	configs.setText(String.valueOf(nconf));
	length.setText(String.valueOf(lenrdf));
	sample.setText(String.valueOf(isampl));
	cutoff.setText(String.valueOf(rcut));
    }
    public void actionPerformed(ActionEvent e)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	String arg = (String)e.getActionCommand();
	if (arg.equals("Run"))
	    {
		name1=atom1.getText();
		name2=atom2.getText();
		fname=history.getText();
		form=format.isSelected();
		nconf=BML.giveInteger(configs.getText(),1);
		lenrdf=BML.giveInteger(length.getText(),1);
		isampl=BML.giveInteger(sample.getText(),1);
		rcut=BML.giveDouble(cutoff.getText(),1);
		gofr();
		rdfFile();
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
    void gofr()
    {
	/*
*********************************************************************

dl_poly/java routine to calculate radial distribution
function for selected atoms from dl_poly HISTORY file 

copyright daresbury laboratory
author  w.smith may 2005

*********************************************************************
*/      

	boolean all1,all2;
	LineNumberReader lnr=null;
	int nat1=0,nat2=0;
	int iconf,natms,k,imcon;
	double rrr,rsq,rcut2,rnrm,tx,ty,tz,ux,uy,uz,f1,f2;
	double cell[]=new double[9];
	double work[]=new double[10];
	double avcell[]=new double[9];
	double info[]=new double[10];

	rnrm=0.0;
	all1=false;
	all2=false;
	name1=name1.toUpperCase();
	name2=name2.toUpperCase();
	if(name1.equals("ALL"))all1=true;
	if(name2.equals("ALL"))all2=true;
	mxrad=Math.max(64,4*(int)lenrdf/4);
	rcut2=rcut*rcut;
	delr=rcut/mxrad;

	// write control variables

	monitor.println("Name of target HISTORY file   : "+fname);
	monitor.println("Label  of atom 1              : "+name1);
	monitor.println("Label  of atom 2              : "+name2);
	monitor.println("Length of RDF array           : "+BML.fmt(lenrdf,8));
	monitor.println("Number of configurations      : "+BML.fmt(nconf,8));
	monitor.println("Cutoff radius                 : "+BML.fmt(rcut,10));
	monitor.println("RDF bin width                 : "+BML.fmt(delr,10));

	// initialize average cell vectors

	for(int i=0;i<9;i++)
	    avcell[i]=0.0;

	// initialise RDF array

	rdf=new double[mxrad];
	
	for(int j=0;j<mxrad;j++)
	    {
		rdf[j]=0.0;
	    }
	
	// process the HISTORY file data

	for(int i=0;i<9;i++)
	    cell[i]=0.0;
	cell[0]=1.0;
	cell[4]=1.0;
	cell[8]=1.0;
	
	// initialise control parameters for HISTORY file reader
	
	info[0]=0.0;
	info[1]=999999;
	info[2]=0.0;
	info[3]=0.0;
	info[4]=0.0;
	info[5]=0.0;
	if(form)
	    {
		lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
		if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1)
		    {
			monitor.println("Error - HISTORY file data error");
			return;
		    }
	    }
	else
	    {
		monitor.println("Error - unformatted read option not active");
		return;
	    }
	natms=BML.nint(info[7]);
	
	name=new String[natms];
	chge=new double[natms];
	weight=new double[natms];
	xyz=new double[3][natms];
	vel=new double[3][natms];
	frc=new double[3][natms];

	OUT:
	for(iconf=0;iconf<nconf;iconf++)
	    {
		lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
		if(BML.nint(info[3])<0 && BML.nint(info[3])!=-1)
		    {
			monitor.println("Error - HISTORY file data error");
			info[0]=-1.0;
			lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
			return;
		    }
		if(lnr == null) break OUT;
		if(iconf==0)
		    {
			info[9]=info[6];
			for(int i=0;i<natms;i++)
			    {
				if(all1 || name1.equals(name[i]))nat1++;
				if(all2 || name2.equals(name[i]))nat2++;
			    }
		    }
		double det=invert(cell,work);
		for(int i=0;i<natms;i++)
		    {
			tx=xyz[0][i];
			ty=xyz[1][i];
			tz=xyz[2][i];
			xyz[0][i]=work[0]*tx+work[3]*ty+work[6]*tz;
			xyz[1][i]=work[1]*tx+work[4]*ty+work[7]*tz;
			xyz[2][i]=work[2]*tx+work[5]*ty+work[8]*tz;
		    }
		if(BML.nint(info[3])==-1)break OUT;
		
		imcon=BML.nint(info[5]);			
		if(imcon < 1 || imcon > 3)
		    {
			monitor.println("Error - incorrect periodic boundary condition");
			info[0]=-1.0;
			lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
			return;
		    }
		
		// running average of cell vectors
		
		f1=((double)iconf)/(iconf+1);
		f2=1.0/(iconf+1);
		for(int i=0;i<9;i++)
		    {
			avcell[i]=f1*avcell[i]+f2*cell[i];
		    }
		
		// calculate radial distribution function
		
		for(int i=1;i<natms;i++)
		    {
			if(all1 || name1.equals(name[i]))
			    {
				for(int j=0;j<i;j++)
				    {
					if(all2 || name2.equals(name[j]))
					    {
						tx=xyz[0][i]-xyz[0][j];
						ty=xyz[1][i]-xyz[1][j];
						tz=xyz[2][i]-xyz[2][j];
						tx-=BML.nint(tx);
						ty-=BML.nint(ty);
						tz-=BML.nint(tz);
						ux=tx*avcell[0]+ty*avcell[3]+tz*avcell[6];
						uy=tx*avcell[1]+ty*avcell[4]+tz*avcell[7];
						uz=tx*avcell[2]+ty*avcell[5]+tz*avcell[8];
						rsq=ux*ux+uy*uy+uz*uz;
						if(rsq < rcut2)
						    {
							rrr=Math.sqrt(rsq);
							k=(int)(rrr/delr);
							if(name[i].equals(name[j]) || name1.equals(name2))
							    {
								rdf[k]+=2.0;
							    }
							else
							    {
								rdf[k]+=1.0;
							    }
						    }
					    }
				    }
			    }
		    }
	    }
	
	if(iconf==nconf-1)
	    {
		info[0]=-1.0;
		lnr=hread(fname,name,lnr,info,cell,chge,weight,xyz,vel,frc);
	    }	
	
	if(BML.nint(info[3])==-1) iconf--;
	
	monitor.println("Number of configurations read: "+BML.fmt(iconf,8));
	
	// normalise radial distribution function
	
	dcell(avcell,work);
	rnrm=0.25*work[9]/(Math.PI*(double)(nat1*nat2)*(double)iconf);
	
	for(int i=0;i<mxrad;i++)
	    {
		rrr=(i+0.5)*delr;
		rdf[i]=rnrm*rdf[i]/(Math.pow(delr,3)*(Math.pow(i+0.5,2)+1.0/12.0));
	    }
    }
    void rdfFile()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to create a RDFDAT file

copyright - daresbury laboratory
author    - w.smith 2005

*********************************************************************
*/
	double rrr;
	try
	    {
		DataOutputStream out = new DataOutputStream(new FileOutputStream("RDFDAT"));
		out.writeBytes("Radial distribution function for atoms "+name1+" "+name2+"\n");
		out.writeBytes(BML.fmt(1,10)+BML.fmt(mxrad,10)+"\n");
		out.writeBytes(BML.fmt(name1,8)+BML.fmt(name2,8)+"\n");

		for(int i=0;i<mxrad;i++)
		    {
			rrr=delr*(i+0.5);
			out.writeBytes(BML.fmt(rrr,14)+BML.fmt(rdf[i],14)+"\n");
		    }

		out.close();
	    }
        catch(Exception e)
	    {
		monitor.println("Error file: " + fname);
	    }
	monitor.println("File RDFDAT created");
    }

}
