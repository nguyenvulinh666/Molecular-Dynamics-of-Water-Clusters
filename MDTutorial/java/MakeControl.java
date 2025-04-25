import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeControl extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to make CONTROL files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static GUI home=null;
    static MakeControl job=null;
    static ProConVar pcv=null;
    static AddOptions ado=null;
    static String title,fname;
    static JComboBox ensemble,electro;
    static double temp,press,tstep,rcut,delr,rvdw,rprim,epsq,taut,taup,ewltol,shktol,qtntol;
    static double jobtim,tclose,fcap;
    static JTextField ttitle,ttemp,tpress,ttstep,trcut,tdelr,trvdw,trprim,tepsq,ttaut,ttaup;
    static JButton make,edit,more,close;
    static JLabel lab1,lab2,lab3,lab4,lab5,lab6,lab7,lab8,lab9,lab10,lab11,lab12,lab13,lab14;
    static int numctr,nstrun,nsteql,mult,nstbpo,nstack,intsta,keyens,keyres,nstrdf;
    static int nstraj,istraj,levcon,nstbts,keyfce;
    static boolean allpairs,lcap,lzeql,lrdf,lprdf,ltraj,ltscal,lzden,lzero,lvdw;
    static Font fontMain;

    // Define the Graphical User Interface

    public MakeControl()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Make CONTROL File");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
	setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
	setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	//        Panel label
	
        lab1 = new JLabel("System Control Variables",JLabel.LEFT);
        fix(lab1,grd,gbc,0,0,2,1);
	
	//        Define the Make button
	
        make = new JButton("Make");
        make.setBackground(butn);
        make.setForeground(butf);
        fix(make,grd,gbc,0,1,1,1);
	
	//        Define the Edit button
	
        edit = new JButton("Edit");
        edit.setBackground(butn);
        edit.setForeground(butf);
        fix(edit,grd,gbc,2,1,1,1);
	
	//        File header
	
        lab2 = new JLabel("File Header",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,1,1);
        ttitle = new JTextField(30);
        ttitle.setBackground(scrn);
        ttitle.setForeground(scrf);
        fix(ttitle,grd,gbc,0,3,3,1);
        
	//        Temperature
	
	ttemp = new JTextField(8);
        ttemp.setBackground(scrn);
        ttemp.setForeground(scrf);
        fix(ttemp,grd,gbc,0,4,1,1);
        lab3 = new JLabel("Temperature (K)",JLabel.LEFT);
        fix(lab3,grd,gbc,1,4,1,1);
        
	//        Pressure
	
        tpress = new JTextField(8);
        tpress.setBackground(scrn);
        tpress.setForeground(scrf);
        fix(tpress,grd,gbc,0,5,1,1);
	lab4 = new JLabel("Pressure (kBars)",JLabel.LEFT);
	fix(lab4,grd,gbc,1,5,1,1);
	
	//        Time step
	
        ttstep = new JTextField(8);
        ttstep.setBackground(scrn);
        ttstep.setForeground(scrf);
        fix(ttstep,grd,gbc,0,6,1,1);
	lab5 = new JLabel("Time step (ps)",JLabel.LEFT);
	fix(lab5,grd,gbc,1,6,1,1);
	
	//        Forces cut off
	
        trcut = new JTextField(8);
        trcut.setBackground(scrn);
        trcut.setForeground(scrf);
        fix(trcut,grd,gbc,0,7,1,1);
	lab6 = new JLabel("Cut off (A)",JLabel.LEFT);
	fix(lab6,grd,gbc,1,7,1,1);
	
	//        Verlet shell width
	
        tdelr = new JTextField(8);
        tdelr.setBackground(scrn);
        tdelr.setForeground(scrf);
        fix(tdelr,grd,gbc,0,8,1,1);
	lab7 = new JLabel("Verlet shell width (A)",JLabel.LEFT);
	fix(lab7,grd,gbc,1,8,1,1);
	
	//        VDW cut off
	
        trvdw = new JTextField(8);
        trvdw.setBackground(scrn);
        trvdw.setForeground(scrf);
        fix(trvdw,grd,gbc,0,9,1,1);
	lab8 = new JLabel("VDW cut off (A)",JLabel.LEFT);
	fix(lab8,grd,gbc,1,9,1,1);
	
	//        Primary cut off
	
        trprim = new JTextField(8);
        trprim.setBackground(scrn);
        trprim.setForeground(scrf);
        fix(trprim,grd,gbc,0,10,1,1);
	lab9 = new JLabel("Primary cut off (A)",JLabel.LEFT);
	fix(lab9,grd,gbc,1,10,1,1);
	
	//        Dielectric Constant
	
        tepsq = new JTextField(8);
        tepsq.setBackground(scrn);
        tepsq.setForeground(scrf);
        fix(tepsq,grd,gbc,0,11,1,1);
	lab10 = new JLabel("Dielectric constant",JLabel.LEFT);
	fix(lab10,grd,gbc,1,11,1,1);
	
	//        Temperature relaxation constant
	
        ttaut = new JTextField(8);
        ttaut.setBackground(scrn);
        ttaut.setForeground(scrf);
        fix(ttaut,grd,gbc,0,12,1,1);
	lab11 = new JLabel("Temp. relaxation (ps)",JLabel.LEFT);
	fix(lab11,grd,gbc,1,12,1,1);
	
	//        Pressure relaxation constant
	
        ttaup = new JTextField(8);
        ttaup.setBackground(scrn);
        ttaup.setForeground(scrf);
        fix(ttaup,grd,gbc,0,13,1,1);
	lab12 = new JLabel("Press. relaxation (ps)",JLabel.LEFT);
	fix(lab12,grd,gbc,1,13,1,1);
	
	//        Choice of ensemble
	
	lab13 = new JLabel("Ensemble.............",JLabel.LEFT);
	fix(lab13,grd,gbc,0,14,1,1);
        ensemble = new JComboBox();
        ensemble.setBackground(scrn);
        ensemble.setForeground(scrf);
	ensemble.addItem("NVE");
	ensemble.addItem("E-NVT");
	ensemble.addItem("B-NVT");
	ensemble.addItem("H-NVT");
	ensemble.addItem("B-NPT");
	ensemble.addItem("H-NPT");
	ensemble.addItem("B-NST");
	ensemble.addItem("H-NST");
	ensemble.addItem("PMF");
        fix(ensemble,grd,gbc,2,14,1,1);
	
	//        Choice of electrostatics
	
	lab14 = new JLabel("Electrostatics.......",JLabel.LEFT);
	fix(lab14,grd,gbc,0,15,1,1);
        electro = new JComboBox();
        electro.setBackground(scrn);
        electro.setForeground(scrf);
	electro.addItem("NONE");
	electro.addItem("EWALD");
	electro.addItem("DISTAN");
	electro.addItem("T-COUL");
	electro.addItem("S-COUL");
	electro.addItem("R-FIELD");
	electro.addItem("SPME");
	electro.addItem("HK-EWALD");
        fix(electro,grd,gbc,2,15,1,1);

	//        Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,16,1,1);
	
	//        Define the More button
	
        more = new JButton("More");
        more.setBackground(butn);
        more.setForeground(butf);
        fix(more,grd,gbc,2,16,1,1);
	
	// Register action buttons
	
	make.addActionListener(this);
	edit.addActionListener(this);
	close.addActionListener(this);
	more.addActionListener(this);
        
    }

    public MakeControl(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for making CONTROL files");
	job=new MakeControl();
	job.pack();
	job.show();
	setParams();
	setPanel();
    }
    void setParams()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	// set default values

	title="CONTROL file generated by DL_POLY/java utility";
	temp=295.0;
	press=0.0;
	tstep=0.001;
	rcut=0.0;
	delr=0.0;
	rvdw=0.0;
	rprim=0.0;
	epsq=1.0;
	taut=1.0;
	taup=1.0;
	fcap=10000.0;
	ewltol=1.E-5;
	qtntol=1.E-5;
	shktol=1.E-5;
	jobtim=100.0;
	tclose=10.0;
	nstrun=10;
	nsteql=10;
	mult=1;
	nstbts=10;
	nstbpo=10;
	nstack=100;
	intsta=10;
	keyens=0;
	keyres=0;
	keyfce=0;
	nstrdf=10;
	nstraj=100;
	istraj=10;
	levcon=0;
	allpairs=false;
	lcap=false;
	lzeql=false;
	lrdf=false;
	lprdf=false;
	ltraj=false;
	ltscal=false;
	lzden=false;
	lzero=false;
    }
    void setPanel()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	//set first panel contents

        ttitle.setText(title);
        ttemp.setText(String.valueOf(temp));
        tpress.setText(String.valueOf(press));
        ttstep.setText(String.valueOf(tstep));
        trcut.setText(String.valueOf(rcut));
	tdelr.setText(String.valueOf(delr));
	trvdw.setText(String.valueOf(rvdw));
	trprim.setText(String.valueOf(rprim));
	tepsq.setText(String.valueOf(epsq));
	ttaut.setText(String.valueOf(taut));
	ttaup.setText(String.valueOf(taup));
	ensemble.setSelectedIndex(keyens);
	electro.setSelectedIndex(keyfce/2);
    }
    void getParams()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	title=ttitle.getText();
	temp=BML.giveDouble(ttemp.getText(),1);
	press=BML.giveDouble(tpress.getText(),1);
	tstep=BML.giveDouble(ttstep.getText(),1);
	rcut=BML.giveDouble(trcut.getText(),1);
	delr=BML.giveDouble(tdelr.getText(),1);
	rvdw=BML.giveDouble(trvdw.getText(),1);
	rprim=BML.giveDouble(trprim.getText(),1);
	epsq=BML.giveDouble(tepsq.getText(),1);
	taut=BML.giveDouble(ttaut.getText(),1);
	taup=BML.giveDouble(ttaup.getText(),1);
	keyens=ensemble.getSelectedIndex();
	keyfce=2*electro.getSelectedIndex();
	if(pcv != null)
	    {
		ewltol=BML.giveDouble(pcv.tewltol.getText(),1);
		qtntol=BML.giveDouble(pcv.tqtntol.getText(),1);
		shktol=BML.giveDouble(pcv.tshktol.getText(),1);
		jobtim=BML.giveDouble(pcv.tjobtim.getText(),1);
		tclose=BML.giveDouble(pcv.ttclose.getText(),1);
		nstrun=BML.giveInteger(pcv.tnstrun.getText(),1);
		nsteql=BML.giveInteger(pcv.tnsteql.getText(),1);
		mult=BML.giveInteger(pcv.tmult.getText(),1);
		nstbpo=BML.giveInteger(pcv.tnstbpo.getText(),1);
		nstack=BML.giveInteger(pcv.tnstack.getText(),1);
		intsta=BML.giveInteger(pcv.tintsta.getText(),1);
		keyres=pcv.restopt.getSelectedIndex();
	    }
	if(ado != null)
	    {
		nstrdf=BML.giveInteger(ado.tnstrdf.getText(),1);
		nstbts=BML.giveInteger(ado.tnstbts.getText(),1);
		nstraj=BML.giveInteger(ado.tnstraj.getText(),1);
		istraj=BML.giveInteger(ado.tistraj.getText(),1);
		levcon=BML.giveInteger(ado.tlevcon.getText(),1);
		fcap=BML.giveDouble(ado.tfcap.getText(),1);
		allpairs=ado.ballpairs.isSelected();
		lcap=ado.blcap.isSelected();
		lvdw=ado.blvdw.isSelected();
		lzeql=ado.blzeql.isSelected();
		lrdf=ado.blrdf.isSelected();
		lprdf=ado.blprdf.isSelected();
		ltraj=ado.bltraj.isSelected();
		ltscal=ado.bltscal.isSelected();
		lzden=ado.blzden.isSelected();
		lzero=ado.blzero.isSelected();
	    }
    }
    int ctrMake()
    {
	/*
***********************************************************************

dl_poly/java program for constructing the dl_poly CONTROL file
defining the simulation control parameters

copyright - daresbury laboratory
author    - w. smith december 2000

***********************************************************************
*/
	boolean kill,lcut,ldelr,lprim,lrvdw;

	// intitialize system variables

	kill=false;
	lcut   = (rcut>0.0);
	ldelr  = (delr>0.0);
	lprim  = (rprim>0.0);
	lrvdw  = (rvdw>0.0);

	// check internal consistency of data

	if(lvdw)
	    keyfce = 2*(keyfce/2);
	else
	    keyfce = 2*(keyfce/2)+1;

	if(!lcut)
	    {
		kill = true;
		monitor.println("Error - no cutoff specified");
	    }

	if(!lrvdw && keyfce%2==1)
	    {
		if(lcut)
		    {
			rvdw=rcut;
		    }
		else
		    {
			kill = true;
			monitor.println("Error - no vdw cutoff set");
		    }
	    }

	if(!ldelr) 
	    {
		kill = true;
		monitor.println("Error - no verlet shell width set");
	    }

	if(mult>1) 
	    {
		if(!lprim) 
		    {
			kill = true;
			monitor.println("Error - no primary cutoff set");
		    }
	    }
	else if(rprim>rcut) 
	    {
		kill = true;
		monitor.println("Error - primary cutoff too large");
	    }

	if(keyens >= 2 && keyens <= 3) 
	    {
		if(taut <= 0.0) 
		    {
			kill=true;
			monitor.println("Error - temperature relaxation time incorrect");
		    }
	    }

	if(keyens >= 4 && keyens <= 7) 
	    {
		if(taut <= 0.0)
		    {
			kill=true;
			monitor.println("Error - temperature relaxation time incorrect");
		    }
		if(taup <= 0.0)
		    {
			kill=true;
			monitor.println("Error - pressure relaxation time incorrect");
		    }
	    }

	if(mult>1) 
	    {
		if(rcut-rprim < delr) 
		    {
			kill = true;
			monitor.println("Error - primary and secondary cutoffs incorrect");
		    }
	    }

	if(rcut < rvdw) 
	    {
		kill = true;
		monitor.println("Error - VdW cutoff exceeds general cutoff");
	    }

	if(allpairs) 
	    {
		if(mult == 1) 
		    {
			kill = true;
			monitor.println("All pairs must use multiple timestep");
		    }
		if(keyfce/2 < 2 || keyfce/2>3) 
		    {
			kill = true;
			monitor.println("Error - electrostatics incorrect for all pairs");
		    }
	    }

	if (kill) 
	    {
		monitor.println("CONTROL file NOT created - see above errors");
		return -1;
	    }

	//  open the CONTROL file for output

	fname="CNTROL."+numctr;

	try
	    {
		DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));

		outStream.writeBytes(title+"\n"+"\n");
	      
		outStream.writeBytes("temperature    "+BML.fmt(temp,8)+"\n");

		outStream.writeBytes("pressure       "+BML.fmt(press,8)+"\n");

		if(keyens ==0)
		    outStream.writeBytes("ensemble nve \n");
		else if(keyens ==1)
		    outStream.writeBytes("ensemble nvt evans \n");
		else if(keyens ==2)
		    outStream.writeBytes("ensemble nvt berendsen "+BML.fmt(taut,8)+"\n");
		else if(keyens ==3)
		    outStream.writeBytes("ensemble nvt hoover "+BML.fmt(taut,8)+"\n");
		else if(keyens ==4)
		    outStream.writeBytes("ensemble npt berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
		else if(keyens ==5)
		    outStream.writeBytes("ensemble npt hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
		else if(keyens ==6)
		    outStream.writeBytes("ensemble nst berendsen "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
		else if(keyens ==7)
		    outStream.writeBytes("ensemble nst hoover "+BML.fmt(taut,8)+BML.fmt(taup,8)+"\n");
		else if(keyens ==8)
		    outStream.writeBytes("ensemble pmf \n");
		outStream.writeBytes("\n");

		outStream.writeBytes("steps          "+BML.fmt(nstrun,8)+"\n");

		outStream.writeBytes("equilibration  "+BML.fmt(nsteql,8)+"\n");

		outStream.writeBytes("multiple step  "+BML.fmt(mult,8)+"\n");

		if(keyres ==1)
		    outStream.writeBytes("restart \n");
		else if(keyres ==2)
		    outStream.writeBytes("restart scale \n");

		if(nstbts>0)
		    outStream.writeBytes("scale          "+BML.fmt(nstbts,8)+"\n");

		if(lzeql)
		    outStream.writeBytes("collect \n");

		if(lzero)
		    outStream.writeBytes("zero temperature optimisation");

		outStream.writeBytes("print          "+BML.fmt(nstbpo,8)+"\n");

		outStream.writeBytes("stack          "+BML.fmt(nstack,8)+"\n");

		outStream.writeBytes("stats          "+BML.fmt(intsta,8)+"\n");

		if(ltraj)
		    outStream.writeBytes("trajectory     "+BML.fmt(nstraj,8)+BML.fmt(istraj,8)+BML.fmt(levcon,8)+"\n");

		if(lrdf)
		    outStream.writeBytes("rdf            "+BML.fmt(nstrdf,8)+"\n");

		outStream.writeBytes	      ("\n");

		outStream.writeBytes("timestep       "+BML.fmt(tstep,8)+"\n");

		if(rprim>0.0)
		    outStream.writeBytes("primary cutoff "+BML.fmt(rprim,8)+"\n");

		outStream.writeBytes("cutoff         "+BML.fmt(rcut,8)+"\n");

		outStream.writeBytes("delr width     "+BML.fmt(delr,8)+"\n");

		outStream.writeBytes("rvdw cutoff    "+BML.fmt(rvdw,8)+"\n");

		if(keyfce/2 ==0)
		    outStream.writeBytes("no electrostatics \n");
		else if(keyfce/2 ==1)
		    outStream.writeBytes("ewald precision"+BML.fmt(ewltol,8)+"\n");
		else if(keyfce/2 ==2) 
		    outStream.writeBytes("distan \n");
		else if(keyfce/2 ==3)
		    outStream.writeBytes("coulomb \n");
		else if(keyfce/2 ==4)
		    outStream.writeBytes("shift \n");
		else if(keyfce/2 == 5) 
		    {
			outStream.writeBytes("reaction field \n");
			outStream.writeBytes("eps constant   "+BML.fmt(epsq,8)+"\n");
		    }
		else if(keyfce/2 ==6)
		    outStream.writeBytes("spme precision"+BML.fmt(ewltol,8)+"\n");
		else if(keyfce/2 ==7)
		    outStream.writeBytes("hke precision"+BML.fmt(ewltol,8)+"\n");

		if(lcap)
		    outStream.writeBytes("cap forces     "+BML.fmt(fcap,8)+"\n");

		if(lvdw)
		    outStream.writeBytes("no vdw forces");

		outStream.writeBytes("shake tolerance"+BML.fmt(shktol,8)+"\n");

		outStream.writeBytes("quaternion tolerance"+BML.fmt(qtntol,8)+"\n");

		if(lzden)
		    outStream.writeBytes("zdensity \n");

		if(lrdf && lprdf)
		    outStream.writeBytes("print rdf \n");

		outStream.writeBytes("\n");

		outStream.writeBytes("job time            "+BML.fmt(jobtim,8)+"\n");

		outStream.writeBytes("close time          "+BML.fmt(tclose,8)+"\n");

		if(allpairs)
		    outStream.writeBytes("all pairs \n");

		outStream.writeBytes("\n");

		outStream.writeBytes("finish \n");

		// close CONTROL file

		outStream.close();
	    }
	catch(Exception e)
	    {
		monitor.println("error - writing file: "+fname);
		return -2;
	    }
	monitor.println(fname+" file created");
	numctr++;
	return 0;
    }
    int ctrLoad()
    {
	/*
***********************************************************************

dl_poly/java program for loading the dl_poly CONTROL file

copyright - daresbury laboratory
author    - w. smith december 2000

**********************************************************************
*/

	int pass=0,n=0;
        String word="",record="";
        LineNumberReader lnr;

	title="";

	// open CONTROL file

        try
	    {
		lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		title=lnr.readLine();
		monitor.println("File header record: "+title);
		pass++;

		// scan through the CONTROL file
		
		while((record=lnr.readLine()) != null)
		    {
			n=BML.countWords(record);
			if(n > 0 && record.charAt(0) != '#')
			    {
				word=BML.giveWord(record,1).toLowerCase();
				
				if(word.indexOf("steps")>=0)
				    {
					nstrun=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("equil")>=0)
				    {
					nsteql=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("restart")>=0)
				    {
					if(record.indexOf("scale")>=0)
					    keyres = 2;
					else
					    keyres = 1;
				    }
				else if(word.indexOf("ensemble")>=0)
				    {
					word=BML.giveWord(record,2).toLowerCase();
					if(word.indexOf("nve")>=0)
					    keyens=0;
					else if(word.indexOf("nvt")>=0)
					    {
						word=BML.giveWord(record,3).toLowerCase();
						if(word.indexOf("evans")>=0)
						    keyens = 1;
						else if(word.indexOf("berendsen")>=0)
						    {                 
							keyens = 2;
							taut=BML.giveDouble(record,n);
						    }
						else if(word.indexOf("hoover")>=0)
						    {
							keyens = 3;
							taut=BML.giveDouble(record,n);
						    }
					    }
					else if(word.indexOf("npt")>=0)
					    {
						word=BML.giveWord(record,3).toLowerCase();
						if(word.indexOf("berendsen")>=0)
						    {
							keyens = 4;
							taut = BML.giveDouble(record,n-1);
							taup = BML.giveDouble(record,n);
						    }
						else if(word.indexOf("hoover")>=0)
						    {
							keyens = 5;
							taut = BML.giveDouble(record,n-1);
							taup = BML.giveDouble(record,n);
						    }
					    }
					else if(word.indexOf("nst")>=0)
					    {
						word=BML.giveWord(record,3).toLowerCase();
						if(word.indexOf("berendsen")>=0)
						    {
							keyens = 6;
							taut = BML.giveDouble(record,n-1);
							taup = BML.giveDouble(record,n);
						    }
						else if(word.indexOf("hoover")>=0)
						    {
							keyens = 7;
							taut = BML.giveDouble(record,n-1);
							taup = BML.giveDouble(record,n);
						    }
					    }
					else if(word.indexOf("pmf")>=0)
					    {
						keyens = 8;
					    }
				    }
				else if(word.indexOf("scale")>=0)
				    {
					ltscal=true;
					nstbts=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("rdf")>=0)
				    {
					lrdf=true;
					nstrdf=BML.giveInteger(record,n);
				    }            
				else if(word.indexOf("zden")>=0)
				    {
					lzden=true;
				    }
				else if(word.indexOf("collect")>=0)
				    {
					lzeql=true;
				    }
				else if(word.indexOf("zero")>=0)
				    {
					lzero=true;
				    }
				else if(word.indexOf("print")>=0)
				    {
					word=BML.giveWord(record,2).toLowerCase();
					if(word.indexOf("rdf")>=0)
					    {
						    lprdf=true;
					    }
					else
					    {
						nstbpo=BML.giveInteger(record,n);
						nstbpo = Math.max(nstbpo,1);
					    }
				    }
				else if(word.indexOf("stack")>=0)
				    {
					nstack=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("stats")>=0)
				    {
					intsta=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("traj")>=0)
				    {
					ltraj=true;
					nstraj=BML.giveInteger(record,2);
					istraj=BML.giveInteger(record,3);
					levcon=BML.giveInteger(record,4);
				    }
				else if(word.indexOf("ewald")>=0)
				    {
					keyfce = 2;
					ewltol=BML.giveDouble(record,n);
			    }            
				else if(word.indexOf("distan")>=0)
				    {
					keyfce = 4;
				    }            
				else if(word.indexOf("coul")>=0)
				    {
					keyfce = 6;
				    }            
				else if(word.indexOf("shift")>=0)
				    {
					keyfce = 8;
				    }
				else if(word.indexOf("reaction")>=0)
				    {            
					keyfce = 10;
				    }            
				else if(word.indexOf("spme")>=0)
				    {            
					keyfce = 12;
					ewltol=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("hke")>=0)
				    {            
					keyfce = 14;
					ewltol=BML.giveDouble(record,n);
			    }            
				else if(word.indexOf("cap")>=0)
				    {
					lcap = true;
					fcap=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("no")>=0)
				    {
					word=BML.giveWord(record,2).toLowerCase();
					if (word.indexOf("vdw")>=0)
					    lvdw=true;
					else if(word.indexOf("elec")>=0)
					    {
						keyfce=0;
					    }
				    }
				else if(word.indexOf("mult")>=0)
				    {
					mult=BML.giveInteger(record,n);
				    }
				else if(word.indexOf("timestep")>=0)
				    {
					tstep=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("temp")>=0)
				    {
					temp=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("pres")>=0)
				    {
					press=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("prim")>=0)
				    {
					rprim=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("cut")>=0)
				    {
					rcut=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("rvdw")>=0)
				    {
					rvdw=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("delr")>=0)
				    {
					delr=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("eps")>=0)
				    {
					epsq=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("shake")>=0)
				    {
					shktol=BML.giveDouble(record,n);
				    }
				else if(word.indexOf("quat")>=0)
				    {
					qtntol=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("job")>=0)
				    {
					jobtim=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("close")>=0)
				    {
					tclose=BML.giveDouble(record,n);
				    }            
				else if(word.indexOf("all")>=0)
				    {
					allpairs = true;
				    }
				else if(word.indexOf("finish")>=0)
				    {
					monitor.println("Finished reading file "+fname);
				    }
			    }
			pass++;
		    }
		monitor.println("File "+fname+" processed lines "+pass);
		lnr.close();
	    }		    
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return -1;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return -2;
	    }
	return 0;
    }

    public void actionPerformed(ActionEvent e)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	int call;
	String arg = (String)e.getActionCommand();
	if (arg.indexOf("Make")>=0)
	    {
		getParams();
		call=ctrMake();
	    }
	else if (arg.equals("Edit"))
	    {
		if((fname=selectFileNameBegins(job,"CNT"))!=null)
		    {
			call=ctrLoad();
			setPanel();
			if(pcv != null) pcv.setParams();
			if(ado != null) ado.setParams();
		    }
		else
		    {
			monitor.println("File selection cancelled");
		    }
	    }
	else if (arg.equals("Close"))
	    {
		byebye();
	    }
	else if (arg.equals("More"))
	    {
		if(pcv == null)
		    pcv=new ProConVar(job);
	    }
    }
    static void byebye()
    {
	if(pcv != null) pcv.byebye();
	job.hide();
    }
}

