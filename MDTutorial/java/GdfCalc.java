import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class GdfCalc extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to calculate van Hove distinct correlation

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    static GUI home;
    static GdfCalc job;
    static HovePlot hovplt=null;
    static double time0,time1,rcut;
    static String name1,name2;
    static int nconf,lencor,isampl,iorig,mxstr;
    static JTextField atom1,atom2,history,configs,length,sample,origin,cutoff;
    static JCheckBox format;
    static boolean form;
    static JButton run,close,status,kill,plot;

    // Define the Graphical User Interface

    public GdfCalc()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        setTitle("Gd(r,t) Calculator");
	
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

	// Instruction label 1
	
        JLabel lab1 = new JLabel("Required HISTORY file:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
	
	// Name of HISTORY file
	
        history = new JTextField(18);
        history.setBackground(scrn);
        history.setForeground(scrf);
        fix(history,grd,gbc,0,2,3,1);
        
	// History file format

        JLabel lab2 = new JLabel("File is formatted?",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);
	format=new JCheckBox("    ");
	format.setBackground(back);
	format.setForeground(fore);
	fix(format,grd,gbc,2,3,1,1);
	
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
        
	// GDF array lengths
	
        JLabel lab5 = new JLabel("GDF array lengths:",JLabel.LEFT);
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
        
	// Origin interval
	
        JLabel lab7 = new JLabel("Origin interval:",JLabel.LEFT);
        fix(lab7,grd,gbc,0,9,2,1);
        origin = new JTextField(8);
        origin.setBackground(scrn);
        origin.setForeground(scrf);
        fix(origin,grd,gbc,2,9,1,1);
        
	// Cutoff radius
	
        JLabel lab8 = new JLabel("Cutoff radius (A):",JLabel.LEFT);
        fix(lab8,grd,gbc,0,10,2,1);
        cutoff = new JTextField(8);
        cutoff.setBackground(scrn);
        cutoff.setForeground(scrf);
        fix(cutoff,grd,gbc,2,10,1,1);
        
	// Define the kill button
	
        kill = new JButton("Kill");
        kill.setBackground(butn);
        kill.setForeground(butf);
        fix(kill,grd,gbc,0,11,1,1);
	
	fix(new JLabel("  "),grd,gbc,1,11,1,1);

	// Define the Status button
	
        status = new JButton("Status");
        status.setBackground(butn);
        status.setForeground(butf);
        fix(status,grd,gbc,2,11,1,1);
	
	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,12,1,1);
	
	// Define the Plot button
	
        plot = new JButton("Plot");
        plot.setBackground(butn);
        plot.setForeground(butf);
        fix(plot,grd,gbc,2,12,1,1);
	
	// Register action buttons
	
	run.addActionListener(this);
	close.addActionListener(this);
	kill.addActionListener(this);
	status.addActionListener(this);
	plot.addActionListener(this);
        
    }

    public GdfCalc(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for calculating GDFs");
	job=new GdfCalc();
	job.pack();
	job.show();
	name1="ALL";
	name2="ALL";
	fname="HISTORY";
        nconf=1000;
	lencor=512;
	isampl=1;
	iorig=1;
	rcut=7.5;
	form=true;
	atom1.setText(name1);
        atom2.setText(name2);
	history.setText(fname);
	format.setSelected(form);
	configs.setText(String.valueOf(nconf));
	length.setText(String.valueOf(lencor));
	sample.setText(String.valueOf(isampl));
	origin.setText(String.valueOf(iorig));
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
		lencor=BML.giveInteger(length.getText(),1);
		isampl=BML.giveInteger(sample.getText(),1);
		iorig=BML.giveInteger(origin.getText(),1);
		rcut=BML.giveDouble(cutoff.getText(),1);
		try
		    {
		      if(pgdf == null)
			{
			    if(!(new File("../java/GDF/GDFPRG.X")).exists())
				{
				    monitor.println("Error - GDFPRG.X executable not available");
				}
			    else
				{
				    int frm=form?1:0;
				    config=new Config();
				    double natms=5000.0;
				    if(config.rdCFG("CONFIG"))natms=(double)config.natms;
				    mxstr=config.natms;
				    monitor.println("natms "+config.natms);
				    time0=(double)System.currentTimeMillis();
				    pgdf=Runtime.getRuntime().exec
				      ("../java/GDF/GDFPRG.X " +frm+" "+name1+" "+name2+
				       " "+fname+" "+natms+" "+nconf+" "+lencor+" "+
				       isampl+" "+iorig+" "+numgdf+" "+mxstr+" "+rcut+" &");
				    monitor.println("../java/GDF/GDFPRG.X " +frm+" "+name1+" "+name2+
				       " "+fname+" "+natms+" "+nconf+" "+lencor+" "+
				       isampl+" "+iorig+" "+numgdf+" "+mxstr+" "+rcut+" &");
				    monitor.println("GDFPRG job submitted: "+String.valueOf(pgdf));
				    numgdf+=1;
				}
			}
		      else
			{
			  monitor.println("Error - GDFPRG job already running");
			}
		    }
		catch(IOException ee)
		  {
		    monitor.println("Error GDFPRG job submission failure");
		    }
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
	else if (arg.equals("Kill") && pgdf != null)
	    {
		monitor.println("Cancelling GDFPRG job: "+String.valueOf(pgdf));
		pgdf.destroy();
		pgdf=null;
		monitor.println("GDFPRG job cancelled");
	    }
	else if (arg.equals("Status") && pgdf != null)
	    {
		try
		    {
			int state=pgdf.exitValue();
			if(state==0)
			    {
				viewFile("GDFPRG.rep");
				monitor.println("GDFPRG has terminated normally");
			    }
			else
			    monitor.println("GDFPRG has terminated abnormally");
			pgdf=null;
		    }
		catch(IllegalThreadStateException ee)
		    {
			time1=((double)System.currentTimeMillis()-time0)/1000.;

			monitor.println("GDFPRG job is running. Elapsed time (s)= "+BML.fmt(time1,9));
		    }
	    }
	else if (arg.equals("Plot"))
	    {
		if(hovplt != null)
		    hovplt.job.hide();
		hovplt=new HovePlot(home);
	    }
    }
}
