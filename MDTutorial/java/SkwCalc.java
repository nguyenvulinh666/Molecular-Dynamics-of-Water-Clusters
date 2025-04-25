import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class SkwCalc extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to calculate dynamic structure factors

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    static GUI home;
    static SkwCalc job;
    static SkwPlot skwplt=null;
    static double time0,time1;
    static int nconf,lencor,iorig,kmax;
    static JTextField history,configs,length,origin,kvect;
    static JCheckBox charge;
    static boolean lchg;
    static JButton run,close,status,kill,plot;

    // Define the Graphical User Interface

    public SkwCalc()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        setTitle("S(k,w) Calculator");
	
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
        
	// Number of configurations
	
        JLabel lab2 = new JLabel("No. configurations:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);
        configs = new JTextField(8);
        configs.setBackground(scrn);
        configs.setForeground(scrf);
        fix(configs,grd,gbc,2,3,1,1);
        
	// SKW array lengths
	
        JLabel lab3 = new JLabel("SKW array lengths:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,2,1);
        length = new JTextField(8);
        length.setBackground(scrn);
        length.setForeground(scrf);
        fix(length,grd,gbc,2,4,1,1);
        
	// Origin interval
	
        JLabel lab4 = new JLabel("Origin interval:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,5,2,1);
        origin = new JTextField(8);
        origin.setBackground(scrn);
        origin.setForeground(scrf);
        fix(origin,grd,gbc,2,5,1,1);
        
	// Maximum k-vector index
	
        JLabel lab5 = new JLabel("Max k-vector index:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,6,2,1);
        kvect = new JTextField(8);
        kvect.setBackground(scrn);
        kvect.setForeground(scrf);
        fix(kvect,grd,gbc,2,6,1,1);
        
	// Use charges option

        JLabel lab6 = new JLabel("Use charges?",JLabel.LEFT);
        fix(lab6,grd,gbc,0,7,2,1);
	charge=new JCheckBox("    ");
	charge.setBackground(back);
	charge.setForeground(fore);
	fix(charge,grd,gbc,2,7,1,1);
	
	// Define the kill button
	
        kill = new JButton("Kill");
        kill.setBackground(butn);
        kill.setForeground(butf);
        fix(kill,grd,gbc,0,8,1,1);
	
	fix(new JLabel("  "),grd,gbc,1,8,1,1);

	// Define the Status button
	
        status = new JButton("Status");
        status.setBackground(butn);
        status.setForeground(butf);
        fix(status,grd,gbc,2,8,1,1);
	
	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,9,1,1);
	
	// Define the Plot button
	
        plot = new JButton("Plot");
        plot.setBackground(butn);
        plot.setForeground(butf);
        fix(plot,grd,gbc,2,9,1,1);
	
	// Register action buttons
	
	run.addActionListener(this);
	close.addActionListener(this);
	kill.addActionListener(this);
	status.addActionListener(this);
	plot.addActionListener(this);
        
    }

    public SkwCalc(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for calculating S(k,w)");
	job=new SkwCalc();
	job.pack();
	job.show();
	fname="HISTORY";
        nconf=1000;
	lencor=512;
	iorig=1;
	kmax=1;
	lchg=false;
	history.setText(fname);
	charge.setSelected(lchg);
	configs.setText(String.valueOf(nconf));
	length.setText(String.valueOf(lencor));
	origin.setText(String.valueOf(iorig));
	kvect.setText(String.valueOf(kmax));
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
		fname=history.getText();
		lchg=charge.isSelected();
		nconf=BML.giveInteger(configs.getText(),1);
		lencor=BML.giveInteger(length.getText(),1);
		iorig=BML.giveInteger(origin.getText(),1);
		kmax=BML.giveInteger(kvect.getText(),1);
		try
		    {
		      if(pskw == null)
			{
			    if(!(new File("../java/SKW/SKWPRG.X")).exists())
				{
				    monitor.println("Error - SKWPRG.X executable not available");
				}
			    else
				{
				    int chg=lchg?1:0;
				    config=new Config();
				    double natms=5000.0;
				    if(config.rdCFG("CONFIG"))natms=(double)config.natms;
				    time0=(double)System.currentTimeMillis();
				    pskw=Runtime.getRuntime().exec
				      ("../java/SKW/SKWPRG.X "+fname+" 1 1 1 "+
				       chg+" "+natms+" "+nconf+" "+kmax+" "+lencor+" "+
				       iorig+" "+numskw+" &");
				    monitor.println("SKWPRG job submitted: "+String.valueOf(pskw));
				    numskw++;
				}
			}
		      else
			{
			  monitor.println("Error - SKWPRG job already running");
			}
		    }
		catch(IOException ee)
		  {
		    monitor.println("Error SKWPRG job submission failure");
		    }
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
	else if (arg.equals("Kill") && pskw != null)
	    {
		monitor.println("Cancelling SKWPRG job: "+String.valueOf(pskw));
		pskw.destroy();
		pskw=null;
		monitor.println("SKWPRG job cancelled");
	    }
	else if (arg.equals("Status") && pskw != null)
	    {
		try
		    {
			int state=pskw.exitValue();
			if(state==0)
			    {
				viewFile("SKWPRG.rep");
				monitor.println("SKWPRG has terminated normally");
			    }
			else
			    monitor.println("SKWPRG has terminated abnormally");
			pskw=null;
		    }
		catch(IllegalThreadStateException ee)
		    {
			time1=((double)System.currentTimeMillis()-time0)/1000.;

			monitor.println("SKWPRG job is running. Elapsed time (s)= "+BML.fmt(time1,9));
		    }
	    }
	else if (arg.equals("Plot"))
	    {
		if(skwplt != null)
		    skwplt.job.hide();
		skwplt=new SkwPlot(home);
	    }
    }
}
