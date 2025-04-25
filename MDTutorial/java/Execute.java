import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class Execute extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to control DL_POLY execution

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static GUI home;
    static double time0,time1;
    static Execute job;
    static JButton exec,stat,zapp,updt,dlte,ctrl,fild,cnfg,tabl,close;
    static WarningBox danger=null;

    // Define the Graphical User Interface

    public Execute()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Run DL_POLY");
	
        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;

	// Run button
	
        exec = new JButton("Run");
        exec.setBackground(butn);
        exec.setForeground(butf);
        fix(exec,grd,gbc,0,0,1,1);
	
	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,1,0,1,1);
	
	// Input file selection

	JLabel lab1=new JLabel("Select required input files:");
	fix(lab1,grd,gbc,0,1,2,1);

	// Select the CONTROL file
	
        ctrl = new JButton("CONTROL");
        ctrl.setBackground(butn);
        ctrl.setForeground(butf);
        fix(ctrl,grd,gbc,0,2,1,1);
	
	// Select the CONFIG file
	
        cnfg = new JButton("CONFIG");
        cnfg.setBackground(butn);
        cnfg.setForeground(butf);
        fix(cnfg,grd,gbc,1,2,1,1);
	
	// Select the FIELD file
	
        fild = new JButton("FIELD");
        fild.setBackground(butn);
        fild.setForeground(butf);
        fix(fild,grd,gbc,0,3,1,1);
	
	// Select the TABLE file
	
        tabl = new JButton("TABLE");
        tabl.setBackground(butn);
        tabl.setForeground(butf);
        fix(tabl,grd,gbc,1,3,1,1);
	
	// Job monitoring options
	
        JLabel lab2 = new JLabel("Job monitoring options:");
        fix(lab2,grd,gbc,0,4,2,1);
	
	// Kill job
	
        zapp = new JButton("Kill");
        zapp.setBackground(butn);
        zapp.setForeground(butf);
        fix(zapp,grd,gbc,0,5,1,1);
        
	// Job status
	
        stat = new JButton("Status");
        stat.setBackground(butn);
        stat.setForeground(butf);
        fix(stat,grd,gbc,1,5,1,1);
        
	// File handling options
	
        JLabel lab3 = new JLabel("File handling options:");
        fix(lab3,grd,gbc,0,6,2,1);
        
	// Clear data files
	
        dlte = new JButton("Clear");
        dlte.setBackground(butn);
        dlte.setForeground(butf);
        fix(dlte,grd,gbc,0,7,1,1);
        
	// Update data files
	
        updt = new JButton("Update");
        updt.setBackground(butn);
        updt.setForeground(butf);
        fix(updt,grd,gbc,1,7,1,1);

	// Register action buttons
	
	exec.addActionListener(this);
	ctrl.addActionListener(this);
	cnfg.addActionListener(this);
	fild.addActionListener(this);
	tabl.addActionListener(this);
	updt.addActionListener(this);
	zapp.addActionListener(this);
	dlte.addActionListener(this);
	stat.addActionListener(this);
	close.addActionListener(this);
        
    }

    public Execute(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	monitor.println("Activated DL_POLY Execute");
	job=new Execute();
	job.pack();
	job.show();
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
	String fname;
	String arg = (String)e.getActionCommand();
	if (arg.equals("Run"))
	    {
		try
		    {
			if(proc == null)
			    {
				if(!(new File("DLPOLY.X")).exists())
				    {
					monitor.println("Error - DLPOLY.X executable not available");
				    }
				else
				    {
					time0=(double)System.currentTimeMillis();
					proc=Runtime.getRuntime().exec("./DLPOLY.X &");
					monitor.println("DL_POLY job submitted: "+String.valueOf(proc));
				    }
			    }
			else
			    {
				monitor.println("Error - DL_POLY job already running");
			    }
		    }
		catch(IOException ee)
		    {
			monitor.println("Error DL_POLY job submission failure");
		    }
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
	else if (arg.equals("CONTROL"))
	    {
		monitor.println("Select required CONTROL file");
		if((fname=selectFileNameBegins(job,"CNT"))!=null)
		    {
			if(copyFile(fname,"CONTROL"))
			    monitor.println("Control file selected");
		    }
		else
		    monitor.println("No file selected");
	    }
	else if (arg.equals("CONFIG"))
	    {
		monitor.println("Select required CONFIG file");
		if((fname=selectFileNameBegins(job,"CFG"))!=null)
		    {
			if(copyFile(fname,"CONFIG"))
			    monitor.println("CONFIG file selected");
		    }
		else
		    monitor.println("No file selected");
	    }
	else if (arg.equals("FIELD"))
	    {
		monitor.println("Select required FIELD file");
		if((fname=selectFileNameBegins(job,"FLD"))!=null)
		    {
			if(copyFile(fname,"FIELD"))
			    monitor.println("FIELD file selected");
		    }
		else
		    monitor.println("No file selected");
	    }
	else if (arg.equals("TABLE"))
	    {
		monitor.println("Select required TABLE file");
		if((fname=selectFileNameBegins(job,"TAB"))!=null)
		    {
			if(copyFile(fname,"TABLE"))
			    monitor.println("TABLE file selected");
		    }
		else
		    monitor.println("No file selected");
	    }
	else if (arg.equals("Kill") && proc != null)
	    {
		monitor.println("Cancelling DL_POLY job: "+String.valueOf(proc));
		proc.destroy();
		proc=null;
		monitor.println("DL_POLY job cancelled");
	    }
	else if (arg.equals("Status") && proc != null)
	    {
		try
		    {
			int state=proc.exitValue();
			if(state==0)
			    monitor.println("DL_POLY has terminated normally");
			else
			    monitor.println("DL_POLY has terminated abnormally");
			proc=null;
		    }
		catch(IllegalThreadStateException ee)
		    {
			time1=((double)System.currentTimeMillis()-time0)/1000.;

			monitor.println("DL_POLY job is running. Elapsed time (s)= "+BML.fmt(time1,9));
		    }
	    }
	else if (arg.equals("Clear"))
	    {
		monitor.println("About to delete all current DL_POLY I/O files");
		danger=new WarningBox(home,"Warning!",true);
		danger.show();
		if(alert)
		    wipeOutFiles();
		else
		    monitor.println("Operation cancelled");
	    }
	else if (arg.equals("Update"))
	    {
		monitor.println("About to overwrite some current DL_POLY I/O files");
		danger=new WarningBox(home,"Warning!",true);
		danger.show();
		if(alert)
		    updateFiles();
		else
		    monitor.println("Operation cancelled");
	    }
    }
    void wipeOutFiles()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	File f=null;
	monitor.println("Deleting all current DL_POLY I/O files .....");
	if((f=new File("CONTROL")).exists())
	    {
		f.delete();
		monitor.println("File CONTROL deleted");
	    }
	if((f=new File("CONFIG")).exists())
	    {
		f.delete();
		monitor.println("File CONFIG deleted");
	    }
	if((f=new File("FIELD")).exists())
	    {
		f.delete();
		monitor.println("File FIELD deleted");
	    }
	if((f=new File("TABLE")).exists())
	    {
		f.delete();
		monitor.println("File TABLE deleted");
	    }
	if((f=new File("OUTPUT")).exists())
	    {
		f.delete();
		monitor.println("File OUTPUT deleted");
	    }
	if((f=new File("REVCON")).exists())
	    {
		f.delete();
		monitor.println("File REVCON deleted");
	    }
	if((f=new File("REVIVE")).exists())
	    {
		f.delete();
		monitor.println("File REVIVE deleted");
	    }
	if((f=new File("REVOLD")).exists())
	    {
		f.delete();
		monitor.println("File REVOLD deleted");
	    }
	if((f=new File("STATIS")).exists())
	    {
		f.delete();
		monitor.println("File STATIS deleted");
	    }
	if((f=new File("HISTORY")).exists())
	    {
		f.delete();
		monitor.println("File HISTORY deleted");
	    }
	monitor.println("All current DL_POLY I/O files deleted");
    }
    void updateFiles()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	File f=null;
	monitor.println("Updating DL_POLY I/O files .....");
	if((f=new File("CONFIG")).exists())
	    if(copyFile("CONFIG","CONFIG.BAK"))
		monitor.println("CONFIG file backup taken: CONFIG.BAK");
	else
	    monitor.println("Error CONFIG file not found");
	if((f=new File("REVCON")).exists())
	    if(copyFile("REVCON","CONFIG"))
		monitor.println("REVCON file renamed as CONFIG file");
	    else
	    monitor.println("Error REVCON file not found");
	if((f=new File("REVOLD")).exists())
	    if(copyData("REVOLD","REVOLD.BAK"))
		monitor.println("REVOLD file backup taken: REVOLD.BAK");
	else
	    monitor.println("Error REVOLD file not found");
	if((f=new File("REVIVE")).exists())
	    if(copyData("REVIVE","REVOLD"))
		monitor.println("REVIVE file renamed as REVOLD file");
	else
	    monitor.println("Error REVIVE file not found");
	monitor.println("DL_POLY I/O files updated.");
    }
}
