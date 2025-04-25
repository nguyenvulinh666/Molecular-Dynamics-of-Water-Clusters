import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class DataArchiver extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class for data archive methods

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static GUI home;
    static DataArchiver job;
    static JButton select,store,fetch,info,close;
    static JTextField dirold,dirnew;
    static JComboBox test;
    static String dname="DEFAULT";
    static WarningBox danger;

    // Define the Graphical User Interface

    public DataArchiver()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Data Archive");
	
        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;

	// Input file selection

	JLabel lab1=new JLabel("Select standard test case.....");
	fix(lab1,grd,gbc,0,0,3,1);

	// Select button
	
        select = new JButton("Select");
        select.setBackground(butn);
        select.setForeground(butf);
        fix(select,grd,gbc,0,1,1,1);
	
	// Test case choice
	
        test = new JComboBox();
        test.setBackground(scrn);
        test.setForeground(scrf);
	test.addItem("TEST1");
	test.addItem("TEST2");
	test.addItem("TEST3");
	test.addItem("TEST4");
	test.addItem("TEST5");
	test.addItem("TEST6");
	test.addItem("TEST7");
	test.addItem("TEST8");
	test.addItem("TEST9");
	test.addItem("TEST10");
	test.addItem("TEST11");
	test.addItem("TEST12");
	test.addItem("TEST13");
	test.addItem("TEST14");
	test.addItem("TEST15");
	test.addItem("TEST16");
	test.addItem("TEST17");
	test.addItem("TEST18");
	test.addItem("TEST19");
	test.addItem("TEST20");
        fix(test,grd,gbc,2,1,1,1);
	
	// Copy files from archive

	JLabel lab2=new JLabel("Data retrieval:");
	fix(lab2,grd,gbc,0,2,1,1);
	JLabel lab4=new JLabel("Directory:    ",JLabel.RIGHT);
	fix(lab4,grd,gbc,2,2,1,1);

	// Fetch files from archive
	
        fetch = new JButton("Fetch");
        fetch.setBackground(butn);
        fetch.setForeground(butf);
        fix(fetch,grd,gbc,0,3,1,1);
        
        dirold = new JTextField(dname);
        dirold.setBackground(scrn);
        dirold.setForeground(scrf);
        fix(dirold,grd,gbc,2,3,1,1);
	
	// Move files to archive

	JLabel lab3=new JLabel("Data storage:");
	fix(lab3,grd,gbc,0,4,1,1);
	JLabel lab5=new JLabel("Directory:    ",JLabel.RIGHT);
	fix(lab5,grd,gbc,2,4,1,1);

        store = new JButton("Store");
        store.setBackground(butn);
        store.setForeground(butf);
        fix(store,grd,gbc,0,5,1,1);
        
        dirnew = new JTextField(dname);
        dirnew.setBackground(scrn);
        dirnew.setForeground(scrf);
        fix(dirnew,grd,gbc,2,5,1,1);
	
	// pad out

	JLabel lab6=new JLabel("    ");
	fix(lab6,grd,gbc,1,6,1,1);

	// Information button

        info = new JButton("Info");
        info.setBackground(butn);
        info.setForeground(butf);
        fix(info,grd,gbc,0,7,1,1);

	// Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,2,7,1,1);
	
	// Register action buttons
	
	select.addActionListener(this);
	info.addActionListener(this);
	fetch.addActionListener(this);
	store.addActionListener(this);
	close.addActionListener(this);
        
    }

    public DataArchiver(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	monitor.println("Activated DL_POLY DataArchiver");
	job=new DataArchiver();
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

	String arg = (String)e.getActionCommand();
	if (arg.equals("Select"))
	    {
		dname=test.getSelectedItem().toString()+"/LF";
		monitor.println("About to overwrite current DL_POLY I/O files");
		danger=new WarningBox(home,"Warning!",true);
		danger.show();
		if(alert)
		    fetchFiles();
		else
		    monitor.println("Operation cancelled");
	    }
	else if (arg.equals("Fetch"))
	    {
		dname=dirold.getText();
		monitor.println("About to overwrite current DL_POLY I/O files");
		danger=new WarningBox(home,"Warning!",true);
		danger.show();
		if(alert)
		    fetchFiles();
		else
		    monitor.println("Operation cancelled");
	    }
	else if (arg.equals("Store"))
	    {
		dname=dirnew.getText();
		storeFiles();
	    }
	else if (arg.equals("Info"))
	    {
		home.viewResource("TestInfo");
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
    void fetchFiles()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	File f=null;
	String fname="";
	monitor.println("Fetching DL_POLY files from directory "+dname+".....");

	if((f=new File(fname="../data/"+dname.trim()+"/CONTROL")).exists())
	    {
		if(copyFile(fname,"CONTROL"))
		    monitor.println("File CONTROL copied");
	    }
	else
	    {
		monitor.println("CONTROL file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/CONFIG")).exists())
	    {
		if(copyFile(fname,"CONFIG"))
		    monitor.println("File CONFIG copied");
	    }
	else
	    {
		monitor.println("CONFIG file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/FIELD")).exists())
	    {
		if(copyFile(fname,"FIELD"))
		    monitor.println("File FIELD copied");
	    }
	else
	    {
		monitor.println("FIELD file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/TABLE")).exists())
	    {
		if(copyFile(fname,"TABLE"))
		    monitor.println("File TABLE copied");
	    }
	else
	    {
		monitor.println("TABLE file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/REVIVE")).exists())
	    {
		if(copyFile(fname,"REVIVE"))
		    monitor.println("File REVIVE copied");
	    }
	else
	    {
		monitor.println("REVIVE file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/REVCON")).exists())
	    {
		if(copyFile(fname,"REVCON"))
		    monitor.println("File REVCON copied");
	    }
	else
	    {
		monitor.println("REVCON file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/OUTPUT")).exists())
	    {
		if(copyFile(fname,"OUTPUT"))
		    monitor.println("File OUTPUT copied");
	    }
	else
	    {
		monitor.println("OUTPUT file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/STATIS")).exists())
	    {
		if(copyFile(fname,"STATIS"))
		    monitor.println("File STATIS copied");
	    }
	else
	    {
		monitor.println("STATIS file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/HISTORY")).exists())
	    {
		if(copyFile(fname,"HISTORY"))
		    monitor.println("File HISTORY copied");
	    }
	else
	    {
		monitor.println("HISTORY file not found");
	    }
	if((f=new File(fname="../data/"+dname.trim()+"/REVOLD")).exists())
	    {
		if(copyFile(fname,"REVOLD"))
		    monitor.println("File REVOLD copied");
	    }
	else
	    {
		monitor.println("REVOLD file not found");
	    }
	monitor.println("DL_POLY file copy completed");
    }
    void storeFiles()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	File f=null;
	String fname="";
	monitor.println("Storing DL_POLY files in directory "+dname+".....");
	
	if(!(f=new File("../data/"+dname.trim())).exists())
	    {

		(new File("../data/"+dname.trim())).mkdir();

		if((f=new File("CONTROL")).exists())
		    {
			fname="../data/"+dname.trim()+"/CONTROL";
			copyFile("CONTROL",fname);
			monitor.println("File CONTROL stored");
			f.delete();
		    }
		if((f=new File("CONFIG")).exists())
		    {
			fname="../data/"+dname.trim()+"/CONFIG";
			copyFile("CONFIG",fname);
			monitor.println("File CONFIG stored");
			f.delete();
		    }
		if((f=new File("FIELD")).exists())
		    {
			fname="../data/"+dname.trim()+"/FIELD";
			copyFile("FIELD",fname);
			monitor.println("File FIELD stored");
			f.delete();
		    }
		if((f=new File("TABLE")).exists())
		    {
			fname="../data/"+dname.trim()+"/TABLE";
			copyFile("TABLE",fname);
			monitor.println("File TABLE stored");
			f.delete();
		    }
		if((f=new File("OUTPUT")).exists())
		    {
			fname="../data/"+dname.trim()+"/OUTPUT";
			copyFile("OUTPUT",fname);
			monitor.println("File OUTPUT stored");
			f.delete();
		    }
		if((f=new File("REVIVE")).exists())
		    {
			fname="../data/"+dname.trim()+"/REVIVE";
			copyFile("REVIVE",fname);
			monitor.println("File REVIVE stored");
			f.delete();
		    }
		if((f=new File("REVOLD")).exists())
		    {
			fname="../data/"+dname.trim()+"/REVOLD";
			copyFile("REVOLD",fname);
			monitor.println("File REVOLD stored");
			f.delete();
		    }
		if((f=new File("REVCON")).exists())
		    {
			fname="../data/"+dname.trim()+"/REVCON";
			copyFile("REVCON",fname);
			monitor.println("File REVCON stored");
			f.delete();
		    }
		if((f=new File("STATIS")).exists())
		    {
			fname="../data/"+dname.trim()+"/STATIS";
			copyFile("STATIS",fname);
			monitor.println("File STATIS stored");
			f.delete();
		    }
		if((f=new File("HISTORY")).exists())
		    {
			fname="../data/"+dname.trim()+"/HISTORY";
			copyFile("HISTORY",fname);
			monitor.println("File HISTORY stored");
			f.delete();
		    }
		monitor.println("DL_POLY files stored");
	    }
	else
	    {
		monitor.println("Error - nominated storage directory already exists");
	    }
    }
}
