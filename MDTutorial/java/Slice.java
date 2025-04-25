import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class Slice extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to take a slice from a CONFIG/REVCON file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    static double slx,sly,slz,top,bot;
    static JButton load,make,close;
    static JTextField dlx,dly,dlz,ubd,lbd;
    static Slice job;
    static GUI home=null;
    
    public Slice()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        setTitle("Slice a REVCON file");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
        gbc.fill=GridBagConstraints.BOTH;
        
	// Define the Load button
	
        load = new JButton("Load");
        load.setBackground(butn);
        load.setForeground(butf);
        fix(load,grd,gbc,0,0,1,1);
	
	// Define the Make button
	
        make = new JButton("Make");
        make.setBackground(butn);
        make.setForeground(butf);
        fix(make,grd,gbc,2,0,1,1);
	
	// Slice direction vector
	
        JLabel lab1 = new JLabel("Slice direction vector:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
        dlx = new JTextField(8);
        dlx.setBackground(scrn);
        dlx.setForeground(scrf);
        fix(dlx,grd,gbc,0,2,1,1);
        
        dly = new JTextField(8);
        dly.setBackground(scrn);
        dly.setForeground(scrf);
        fix(dly,grd,gbc,1,2,1,1);
        
        dlz = new JTextField(8);
        dlz.setBackground(scrn);
        dlz.setForeground(scrf);
        fix(dlz,grd,gbc,2,2,1,1);
	fix(new JLabel("        "),grd,gbc,1,3,1,1);
        
	// Upper bound of slice
	
        JLabel lab2 = new JLabel("Upper bound:",JLabel.RIGHT);
        fix(lab2,grd,gbc,0,4,2,1);
        ubd = new JTextField(8);
        ubd.setBackground(scrn);
        ubd.setForeground(scrf);
        fix(ubd,grd,gbc,2,4,1,1);

	// Lower bound of slice
	
        JLabel lab3 = new JLabel("Lower bound:",JLabel.RIGHT);
        fix(lab3,grd,gbc,0,5,2,1);
        lbd = new JTextField(8);
        lbd.setBackground(scrn);
        lbd.setForeground(scrf);
        fix(lbd,grd,gbc,2,5,1,1);
        
	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,6,1,1);
	
	// Register action buttons
	
	load.addActionListener(this);
	make.addActionListener(this);
	close.addActionListener(this);
	
    }

    // Constructor method
    
    public Slice(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	monitor.println("Activated panel for slicing a REVCON file");
	home=here;
	
	// Set up Graphical User interface
	
	job = new Slice();
	job.pack();
	job.show();
	setValues();
    }
    
    // Set default values
    
    static void setValues()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        slx=0.0;
	sly=0.0;
	slz=1.0;
	top=3.0;
	bot=-3.0;
        dlx.setText(String.valueOf(slx));
        dly.setText(String.valueOf(sly));
        dlz.setText(String.valueOf(slz));
        ubd.setText(String.valueOf(top));
        lbd.setText(String.valueOf(bot));
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
	int call;
	String arg = (String)e.getActionCommand();
	if (arg.equals("Load"))
	    {
		slx=BML.giveDouble(dlx.getText(),1);
		sly=BML.giveDouble(dly.getText(),1);
		slz=BML.giveDouble(dlz.getText(),1);
		top=BML.giveDouble(ubd.getText(),1);
		bot=BML.giveDouble(lbd.getText(),1);
		config=new Config(home,"CFG");
		call=sliceFile();
	    }
	if (arg.equals("Make"))
	    {
		slx=BML.giveDouble(dlx.getText(),1);
		sly=BML.giveDouble(dly.getText(),1);
		slz=BML.giveDouble(dlz.getText(),1);
		top=BML.giveDouble(ubd.getText(),1);
		bot=BML.giveDouble(lbd.getText(),1);
		call=sliceFile();
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
    int sliceFile()
    {
	/*
**********************************************************************

dl_poly/java utility to cut a slice from a REVCON file

copyright daresbury laboratory

author w. smith march 2001

**********************************************************************
*/
	String title="";
	String aname="";
	String record="";
	int n,npt,natms,imcon;
	double ddd,sss;
	Element[] atoms;
	double[] cell;
	double[][] xyz;
		
	natms=config.natms;
	if(natms==0)return -1;
	imcon=config.imcon;
	title=config.title;
	cell=config.cell;
	xyz=config.xyz;
	atoms=config.atoms;

	npt=fname.indexOf(".");
	if(npt > 0) 
	    fname="CFGSLC"+fname.substring(npt);
	else
	    fname="CFGSLC";

	n=0;
	sss=Math.sqrt(slx*slx+sly*sly+slz*slz);
	for(int i=0;i<natms;i++)
	    {
		ddd=(slx*xyz[0][i]+sly*xyz[1][i]+slz*xyz[2][i])/sss;

		if(ddd > bot && ddd < top)
		    {
			atoms[n]=new Element(atoms[i].zsym);
			xyz[0][n]=xyz[0][i];
			xyz[1][n]=xyz[1][i];
			xyz[2][n]=xyz[2][i];
			n++;
		    }
		}
	natms=n;

	// Define new Config object

	config=new Config();
	config.title=title;
	config.mxatms=natms;
	config.natms=natms;
	config.imcon=imcon;
	config.atoms=atoms;
	config.xyz=xyz;
	config.cell=cell;
	config.boxCELL();

	// write new CONFIG file

	if(!config.configWrite(fname)) return -2;
	monitor.println("File "+fname+" created");
	monitor.println("Number of atoms in "+fname+" : "+natms);

	// Draw sliced structure

	pane.restore();

	return 0;
    }
}
