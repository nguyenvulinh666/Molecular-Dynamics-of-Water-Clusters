import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class AddOptions extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java Addoptions class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static MakeControl home;
    static ProConVar pcv;
    static AddOptions job;
    static JTextField tfcap,tnstrdf,tnstraj,tistraj,tlevcon,ttscal,tnstbts;
    static int nstrdf,nstraj,istraj,levcon,tscal,nstbts;
    static double fcap;
    static JButton close;
    static JLabel lab1,lab2,lab3,lab4,lab5,lab6,lab7,lab8,lab9,lab10,lab11;
    static JLabel lab12,lab13,lab14,lab15;
    static JCheckBox ballpairs,blcap,blzeql,blvdw,blrdf,blprdf;
    static JCheckBox bltraj,bltscal,blzden,blzero;
    static boolean allpairs,lcap,lzeql,lvdw,lrdf,lprdf;
    static boolean ltraj,ltscal,lzden,lzero;

    // Define the Graphical User Interface

    public AddOptions()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("More Options");

        getContentPane().setForeground(fore);
        getContentPane().setBackground(back);
	setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
	setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	//        Panel label
	
        JLabel lab1 = new JLabel("Select options:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,0,3,1);
	
	//        All pairs
	
	ballpairs = new JCheckBox("All pairs");
        ballpairs.setForeground(fore);
        ballpairs.setBackground(back);
        fix(ballpairs,grd,gbc,0,1,1,1);
        
	//        Cap forces
	
        blcap = new JCheckBox("Cap Forces");
        blcap.setForeground(fore);
        blcap.setBackground(back);
        fix(blcap,grd,gbc,0,2,1,1);
	
	//        Force cap
	
        tfcap = new JTextField(8);
        tfcap.setBackground(scrn);
        tfcap.setForeground(scrf);
        fix(tfcap,grd,gbc,0,3,1,1);
	JLabel lab2 = new JLabel("Force cap)",JLabel.LEFT);
	fix(lab2,grd,gbc,1,3,2,1);
	
	//        Collect stats during equilibration
	
	blzeql = new JCheckBox("Collect");
        blzeql.setForeground(fore);
        blzeql.setBackground(back);
        fix(blzeql,grd,gbc,0,4,1,1);
	
	//        Disable VDW forces
	
	blvdw = new JCheckBox("Disable VDW forces");
        blvdw.setForeground(fore);
        blvdw.setBackground(back);
        fix(blvdw,grd,gbc,0,5,2,1);
	
	//        Calculate RDF
	
	blrdf = new JCheckBox("Calculate RDF");
        blrdf.setForeground(fore);
        blrdf.setBackground(back);
        fix(blrdf,grd,gbc,0,6,1,1);
	
	//        RDF start time step
	
        tnstrdf = new JTextField(8);
        tnstrdf.setBackground(scrn);
        tnstrdf.setForeground(scrf);
        fix(tnstrdf,grd,gbc,0,7,1,1);
	JLabel lab5 = new JLabel("RDF interval",JLabel.LEFT);
	fix(lab5,grd,gbc,1,7,2,1);
	
	//        Print RDF
	
	blprdf = new JCheckBox("print RDF");
        blprdf.setForeground(fore);
        blprdf.setBackground(back);
        fix(blprdf,grd,gbc,0,8,1,1);
	
	//        Produce History file
	
	bltraj = new JCheckBox("Produce History file");
        bltraj.setForeground(fore);
        bltraj.setBackground(back);
        fix(bltraj,grd,gbc,0,9,2,1);
	
	//        History file start,increment and data level
	
	tnstraj = new JTextField(8);
        tnstraj.setBackground(scrn);
        tnstraj.setForeground(scrf);
	fix(tnstraj,grd,gbc,0,10,1,1);
	tistraj = new JTextField(8);
        tistraj.setBackground(scrn);
        tistraj.setForeground(scrf);
	fix(tistraj,grd,gbc,1,10,1,1);
	tlevcon = new JTextField(8);
        tlevcon.setBackground(scrn);
        tlevcon.setForeground(scrf);
	fix(tlevcon,grd,gbc,2,10,1,1);
	
	//        Enable Temp scaling
	
	bltscal = new JCheckBox("Enable T scaling");
        bltscal.setForeground(fore);
        bltscal.setBackground(back);
        fix(bltscal,grd,gbc,0,11,2,1);
	
	//        Temp scaling interval
	
	tnstbts = new JTextField(8);
        tnstbts.setBackground(scrn);
        tnstbts.setForeground(scrf);
	fix(tnstbts,grd,gbc,0,12,1,1);
	JLabel lab8 = new JLabel("T scaling interval",JLabel.LEFT);
	fix(lab8,grd,gbc,1,12,2,1);

	//        Z density calculation
	
	blzden = new JCheckBox("Z density");
        blzden.setForeground(fore);
        blzden.setBackground(back);
        fix(blzden,grd,gbc,0,13,1,1);

	//        Zero K MD option
	
	blzero = new JCheckBox("Zero K MD");
        blzero.setForeground(fore);
        blzero.setBackground(back);
        fix(blzero,grd,gbc,0,14,1,1);

	//        Define the Close button
	
        close = new JButton("Close");
	close.setBackground(butn);
	close.setForeground(butf);
        fix(close,grd,gbc,0,15,1,1);
	
	// Register action buttons
	
	close.addActionListener(this);
        
    }

    public AddOptions(MakeControl here,ProConVar where)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	pcv=where;
	job=new AddOptions();
	job.pack();
	job.show();
	setParams();
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

	home.ado=job;

	// set parameters values

	allpairs=home.allpairs;
	lcap=home.lcap;
	fcap=home.fcap;
	lzeql=home.lzeql;
	lvdw=home.lvdw;
	lprdf=home.lprdf;
	lrdf=home.lrdf;
	nstrdf=home.nstrdf;
	ltraj=home.ltraj;
	nstraj=home.nstraj;
	istraj=home.istraj;
	levcon=home.levcon;
	ltscal=home.ltscal;
	nstbts=home.nstbts;
	lzden=home.lzden;
	lzero=home.lzero;

	// set third panel

        ballpairs.setSelected(allpairs);
        blcap.setSelected(lcap);
        tfcap.setText(String.valueOf(fcap));
        blzeql.setSelected(lzeql);
        blvdw.setSelected(lvdw);
	blprdf.setSelected(lprdf);
        blrdf.setSelected(lrdf);
        tnstrdf.setText(String.valueOf(nstrdf));
	bltraj.setSelected(ltraj);
	tnstraj.setText(String.valueOf(nstraj));
	tistraj.setText(String.valueOf(istraj));
	tlevcon.setText(String.valueOf(levcon));
	bltscal.setSelected(ltscal);
	tnstbts.setText(String.valueOf(nstbts));
	blzden.setSelected(lzden);
	blzero.setSelected(lzero);

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
	if (arg.equals("Close"))
	    {
		byebye();
	    }
    }
    void byebye()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	//shut down panel

	home.nstrdf=BML.giveInteger(tnstrdf.getText(),1);
	home.nstbts=BML.giveInteger(tnstbts.getText(),1);
	home.nstraj=BML.giveInteger(tnstraj.getText(),1);
	home.istraj=BML.giveInteger(tistraj.getText(),1);
	home.levcon=BML.giveInteger(tlevcon.getText(),1);
	home.fcap=BML.giveDouble(tfcap.getText(),1);
	home.allpairs=ballpairs.isSelected();
	home.lcap=blcap.isSelected();
	home.lvdw=blvdw.isSelected();
	home.lzeql=blzeql.isSelected();
	home.lrdf=blrdf.isSelected();
	home.lprdf=blprdf.isSelected();
	home.ltraj=bltraj.isSelected();
	home.ltscal=bltscal.isSelected();
	home.lzden=blzden.isSelected();
	home.lzero=blzero.isSelected();
	home.ado=null;
	pcv.ado=null;
	job.hide();
    }
}
