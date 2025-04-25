import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class ProConVar extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI to assist making CONTROL file

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static MakeControl home;
    static ProConVar job;
    static AddOptions ado;
    static JComboBox restopt;
    static JTextField tnstrun,tnsteql,tmult,tnstbpo,tnstack,tintsta,tewltol,tshktol,tqtntol;
    static JTextField tjobtim,ttclose;
    static int nstrun,nsteql,mult,nstbpo,nstack,intsta,keyres;
    static double jobtim,tclose,ewltol,shktol,qtntol;
    static JButton more,close;
    static JLabel lab1,lab2,lab3,lab4,lab5,lab6,lab7,lab8,lab9,lab10,lab11;

    // Define the Graphical User Interface

    public ProConVar()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Program Controls");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
	setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	//        Panel label
	
        lab1 = new JLabel("Program Control Variables",JLabel.LEFT);
        fix(lab1,grd,gbc,0,0,3,1);
	
	//        Number of steps
	
	tnstrun = new JTextField(8);
        tnstrun.setBackground(scrn);
        tnstrun.setForeground(scrf);
        fix(tnstrun,grd,gbc,0,1,1,1);
        lab1 = new JLabel("Number of steps",JLabel.LEFT);
        fix(lab1,grd,gbc,1,1,2,1);
        
	//        Equilibration steps
	
        tnsteql = new JTextField(8);
        tnsteql.setBackground(scrn);
        tnsteql.setForeground(scrf);
        fix(tnsteql,grd,gbc,0,2,1,1);
	lab2 = new JLabel("Equilibration steps",JLabel.LEFT);
	fix(lab2,grd,gbc,1,2,2,1);
	
	//        Multiple time step
	
        tmult = new JTextField(8);
        tmult.setBackground(scrn);
        tmult.setForeground(scrf);
        fix(tmult,grd,gbc,0,3,1,1);
	lab3 = new JLabel("Multiple time step",JLabel.LEFT);
	fix(lab3,grd,gbc,1,3,2,1);
	
	//        Print interval
	
        tnstbpo = new JTextField(8);
        tnstbpo.setBackground(scrn);
        tnstbpo.setForeground(scrf);
        fix(tnstbpo,grd,gbc,0,4,1,1);
	lab4 = new JLabel("Print interval",JLabel.LEFT);
	fix(lab4,grd,gbc,1,4,2,1);
	
	//        Stack interval
	
        tnstack = new JTextField(8);
        tnstack.setBackground(scrn);
        tnstack.setForeground(scrf);
        fix(tnstack,grd,gbc,0,5,1,1);
	lab5 = new JLabel("Stack interval",JLabel.LEFT);
	fix(lab5,grd,gbc,1,5,2,1);
	
	//        Stats interval
	
        tintsta = new JTextField(8);
        tintsta.setBackground(scrn);
        tintsta.setForeground(scrf);
        fix(tintsta,grd,gbc,0,6,1,1);
	lab5 = new JLabel("Stats interval",JLabel.LEFT);
	fix(lab5,grd,gbc,1,6,2,1);
	
	//        Ewald precision
	
        tewltol = new JTextField(8);
        tewltol.setBackground(scrn);
        tewltol.setForeground(scrf);
        fix(tewltol,grd,gbc,0,7,1,1);
	lab6 = new JLabel("Ewald precision",JLabel.LEFT);
	fix(lab6,grd,gbc,1,7,2,1);
	
	//        SHAKE tolerance
	
        tshktol = new JTextField(8);
        tshktol.setBackground(scrn);
        tshktol.setForeground(scrf);
        fix(tshktol,grd,gbc,0,8,1,1);
	lab7 = new JLabel("SHAKE tolerance",JLabel.LEFT);
	fix(lab7,grd,gbc,1,8,2,1);
	
	//        Quaternion tolerance
	
        tqtntol = new JTextField(8);
        tqtntol.setBackground(scrn);
        tqtntol.setForeground(scrf);
        fix(tqtntol,grd,gbc,0,9,1,1);
	lab8 = new JLabel("Quaternion tolerance",JLabel.LEFT);
	fix(lab8,grd,gbc,1,9,2,1);
	
	//        Job time
	
        tjobtim = new JTextField(8);
        tjobtim.setBackground(scrn);
        tjobtim.setForeground(scrf);
        fix(tjobtim,grd,gbc,0,10,1,1);
	lab9 = new JLabel("Job time (s)",JLabel.LEFT);
	fix(lab9,grd,gbc,1,10,2,1);
	
	//        Close time
	
        ttclose = new JTextField(8);
        ttclose.setBackground(scrn);
        ttclose.setForeground(scrf);
        fix(ttclose,grd,gbc,0,11,1,1);
	lab10 = new JLabel("Close time (s)",JLabel.LEFT);
	fix(lab10,grd,gbc,1,11,2,1);
	
	//        Restart option
	
	lab11 = new JLabel("Restart option",JLabel.LEFT);
	fix(lab11,grd,gbc,0,12,2,1);
        restopt = new JComboBox();
        restopt.setBackground(scrn);
        restopt.setForeground(scrf);
	restopt.addItem("NONE");
	restopt.addItem("RESTART");
	restopt.addItem("T-SCALE");
        fix(restopt,grd,gbc,2,12,1,1);
	
	//        Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,13,1,1);
	
	//        Define the More button
	
        more = new JButton("More");
        more.setBackground(butn);
        more.setForeground(butf);
        fix(more,grd,gbc,2,13,1,1);
	
	// Register action buttons
	
	close.addActionListener(this);
	more.addActionListener(this);
        
    }

    public ProConVar(MakeControl here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	job=new ProConVar();
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
	// set parameters

	mult=home.mult;
	nstbpo=home.nstbpo;
	nstrun=home.nstrun;
	nsteql=home.nsteql;
	nstack=home.nstack;
	intsta=home.intsta;
	ewltol=home.ewltol;
	shktol=home.shktol;
	qtntol=home.qtntol;
	jobtim=home.jobtim;
	tclose=home.tclose;
	keyres=home.keyres;

	// set second panel contents

        tmult.setText(String.valueOf(mult));
        tnstrun.setText(String.valueOf(nstrun));
        tnsteql.setText(String.valueOf(nsteql));
        tnstbpo.setText(String.valueOf(nstbpo));
	tnstack.setText(String.valueOf(nstack));
	tintsta.setText(String.valueOf(intsta));
	tewltol.setText(String.valueOf(ewltol));
	tshktol.setText(String.valueOf(shktol));
	tqtntol.setText(String.valueOf(qtntol));
	tjobtim.setText(String.valueOf(jobtim));
	ttclose.setText(String.valueOf(tclose));
	restopt.setSelectedIndex(keyres);
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
	else if (arg.equals("More"))
	    {
		if(ado == null)
		    ado=new AddOptions(home,job);
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
	home.ewltol=BML.giveDouble(tewltol.getText(),1);
	home.qtntol=BML.giveDouble(tqtntol.getText(),1);
	home.shktol=BML.giveDouble(tshktol.getText(),1);
	home.jobtim=BML.giveDouble(tjobtim.getText(),1);
	home.tclose=BML.giveDouble(ttclose.getText(),1);
	home.nstrun=BML.giveInteger(tnstrun.getText(),1);
	home.nsteql=BML.giveInteger(tnsteql.getText(),1);
	home.mult=BML.giveInteger(tmult.getText(),1);
	home.nstbpo=BML.giveInteger(tnstbpo.getText(),1);
	home.nstack=BML.giveInteger(tnstack.getText(),1);
	home.intsta=BML.giveInteger(tintsta.getText(),1);
	home.keyres=restopt.getSelectedIndex();
	if(ado !=null) ado.byebye();
	home.pcv=null;
	job.hide();
    }
}
