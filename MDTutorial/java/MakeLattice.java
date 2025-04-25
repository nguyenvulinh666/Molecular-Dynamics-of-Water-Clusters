import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class MakeLattice extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to construct lattice CONFIG files 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static final int mxbas=20;
    static int nbas,nxnum,nynum,nznum,atnum;
    static int status,imcon;
    static double sxnum,synum,sznum;
    static String atnam,aname[];
    static int no[];
    static double unit[],uvw[][];
    static JButton make,clear,close,enter;
    static JTextField ax,ay,az,bx,by,bz,cx,cy,cz,sx,sy,sz;
    static JTextField nx,ny,nz,atna,atco;
    static MakeLattice job;
    static boolean kill=true;
    static GUI home=null;
    
    public MakeLattice()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Make Lattice");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
        gbc.fill=GridBagConstraints.BOTH;
        
	//        Define the Make button
	
        make = new JButton("Make");
        make.setBackground(butn);
        make.setForeground(butf);
        fix(make,grd,gbc,0,0,1,1);
	
	//        Define the Clear button
	
        clear = new JButton("Clear");
        clear.setBackground(butn);
        clear.setForeground(butf);
        fix(clear,grd,gbc,1,0,1,1);
	
	//        Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,2,0,1,1);
	
	//        Instruction label 1
	
        JLabel lab1 = new JLabel("Enter unit cell vectors:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,3,1);
	
	//        Instruction label 2
	
        JLabel lab2 = new JLabel("A vector:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,1,1);
        
	//        Unit cell component Ax
	
        ax = new JTextField(8);
        ax.setBackground(scrn);
        ax.setForeground(scrf);
        fix(ax,grd,gbc,0,3,1,1);
        
	//        Unit cell component Ay
	
        ay = new JTextField(8);
        ay.setBackground(scrn);
        ay.setForeground(scrf);
        fix(ay,grd,gbc,1,3,1,1);
        
	//        Unit cell component Az
	
        az = new JTextField(8);
        az.setBackground(scrn);
        az.setForeground(scrf);
        fix(az,grd,gbc,2,3,1,1);
        
	//        Instruction label 3
	
        JLabel lab3 = new JLabel("B vector:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,1,1);
        
	//        Unit cell component Bx
	
        bx = new JTextField(8);
        bx.setBackground(scrn);
        bx.setForeground(scrf);
        fix(bx,grd,gbc,0,5,1,1);
        
	//        Unit cell component By
	
        by = new JTextField(8);
        by.setBackground(scrn);
        by.setForeground(scrf);
        fix(by,grd,gbc,1,5,1,1);
        
	//        Unit cell component Bz
	
        bz = new JTextField(8);
        bz.setBackground(scrn);
        bz.setForeground(scrf);
        fix(bz,grd,gbc,2,5,1,1);
        
	//        Instruction label 4
	
        JLabel lab4 = new JLabel("C vector:",JLabel.LEFT);
        fix(lab4,grd,gbc,0,6,1,1);
        
	//        Unit cell component Cx
	
        cx = new JTextField(8);
        cx.setBackground(scrn);
        cx.setForeground(scrf);
        fix(cx,grd,gbc,0,7,1,1);
        
	//        Unit cell component Cy
	
        cy = new JTextField(8);
        cy.setBackground(scrn);
        cy.setForeground(scrf);
        fix(cy,grd,gbc,1,7,1,1);
        
	//        Unit cell component Cz
	
        cz = new JTextField(8);
        cz.setBackground(scrn);
        cz.setForeground(scrf);
        fix(cz,grd,gbc,2,7,1,1);
        
	//        Instruction label 5
	
        JLabel lab5 = new JLabel("Replication in A B C directions:",JLabel.LEFT);
        fix(lab5,grd,gbc,0,8,3,1);
        
	//        Replication in direction A
	
        nx = new JTextField(5);
        nx.setBackground(scrn);
        nx.setForeground(scrf);
        fix(nx,grd,gbc,0,9,1,1);
        
	//        Replication in direction B
	
        ny = new JTextField(5);
        ny.setBackground(scrn);
        ny.setForeground(scrf);
        fix(ny,grd,gbc,1,9,1,1);
        
	//        Replication in direction C
	
        nz = new JTextField(5);
        nz.setBackground(scrn);
        nz.setForeground(scrf);
        fix(nz,grd,gbc,2,9,1,1);
        
	//        Instruction label 6
	
        JLabel lab6 = new JLabel("Enter unit cell contents:",JLabel.LEFT);
        fix(lab6,grd,gbc,0,10,3,1);
        
	//        Instruction label 7
	
        JLabel lab7 = new JLabel("Atom name:",JLabel.RIGHT);
        fix(lab7,grd,gbc,0,11,2,1);
        
	//        Atom name field
	
        atna = new JTextField(8);
        atna.setBackground(scrn);
        atna.setForeground(scrf);
        fix(atna,grd,gbc,2,11,1,1);
        
	//        Instruction label 9
	
        JLabel lab9 = new JLabel("Fractional coordinates:",JLabel.LEFT);
        fix(lab9,grd,gbc,0,12,3,1);
        
	//        Fractional coordinate Sx
	
        sx = new JTextField(8);
        sx.setBackground(scrn);
        sx.setForeground(scrf);
        fix(sx,grd,gbc,0,13,1,1);
        
	//        Fractional coordinate Sy
	
        sy = new JTextField(8);
        sy.setBackground(scrn);
        sy.setForeground(scrf);
        fix(sy,grd,gbc,1,13,1,1);
        
	//        Fractional coordinate Sz
	
        sz = new JTextField(8);
        sz.setBackground(scrn);
        sz.setForeground(scrf);
        fix(sz,grd,gbc,2,13,1,1);
        
	//      Enter atomic data button
        
        enter = new JButton("Enter");
        enter.setBackground(butn);
        enter.setForeground(butf);
        fix(enter,grd,gbc,0,14,1,1);
	
	//        Instruction label 10
	
        JLabel lab10 = new JLabel("Atom count:",JLabel.RIGHT);
        fix(lab10,grd,gbc,1,14,1,1);
        
	//        Atom count field
	
        atco = new JTextField(8);
        atco.setBackground(back);
        atco.setForeground(fore);
        fix(atco,grd,gbc,2,14,1,1);
	
	// Register action buttons
	
	make.addActionListener(this);
	close.addActionListener(this);
	enter.addActionListener(this);
	clear.addActionListener(this);
	
    }

    // Constructor method
    
    public MakeLattice(GUI home)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	
	monitor.println("Activated panel for making lattice CONFIG files");
	kill=false;
	this.home=home;
	
	// Define arrays
	
	no = new int[mxbas];
	unit = new double[9];
	aname = new String[mxbas];
	uvw = new double[3][mxbas];
	
	// Set up Graphical User interface
	
	job = new MakeLattice();
	job.pack();
	job.show();
	setValues();
	
    }
    
    // Set initial values
    
    static void setValues()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        nbas=0;
        nxnum=1;
        nynum=1;
        nznum=1;
        atnum=0;
	imcon=3;
        atnam="";
        status=0;
        sxnum=0.0;
        synum=0.0;
        sznum=0.0;
        for (int i=0;i<9;i++)
	    unit[i]=0.0;
        atco.setText(String.valueOf(nbas));
        ax.setText(String.valueOf(unit[0]));
        ay.setText(String.valueOf(unit[1]));
        az.setText(String.valueOf(unit[2]));
        bx.setText(String.valueOf(unit[3]));
        by.setText(String.valueOf(unit[4]));
        bz.setText(String.valueOf(unit[5]));
        cx.setText(String.valueOf(unit[6]));
        cy.setText(String.valueOf(unit[7]));
        cz.setText(String.valueOf(unit[8]));
        nx.setText(String.valueOf(nxnum));
        ny.setText(String.valueOf(nynum));
        nz.setText(String.valueOf(nznum));
        sx.setText(String.valueOf(sxnum));
        sy.setText(String.valueOf(synum));
        sz.setText(String.valueOf(sznum));
        atna.setText(atnam);
    }
    
    // Add atom to unit cell
    
    static int addatom()
    {
	/*
**********************************************************************

dl_poly/java utility to add atoms to unit cell

author    -  w.smith october 2000
copyright - daresbury laboratory 2000

*********************************************************************
*/
        boolean safe;
        safe=true;
	
        if(nbas < mxbas)
	    {
		if(atnam.equals(""))
		    {
			println("Error - atom type unrecognised");
			safe=false;
		    }
		else
		    {
			atnum=atmnum(atnam);
			safe=true;
		    }
		if (safe)
		    {
			nbas++;
			no[nbas-1]=atnum;
			aname[nbas-1]=BML.fmt(atnam,8);
			uvw[0][nbas-1]=sxnum;
			uvw[1][nbas-1]=synum;
			uvw[2][nbas-1]=sznum;
		    }
	    }
        else
	    {
		println("Error - Maximum number of basis atoms reached");
	    }
        return nbas;
    }
    
    // Generate the super cell
    
    static int genlat()
    {
	/*
***********************************************************************

dl_poly/java routine to generate a perfect lattice with a
general unit cell

copyright - daresbury laboratory
author    - w.smith october 2000

***********************************************************************
*/
        int kkk;
        fname="CFGLAT.";
        double xx,yy,zz,xs,ys,zs,wdth;
	double cell[]=new double[9];
        double cprp[]=new double[10];
        int status=0;
        int natms=nbas*nxnum*nynum*nznum;
	Element atoms[]=new Element[natms];
	double xyz[][]=new double[3][natms];

        cell[0]=nxnum*unit[0];
        cell[1]=nxnum*unit[1];
        cell[2]=nxnum*unit[2];
        cell[3]=nynum*unit[3];
        cell[4]=nynum*unit[4];
        cell[5]=nynum*unit[5];
        cell[6]=nznum*unit[6];
        cell[7]=nznum*unit[7];
        cell[8]=nznum*unit[8];
        
	// set up CONFIG file for writing
	
        try
	    {
		fname+=String.valueOf(numlat);
		DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
		outStream.writeBytes("Lattice file generated by FILE_MAKER utility\n");
		outStream.writeBytes(BML.fmt(0,10)+BML.fmt(3,10)+"\n");
		outStream.writeBytes(BML.fmt(cell[0],20)+BML.fmt(cell[1],20)+BML.fmt(cell[2],20)+"\n");
		outStream.writeBytes(BML.fmt(cell[3],20)+BML.fmt(cell[4],20)+BML.fmt(cell[5],20)+"\n");
		outStream.writeBytes(BML.fmt(cell[6],20)+BML.fmt(cell[7],20)+BML.fmt(cell[8],20)+"\n");
		
		// set up lattice
		
		kkk=0;
		for (int k=0;k<nznum;k++)
		    {
			for (int j=0;j<nynum;j++)
			    {
				for (int i=0;i<nxnum;i++)
				    {
					for (int n=0;n<nbas;n++)
					    {
						xs=i+uvw[0][n]-0.5*nxnum;
						ys=j+uvw[1][n]-0.5*nynum;
						zs=k+uvw[2][n]-0.5*nznum;
						xx=unit[0]*xs+unit[3]*ys+unit[6]*zs;
						yy=unit[1]*xs+unit[4]*ys+unit[7]*zs;
						zz=unit[2]*xs+unit[5]*ys+unit[8]*zs;
						outStream.writeBytes(aname[n]+BML.fmt(kkk+1,10)+BML.fmt(no[n],10)+"\n");
						outStream.writeBytes(BML.fmt(xx,20)+BML.fmt(yy,20)+BML.fmt(zz,20)+"\n");
						atoms[kkk]=new Element(aname[n]);
						xyz[0][kkk]=xx;
						xyz[1][kkk]=yy;
						xyz[2][kkk]=zz;
						kkk++;
					    }
				    }
			    }
		    }
		outStream.close();
	    }
        catch(Exception e)
	    {
		println("error - writing file: "+fname);
		status=-3;
		return status;
	    }
        dcell(unit,cprp);
        wdth=0.5*BML.min(nxnum*cprp[7],nynum*cprp[8],nznum*cprp[9]);
        println("file "+fname+" created");
        println("number of atoms in system = "+BML.fmt(nbas*nxnum*nynum*nznum,10));
        println("maximum radius of cutoff  = "+BML.fmt(wdth,10));
	numlat++;

    	// Create Config object

	config =new Config();
	config.mxatms=natms;
	config.natms=natms;
	config.imcon=imcon;
	config.atoms=atoms;
	config.cell=cell;
	config.xyz=xyz;
	config.title="Lattice file generated by FILE_MAKER utility";
	config.boxCELL();

	// draw structure

	pane.restore();

        return status;
    }
    static void println(String s)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    	if(home != null)
	    monitor.println(s);
	else
	    System.out.println(s);
    }

    // Interpret the button and textfield actions
	
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
	if (arg.equals("Make"))
	    {
		unit[0]=BML.giveDouble(ax.getText(),1);
		unit[1]=BML.giveDouble(ay.getText(),1);
		unit[2]=BML.giveDouble(az.getText(),1);
		unit[3]=BML.giveDouble(bx.getText(),1);
		unit[4]=BML.giveDouble(by.getText(),1);
		unit[5]=BML.giveDouble(bz.getText(),1);
		unit[6]=BML.giveDouble(cx.getText(),1);
		unit[7]=BML.giveDouble(cy.getText(),1);
		unit[8]=BML.giveDouble(cz.getText(),1);
		nxnum=BML.giveInteger(nx.getText(),1);
		nynum=BML.giveInteger(ny.getText(),1);
		nznum=BML.giveInteger(nz.getText(),1);
		status=genlat();
	    }
	else if (arg.equals("Enter"))
	    {
		atnam=atna.getText();
		sxnum=BML.giveDouble(sx.getText(),1);
		synum=BML.giveDouble(sy.getText(),1);
		sznum=BML.giveDouble(sz.getText(),1);
		nbas=addatom();
		atnam="";
		sxnum=0.0;
		synum=0.0;
		sznum=0.0;
		atna.setText(atnam);
		atco.setText(String.valueOf(nbas));
		sx.setText(String.valueOf(sxnum));
		sy.setText(String.valueOf(synum));
		sz.setText(String.valueOf(sznum));
	    }
	else if (arg.equals("Clear"))
	    {
		setValues();
	    }
	else if (arg.equals("Close"))
	    {
		if(kill)
		    System.exit(0);
		else
		    job.hide();
	    }
    }
}
