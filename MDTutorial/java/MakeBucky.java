import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class MakeBucky extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to construct fullerene CONFIG files

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static GUI home;
    static MakeBucky job;
    static double ccbond;
    static int natms,numx,numz,imcon;
    static double[] cell;
    static double[][] xyz;
    static Element[] atoms;
    static JTextField bond,nrx,nrz;
    static JButton make1,make2,close;

    // Define the Graphical User Interface

    public MakeBucky()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
        setTitle("Make Fullerene");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	//        Define the Make C60 button
	
        make1 = new JButton("C60");
        make1.setBackground(butn);
        make1.setForeground(butf);
        fix(make1,grd,gbc,0,0,1,1);
	
	//        Define the Make Tube button
	
        make2 = new JButton("Tube");
        make2.setBackground(butn);
        make2.setForeground(butf);
        fix(make2,grd,gbc,2,0,1,1);
        fix(new JLabel(" "),grd,gbc,1,1,1,1);
	
	//        Bond length
	
        JLabel lab1 = new JLabel("C-C Bond (A):",JLabel.LEFT);
        fix(lab1,grd,gbc,0,2,1,1);
        bond = new JTextField(8);
        bond.setBackground(scrn);
        bond.setForeground(scrf);
        fix(bond,grd,gbc,2,2,1,1);
        
	//        Tube size - rings in x direction
	
        JLabel lab2 = new JLabel("Tube size : rings X Y",JLabel.LEFT);
        fix(lab2,grd,gbc,0,3,2,1);
        nrx = new JTextField(8);
        nrx.setBackground(scrn);
        nrx.setForeground(scrf);
        fix(nrx,grd,gbc,0,4,1,1);
        
	//        Tube size - rings in y direction
	
        nrz = new JTextField(8);
        nrz.setBackground(scrn);
        nrz.setForeground(scrf);
        fix(nrz,grd,gbc,2,4,1,1);

	//        Define the Close button
	
	fix(new JLabel("  "),grd,gbc,1,5,1,1);
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,2,6,1,1);
	
	// Register action buttons
	
	make1.addActionListener(this);
	make2.addActionListener(this);
	close.addActionListener(this);
        
    }

    public MakeBucky(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for making fullerene CONFIG files");
	job=new MakeBucky();
	job.pack();
	job.show();
	ccbond=ccab;
	numx=8;
	numz=12;
	imcon=2;
        bond.setText(String.valueOf(ccbond));
        nrx.setText(String.valueOf(numx));
        nrz.setText(String.valueOf(numz));
    }

     int bucky()
    {
	/*
*********************************************************************

construction of buckminster fullerene molecule

copyright daresbury laboratory
author w. smith     november 2000

*********************************************************************
*/     
	
	boolean same;
	int call,i,j,k,n,m;
	double dr,gam,xhi,xlo,base;
	double[][] o;
	o=new double[9][2];
	cell=new double[9];
	atoms=new Element[61];
	xyz=new double[3][61];

	xhi=0.0;
	xlo=0.0;
	dr=Math.PI/180.0;
	xyz[0][0]=0.0;
	xyz[1][0]=0.850650808*ccbond;
	xyz[2][0]=2.327438436*ccbond;
	o[0][0]=Math.cos(dr*72.0);
	o[1][0]=Math.sin(dr*72.0);
	o[2][0]=0.0;
	o[3][0]=-Math.sin(dr*72.0);
	o[4][0]=Math.cos(dr*72.0);
	o[5][0]=0.0;
	o[6][0]=0.0;
	o[7][0]=0.0;
	o[8][0]=1.0;
	gam=dr*37.37736815;
	o[0][1]=Math.cos(dr*120.0);
	o[1][1]=Math.sin(dr*120.0)*Math.cos(gam);
	o[2][1]=Math.sin(dr*120.0)*Math.sin(gam);
	o[3][1]=-Math.sin(dr*120.0)*Math.cos(gam);
	o[4][1]=Math.cos(dr*120.0)*Math.pow(Math.cos(gam),2)+Math.pow(Math.sin(gam),2);
	o[5][1]=Math.cos(gam)*Math.sin(gam)*(Math.cos(dr*120.0)-1.0);
	o[6][1]=-Math.sin(dr*120.0)*Math.sin(gam);
	o[7][1]=Math.cos(gam)*Math.sin(gam)*(Math.cos(dr*120.0)-1.0);
	o[8][1]=Math.cos(dr*120.0)*Math.pow(Math.sin(gam),2)+Math.pow(Math.cos(gam),2);
	natms=1;
	atoms[0]=new Element("C_R     ");
	for(m=0;m<100;m++)
	    {
		if(natms<60)
		    {
			for(i=0;i<2;i++)
			    {
				for (j=0;j<60;j++)
				    {
					if(j < natms)
					    {
						n=natms;
						xyz[0][n]=o[0][i]*xyz[0][j]+o[3][i]*xyz[1][j]+o[6][i]*xyz[2][j];
						xyz[1][n]=o[1][i]*xyz[0][j]+o[4][i]*xyz[1][j]+o[7][i]*xyz[2][j];
						xyz[2][n]=o[2][i]*xyz[0][j]+o[5][i]*xyz[1][j]+o[8][i]*xyz[2][j];
						same=false;
						for (k=0;k<natms;k++)
						    {
							if(!same)
							    {
								same=true;
								if(Math.abs(xyz[0][n]-xyz[0][k]) > 1.e-4)
								    same=false;
								if(Math.abs(xyz[1][n]-xyz[1][k]) > 1.e-4)
								    same=false;
								if(Math.abs(xyz[2][n]-xyz[2][k]) > 1.e-4)
								    same=false;
							    }
						    }
						if(!same)
						    {
							natms++;
							atoms[n]=new Element("C_R     ");
							xhi=Math.max(xhi,xyz[0][n]);
							xlo=Math.min(xlo,xyz[0][n]);
						    }
					    }
				    }
			    }
		    }
	    }
	monitor.println("Number of atoms created : "+natms);

	// define cell vectors

	base=xhi-xlo+2.0*ccbond;
	for(i=0;i<9;i++)
	    cell[i]=0.0;
	cell[0]=base;
	cell[4]=base;
	cell[8]=base;

	// Create Config object

	config =new Config();
	config.mxatms=natms;
	config.natms=natms;
	config.imcon=imcon;
	config.atoms=atoms;
	config.cell=cell;
	config.xyz=xyz;
	config.title="Buckminster Fullerene C60";
	config.boxCELL();

	// write CONFIG file

	fname="CFGBUK."+String.valueOf(numbuk);
	if(!config.configWrite(fname)) return -1;
	monitor.println("File "+fname+" created");
	numbuk++;


	// draw structure

	pane.restore();

	return 0;
    }      
    int tube()
    {
	/*
*********************************************************************

construction of buckminster fullerene tube

copyright daresbury laboratory
author w. smith     november 2000

*********************************************************************
*/     

	int n,call,levels;
	double height,alp0,alp2,alp4,rad,xxx,yyy,zzz,ang;

	n=0;
	height=ccbond*(1.5*numz+0.5);
	levels=2*(numz+1);
	natms=numx*levels;
	cell=new double[9];
	atoms=new Element[natms];
	xyz=new double[3][natms];

	alp0=2.0*Math.PI/numx;
	alp2=alp0/2.0;
	alp4=alp0/4.0;
	rad=0.25*Math.sqrt(3.0)*ccbond/Math.sin(alp4);

	for(int i=0;i<9;i++)
	    cell[i]=0.0;
	cell[0]=4.0*rad;
	cell[4]=4.0*rad;
	cell[8]=height+ccbond;
	zzz=0.5*(ccbond-height-ccbond);
	ang=0.0;
	for(int k=0;k<levels;k++)
	    {
		for(int i=0;i<numx;i++)
		    {
			atoms[n]=new Element("C_R     ");
			xyz[0][n]=rad*Math.cos(ang+i*alp0);
			xyz[1][n]=rad*Math.sin(ang+i*alp0);
			xyz[2][n]=zzz;
			n++;
		    }
		if(k%2 == 0)
		    {
			ang+=alp2;
			zzz+=0.5*ccbond;
		    }
		else
		    {
			zzz+=ccbond;
		    }
	    }
	natms=n;
	monitor.println("Number of atoms created : "+natms);

	// Create Config object

	config =new Config();
	config.mxatms=natms;
	config.natms=natms;
	config.imcon=imcon;
	config.atoms=atoms;
	config.cell=cell;
	config.xyz=xyz;
	config.title="Buckminsterfuller tube "+BML.fmt(numx,5)+"  x"+BML.fmt(numz,5);
	config.boxCELL();

	// write CONFIG file

	fname="CFGBUK."+String.valueOf(numbuk);
	if(!config.configWrite(fname)) return -1;
	monitor.println("File "+fname+" created");
	monitor.println("radius of tube (A) is: "+BML.fmt(rad,20));
	numbuk++;

	// Draw structure

	pane.restore();

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
	if (arg.equals("C60"))
	    {
		ccbond=BML.giveDouble(bond.getText(),1);
		call=bucky();
	    }
	else if (arg.equals("Tube"))
	    {
		ccbond=BML.giveDouble(bond.getText(),1);
		numx=BML.giveInteger(nrx.getText(),1);
		numz=BML.giveInteger(nrz.getText(),1);
		call=tube();
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
}
