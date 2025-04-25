import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

public class SkwPlot extends Super implements ActionListener
{
	/*
*********************************************************************

dl_poly/java GUI class to plot dynamic structure factors

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    static GUI home;
    static SkwPlot job;
    static String title;
    static int npnts,kv1,kv2,kv3;
    static JTextField kvec1,kvec2,kvec3;
    static JButton plot,close;
    static double time;
    static double[] xx,yy;

    // Define the Graphical User Interface

    public SkwPlot()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        setTitle("S(k,w) Plotter");

        getContentPane().setBackground(back);
        getContentPane().setForeground(fore);
        setFont(fontMain);
        GridBagLayout grd = new GridBagLayout();
        GridBagConstraints gbc = new GridBagConstraints();
        getContentPane().setLayout(grd);
        
	gbc.fill=GridBagConstraints.BOTH;
        
	// Define the Plot button
	
        plot = new JButton("  Plot  ");
        plot.setBackground(butn);
        plot.setForeground(butf);
        fix(plot,grd,gbc,0,0,1,1);
	
	fix(new JLabel("  "),grd,gbc,1,0,1,1);

	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,2,0,1,1);
	
	// K vector X component
	
        JLabel lab1 = new JLabel("K-vector X index:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,2,1);
        kvec1 = new JTextField(5);
        kvec1.setBackground(scrn);
        kvec1.setForeground(scrf);
        fix(kvec1,grd,gbc,2,1,1,1);

	// K vector Y component
	
        JLabel lab2 = new JLabel("K-vector Y index:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,2,1);
        kvec2 = new JTextField(5);
        kvec2.setBackground(scrn);
        kvec2.setForeground(scrf);
        fix(kvec2,grd,gbc,2,2,1,1);

	// K vector Z component
	
        JLabel lab3 = new JLabel("K-vector Z index:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,3,2,1);
        kvec3 = new JTextField(5);
        kvec3.setBackground(scrn);
        kvec3.setForeground(scrf);
        fix(kvec3,grd,gbc,2,3,1,1);
        
	// Register action buttons
	
	plot.addActionListener(this);
	close.addActionListener(this);
        
    }

    public SkwPlot(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for plotting S(k,w) functions");
	job=new SkwPlot();
	job.pack();
	job.show();
	npnts=0;
	kv1=1;
	kv2=0;
	kv3=0;
	fname="";
	kvec1.setText(String.valueOf(kv1));
	kvec2.setText(String.valueOf(kv2));
	kvec3.setText(String.valueOf(kv3));
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
	if (arg.equals("  Plot  "))
	    {
		kv1=BML.giveInteger(kvec1.getText(),1);
		kv2=BML.giveInteger(kvec2.getText(),1);
		kv3=BML.giveInteger(kvec3.getText(),1);
		fname=selectFileNameBegins(job,"DEN");
		if(fname != null)
		    {
			npnts=rdskw(fname,kv1,kv2,kv3);
			if(npnts>0)
			    {
				skwXY(npnts,kv1,kv2,kv3);
				if(graf != null)
				    graf.job.hide();
				graf=new GraphDraw(home);
				if(fname.toUpperCase().indexOf("DENFKT")>=0)
				    {
					graf.xlabel.setText("Time (ps)");
					graf.ylabel.setText("F(k,t)");
					graf.plabel.setText("F(k,t) of ("+
					BML.fmt(kv1,2)+","+BML.fmt(kv2,2)+","+
                                        BML.fmt(kv3,2)+")");
					graf.extraPlot(npnts,xx,yy);
				    }
				else if(fname.toUpperCase().indexOf("DENSKW")>=0)
				    {
					graf.xlabel.setText("Freq (/ps)");
					graf.ylabel.setText("S(k,w)");
					graf.plabel.setText("S(k,w) of ("+
					BML.fmt(kv1,2)+","+BML.fmt(kv2,2)+","+
                                        BML.fmt(kv3,2)+")");
					graf.extraPlot(npnts,xx,yy);
				    }
			    }
		    }
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
    int rdskw(String fname,int kk1,int kk2,int kk3)
    {
	/*
*********************************************************************
     
dl_poly/java routine to read a DL_POLY SKW Data file
     
copyright - daresbury laboratory
author    - w.smith march 2001

*********************************************************************
*/      
	int nskw,npts,k1,k2,k3;
	String record;
	boolean found=false;

        try
	    {
		LineNumberReader lnr = new LineNumberReader(new FileReader(fname));
		monitor.println("Reading file: "+fname);
		title = lnr.readLine();
		monitor.println("File header record: "+title);
		record = lnr.readLine();
		monitor.println("File subheader    : "+record);
		record = lnr.readLine();
		nskw=BML.giveInteger(record,1);
		npts=BML.giveInteger(record,2);
		xx=new double[npts];
		yy=new double[npts];
	    OUT:
		for(int n=0;n<nskw;n++)
		    {
			record=lnr.readLine();
			k1=BML.giveInteger(record,2);
			k2=BML.giveInteger(record,3);
			k3=BML.giveInteger(record,4);
			for(int i=0;i<npts;i++)
			    {
				record=lnr.readLine();
				xx[i]=BML.giveDouble(record,1);
				yy[i]=BML.giveDouble(record,2);
			    }

			if(k1==kk1 && k2==kk2 && k3==kk3)
			    {
				found=true;
				break OUT;
			    }
		    }
		if(!found)
		    {
			monitor.println("Error - requested function not present in file");
			lnr.close();
			return -1;
		    }
		lnr.close();
	    }
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: " + fname);
		return -2;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: " + fname + " "+e);
		return -3;
	    }
	monitor.println("Number of points ploted:"+BML.fmt(npts,6));
	return npts;
    }
    void skwXY(int npts,int kk1,int kk2,int kk3)
    {
	/*
*********************************************************************

dl_poly/java GUI routine to create a SKW XY file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	String gname;
	String[] header;
        int nhead=4,call;

	header=new String[nhead];
	gname="SKW"+String.valueOf(numskw)+".XY";
	numskw++;
	header[0]=" Skw Plotting Program";
	if(fname.toUpperCase().indexOf("DENSKW")>=0)
	    {
		header[2]="Freq (/ps)";
		header[3]="S(k,w)";
		header[1]="S(k,w) of ("+BML.fmt(kv1,2)+","+BML.fmt(kv2,2)+","+BML.fmt(kv3,2)+")";
	    }
	else
	    {
		header[2]="Time (ps)";
		header[3]="F(k,t)";
		header[1]="F(k,t) of ("+BML.fmt(kv1,2)+","+BML.fmt(kv2,2)+","+BML.fmt(kv3,2)+")";
	    }
	call=putXY(gname,header,nhead,npts,xx,yy);
	if(call==0)
	    monitor.println("PLOT file "+gname+" created");
    }
}

