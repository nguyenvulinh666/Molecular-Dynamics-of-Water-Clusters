import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class ChgDefaults extends Super implements ActionListener
{
    /*
**********************************************************************

dl_poly/java routine to change the defaults for the GUI

copyright - daresbury laboratory
author    - w. smith january 2001

**********************************************************************
*/
    static GUI home;
    static ChgDefaults job;
    static JTextField rotang,tradis,bondtol;
    static JButton set,close;

    public ChgDefaults()
    {
	setTitle("Change Defaults");
	getContentPane().setForeground(fore);
	getContentPane().setBackground(back);
	setFont(fontMain);
	getContentPane().setLayout(new GridLayout(4,2));

	// define the buttons
	    
	set=new JButton("Set");
	set.setForeground(butf);
	set.setBackground(butn);
	close=new JButton("Close");
	close.setForeground(butf);
	close.setBackground(butn);
	getContentPane().add(set);
	getContentPane().add(close);

	// Fixed rotation angle

	JLabel l07=new JLabel("Rotation angle (deg)");
	getContentPane().add(l07);
	rotang=new JTextField(8);
	rotang.setBackground(scrn);
	rotang.setForeground(scrf);
	getContentPane().add(rotang);

	// Fixed translation distance

	JLabel l08=new JLabel("Translation dist (A)");
	getContentPane().add(l08);
	tradis=new JTextField(8);
	tradis.setBackground(scrn);
	tradis.setForeground(scrf);
	getContentPane().add(tradis);

	// Bond acceptance tolerance (percent)

	JLabel l09=new JLabel("Bond Tolerance (%)");
	getContentPane().add(l09);
	bondtol=new JTextField(8);
	bondtol.setBackground(scrn);
	bondtol.setForeground(scrf);
	getContentPane().add(bondtol);

	// Register action buttons
	
	set.addActionListener(this);
	close.addActionListener(this);

        
    }
    public ChgDefaults(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for changing defaults");
	job=new ChgDefaults();
	job.pack();
	job.show();
	defaultSet();
    }
    public void defaultSet()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	rotang.setText(String.valueOf(rotdef));
	tradis.setText(String.valueOf(tradef));
	bondtol.setText(String.valueOf(bondpc));
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
	if (arg.equals("Set"))
	    {
		rotdef=BML.giveDouble(rotang.getText(),1);
		tradef=BML.giveDouble(tradis.getText(),1);
		bondpc=BML.giveDouble(bondtol.getText(),1);
		home.cal=Math.cos(Math.PI*rotdef/180.0);
		home.sal=Math.sin(Math.PI*rotdef/180.0);
		home.incx=tradef;
		home.incy=tradef;
		home.incz=tradef;
		home.setPaneParams();
		monitor.println("GUI defaults now changed");
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
}

