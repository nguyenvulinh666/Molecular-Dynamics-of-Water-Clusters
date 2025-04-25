import java.awt.*;
import javax.swing.*;
import java.awt.event.*;

public class Monitor extends Super
{
	/*
*********************************************************************

dl_poly/java GUI class to open a monitor window

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    static JTextArea board;
    static Monitor job;
    static GUI home;
    static JScrollPane scroll;
    
    public Monitor()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	// Set up frame
	
	setTitle("Monitor");
	setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
        Font fontMain=new Font("Verdana",Font.PLAIN,14);
	GridLayout grd=new GridLayout(1,1);
	getContentPane().setLayout(grd);
	board = new JTextArea();
	board.setBackground(scrn);
	board.setForeground(scrf);
	board.setFont(fontMain);
	board.setEditable(false);
	board.setLineWrap(true);
	board.setWrapStyleWord(true);
	scroll = new JScrollPane(board);
	scroll.setPreferredSize(new Dimension(home.tilex,home.tiley/4));
	getContentPane().add(scroll);
	
    }

    // Monitor constructor

    public Monitor(GUI home)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	this.home=home;
	job=new Monitor();
	job.pack();
	job.show();
    }
    
    // Clear the text area

    static void clearScreen()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	board.setText("");
    }
    
    // Write comments to text area
    
    static void println(String s)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	try
	    {
		board.append(s+"\n");
		board.scrollRectToVisible(new Rectangle(0,board.getHeight()-2,1,1));
            }
	catch(Exception e)
            {
                System.out.println(s);
            }
    }
    static void print(String s)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000
v
*********************************************************************
*/
	try
	    {
		board.append(s);
            }
	catch(Exception e)
            {
                System.out.println(s);
            }
    }
}
