import java.awt.*;
import java.io.*;
import java.awt.event.*;
import java.awt.print.*;
import java.awt.geom.*;
import javax.swing.*;

// Define the DL_POLY Graphical User Interface

public class GUI extends Super
{
	/*
*********************************************************************

main dl_poly/java GUI class

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/

    private static GUI job;
    static BondLengths bondlen=null;
    static ChgDefaults newdefs=null;
    static MakeControl makctr=null;
    static ColorScheme art;
    static int tilex=600,tiley=525;
    static double cal,sal,incx,incy,incz;
    static MakeLattice maklat=null;
    static MakeBucky makbuk=null;
    static MakePoly makpol=null;
    static MakeChain makchain=null;
    static Nfold enfold=null;
    static WaterAdd addh2o=null;
    static Execute runjob=null;
    static DataArchiver datarc=null;
    static RDFPlot rdfplt=null;
    static RDFCalc rdfcal=null;
    static SokPlot sokplt=null;
    static ZdenPlot zdnplt=null;
    static Slice slcrev=null;
    static RunMSD msdrun=null;
    static RunVAF vafrun=null;
    static RunFAF fafrun=null;
    static GslCalc gslcal=null;
    static GdfCalc gdfcal=null;
    static SkwCalc skwcal=null;
    static StatProp staprp=null;
    static MakeBlankField makblank=null;
    static MakeDreiField makdrei=null;
    static MakeOPLSField makopls=null;
    static MakeTable maktable=null;
    static MakeCeramField makceram=null;
    GridBagLayout grd;
    GridBagConstraints gbc;
    JButton AC,BC,CC,DC,EC,FC,GC,HC,IC,JC;
    JMenuBar top = new JMenuBar();
    JMenu Editor = MyMenu("Editor");

    boolean edit=false;

    public GUI()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/

	// Set up Frame

	setTitle("DL_POLY GUI");

	fontMain=new Font("Verdana",Font.BOLD,14);
	setFont(fontMain);
	grd = new GridBagLayout();
	gbc = new GridBagConstraints();
	getContentPane().setLayout(grd);
	gbc.fill=GridBagConstraints.BOTH;

	// Define the graphics window

        pane=new Editor(this);
        pane.setPreferredSize(new Dimension(tilex,tiley));
	pane.setBackground(scrn);
	pane.setForeground(scrf);
        fix(pane,grd,gbc,0,0,1,14,tilex,tiley);
	setPaneParams();

	// Define the Graphical buttons

	defGraphButtons();

	// Define menu bar

	defMenuBar();

	// register object to receive and handle window events

	addWindowListener(new GUIWindowAdapter());

    }
    
    // Define the menu bar

    void defMenuBar()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to define the graphics buttons

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	//JMenuBar top = new JMenuBar();
	top.setForeground(fore);
	top.setBackground(back);
	top.setBorderPainted(false);

	// File menu A

	JMenuItem itemAA,itemAB,itemAC,itemAD,itemAE,itemAF;
	JMenu File = MyMenu("File");
	File.add(itemAA = MyMenuItem("Quit"));
	File.add(itemAB = MyMenuItem("View File"));
	File.add(itemAC = MyMenuItem("Delete file"));
	File.add(itemAD = MyMenuItem("Defaults"));
	File.add(itemAE = MyMenuItem("Reset"));
	File.add(itemAF = MyMenuItem("Print"));
	File.setFont(fontMain);
	File.setForeground(fore);
	File.setBackground(back);
	top.add(File);
	itemAA.addActionListener(new GUIMenuHandler());
	itemAB.addActionListener(new GUIMenuHandler());
	itemAC.addActionListener(new GUIMenuHandler());
	itemAD.addActionListener(new GUIMenuHandler());
	itemAE.addActionListener(new GUIMenuHandler());
	itemAF.addActionListener(new GUIMenuHandler());

	// FileMaker menu B

	JMenuItem itemBA;
	JMenuItem itemBBa,itemBBb,itemBBc,itemBBd;
	JMenuItem itemBCa,itemBCb,itemBCc,itemBCd,itemBCe;
	JMenuItem itemBDa,itemBDb,itemBDc,itemBDd,itemBDe;
	JMenuItem itemBEa,itemBEb,itemBEc;
	JMenu FileMaker  = MyMenu("FileMaker");
	FileMaker.add(itemBA = MyMenuItem("CONTROL"));
	JMenu subFileMakerBB = MyMenu("CONFIG");
	subFileMakerBB.setBackground(new Color(175,175,35));
	subFileMakerBB.add(itemBBa = MyMenuItem("Lattice"));
	subFileMakerBB.add(itemBBb = MyMenuItem("Chain"));
	subFileMakerBB.add(itemBBc = MyMenuItem("Polymer"));
	subFileMakerBB.add(itemBBd = MyMenuItem("Bucky"));
	subFileMakerBB.setForeground(Color.black);
        FileMaker.add(subFileMakerBB);
	JMenu subFileMakerBC = MyMenu("FIELD");
	subFileMakerBC.add(itemBCa = MyMenuItem("Blank"));
	subFileMakerBC.add(itemBCb = MyMenuItem("Dreiding"));
	subFileMakerBC.add(itemBCc = MyMenuItem("OPLS"));
	subFileMakerBC.add(itemBCd = MyMenuItem("Ceramics"));
	subFileMakerBC.add(itemBCe = MyMenuItem("Table"));
	subFileMakerBC.setForeground(Color.black);
        FileMaker.add(subFileMakerBC);
	JMenu subFileMakerBD = MyMenu("File Type");
	subFileMakerBD.add(itemBDa = MyMenuItem("CFG"));
	subFileMakerBD.add(itemBDb = MyMenuItem("XYZ"));
	subFileMakerBD.add(itemBDc = MyMenuItem("SEQ"));
	subFileMakerBD.add(itemBDd = MyMenuItem("MSI"));
	subFileMakerBD.add(itemBDe = MyMenuItem("MDR"));
	subFileMakerBD.setForeground(Color.black);
        FileMaker.add(subFileMakerBD);
	JMenu subFileMakerBE = MyMenu("Tools");
	subFileMakerBE.add(itemBEa = MyMenuItem("N_fold"));
	subFileMakerBE.add(itemBEb = MyMenuItem("BondLengths"));
	subFileMakerBE.add(itemBEc = MyMenuItem("Add Water"));
	subFileMakerBE.setForeground(Color.black);
        FileMaker.add(subFileMakerBE);
	FileMaker.setForeground(fore);
	FileMaker.setBackground(back);
	FileMaker.setFont(fontMain);
	top.add(FileMaker);
	itemBA.addActionListener(new GUIMenuHandler());
	itemBBa.addActionListener(new GUIMenuHandler());
	itemBBb.addActionListener(new GUIMenuHandler());
	itemBBc.addActionListener(new GUIMenuHandler());
	itemBBd.addActionListener(new GUIMenuHandler());
	itemBCa.addActionListener(new GUIMenuHandler());
	itemBCb.addActionListener(new GUIMenuHandler());
	itemBCc.addActionListener(new GUIMenuHandler());
	itemBCd.addActionListener(new GUIMenuHandler());
	itemBCe.addActionListener(new GUIMenuHandler());
	itemBDa.addActionListener(new GUIMenuHandler());
	itemBDb.addActionListener(new GUIMenuHandler());
	itemBDc.addActionListener(new GUIMenuHandler());
	itemBDd.addActionListener(new GUIMenuHandler());
	itemBDe.addActionListener(new GUIMenuHandler());
	itemBEa.addActionListener(new GUIMenuHandler());
	itemBEb.addActionListener(new GUIMenuHandler());
	itemBEc.addActionListener(new GUIMenuHandler());

	// Execute menu C

	JMenuItem itemCA,itemCB;
	JMenu Execute  = MyMenu("Execute");
	Execute.add(itemCA = MyMenuItem("Run DL_POLY"));
	Execute.add(itemCB = MyMenuItem("Store/Fetch Data"));
	Execute.setForeground(fore);
	Execute.setBackground(back);
	Execute.setFont(fontMain);
	top.add(Execute);
	itemCA.addActionListener(new GUIMenuHandler());
	itemCB.addActionListener(new GUIMenuHandler());

	// Analysis menu D

	JMenuItem itemD0;
	JMenuItem itemDAa,itemDAb,itemDAc,itemDAd,itemDAe;
	JMenuItem itemDBa,itemDBb,itemDBc;
	JMenuItem itemDCa,itemDCb,itemDCc;
	JMenuItem itemDDa,itemDDb,itemDDc;
	JMenuItem itemDEa,itemDEb;
	JMenu Analysis  = MyMenu("Analysis");
	Analysis.add(itemD0 = MyMenuItem("Statistics"));
	JMenu subAnalysisDA = MyMenu("Structure");
	subAnalysisDA.add(itemDAa = MyMenuItem("RDF_Plot"));
	subAnalysisDA.add(itemDAb = MyMenuItem("RDF_Calc"));
	subAnalysisDA.add(itemDAc = MyMenuItem("S(k)"));
	subAnalysisDA.add(itemDAd = MyMenuItem("Z_Density"));
	subAnalysisDA.add(itemDAe = MyMenuItem("Slice"));
	subAnalysisDA.setForeground(Color.black);
        Analysis.add(subAnalysisDA);
	JMenu subAnalysisDB = MyMenu("Dynamics");
	subAnalysisDB.add(itemDBa = MyMenuItem("MSD"));
	subAnalysisDB.add(itemDBb = MyMenuItem("VAF"));
	subAnalysisDB.add(itemDBc = MyMenuItem("FAF"));
	subAnalysisDB.setForeground(Color.black);
        Analysis.add(subAnalysisDB);
	JMenu subAnalysisDC = MyMenu("van Hove");
	subAnalysisDC.add(itemDCa = MyMenuItem("Gs(r,t)"));
	subAnalysisDC.add(itemDCb = MyMenuItem("Gd(r,t)"));
	subAnalysisDC.add(itemDCc = MyMenuItem("S(k,w)"));
	subAnalysisDC.setForeground(Color.black);
        Analysis.add(subAnalysisDC);
	JMenu subAnalysisDD = MyMenu("Display");
	subAnalysisDD.add(itemDDa = MyMenuItem("CONFIG"));
	subAnalysisDD.add(itemDDb = MyMenuItem("REVCON"));
	subAnalysisDD.setForeground(Color.black);
        Analysis.add(subAnalysisDD);
	JMenu subAnalysisDE = MyMenu("Tools");
	subAnalysisDE.add(itemDEa = MyMenuItem("What Atoms?"));
	subAnalysisDE.add(itemDEb = MyMenuItem("Plot"));
	subAnalysisDE.setForeground(Color.black);
        Analysis.add(subAnalysisDE);
	Analysis.setForeground(fore);
	Analysis.setBackground(back);
	Analysis.setFont(fontMain);
	top.add(Analysis);
	itemD0.addActionListener(new GUIMenuHandler());
	itemDAa.addActionListener(new GUIMenuHandler());
	itemDAb.addActionListener(new GUIMenuHandler());
	itemDAc.addActionListener(new GUIMenuHandler());
	itemDAd.addActionListener(new GUIMenuHandler());
	itemDAe.addActionListener(new GUIMenuHandler());
	itemDBa.addActionListener(new GUIMenuHandler());
	itemDBb.addActionListener(new GUIMenuHandler());
	itemDBc.addActionListener(new GUIMenuHandler());
	itemDCa.addActionListener(new GUIMenuHandler());
	itemDCb.addActionListener(new GUIMenuHandler());
	itemDCc.addActionListener(new GUIMenuHandler());
	itemDDa.addActionListener(new GUIMenuHandler());
	itemDDb.addActionListener(new GUIMenuHandler());
	itemDEa.addActionListener(new GUIMenuHandler());
	itemDEb.addActionListener(new GUIMenuHandler());

	// Information menu F

	JMenuItem itemEA,itemEB,itemEC,itemED,itemEE,itemEF;
	JMenuItem itemEG,itemEH;
	JMenu Help  = MyMenu("Information");
	Help.add(itemEA = MyMenuItem("About DL_POLY"));
	Help.add(itemEB = MyMenuItem("Disclaimer"));
	Help.add(itemEC = MyMenuItem("Licence"));
	Help.add(itemED = MyMenuItem("Acknowledgements"));
	Help.add(itemEE = MyMenuItem("MINIDREI"));
	Help.add(itemEF = MyMenuItem("MINIOPLS"));
	Help.add(itemEG = MyMenuItem("CERAMICS"));
	Help.add(itemEH = MyMenuItem("Clear Text"));
	Help.setForeground(fore);
	Help.setBackground(back);
	Help.setFont(fontMain);
	top.add(Help);
	itemEA.addActionListener(new GUIMenuHandler());
	itemEB.addActionListener(new GUIMenuHandler());
	itemEC.addActionListener(new GUIMenuHandler());
	itemED.addActionListener(new GUIMenuHandler());
	itemEE.addActionListener(new GUIMenuHandler());
	itemEF.addActionListener(new GUIMenuHandler());
	itemEG.addActionListener(new GUIMenuHandler());
	itemEH.addActionListener(new GUIMenuHandler());

	// Edit menu G

	JMenuItem itemGAa,itemGAb,itemGAc,itemGAd,itemGAe;
	JMenuItem itemGAf,itemGAg,itemGAh,itemGAi,itemGAj;
	JMenuItem itemGAk,itemGAl,itemGAm,itemGAn,itemGAo;
	JMenuItem itemGAp,itemGAq;
	JMenuItem itemGBa,itemGBb,itemGBc,itemGBd,itemGBe;
	JMenuItem itemGCa,itemGCb,itemGCc,itemGCd,itemGCe;
	JMenuItem itemGCf,itemGCg,itemGCh,itemGCi,itemGCj;
	JMenuItem itemGCk,itemGBf;
	JMenu subEditorGA = MyMenu("Atoms");
	subEditorGA.add(itemGAa = MyMenuItem("H_"));
	subEditorGA.add(itemGAb = MyMenuItem("C_3"));
	subEditorGA.add(itemGAc = MyMenuItem("C_2"));
	subEditorGA.add(itemGAd = MyMenuItem("C_1"));
	subEditorGA.add(itemGAe = MyMenuItem("C_R"));
	subEditorGA.add(itemGAf = MyMenuItem("O_3"));
	subEditorGA.add(itemGAg = MyMenuItem("O_2"));
	subEditorGA.add(itemGAh = MyMenuItem("N_3"));
	subEditorGA.add(itemGAi = MyMenuItem("N_2"));
	subEditorGA.add(itemGAj = MyMenuItem("N_1"));
	subEditorGA.add(itemGAk = MyMenuItem("P_3"));
	subEditorGA.add(itemGAl = MyMenuItem("P_2"));
	subEditorGA.add(itemGAm = MyMenuItem("S_3"));
	subEditorGA.add(itemGAn = MyMenuItem("S_2"));
	subEditorGA.add(itemGAo = MyMenuItem("OW"));
	subEditorGA.add(itemGAp = MyMenuItem("HW"));
	subEditorGA.add(itemGAq = MyMenuItem("H__HB"));
	subEditorGA.setForeground(Color.black);
	Editor.add(subEditorGA);
	JMenu subEditorGB = MyMenu("Box");
	subEditorGB.add(itemGBa = MyMenuItem("None"));
	subEditorGB.add(itemGBb = MyMenuItem("Cubic"));
	subEditorGB.add(itemGBc = MyMenuItem("O-Rhombic"));
	subEditorGB.add(itemGBd = MyMenuItem("T-Octahedral"));
	subEditorGB.add(itemGBe = MyMenuItem("R-Dodecahedral"));
	subEditorGB.add(itemGBf = MyMenuItem("Hexagonal"));
	subEditorGB.setForeground(Color.black);
	Editor.add(subEditorGB);
	JMenu subEditorGC = MyMenu("Fragment");
	subEditorGC.add(itemGCa = MyMenuItem("Search"));
	subEditorGC.add(itemGCb = MyMenuItem("Alanine"));
	subEditorGC.add(itemGCc = MyMenuItem("Benzene"));
	subEditorGC.add(itemGCd = MyMenuItem("Glucose"));
	subEditorGC.add(itemGCe = MyMenuItem("i-Butane"));
	subEditorGC.add(itemGCf = MyMenuItem("Naphthalene"));
	subEditorGC.add(itemGCg = MyMenuItem("Styrene"));
	subEditorGC.add(itemGCh = MyMenuItem("c-Hexane"));
	subEditorGC.add(itemGCi = MyMenuItem("n-Butane"));
	subEditorGC.add(itemGCj = MyMenuItem("n-Decane"));
	subEditorGC.add(itemGCk = MyMenuItem("n-Hexane"));
	subEditorGC.setForeground(Color.black);
	Editor.add(subEditorGC);
	Editor.setForeground(fore);
	Editor.setBackground(back);
	Editor.setFont(fontMain);
	itemGAa.addActionListener(new GUIMenuHandler());
	itemGAb.addActionListener(new GUIMenuHandler());
	itemGAc.addActionListener(new GUIMenuHandler());
	itemGAd.addActionListener(new GUIMenuHandler());
	itemGAe.addActionListener(new GUIMenuHandler());
	itemGAf.addActionListener(new GUIMenuHandler());
	itemGAg.addActionListener(new GUIMenuHandler());
	itemGAh.addActionListener(new GUIMenuHandler());
	itemGAi.addActionListener(new GUIMenuHandler());
	itemGAj.addActionListener(new GUIMenuHandler());
	itemGAk.addActionListener(new GUIMenuHandler());
	itemGAl.addActionListener(new GUIMenuHandler());
	itemGAm.addActionListener(new GUIMenuHandler());
	itemGAn.addActionListener(new GUIMenuHandler());
	itemGAo.addActionListener(new GUIMenuHandler());
	itemGBa.addActionListener(new GUIMenuHandler());
	itemGBb.addActionListener(new GUIMenuHandler());
	itemGBc.addActionListener(new GUIMenuHandler());
	itemGBd.addActionListener(new GUIMenuHandler());
	itemGBe.addActionListener(new GUIMenuHandler());
	itemGBf.addActionListener(new GUIMenuHandler());
	itemGCa.addActionListener(new GUIMenuHandler());
	itemGCb.addActionListener(new GUIMenuHandler());
	itemGCc.addActionListener(new GUIMenuHandler());
	itemGCd.addActionListener(new GUIMenuHandler());
	itemGCe.addActionListener(new GUIMenuHandler());
	itemGCf.addActionListener(new GUIMenuHandler());
	itemGCg.addActionListener(new GUIMenuHandler());
	itemGCh.addActionListener(new GUIMenuHandler());
	itemGCi.addActionListener(new GUIMenuHandler());
	itemGCj.addActionListener(new GUIMenuHandler());
	itemGCk.addActionListener(new GUIMenuHandler());

	// Invoke menu bar

	setJMenuBar(top);
    }

    // define the graphics buttons

    void defGraphButtons()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to define the graphics buttons

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/	
	JButton AA,AB,BA,BB,CA,CB,DA,DB,EA,EB;
	JButton FA,FB,GA,GB,HA,HB,IA,IB,JA,JB;
        AA=MyButton("New");
	fix(AA,grd,gbc,2,0,1,1,20,10);
        AB=MyButton("Clr");
	fix(AB,grd,gbc,3,0,1,1,20,10);
        BA=MyButton("Rst");
	fix(BA,grd,gbc,2,1,1,1,20,10);
        BB=MyButton("Edt");
	fix(BB,grd,gbc,3,1,1,1,20,10);

	// second block of buttons

        CA=MyButton("Tx-");
	fix(CA,grd,gbc,2,3,1,1,20,10);
        CB=MyButton("Tx+");
	fix(CB,grd,gbc,3,3,1,1,20,10);
	DA=MyButton("Ty-");
	fix(DA,grd,gbc,2,4,1,1,20,10);
        DB=MyButton("Ty+");
	fix(DB,grd,gbc,3,4,1,1,20,10);

	// third block of buttons

        EA=MyButton("Tz-");
	fix(EA,grd,gbc,2,6,1,1,20,10);
        EB=MyButton("Tz+");
	fix(EB,grd,gbc,3,6,1,1,20,10);
        FA=MyButton("Rot");
	fix(FA,grd,gbc,2,7,1,1,20,10);
        FB=MyButton("Tra");
	fix(FB,grd,gbc,3,7,1,1,20,10);

	// fourth block of buttons

        GA=MyButton("Rx-");
	fix(GA,grd,gbc,2,9,1,1,20,10);
        GB=MyButton("Rx+");
	fix(GB,grd,gbc,3,9,1,1,20,10);
        HA=MyButton("Ry-");
	fix(HA,grd,gbc,2,10,1,1,20,10);
        HB=MyButton("Ry+");
	fix(HB,grd,gbc,3,10,1,1,20,10);

	// fifth block of buttons

        IA=MyButton("Rz-");
	fix(IA,grd,gbc,2,12,1,1,20,10);
        IB=MyButton("Rz+");
	fix(IB,grd,gbc,3,12,1,1,20,10);
        JA=MyButton("H2O");
	fix(JA,grd,gbc,2,13,1,1,20,10);
        JB=MyButton("Bnd");
	fix(JB,grd,gbc,3,13,1,1,20,10);
                
	// Register button events listeners

	AA.addActionListener(new GUIButtonHandler());
	AB.addActionListener(new GUIButtonHandler());
	BA.addActionListener(new GUIButtonHandler());
	BB.addActionListener(new GUIButtonHandler());
	CA.addActionListener(new GUIButtonHandler());
	CB.addActionListener(new GUIButtonHandler());
	DA.addActionListener(new GUIButtonHandler());
	DB.addActionListener(new GUIButtonHandler());
	EA.addActionListener(new GUIButtonHandler());
	EB.addActionListener(new GUIButtonHandler());
	FA.addActionListener(new GUIButtonHandler());
	FB.addActionListener(new GUIButtonHandler());
	GA.addActionListener(new GUIButtonHandler());
	GB.addActionListener(new GUIButtonHandler());
	HA.addActionListener(new GUIButtonHandler());
	HB.addActionListener(new GUIButtonHandler());
	IA.addActionListener(new GUIButtonHandler());
	IB.addActionListener(new GUIButtonHandler());
	JA.addActionListener(new GUIButtonHandler());
	JB.addActionListener(new GUIButtonHandler());
    }

    void displayEditor()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to define the drawing buttons

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/	
	AC=MyButton("Drw");
	fix(AC,grd,gbc,4,0,1,1,20,10);
        BC=MyButton("Lnk");
	fix(BC,grd,gbc,4,1,1,1,20,10);
        CC=MyButton("Del");
	fix(CC,grd,gbc,4,3,1,1,20,10);
	DC=MyButton("ADH");
	fix(DC,grd,gbc,4,4,1,1,20,10);
        EC=MyButton("Grp");
	fix(EC,grd,gbc,4,6,1,1,20,10);
        FC=MyButton("Opt");
	fix(FC,grd,gbc,4,7,1,1,20,10);
        GC=MyButton("Sav");
	fix(GC,grd,gbc,4,9,1,1,20,10);
        HC=MyButton("Dup");
	fix(HC,grd,gbc,4,10,1,1,20,10);
        IC=MyButton("Box");
	fix(IC,grd,gbc,4,12,1,1,20,10);
        JC=MyButton("Frg");
	fix(JC,grd,gbc,4,13,1,1,20,10);

	// Register button events listeners

	AC.addActionListener(new GUIButtonHandler());
	BC.addActionListener(new GUIButtonHandler());
	CC.addActionListener(new GUIButtonHandler());
	DC.addActionListener(new GUIButtonHandler());
	EC.addActionListener(new GUIButtonHandler());
	FC.addActionListener(new GUIButtonHandler());
	GC.addActionListener(new GUIButtonHandler());
	HC.addActionListener(new GUIButtonHandler());
	IC.addActionListener(new GUIButtonHandler());
	JC.addActionListener(new GUIButtonHandler());

	// Add Editor menu to menu bar

	top.add(Editor);
	setJMenuBar(top);

    }	

    void disposeEditor()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to remove the drawing buttons

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/	
	// Activate safety backup of edited molecule

	pane.backupEditFile();

	// Remove editor buttons from panel

	getContentPane().remove(AC);
	getContentPane().remove(BC);
	getContentPane().remove(CC);
	getContentPane().remove(DC);
	getContentPane().remove(EC);
	getContentPane().remove(FC);
	getContentPane().remove(GC);
	getContentPane().remove(HC);
	getContentPane().remove(IC);
	getContentPane().remove(JC);

	// Remove Editor menu from menu bar

	top.remove(Editor);
	setJMenuBar(top);

    }

    // Handle Button events

    class GUIButtonHandler implements ActionListener
    {
	/*
*********************************************************************

dl_poly/java GUI class to handle button events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	public void actionPerformed(ActionEvent e)
	{
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    String arg = e.getActionCommand();

	    if (e.getSource() instanceof JButton)
		{
		    if(arg.equals("New"))
			{
			    config=new Config(job,ftype);
			    pane.restore();
			}
		    else if(arg.equals("Clr"))
			{
			    config=null;
			    pane.restore();
			}       
		    else if(arg.equals("Tx+"))
			switchTxp();
		    else if(arg.equals("Tx-"))
			switchTxm();
		    else if(arg.equals("Ty+"))
			switchTyp();
		    else if(arg.equals("Ty-"))
			switchTym();
		    else if(arg.equals("Tz+"))
			switchTzp();
		    else if(arg.equals("Tz-"))
			switchTzm();
		    else if(arg.equals("Tra"))
			switchTranslate();
		    else if(arg.equals("Rx+"))
			switchRxp();
		    else if(arg.equals("Rx-"))
			switchRxm();
		    else if(arg.equals("Ry+"))
			switchRyp();
		    else if(arg.equals("Ry-"))
			switchRym();
		    else if(arg.equals("Rz+"))
			switchRzp();
		    else if(arg.equals("Rz-"))
			switchRzm();
		    else if(arg.equals("Rot"))
			switchRotate();
		    else if(arg.equals("Rst"))
			switchRestore();
		    else if(arg.equals("Edt"))
			switchEdit();
		    else if(arg.equals("Drw"))
			switchDraw();
		    else if(arg.equals("Bnd"))
			switchBonds();
		    else if(arg.equals("Grp"))
			switchGroup();
		    else if(arg.equals("Opt"))
			switchOpt();
		    else if(arg.equals("Sav"))
			pane.saveEdit();
		    else if(arg.equals("Del"))
			switchDelete();
		    else if(arg.equals("ADH"))
			switchAddHydrogen();
		    else if(arg.equals("Lnk"))
			switchAddLinks();
		    else if(arg.equals("H2O"))
			pane.showWater();
		    else if(arg.equals("Dup"))
			switchDuplicate();
		    else if(arg.equals("Box"))
			switchBoxEditor();
		    else if(arg.equals("Frg"))
			switchAddFragment();
		}
	}
    }

    // Handle Menu events

    class GUIMenuHandler implements ActionListener
    {
	/*
*********************************************************************

dl_poly/java GUI class to handle menu and button events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	public void actionPerformed(ActionEvent e)
	{
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    String arg = e.getActionCommand();

	    if (e.getSource() instanceof JMenuItem)
		{
		    if(arg.equals("Quit"))
			{
			    System.exit(0);
			}
		    else if(arg.equals("View File"))
			{
			    getViewFile();
			}
		    else if(arg.equals("Delete file"))
			{
			    zappFile();
			}
		    else if(arg.equals("Defaults"))
			{
			    if(newdefs != null)
				newdefs.job.hide();
			    newdefs=new ChgDefaults(job);
			}
		    else if(arg.equals("Reset"))
			{
			    startGUI();
			    pane.repaint();
			    monitor.clearScreen();
			    monitor.println("GUI variables reinitialized");
			    }
		    else if(arg.equals("Print"))
			{
			    pane.printOut();
			}
		    else if(arg.equals("CONTROL"))
			{
			    if(makctr != null)
				makctr.byebye();
			    makctr=new MakeControl(job);
			}
		    else if(arg.equals("CFG"))
			{
			    ftype="CFG";
			    monitor.println("CFG file input selected");
			}
		    else if(arg.equals("XYZ"))
			{
			    ftype="XYZ";
			    monitor.println("XYZ file input selected");
			}
		    else if(arg.equals("SEQ"))
			{
			    ftype="SEQ";
			    monitor.println("SEQ file input selected");
			}
		    else if(arg.equals("MSI"))
			{
			    ftype="MSI";
			    monitor.println("MSI file input selected");
			}
		    else if(arg.equals("MDR"))
			{
			    ftype="MDR";
			    monitor.println("MDR file input selected");
			}
		    else if(arg.equals("Lattice"))
			{
			    if(maklat != null)
				maklat.job.hide();
				maklat=new MakeLattice(job);
				}
		    else if(arg.equals("Bucky"))
			{
			    if(makbuk != null)
				makbuk.job.hide();
			    makbuk=new MakeBucky(job);
			}
		    else if(arg.equals("Polymer"))
			{
			    if(makpol != null)
				makpol.job.hide();
			    makpol=new MakePoly(job);
			}
		    else if(arg.equals("Chain"))
			{
			    if(makchain != null)
				makchain.job.hide();
			    makchain=new MakeChain(job);
			}
		    else if(arg.equals("N_fold"))
			{
			    if(enfold != null)
				enfold.job.hide();
			    enfold=new Nfold(job);
			}
		    else if(arg.equals("BondLengths"))
			{
			    if(bondlen != null)
				bondlen.job.hide();
			    bondlen=new BondLengths(job);
			}
		    else if(arg.equals("Add Water"))
			{
			    if(addh2o != null)
				addh2o.job.hide();
			    addh2o=new WaterAdd(job);
			}
		    else if(arg.equals("Blank"))
			{
			    if(makblank != null)
				makblank.job.hide();
			    if(edit && !pane.safe)
				pane.backupEditFile();
			    makblank=new MakeBlankField(job);
			}
		    else if(arg.equals("Dreiding"))
			{
			    if(makdrei != null)
				makdrei.job.hide();
			    if(edit && !pane.safe)
				pane.backupEditFile();
			    makdrei=new MakeDreiField(job);
			}
		    else if(arg.equals("OPLS"))
			{
			    if(makopls != null)
				makopls.job.hide();
			    if(edit && !pane.safe)
				pane.backupEditFile();
			    makopls=new MakeOPLSField(job);
			}
		    else if(arg.equals("Ceramics"))
			{
			    if(makceram != null)
				makceram.job.hide();
			    if(edit && !pane.safe)
				pane.backupEditFile();
			    makceram=new MakeCeramField(job);
			    }
		    else if(arg.equals("Table"))
			{
			    if(maktable != null)
				maktable.job.hide();
			    maktable=new MakeTable(job);
			    }
		    else if(arg.equals("Run DL_POLY"))
			{
			    if(runjob != null)
				runjob.job.hide();
			    runjob=new Execute(job);
			}
		    else if(arg.equals("Store/Fetch Data"))
			{
			    if(datarc != null)
				datarc.job.hide();
			    datarc=new DataArchiver(job);
			}
		    else if(arg.equals("RDF_Plot"))
			{
			    if(rdfplt != null)
				rdfplt.job.hide();
			    rdfplt=new RDFPlot(job);
			}
		    else if(arg.equals("RDF_Calc"))
			{
			    if(rdfcal != null)
				rdfcal.job.hide();
			    rdfcal=new RDFCalc(job);
			}
		    else if(arg.equals("S(k)"))
			{
			    if(sokplt != null)
				sokplt.job.hide();
			    sokplt=new SokPlot(job);
			}
		    else if(arg.equals("Z_Density"))
			{
			    if(zdnplt != null)
				zdnplt.job.hide();
			    zdnplt=new ZdenPlot(job);
			}
		    else if(arg.equals("Slice"))
			{
			    if(slcrev != null)
				slcrev.job.hide();
			    slcrev=new Slice(job);
			}
		    else if(arg.equals("MSD"))
			{
			    if(msdrun != null)
				msdrun.job.hide();
			    msdrun=new RunMSD(job);
			}
		    else if(arg.equals("VAF"))
			{
			    if(vafrun != null)
				vafrun.job.hide();
			    vafrun=new RunVAF(job);
			}
		    else if(arg.equals("FAF"))
			{
			    if(fafrun != null)
				fafrun.job.hide();
			    fafrun=new RunFAF(job);
			}
		    else if(arg.equals("Gs(r,t)"))
			{
			    if(gslcal != null)
				gslcal.job.hide();
			    gslcal=new GslCalc(job);
			}
		    else if(arg.equals("Gd(r,t)"))
			{
			    if(gdfcal != null)
				gdfcal.job.hide();
			    gdfcal=new GdfCalc(job);
			}
		    else if(arg.equals("S(k,w)"))
			{
			    if(skwcal != null)
				skwcal.job.hide();
			    skwcal=new SkwCalc(job);
			}
		    else if(arg.equals("CONFIG"))
			{
			    config=new Config(job,"CONFIG");
			    pane.restore();
			}
		    else if(arg.equals("REVCON"))
			{
			    config=new Config(job,"REVCON");
			    pane.restore();
			}
		    else if(arg.equals("Plot"))
			{
			    if(graf != null)
				graf.job.hide();
			    graf=new GraphDraw(job);
			}
		    else if(arg.equals("What Atoms?"))
			{
			    monitor.println("Select required CONFIG file for input");
			    if((fname=selectFileNameContains(job,"CFG"))!=null)
				{
				    whatAtoms(fname);
				}
			    else
				{
				    monitor.println("File selection cancelled");
				}
			}
		    else if(arg.equals("Statistics"))
			{
			    if(graf != null)
				graf.job.hide();
			    graf=new GraphDraw(job);
			    if(staprp != null)
				staprp.job.hide();
			    staprp=new StatProp(job);
			    }
		    else if(arg.equals("About DL_POLY"))
			{
			    viewResource("About_DL_POLY");
			}
		    else if(arg.equals("Disclaimer"))
			{
			    viewResource("Disclaimer");
			}
		    else if(arg.equals("Licence"))
			{
			    viewResource("Licence");
			}
		    else if(arg.equals("Acknowledgements"))
			{
			    viewResource("Acknowledge");
			}
		    else if(arg.equals("MINIDREI"))
			{
			    viewResource("MINIDREI");
			}
		    else if(arg.equals("MINIOPLS"))
			{
			    viewResource("MINIOPLS");
			}
		    else if(arg.equals("CERAMICS"))
			{
			    viewResource("CERAMICS");
			}
		    else if(arg.equals("Clear Text"))
			{
			    monitor.clearScreen();
			}
		    else if (arg.equals("C_1"))
			{
			    pane.setDrawAtom("C_1");
			}
		    else if (arg.equals("C_2"))
			{
			    pane.setDrawAtom("C_2");
			}
		    else if (arg.equals("C_R"))
			{
			    pane.setDrawAtom("C_R");
			}
		    else if (arg.equals("C_3"))
			{
			    pane.setDrawAtom("C_3");
			}
		    else if (arg.equals("O_2"))
			{
			    pane.setDrawAtom("O_2");
			}
		    else if (arg.equals("O_3"))
			{
			    pane.setDrawAtom("O_3");
			}
		    else if (arg.equals("H_"))
			{
			    pane.setDrawAtom("H_");
			}
		    else if (arg.equals("N_1"))
			{
			    pane.setDrawAtom("N_1");
			}
		    else if (arg.equals("N_2"))
			{
			    pane.setDrawAtom("N_2");
			}
		    else if (arg.equals("N_3"))
			{
			    pane.setDrawAtom("N_3");
			}
		    else if (arg.equals("S_2"))
			{
			    pane.setDrawAtom("S_2");
			}
		    else if (arg.equals("S_3"))
			{
			    pane.setDrawAtom("S_3");
			}
		    else if (arg.equals("P_2"))
			{
			    pane.setDrawAtom("P_2");
			}
		    else if (arg.equals("P_3"))
			{
			    pane.setDrawAtom("P_3");
			}
		    else if (arg.equals("OW"))
			{
			    pane.setDrawAtom("OW");
			}
		    else if (arg.equals("HW"))
			{
			    pane.setDrawAtom("HW");
			}
		    else if (arg.equals("H__HB"))
			{
			    pane.setDrawAtom("H__HB");
			}
		    else if (arg.equals("None"))
			{
			    pane.setBoxType("NON");
			}
		    else if (arg.equals("Cubic"))
			{
			    pane.setBoxType("CUB");
			}
		    else if (arg.equals("O-Rhombic"))
			{
			    pane.setBoxType("ORH");
			}
		    else if (arg.equals("T-Octahedral"))
			{
			    pane.setBoxType("OCT");
			}
		    else if (arg.equals("R-Dodecahedral"))
			{
			    pane.setBoxType("DEC");
			}
		    else if (arg.equals("Hexagonal"))
			{
			    pane.setBoxType("HEX");
			}
		    else if (arg.equals("Search"))
			{
			    pane.setFragmentType("SELECTED");
			}
		    else if (arg.equals("Alanine"))
			{
			    pane.setFragmentType("ALANINE");
			}
		    else if (arg.equals("Benzene"))
			{
			    pane.setFragmentType("BENZENE");
			}
		    else if (arg.equals("Naphthalene"))
			{
			    pane.setFragmentType("NAPHTHALENE");
			}
		    else if (arg.equals("n-Hexane"))
			{
			    pane.setFragmentType("n_HEXANE");
			}
		    else if (arg.equals("n-Decane"))
			{
			    pane.setFragmentType("n_DECANE");
			}
		    else if (arg.equals("n-Butane"))
			{
			    pane.setFragmentType("n_BUTANE");
			}
		    else if (arg.equals("c-Hexane"))
			{
			    pane.setFragmentType("c_HEXANE");
			}
		    else if (arg.equals("i-Butane"))
			{
			    pane.setFragmentType("i_BUTANE");
			}
		    else if (arg.equals("Styrene"))
			{
			    pane.setFragmentType("STYRENE");
			}
		    else if (arg.equals("Glucose"))
			{
			    pane.setFragmentType("GLUCOSE");
			}
		}
	}

	// Handle item events

	public void itemStateChanged(ItemEvent ie)
	{
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    pane.repaint();
	}
    }
    
    // Main method

    public static void main(String args[])
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	// Select colour scheme

	String scheme="picasso";
	if(args.length>0)scheme=args[0].toLowerCase();
	art=new ColorScheme(scheme);

	// Set up Graphical User interface

	job = new GUI();
	job.pack();
	job.show();
	monitor = new Monitor(job);
	monitor.println(page);
	startGUI();
    }
    static void startGUI()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/

	// Initialise GUI variables

	ftype="CFG";
	bondpc=BONDPERCENT;
	rotdef=DEFINE_ROTATION;
	tradef=DEFINE_TRANSLATION;
	incx=tradef;
	incy=tradef;
	incz=tradef;
	cal=Math.cos(Math.PI*rotdef/180.0);
	sal=Math.sin(Math.PI*rotdef/180.0);

	// Define default bondlengths

	cc1b=1.54;
	cc2b=1.34;
	cc3b=1.20;
	ccab=1.39;
	ch1b=1.09;
	cn1b=1.47;
	cn2b=1.27;
	cn3b=1.16;
	cnab=1.35;
	co1b=1.43;
	co2b=1.22;
	coab=1.36;
	nh1b=1.01;
	oh1b=0.96;
	if (bondlen != null)
	    bondlen.bondLengthSet();

	// File numbers

	numlat=0;
	numpol=0;
	numbuk=0;
	numblk=0;
	numchn=0;
	numctr=0;
	numtab=0;
	numxyz=0;
	nummsi=0;
	nummdr=0;
	numseq=0;
	numrdf=0;
	numsok=0;
	numzdn=0;
	nummsd=0;
	numvaf=0;
	numfaf=0;
	numstat=0;
	numgdf=0;
	numhov=0;
	numgsl=0;
	numskw=0;
	numh2o=0;
	numdre=0;
	numopl=0;
	numcer=0;
	numsav=0;

	// nullify processes

	proc=null;
	prdf=null;
	pgdf=null;
	pskw=null;
    }

    void setPaneParams()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.cal=cal;
	pane.sal=sal;
	pane.incx=incx;;
	pane.incy=incy;
	pane.incz=incz;
	pane.edit=edit;
	pane.fragment=null;
    }

    void switchEdit()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		edit=false;
		disposeEditor();
		setPaneParams();
	    }
	else
	    {
		edit=true;
		displayEditor();
		setPaneParams();
		pane.oper=0;
		pane.news="NULL";
		pane.config=null;
		pane.setDrawAtom("C_3");
	    }
	job.pack();
	job.show();
	pane.restore();
	pane.redraw();
    }

    void switchDuplicate()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 11)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=11;
		pane.news="DUPLICATE";
	    }
	pane.redraw();
    }

    void switchBoxEditor()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 12)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=12;
		pane.defineBox();
		pane.news="EDIT BOX";
	    }
	pane.redraw();
    }

    void switchDraw()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 1)
	    {
		pane.oper=0;
		pane.lpen=false;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=1;
		pane.lpen=true;
		pane.link=false;
		pane.news="DRAW";
	    }
	pane.redraw();
    }

    void switchBonds()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/

	if(edit)
	    {
		if(pane.oper == 14)
		    {
			pane.oper=0;
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=14;
			pane.news="BONDS";
		    }
	    }
	else
	    {
		config.switchBonds();
		pane.makeBonds();
	    }
	pane.redraw();
    }

    void switchAddFragment()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(pane.oper == 13)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=13;
		pane.news="FRAGMENT";
	    }
	pane.redraw();
    }

    void switchAddLinks()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 2)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=2;
		pane.news="LINK";
	    }
	pane.redraw();
    }

    void switchAddHydrogen()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 10)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=10;
		pane.news="ADD/DELETE HYDROGEN";
	    }
	pane.redraw();
    }

    void switchGroup()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.mark0=-1;
	pane.mark1=-1;
	pane.mark2=-1;
	if(pane.oper == 5)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=5;
		pane.news="GROUP";
	    }
	pane.redraw();
    }

    void switchOpt()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(pane.oper == 9)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=9;
		pane.news="OPTIMIZATION";
	    }
	pane.redraw();
    }

    void switchDelete()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(pane.oper == 3)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=3;
		pane.mark0=-1;
		pane.mark1=-1;
		pane.mark2=-1;
		pane.news="DELETE";
	    }
	pane.redraw();
    }

    void switchTzm()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 40)
		    {
			pane.oper=0;
			pane.act="";
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=40;
			pane.act="Tz-";
			pane.news="SHIFT Z-";
		    }
	    }
	else
	    {
		pane.displace(2,-incz);
	    }
	pane.redraw();
    }

    void switchTzp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 41)
		    {
			pane.oper=0;
			pane.act="";
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=41;
			pane.act="Tz+";
			pane.news="SHIFT Z+";
		    }
	    }
	else
	    {
		pane.displace(2,incz);
	    }
	pane.redraw();
    }

    void switchTxm()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 42)
		    {
			pane.oper=0;
			pane.act="";
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=42;
			pane.act="Tx-";
			pane.news="SHIFT X-";
		    }
	    }
	else
	    {
		pane.displace(0,-incx);
	    }
	pane.redraw();
    }

    void switchTxp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 43)
		    {
			pane.oper=0;
			pane.act="";
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=43;
			pane.act="Tx+";
			pane.news="SHIFT X+";
		    }
	    }
	else
	    {
		pane.displace(0,incx);
	    }
	pane.redraw();
    }

    void switchTym()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 44)
		    {
			pane.oper=0;
			pane.act="";
			pane.news="NULL";
		    }
		else
		    {
			pane.oper=44;
			pane.act="Ty-";
			pane.news="SHIFT Y-";
		    }
	    }
	else
	    {
		pane.zoom(-1);
	    }
	pane.redraw();
    }

    void switchTyp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	{
	    if(pane.oper == 45)
		{
		    pane.oper=0;
		    pane.act="";
		    pane.news="NULL";
		}
	    else
		{
		    pane.oper=45;
		    pane.act="Ty+";
		    pane.news="SHIFT Y+";
		}
	}
	else
	    {
		pane.zoom(1);
	    }
	pane.redraw();
    }

    void switchTranslate()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(pane.oper == 6)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=6;
		pane.news="MOVE";
	    }
	pane.redraw();
    }

    void switchRzm()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Rz-"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Rz-";
				pane.news="ROTATE Z";
			    }
		    }
		else
		    {
			if(pane.oper == 80)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=80;
				pane.act="Rz-";
				pane.news="ROTATE Z-";
			    }
		    }
	    }
	else
	    {
		pane.rotate(0,1,cal,-sal);
	    }
	pane.redraw();
    }

    void switchRzp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Rz+"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Rz+";
				pane.news="ROTATE Z";
			    }
		    }
		else
		    {
			if(pane.oper == 81)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=81;
				pane.act="Rz+";
				pane.news="ROTATE Z+";
			    }
		    }
	    }
	else
	    {
		pane.rotate(0,1,cal,sal);
	    }
	pane.redraw();
    }

    void switchRxm()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Rx-"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Rx-";
				pane.news="ROTATE X";
			    }
		    }
		else
		    {
			if(pane.oper == 82)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=82;
				pane.act="Rx-";
				pane.news="ROTATE X-";
			    }
		    }
	    }
	else
	    {
		pane.rotate(1,2,cal,-sal);
	    }
	pane.redraw();
    }

    void switchRxp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Rx+"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Rx+";
				pane.news="ROTATE X";
			    }
		    }
		else
		    {
			if(pane.oper == 83)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=83;
				pane.act="Rx+";
				pane.news="ROTATE X+";
			    }
		    }
	    }
	else
	    {
		pane.rotate(1,2,cal,sal);
	    }
	pane.redraw();
    }

    void switchRym()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Ry-"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Ry-";
				pane.news="ROTATE Y";
			    }
		    }
		else
		    {
			if(pane.oper == 84)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=84;
				pane.act="Ry-";
				pane.news="ROTATE Y-";
			    }
		    }
	    }
	else
	    {
		pane.rotate(2,0,cal,-sal);
	    }
	pane.redraw();
    }
    
    void switchRyp()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    {
		if(pane.oper == 7)
		    {
			if(pane.act.equals("Ry+"))
			    {
				pane.act="";
				pane.news="ROTATE";
			    }
			else
			    {
				pane.act="Ry+";
				pane.news="ROTATE Y";
			    }
		    }
		else
		    {
			if(pane.oper == 85)
			    {
				pane.oper=0;
				pane.act="";
				pane.news="NULL";
			    }
			else
			    {
				pane.oper=85;
				pane.act="Ry+";
				pane.news="ROTATE Y+";
			    }
		    }
	    }
	else
	    {
		pane.rotate(2,0,cal,sal);
	    }
	pane.redraw();
    }

    void switchRotate()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	pane.act="";
	if(pane.oper == 7)
	    {
		pane.oper=0;
		pane.news="NULL";
	    }
	else
	    {
		pane.oper=7;
		pane.news="ROTATE";
	    }
	pane.redraw();
    }
    
    void switchRestore()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(edit)
	    pane.editRestore();
	else
	    pane.restore();
    }

    void getViewFile()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	
	// View a text file
	
	JFileChooser fc = new JFileChooser(new File("./"));
	fc.setDialogTitle("Select file for viewing");
	fc.setForeground(fore);
	fc.setBackground(back);
	int key=fc.showDialog(job,"View");
	if(key == JFileChooser.APPROVE_OPTION)
	    {
		fname=fc.getSelectedFile().getPath();
		viewFile(fname);
	    }
	else if(key == JFileChooser.CANCEL_OPTION)
	    {
		monitor.println("File selection cancelled");
	    }
	else
	    {
		monitor.println("Error - file viewer failure");
	    }
    }
    void viewResource(String fileName)
    {
	/*
*********************************************************************

dl_poly/java routine to display a Jar resource file

copyright - daresbury laboratory
author    - w.smith may 2005

*********************************************************************
 */
	monitor.clearScreen();

	try
	    {
		String text=""; 
		//ClassLoader c = getClass().getClassLoader();
		//InputStream instream = c.getResourceAsStream(fileName);
		File stuff = new File("../java/"+fileName);
		FileInputStream instream = new FileInputStream(stuff);
		InputStreamReader isr = new InputStreamReader(instream);
		BufferedReader reader = new BufferedReader(isr);
		while((text=reader.readLine()) != null)
		    {
			monitor.println(text);
		    }
		reader.close();
	    }
	catch(Exception e)
	    {
		System.err.println("Error: reading file " + fileName);
	    }
    }
    void zappFile()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	
	// Delete a file selected from browser
	
	monitor.println("Select file for deletion");
	JFileChooser fc = new JFileChooser(new File("./"));
	fc.setDialogTitle("Select file for deletion");
	fc.setForeground(fore);
	fc.setBackground(back);
	int key=fc.showDialog(job,"Delete");
	if(key == JFileChooser.APPROVE_OPTION)
	    {
		WarningBox danger=new WarningBox(job,"Warning!",true);
		danger.show();
		if(alert)
		    {
			fname=fc.getSelectedFile().getPath();
			(new File(fname)).delete();
			monitor.println("File "+fname+" has been deleted");
		    }
		else
		    {
			monitor.println("File deletion aborted");
		    }
	    }
	else if(key == JFileChooser.CANCEL_OPTION) 
	    {
		monitor.println("File selection cancelled");
	    }
	else
	    {
		monitor.println("Error - file deletion failure");
	    }
    }
    void draggedResponse(int a,int b){}
    JMenu MyMenu(String s)
    {
	    JMenu  mine = new JMenu(s);  
	    mine.setFont(fontMain);
	    mine.setForeground(fore);
	    mine.setBackground(back);
	    return mine;
    }
    JMenuItem MyMenuItem(String s)
    {
	    JMenuItem  mine = new JMenuItem(s);  
	    mine.setFont(fontMain);
	    mine.setForeground(fore);
	    mine.setBackground(back);
	    return mine;
    }
    JButton MyButton(String s)
    {
	    JButton mine=new JButton(s);
	    mine.setBackground(butn);
	    mine.setForeground(butf);
	    return mine;
    }


    // Support classes

    class GUIWindowAdapter implements WindowListener
    {
	/*
*********************************************************************

dl_poly/java GUI class for window events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	public void windowClosing(WindowEvent e)
        {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    if(edit && !pane.safe)
		pane.backupEditFile();
	    System.exit(0);
        }
	public void windowActivated(WindowEvent e){}
	public void windowClosed(WindowEvent e){}
	public void windowDeactivated(WindowEvent e){}
	public void windowDeiconified(WindowEvent e){}
	public void windowIconified(WindowEvent e){}
	public void windowOpened(WindowEvent e){}
    }

}
