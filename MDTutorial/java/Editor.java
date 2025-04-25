import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.print.*;
import javax.swing.*;
import java.util.*;

public class Editor extends JComponent implements Printable
{
	/*
*********************************************************************

dl_poly/java GUI class for rendering molecular pictures

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
    GUI home;
    Monitor monitor;
    Image offscreenImg;
    Graphics offscreenG;
    Config config=null,cfgsav=null,fragment=null;
    boolean edit=false,link,lpen,water=true,lmove=false,safe=true;
    int oper=0,nangs,tilex,tiley,boxtyp=1,nhdel=0,pen=0;
    int mxgroup,ngroup,sx,sy,fx,fy,mark0=-1,mark1=-1,mark2=-1;
    int[] lst,xpos,ypos,zpos,dia,xvrt,yvrt,zvrt,group;
    double scale,cal,sal,xx2,yy2,zz2,dd2,fac,incx,incy,incz;
    boolean[] led;
    String news,act,frg;
    Element atom;

    public Editor(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	home=here;
	monitor=home.monitor;

	// Register Mouse Listeners

	addMouseListener(new MousePoints());
	addMouseMotionListener(new MouseMoves());

	// Register object for screen resize
	
       	addComponentListener(new MyScreenSizer());
	
    }

    public void printOut()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to print the GUI

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	PrinterJob prnjob=PrinterJob.getPrinterJob();
	PageFormat pf=prnjob.defaultPage();
	pf.setOrientation(PageFormat.LANDSCAPE);
	prnjob.setPrintable(this,pf);
	try
	    {
		if(prnjob.printDialog())
		    {
			monitor.println("Initiating print .......");
			prnjob.print();
			monitor.println("Print complete");
		    }
		else
		    {
			monitor.println("No print initiated");
		    }
	    }
	catch(PrinterException pe)
	    {
		monitor.println("Error - problem with print "+pe.toString());
	    }

    }

    public int print(Graphics g,PageFormat pf,int pageIndex)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	if(pageIndex > 0)
	    {
		return Printable.NO_SUCH_PAGE;
	    }
	Graphics2D g2d=(Graphics2D)g;
	g2d.translate(pf.getImageableX(),pf.getImageableY());
	paint(g2d);
	return Printable.PAGE_EXISTS;
    }

    public void paint(Graphics g) 
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	boolean showbond;
	int j,k,u,w,rr,xd,yd,zd,ka,kb,ja,jb;

	tilex=getSize().width;
	tiley=getSize().height;
	offscreenImg = createImage(tilex,tiley);
	offscreenG = offscreenImg.getGraphics();
	offscreenG.setColor(home.scrn);
	offscreenG.fillRect(0,0,tilex,tiley);
	offscreenG.setColor(Color.black);
	g.setFont(new Font("Monaco",Font.PLAIN,10));
	if(edit)
	    {
		offscreenG.drawString("EDIT: "+news,10,10);
		offscreenG.drawString(atom.zsym,tilex-50,10);
		offscreenG.drawLine(tilex/2,tiley/2+10,tilex/2,tiley/2-10);
		offscreenG.drawLine(tilex/2+10,tiley/2,tilex/2-10,tiley/2);
		if(oper == 13 && fragment != null)
		    offscreenG.drawString("FRAG: "+frg,10,tiley-10);
	    }
	else if(oper == 6 || oper == 7)
	    offscreenG.drawString("VIEW: "+news,10,10);
	else
	    offscreenG.drawString("VIEW",10,10);
	if(config != null && config.natms > 0)
	    {
		//      order in z direction
	
		Super.ShellSort(config.natms,lst,zpos);
	
		//      Mark clicked atom
		
		if(mark0>=0)
		    {
			w=dia[mark0];
			offscreenG.setColor(Color.red);
			offscreenG.fillOval(xpos[mark0]-4,ypos[mark0]-4,w+8,w+8);
		    }
		offscreenG.setColor(Color.black);
		
		if(config.imcon > 0)
		    {
			//      check visibility of cell edges
			
			checkEdges();
			
			//      Draw hidden box edges
			
			for (int i=0;i<config.nedge;i++)
			    {
				if(!led[i])
				    {
					j=config.edge[0][i];
					k=config.edge[1][i];
					offscreenG.drawLine(xvrt[j],yvrt[j],xvrt[k],yvrt[k]);
				    }
			    }
			
		    }

		//      Draw atoms
		
		for (int i=0;i<config.natms;i++)
		    {
			if((water ||!(config.atoms[lst[i]].zsym.indexOf("OW") == 0 ||
				     config.atoms[lst[i]].zsym.indexOf("HW") == 0)) &&
			             !(config.atoms[lst[i]].zsym.indexOf("QW") == 0))
			    {
				j=lst[i];
				w=dia[j];
				if(j == mark0 || j == mark1 || j == mark2 || (ngroup > 0 && group[j] > -1))
				    {
					offscreenG.setColor(Color.red);
					offscreenG.fillOval(xpos[j]-4,ypos[j]-4,w+8,w+8);
				    }
				offscreenG.fillOval(xpos[j]-2,ypos[j]-2,w+4,w+4);
				offscreenG.setColor(config.atoms[j].zcol);
				offscreenG.fillOval(xpos[j],ypos[j],w,w);
				offscreenG.setColor(Color.black);
				w=w/2;
				if(config.nbnds>0)
				    {
					for(int m=0;m<config.lbnd[j];m++)
					    {	
						showbond=true;
						k=config.bond[m][j];
						if(config.imcon > 0)
						    {
							if(Math.abs((double)(xpos[k]-xpos[j])) > scale*4)showbond=false;
							if(Math.abs((double)(ypos[k]-ypos[j])) > scale*4)showbond=false;
							if(Math.abs((double)(zpos[k]-zpos[j])) > scale*4)showbond=false;
						    }
						if(config.atoms[k].zsym.indexOf("QW") == 0)showbond=false;
						if(showbond)
						    {
							u=dia[k]/2;
							if(j>k && zpos[k]==zpos[j])
							    {
								xd=xpos[k]-xpos[j]+u-w;
								yd=ypos[k]-ypos[j]+u-w;
								rr=(int)Math.sqrt(xd*xd+yd*yd);
								ja=xpos[j]+w+(int)((double)(xd*w)/rr);
								jb=ypos[j]+w+(int)((double)(yd*w)/rr);
								ka=xpos[k]+u-(int)((double)(xd*u)/rr);
								kb=ypos[k]+u-(int)((double)(yd*u)/rr);
								offscreenG.drawLine(ja,jb+1,ka,kb+1);
								offscreenG.drawLine(ja+1,jb,ka+1,kb);
								offscreenG.drawLine(ja,jb,ka,kb);
								offscreenG.drawLine(ja-1,jb,ka-1,kb);
								offscreenG.drawLine(ja,jb-1,ka,kb-1);
							    }
							else if(zpos[k]>zpos[j])
							    {
								xd=xpos[k]-xpos[j]+u-w;
								yd=ypos[k]-ypos[j]+u-w;
								zd=zpos[k]-zpos[j]+u-w;
								rr=(int)Math.sqrt(xd*xd+yd*yd+zd*zd);
								ja=xpos[j]+w+(int)((double)(xd*w)/rr);
								jb=ypos[j]+w+(int)((double)(yd*w)/rr);
								ka=xpos[k]+u;
								kb=ypos[k]+u;
								offscreenG.drawLine(ja,jb+1,ka,kb+1);
								offscreenG.drawLine(ja+1,jb,ka+1,kb);
								offscreenG.drawLine(ja,jb,ka,kb);
								offscreenG.drawLine(ja-1,jb,ka-1,kb);
								offscreenG.drawLine(ja,jb-1,ka,kb-1);
							    }
						    }
					    }
				    }
			    }
		    }

		if(config.imcon > 0)
		    {
			
			//      Draw non hidden box edges
			
			for (int i=0;i<config.nedge;i++)
			    {
				if(led[i])
				    {
					j=config.edge[0][i];
					k=config.edge[1][i];
					offscreenG.drawLine(xvrt[j],yvrt[j],xvrt[k],yvrt[k]);
				    }
			    }
		    }
		if(lmove && oper == 5)
		    {
			offscreenG.drawLine(sx,sy,fx,sy);
			offscreenG.drawLine(sx,sy,sx,fy);
			offscreenG.drawLine(sx,fy,fx,fy);
			offscreenG.drawLine(fx,sy,fx,fy);
		    }
	    }
	g.drawImage(offscreenImg,0,0,this);
    }

    void redraw()
    {	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	double rad;

	if(config != null)
	    {
		if(edit)
		    fac=0.3;
		else
		    fac=0.5;
		lst=new int[config.natms];
		xvrt=new int[config.nvrt];
		yvrt=new int[config.nvrt];
		zvrt=new int[config.nvrt];
		led=new boolean[config.nedge];
		xpos=new int[config.mxatms];
		ypos=new int[config.mxatms];
		zpos=new int[config.mxatms];
		dia=new int[config.mxatms];
		if(config.imcon > 0)
		    {
			if(edit)config.setBox();
			double size=BML.max(config.cell[0]+config.cell[3]+config.cell[6],
					    config.cell[1]+config.cell[4]+config.cell[7],
					    config.cell[2]+config.cell[5]+config.cell[8]);
			//scale=0.75*Math.min(tilex,tiley)/Math.max(size,1.0);
		    }
		for (int i=0;i<config.natms;i++)
		    {
			lst[i]=i;
			rad=fac*config.atoms[i].zrad;
			dia[i]=(int)(2*scale*rad);
			xpos[i]=(int)(scale*(config.xyz[0][i]-rad))+tilex/2;
			ypos[i]=-(int)(scale*(config.xyz[2][i]+rad))+tiley/2;
			zpos[i]=-(int)(scale*config.xyz[1][i]);
		    }
		for (int i=0;i<config.nvrt;i++)
		    {
			xvrt[i]=(int)(scale*config.vrt[0][i])+tilex/2;
			yvrt[i]=-(int)(scale*config.vrt[2][i])+tiley/2;
			zvrt[i]=-(int)(scale*config.vrt[1][i]);
		    }
	    }
	repaint();
    }

    void setDrawAtom(String atm)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/

	atom=new Element(atm);
	redraw();
    }

    void setFragmentType(String fra)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	boolean chk;
	frg=fra;

	if(fra.equals("SELECTED"))
	    {
		fragment=new Config(home,"CFG");
		if(fragment != null)
		    {
			if(config != null && config.nbnds > 0)
			    fragment.switchBonds();
		    }
	    }
	else
	    {
		fragment=new Config();
		if(fragment.rdCFG("../java/FRGLIB/"+fra+".FRG"))
		    {
			if(config != null && config.nbnds > 0)
			    fragment.switchBonds();
		    }
		else
		    fragment=null;
	    }
	redraw();
    }

    void setBoxType(String box)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/

	if(box.equals("NON"))
	    boxtyp=0;
	else if(box.equals("CUB"))
	    boxtyp=1;	
	else if(box.equals("ORH"))
	    boxtyp=2;
	else if(box.equals("OCT"))
	    boxtyp=4;
	else if(box.equals("DEC"))
	    boxtyp=5;
	else if(box.equals("HEX"))
	    boxtyp=7;
	redraw();
    }

    void makeBonds()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(config != null)
	    {
		config.switchBonds();
		if(config.nbnds > 0)
		    {
			getJoinArray();
		    }
	    }
	redraw();
    }

    void showWater()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	if(water)
	    {
		monitor.println("Water visibility is OFF");
		water=false;
	    }
	else
	    {
		monitor.println("Water visibility is ON");
		water=true;
	    }
	redraw();
    }

    void saveEdit()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	String fname="CFGEDT";

	cfgsav=copyConfig(config);
	monitor.println("Configuration saved");
	home.configWrite(fname,cfgsav);
	home.config=cfgsav;
	safe=true;
	redraw();
    }

    Config copyConfig(Config cfg)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	Config cfgcpy;
	int[] key,lok;
	int c,m,n,o,p,q;

	key=new int[cfg.natms];
	lok=new int[cfg.natms];

	cfgcpy=new Config(cfg.mxatms);

	cfgcpy.title=cfg.title;
	cfgcpy.natms=cfg.natms;
	cfgcpy.nbnds=cfg.nbnds;
	cfgcpy.imcon=cfg.imcon;
	cfgcpy.atoms=new Element[cfg.mxatms];
	cfgcpy.xyz=new double[3][cfg.mxatms];
	if(cfg.imcon > 0)
	    {
		for (int i=0;i<9;i++)
		    cfgcpy.cell[i]=cfg.cell[i];
		cfgcpy.setBox();
	    }
	for (int i=0;i<cfg.natms;i++)
	    {
		cfgcpy.atoms[i]=new Element(cfg.atoms[i].zsym);
		cfgcpy.xyz[0][i]=cfg.xyz[0][i];
		cfgcpy.xyz[1][i]=cfg.xyz[1][i];
		cfgcpy.xyz[2][i]=cfg.xyz[2][i];
	    }

	// construct contiguous version

	if(cfg.nbnds > 0)
	    {
		c=cfgcpy.getBonds();
		
		for (int i=0;i<cfgcpy.natms;i++)
		    {
			lok[i]=i;
			key[i]=-1;
		    }
		for (int i=0;i<cfgcpy.nbnds;i++)
		    {
			m=Math.min(cfgcpy.join[0][i],cfgcpy.join[1][i]);
			n=Math.max(cfgcpy.join[0][i],cfgcpy.join[1][i]);
			
			o=key[m];
			if(o < 0)
			    {
				o=m;
				key[m]=m;
			    }
			if(key[n] < 0)
			    {
				key[n]=o;
			    }
			else if (key[n] != o)
			    {
				p=Math.min(o,key[n]);
				q=Math.max(o,key[n]);
				
				for (int j=0;j<cfgcpy.natms;j++)
				    {
					if(key[j] == q) key[j]=p;
				    }
				o=p;
			    }
		    }
		
		Super.ShellSort(cfgcpy.natms,lok,key);
		
		for (int i=0;i<cfgcpy.natms;i++)
		    {
			p=lok[i];
			cfgcpy.atoms[i]=new Element(cfg.atoms[p].zsym);
			cfgcpy.xyz[0][i]=cfg.xyz[0][p];
			cfgcpy.xyz[1][i]=cfg.xyz[1][p];
			cfgcpy.xyz[2][i]=cfg.xyz[2][p];
		    }
		c=cfgcpy.getBonds();
	    }
	return cfgcpy;
    }

    void backupEditFile()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	String fname="CFGEDT";

	// Save backup data file

	if(!safe)
	    {
		cfgsav=copyConfig(config);
		home.config=cfgsav;
		try
		    {
			DataOutputStream outStream = new DataOutputStream(new FileOutputStream(fname));
			outStream.writeBytes(cfgsav.title+"\n");
			outStream.writeBytes(BML.fmt(0,10)+BML.fmt(cfgsav.imcon,10)+"\n");
			if(cfgsav.imcon > 0)
			    {
				outStream.writeBytes(BML.fmt(cfgsav.cell[0],20)+BML.fmt(cfgsav.cell[1],20)+BML.fmt(cfgsav.cell[2],20)+"\n");
				outStream.writeBytes(BML.fmt(cfgsav.cell[3],20)+BML.fmt(cfgsav.cell[4],20)+BML.fmt(cfgsav.cell[5],20)+"\n");
				outStream.writeBytes(BML.fmt(cfgsav.cell[6],20)+BML.fmt(cfgsav.cell[7],20)+BML.fmt(cfgsav.cell[8],20)+"\n");
			    }
			for (int k=0;k<cfgsav.natms;k++)
			    {
				outStream.writeBytes(BML.fmt(cfgsav.atoms[k].zsym,8)+BML.fmt(k+1,10)+BML.fmt(cfgsav.atoms[k].znum,10)+"\n");
				outStream.writeBytes(BML.fmt(cfgsav.xyz[0][k],20)+BML.fmt(cfgsav.xyz[1][k],20)+BML.fmt(cfgsav.xyz[2][k],20)+"\n");
			    }
			outStream.close();
		    }
		catch(Exception e)
		    {
			monitor.println("Error - writing backup file: "+fname);
		    }
		monitor.println("Edit backup file written: CFGEDT");
		safe=true;
	    }
    }

   void restore()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	oper=0;
	fac=0.5;
	ngroup=0;
	mark0=-1;
	mark1=-1;
	mark2=-1;
	safe=true;
	scale=37.5;

	if(home.config == null)
	    config=null;
	else
	    config=copyConfig(home.config);

	if(config != null)
	    {
		if(config.natms > 0)
		    getAtomicArrays();
		if(config.nbnds > 0)
		    {
			getJoinArray();
		    }
		if(config.imcon > 0)
		    {
			config.setBox();
			double size=BML.max(config.cell[0]+config.cell[3]+config.cell[6],
					    config.cell[1]+config.cell[4]+config.cell[7],
					    config.cell[2]+config.cell[5]+config.cell[8]);
			scale=0.75*Math.min(tilex,tiley)/Math.max(size,1.0);
		    }
	    }
	redraw();
    }

    void editRestore()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	fac=0.3;
	oper=0;
	ngroup=0;
	mark0=-1;
	mark1=-1;
	mark2=-1;
	scale=37.5;
	if(cfgsav == null)
	    {
		monitor.println("No edit backup available");
	    }
	else
	    {
		config=cfgsav;
		home.config=cfgsav;
		monitor.println("Edit backup restored");
		safe=true;
		lpen=false;
		link=false;
		news="NULL";
		if(config.natms > 0)
		    getAtomicArrays();
		if(config.nbnds > 0)
		    {
			getJoinArray();
		    }
		if(config.imcon > 0)
		    config.setBox();
		redraw();
	    }
    }

    void checkEdges()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to check visibility of box edges

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	int w,n1,n2,xa,xb,ya,yb;
	double px,py,pp,sx,sy,ss,uu,dd,ww,rr,zz,tt;

	if(config != null)
	    {
		for (int n=0;n<config.nedge;n++)
		    {
			led[n]=true;
			n1=config.edge[0][n];
			n2=config.edge[1][n];
			xa=Math.min(xvrt[n1],xvrt[n2]);
			xb=Math.max(xvrt[n1],xvrt[n2]);
			ya=Math.min(yvrt[n1],yvrt[n2]);
			yb=Math.max(yvrt[n1],yvrt[n2]);
		    OUT:
			for(int i=0;i<config.natms;i++)
			    {
				if(xpos[i]>xa-dia[i] && xpos[i]<xb && ypos[i]>ya-dia[i] && ypos[i]<yb)
				    {
					w=dia[i]/2;
					ww=0.25*(double)(dia[i]*dia[i]);
					px=(double)(xvrt[n2]-xvrt[n1]);
					py=(double)(yvrt[n2]-yvrt[n1]);
					pp=Math.sqrt(px*px+py*py);
					sx=(double)(xpos[i]+w-xvrt[n1]);
					sy=(double)(ypos[i]+w-yvrt[n1]);
					uu=(sx*px+sy*py)/pp;
					ss=sx*sx+sy*sy-uu*uu;
					if(ss < ww)
					    {
						dd=(double)(zpos[i]-zvrt[n1]);
						zz=(double)(zvrt[n2]-zvrt[n1]);
						tt=uu*zz/pp;
						if (tt < dd)
						    {
							led[n]=false;
							break OUT;
						    }
					    }
				    }
			    }
		    }
	    }
    }

    void zoom(int k)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double rad;

	if(edit)return;

	if (config != null)
	    {
		if(scale==0.0)
		    {
			if(config.imcon > 0)
			    {
				double size=BML.max(config.cell[0]+config.cell[3]+config.cell[6],
						    config.cell[1]+config.cell[4]+config.cell[7],
						    config.cell[2]+config.cell[5]+config.cell[8]);
				scale=0.75*Math.min(tilex,tiley)/Math.max(size,1.0);
			    }
			else
			    scale=37.5;
		    }
		else
		    {
			if(k > 0)
			    scale/=1.5;
			else
			    scale*=1.5;
		    }
		redraw();
	    }
    }
    
    void rescale()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double rad;

	if(edit)return;

	if(config!=null)
	    {
		if(config.imcon > 0)
		    {
			double size=BML.max(config.cell[0]+config.cell[3]+config.cell[6],
					    config.cell[1]+config.cell[4]+config.cell[7],
					    config.cell[2]+config.cell[5]+config.cell[8]);
			scale=0.75*Math.min(tilex,tiley)/Math.max(size,1.0);
		    }
		else
		    scale=37.5;
		redraw();
	    }
    }

    public void displace(int k,double inc)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double rad;

	if(edit)return;

	if(config != null)
	    {
		for (int i=0;i<config.natms;i++)
		    config.xyz[k][i]+=inc;
		for (int i=0;i<config.nvrt;i++)
		    config.vrt[k][i]+=inc;
		if (k == 0)
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				xpos[i]=(int)(scale*(config.xyz[0][i]-rad))+tilex/2;
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				xvrt[i]=(int)(scale*config.vrt[0][i])+tilex/2;
			    }
		    }
		else if (k == 1)
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				ypos[i]=-(int)(scale*(config.xyz[2][i]+rad))+tiley/2;
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				yvrt[i]=-(int)(scale*config.vrt[2][i])+tiley/2;
			    }
		    }
		else if (k == 2)
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				zpos[i]=-(int)(scale*config.xyz[1][i]);
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				zvrt[i]=-(int)(scale*config.vrt[1][i]);
			    }
		    }
		repaint();
	    }
    }

    public void rotate(int j,int k,double c,double s)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double u,v,rad;

	if(config != null)
	    {
		for (int i=0;i<config.natms;i++)
		    {
			u=config.xyz[j][i];
			v=config.xyz[k][i];
			config.xyz[j][i]=(c*u-s*v);
			config.xyz[k][i]=(s*u+c*v);
		    }
		for (int i=0;i<config.nvrt;i++)
		    {
			u=config.vrt[j][i];
			v=config.vrt[k][i];
			config.vrt[j][i]=(c*u-s*v);
			config.vrt[k][i]=(s*u+c*v);
		    }
		if ((j == 0)&&(k == 1))
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				xpos[i]=(int)(scale*(config.xyz[0][i]-rad))+tilex/2;
				ypos[i]=-(int)(scale*(config.xyz[2][i]+rad))+tiley/2;
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				xvrt[i]=(int)(scale*config.vrt[0][i])+tilex/2;
				yvrt[i]=-(int)(scale*config.vrt[2][i])+tiley/2;
			    }
		    }
		else if ((j == 2)&&(k == 0))
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				xpos[i]=(int)(scale*(config.xyz[0][i]-rad))+tilex/2;
				zpos[i]=-(int)(scale*config.xyz[1][i]);
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				xvrt[i]=(int)(scale*config.vrt[0][i])+tilex/2;
				zvrt[i]=-(int)(scale*config.vrt[1][i]);
			    }
		    }
		else if ((j == 1)&&(k == 2))
		    {
			for (int i=0;i<config.natms;i++)
			    {
				rad=fac*config.atoms[i].zrad;
				ypos[i]=-(int)(scale*(config.xyz[2][i]+rad))+tiley/2;
				zpos[i]=-(int)(scale*config.xyz[1][i]);
			    }
			for (int i=0;i<config.nvrt;i++)
			    {
				yvrt[i]=-(int)(scale*config.vrt[2][i])+tiley/2;
				zvrt[i]=-(int)(scale*config.vrt[1][i]);
			    }
		    }
		repaint();
	    }
    }

    void shiftAtoms()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double rad;

	if (config != null)
	    {

		for (int i=0;i<config.natms;i++)
		    {
			config.xyz[0][i]+=(fx-sx)/scale;
			config.xyz[2][i]-=(fy-sy)/scale;
		    }
		for (int i=0;i<config.nvrt;i++)
		    {
			config.vrt[0][i]+=(fx-sx)/scale;
			config.vrt[2][i]-=(fy-sy)/scale;
		    }
		redraw();
	    }
    }

    void turnAtoms()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double rad,ca,cb,sa,sb,tx,ty,tz;

	if(edit)return;

	if (config != null)
	    {
		ca=Math.cos((sx-fx)*Math.PI/tilex);
		sa=Math.sin((sx-fx)*Math.PI/tilex);
		cb=Math.cos((sy-fy)*Math.PI/tiley);
		sb=Math.sin((sy-fy)*Math.PI/tiley);
		for (int i=0;i<config.natms;i++)
		    {
			tx=config.xyz[0][i];
			ty=config.xyz[2][i];
			tz=config.xyz[1][i];
			config.xyz[0][i]=ca*tx+sa*tz;
			config.xyz[1][i]=-sa*cb*tx+sb*ty+ca*cb*tz;
			config.xyz[2][i]=sa*sb*tx+cb*ty-ca*sb*tz;
		    }
		for (int i=0;i<config.nvrt;i++)
		    {
			tx=config.vrt[0][i];
			ty=config.vrt[2][i];
			tz=config.vrt[1][i];
			config.vrt[0][i]=ca*tx+sa*tz;
			config.vrt[1][i]=-sa*cb*tx+sb*ty+ca*cb*tz;
			config.vrt[2][i]=sa*sb*tx+cb*ty-ca*sb*tz;
		    }
		redraw();
	    }
    }

    void identifyAtom()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;
	double aaa,xx1,yy1,zz1,dd1;
	
	if(oper == 0)
	    {

		k=getAtom();

		if(k < 0)
		    {
			mark0=-1;
			mark1=-1;
			mark2=-1;
			repaint();
			return;
		    }

		// Identify clicked atoms

		if(mark0==k)
		    {
			if(edit)
			    {
				// Swap atom identities if different

				if(!atom.zsym.equals(config.atoms[k].zsym))
				    {
					config.atoms[k]=new Element(atom.zsym);
					monitor.println("Atom substitution:");
					monitor.println("ATOM: "+BML.fmt(k+1,10)+"     "+BML.fmt(config.atoms[k].zsym,8));
				    }
			    }
			else
			    {
				mark0=-1;
				mark1=-1;
				mark2=-1;
			    }
			redraw();
		    }
		else
		    {
			mark2=mark1;
			mark1=mark0;
			mark0=k;
			if(mark1 < 0)
			    monitor.println("ATOM: "+BML.fmt(k+1,10)+"     "+BML.fmt(config.atoms[k].zsym,8));
			else
			    {
				xx1=config.xyz[0][mark1]-config.xyz[0][mark0];
				yy1=config.xyz[1][mark1]-config.xyz[1][mark0];
				zz1=config.xyz[2][mark1]-config.xyz[2][mark0];
				dd1=Math.sqrt(Math.pow(xx1,2)+Math.pow(yy1,2)+Math.pow(zz1,2));
				if(mark2 < 0)
				    {
					monitor.println("ATOM: "+BML.fmt(k+1,10)+"     "+
							BML.fmt(config.atoms[k].zsym,8)+
							"     "+BML.fmt(dd1,10));
				    }
				else
				    {
					aaa=-(xx1*xx2+yy1*yy2+zz1*zz2)/(dd1*dd2);
					if(aaa > 1)aaa=1;
					if(aaa < -1)aaa=-1;
					aaa=(180/Math.PI)*Math.acos(aaa);
					monitor.println("ATOM: "+BML.fmt(k+1,10)+"     "+
							   BML.fmt(config.atoms[k].zsym,8)+
							   "     "+BML.fmt(dd1,10)+"     "+
							   BML.fmt(aaa,10));
				    }
				xx2=xx1;
				yy2=yy1;
				zz2=zz1;
				dd2=dd1;
			    }
			redraw();
		    }
	    }
    }

    int getAtom()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;
	double ddd,rscale;

	// Get identity of clicked atom

	k=-1;
	if(config != null)
	    {
		rscale=1.0/scale;
		for(int i=0;i<config.natms;i++)
		    {
			ddd=Math.sqrt(Math.pow(config.xyz[0][i]-sx*rscale,2)+
				      Math.pow(config.xyz[2][i]-sy*rscale,2));
			if(ddd < fac*config.atoms[i].zrad)
			    {
				if(k < 0)
				    k=i;
				else if(config.xyz[1][k] > config.xyz[1][i])
				    k=i;
			    }
		    }
	    }
	return k;
    }

    int cgropt(int keyopt,int natms,double step,int ida[],double hnrm[],
	       double fff[],double grad[],double hhh[][],double dxyz[][])
    {
	/*
*********************************************************************

dl-poly/java GUI conjugate gradient routine

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int j;
	double ggg,stride,gam2;

	// Magnitude of current gradient vector

	ggg=0.0;
	for(int i=0;i<natms;i++)
	    {
		ggg+=(dxyz[0][i]*dxyz[0][i]+dxyz[1][i]*dxyz[1][i]+
		      dxyz[2][i]*dxyz[2][i]);
	    }
	ggg=Math.sqrt(ggg);

	if(keyopt == 0)
	    {
		// Set original search direction (vector hhh)

		keyopt=1;
		hnrm[0]=ggg;
		grad[0]=ggg;
		grad[2]=ggg;
		fff[2]=fff[0];

		for(int i=0;i<natms;i++)	    
		    {
			j=ida[i];
			hhh[0][i]=dxyz[0][i];
			hhh[1][i]=dxyz[1][i];
			hhh[2][i]=dxyz[2][i];
			config.xyz[0][j]+=(step*hhh[0][i]);
			config.xyz[1][j]+=(step*hhh[1][i]);
			config.xyz[2][j]+=(step*hhh[2][i]);
		    }
	    }
	else if(keyopt == 1)
	    {
		// Line search along chosen direction

		stride=step;
		fff[1]=fff[2];
		fff[2]=fff[0];
		grad[1]=grad[2];

		grad[2]=0.0;
		for(int i=0;i<natms;i++)
		    {
			grad[2]+=(hhh[0][i]*dxyz[0][i]+hhh[1][i]*dxyz[1][i]+
				  hhh[2][i]*dxyz[2][i]);
		    }
		grad[2]=grad[2]/hnrm[0];

		    // Linear extrapolation to minimum

		if(grad[2] < 0)
		    {
			stride=step*grad[2]/(grad[1]-grad[2]);
			keyopt=2;
		    }
		for(int i=0;i<natms;i++)        
		    {
			j=ida[i];
			config.xyz[0][j]+=(stride*hhh[0][i]);
			config.xyz[1][j]+=(stride*hhh[1][i]);
			config.xyz[2][j]+=(stride*hhh[2][i]);
		    }
	    }
	else if(keyopt == 2)
	    {
		fff[1]=fff[2];
		fff[2]=fff[0];

		// Check for global convergence

		if(Math.abs(ggg/natms) < 0.0000001)
		    {
			return 999;
		    }

		// Construct conjugate search vector
        
		gam2=Math.pow((ggg/grad[0]),2);

		grad[0]=ggg;
		hnrm[0]=0.0;
		grad[2]=0.0;
		for(int i=0;i<natms;i++)
		    {
			hhh[0][i]=dxyz[0][i]+gam2*hhh[0][i];
			hhh[1][i]=dxyz[1][i]+gam2*hhh[1][i];
			hhh[2][i]=dxyz[2][i]+gam2*hhh[2][i];
			hnrm[0]+=(hhh[0][i]*hhh[0][i]+hhh[1][i]*hhh[1][i]+
				  hhh[2][i]*hhh[2][i]);
			grad[2]+=(hhh[0][i]*dxyz[0][i]+hhh[1][i]*dxyz[1][i]+
				  hhh[2][i]*dxyz[2][i]);
		    }
		hnrm[0]=Math.sqrt(hnrm[0]);
		grad[2]/=hnrm[0];
		for(int i=0;i<natms;i++)
		    {
			j=ida[i];
			config.xyz[0][j]+=(step*hhh[0][i]);
			config.xyz[1][j]+=(step*hhh[1][i]);
			config.xyz[2][j]+=(step*hhh[2][i]);
		    }
		keyopt=1;
	    }
	return keyopt;
    }

    void moleculeBuilder()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int j,k,m,mnth;
	double ddd,rscale;

	rscale=1.0/scale;

	if(oper == 1)
	    {

		// set up initial arrays
		
		if(config == null)
		    {
			link=false;
			config=new Config();
			Calendar datum=Calendar.getInstance();
			mnth=datum.get(Calendar.MONTH)+1;
			config.title="GUI Reference : "+
			    datum.get(Calendar.YEAR)+"/"+mnth+
			    "/"+datum.get(Calendar.DAY_OF_MONTH)+"/"
			    +datum.get(Calendar.HOUR_OF_DAY)+":"
			    +datum.get(Calendar.MINUTE);
		    }

		// Build molecular structure

	    OUT:
		if(lpen)
		    {
			// Check identity of new atom

			k=-1;
			for(int i=0;i<config.natms;i++)
			    {
				ddd=Math.sqrt(Math.pow(config.xyz[0][i]-sx*rscale,2)+
					      Math.pow(config.xyz[2][i]-sy*rscale,2));
				if(ddd < fac*config.atoms[i].zrad) 
				    {
					if(k == -1)
					    k=i+1;
					else if(config.xyz[1][i] < config.xyz[1][k])
					    k=i+1;
				    }
			    }

			// Same atom as last time - finish current fragment

			if(k == config.natms) 
			    {
				link=false;
				lpen=false;
				break OUT;
			    }

			if(k < 0)k=config.natms;

			// Add bond to structure

			if(link)
			    {
				j=Math.min(config.natms,k)-1;
				if(k == config.natms)
				    {

				// New atom - bonds to previous entered atom
				    
					m=config.natms;
					config.bond[config.lbnd[m]][m]=j;
					config.bond[config.lbnd[j]][j]=m;
					config.lbnd[m]++;
					config.lbnd[j]++;
					config.join[0][config.nbnds]=Math.min(j,m);
					config.join[1][config.nbnds]=Math.max(j,m);
					config.nbnds++;
				    }
				else
				    {
				    
					// Old atom - connect to existing fragment
				    
					m=config.natms-1;
					config.bond[config.lbnd[m]][m]=j;
					config.bond[config.lbnd[j]][j]=m;
					config.lbnd[m]++;
					config.lbnd[j]++;
					config.join[0][config.nbnds]=Math.min(j,m);
					config.join[1][config.nbnds]=Math.max(j,m);
					config.nbnds++;
				    }
			    }
			if(k == config.natms)
			    {
				// Add new atom to fragment
			    
				config.atoms[config.natms]=new Element(atom.zsym);
				config.xyz[0][config.natms]=sx*rscale;
				config.xyz[2][config.natms]=sy*rscale;
				config.xyz[1][config.natms]=0.0;
				link=true;
				safe=false;
				config.natms++;
			    }
			else
			    {

				// Old atom - finish current fragment

				link=false;
				lpen=false;
				config.xyz[1][config.natms-1]=config.xyz[1][k-1];
				break OUT;
			    }
		    }
		else
		    lpen=true;

		// extend arrays if atom count equals current array sizes

		if(config.natms == config.mxatms) 
		    getAtomicArrays();
		if(config.nbnds == config.mxjoin) 
		    getJoinArray();
		redraw();
	    }
    }

    void linkMaker()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	if(oper == 2)
	    {
		if(mark0 < 0)
		    {
			mark0=getAtom();
			repaint();
		    }
		else if(mark1 < 0)
		    {
			mark1=getAtom();
			if(mark0 == mark1)
			    {
				mark1=-1;
				repaint();
			    }
			else if(mark1 >= 0)
			    {
				if(!duplicateBond())
				    {
					config.bond[config.lbnd[mark0]][mark0]=mark1;
					config.bond[config.lbnd[mark1]][mark1]=mark0;
					config.lbnd[mark0]++;
					config.lbnd[mark1]++;
					config.join[0][config.nbnds]=Math.min(mark0,mark1);
					config.join[1][config.nbnds]=Math.max(mark0,mark1);
					config.nbnds++;
					mark0=-1;
					mark1=-1;
					repaint();
				    }
			    }
		    }
		if(config.nbnds == config.mxjoin)
		    getJoinArray();
	    }
    }

    void makeGroup()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;

	k=getAtom();
	if(k < 0) return;

	if(ngroup == 0)
	    {
		mxgroup=config.mxatms;
		group=new int[mxgroup];
		for(int i=0;i<mxgroup;i++)
		    group[i]=-1;
		group[k]=1;
		ngroup++;
	    }
	else if(group[k] > 0)
	    {
		group=null;
		ngroup=0;
	    }
	else
	    {
		group[k]=1;
		ngroup++;
	    }
	repaint();
    }

    void frameMolecule()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	double rscale;
	int bx,ex,by,ey;

	sx=sx-tilex/2;
	sy=tiley/2-sy;
	fx=fx-tilex/2;
	fy=tiley/2-fy;
	rscale=1.0/scale;
	bx=Math.min(sx,fx);
	ex=Math.max(sx,fx);
	by=Math.min(sy,fy);
	ey=Math.max(sy,fy);

	if(ngroup == 0)
	    {
		mxgroup=config.mxatms;
		group=new int[mxgroup];
		for(int i=0;i<config.natms;i++)
		    group[i]=-1;
	    }

	for(int i=0;i<config.natms;i++)
	    {
		if((config.xyz[0][i] >= bx*rscale) && (config.xyz[0][i] <= ex*rscale) &&
		   (config.xyz[2][i] >= by*rscale) && (config.xyz[2][i] <= ey*rscale))
		    {
			group[i]=1;
			ngroup++;
		    }
	    }

	repaint();
    }

    void addDeleteHydrogen()
    {	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k,n;

	if(oper == 10)
	    {
		if(ngroup == 0)
		    {
			// Get identity of clicked atom
			    
			k=getAtom();
			if(k >= 0) 
			    {
				hydrogenAtom(k);
			    }
			else
			    {
				nhdel=0;
				n=config.natms;
				for(int i=0;i<n;i++)
				    if(config.atoms[i].znum == 1) nhdel++;
				if(nhdel > 0)
				    {
					delHAtoms();
				    }
				else
				    {
					for(int i=0;i<n;i++)
					    {
						hydrogenAtom(i);
					    }
				    }
			    }
		    }
		else
		    {
			nhdel=0;
			n=config.natms;
			for(int i=0;i<n;i++)
			    if(group[i] > 0 && config.atoms[i].znum == 1) nhdel++;
			if(nhdel > 0)
			    {
				delHAtoms();
			    }
			else
			    {
				for (int i=0;i<n;i++)
				    {
					if(group[i] > 0)
					    {
						hydrogenAtom(i);
					    }
				    }
			    }
			ngroup=0;
		    }
		safe=false;
		redraw();
	    }
    }

    void hydrogenAtom(int k)
    {	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	String str;
	int nhy,hyb,val;
	double stick;

	nhy=0;
	hyb=0;
	val=0;
	str=config.atoms[k].zsym;
	stick=config.atoms[k].zrad+(new Element("H_")).zrad;

	if(str.indexOf("C_1") >= 0)
	    {
		hyb=1;
		val=2;
		nhy=2-config.lbnd[k];
	    }
	else if (str.indexOf("C_2") >= 0)
	    {
		hyb=2;
		val=3;
		nhy=3-config.lbnd[k];
	    }
	else if (str.indexOf("C_R") >= 0)
	    {
		hyb=2;
		val=3;
		nhy=3-config.lbnd[k];
	    }
	else if (str.indexOf("C_3") >= 0)
	    {
		hyb=3;
		val=4;
		nhy=4-config.lbnd[k];
	    }
	else if (str.indexOf("O_3") >= 0)
	    {
		hyb=3;
		val=2;
		nhy=2-config.lbnd[k];
	    }
	else if (str.indexOf("OW") >= 0)
	    {
		hyb=3;
		val=2;
		nhy=2-config.lbnd[k];
	    }
	else if (str.indexOf("N_2") >= 0)
	    {
		hyb=2;
		val=2;
		nhy=2-config.lbnd[k];
	    }
	else if (str.indexOf("N_3") >= 0)
	    {
		hyb=3;
		val=3;
		nhy=3-config.lbnd[k];
	    }
	else if (str.indexOf("S_3") >= 0)
	    {
		hyb=3;
		val=2;
		nhy=2-config.lbnd[k];
	    }
	else if (str.indexOf("P_2") >= 0)
	    {
		hyb=2;
		val=3;
		nhy=3-config.lbnd[k];
	    }
	else if (str.indexOf("P_3") >= 0)
	    {
		hyb=3;
		val=4;
		nhy=4-config.lbnd[k];
	    }
	if(nhy > 0)
	    addHAtom(k,nhy,hyb,val,stick);
    }

    void addHAtom(int k,int nhy,int hyb,int val,double stick)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int a,b,c;
	double[] vv,va,vb,vc,vd;
	double xx,yy,zz,ct,st,cp,sp,atn;

	if(hyb == 1 && nhy == 1 && val == 2)
	    {
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		addVH(k,stick,va);
	    }
	else if(hyb == 2 && nhy == 1 && val == 2)
	    {
		vv=new double[3];
		xx=0.5;
		yy=0.86602540375;
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		atn=Math.atan2(va[1],va[0]);
		cp=Math.cos(atn);
		sp=Math.sin(atn);
		vv[0]=va[0]*xx-sp*yy;
		vv[1]=va[1]*xx+cp*yy;
		vv[2]=va[2]*xx;
		addVH(k,stick,vv);
	    }
	else if(hyb == 2 && nhy == 1 && val == 3)
	    {
		vv=new double[3];
		a=config.bond[0][k];
		b=config.bond[1][k];
		va=BML.ndxyz(k,a,config.xyz);
		vb=BML.ndxyz(k,b,config.xyz);
		vv[0]=va[0]+vb[0];
		vv[1]=va[1]+vb[1];
		vv[2]=va[2]+vb[2];
		BML.vnorm(vv);
		addVH(k,stick,vv);
	    }
	else if(hyb == 2 && nhy == 2 && val == 2)
	    {
		vv=new double[3];
		vv[0]=0.5;
		vv[1]=0.86602540384;
		vv[2]=0.0;
		addVH(k,stick,vv);
		vv[1]=-0.86602540384;
		addVH(k,stick,vv);
	    }
	else if(hyb == 2 && nhy == 2 && val == 3)
	    {
		vv=new double[3];
		xx=0.5;
		yy=0.86602540375;
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		atn=Math.atan2(va[1],va[0]);
		cp=Math.cos(atn);
		sp=Math.sin(atn);
		vv[0]=va[0]*xx-sp*yy;
		vv[1]=va[1]*xx+cp*yy;
		vv[2]=va[2]*xx;
		addVH(k,stick,vv);
		vv[0]=va[0]*xx+sp*yy;
		vv[1]=va[1]*xx-cp*yy;
		vv[2]=va[2]*xx;
		addVH(k,stick,vv);
	    }
	else if(hyb == 2 && nhy == 3 && val == 3)
	    {
		vv=new double[3];
		vv[0]=1.0;
		vv[1]=0.0;
		vv[2]=0.0;
		addVH(k,stick,vv);
		vv[0]=0.5;
		vv[1]=0.86602540384;
		addVH(k,stick,vv);
		vv[2]=-0.86602540384;
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 1 && val == 2)
	    {
		vv=new double[3];
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		atn=Math.atan2(va[1],va[0]);
		cp=Math.cos(atn);
		sp=Math.sin(atn);
		xx=1.0/3.0;
		yy=0.94280904165;
		vv[0]=xx*va[0]-yy*sp;
		vv[1]=xx*va[1]+yy*cp;
		vv[2]=xx*va[2];
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 1 && val == 3)
	    {
		vv=new double[3];
		vc=new double[3];
		a=config.bond[0][k];
		b=config.bond[1][k];
		va=BML.ndxyz(k,a,config.xyz);
		vb=BML.ndxyz(k,b,config.xyz);
		vc[0]=va[0]+vb[0];
		vc[1]=va[1]+vb[1];
		vc[2]=va[2]+vb[2];
		BML.vnorm(vc);
		vd=BML.cross(va,vb);
		BML.vnorm(vd);
		xx=0.57735026924;
		yy=0.81649658092;
		vv[0]=xx*vc[0]+yy*vd[0];
		vv[1]=xx*vc[1]+yy*vd[1];
		vv[2]=xx*vc[2]+yy*vd[2];
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 1 &&  val == 4)
	    {
		vv=new double[3];
		a=config.bond[0][k];
		b=config.bond[1][k];
		c=config.bond[2][k];
		va=BML.ndxyz(k,a,config.xyz);
		vb=BML.ndxyz(k,b,config.xyz);
		vc=BML.ndxyz(k,c,config.xyz);
		vv[0]=va[0]+vb[0]+vc[0];
		vv[1]=va[1]+vb[1]+vc[1];
		vv[2]=va[2]+vb[2]+vc[2];
		BML.vnorm(vv);
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 2 && val == 2)
	    {
		vv=new double[3];
		vv[0]=1.0/3.0;
		vv[1]=0.81649658092;
                vv[2]=-0.47140452082;
		addVH(k,stick,vv);
		vv[1]=-0.81649658092;
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 2 && val == 3)
	    {
		vv=new double[3];
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		atn=Math.atan2(va[1],va[0]);
		cp=Math.cos(atn);
		sp=Math.sin(atn);
		st=va[2];
		ct=Math.sqrt(va[0]*va[0]+va[1]*va[1]);
		xx=1.0/3.0;
		yy=0.81649658092;
                zz=-0.47140452082;
		vv[0]=xx*va[0]-yy*sp-zz*cp*st;
		vv[1]=xx*va[1]+yy*cp-zz*sp*st;
		vv[2]=xx*st+zz*ct;
		addVH(k,stick,vv);
		vv[0]=xx*va[0]+yy*sp-zz*cp*st;
		vv[1]=xx*va[1]-yy*cp-zz*sp*st;
		vv[2]=xx*st+zz*ct;
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 2 && val == 4)
	    {
		vv=new double[3];
		vc=new double[3];
		a=config.bond[0][k];
		b=config.bond[1][k];
		va=BML.ndxyz(k,a,config.xyz);
		vb=BML.ndxyz(k,b,config.xyz);
		vc[0]=va[0]+vb[0];
		vc[1]=va[1]+vb[1];
		vc[2]=va[2]+vb[2];
		BML.vnorm(vc);
		vd=BML.cross(va,vb);
		BML.vnorm(vd);
		xx=0.57735026924;
		yy=0.81649658092;
		vv[0]=xx*vc[0]+yy*vd[0];
		vv[1]=xx*vc[1]+yy*vd[1];
		vv[2]=xx*vc[2]+yy*vd[2];
		addVH(k,stick,vv);
		vv[0]=xx*vc[0]-yy*vd[0];
		vv[1]=xx*vc[1]-yy*vd[1];
		vv[2]=xx*vc[2]-yy*vd[2];
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 3 && val == 3)
	    {
		vv=new double[3];
		vv[0]=1.0/3.0;
		vv[1]=0.0;
		vv[2]=0.94280904165;
		addVH(k,stick,vv);
		vv[2]=-0.47140452076;
		vv[1]=0.81649658092;
		addVH(k,stick,vv);
		vv[1]=-0.81649658092;
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 3 && val == 4)
	    {
		vv=new double[3];
		a=config.bond[0][k];
		va=BML.ndxyz(k,a,config.xyz);
		atn=Math.atan2(va[1],va[0]);
		cp=Math.cos(atn);
		sp=Math.sin(atn);
		ct=va[2];
		st=Math.sqrt(va[0]*va[0]+va[1]*va[1]);
		xx=0.94280904165;
		yy=0.0;
		zz=1.0/3.0;
		vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
		vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
		vv[2]=ct*zz-st*xx;
		addVH(k,stick,vv);
		xx=-0.47140452076;
		yy=0.81649658092;
		vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
		vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
		vv[2]=ct*zz-st*xx;
		addVH(k,stick,vv);
		yy=-0.81649658092;
		vv[0]=ct*cp*xx-sp*yy+cp*st*zz;
		vv[1]=ct*sp*xx+cp*yy+sp*st*zz;
		vv[2]=ct*zz-st*xx;
		addVH(k,stick,vv);
	    }
	else if(hyb == 3 && nhy == 4 && val == 4)
	    {
		vv=new double[3];
		vv[0]=0.57735026918;
		vv[1]=0.57735026918;
		vv[2]=0.57735026918;
		addVH(k,stick,vv);
		vv[0]=-0.57735026918;
		vv[1]=-0.57735026918;
		vv[2]=0.57735026918;
		addVH(k,stick,vv);
		vv[0]=0.57735026918;
		vv[1]=-0.57735026918;
		vv[2]=-0.57735026918;
		addVH(k,stick,vv);
		vv[0]=-0.57735026918;
		vv[1]=0.57735026918;
		vv[2]=-0.57735026918;
		addVH(k,stick,vv);
	    }
    }

    void addVH(int k,double aaa,double vvv[])
    {
	/*
*********************************************************************

dl_poly/java GUI routine to add hydrogen atom and bond vector

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	// create bond arrays if nonexistent

	if(config.nbnds == 0)
	    {
		int[] lll=new int[config.mxatms];
		int[][] bbb=new int[config.mxcon][config.mxatms];
		config.bond=bbb;
		config.lbnd=lll;
	    }

	if (config.atoms[k].zsym.indexOf("OW") >=0 )
	    config.atoms[config.natms]=new Element("HW");
	else
	    config.atoms[config.natms]=new Element("H_");
	config.xyz[0][config.natms]=config.xyz[0][k]+aaa*vvv[0];
	config.xyz[1][config.natms]=config.xyz[1][k]+aaa*vvv[1];
	config.xyz[2][config.natms]=config.xyz[2][k]+aaa*vvv[2];
	config.join[0][config.nbnds]=Math.min(k,config.natms);
	config.join[1][config.nbnds]=Math.max(k,config.natms);
	config.bond[config.lbnd[k]][k]=config.natms;
	config.bond[0][config.natms]=k;
	config.lbnd[config.natms]=1;
	config.lbnd[k]++;
	config.nbnds++;
	config.natms++;

	// extend arrays if atom count equals current array sizes

	if(config.natms == config.mxatms)
	    getAtomicArrays();
	if(config.nbnds == config.mxjoin)
	    getJoinArray();
    }

    void delHAtoms()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k=0;

	if(ngroup == 0)
	    {
		for (int i=0;i<nhdel;i++)
		    {
		    OUT:
			for(int j=k;j<config.natms;j++)
			    {
				if(config.atoms[j].znum == 1)
				    {
					deleteAtom(j);
					k=j;
					break OUT;
				    }
			    }
		    }
	    }
	else
	    {
		for (int i=0;i<nhdel;i++)
		    {
		    OUT:
			for(int j=k;j<config.natms;j++)
			    {
				if(group[j] > 0 && config.atoms[j].znum == 1)
				    {
					deleteAtom(j);
					k=j;
					break OUT;
				    }
			    }
		    }
	    }
	safe=false;
	redraw();
    }

    void getAtomicArrays()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	config.mxatms*=2;
	Element[] aaa=new Element[config.mxatms];
	double[][] ccc=new double[3][config.mxatms];

	for(int i=0;i<config.natms;i++)
	    {
		aaa[i]=new Element(config.atoms[i].zsym);
		ccc[0][i]=config.xyz[0][i];
		ccc[1][i]=config.xyz[1][i];
		ccc[2][i]=config.xyz[2][i];
	    }
	config.atoms=aaa;
	config.xyz=ccc;

	if(config.nbnds >0)
	    {
		int[] lll=new int[config.mxatms];
		int[][] bbb=new int[config.mxcon][config.mxatms];
		
		for(int i=0;i<config.natms;i++)
		    {
			lll[i]=config.lbnd[i];
			for(int n=0;n<config.lbnd[i];n++)
			    bbb[n][i]=config.bond[n][i];
		    }
		for(int i=config.natms;i<config.mxatms;i++)
		    {
			lll[i]=0;
		    }
		config.bond=bbb;
		config.lbnd=lll;
	    }
    }

    void getJoinArray()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	config.mxjoin*=2;
	boolean[] aaa=new boolean[config.mxjoin];
	int[][] bbb=new int[2][config.mxjoin];
	for(int i=0;i<config.nbnds;i++)
	    {
		bbb[0][i]=config.join[0][i];
		bbb[1][i]=config.join[1][i];
	    }
	config.join=bbb;
    }

    void getGroupArray()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2003

*********************************************************************
*/
	mxgroup*=2;
	int[] bbb=new int[mxgroup];
	for(int i=0;i<config.natms;i++)
	    {
		bbb[i]=group[i];
	    }
	group=bbb;
    }

    boolean duplicateBond()
    {
	/*
*********************************************************************

dl_poly/java GUI routine to remove duplicated bonds

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int n,m,k;
	    
	k=-1;
	n=Math.min(mark0,mark1);
	m=Math.max(mark0,mark1);

    OUT:
	for(int i=0;i<config.nbnds;i++)
	    {
		if(n == config.join[0][i] && m == config.join[1][i])
		    {
			k=i;
			config.join[0][i]=-1;
			break OUT;
		    }
	    }
	if(k < 0) return false;

	// Delete duplicated bond
	    
	config.nbnds--;
	for(int i=k;i<config.nbnds;i++)
	    {
		config.join[0][i]=config.join[0][i+1];
		config.join[1][i]=config.join[1][i+1];
	    }
	mark0=-1;
	mark1=-1;

	// Rebuild connection table
	    
	rebuildConnections();
	redraw();

	return true;
    }

    void deleteAtom(int k)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int j,m,n;

	// Remove atom
	    
	config.natms--;
	safe=false;
	for(int i=k;i<config.natms;i++)
	    {
		config.atoms[i]=new Element(config.atoms[i+1].zsym);
		config.xyz[0][i]=config.xyz[0][i+1];
		config.xyz[1][i]=config.xyz[1][i+1];
		config.xyz[2][i]=config.xyz[2][i+1];
		if(ngroup > 0)group[i]=group[i+1];
	    }
	    
	// Mark bonds for removal
	    
	for(int i=0;i<config.nbnds;i++)
	    {
		if(config.join[0][i] == k || config.join[1][i] == k)
		    {
			config.join[0][i]=-1;
		    }
	    }
	// Delete cancelled bonds
	    
	j=0;
	for(int i=0;i<config.nbnds;i++)
	    {
		if(config.join[0][i] >= 0)
		    {
			m=config.join[0][i];
			n=config.join[1][i];
			if(m > k) m--;
			if(n > k) n--;
			config.join[0][j]=m;
			config.join[1][j]=n;
			j++;
		    }
	    }
	config.nbnds=j;
	    
	// Rebuild connection table
	    
	rebuildConnections();
    }

    void deleteGroup()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;
	int n=0;

	if(oper == 3)
	    {
		if(ngroup == 0)
		    {
			// Get identity of clicked atom
			    
			k=getAtom();
			if(k >= 0) deleteAtom(k);
		    }
		else
		    {
			for (int i=0;i<ngroup;i++)
			    {
			    OUT:
				for(int j=n;j<config.natms;j++)
				    {
					if(group[j] > 0)
					    {
						deleteAtom(j);
						n=j;
						break OUT;
					    }
				    }
			    }
			ngroup=0;
		    }
		safe=false;
		redraw();
	    }
    }

    void duplicateGroup()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k,n,m;
	int iab[];

	if(oper == 11 && ngroup > 0)
	    {
		m=config.natms+ngroup;
		if(m >= mxgroup)
		    getGroupArray();
		if(m >= config.mxatms)
		    getAtomicArrays();
		iab=new int[config.mxatms];
		k=config.natms;
		for (int i=0;i<config.natms;i++)
		    {
			iab[i]=-1;
			if(group[i] > 0)
			    {
				iab[i]=k;
				config.atoms[k]=new Element(config.atoms[i].zsym);
				config.xyz[0][k]=config.xyz[0][i]+incx;
				config.xyz[1][k]=config.xyz[1][i];
				config.xyz[2][k]=config.xyz[2][i]+incz;
				group[k]=1;
				k++;
			    }
		    }
		k=config.nbnds;
		for (int i=0;i<config.nbnds;i++)
		    {
			n=config.join[0][i];
			m=config.join[1][i];
			if(group[n] > 0 && group[m] > 0)
			    {
				n=iab[n];
				m=iab[m];
				if(k == config.mxjoin)
				    getJoinArray();
				config.join[0][k]=n;
				config.join[1][k]=m;
				config.bond[config.lbnd[n]][n]=m;
				config.bond[config.lbnd[m]][m]=n;
				config.lbnd[n]++;
				config.lbnd[m]++;
				k++;
			    }
		    }
		config.nbnds=k;
		for(int i=0;i<config.natms;i++)
		    group[i]=-1;
		config.natms+=ngroup;
		safe=false;
		oper=6;
		news="MOVE";
	    }
	redraw();
    }
    
    void insertFragment()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k,n,m;

	if(fragment != null)
	    {
		if(config == null)
		    config=new Config();
		if(home.config == null) home.config=fragment;
		m=config.natms+fragment.natms;
		if(ngroup == 0)
		    group=new int[config.mxatms];
		if(m >= mxgroup)
		    {
			mxgroup=m;
			getGroupArray();
		    }
		if(m >= config.mxatms)
		    {
			config.mxatms=m;
			getAtomicArrays();
		    }
		k=config.natms;
		for (int i=0;i<config.natms;i++)
		    group[i]=-1;
		for (int i=0;i<fragment.natms;i++)
		    {
			config.atoms[k]=new Element(fragment.atoms[i].zsym);
			config.xyz[0][k]=fragment.xyz[0][i];
			config.xyz[1][k]=fragment.xyz[1][i];
			config.xyz[2][k]=fragment.xyz[2][i];
			group[k]=1;
			k++;
		    }
		config.natms=k;
		k=config.nbnds;
		for (int i=0;i<fragment.nbnds;i++)
		    {
			n=fragment.join[0][i]+config.nbnds;
			m=fragment.join[1][i]+config.nbnds;
			if(k == config.mxjoin)
			    getJoinArray();
			config.join[0][k]=n;
			config.join[1][k]=m;
			config.bond[config.lbnd[n]][n]=m;
			config.bond[config.lbnd[m]][m]=n;
			config.lbnd[n]++;
			config.lbnd[m]++;
			k++;
		    }
		config.nbnds=k;
		ngroup=fragment.natms;
		safe=false;
		oper=6;
		news="MOVE";
	    }
	redraw();
	//restore();
    }
    
    void rebuildConnections()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int j,k;

	for(int i=0;i<config.natms;i++)
	    config.lbnd[i]=0;
		
	if(config.nbnds > 0)
	    {
		for(int i=0;i<config.nbnds;i++)
		    {
			j=config.join[0][i];
			k=config.join[1][i];
			config.bond[config.lbnd[j]][j]=k;
			config.bond[config.lbnd[k]][k]=j;
			config.lbnd[j]++;
			config.lbnd[k]++;
		    }
	    }
    }

    void shiftMolecule(String act)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;
	double rad;

	k=getAtom();

	if(ngroup == 0)
	    {
		if(k >= 0)
		    {
			if(act.equals("Tz+"))
			    config.xyz[2][k]+=0.1;
			else if(act.equals("Tz-"))
			    config.xyz[2][k]-=0.1;
			else if(act.equals("Tx+"))
			    config.xyz[0][k]+=0.1;
			else if(act.equals("Tx-"))
			    config.xyz[0][k]-=0.1;
			else if(act.equals("Ty+"))
			    config.xyz[1][k]+=0.1;
			else if(act.equals("Ty-"))
			    config.xyz[1][k]-=0.1;
		    }
		else
		    {
			for(int i=0;i<config.natms;i++)
			    {
				if(act.equals("Tz+"))
				    config.xyz[2][i]+=0.1;
				else if(act.equals("Tz-"))
				    config.xyz[2][i]-=0.1;
				else if(act.equals("Tx+"))
				    config.xyz[0][i]+=0.1;
				else if(act.equals("Tx-"))
				    config.xyz[0][i]-=0.1;
				else if(act.equals("Ty+"))
				    config.xyz[1][i]+=0.1;
				else if(act.equals("Ty-"))
				    config.xyz[1][i]-=0.1;
			    }
		    }
	    }
	else
	    {		    
		if(k >= 0 && group[k] < 0)
		    {
			if(act.equals("Tz+"))
			    config.xyz[2][k]+=0.1;
			else if(act.equals("Tz-"))
			    config.xyz[2][k]-=0.1;
			else if(act.equals("Tx+"))
			    config.xyz[0][k]+=0.1;
			else if(act.equals("Tx-"))
			    config.xyz[0][k]-=0.1;
			else if(act.equals("Ty+"))
			    config.xyz[1][k]+=0.1;
			else if(act.equals("Ty-"))
			    config.xyz[1][k]-=0.1;
		    }
		else
		    {
			for (int i=0;i<config.natms;i++)
			    {
				if(group[i] > 0)
				    {
					if(act.equals("Tz+"))
					    config.xyz[2][i]+=0.1;
					else if(act.equals("Tz-"))
					    config.xyz[2][i]-=0.1;
					else if(act.equals("Tx+"))
					    config.xyz[0][i]+=0.1;
					else if(act.equals("Tx-"))
					    config.xyz[0][i]-=0.1;
					else if(act.equals("Ty+"))
					    config.xyz[1][i]+=0.1;
					else if(act.equals("Ty-"))
					    config.xyz[1][i]-=0.1;
				    }
			    }
		    }
	    }
	safe=false;
	redraw();
    }

    void moveMolecule()
    {
	/*
*********************************************************************

dl_Poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k;
	double dx,dz,rad;

	dx=(fx-sx)/scale;
	dz=(sy-fy)/scale;

	if(ngroup == 0)
	    {
		sx=sx-tilex/2;
		sy=tiley/2-sy;
		k=getAtom();
		if(k >= 0)
		    {
			config.xyz[0][k]+=dx;
			config.xyz[2][k]+=dz;
		    }
		else
		    {
			for(int i=0;i<config.natms;i++)
			    {
				config.xyz[0][i]+=dx;
				config.xyz[2][i]+=dz;
			    }
		    }
	    }
	else
	    {
		for(int i=0;i<config.natms;i++)
		    {
			if(group[i] > 0)
			    {
				config.xyz[0][i]+=dx;
				config.xyz[2][i]+=dz;
			    }
		    }
	    }
	safe=false;
	redraw();
    }

    void turnMolecule(String act)
    {
	/*
*********************************************************************

dl_Poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int n,j,k;
	double u,v,ca,sa,rad;

	j=-1;
	k=-1;
	sa=sal;
	ca=cal;

	if(act.equals("Rx+"))
	    {
		j=1;
		k=2;
		sa=sal;
	    }
	else if(act.equals("Rx-"))
	    {
		j=1;
		k=2;
		sa=-sal;
	    }
	else if(act.equals("Ry+"))
	    {
		j=2;
		k=0;
		sa=sal;
	    }
	else if(act.equals("Ry-"))
	    {
		j=2;
		k=0;
		sa=-sal;
	    }
	else if(act.equals("Rz+"))
	    {
		j=0;
		k=1;
		sa=sal;
	    }
	else if(act.equals("Rz-"))
	    {
		j=0;
		k=1;
		sa=-sal;
	    }

	if(ngroup == 0)
	    {
		n=getAtom();
		if(n >= 0)
		    {
			u=config.xyz[j][n];
			v=config.xyz[k][n];
			config.xyz[j][n]=(ca*u-sa*v);
			config.xyz[k][n]=(sa*u+ca*v);
		    }
		else
		    {
			for (int i=0;i<config.natms;i++)
			    {
				u=config.xyz[j][i];
				v=config.xyz[k][i];
				config.xyz[j][i]=(ca*u-sa*v);
				config.xyz[k][i]=(sa*u+ca*v);
			    }
		    }
	    }
	else
	    {
		for (int i=0;i<config.natms;i++)
		    {
			if(group[i] > 0)
			    {
				u=config.xyz[j][i];
				v=config.xyz[k][i];
				config.xyz[j][i]=(ca*u-sa*v);
				config.xyz[k][i]=(sa*u+ca*v);
			    }
		    }
	    }
	safe=false;
	redraw();
    }

    void rotateMolecule()
    {
	/*
*********************************************************************

dl_Poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	double[] r={1,0,0,0,1,0,0,0,1};
	double rad,r1,r2,ca,cb,sa,sb,tx,ty,tz;

	if(act.equals("Rx+") || act.equals("Rx-"))
	    {
		ca=Math.cos((fy-sy)*Math.PI/tiley);
		sa=Math.sin((fy-sy)*Math.PI/tiley);
		r[4]= ca;
		r[5]= sa;
		r[7]=-sa;
		r[8]= ca;
	    }
	else if(act.equals("Rz+") || act.equals("Rz-"))
	    {
		cb=Math.cos((fx-sx)*Math.PI/tilex);
		sb=Math.sin((fx-sx)*Math.PI/tilex);
		r[0]= cb;
		r[1]= sb;
		r[3]=-sb;
		r[4]= cb;
	    }
	else if(act.equals("Ry+") || act.equals("Ry-"))
	    {
		if(sx > tilex/2)
		    rad=(fy-sy)*Math.PI/tiley;
		else
		    rad=-(fy-sy)*Math.PI/tiley;
		ca=Math.cos(rad);
		sa=Math.sin(rad);
		r[0]= ca;
		r[2]=-sa;
		r[6]= sa;
		r[8]= ca;
	    }
	else
	    {
		ca=Math.cos((sx-fx)*Math.PI/tilex);
		sa=Math.sin((sx-fx)*Math.PI/tilex);
		cb=Math.cos((sy-fy)*Math.PI/tiley);
		sb=Math.sin((sy-fy)*Math.PI/tiley);
		r[0]= ca;
		r[2]= sa*sb;
		r[1]=-sa*cb;
		r[6]= 0;
		r[8]= cb;
		r[7]= sb;
		r[3]= sa;
		r[5]=-ca*sb;
		r[4]= ca*cb;
	    }
	if(ngroup > 0)
	    {
		for (int i=0;i<config.natms;i++)
		    {
			if(group[i] > 0)
			    {
				tx=config.xyz[0][i];
				ty=config.xyz[1][i];
				tz=config.xyz[2][i];
				config.xyz[0][i]=r[0]*tx+r[3]*ty+r[6]*tz;
				config.xyz[1][i]=r[1]*tx+r[4]*ty+r[7]*tz;
				config.xyz[2][i]=r[2]*tx+r[5]*ty+r[8]*tz;
			    }
		    }
	    }
	else
	    {
		for (int i=0;i<config.natms;i++)
		    {
			tx=config.xyz[0][i];
			ty=config.xyz[1][i];
			tz=config.xyz[2][i];
			config.xyz[0][i]=r[0]*tx+r[3]*ty+r[6]*tz;
			config.xyz[1][i]=r[1]*tx+r[4]*ty+r[7]*tz;
			config.xyz[2][i]=r[2]*tx+r[5]*ty+r[8]*tz;
		    }
	    }
	safe=false;
	redraw();
    }

    void Optimize()
    {
	/*
*********************************************************************

dl_Poly/java GUI routine to optimize structures

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/

	double step;
	int k,keyopt,matms,mbnds,mxangs,mxcon;
	double[] fff,grad,hnrm,length,angcon;
	double[][] hhh,dxyz;
	int[][] jopt,tbond,angle;
	int[] ida,idb,tlbnd;

	if(config.natms < 2) return;

	if(oper == 9)
	    {
		if(ngroup == 0)
		    {
			jopt=config.join;
			tbond=config.bond;
			tlbnd=config.lbnd;
			mbnds=config.nbnds;
			matms=config.natms;
			ida=new int[config.natms];
			idb=new int[config.natms];
			for(int i=0;i<config.natms;i++)
			    {
				ida[i]=i;
				idb[i]=i;
			    }
		    }
		else
		    {
			mbnds=0;
			matms=ngroup;
			mxcon=config.mxcon;
			jopt=new int[2][config.nbnds];
			for (int i=0;i<config.nbnds;i++)
			    {
				if(group[config.join[0][i]] > 0 || group[config.join[1][i]] > 0)
				    {
					jopt[0][mbnds]=config.join[0][i];
					jopt[1][mbnds]=config.join[1][i];
					mbnds++;
				    }
			    }
			k=0;
			ida=new int[matms];
			idb=new int[config.natms];
			tlbnd=new int[matms];
			tbond=new int[mxcon][matms];
			for (int i=0;i<config.natms;i++)
			    {
				idb[i]=-1;
				if(group[i] > 0)
				    {
					ida[k]=i;
					idb[i]=k;
					tlbnd[k]=config.lbnd[i];
					for (int j=0;j<config.lbnd[i];j++)
					    {
						tbond[j][k]=config.bond[j][i];
					    }
					k++;
				    }
			    }
		    }
		keyopt=0;
		step=0.01;
		fff=new double[3];
		grad=new double[3];
		hnrm=new double[1];
		hhh=new double[3][matms];
		dxyz=new double[3][matms];
		length=new double[mbnds];
		mxangs=matms*12;
		angle=new int[3][mxangs];
		angcon=new double[mxangs];

		    
		// Define bond properties
		    
		bondProperties(matms,mbnds,mxangs,ida,jopt,angle,tlbnd,tbond,length,angcon);
		    
		k=0;
		while (k < 1000 && keyopt != 999)
		    {
			fff[0]=energyGradient(matms,mbnds,idb,jopt,tlbnd,tbond,angle,
					      length,angcon,dxyz);
			keyopt=cgropt(keyopt,matms,step,ida,hnrm,fff,grad,hhh,dxyz);
			k++;
		    }
		if(keyopt == 999)
		    news="OPTIMIZATION: CONVERGED";
		else
		    news="OPTIMIZATION: ENERGY ="+BML.fmt(fff[0],12);
		safe=false;
		redraw();
	    }
    }

    void bondProperties(int natms,int nbnds,int mxangs,int ida[],int join[][],int angle[][],
			int lbnd[],int bond[][],double length[],double angcon[])
    {
	/*
*********************************************************************

dl_Poly/java GUI routine to define bond and angle properties

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int n,m;

	// Define target bond lengths

	for(int i=0;i<nbnds;i++)
	    {
		n=join[0][i];
		m=join[1][i];
		length[i]=config.atoms[n].zrad+config.atoms[m].zrad;
	    }

	// Define target bond angles

	n=0;
	for(int i=0;i<natms;i++)
	    {
		if(lbnd[i] > 1)
		    {
			for(int j=1;j<lbnd[i];j++)
			    {
				for(int k=0;k<j;k++)
				    {
					m=ida[i];
					angle[0][n]=bond[j][i];
					angle[1][n]=m;
					angle[2][n]=bond[k][i];
					if(config.atoms[m].zsym.charAt(2) == '3')
					    angcon[n]=-1.0/3.0;
					else if(config.atoms[m].zsym.charAt(2) == '2')
					    angcon[n]=-1.0/2.0;
					else if(config.atoms[m].zsym.charAt(2) == 'R')
					    angcon[n]=-1.0/2.0;
					else if(config.atoms[m].zsym.charAt(2) == '1')
					    angcon[n]=-1.0;
					n++;
					if(n == mxangs)
					    {
						mxangs*=2;
						int aaa[][]=new int[3][mxangs];
						double bbb[]=new double[mxangs];
						for(int o=0;o<n;o++)
						    {
							bbb[o]=angcon[o];
							aaa[0][o]=angle[0][o];
							aaa[1][o]=angle[1][o];
							aaa[2][o]=angle[2][o];
						    }
						angle=aaa;
						angcon=bbb;
					    }
				    }
			    }
		    }
	    }
	nangs=n;
    }

    double energyGradient(int natms,int nbnds,int idb[],int join[][],int lbnd[],int bond[][],
			  int angle[][],double length[],double angcon[],double dxyz[][])
    {
	/*
*********************************************************************

dl_Poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	int k,l,n,m,q,b,c;
	double energy,harm,hcos,dx,dy,dz,dxa,dya,dza,dxb,dyb,dzb;
	double rr,ra,rb,coz,gamma,gxa,gya,gza,gxb,gyb,gzb,rt6,sig;
	double sig2,rsig,rs3,rs6,rsq,gmax;

	harm=1.0;
	hcos=0.5;
	gmax=0.0;
	energy=0.0;

	// Initialise the gradient array

	for(int i=0;i<natms;i++)
	    {
		dxyz[0][i]=0.0;
		dxyz[1][i]=0.0;
		dxyz[2][i]=0.0;
	    }

	// Contribution from harmonic bonds

	for(int i=0;i<nbnds;i++)
	    {
		n=join[0][i];
		m=join[1][i];
		dx=config.xyz[0][m]-config.xyz[0][n];
		dy=config.xyz[1][m]-config.xyz[1][n];
		dz=config.xyz[2][m]-config.xyz[2][n];
		rr=Math.sqrt(dx*dx+dy*dy+dz*dz);
		energy+=0.5*harm*Math.pow((rr-length[i]),2);
		gamma=-harm*(rr-length[i])/rr;
		gmax=Math.max(gmax,Math.abs(gamma));
		m=idb[m];
		if(m >= 0)
		    {
			dxyz[0][m]+=(dx*gamma);
			dxyz[1][m]+=(dy*gamma);
			dxyz[2][m]+=(dz*gamma);
		    }
		n=idb[n];
		if(n >= 0)
		    {
			dxyz[0][n]-=(dx*gamma);
			dxyz[1][n]-=(dy*gamma);
			dxyz[2][n]-=(dz*gamma);
		    }
	    }

	// Contribution from cosine angle potentials

	for(int i=0;i<nangs;i++)
	    {
		k=angle[0][i];
		n=angle[1][i];
		m=angle[2][i];
		dxa=config.xyz[0][k]-config.xyz[0][n];
		dya=config.xyz[1][k]-config.xyz[1][n];
		dza=config.xyz[2][k]-config.xyz[2][n];
		ra=Math.sqrt(dxa*dxa+dya*dya+dza*dza);
		dxa/=ra;
		dya/=ra;
		dza/=ra;
		dxb=config.xyz[0][m]-config.xyz[0][n];
		dyb=config.xyz[1][m]-config.xyz[1][n];
		dzb=config.xyz[2][m]-config.xyz[2][n];
		rb=Math.sqrt(dxb*dxb+dyb*dyb+dzb*dzb);
		dxb/=rb;
		dyb/=rb;
		dzb/=rb;
		coz=dxa*dxb+dya*dyb+dza*dzb;
		energy+=(0.5*hcos*Math.pow((coz-angcon[i]),2));
		gamma=-hcos*(coz-angcon[i]);
		gmax=Math.max(gmax,Math.abs(gamma));
		gxa=gamma*(dxb-dxa*coz)/ra;
		gya=gamma*(dyb-dya*coz)/ra;
		gza=gamma*(dzb-dza*coz)/ra;
		gxb=gamma*(dxa-dxb*coz)/rb;
		gyb=gamma*(dya-dyb*coz)/rb;
		gzb=gamma*(dza-dzb*coz)/rb;
		k=idb[k];
		if(k >= 0)
		    {
			dxyz[0][k]+=gxa;
			dxyz[1][k]+=gya;
			dxyz[2][k]+=gza;
		    }
		n=idb[n];
		if(n >= 0)
		    {
			dxyz[0][n]-=(gxa+gxb);
			dxyz[1][n]-=(gya+gyb);
			dxyz[2][n]-=(gza+gzb);
		    }
		m=idb[m];
		if(m >= 0)
		    {
			dxyz[0][m]+=gxb;
			dxyz[1][m]+=gyb;
			dxyz[2][m]+=gzb;
		    }
	    }
	return energy;
    }

    void defineBox()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	if(config != null)
	    {
		if(config.imcon != boxtyp)
		    {
			config.imcon=boxtyp;
			config.cell[0]=10;
			config.cell[4]=10;
			config.cell[8]=10;
			config.setBox();
			safe=false;
		    }
	    }
    }

    void boxResize()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	double scy=Math.max(0.1,1+(fy-sy)/(double)tilex);
	if(config != null)
	    {
		if(config.imcon == 2 || config.imcon == 7)
		    {
			double scx=Math.max(0.1,1+(fx-sx)/(double)tilex);
			config.cell[0]*=scx;
			config.cell[4]*=scx;
			config.cell[8]*=scy;
			for (int i=0;i<config.nvrt;i++)
			    {
				config.vrt[0][i]*=scx;
				config.vrt[1][i]*=scx;
				config.vrt[2][i]*=scy;
			    }
		    }
		else
		    {
			config.cell[0]*=scy;
			config.cell[4]*=scy;
			config.cell[8]*=scy;
			for (int i=0;i<config.nvrt;i++)
			    {
				config.vrt[0][i]*=scy;
				config.vrt[1][i]*=scy;
				config.vrt[2][i]*=scy;
			    }
		    }
		pen++;
		safe=false;
		if(pen == 5)
		    {
			pen=0;
			sx=fx;
			sy=fy;
			redraw();
		    }
	    }
    }

    class MyScreenSizer implements ComponentListener
    {
	/*
*********************************************************************

dl_poly/java GUI class to handle screen events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	public void componentResized(ComponentEvent e)
	{
	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    Dimension arg = ((Component)e.getSource()).getSize();

	    tilex = arg.width;
	    tiley = arg.height;
	    if(!edit) rescale();
	}
	public void componentMoved(ComponentEvent e){}
	public void componentShown(ComponentEvent e){}
	public void componentHidden(ComponentEvent e){}
    }

    class MousePoints implements MouseListener
    {
	/*
*********************************************************************

dl_poly/java GUI class to handle mouse events

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	public void mouseClicked(MouseEvent e)
	{
	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	    sx=e.getX()-tilex/2;
	    sy=tiley/2-e.getY();

	    if(edit)
		{		    
		    if(oper == 0)
			{
			    identifyAtom();
			}
		    else if(oper == 1)
			{
			    moleculeBuilder();
			}
		    else if(oper == 2)
			{
			    linkMaker();
			}
		    else if(oper == 3)
			{
			    deleteGroup();
			}
		    else if(oper/10 == 4)
			{
			    shiftMolecule(act);
			}
		    else if(oper == 5)
			{
			    makeGroup();
			}
		    else if(oper/10 == 8)
			{
			    turnMolecule(act);
			}
		    else if(oper == 9)
			{
			    Optimize();
			}
		    else if(oper == 10)
			{
			    addDeleteHydrogen();
			}
		    else if(oper == 11)
			{
			    duplicateGroup();
			}
		    else if(oper == 13)
			{
			    insertFragment();
			}
		    else if(oper == 14)
			{
			    makeBonds();
			}
		}
	    else 
		{
		    oper=0;
		    identifyAtom();
		}
	}

	public void mousePressed(MouseEvent e)
	{
	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    sx=e.getX();
	    sy=e.getY();
	    lmove=false;
	}

	public void mouseReleased(MouseEvent e)
	{
	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	    int j;
	    double aaa,dd1,xx1,yy1,zz1,w;
	    
	    fx=e.getX();
	    fy=e.getY();

	    if(sx != fx && sy != fy)
		{
		    if(edit)
			{
			    lmove=false;
			    if(oper == 5)
				{
				    frameMolecule();
				}
			    else if(oper == 6)
				{
				    moveMolecule();
				}
			    else if(oper == 7)
				{
				    rotateMolecule();
				}
			}
		    else
			{
			    if(oper == 6)
				shiftAtoms();
			    else if(oper == 7)
				turnAtoms();
			}
		}
	}
	public void mouseEntered(MouseEvent e){}
	public void mouseExited(MouseEvent e){}
    }
    
    class MouseMoves implements MouseMotionListener
    {
	public void mouseDragged(MouseEvent e)
	{
	    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2002

*********************************************************************
*/
	    fx=e.getX();
	    fy=e.getY();
	    if(edit)
		{
		    if(oper == 5)
			{
			    lmove=true;
			    repaint();
			}
		    else if(oper == 6)
			{
			    moveMolecule();
			    sx=fx;
			    sy=fy;
			}
		    else if(oper == 7)
			{
			    rotateMolecule();
			    sx=fx;
			    sy=fy;
			}
		    else if(oper == 12)
			{
			    boxResize();
			}
		}
	}
	public void mouseMoved(MouseEvent e){}	
    }
    
}
