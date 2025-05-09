import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;

// Define the Graphical User Interface

public class WaterAdd extends Super implements ActionListener
{
    /*
*********************************************************************

dl_poly/java GUI class to add water to a CONFIG file

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    static double slx,sly,slz,top,bot,dws,dww;
    static JButton load,make,close;
    static JCheckBox slab;
    static JTextField dlx,dly,dlz,ubd,lbd,wsd,wwd;
    static WaterAdd job;
    static boolean lslab,lfetch;
    static GUI home=null;
    
    public WaterAdd()
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
        setTitle("Add Water Panel");

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
	fix(new JLabel("        "),grd,gbc,1,0,1,1);
	
	// Define the Make button
	
        make = new JButton("Make");
        make.setBackground(butn);
        make.setForeground(butf);
        fix(make,grd,gbc,2,0,1,1);
	
	// Water-solute minimum distance
	
        JLabel lab1 = new JLabel("Min water-solute distance:",JLabel.LEFT);
        fix(lab1,grd,gbc,0,1,2,1);
        wsd = new JTextField(8);
        wsd.setBackground(scrn);
        wsd.setForeground(scrf);
        fix(wsd,grd,gbc,2,1,1,1);

	// Water-water minimum distance
	
        JLabel lab2 = new JLabel("Min water-water distance:",JLabel.LEFT);
        fix(lab2,grd,gbc,0,2,2,1);
        wwd = new JTextField(8);
        wwd.setBackground(scrn);
        wwd.setForeground(scrf);
        fix(wwd,grd,gbc,2,2,1,1);

	// Water slab option

	slab=new JCheckBox("Slab?");
        slab.setBackground(back);
        slab.setForeground(fore);
	fix(slab,grd,gbc,0,3,1,1);
	
	// Slab direction vector
	
        JLabel lab3 = new JLabel("Slab direction vector:",JLabel.LEFT);
        fix(lab3,grd,gbc,0,4,3,1);
        dlx = new JTextField(8);
        dlx.setBackground(scrn);
        dlx.setForeground(scrf);
        fix(dlx,grd,gbc,0,5,1,1);
        
        dly = new JTextField(8);
        dly.setBackground(scrn);
        dly.setForeground(scrf);
        fix(dly,grd,gbc,1,5,1,1);
        
        dlz = new JTextField(8);
        dlz.setBackground(scrn);
        dlz.setForeground(scrf);
        fix(dlz,grd,gbc,2,5,1,1);
        
	// Upper bound of slab
	
        JLabel lab4 = new JLabel("Upper bound:",JLabel.RIGHT);
        fix(lab4,grd,gbc,0,6,2,1);
        ubd = new JTextField(8);
        ubd.setBackground(scrn);
        ubd.setForeground(scrf);
        fix(ubd,grd,gbc,2,6,1,1);

	// Lower bound of slab
	
        JLabel lab5 = new JLabel("Lower bound:",JLabel.RIGHT);
        fix(lab5,grd,gbc,0,7,2,1);
        lbd = new JTextField(8);
        lbd.setBackground(scrn);
        lbd.setForeground(scrf);
        fix(lbd,grd,gbc,2,7,1,1);
        
	// Define the Close button
	
        close = new JButton("Close");
        close.setBackground(butn);
        close.setForeground(butf);
        fix(close,grd,gbc,0,8,1,1);
	
	// Register action buttons
	
	load.addActionListener(this);
	make.addActionListener(this);
	close.addActionListener(this);
	
    }

    // Constructor method
    
    public WaterAdd(GUI here)
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
	home=here;
	monitor.println("Activated panel for adding water to a CONFIG file");
	
	// Set up Graphical User interface
	
	job = new WaterAdd();
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
	dww=3.0;
	dws=3.0;
        slx=0.0;
	sly=0.0;
	slz=1.0;
	top=3.0;
	bot=-3.0;
	lslab=false;
	lfetch=false;
	slab.setSelected(lslab);
        wsd.setText(String.valueOf(dws));
        wwd.setText(String.valueOf(dww));
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
		lfetch=true;
		lslab=slab.isSelected();
		dws=BML.giveDouble(wsd.getText(),1);
		dww=BML.giveDouble(wwd.getText(),1);
		slx=BML.giveDouble(dlx.getText(),1);
		sly=BML.giveDouble(dly.getText(),1);
		slz=BML.giveDouble(dlz.getText(),1);
		top=BML.giveDouble(ubd.getText(),1);
		bot=BML.giveDouble(lbd.getText(),1);
		call=addwater(lslab,lfetch,dws,dww,slx,sly,slz,bot,top);
	    }
	else if (arg.equals("Make"))
	    {
		lfetch=false;
		if(config == null)lfetch=true;
		lslab=slab.isSelected();
		dws=BML.giveDouble(wsd.getText(),1);
		dww=BML.giveDouble(wwd.getText(),1);
		slx=BML.giveDouble(dlx.getText(),1);
		sly=BML.giveDouble(dly.getText(),1);
		slz=BML.giveDouble(dlz.getText(),1);
		top=BML.giveDouble(ubd.getText(),1);
		bot=BML.giveDouble(lbd.getText(),1);
		call=addwater(lslab,lfetch,dws,dww,slx,sly,slz,bot,top);
	    }
	else if (arg.equals("Close"))
	    {
		job.hide();
	    }
    }
    int addwater (boolean lslab,boolean lfetch,double tolw_s,double tolw_w,
		  double  slx,double sly,double slz,double slb,double slt)
    {
	/*      
***********************************************************************
     
dl_poly/java utility to add SPC water molecules to a structure to fill
out the MD cell.
Assumes atomic positions are in a form compatible
with the CONFIG file used in DL_POLY with periodic boundary
conditions.

Water is added from 'WATER.300K' file

copyright daresbury laboratory
author  of original fortran code - t. forester  feb 1993
adaptation to dl_poly/java gui   - w. smith march 2001

***********************************************************************
    */

	String newfile,title,namstr;
	int i,j,k,m,imax,imcon1,iplus,ird,irem,iwr,iwt,jmax,kmax;
	int levcfg,naa,natm1,ncells,mxwater,status,natms,imcon,newnum;
	int ilx,ily,ilz;
	double det,rrr,tolnce,tolow,rlink;

	status=0;
	mxwater=256;

	double[] cell;
	int[] lct=null;
	double ddc[]=new double[3];
	double axcub[]=new double[3];
	double ox[]=new double[mxwater];
	double oy[]=new double[mxwater];
	double oz[]=new double[mxwater];
	double h1x[]=new double[mxwater];
	double h1y[]=new double[mxwater];
	double h1z[]=new double[mxwater];
	double h2x[]=new double[mxwater];
	double h2y[]=new double[mxwater];
	double h2z[]=new double[mxwater];
	double rcell[]=new double[9];
	double celprp[]=new double[10];
      
	int[] nix={  0, 1, 0, 0,-1, 1, 0,-1, 1, 0,-1, 1,-1, 1,-1, 0, 0, 1,-1, 0, 1,-1, 0, 1,-1, 1,-1};
	int[] niy={  0, 0, 1, 0, 1, 1,-1, 0, 0, 1,-1,-1, 1, 1, 0,-1, 0,-1,-1, 1, 0, 0,-1, 1, 1,-1,-1};
	int[] niz={  0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,-1, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1};

	// tolerance for O / solute distance

	tolnce=tolw_s*tolw_s;

	// tolerance for water-site / water-site distance

	tolow=tolw_w*tolw_w;

	// read the water configuration file

	status=getWaterFile("WATER.300K",mxwater,axcub,ox,oy,oz,h1x,h1y,h1z,
			    h2x,h2y,h2z);
	if(status<0)
	    {
		monitor.println("Error - processing water file");
		return -1;
	    }
	naa=status;

	// read input CONFIG file

	if(lfetch)
	    {
		config=new Config(home,"CFG");
	    }
	if(config==null)return 0;
	title=config.title;
	cell=config.cell;
	dcell(cell,celprp);
	natms=config.natms;
	imcon=config.imcon;
	newnum=BML.nint(1536*celprp[9]/Math.pow(axcub[0],3));
	int lcells[]=new int[newnum];
	int lstrem[]=new int[newnum];
	Element atoms[]=new Element[newnum];
	double xyz[][]=new double[3][newnum];
	double ddd[][]=new double[3][newnum];
	double sss[][]=new double[3][newnum];
	det=invert(cell,rcell);
	imcon1=imcon;
	imcon=Math.min(imcon,3);

	// place solute atoms into linklist

	for(i=0;i<natms;i++)
	    {
		atoms[i]=new Element(config.atoms[i].zsym);
		xyz[0][i]=config.xyz[0][i];
		xyz[1][i]=config.xyz[1][i];
		xyz[2][i]=config.xyz[2][i];
	    }

	// Construct link cells

	rlink=Math.sqrt(tolnce);
	ilx=(int)(celprp[6]/rlink);
	ily=(int)(celprp[7]/rlink);
	ilz=(int)(celprp[8]/rlink);
	ncells=ilx*ily*ilz;
	lct=new int[ncells];
	status=cells(natms,imcon,tolnce,lct,lcells,ddc,cell,rcell,sss,xyz);
	if(status<0)
	    {
		monitor.println("Error - method CELLS failure: first call");
		return -2;
	    }

	// define number of water blocks to use

	imax=(int)((Math.abs(cell[0])+Math.abs(cell[3])+Math.abs(cell[6]))/axcub[0]+.99)/2+1 ;
	jmax=(int)((Math.abs(cell[1])+Math.abs(cell[4])+Math.abs(cell[7]))/axcub[0]+.99)/2+1;
	kmax=(int)((Math.abs(cell[2])+Math.abs(cell[5])+Math.abs(cell[8]))/axcub[0]+.99)/2+1 ;
      
	monitor.println("CONFIG and WATER.300K files read successfully");
	monitor.println("Inserting SPC water into CONFIG file");

	iplus=0;
	status=spcins(atoms,h1x,h1y,h1z,h2x,h2y,h2z,ox,oy,oz,imax,jmax,kmax,
		      naa,natms,newnum,axcub[0],tolnce,lct,lcells,
		      nix,niy,niz,xyz,sss,0.5*cell[8],-0.5*cell[8],
		      ddc,cell,rcell);
	if(status<0)
	    {
		monitor.println("Error - method SPCINS failure");
		return -3;
	    }
	iplus=status;

	// set up list for removals

	natm1=natms+3*iplus;
	for(i=0;i<natm1;i++)
	    lstrem[i]=0;

	// check for bad interactions (less than tolw_w apart)
      
	tolnce=Math.max(tolow,9.0);

	// place all atoms in link cells

	rlink=Math.sqrt(tolnce);
	ilx=(int)(celprp[6]/rlink);
	ily=(int)(celprp[7]/rlink);
	ilz=(int)(celprp[8]/rlink);
	ncells=ilx*ily*ilz;
	lct=new int[ncells];
	status=cells(natm1,imcon,tolnce,lct,lcells,ddc,cell,rcell,sss,xyz);
	if(status<0)
	    {
		monitor.println("Error - method CELLS failure: second call");
		return -4;
	    }
	ncells=status;
	tolnce=tolow;

	// additional checks on non standard periodic boundaries

	if(imcon1>3 && imcon1!=6)
	    {
		imcon=imcon1;
		status=imgchk(imcon,natms,iplus,tolnce,lstrem,sss,xyz,cell,rcell,ddd);
		if(status<0)
		    {
			monitor.println("Error - method IMGCHK failed");
			return -5;
		    }
		iplus=status;
	    }
	else
	    {
		status=spcchk(iplus,ncells,natms,tolnce,lct,lcells,lstrem,
			      nix,niy,niz,sss,xyz,ddc,cell);
		if(status<0)
		    {
			monitor.println("Error - method SPCCHK failed");
			return -6;
		    }
		iplus=status;
	    }

	// reduce water volume to slab if requested

	if(lslab)
	    {

		// normalise the projection vector
      
		rrr=Math.sqrt(slx*slx+sly*sly+slz*slz);
		slx=slx/rrr;
		sly=sly/rrr;
		slz=slz/rrr;

		// projections along slab vector

		for(i=natms;i<natms+3*iplus;i++)
		    {
			rrr=slx*xyz[0][i]+sly*xyz[1][i]+slz*xyz[2][i];
			if(rrr<slb || rrr>slt)
			    {
				irem=natms-1+3*((i-natms)/3);
				lstrem[irem+1]=1;
				lstrem[irem+2]=1;
				lstrem[irem+3]=1;
			    }
		    }

		// remove unwanted water molecules

		iplus=cleanup(natms,iplus,lstrem,xyz);

	    }

	// write out final structure

	newfile="CFGH2O."+String.valueOf(numh2o);
	numh2o++;

	try
	    {
		DataOutputStream outStream=new DataOutputStream(new FileOutputStream(newfile));
		outStream.writeBytes(title+"\n");

		// imcon=10 for ellipsoidal boundary

		if(imcon!=10)
		    {
			outStream.writeBytes(BML.fmt(0,10)+BML.fmt(imcon,10)+"\n");
			outStream.writeBytes(BML.fmt(cell[0],20)+BML.fmt(cell[1],20)+BML.fmt(cell[2],20)+"\n");
			outStream.writeBytes(BML.fmt(cell[3],20)+BML.fmt(cell[4],20)+BML.fmt(cell[5],20)+"\n");
			outStream.writeBytes(BML.fmt(cell[6],20)+BML.fmt(cell[7],20)+BML.fmt(cell[8],20)+"\n");
		    }
		else
		    {
			outStream.writeBytes(BML.fmt(0,10)+BML.fmt(0,10)+"\n");
		    }
		for (k=0;k<natms+3*iplus;k++)
		    {
			outStream.writeBytes(BML.fmt(atoms[k].zsym,8)+BML.fmt(k+1,10)+"\n");
			outStream.writeBytes(BML.fmt(xyz[0][k],20)+BML.fmt(xyz[1][k],20)+BML.fmt(xyz[2][k],20)+"\n");
		    }
		outStream.close();
	    }
	catch(Exception e)
	    {
		monitor.println("Error - writing file: "+newfile);
		return -7;
	    }

	monitor.println("All done: added "+iplus+" water molecules");
	monitor.println(newfile+" file created");

	natms+=3*iplus;

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

	// Draw structure

	pane.restore();

	return 0;
    
    }
    int  cleanup(int natms,int iplus,int lstrem[],double xyz[][])
    {
	/*
*********************************************************************

dl_poly/java routine to remove unwanted water molecules from a CONFIG 
file

copyright daresbury laboratory
author: t. forester feb 1993
java  : w. smith march 2001

*********************************************************************
*/

	int i,j;

	j=natms;

	for(i=natms;i<natms+3*iplus;i++)
	    {
		if(lstrem[i]==0)
		    {
			lstrem[j]=0;
			xyz[0][j]=xyz[0][i];
			xyz[1][j]=xyz[1][i];
			xyz[2][j]=xyz[2][i];
			j++;
		    }
	    }
	iplus=(j-natms)/3;

	return iplus;
    }
    int cells(int natms,int imcon,double tolnce,int lct[],int lcells[],double ddc[],
	      double cell[],double rcell[],double sss[][],double xyz[][])
    {
	/*
**********************************************************************

dl_poly/java routine to construct link cells

copyright daresbury laboratory
author: t. forester feb 1993
java  : w.smith march 2001

**********************************************************************
*/

	boolean linc;
	double rlink,xdc,ydc,zdc;
	int i,icell,ilx,ily,ilz,ix,iy,iz,j,status,ncells;
	double celprp[]=new double[10];

	status=0;
	dcell(cell,celprp);
	if(celprp[9]<1.e-6)
	    {
		monitor.println("Error - system volume too small");
		return -1;
	    }

	rlink=Math.sqrt(tolnce);
	ilx=(int)(celprp[6]/rlink);
	ily=(int)(celprp[7]/rlink);
	ilz=(int)(celprp[8]/rlink);
	ncells=ilx*ily*ilz;

	// check system is big enough
        
	linc =false;
	if(ilx<3) linc=true;
	if(ily<3) linc=true;
	if(ilz<3) linc=true;

	if(linc)
	    {
		monitor.println("Error - system too small for link cells");
		return -3;
	    }

	// calculate link cell indices
        
	for(i=0;i<ncells;i++)
	    lct[i]=-1;

// link-cell cutoff for reduced space
        
	xdc=ddc[0]=(double)ilx;
	ydc=ddc[1]=(double)ily;
	zdc=ddc[2]=(double)ilz;

	// reduced space coordinates
      
	for(i=0;i<natms;i++)
	    {
		sss[0][i]=xyz[0][i];
		sss[1][i]=xyz[1][i];
		sss[2][i]=xyz[2][i];
	    }

	images(natms,imcon,cell,sss);

	for(i=0;i<natms;i++)
	    {
		xyz[0][i]=sss[0][i];
		xyz[1][i]=sss[1][i];
		xyz[2][i]=sss[2][i];
	    }
        
	for(i=0;i<natms;i++)
	    {
		sss[0][i]=(rcell[0]*xyz[0][i]+rcell[3]*xyz[1][i]+rcell[6]*xyz[2][i])+0.5;
		sss[1][i]=(rcell[1]*xyz[0][i]+rcell[4]*xyz[1][i]+rcell[7]*xyz[2][i])+0.5;
		sss[2][i]=(rcell[2]*xyz[0][i]+rcell[5]*xyz[1][i]+rcell[8]*xyz[2][i])+0.5;
	    }

	// link neighbours 
        
	for(i=0;i<natms;i++)
	    {
		ix=(int)(xdc*sss[0][i]);
		iy=(int)(ydc*sss[1][i]);
		iz=(int)(zdc*sss[2][i]);
          
		icell=ix+ilx*(iy+ily*iz);
          
		j=lct[icell];
		lct[icell]=i;
		lcells[i]=j;
	    }

	return  ncells;
    }

    int spcins(Element atoms[],double h1x[],double h1y[],double h1z[],double h2x[],double h2y[],
	       double h2z[],double ox[],double oy[],double oz[],int imax,int jmax,int kmax,
	       int naa,int natms,int newnum,double axcub,double tolnce,int lct[],
	       int lcells[],int nix[],int niy[],int niz[],double xyz[][],double sss[][],
	       double slt,double slb,double ddc[],double cell[],double rcell[])
    {
	/*
**********************************************************************

dl_poly/java routine to insert SPC water molecules into a configuration
file

copyright daresbury laboratory
author: t. forester feb 1993
java  : w. smith march 2001

**********************************************************************
*/
      
	boolean ladd;
	String ow="OW      ";
	String hw="HW      ";
	double cx,cy,cz,rsq,sx,sxd,sy,syd,sz,szd,xd,xi,xo,yd,yi;
	double yo,zd,zi,zo,xdc,ydc,zdc;
	int i,ilx,ily,ilz,iplus,ix,iy,iz,j,jc,jj,jx,jy,jz;
	int k,k1,k2,k3,kc,l,status,nsbcll;

	iplus=0;
	status=0;

	xdc=ddc[0];
	ydc=ddc[1];
	zdc=ddc[2];

	ilx=(int)(xdc+0.5);
	ily=(int)(ydc+0.5);
	ilz=(int)(zdc+0.5);

	for(i=-imax;i<=imax;i++)
	    {
		xi=i*axcub;
		for(j=-jmax;j<=jmax;j++)
		    {
			yi=j*axcub;
			for(k=-kmax;k<=kmax;k++)
			    {
				zi=k*axcub ;
				for(l=0;l<naa;l++)
				    {
					xo=xi+ox[l];
					yo=yi+oy[l];
					zo=zi+oz[l];
					if(zo<slt && zo>slb)
					    {
						// check if water is in basic MD cell
              
						sx=(rcell[0]*xo+rcell[3]*yo+rcell[6]*zo);
						sy=(rcell[1]*xo+rcell[4]*yo+rcell[7]*zo);
						sz=(rcell[2]*xo+rcell[5]*yo+rcell[8]*zo);
              
						if(Math.abs(sx)<=0.5)
						    {
							if(Math.abs(sy)<=0.5)
							    {
								if(Math.abs(sz)<=0.5)
								    {
									// compute link cell index
                
									sx=sx+0.5;
									sy=sy+0.5;
									sz=sz+0.5;

									ix=(int)(xdc*sx);
									iy=(int)(ydc*sy);
									iz=(int)(zdc*sz);
								    
									// flag to add water

									ladd=true;

									// loop over nearby link cells of solute
                
									nsbcll=27;
									for(kc=0;kc<nsbcll;kc++)
									    {
										jx=ix+nix[kc];
										jy=iy+niy[kc];
										jz=iz+niz[kc];

										// minimum image convention for link cells

										cx=lnkimg(jx,ilx);
										cy=lnkimg(jy,ily);
										cz=lnkimg(jz,ilz);
										jx-=(cx*ilx);
										jy-=(cy*ily);
										jz-=(cz*ilz);

										// index of neighbouring cell
              
										jc=jx+ilx*(jy+ily*jz);
										jj=lct[jc];

										// ignore if empty

										while(jj>=0)
										    {
											// distance in real space : minimum image applied
                      
											sxd=sss[0][jj]-sx+cx;
											syd=sss[1][jj]-sy+cy;
											szd=sss[2][jj]-sz+cz;
                      
											xd=cell[0]*sxd+cell[3]*syd+cell[6]*szd;
											yd=cell[1]*sxd+cell[4]*syd+cell[7]*szd;
											zd=cell[2]*sxd+cell[5]*syd+cell[8]*szd;

											rsq=xd*xd+yd*yd+zd*zd;

											// test of distance

											if(rsq<tolnce)
											    ladd=false;

											jj=lcells[jj];
										    }
									    }

									// add O to structure
                
									if(ladd)
									    {
										iplus++;
										k1=natms+(iplus-1)*3;
										if(k1+3>=newnum)
										    {
											monitor.println("Error - too many water molecules added");
											return -1;
										    }
									    
										xyz[0][k1]=xo;
										xyz[1][k1]=yo;
										xyz[2][k1]=zo;
										atoms[k1]=new Element(ow);
                
										k2=k1+1;
										atoms[k2]=new Element(hw);
										xyz[0][k2]=xi+h1x[l];
										xyz[1][k2]=yi+h1y[l];
										xyz[2][k2]=zi+h1z[l];
									    
										k3=k1+2;
										atoms[k3]=new Element(hw);
										xyz[0][k3]=xi+h2x[l];
										xyz[1][k3]=yi+h2y[l];
										xyz[2][k3]=zi+h2z[l];
									    }
								    }
							    }
						    }
					    }
				    }
			    }
		    }
	    }
	return iplus;
    }

    int spcchk(int iplus,int ncells,int natms,double tolnce,int lct[],
	       int lcells[],int lstrem[],int nix[],int niy[],int niz[],
	       double sss[][],double xyz[][],double ddc[],double cell[])
    {
	/*
**********************************************************************

dl_poly/java routine to check contacts between SPC water molecules
and other species

copyright daresbury laboratory
author: t.forester feb 1993
java  : w. smith march 2001

**********************************************************************
*/
	double cx,cy,cz,rsq,sxd,syd,szd,xd,yd,zd,xdc,ydc,zdc;
	int ic,ii,ik,ilx,ily,ilz,irem,ix,iy,iz,jc,jj,jk,jx,jy,jz,kc,nsbcll;

	xdc=ddc[0];
	ydc=ddc[1];
	zdc=ddc[2];

	ilx=(int)(xdc+0.5);
	ily=(int)(ydc+0.5);
	ilz=(int)(zdc+0.5);

	ix=0;
	iy=0;
	iz=0;

	for(ic=0;ic<ncells;ic++)
	    {
		ii=lct[ic];
		if(ii>=0)
		    {
		    
			// secondary loop over subcells
		    
			nsbcll=14;
			for(kc=0;kc<nsbcll;kc++)
			    {
				ii=lct[ic];
			    
				jx=ix+nix[kc];
				jy=iy+niy[kc];
				jz=iz+niz[kc];
			    
				// minimum image convention for link cells
			    
				cx=lnkimg(jx,ilx);
				cy=lnkimg(jy,ily);
				cz=lnkimg(jz,ilz);
				jx-=(cx*ilx);
				jy-=(cy*ily);
				jz-=(cz*ilz);
			    
				// index of neighbouring cell
			    
				jc=jx+ilx*(jy+ily*jz);
				jj=lct[jc];
			    
				// ignore if no water present
			    
				if(jj>=0)
				    {
					while(ii>=natms)
					    {
						// water molecule id for site ii
					    
						ik=(ii-natms)/3;
					    
						if(ic==jc) jj=lcells[ii];
						if(jj>=0)
						    {
							while(jj>=natms)
							    {
								// distance in real space : minimum image applied
							    
								sxd=sss[0][jj]-sss[0][ii]+cx;
								syd=sss[1][jj]-sss[1][ii]+cy;
								szd=sss[2][jj]-sss[2][ii]+cz;
							    
								xd=cell[0]*sxd+cell[3]*syd+cell[6]*szd;
								yd=cell[1]*sxd+cell[4]*syd+cell[7]*szd;
								zd=cell[2]*sxd+cell[5]*syd+cell[8]*szd;
							    
								rsq=xd*xd+yd*yd+zd*zd;
							    
								// test of distance
							    
								if(rsq<tolnce)
								    {
									// make sure sites are not on same water molecule!
								    
									jk=(jj-natms)/3;
								    
									if(ik>=0 && jk>=0 && ik!=jk)
									    {
										irem=natms-1+3*(Math.max(ik,jk));
										lstrem[irem+1]=1;
										lstrem[irem+2]=1;
										lstrem[irem+3]=1;
									    }
								    }
								jj=lcells[jj];
							    }
						    }
						jj=lct[jc];
						ii=lcells[ii];
					    }
				    }
			    }
		    }
		ix++;
		if(ix>=ilx)
		    {
			ix=0;
			iy++;
			if(iy>=ily)
			    {
				iy=0;
				iz++;
			    }
		    }
	    }
    
	// remove unwanted water molecules
    
	iplus=cleanup(natms,iplus,lstrem,xyz);
    
	return iplus;
    }

    int imgchk(int imcon,int natms,int iplus,double tolnce,int lstrem[],
	       double sss[][],double xyz[][],double cell[],double rcell[],
	       double ddd[][])
    {
	/*
**********************************************************************

dlpoly utility subroutine to remove additional spc waters from
unusual (periodic) bondaries: TO, RD, HP ellipsoidal

copyright daresbury laboratory
author: t.forester october 1994
java  : w. smith march 2001

**********************************************************************
*/
	boolean check;
	double rsq,sx,sy,sz;
	int i,ii,ik,irem,j,jk,natm1,status;
      
	status=0;
	natm1=natms+3*iplus;

	// Truncated octahedral (TO)

	if(imcon==4)
	    {
		// check that Ow is in central box
	      
		for(i=natms;i<natm1;i+=3)
		    {
			sx=Math.abs(sss[0][i]-0.5);
			sy=Math.abs(sss[1][i]-0.5);
			sz=Math.abs(sss[2][i]-0.5);

			if(sx+sy+sz>0.75)
			    {
				irem=natms-1+3*((i-natms)/3);
				lstrem[irem+1]=1;
				lstrem[irem+2]=1;
				lstrem[irem+3]=1;
			    }

		    }
	    }
	else if(imcon==5)
	    {
		// Rhombic dodecahedral

		for(i=natms;i<natm1;i+=3)
		    {
			sx=Math.abs(sss[0][i]-0.5);
			sy=Math.abs(sss[1][i]-0.5);
			sz=Math.abs(sss[2][i]-0.5);

			if(sx+sy+2.0*sz>1.0)
			    {
				irem=natms-1+3*((i-natms)/3);
				lstrem[irem+1]=1;
				lstrem[irem+2]=1;
				lstrem[irem+3]=1;
			    }
		    }
	    }
	else if(imcon==7)
	    {
		// Hexagonal Prism

		for(i=natms;i<natm1;i+=3)
		    {
			sx=Math.abs(sss[0][i]-0.5);
			sy=Math.abs(sss[1][i]-0.5);
          
			if(3.0*sx+sy>=1.0)
			    {            
				irem=natms-1+3*((i-natms)/3);
				lstrem[irem+1]=1;
				lstrem[irem+2]=1;
				lstrem[irem+3]=1;
			    }
		    }
	    }
	else if(imcon==10)
	    {
		// ellipsoidal shell

		for(i=natms;i<natm1;i+=3)
		    {
			sx=Math.abs(sss[0][i]-0.5);
			sy=Math.abs(sss[1][i]-0.5);
			sz=Math.abs(sss[2][i]-0.5);

			if(sx*sx+sy*sy+sz*sz>0.25)
			    {
				irem=natms-1+3*((i-natms)/3);
				lstrem[irem+1]=1;
				lstrem[irem+2]=1;
				lstrem[irem+3]=1;
			    }
		    }
	    }

	// remove unwanted waters

	iplus=cleanup(natms,iplus,lstrem,xyz);
	natm1=natms+3*iplus;

	for(i=natms;i<natm1-3;i++)
	    {
		// water molecule number

		ik=(i-natms)/3;
		ii=0;

		for(j=i+1;j<natm1;j++)
		    {
			jk=(j-natms)/3;
		      
			if(jk!=ik)
			    {
				ddd[0][ii]=xyz[0][i]-xyz[0][j];
				ddd[1][ii]=xyz[1][i]-xyz[1][j];
				ddd[2][ii]=xyz[2][i]-xyz[2][j];
				ii++;
			    }
		    }

		// periodic image
	    
		images(ii,imcon,cell,ddd);
        
		check=true;

		for(j=0;j<ii;j++)
		    {
			rsq=ddd[0][j]*ddd[0][j]+ddd[1][j]*ddd[1][j]+ddd[2][j]*ddd[2][j];
			if(rsq<tolnce) check=false;
		    }
		if(!check)
		    {
			irem=natms-1+3*ik;
			lstrem[irem+1]=1;
			lstrem[irem+2]=1;
			lstrem[irem+3]=1;
		    }
	    }

	// remove unwanted waters

	iplus=cleanup(natms,iplus,lstrem,xyz);
	return iplus;
    }
    double lnkimg(int jy,int ily)
    /*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2001

*********************************************************************
*/
    {
	if(jy>=ily) 
	    {
		return 1.0;
	    }
	else if(jy<0)
	    {
		return -1.0;      
	    }
	return 0.0;
    }
    int getWaterFile(String waterFile,int mxwater,double axcub[],double ox[],double oy[],
		     double oz[],double h1x[],double h1y[],double h1z[],double h2x[],
		     double h2y[],double h2z[])
    {
	/*
*********************************************************************

dl_poly/java GUI routine 

copyright - daresbury laboratory
author    - w.smith 2000

*********************************************************************
*/
	int naa;

	naa=0;
	axcub[0]=0.0;
	axcub[1]=0.0;
	axcub[2]=0.0;

	try
	    {
		String record=""; 
		//ClassLoader c = getClass().getClassLoader();
		//InputStream instream = c.getResourceAsStream(waterFile);
		File stuff = new File("../java/"+waterFile);
		FileInputStream instream = new FileInputStream(stuff);
		InputStreamReader isr = new InputStreamReader(instream);
		BufferedReader reader = new BufferedReader(isr);
		monitor.println("Reading file: "+waterFile);
		record=reader.readLine();
		monitor.println("File header record: "+record);
		record=reader.readLine();
		axcub[0]=BML.giveDouble(record,1);
		naa=BML.giveInteger(record,2);
		if(naa>mxwater)
		    {
			monitor.println("Error - too many water molecules in WATER.300K file");
			reader.close();
			return -1;
		    }
		for(int i=0;i<naa;i++)
		    {
			record=reader.readLine();
			record=reader.readLine();
			ox[i]=BML.giveDouble(record,1);
			oy[i]=BML.giveDouble(record,2);
			oz[i]=BML.giveDouble(record,3);
			record=reader.readLine();
			record=reader.readLine();
			h1x[i]=BML.giveDouble(record,1);
			h1y[i]=BML.giveDouble(record,2);
			h1z[i]=BML.giveDouble(record,3);
			record=reader.readLine();
			record=reader.readLine();
			h2x[i]=BML.giveDouble(record,1);
			h2y[i]=BML.giveDouble(record,2);
			h2z[i]=BML.giveDouble(record,3);

		    }
		reader.close();
	    }
	catch(FileNotFoundException e)
	    {
		monitor.println("Error - file not found: "+waterFile);
		return -2;
	    }
	catch(Exception e)
	    {
		monitor.println("Error reading file: "+waterFile+" "+e);
		return -3;
	    }
	monitor.println("Water configuration read");

	return naa;
    }
}
