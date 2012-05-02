package ch.psi.tomcat.tipl;

import java.io.*;
import java.net.*;
import java.security.Principal;
import java.util.*;
import java.io.IOException;
import java.io.PrintStream;
import java.net.InetAddress;
import java.lang.Math;
import java.awt.event.*;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;


/** DistLabel is the class used for labeling bubbles based on a distance map and a mask */
public class DistLabel extends TIPLPlugin {
	public String PluginName() { return "DistLabel";}
	
	public int[] distmap;
	public volatile int[] labels;
	public boolean[] mask;
	public boolean[] diffmask;
	
    
	public int maxlabel;
	public int unfilledVox=0;
	
	private final int MAXDIST=4000;
	private final int MAXDISTVAL=32765;
	private final int MAXLABEL=65543;
	private final int OUTERSHELL=2;
	private final double	MINWALLDIST=4;
	private final double	FLATCRIT=3;//0.401;
	private final int MAXOVERLAP=40;
	private final double	ABSMINWALLDIST=2;
    
	private final boolean DISTGROW=true;
	
	private double MINNSCORE=-0.9;
	private int marchMode=0;
	/** Scalar used for scaling the distance map into voxel space */
	public double distScalar=(MAXDIST+0.0)/(MAXDISTVAL+0.0);
	
	public DistLabel(short[] inputmap,D3int idim,D3int ioffset) {
		aimLength=inputmap.length;
		distmap=new int[aimLength];
		mask=new boolean[aimLength];
		for (int i=0; i<aimLength; i++) {
			distmap[i]=(int) inputmap[i];
			mask[i]=distmap[i]>0;
		}
		Init(idim,ioffset);
	}
	/** The constructor 
     * @param inputmap Is the distance map (scaled into a short) against which the bubbles are grown
     * @param inputmask Is the mask containing the volume where the bubbles are
     */
	public DistLabel(short[] inputmap,boolean[] inputmask,D3int idim,D3int ioffset) {
		aimLength=inputmap.length;
		distmap=new int[aimLength];
		mask=inputmask;
		for (int i=0; i<aimLength; i++) {
			distmap[i]=(int) inputmap[i];
		}
		Init(idim,ioffset);
	}
    /** Constructor based on the distance map and mask aim images
     @param imap        Distance map image
     @param imask       Mask to be filld and identified
     */
    public DistLabel(VirtualAim imap, VirtualAim imask) {
        LoadAimData(imap.getIntAim(),imask.getBoolAim(),imap.dim,imap.offset);
    }
    /** Constructor based on just the distance map, assume everything is available mask 
     @param imap        Distance map image
     */
    public DistLabel(VirtualAim imap) {
        LoadAimData(imap.getIntAim(),imap.dim,imap.offset);
    }
	public DistLabel(int[] inputmap,D3int idim,D3int ioffset) {
        LoadAimData(inputmap,idim,ioffset);
	}
	public DistLabel(int[] inputmap,boolean[] inputmask,D3int idim,D3int ioffset) {
		LoadAimData(inputmap,inputmask,idim,ioffset);
	}
    protected void LoadAimData(int[] inputmap,D3int idim,D3int ioffset) {
        aimLength=inputmap.length;
		distmap=inputmap;
		for (int i=0; i<aimLength; i++) mask[i]=distmap[i]>0;
		Init(idim,ioffset);
    }
    protected void LoadAimData(int[] inputmap,boolean[] inputmask,D3int idim,D3int ioffset) {
        aimLength=inputmap.length;
		distmap=inputmap;
		mask=inputmask;
		Init(idim,ioffset);
    }
	private void Init(D3int idim,D3int ioffset) {
		if (distmap.length!=mask.length) {
			System.out.println("SIZES DO NOT MATCH!!!!!!!!");
			return;
		}
		
		labels=new int[aimLength];
		diffmask=new boolean[aimLength];
        
		
		InitDims(idim,ioffset);
		
		InitDiffmask();
		
		
		
	}
    final private boolean useLaplacian=false;
    /** temporary solution until a more robust stack data managament gets integrated, since java cannot handle more than 8e9 array elements AIM style referencing from IPL will not work **/
    private int dmapGet(int x, int y, int z) {
        return distmap[(z*dim.y + y)*dim.x + x];
    }
    /** temporary solution until a more robust stack data managament gets integrated, since java cannot handle more than 8e9 array elements AIM style referencing from IPL will not work **/
    private int[][][] dmapGet(int x, int y, int z, D3int nSize) {
        int off = (z*dim.y + y)*dim.x + x;
        int[][][] outMap=new int[nSize.x*2+1][nSize.z*2+1][nSize.z*2+1];
        int dz=0;
        for(int z2=max(z-nSize.z,lowz); z2<=min(z+nSize.z,uppz-1); z2++,dz++) {
            int dy=0;
            for(int y2=max(y-nSize.y,lowy); y2<=min(y+nSize.y,uppy-1); y2++,dy++) {
                int off2 = (z2*dim.y + y2)*dim.x + max(x-nSize.x,lowx);
                int dx=0;
                for(int x2=max(x-nSize.x,lowx); x2<=min(x+nSize.x,uppx-1); x2++, off2++,dx++) {
                    outMap[dx][dy][dz]=distmap[off2];
                }
            }
        }
        return outMap;
    }
    
    /** Run the distance label initialization routines in parallel*/
    private class dlRunner extends Thread {
        int sslice,fslice,asType;
        volatile DistLabel parent;
        public boolean isFinished=false;
        public boolean hasStarted=false;
        public SeedList threadSeedList=null;
        public dlRunner(DistLabel iparent,int isslice,int ifslice) {
            super("dlRunner:<"+isslice+", "+ifslice+">");
            sslice=isslice;
            fslice=ifslice;
            parent=iparent;
        }
        /** Distribute (using divideThreadWork) and process (using processWork) the work across the various threads */
        public void run() {
            hasStarted=true;
            isFinished=false;
            
            threadSeedList=parent.locateSeeds(sslice,fslice);
            isFinished=true;
            System.out.println("DlRunner Finished:, <"+sslice+", "+fslice+">"+threadSeedList+", k-"+threadSeedList.killedLabels);
        }
    }
    
    /** Object to divide the thread work into supportCores equal parts, default is z-slices */
    public int[] divideSlices(int cThread) {
        int minSlice=lowz+OUTERSHELL;
        int maxSlice=(uppz-OUTERSHELL);
        
        int range=(maxSlice-minSlice)/neededCores;
        
        int startSlice=minSlice;
        int endSlice=startSlice+range;
        
        for (int i=0;i<cThread;i++) {
            startSlice=endSlice; // must overlap since i<endSlice is always used, endslice is never run
            endSlice=startSlice+range;
        }
        if (cThread==(neededCores-1)) endSlice=maxSlice;
        if (cThread>=neededCores) return null;
        return (new int[] {startSlice,endSlice});    
    }
    SeedList startSeedList;
    
    public void InitDiffmask() {
        
        Thread myThread = Thread.currentThread();
        ArrayList<dlRunner> threadList=new ArrayList<dlRunner>();
        startSeedList=new SeedList(MAXOVERLAP);
        
        jStartTime=System.currentTimeMillis();
        // Call the other threads
        for (int i=0; i<neededCores; i++) {             // setup the background threads
            int[] mySlices=divideSlices(i);
            dlRunner bgThread = new dlRunner(this,mySlices[0],mySlices[1]);
            threadList.add(bgThread);
            bgThread.start();
        }
        
        for (dlRunner theThread : threadList) {    // for all other threads:
            try {
                theThread.join();               // wait until thread has finished
            } catch (InterruptedException e) {
                System.out.println("ERROR - Thread : "+theThread+" was interrupted, proceed carefully!");
            }
            startSeedList.killedLabels+=theThread.threadSeedList.killedLabels;
            for (SeedLabel cBubble : theThread.threadSeedList.sl) if (cBubble.isValid()) startSeedList.addbubble(cBubble.cx,cBubble.cy,cBubble.cz,cBubble.rad);
            
        }
        startSeedList.finalize();
        
        String outString="BubbleSeeds: Ran in "+StrRatio(System.currentTimeMillis()-jStartTime,1000)+" seconds on "+neededCores+" cores and found:"+startSeedList.toString();
        System.out.println(outString);
        procLog+=outString+"\n";
    }
    
    
    
	private SeedList locateSeeds(int startSlice, int finalSlice) {
        SeedList tempSeedList=new SeedList(MAXOVERLAP);
		D3int iNeighborSize=new D3int(2);
		int off=0;
		double avgGrad=0;
		double avgSGrad=0;
		double avgLap=0;
		double avgSLap=0;
		int inVox=0;
		unfilledVox=0;
		for(int z=startSlice; z<finalSlice; z++) {
            for(int y=lowy+OUTERSHELL; y<(uppy-OUTERSHELL); y++) {
                off = (z*dim.y + y)*dim.x + lowx+OUTERSHELL;
                for(int x=lowx+OUTERSHELL; x<(uppx-OUTERSHELL); x++, off++) {
					// The code is optimized so the least number of voxels make it past the first check
                    float cVDist=(float) distScalar*distmap[off];
                    
					if ((mask[off]) && ((cVDist)>MINWALLDIST)) {
						unfilledVox++;
						double gradX=0.0,gradY=0.0,gradZ=0.0;
						double lapVal=0.0;
						int lapCount=0;
						int gradCount=0;
                        double[][] hessianMatrix = new double[3][3];
                        double minCriteria=0; // variable for the minimum criteria (either laplacian or 
						if (useLaplacian) {
                            for(int z2=max(z-iNeighborSize.z,lowz); z2<=min(z+iNeighborSize.z,uppz-1); z2++) {
                                for(int y2=max(y-iNeighborSize.y,lowy); y2<=min(y+iNeighborSize.y,uppy-1); y2++) {
                                    int off2 = (z2*dim.y + y2)*dim.x + max(x-iNeighborSize.x,lowx);
                                    for(int x2=max(x-iNeighborSize.x,lowx); x2<=min(x+iNeighborSize.x,uppx-1); x2++, off2++) {
                                        if (off!=off2) {
                                            if ((Math.abs(x2-x)<=1) && (Math.abs(x2-x)<=1) && (Math.abs(x2-x)<=1)) { //Local gradient
                                                gradX+=(x2-x)*(distmap[off2]);
                                                gradY+=(y2-y)*(distmap[off2]);
                                                gradZ+=(z2-z)*(distmap[off2]);
                                                gradCount++;
                                            }
                                            
                                            if (((x2!=x) ? 1 : 0)+((y2!=y) ? 1 : 0)+((z2!=z) ? 1 : 0)<=1) { // slightly larger laplacian
                                                lapVal+=(distScalar*(distmap[off2]));
                                                lapCount++;
                                            }
                                            
                                            // First derivative is 0 and second derivative is less than 0 (local maxima)
                                            
                                        }
                                    }
                                }
                            }
                            lapVal+=-lapCount*(distScalar*(distmap[off]));
                            lapVal/=lapCount;
                            minCriteria=lapVal;
                            gradX/=gradCount;
                            gradY/=gradCount;
                            gradZ/=gradCount;
                        } else {
                            double curVal = 2 * distmap[off];
                            
                            int[][][] localValues=dmapGet(x,y,z,new D3int(1));
                            
                            // xx
                            hessianMatrix[0][0] = localValues[2][1][1]-curVal+localValues[0][1][1]; //dmapGet(x+1,y,z)-curVal+dmapGet(x-1,y,z); 
                            
                            // yy
                            hessianMatrix[1][1] = localValues[1][2][1]-curVal+localValues[1][0][1]; //dmapGet(x,y+1,z)-curVal+dmapGet(x,y-1,z); //img.get(x, y + 1, z) - temp + img.get(x, y - 1, z);
                            
                            // zz
                            hessianMatrix[2][2] = localValues[1][1][2]-curVal+localValues[1][1][0]; // dmapGet(x, y, z + 1) - curVal +dmapGet(x, y, z - 1);
                            
                            // xy
                            hessianMatrix[0][1] = hessianMatrix[1][0] = 
                            (
                             (localValues[2][2][1] - localValues[0][2][1]) / 2
                             -
                             (localValues[2][0][1]  - localValues[0][0][1]) / 2
                             ) / 2;
                            
                            
                            // xz
                            hessianMatrix[0][2] = hessianMatrix[2][0] =
                            (
                             (localValues[2][1][2] - localValues[0][1][2]) / 2
                             -
                             (localValues[2][1][0] - localValues[0][1][0]) / 2
                             ) / 2;
                            
                            // yz
                            hessianMatrix[1][2] = hessianMatrix[2][1] =
                            (
                             (localValues[1][2][2]- localValues[1][0][2]) / 2
                             -
                             (localValues[1][2][0] - localValues[1][0][0]) / 2
                             ) / 2;
                            
                            Matrix M = new Matrix(hessianMatrix);
                            EigenvalueDecomposition E = new EigenvalueDecomposition(M);
                            
                            double[] result = E.getRealEigenvalues();
                            for(int iz=-1;iz<=1;iz++) {
                                for(int iy=-1;iy<=1;iy++) {
                                    for(int ix=-1;ix<=1;ix++) {
                                        gradX+=ix*localValues[ix+1][iy+1][iz+1];
                                        gradY+=iy*localValues[ix+1][iy+1][iz+1];
                                        gradZ+=iz*localValues[ix+1][iy+1][iz+1];
                                        gradCount++;
                                    }
                                }
                            }
                            if ((result[0]<=0) && (result[1]<=0) && (result[2]<=0)) minCriteria=-5; // Negative Definite
                            else lapVal=1;
                        }
                        
						gradX*=distScalar;
                        gradY*=distScalar;
                        gradZ*=distScalar;
						
						double cGrad=Math.sqrt(gradX*gradX+gradY*gradY+gradZ*gradZ);
						
                        avgGrad+=cGrad;
						avgSGrad+=cGrad*cGrad;
						avgLap+=lapVal;
						avgSLap+=lapVal*lapVal;
						inVox++;
						if ((cGrad<=FLATCRIT) && (minCriteria<-0.1)) {
							//System.out.println("GradVal:"+cGrad+"->("+gradX+", "+gradY+", "+gradZ+"), "+gradCount+", Lap:"+lapVal+", LC"+lapCount);
							diffmask[off]=true;
                            tempSeedList.addbubble(x,y,z,cVDist);
						}
					}	// End mask and dist-check
                } // End x
            } // End y
		} // End z
		avgGrad/=inVox;
		avgLap/=inVox;
		System.out.println("Average GradVal:"+avgGrad+" - STD:"+Math.sqrt(avgSGrad/inVox-avgGrad*avgGrad)+", Lap:"+avgLap+" - STD:"+Math.sqrt(avgSLap/inVox-avgLap*avgLap));
		//procLog+="CMD:LocalMaxima: GradVal:"+avgGrad+" - STD:"+Math.sqrt(avgSGrad/inVox-avgGrad*avgGrad)+", Lap:"+avgLap+" - STD:"+Math.sqrt(avgSLap/inVox-avgLap*avgLap)+"\n";
        return tempSeedList;
        
	}
    private class bubbleFiller extends Thread {
        int cLabel;
        SeedLabel cBubble;
        volatile DistLabel parent;
        public boolean isFinished=false;
        public boolean hasStarted=false;
        public bubbleFiller(DistLabel iparent,int nLabel, SeedLabel nBubble, int core) {
            super("BubbleFiller["+core+"]:<"+nLabel+", "+nBubble+">");
            cBubble=nBubble;
            cLabel=nLabel;
            parent=iparent;
        }
        /** Distribute (using divideThreadWork) and process (using processWork) the work across the various threads */
        public void run() {
            hasStarted=true;
            isFinished=false;
            long ist=System.currentTimeMillis();
            try {
                parent.fillBubble(cLabel,cBubble.cx,cBubble.cy,cBubble.cz,cBubble.rad);
            } catch (Exception e) {
                System.out.println("ERROR - Thread : "+this+" has crashed, proceed carefully!");
                e.printStackTrace();
            }
            isFinished=true;
            System.out.println("BubbleFiller Finished:"+StrRatio(System.currentTimeMillis()-ist,1000)+":, <"+this+">");
        }
    }
    int remVoxels=aimLength;
    int totalVoxels=aimLength;
    protected boolean isCoreFree(bubbleFiller[] coreArray,int coreIndex) {
        if (coreArray[coreIndex]!=null) {
            return (coreArray[coreIndex].isFinished);
        } else
            return true;
    }
    public void run() {
		
		int curLabel=1;
		long cVox=0;
		long sVox=0;
		int off=0;
		short val;
		
		int maxVal=-1;
		double emptyVoxels=0;
		double outsideMask=0;
		double fullVoxels=0;
        
		
		
		// Identify local maxima regions
        
		
		// Run loop to make bubbles
		int clabel=1;
		
        Thread myThread = Thread.currentThread();
        jStartTime=System.currentTimeMillis();
        bubbleFiller[] bfArray=new bubbleFiller[neededCores];
        
        // Call the other threads
        int curCore=0;
        for (SeedLabel cBubble : startSeedList.sl) {
            int toff = (cBubble.cz*dim.y + cBubble.cy)*dim.x + cBubble.cx;
            if (diffmask[toff]) {
                while (!isCoreFree(bfArray,curCore)) {
                    curCore++;
                    if (curCore>=neededCores) {
                        curCore=0;
                        //myThread.sleep(10);
                        myThread.yield();
                    }
                }
                bfArray[curCore]=new bubbleFiller(this,clabel,cBubble,curCore);
                bfArray[curCore].start();
                clabel++;
            }
        }
        // Wind down
        for(int i=0;i<neededCores;i++) {
            if (bfArray[i]!=null) {
                try {
                    bfArray[i].join();               // wait until thread has finished
                } catch (InterruptedException e) {
                    System.out.println("ERROR - Thread : "+bfArray[i]+" was interrupted, proceed carefully!");
                }
            }
        }
        
        
        
        String outString="BubbleFiller: Ran in "+StrRatio(System.currentTimeMillis()-jStartTime,1000)+" seconds on "+neededCores+" cores";
        System.out.println(outString);
        procLog+=outString+"\n";
        
		procLog+="CMD:DistLabel: Max Label:"+clabel+"\n";
		runCount++;
	}
    protected void fillBubble(int clabel,int xmax,int ymax, int zmax, double cMaxVal) {
        
        int off;
        int bubbleSize=0;
        
        // Fill in the bubble based on the distance map
        
        // Check for Overlap with Other Bubbles!
        int overlapVox=0;
        int selfVox=0;
        boolean changedPos=true;
        
        double nDist=cMaxVal;
        int nxi,nyi,nzi;
        int	tlowx,tlowy,tlowz,tuppx,tuppy,tuppz;
        // Initial Values
        tlowx=max(xmax-cMaxVal,lowx);
        tlowy=max(ymax-cMaxVal,lowy);
        tlowz=max(zmax-cMaxVal,lowz);
        tuppx=min(xmax+cMaxVal,uppx);
        tuppy=min(ymax+cMaxVal,uppy);
        tuppz=min(zmax+cMaxVal,uppz);
        // Movable Values
        double startX,startY,startZ,startR;
        startX=xmax;
        startY=ymax;
        startZ=zmax;
        startR=cMaxVal;
        // Initialize Variables
        nxi=xmax;
        nyi=ymax;
        nzi=zmax;
        
        // Find maximum and verify overlap
        while (changedPos) {
            
            overlapVox=0;
            selfVox=0;
            changedPos=false;
            
            tlowx=max(xmax-cMaxVal,lowx);
            tlowy=max(ymax-cMaxVal,lowy);
            tlowz=max(zmax-cMaxVal,lowz);
            tuppx=min(xmax+cMaxVal,uppx);
            tuppy=min(ymax+cMaxVal,uppy);
            tuppz=min(zmax+cMaxVal,uppz);
            
            
            for(int z=tlowz; z<tuppz; z++) {
                double zd=(zmax-z)*(zmax-z);
                for(int y=tlowy; y<tuppy; y++) {
                    double yd=(ymax-y)*(ymax-y);
                    off = (z*dim.y + y)*dim.x + lowx;
                    off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
                    for(int x=tlowx; x<tuppx; x++, off++) {
                        if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
                            if (labels[off]>0) overlapVox++;
                            if (mask[off]) selfVox++;
                            if ((distScalar*distmap[off])>nDist) {
                                nDist=(distScalar*distmap[off]);
                                nxi=x;
                                nyi=y;
                                nzi=z;
                                changedPos=true;
                            }
                        }
                    }
                }
            }
            if (changedPos) {
                xmax=nxi;
                ymax=nyi;
                zmax=nzi;
                cMaxVal=nDist;
            }
        }
        if ((overlapVox*100)/(selfVox+overlapVox+0.1)>MAXOVERLAP) {
            int silenced=0;
            
            for(int cleanVal=0;cleanVal<2;cleanVal++) {
                if (cleanVal==1) {
                    xmax=(int) startX;
                    ymax=(int) startY;
                    zmax=(int) startZ;
                    cMaxVal=startR;
                }
                tlowx=max(xmax-cMaxVal,lowx);
                tlowy=max(ymax-cMaxVal,lowy);
                tlowz=max(zmax-cMaxVal,lowz);
                tuppx=min(xmax+cMaxVal,uppx);
                tuppy=min(ymax+cMaxVal,uppy);
                tuppz=min(zmax+cMaxVal,uppz);
                for(int z=tlowz; z<tuppz; z++) {
                    double zd=(zmax-z)*(zmax-z);
                    for(int y=tlowy; y<tuppy; y++) {
                        double yd=(ymax-y)*(ymax-y);
                        off = (z*dim.y + y)*dim.x + lowx;
                        off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
                        for(int x=tlowx; x<tuppx; x++, off++) {
                            if (mask[off]) {
                                if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
                                    diffmask[off]=false;
                                    silenced++;
                                }
                            }
                        }
                    }
                }
            }
            System.out.println("excessive overlap: "+(overlapVox*100.0)/(selfVox+overlapVox)+"%, removing DiffMask Eligibility, "+silenced+" silenced vx");
        } else {
            
            tlowx=max(xmax-cMaxVal,lowx);
            tlowy=max(ymax-cMaxVal,lowy);
            tlowz=max(zmax-cMaxVal,lowz);
            tuppx=min(xmax+cMaxVal,uppx);
            tuppy=min(ymax+cMaxVal,uppy);
            tuppz=min(zmax+cMaxVal,uppz);
            
            for(int z=tlowz; z<tuppz; z++) {
                double zd=(zmax-z)*(zmax-z);
                for(int y=tlowy; y<tuppy; y++) {
                    double yd=(ymax-y)*(ymax-y);
                    off = (z*dim.y + y)*dim.x + lowx;
                    off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
                    for(int x=tlowx; x<tuppx; x++, off++) {
                        if ((distScalar*distmap[off])>ABSMINWALLDIST) {
                            if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
                                if (mask[off]) {
                                    mask[off]=false;
                                    diffmask[off]=false;
                                    labels[off]=clabel;
                                    unfilledVox--;
                                    bubbleSize++;
                                    
                                    
                                } else {
                                    // Already occupied, empty voxel since it is shared between two (or more) bubbles
                                    labels[off]=0;
                                    
                                }
                            }
                        }
                    }
                    
                }
            }
            double bubRad=Math.pow(bubbleSize/(1.33*3.14159),0.333);
            
            System.out.println("CurBubble("+clabel+") Filled.Rad: "+bubRad+" Voxels: "+bubbleSize+", Overlap: "+(overlapVox*100)/(selfVox+overlapVox+1)+" %");
            if (DISTGROW) {
                int changes=1;
                int iter=0;
                while (changes>0) {
                    
                    changes=0;
                    // Fill in bubble with distgrow
                    for(int z=tlowz; z<tuppz; z++) {
                        //System.out.println+"Slice :: "+z+endl;
                        for(int y=tlowy; y<tuppy; y++) {
                            off = (z*dim.y + y)*dim.x + lowx;
                            off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
                            for(int x=tlowx; x<tuppx; x++, off++) {
                                // [off] is the voxel to be filled
                                // [off2] is the current voxel (filled, the one doing the expanding)
                                
                                boolean distMatch=false;
                                if (mask[off]) {// Eligble empty voxel
                                    double gradX=0.0,gradY=0.0,gradZ=0.0;
                                    
                                    for(int z2=max(z-neighborSize.z,lowz); z2<=min(z+neighborSize.z,uppz-1); z2++) {
                                        for(int y2=max(y-neighborSize.y,lowy); y2<=min(y+neighborSize.y,uppy-1); y2++) {
                                            int off2 = (z2*dim.y + y2)*dim.x + max(x-neighborSize.x,lowx);
                                            for(int x2=max(x-neighborSize.x,lowx); x2<=min(x+neighborSize.x,uppx-1); x2++, off2++) {
                                                if ((off!=off2) && (!distMatch)) {
                                                    totalVoxels++;
                                                    gradX+=(x2-x)*distScalar*(distmap[off2]);
                                                    gradY+=(y2-y)*distScalar*(distmap[off2]);
                                                    gradZ+=(z2-z)*distScalar*(distmap[off2]);
                                                    if ((labels[off2]==clabel))  { // valid seed voxel
                                                        
                                                        // Now check distance map
                                                        double cDist=Math.sqrt(Math.pow(x2-x,2)+Math.pow(y2-y,2)+Math.pow(z2-z,2));
                                                        double slope=((double) distScalar*distmap[off]-(double) distScalar*distmap[off2])/cDist;
                                                        //System.out.println+"Current Slope:"+slope+endl;
                                                        if ((distmap[off]==cMaxVal)) {
                                                            distMatch=true;
                                                            //System.out.println+"Flooding..."+endl;
                                                        }
                                                        else {
                                                            switch (marchMode) {
                                                                case 0: // Expands against distance gradient
                                                                    distMatch=slope<MINNSCORE;
                                                                    break;
                                                                case 1: // Expands against distance gradient
                                                                    distMatch=slope<=MINNSCORE;
                                                                    break;
                                                                case 2: // Expands with distance gradient
                                                                    distMatch=slope>MINNSCORE;
                                                                    break;
                                                                case 3: // Expands with distance gradient
                                                                    distMatch=slope>=MINNSCORE;
                                                                    break;
                                                            }
                                                        }
                                                        if (distMatch) break; // break only breaks out of one loop so added (!distMatch line above)
                                                    }
                                                    
                                                }
                                                
                                            }
                                        }
                                    }
                                    
                                    if (distMatch) {
                                        //System.out.println+"Slope Match: "+slope*dst[ik]+endl;
                                        labels[off]=clabel;
                                        unfilledVox--;
                                        mask[off]=false; 
                                        tlowx=min(tlowx,max(x-1,lowx));
                                        tlowy=min(tlowy,max(y-1,lowy));
                                        tlowz=min(tlowz,max(z-1,lowz));
                                        tuppx=max(tuppx,min(x+1,uppx));
                                        tuppy=max(tuppy,min(y+1,uppy));
                                        tuppz=max(tuppz,min(z+1,uppz));
                                        changes++;
                                        bubbleSize++;
                                    }
                                    
                                }
                                
                            }
                        }
                    }
                    iter++;
                    
                }
                System.out.println("Iters:"+iter+" : "+bubbleSize+" unfilled voxels (%) "+String.format("%.2f",(unfilledVox*100.0)/(aimLength)));
            }
            remVoxels-=bubbleSize;
        }
    }
	public void oldrun() {
		
		int curLabel=1;
		long cVox=0;
		long sVox=0;
		int off=0;
		short val;
		
		int maxVal=-1;
		double emptyVoxels=0;
		double outsideMask=0;
		double fullVoxels=0;
        
		int remVoxels=aimLength;
		int totalVoxels=aimLength;
		int	tlowx,tlowy,tlowz,tuppx,tuppy,tuppz;
		// Identify local maxima regions
        
		
		// Run loop to make bubbles
		int clabel=1;
		while (clabel<MAXLABEL) {
			System.out.println("Scanning Image for largest remaining local maximum...");
			double cMaxVal=0;
			int xmax=-1,ymax=-1,zmax=-1;
			for(int z=lowz+OUTERSHELL; z<(uppz-OUTERSHELL); z++) {
				for(int y=lowy+OUTERSHELL; y<(uppy-OUTERSHELL); y++) {
					off = (z*dim.y + y)*dim.x + lowx+OUTERSHELL;
					for(int x=lowx+OUTERSHELL; x<(uppx-OUTERSHELL); x++, off++) {
						if ((distScalar*distmap[off])>MINWALLDIST) {
							if (mask[off] && diffmask[off]) {
								if ((distScalar*distmap[off])>cMaxVal) {
									cMaxVal=(distScalar*distmap[off]);
									xmax=x;
									ymax=y;
									zmax=z;
								}
							}
						}
					}
				}
			}
			if (cMaxVal<1) {
				System.out.println("No more sufficient seed bubbles found!");
				break;
			}
            
			
			int bubbleSize=0;
			
			// Fill in the bubble based on the distance map
            
			// Check for Overlap with Other Bubbles!
			int overlapVox=0;
			int selfVox=0;
			boolean changedPos=true;
			
			double nDist=cMaxVal;
			int nxi,nyi,nzi;
			// Initial Values
			tlowx=max(xmax-cMaxVal,lowx);
			tlowy=max(ymax-cMaxVal,lowy);
			tlowz=max(zmax-cMaxVal,lowz);
			tuppx=min(xmax+cMaxVal,uppx);
			tuppy=min(ymax+cMaxVal,uppy);
			tuppz=min(zmax+cMaxVal,uppz);
			// Movable Values
			double startX,startY,startZ,startR;
			startX=xmax;
			startY=ymax;
			startZ=zmax;
			startR=cMaxVal;
			// Initialize Variables
			nxi=xmax;
			nyi=ymax;
			nzi=zmax;
			
			// Find maximum and verify overlap
			while (changedPos) {
				
				overlapVox=0;
				selfVox=0;
				changedPos=false;
				
				tlowx=max(xmax-cMaxVal,lowx);
				tlowy=max(ymax-cMaxVal,lowy);
				tlowz=max(zmax-cMaxVal,lowz);
				tuppx=min(xmax+cMaxVal,uppx);
				tuppy=min(ymax+cMaxVal,uppy);
				tuppz=min(zmax+cMaxVal,uppz);
				
				
				for(int z=tlowz; z<tuppz; z++) {
					double zd=(zmax-z)*(zmax-z);
					for(int y=tlowy; y<tuppy; y++) {
						double yd=(ymax-y)*(ymax-y);
						off = (z*dim.y + y)*dim.x + lowx;
						off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
						for(int x=tlowx; x<tuppx; x++, off++) {
							if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
								if (labels[off]>0) overlapVox++;
								if (mask[off]) selfVox++;
								if ((distScalar*distmap[off])>nDist) {
									nDist=(distScalar*distmap[off]);
									nxi=x;
									nyi=y;
									nzi=z;
									changedPos=true;
								}
							}
						}
					}
				}
				if (changedPos) {
					xmax=nxi;
					ymax=nyi;
					zmax=nzi;
					cMaxVal=nDist;
				}
			}
			if ((overlapVox*100)/(selfVox+overlapVox+0.1)>MAXOVERLAP) {
				int silenced=0;
				
				for(int cleanVal=0;cleanVal<2;cleanVal++) {
					if (cleanVal==1) {
						xmax=(int) startX;
						ymax=(int) startY;
						zmax=(int) startZ;
						cMaxVal=startR;
					}
					tlowx=max(xmax-cMaxVal,lowx);
					tlowy=max(ymax-cMaxVal,lowy);
					tlowz=max(zmax-cMaxVal,lowz);
					tuppx=min(xmax+cMaxVal,uppx);
					tuppy=min(ymax+cMaxVal,uppy);
					tuppz=min(zmax+cMaxVal,uppz);
					for(int z=tlowz; z<tuppz; z++) {
						double zd=(zmax-z)*(zmax-z);
						for(int y=tlowy; y<tuppy; y++) {
							double yd=(ymax-y)*(ymax-y);
							off = (z*dim.y + y)*dim.x + lowx;
							off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
							for(int x=tlowx; x<tuppx; x++, off++) {
								if (mask[off]) {
									if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
										diffmask[off]=false;
										silenced++;
									}
								}
							}
						}
					}
				}
				System.out.println("excessive overlap: "+(overlapVox*100.0)/(selfVox+overlapVox)+"%, removing DiffMask Eligibility, "+silenced+" silenced vx");
			} else {
				
				tlowx=max(xmax-cMaxVal,lowx);
				tlowy=max(ymax-cMaxVal,lowy);
				tlowz=max(zmax-cMaxVal,lowz);
				tuppx=min(xmax+cMaxVal,uppx);
				tuppy=min(ymax+cMaxVal,uppy);
				tuppz=min(zmax+cMaxVal,uppz);
				
				for(int z=tlowz; z<tuppz; z++) {
					double zd=(zmax-z)*(zmax-z);
					for(int y=tlowy; y<tuppy; y++) {
						double yd=(ymax-y)*(ymax-y);
						off = (z*dim.y + y)*dim.x + lowx;
						off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
						for(int x=tlowx; x<tuppx; x++, off++) {
							if ((distScalar*distmap[off])>ABSMINWALLDIST) {
								if (Math.sqrt((xmax-x)*(xmax-x)+yd+zd)<cMaxVal) {
									if (mask[off]) {
										mask[off]=false;
										diffmask[off]=false;
										labels[off]=clabel;
										unfilledVox--;
										bubbleSize++;
										
										
									} else {
										// Already occupied, empty voxel since it is shared between two (or more) bubbles
										labels[off]=0;
										
									}
								}
							}
						}
						
					}
				}
				double bubRad=Math.pow(bubbleSize/(1.33*3.14159),0.333);
				
				System.out.println("CurBubble("+clabel+") Filled.Rad: "+bubRad+" Voxels: "+bubbleSize+", Remaining: "+((remVoxels*100)/totalVoxels+1)+" %"+", Overlap: "+(overlapVox*100)/(selfVox+overlapVox+1)+" %");
				if (DISTGROW) {
					int changes=1;
					int iter=0;
					while (changes>0) {
						
						changes=0;
						// Fill in bubble with distgrow
						for(int z=tlowz; z<tuppz; z++) {
							//System.out.println+"Slice :: "+z+endl;
							for(int y=tlowy; y<tuppy; y++) {
								off = (z*dim.y + y)*dim.x + lowx;
								off+=(tlowx-lowx);// since loop in next line starts at tlowx instead of lowx
								for(int x=tlowx; x<tuppx; x++, off++) {
									// [off] is the voxel to be filled
									// [off2] is the current voxel (filled, the one doing the expanding)
									
									boolean distMatch=false;
									if (mask[off]) {// Eligble empty voxel
										double gradX=0.0,gradY=0.0,gradZ=0.0;
										
										for(int z2=max(z-neighborSize.z,lowz); z2<=min(z+neighborSize.z,uppz-1); z2++) {
											for(int y2=max(y-neighborSize.y,lowy); y2<=min(y+neighborSize.y,uppy-1); y2++) {
												int off2 = (z2*dim.y + y2)*dim.x + max(x-neighborSize.x,lowx);
												for(int x2=max(x-neighborSize.x,lowx); x2<=min(x+neighborSize.x,uppx-1); x2++, off2++) {
													if ((off!=off2) && (!distMatch)) {
														totalVoxels++;
														gradX+=(x2-x)*distScalar*(distmap[off2]);
														gradY+=(y2-y)*distScalar*(distmap[off2]);
														gradZ+=(z2-z)*distScalar*(distmap[off2]);
														if ((labels[off2]==clabel))  { // valid seed voxel
															
															// Now check distance map
															double cDist=Math.sqrt(Math.pow(x2-x,2)+Math.pow(y2-y,2)+Math.pow(z2-z,2));
															double slope=((double) distScalar*distmap[off]-(double) distScalar*distmap[off2])/cDist;
															//System.out.println+"Current Slope:"+slope+endl;
															if ((distmap[off]==cMaxVal)) {
																distMatch=true;
																//System.out.println+"Flooding..."+endl;
															}
															else {
																switch (marchMode) {
																	case 0: // Expands against distance gradient
																		distMatch=slope<MINNSCORE;
																		break;
																	case 1: // Expands against distance gradient
																		distMatch=slope<=MINNSCORE;
																		break;
																	case 2: // Expands with distance gradient
																		distMatch=slope>MINNSCORE;
																		break;
																	case 3: // Expands with distance gradient
																		distMatch=slope>=MINNSCORE;
																		break;
																}
															}
															if (distMatch) break; // break only breaks out of one loop so added (!distMatch line above)
														}
														
													}
													
												}
											}
										}
										
										if (distMatch) {
											//System.out.println+"Slope Match: "+slope*dst[ik]+endl;
											labels[off]=clabel;
											unfilledVox--;
											mask[off]=false; 
											tlowx=min(tlowx,max(x-1,lowx));
											tlowy=min(tlowy,max(y-1,lowy));
											tlowz=min(tlowz,max(z-1,lowz));
											tuppx=max(tuppx,min(x+1,uppx));
											tuppy=max(tuppy,min(y+1,uppy));
											tuppz=max(tuppz,min(z+1,uppz));
											changes++;
											bubbleSize++;
										}
										
									}
									
								}
							}
						}
						iter++;
						
					}
					System.out.println("Iters:"+iter+" : "+bubbleSize+" unfilled voxels (%) "+String.format("%.2f",(unfilledVox*100.0)/(aimLength)));
				}
				remVoxels-=bubbleSize;
				clabel++;
			}
			
			//System.out.println("Distance :"+String.format("%.2f",fDist)+"("+curIter+"), valid(%) : "+String.format("%.2f",(100.0*validVoxels)/(emptyVoxels+fullVoxels+1.0))+", changes(pm) : "+String.format("%.2f",(1000.0*changescnt)/(emptyVoxels+fullVoxels+1.0))+", porosity(%) : "+String.format("%.2f",(100.*emptyVoxels)/(emptyVoxels+fullVoxels+1.0)));
			
		}
		procLog+="CMD:DistLabel: Max Label:"+clabel+"\n";
		runCount++;
	}
    /** Write the labeled bubbles to an image based on the template */
	public VirtualAim ExportAim(VirtualAim templateAim) {
		if (isInitialized) {
			if (runCount>0) {
				VirtualAim outAimData=templateAim.inheritedAim(labels,dim,offset);
                outAimData.appendProcLog(procLog);
                return outAimData;
			} else {
				System.err.println("The plug-in : "+PluginName()+", has not yet been run, exported does not exactly make sense, original data will be sent.");
				return templateAim.inheritedAim(labels,dim,offset);
			}
		} else {
			System.err.println("The plug-in : "+PluginName()+", has not yet been initialized, exported does not make any sense");
			return templateAim.inheritedAim(templateAim.getBoolAim(),dim,offset);
			
		}
	}
    /** export the bubble seeds if anyone actually wants them */
	public VirtualAim ExportBubbleseedsAim(VirtualAim templateAim) {
		if (isInitialized) {
            VirtualAim outAimData=templateAim.inheritedAim(diffmask,dim,offset);
            outAimData.appendProcLog(procLog);
            return outAimData;
			
		} else {
			System.err.println("The plug-in : "+PluginName()+", has not yet been initialized, exported does not make any sense");
			return templateAim.inheritedAim(templateAim.getBoolAim(),dim,offset);
			
		}
	}
    
    class SeedLabel implements Comparable {
        //public int label;
        
        public int cx;
        public int cy;
        public int cz;
        
        public float rad;
        public boolean valid;
        
        
        // Scale distance by STD
        
        public SeedLabel(int x, int y, int z, float irad) {
            cx=x;
            cy=y;
            cz=z;
            rad=irad;
            valid=true;
        }
        public double dist(int x,int y, int z) {
            return dist((float) x,(float) y,(float) z);
        }
        
        public double dist(float x,float y, float z) {
            return Math.sqrt(Math.pow(x-cx,2)+Math.pow(y-cy,2)+Math.pow(z-cz,2));
            
        }
        /** percentage the bubble nseed overlaps with this bubble (roughly) **/
        public double overlap(SeedLabel nSeed) {
            double cdist=nSeed.dist(cx,cy,cz);
            if (cdist>(rad+nSeed.rad)) return 0;
            else return 100*(rad+nSeed.rad-cdist)/nSeed.rad;
        }
        public int compareTo(Object otherBubble) //  must be defined if we are implementing //Comparable interface
        {
            if(!(otherBubble instanceof SeedLabel))
            {
                throw new ClassCastException("Not valid SeedLabel object");
            }
            SeedLabel tempName = (SeedLabel)otherBubble;
            // eliminate the duplicates when you sort 
            if(this.rad >tempName.rad)
            {
                return -1;
            }else if (this.rad < tempName.rad){
                return 1;
            }else{
                return 0;
            }
        }
        public void superSeed(int x, int y, int z, float nrad) {
            String curString=toString();
            cx=x;cy=y;cz=z;rad=nrad;
            String newString=toString();
            System.out.println("Too much overlap, bubble being replaced : "+curString+" -> "+newString);
        }
        public boolean isValid() { return valid;}
        public void invalidate() { valid=false;}
        
        public String toString() {
            return "COV:("+cx+","+cy+","+cz+"), Rad:("+rad+")";
        }
    }
    class SeedList {
        public int label;
        ArrayList<SeedLabel> sl;
        public int killedLabels=0;
        public float overlapLimit;
        public boolean scaleddist=false; 
        
        public SeedList(float ioverlapLimit) {
            sl=new ArrayList<SeedLabel>();
            
            overlapLimit=ioverlapLimit;
        }
        public void addbubble(int x, int y, int z, float rad) {
            SeedLabel newSeed=new SeedLabel(x,y,z,rad);
            boolean keepSeed=true;
            boolean hasReplaced=false;
            for(SeedLabel curSeed : sl) {
                if (curSeed.isValid()) {
                    if (newSeed.overlap(curSeed)>=overlapLimit) {
                        if (hasReplaced) {
                            // if it has already replaced another bubble then an excessively overlapping bubble needs to be marked for deletion
                            curSeed.invalidate();
                            killedLabels++;
                        } else {
                            if (newSeed.rad>curSeed.rad) {
                                curSeed.superSeed(x,y,z,rad);
                                killedLabels++;
                                hasReplaced=true;
                            } else {
                                //System.out.println("New Seed Executed:"+newSeed);
                                keepSeed=false;
                                killedLabels++;
                                
                                break;
                            }
                        }
                    }
                }
            }
            if (keepSeed) append(newSeed);
            
        }
        protected void append(SeedLabel newSeed) {
            sl.add(newSeed);
        }
        public void finalize() {
            ArrayList<SeedLabel> filteredSL = new ArrayList<SeedLabel>();
            for (SeedLabel cSeed : sl) {
                if (cSeed.isValid()) {
                    filteredSL.add(cSeed);
                }            
            }
            Collections.sort(filteredSL);
            sl=filteredSL;
        }
        
        public Iterator kiter() {
            return sl.iterator();
        }
        public int length() {
            return sl.size();
        }
        public String toString() {
            double maxRad=-1;
            if (sl.size()>0) maxRad=((SeedLabel) sl.get(0)).rad;
            return " "+length()+" bubbles ("+killedLabels+"), largest with radius:"+maxRad+" ";
        }
        
    }
    
}	





