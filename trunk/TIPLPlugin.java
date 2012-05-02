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


// Used as a replacement for the moment function as it allows much more control over data
// and communication with webservices (potentially?)
/** Abstract Class for performing TIPLPlugin */
abstract public class TIPLPlugin implements Runnable {
	/** Function to return the name of the current plug-in, useful for ProcLog and other measures */
	abstract public String PluginName();
    /** Run function to be executed by the plug-in. This does the majority of the computational work and is thus setup so that it can be run in a seperate thread from the main functions */
	abstract public void run();
    /** All plug-ins have an interface for exporting the main result to an Aim class based on a template aim, many have other methods for exporting the secondary results (distance maps, histograms, shape analyses, but these need to be examined individually
     @param templateAim     TemplateAim is an aim file which will be used in combination with the array of data saved in the plugin to generate a full aim output class (element size, procedural log, etc..)*/
	abstract public VirtualAim ExportAim(VirtualAim templateAim);
    
	/** A generic interface for thresholding operations, typically only one function is needed but the interface provides all anyways */
	public interface TIPLFilter {
        double meanV=0.0;
        double callT=0;
        /** whether or not the group should be accepted based on a boolean */
		boolean accept(boolean voxCount);
        /** whether or not the group should be accepted based on voxel-count */
		boolean accept(int voxCount);
		/** whether or not a voxel/group should be accepted based on a float value */
		boolean accept(float value);
		/** wheter or not a group should be accepted based on a label number (for sorted lists) or voxel count */
		boolean accept(int labelNumber, int voxCount);
		/** A three metric acceptance criteria */
		boolean accept(int labelNumber, int voxCount,int thirdMetric);
	}
	/** The kernel used for filtering operations */
	public interface filterKernel {
		public String filterName();
		/** Given voxel @ x1, y1, z1 the function adds weighting coefficient for voxel ( x2, y2, z2 ), value value. */
		public void addpt(double x1,double x2, double y1,double y2, double z1, double z2,double value);
		/** The final value of voxel x1, y1, z1 */
		public double value();
		/** Reset for use with a new voxel */
		public void reset();
	}
	// Default Kernels
    /** A sobel gradient filter useful for finding edges or performing edge enhancement */
    public static filterKernel gradientFilter() {
        return new filterKernel() {
            double gx=0.0; // average over all voxels
            double gy=0.0;
            double gz=0.0;
            public String filterName() {return "SobelGradient";}
            public void addpt(double oposx,double ox, double oposy,double oy, double oposz, double oz,double dcVox) {
                double cDistXY=Math.max(4-Math.sqrt(Math.pow((oposx-ox+0.0),2.0)+Math.pow((oposy-oy+0.0),2.0)),0.0);
                double cDistYZ=Math.max(4-Math.sqrt(Math.pow((oposz-oz+0.0),2.0)+Math.pow((oposy-oy+0.0),2.0)),0.0);
                double cDistXZ=Math.max(4-Math.sqrt(Math.pow((oposz-oz+0.0),2.0)+Math.pow((oposx-ox+0.0),2.0)),0.0);
                double cDistX=(oposx-ox+0.0);
                double cDistY=(oposy-oy+0.0);
                double cDistZ=(oposz-oz+0.0);
                double cgx=(double) Math.round(cDistX);
                double cgy=(double) Math.round(cDistY);
                double cgz=(double) Math.round(cDistZ);
                if (Math.abs(cgx)>1.5) cgx=0;
                if (Math.abs(cgy)>1.5) cgy=0;
                if (Math.abs(cgz)>1.5) cgz=0;
                gx+=cgx*cDistYZ*dcVox;
                gy+=cgy*cDistXZ*dcVox;
                gz+=cgy*cDistXY*dcVox;
            }
            public double value() {
                return Math.sqrt(gx*gx+gy*gy+gz*gz);
            }
            public void reset() {
                gx=0;
                gy=0;
                gz=0;
            }
        };
    }
    /** A laplacian filter */
    public static filterKernel laplaceFilter() {
        return new filterKernel() {
            double cVoxCent=0.0; // average over all voxels
            double cVoxEdge=0.0;
            int cVoxEdgeCnt=0;
            int cVoxCentCnt=0;
            public String filterName () {return "Laplacian";}
            
            public void addpt(double oposx,double ox, double oposy,double oy, double oposz, double oz,double dcVox) {
                double cDist=Math.sqrt(Math.pow((oposx-ox+0.0),2.0)+Math.pow((oposy-oy+0.0),2.0)+Math.pow((oposz-oz+0.0),2.0));
                if (Math.round(cDist)<1) {
                    cVoxCent+=dcVox;
                    cVoxCentCnt++;
                } else if (Math.round(cDist)<2.1) {
                    cVoxEdge+=dcVox;
                    cVoxEdgeCnt++;
                }
            }
            public double value() {
                if (cVoxCentCnt>0) {
                    if (cVoxEdgeCnt>0) {
                        return Math.sqrt(Math.pow(cVoxEdge/cVoxEdgeCnt-cVoxCent/cVoxCentCnt,2));
                    }
                }
                return 0;
            }
            public void reset() {
                cVoxCent=0.0;
                cVoxEdge=0.0;
                cVoxEdgeCnt=0;
                cVoxCentCnt=0;
            }
        };	
    }
    /** Gaussian spatial filter (smoothing), isotrpic radii radr */
	public static filterKernel gaussFilter(double radr) {
		return gaussFilter(radr,radr,radr);
	}
	
    /** Gaussian spatial filter (smoothing), radii in each direction in voxels */
	public static filterKernel gaussFilter(double radx,double rady, double radz) {
		final double fRadX=radx;
		final double fRadY=rady;
		final double fRadZ=radz;
		filterKernel gaussFilterOut=new filterKernel() {
			double cVox=0.0; // average over all voxels
			double cWeight=0.0; // weight by how much of the voxel is within the field
			public String filterName() {return "Gaussian";}
			public void addpt(double oposx,double ox, double oposy,double oy, double oposz, double oz,double dcVox) {
				double cDist=Math.sqrt(Math.pow((oposx-ox+0.0)/fRadX,2.0)+Math.pow((oposy-oy+0.0)/fRadY,2.0)+Math.pow((oposz-oz+0.0)/fRadZ,2.0));
				cDist=Math.exp(-1*Math.pow(cDist,2));
				cWeight+=cDist;
				cVox+=dcVox*cDist;
			}
			public double value() {
				if (cWeight>0) return (cVox/cWeight);
				return 0;
				
			}
			public void reset() {
				cWeight=0.0;
				cVox=0.0;
			}
		};
		return gaussFilterOut;
	}
    
	/** The kernel used can be custom defined using the morphKernel interface containing just one function */
	public interface morphKernel {
		/** Given voxel (linear position off) @ x1, y1, z1 is (boolean) the voxel ( x2, y2, z2 ) at linear position (off2) inside. 
		 *off and off2 are generally not used can be used for hashtables to speed up operations. 
		 * <li>
		 * <p> Example kernel: Spherical
		 * <p>	sphericalKernel= new morphKernel() {// Semi-axes radx,rady, radz
		 * <p>		public boolean inside(int off,int off2, int x1, int x2, int y1, int y2, int z1, int z2) {
		 * <p>			return (Math.pow((x1-x2)/radx,2)+Math.pow((y1-y2)/rady,2)+Math.pow((z1-z2)/radz,2))<=1;
		 * <p>		}
		 * <p>	};
		 * 
		 */
		public boolean inside(int off,int off2, int x1, int x2, int y1, int y2, int z1, int z2);
	}
	public D3int neighborSize=new D3int(1); // Neighborhood size
	
	/** A ellipsoidal kernel radius radxi, radyi, radzi */
	public morphKernel sphKernel(double radxi,double radyi, double radzi) {
		final double radx=radxi;
		final double rady=radyi;
		final double radz=radzi;
		return new morphKernel() {// Semi-axes radx,rady, radz
			public boolean inside(int off,int off2, int x1, int x2, int y1, int y2, int z1, int z2) {
				return (Math.pow((x1-x2)/radx,2)+Math.pow((y1-y2)/rady,2)+Math.pow((z1-z2)/radz,2))<=1;
			}
		};
	}
	/** A ellipsoidal kernel radius radxi, radyi, radzi */
	public morphKernel sphKernel(double rad) {
		return sphKernel(rad,rad,rad);
	}
	/** The actual kernel to use for TIPLPluginlogical operations, default is null, which means function reverts to its default */
	public morphKernel neighborKernel=null; 
	
	
	/** has the plug-in been properly initialized with data */
	boolean isInitialized=false;
	/** does this plugin support multi-threaded processing */
    public boolean supportsThreading=false;
    /** how many cores are supported (additional threads to launch) */
    public static int supportedCores=Runtime.getRuntime().availableProcessors();
    /** how many cores does the plugin want (-1 = as many as possible)*/
    public int neededCores=-1;
    protected boolean checkMaxCores=true;
    /** which thread ran the script (will normally not work afterwards */
    protected Thread launchThread=null;
    /** The work assigned to each thread */
    protected Hashtable workForThread=null;
    /** Job Start Time*/
    protected long jStartTime;
    
	/** has the plug-in actually been run at least once */
	int runCount=0;
	/** Dimensions of dataset which has been loaded */
	D3int dim;
    /** Dimensions of border around datatset which should not be processed (normally 0) */
	D3int offset;
    /** Position of bottom right corner */
    //D3int pos; 
	
	int aimLength;
	int lowx,lowy,lowz,uppx,uppy,uppz;
	/** Procedure Log for the function, should be added back to the aim-file after function operation is complete */
	public String procLog="";
	public TIPLPlugin() {
		isInitialized=false;
	}
	/** constructor function taking boolean (other castings just convert the array first) linear array and the dimensions */
	public TIPLPlugin(D3int idim,D3int ioffset) {
		InitDims(idim,ioffset);
	}
    /** Object to divide the thread work into supportCores equal parts, default is z-slices */
    public Object divideThreadWork(int cThread) {
        int minSlice=lowz;
        int maxSlice=uppz;
        int myNeededCores=neededCores;
        
        if (2*neededCores>(maxSlice-minSlice)) myNeededCores=(maxSlice-minSlice)/2; // At least 2 slices per core and definitely no overlap
        if (myNeededCores<1) myNeededCores=1;
        int range=(maxSlice-minSlice)/myNeededCores;
        
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
    /** Distribute (using divideThreadWork) and process (using processWork) the work across the various threads, returns true if thread is launch thread */
	public boolean runMulticore() {
        
        Thread myThread = Thread.currentThread();
        if (myThread==launchThread) { // Distribute work load
            workForThread = new Hashtable(neededCores-1);
            jStartTime=System.currentTimeMillis();
            // Call the other threads
            for (int i=1; i<neededCores; i++) {             // setup the background threads
                Object myWork=divideThreadWork(i);
                Thread bgThread = new Thread(this,"BG:"+PluginName()+" : "+myWork);
                workForThread.put(bgThread, myWork);
                bgThread.start();
            }
            
            processWork(divideThreadWork(0)); // Do one share of the work while waiting
            if (workForThread != null) {
                while (workForThread.size()>0) {    // for all other threads:
                    Thread theThread = (Thread)workForThread.keys().nextElement();
                    
                    try {
                        theThread.join();               // wait until thread has finished
                    } catch (InterruptedException e) {
                        System.out.println("ERROR - Thread : "+theThread+" was interrupted, proceed carefully!");
                        
                    }
                    workForThread.remove(theThread);  // and remove it from the list.
                }
            }
            
            String outString="MCJob Ran in "+StrRatio(System.currentTimeMillis()-jStartTime,1000)+" seconds in "+neededCores;
            System.out.println(outString);
            procLog+=outString+"\n";
        } else if (workForThread!=null && workForThread.containsKey(myThread)) { // Process work load
            Object myWork=workForThread.get(myThread);
            if (myWork==null) {
                // gommer go schlofe
                System.out.println("Too many cores for too little work!");
            } else processWork(myWork); 
        } else
            System.out.println("ERROR - TIPLPlugin:"+PluginName()+" internal error:unsolicited background thread:"+myThread);
        return (myThread==launchThread);
    }
    /** Placeholder for the real function which processes the work, basically takes the section of the computational task and processes it 
     @param myWork      myWork (typically int[2] with starting and ending slice) is an object which contains the work which needs to be processed by the given thread*/
    
    protected void processWork(Object myWork) {
        System.out.println("THIS IS AN pseudo-ABSTRACT FUNCTION AND DOES NOTHING, PLEASE EITHER TURN OFF MULTICORE SUPPORT OR REWRITE YOUR PLUGIN!!!!");
        return; // nothing to do
    }
	protected void InitDims(D3int idim,D3int ioffset) {
		dim=new D3int(idim.x,idim.y,idim.z);
		offset=new D3int(ioffset.x,ioffset.y,ioffset.z);

		
		// boundaries for evalution
		lowx = offset.x ;
		lowy = offset.y ;
		lowz = offset.z ;
		
		uppx = dim.x - offset.x ;
		uppy = dim.y - offset.y ;
		uppz = dim.z - offset.z ;
        
        launchThread=Thread.currentThread();
        
        if (neededCores<1) {
            System.out.println(PluginName()+" has "+neededCores+" number of cores, defaulting to :"+supportedCores);
            neededCores=supportedCores;
            
        }
        if (checkMaxCores) neededCores=min(neededCores,supportedCores);
        
		isInitialized=true;
		
		
	}
    /** Provides a series of useful log and screen writing functions */
    public static String StrPctRatio(double a, double b) {return String.format("%.2f",((a*100)/(b)))+"%";}
    /** Provides a series of useful log and screen writing functions */
    public static String StrRatio(double a, double b) {return String.format("%.2f",((a)/(b)));}
    /** Provides a series of useful log and screen writing functions */
    public static String StrMvx(double a) {return String.format("%.2f",(((double) a)/(1.0e6)))+" MVx";}
    
	public static int min(int a,int b) { if (a>b) return b; else return a; }
	public static int max(int a,int b) { if (a>b) return a; else return b; }
	public static int min(double a,int b) { if (a>b) return b; else return ((int) a); }
	public static int max(double a,int b) { if (a>b) return ((int) a); else return b; }
	
	/** stationaryKernel uses a hashtable for stationary ( f(a,b)=f(a+c,b+c) ) kernels to speed up calculations */
	public class stationaryKernel implements TIPLPlugin.morphKernel{ // basically caches stationary kernels
		TIPLPlugin.morphKernel myKernel;
		Hashtable kvalues=new Hashtable();
		String kernelName="Default";
		boolean isFullKernel=false;
		boolean isStationary=true;
		public stationaryKernel(TIPLPlugin.morphKernel inKernel) {
			myKernel=inKernel;
		}
		public stationaryKernel(TIPLPlugin.morphKernel inKernel,String kName) {
			myKernel=inKernel;
			kernelName=kName;
		}
		public stationaryKernel(TIPLPlugin.morphKernel inKernel,String kName, boolean inIsStationary) {
			myKernel=inKernel;
			kernelName=kName;
			isStationary=inIsStationary;
		}
		public stationaryKernel() {
			myKernel=null;
			kernelName="Full";
			isFullKernel=true;
		}
		public boolean inside(int off,int off2,int x,int x2,int y,int y2,int z, int z2) {
			if (isFullKernel) return true;
			int offD=off2-off;
			if (isStationary) if (kvalues.containsKey(offD)) return ((Boolean) kvalues.get(offD)).booleanValue();
			boolean outVal=myKernel.inside(off,off2,x,x2,y,y2,z,z2);
			if (isStationary) kvalues.put(offD,new Boolean(outVal));
			return outVal;
		}
	}
	// Actual kernel definitons at the end because they screw up the formatting in most text editors
	public static morphKernel fullKernel=new morphKernel() {
	public boolean inside(int off,int off2, int x1, int x2, int y1, int y2, int z1, int z2) {
	return true;
}
};
/** A cross like kernel with ones along the x,y,z axes and all other values set to 0 */
public static morphKernel dKernel=new morphKernel() {
public boolean inside(int off, int off2, int x1, int x2, int y1, int y2, int z1, int z2) {
if (((x1==x2 ? 0 : 1)+(y1==y2 ? 0 : 1)+(z1==z2 ? 0 : 1))==1) return true;
else return false;
}
};

}






