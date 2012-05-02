
package ch.psi.tomcat.tipl;

import java.util.*;
import ch.psi.tomcat.tipl.*;

/** 
 The basis for the other image processing tools this class serves as a framework for storing, reading, and writing 3D image data. It is built using the memory model from OpenVMS and the framework of ImageJ so it is in theory 100% compatible with both approaches to image processing. 
 
 <li> Oct 10, 2011 - Fixed support for reading in float arrays from float named images
 <li> Dec 18, 2011 - Recoded to work on normal systems (windows, osx, linux) using tiff stacks
 <li> Jan 25, 2012 - Restructure class as ImageStack from ImageJ and added preview, ImagePlus and image processing capabilities
 */
public interface TImg {
    /** The function the reader uses to initialize an AIM */
    public boolean InitializeImage(D3int dPos, D3int cDim, D3int dOffset, D3float elSize,int imageType);
    /** The aim type of the image (0=char, 1=short, 2=int, 3=float, 10=bool, -1 same as input) */
	public int getImageType();
    /** The aim type of the image (0=char, 1=short, 2=int, 3=float, 10=bool, -1 same as input) */
    public void setImageType(int inData);
    /** The factor to scale bool/short/int/char values by when converting to/from float (for distance maps is  (1000.0/32767.0)) */
    public float getShortScaleFactor(); 
    public void setShortScaleFactor(float ssf);
    public float readShortScaleFactor(); 
	/** Is the image signed (should an offset be added / subtracted when the data is loaded to preserve the sign) */
	public boolean getSigned();
    public void setSigned(boolean inData);
    
    /** The size of the image */
    public D3int getDim();
    /** The size of the image */
    public void setDim(D3int inData);
    public int getSlices();
	/** The position of the bottom leftmost voxel in the image in real space, only needed for ROIs */
	public D3int getPos();
    /** The position of the bottom leftmost voxel in the image in real space, only needed for ROIs */
    public void setPos(D3int inData);
	/** The size of the border around the image which does not contain valid voxel data */
	public D3int getOffset();
    /** The size of the border around the image which does not contain valid voxel data */
    public void setOffset(D3int inData);
	/** The element size (in mm) of a voxel */
	public D3float getElSize();
    /** The element size (in mm) of a voxel */
    public void setElSize(D3float inData);
    public boolean CheckSizes(String otherPath);
	public boolean CheckSizes(TImg otherTImg);
    /** Is the data in good shape */
	public boolean isGood();
    /** The name of the data */
    public String getSampleName();
    /** The path of the data (whatever it might be) */
    public String getPath();
	/** Procedure Log, string containing past operations and information on the aim-file */
	public String getProcLog();
    public String appendProcLog(String inData);
    
    /** Whether or not basic compression (compatible almost everywhere) should be used when writing data */
	public boolean getCompression();
    public void setCompression(boolean inData);
    /** Load the full image into a cache so it can be quickly read (not applicable to all data types) */
    public void CacheFullImage();
    
    public Float[] GetXYZVec(int cIndex,int sliceNumber);
    public boolean[] getBoolArray(int sliceNumber);
    public char[] getByteArray(int sliceNumber);
    public short[] getShortArray(int sliceNumber);
    public boolean setShortArray(int sliceNumber,short[] cSlice);
    public int[] getIntArray(int sliceNumber);
    
    public float[] getFloatArray(int sliceNumber);
    
	
    

}
