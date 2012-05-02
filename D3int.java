/**
 This class represents an array of disk-resident images.
 Oct 10, 2011 - Fixed support for reading in float arrays from float named images
 Dec 18, 2011 - Recoded to work on normal systems (windows, osx, linux) using tiff stacks
 */
package ch.psi.tomcat.tipl;
public class D3int {
	public int x;
	public int y;
	public int z;
	public D3int(int xi,int yi,int zi) {
		setVals(xi,yi,zi);
	}
	public D3int(D3int xi) {
		setVals(xi.x,xi.y,xi.z);
	}
	public D3int(int xi) {
		setVals(xi,xi,xi);
	}
	public D3int() {
		setVals(0,0,0);
	}
	public void setVals(int xi,int yi,int zi) {
		x=xi;
		y=yi;
		z=zi;
	}
	public String toString() {
		return "("+x+","+y+","+z+")";
	}
	public double prod() {
		double out=x;
		out*=y;
		out*=z;
		return out;
	}
	
}
