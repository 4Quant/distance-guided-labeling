/**
 This class represents an array of disk-resident images.
 Oct 10, 2011 - Fixed support for reading in float arrays from float named images
 Dec 18, 2011 - Recoded to work on normal systems (windows, osx, linux) using tiff stacks
 */
package ch.psi.tomcat.tipl;
public class D3float {
	// called a float on VMS but for other systems double is probably more reliable
	public double x=0.0;
	public double y=0.0;
	public double z=0.0;
	public D3float(double xi,double yi,double zi) {
		setVals(xi,yi,zi);
	}
	public D3float(D3float xi) {
		setVals(xi.x,xi.y,xi.z);
	}
	public D3float() {
		setVals(0.0,0.0,0.0);
	}
    public void setVals(double xi,double yi,double zi) {
		x=xi;
		y=yi;
		z=zi;
	}
	public String toString() {
		return "("+String.format("%.4f",x)+","+String.format("%.4f",y)+","+String.format("%.4f",z)+")";
	}
}