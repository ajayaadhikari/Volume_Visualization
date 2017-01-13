/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

import java.io.File;
import java.io.IOException;

/**
 *
 * @author michel
 * @Anna 
 * Volume object: This class contains the object and assumes that the distance between the voxels in x,y and z are 1 
 */
public class Volume {
    
    public Volume(int xd, int yd, int zd) {
        data = new short[xd*yd*zd];
        dimX = xd;
        dimY = yd;
        dimZ = zd;
    }
    
    public Volume(File file) {
        
        try {
            VolumeIO reader = new VolumeIO(file);
            dimX = reader.getXDim();
            dimY = reader.getYDim();
            dimZ = reader.getZDim();
            data = reader.getData().clone();
            computeHistogram();
        } catch (IOException ex) {
            System.out.println("IO exception");
        }
        
    }
    
    
    public short getVoxel(int x, int y, int z) {
        return data[x + dimX*(y + dimY * z)];
    }
    
    public void setVoxel(int x, int y, int z, short value) {
        data[x + dimX*(y + dimY*z)] = value;
    }

    public void setVoxel(int i, short value) {
        data[i] = value;
    }
    
    public short getVoxelInterpolate(double[] coord) {
    /* to be implemented: get the trilinear interpolated value. 
        The current implementation gets the Nearest Neightbour */
        
        if (coord[0] < 0 || coord[0] > (dimX-1) || coord[1] < 0 || coord[1] > (dimY-1)
                || coord[2] < 0 || coord[2] > (dimZ-1)) {
            return 0;
        }
        
        // Compute the ceil and floor of the given coordinates in all dimensions
        double x0 = Math.floor(coord[0]);
        double x1 = Math.ceil(coord[0]);
        double y0 = Math.floor(coord[1]);
        double y1 = Math.ceil(coord[1]);
        double z0 = Math.floor(coord[2]);
        double z1 =  Math.ceil(coord[2]);
        
        // Compute the distance from the floor of the coordinate (d1) to the coordinate in all dimensions
        // The distance between the surrounding voxels (d2) for each dimension is 1.
        // So it is equal to the ratio between d1 and d2,
        double x_ratio = (coord[0] - x0);
        double y_ratio = (coord[1] - y0);
        double z_ratio = (coord[2] - z0);
        
        double c00 = getVoxel((int) x0, (int) y0, (int) z0) * (1 - x_ratio) +
                     getVoxel((int) x1, (int) y0, (int) z0) * x_ratio;
        double c01 = getVoxel((int) x0, (int) y0, (int) z1) * (1 - x_ratio) +
                     getVoxel((int) x1, (int) y0, (int) z1) * x_ratio;
        double c10 = getVoxel((int) x0, (int) y1, (int) z0) * (1 - x_ratio) +
                     getVoxel((int) x1, (int) y1, (int) z0) * x_ratio;
        double c11 = getVoxel((int) x0, (int) y1, (int) z1) * (1 - x_ratio) +
                     getVoxel((int) x1, (int) y1, (int) z1) * x_ratio;
        
        double c0 = c00 * (1 - y_ratio) + c10 * y_ratio;
        double c1 = c01 * (1 - y_ratio) + c11 * y_ratio;
        
        double result = c0 * (1 - z_ratio) + c1 * z_ratio;
        
        return (short)result;
        
    }
    
    public short getVoxel(int i) {
        return data[i];
    }
    
    public int getDimX() {
        return dimX;
    }
    
    public int getDimY() {
        return dimY;
    }
    
    public int getDimZ() {
        return dimZ;
    }

    public short getMinimum() {
        short minimum = data[0];
        for (int i=0; i<data.length; i++) {
            minimum = data[i] < minimum ? data[i] : minimum;
        }
        return minimum;
    }

    public short getMaximum() {
        short maximum = data[0];
        for (int i=0; i<data.length; i++) {
            maximum = data[i] > maximum ? data[i] : maximum;
        }
        return maximum;
    }
 
    public int[] getHistogram() {
        return histogram;
    }
    
    private void computeHistogram() {
        histogram = new int[getMaximum() + 1];
        for (int i=0; i<data.length; i++) {
            histogram[data[i]]++;
        }
    }
    
    private int dimX, dimY, dimZ;
    private short[] data;
    private int[] histogram;
}
