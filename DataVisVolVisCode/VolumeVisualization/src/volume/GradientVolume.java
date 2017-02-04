/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 * @ Anna
 * This class contains the pre-computes gradients of the volume. This means calculates the gradient
 * at all voxel positions, and provides functions
 * to get the gradient at any position in the volume also continuous..
*/
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }
    
    public VoxelGradient getGradientNN(double[] coord) {
        /* Nearest neighbour interpolation applied to provide the gradient */
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
                || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
        }

        int x = (int) Math.round(coord[0]);
        int y = (int) Math.round(coord[1]);
        int z = (int) Math.round(coord[2]);
        
        return getGradient(x, y, z);
    }
    
    private VoxelGradient interpolate(VoxelGradient g0, VoxelGradient g1, float factor) {
        return new VoxelGradient(g1.x * factor + g0.x * (1-factor),
                                 g1.y * factor + g0.y * (1-factor),
                                 g1.z * factor + g0.z * (1-factor));
    }

    
    public VoxelGradient getGradient(double[] coord) {
    /* To be implemented: Returns trilinear interpolated gradient based on the precomputed gradients. 
     *   Use function interpolate. Use getGradientNN as bases */
        if (coord[0] < 0 || coord[0] > (dimX-2) || coord[1] < 0 || coord[1] > (dimY-2)
            || coord[2] < 0 || coord[2] > (dimZ-2)) {
            return zero;
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
        float x_ratio = (float)(coord[0] - x0);
        float y_ratio = (float)(coord[1] - y0);
        float z_ratio = (float)(coord[2] - z0);
        
        VoxelGradient c00 = interpolate(getGradient((int) x0, (int) y0, (int) z0),
                                        getGradient((int) x1, (int) y0, (int) z0),
                                        x_ratio);
        VoxelGradient c01 = interpolate(getGradient((int) x0, (int) y0, (int) z1),
                                        getGradient((int) x1, (int) y0, (int) z1),
                                        x_ratio);
        VoxelGradient c10 = interpolate(getGradient((int) x0, (int) y1, (int) z0),
                                        getGradient((int) x1, (int) y1, (int) z0),
                                        x_ratio);
        VoxelGradient c11 = interpolate(getGradient((int) x0, (int) y1, (int) z1),
                                        getGradient((int) x1, (int) y1, (int) z1),
                                        x_ratio);
        VoxelGradient c0 = interpolate(c00, c10, y_ratio);
        VoxelGradient c1 = interpolate(c01,c11, y_ratio);
        
        return interpolate(c0, c1, z_ratio);
    }
  
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
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

    private void compute() {
        for (int i = 0; i < data.length; i++) {
            data[i] = zero;
        }
        this.maxGradientMagnitude = 0;
        for (int i = 1; i < volume.getDimX() - 1; i++) {
            for (int j = 1; j < volume.getDimY() - 1; j++) {
                for (int k = 1; k < volume.getDimZ() - 1; k++) {
                    float gx = ((volume.getVoxel(i - 1, j, k) - volume.getVoxel(i + 1, j, k)) / 2.0f);
                    float gy = ((volume.getVoxel(i, j - 1, k) - volume.getVoxel(i, j + 1, k)) / 2.0f);
                    float gz = ((volume.getVoxel(i, j, k - 1) - volume.getVoxel(i, j, k + 1)) / 2.0f);

                    VoxelGradient value = new VoxelGradient(gx, gy, gz);
                    if (value.mag > this.maxGradientMagnitude) {
                        this.maxGradientMagnitude = value.mag;
                    }
                    this.setGradient(i, j, k, value);
                }
            }
        }
    }
    
    public double getMaxGradientMagnitude() {
        /* to be implemented: Returns the maximum gradient magnitude*/
        return this.maxGradientMagnitude;
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
    double maxGradientMagnitude;

}
