/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 * @ Anna
 *  This class contains one vector as its x,y,z coordinate and includes the magnitude. 
 *  This will the voxel element in a gradient volume 
 */
public class VoxelGradient {

    public float x, y, z;
    public float mag;
    
    public VoxelGradient() {
        x = y = z = mag = 0.0f;
    }
    
    public VoxelGradient(float gx, float gy, float gz) {
        x = gx;
        y = gy;
        z = gz;
        mag = (float) Math.sqrt(x*x + y*y + z*z);
    }
    
    public double[] normalize(){
        double[] result = new double[] {x,y,z};
        if (mag != 0){
            result[0] = x/mag;
            result[1] = y/mag;
            result[2] = z/mag;
        }
        return result;
    }
}
