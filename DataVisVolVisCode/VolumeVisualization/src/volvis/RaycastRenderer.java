/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 * @Anna
 * This class has the main code that generates the raycasting result image. 
 * The connection with the interface is already given.  
 * The different modes mipMode, slicerMode, etc. are already correctly updated
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private boolean mipMode = false;
    private boolean slicerMode = true;
    private boolean compositingMode = false;
    private boolean tf2dMode = false;
    private boolean shadingMode = false;
    private boolean mipXray = true;
    
    // Shading constants
    private float kAmbient = 0.0f;//0.1f;
    private float kDiff = 0.0f;
    private float kSpec = 1f;
    private float alpha = 10;
    private float iAmbient = 0.5f;
    private float iDiff = 1f;
    private float iSpec = 1f;
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        tFunc.setTestFunc();
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
     
    public void setShadingMode(boolean mode) {
        shadingMode = mode;
        changed();
    }
    
    public void setMIPMode() {
        setMode(false, true, false, false);
    }
    
    public void setSlicerMode() {
        setMode(true, false, false, false);
    }
    
    public void setCompositingMode() {
        setMode(false, false, true, false);
    }
    
    public void setTF2DMode() {
        setMode(false, false, false, true);
    }
    
    private void setMode(boolean slicer, boolean mip, boolean composite, boolean tf2d) {
        slicerMode = slicer;
        mipMode = mip;
        compositingMode = composite;
        tf2dMode = tf2d;        
        changed();
    }
    
        
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }
    
    /**
     * 
     * @param plane_pos
     * @param plane_normal
     * @param line_pos
     * @param line_dir
     * @param intersection
     * @return 
     * 
     * Find out whether the plane perpendicular to @plane_normal and going through @plane_pos intersects with the
     * line given by @line_pos and @line_dir. If true fill in intersection and return true, else return false.
     */
    private boolean intersectLinePlane(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection) {

        double[] tmp = new double[3];

        for (int i = 0; i < 3; i++) {
            tmp[i] = plane_pos[i] - line_pos[i];
        }

        double denom = VectorMath.dotproduct(line_dir, plane_normal);
        if (Math.abs(denom) < 1.0e-8) {
            return false;
        }

        double t = VectorMath.dotproduct(tmp, plane_normal) / denom;

        for (int i = 0; i < 3; i++) {
            intersection[i] = line_pos[i] + t * line_dir[i];
        }

        return true;
    }

    /**
     * 
     * @param intersection
     * @param xb
     * @param xe
     * @param yb
     * @param ye
     * @param zb
     * @param ze
     * @return 
     * 
     * Check whether the point @intersection falls into the cube specified by @xb, @xe, @yb, @ye, @zb and @ze.
     */
    private boolean validIntersection(double[] intersection, double xb, double xe, double yb,
            double ye, double zb, double ze) {

        return (((xb - 0.5) <= intersection[0]) && (intersection[0] <= (xe + 0.5))
                && ((yb - 0.5) <= intersection[1]) && (intersection[1] <= (ye + 0.5))
                && ((zb - 0.5) <= intersection[2]) && (intersection[2] <= (ze + 0.5)));

    }
    
    /**
     * 
     * @param plane_pos
     * @param plane_normal
     * @param line_pos
     * @param line_dir
     * @param intersection
     * @param entryPoint
     * @param exitPoint 
     * 
     * Find out whether the plane perpendicular to @plane_normal and going through @plane_pos intersects with the
     * line given by @line_pos and line_dir. If they intersect, @entryPoint or @exitPoint is filled. @entrypoint is filled
     * if @plane_normal and @line_dir have the same direction and @exitPoint is filled if @plane_normal and @line_dir have 
     * opposite direction.
     */
    private void intersectFace(double[] plane_pos, double[] plane_normal,
            double[] line_pos, double[] line_dir, double[] intersection,
            double[] entryPoint, double[] exitPoint) {

        boolean intersect = intersectLinePlane(plane_pos, plane_normal, line_pos, line_dir,
                intersection);
        if (intersect) {

            //System.out.println("Plane pos: " + plane_pos[0] + " " + plane_pos[1] + " " + plane_pos[2]);
            //System.out.println("Intersection: " + intersection[0] + " " + intersection[1] + " " + intersection[2]);
            //System.out.println("line_dir * intersection: " + VectorMath.dotproduct(line_dir, plane_normal));

            double xpos0 = 0;
            double xpos1 = volume.getDimX();
            double ypos0 = 0;
            double ypos1 = volume.getDimY();
            double zpos0 = 0;
            double zpos1 = volume.getDimZ();

            if (validIntersection(intersection, xpos0, xpos1, ypos0, ypos1,
                    zpos0, zpos1)) {
                if (VectorMath.dotproduct(line_dir, plane_normal) > 0) {
                    entryPoint[0] = intersection[0];
                    entryPoint[1] = intersection[1];
                    entryPoint[2] = intersection[2];
                } else {
                    exitPoint[0] = intersection[0];
                    exitPoint[1] = intersection[1];
                    exitPoint[2] = intersection[2];
                }
            }
        }
    }
    
  

    int traceRayMIP(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {
        VectorMath.normalize(viewVec);
        double totalDistance = VectorMath.distance(entryPoint, exitPoint);
        int maxValue = 0;
        // Raycast through the volume starting from @entryPoint and ending at @exitPoint along the direction
        // given by @viewVec with @sampleStep as sampling rate.
        for(double i = 0.0; true; i+=sampleStep){
            double[] coordinate = VectorMath.add(exitPoint, VectorMath.scale(viewVec, i));
            // If the sample point is outside of the volume break
            double distance = VectorMath.distance(exitPoint, coordinate);
            if (distance > totalDistance){
                break;
            }
            
            int currentValue = volume.getVoxelInterpolate(coordinate);
            if (currentValue > maxValue){
                maxValue = currentValue;
            }
        }
        
        int c_alpha, c_red, c_green, c_blue;
        // If xray mode is set, set the color as white
        if (mipXray){
            c_alpha = maxValue;
            c_red = 255;
            c_green = 255;
            c_blue = 255; 
        }
        else {
            // Get the color of @maxValue
            TFColor voxelColor = new TFColor(); 
            voxelColor = tFunc.getColor(maxValue);

            // BufferedImage expects a pixel color packed as ARGB in an int
            c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
            c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
            c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
            c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
        }

        int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
        return pixelColor;
    }
    
    

   
    void computeEntryAndExit(double[] p, double[] viewVec, double[] entryPoint, double[] exitPoint) {

        for (int i = 0; i < 3; i++) {
            entryPoint[i] = -1;
            exitPoint[i] = -1;
        }

        double[] plane_pos = new double[3];
        double[] plane_normal = new double[3];
        double[] intersection = new double[3];

        VectorMath.setVector(plane_pos, volume.getDimX(), 0, 0);
        VectorMath.setVector(plane_normal, 1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, -1, 0, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, volume.getDimY(), 0);
        VectorMath.setVector(plane_normal, 0, 1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, -1, 0);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, volume.getDimZ());
        VectorMath.setVector(plane_normal, 0, 0, 1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

        VectorMath.setVector(plane_pos, 0, 0, 0);
        VectorMath.setVector(plane_normal, 0, 0, -1);
        intersectFace(plane_pos, plane_normal, p, viewVec, intersection, entryPoint, exitPoint);

    }

    void raycast(double[] viewMatrix) {
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);


        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] entryPoint = new double[3];
        double[] exitPoint = new double[3];
        
        int increment=1;
        float sampleStep=0.2f;
        
        // Through experimentation, we decided to increase the sample step by five and we decrease the number of pixels by four, 
        // to find a good balance between two. 
        if (this.interactiveMode){
            increment = 2;
            sampleStep = 1f;
        }

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        for (int j = 0; j < image.getHeight(); j += increment) {
            for (int i = 0; i < image.getWidth(); i += increment) {
                // compute starting points of rays in a plane shifted backwards to a position behind the data set
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) - viewVec[0] * imageCenter
                        + volume.getDimX() / 2.0;
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) - viewVec[1] * imageCenter
                        + volume.getDimY() / 2.0;
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) - viewVec[2] * imageCenter
                        + volume.getDimZ() / 2.0;

                computeEntryAndExit(pixelCoord, viewVec, entryPoint, exitPoint);
                if ((entryPoint[0] > -1.0) && (exitPoint[0] > -1.0)) {
                    int pixelColor = 0;
                   if(mipMode) {
                        pixelColor= traceRayMIP(entryPoint,exitPoint,viewVec,sampleStep);
                   }
                   else if(compositingMode) {
                       pixelColor = composite(entryPoint,exitPoint,viewVec, sampleStep);
                   }
                   else if(tf2dMode) {
                       pixelColor = tranferFuntion2D(entryPoint,exitPoint,viewVec, sampleStep);
                   }
                      
                    for (int ii = i; ii < i + increment; ii++) {
                        for (int jj = j; jj < j + increment; jj++) {
                            image.setRGB(ii, jj, pixelColor);
                        }
                    }
                }

            }
        }
    }
    
    int tranferFuntion2D(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {
        //Get the user selected values
        int selectedIntensity = tfEditor2D.triangleWidget.baseIntensity;
        double selectedRadius = tfEditor2D.triangleWidget.radius;
        TFColor selectedColor = tfEditor2D.triangleWidget.color;
        
        VectorMath.normalize(viewVec);
        double totalDistance = VectorMath.distance(entryPoint, exitPoint);
        TFColor currentColor = new TFColor(); 
        TFColor previousColor = new TFColor(0,0,0,0);
        // Raycast through the volume starting from @exitPoint and ending at @entryPoint along the direction
        // given by @viewVec with @sampleStep as sampling rate.
        // We are sampeling back to front here.
        for(double i = 0.0; true; i+=sampleStep){
            double[] coordinate = VectorMath.add(exitPoint, VectorMath.scale(viewVec, i));
            
            // If this sample is beyond the exitPoint point break
            double distance = VectorMath.distance(exitPoint, coordinate);
            if (distance > totalDistance){
                break;
            }
            int currentValue = volume.getVoxelInterpolate(coordinate);
            VoxelGradient gradient = gradients.getGradient(coordinate); 
            currentColor.r = selectedColor.r;
            currentColor.g = selectedColor.g;
            currentColor.b = selectedColor.b;
            
            // Use Levoy's formula to set the opacity of the current sample
            if (currentValue == selectedIntensity && gradient.mag == 0) {
                currentColor.a = selectedColor.a * 1.0;
            }
            else if (gradient.mag > 0.0 && 
                    ((currentValue - selectedRadius * gradient.mag) <= selectedIntensity) &&
                    ((currentValue + selectedRadius * gradient.mag) >= selectedIntensity))  {
                currentColor.a = selectedColor.a*(1.0 - (1 / selectedRadius) * (Math.abs((selectedIntensity - currentValue)/ gradient.mag)));
            }
            else {
                currentColor.a = 0.0;
            }
            if (shadingMode) 
                getShade(viewVec, gradient, selectedColor, currentColor);

            previousColor.r = currentColor.a*currentColor.r + (1-currentColor.a) * previousColor.r;
            previousColor.g = currentColor.a*currentColor.g + (1-currentColor.a) * previousColor.g;
            previousColor.b = currentColor.a*currentColor.b + (1-currentColor.a) * previousColor.b;
        }
        previousColor.a = 1.0;
        return getColor(previousColor);
    }  
    
    void getShade(double[] viewVector, VoxelGradient gradient, TFColor selectedColor, TFColor resultingColor){
        // We are taking the direction vector from the point on the suface 
        // toward the light as the same as the viewVector
        double[] lightDir = viewVector;        
        // Normal at the point on the surface
        double[] normal = gradient.normalize();
        // Compute the direction that a perfectly reflected ray of light would take form the surface
        double[] perfectReflectionDir = VectorMath.sub(VectorMath.scale(normal, 2 * (VectorMath.dotproduct(lightDir, normal))), lightDir) ;
        
        double normalLightDot = VectorMath.dotproduct(lightDir, normal);
        
        // The constants are defined at the beginning of the file
        double K_a = kAmbient*iAmbient;
        double K_d = kDiff*iDiff* normalLightDot;
        double K_s = kSpec*iSpec* Math.pow(VectorMath.dotproduct(perfectReflectionDir, viewVector), alpha);
        
        // If the angle between the normal and the direction of the light greater than 180 degrees, then no light falls on the surface
        // In this case the opacity of this point is set to zero.
        if (normalLightDot <= 0){
            resultingColor.a = 0;
        }
        else{
            resultingColor.r = selectedColor.r * K_a + selectedColor.r * K_d + K_s;
            resultingColor.g = selectedColor.g * K_a + selectedColor.g * K_d + K_s;
            resultingColor.b = selectedColor.b * K_a + selectedColor.b * K_d + K_s; 
  
        }
   }
    
    int getColor(TFColor color){
        // BufferedImage expects a pixel color packed as ARGB in an int
        int c_alpha = color.a <= 1.0 ? (int) Math.floor(color.a * 255) : 255;
        int c_red = color.r <= 1.0 ? (int) Math.floor(color.r * 255) : 255;
        int c_green = color.g <= 1.0 ? (int) Math.floor(color.g * 255) : 255;
        int c_blue = color.b <= 1.0 ? (int) Math.floor(color.b * 255) : 255;
        return (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
    }
    
    int composite(double[] entryPoint, double[] exitPoint, double[] viewVec, double sampleStep) {

        VectorMath.normalize(viewVec);
        double totalDistance = VectorMath.distance(entryPoint, exitPoint);
        TFColor currentColor; 
        TFColor previousColor = new TFColor();
        // Raycast through the volume starting from @exitPoint and ending at @entryPoint along the direction
        // given by @viewVec with @sampleStep as sampling rate.
        // We are sampeling back to front here.
        for(double i = 0.0; true; i+=sampleStep){
            double[] coordinate = VectorMath.add(exitPoint, VectorMath.scale(viewVec, i));
            double distance = VectorMath.distance(exitPoint, coordinate);
            
            if (distance > totalDistance){
                break;
            }
            
            int currentValue = volume.getVoxelInterpolate(coordinate);
            currentColor = tFunc.getColor(currentValue);
            previousColor.r = currentColor.a*currentColor.r + (1-currentColor.a) * previousColor.r;
            previousColor.g = currentColor.a*currentColor.g + (1-currentColor.a) * previousColor.g;
            previousColor.b = currentColor.a*currentColor.b + (1-currentColor.a) * previousColor.b;
        }
        previousColor.a = 1.0;

        return getColor(previousColor);
    }    
    
    /**
     * 
     * @param viewMatrix: This varaible contains the direction of the viewpoint, and the direction of the plane
     *                    perpendicular to it.
     * Visualize a slice throught the middelpoint of the volume. This plane has the same direction as the plane
     * perpendular to the direction of the viewpoint.
     * 
     */
    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = volume.getVoxelInterpolate(pixelCoord);
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                
                // Alternatively, apply the transfer function to obtain a color
                TFColor auxColor = new TFColor(); 
                auxColor = tFunc.getColor(val);

                image.setRGB(i, j, getColor(voxelColor));
            }
        }


    }


    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        if (slicerMode) {
            slicer(viewMatrix);    
        } else {
            raycast(viewMatrix);
        }
        
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
