import ij.*;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import java.awt.*;
import java.util.Objects;

/**
 * This plugin is the implemenation of "SynQuant: An Automatic Tool to Quantify Synapses from Microscopy Images". It 
 * 1. segments synapse/puncta on 2d or 3d image(s) using Probability principled method
 * 2. extracts dendrite from 2d image and use steerable filter
 * 3. quantifies the synapse detection and dendrite extraction results if input data are all 2d
 * 
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @version 1.2
 * @date 2019-06-17
 *
 *
 */

public class SynQuantSimple implements PlugIn, DialogListener{
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp = null; // Current Image
	protected ImagePlus[] impVec; // Current Image
	protected int width; //image width
	protected int height; //image height
	protected int timePts; //image height
	protected int zSlice; //image height
	protected int numChannels; //image height
	protected double fdr = 0.05; // fdr control threshold 
	protected double zscore_thres = 0;
	protected int MinSize,MaxSize; // synapse size range
	protected double minFill, maxWHRatio;
	protected ImagePlus outputImp=null; // detection results output
	protected ImagePlus outputImp1=null; // detection results output
	int [][][][] synIdx;
	double [][][][] synZscore;
	double slideThrZ;
	int [][][] sliderSynMap; 	// synapse map after post-processing
	int way2combinePostPre = 1; // the way to combine pre- and post-channel, 0:intersect, 1: add pre-zscore to post-zscore
	double ExtendedDistance=0; // extended distance
	double zAxisMultiplier=1; // z axis extended distance multiplier
	
	/***Show the dialog for parameter input, then start synQuant***/
	public void run(String arg) {
		try {
			if (showDialog()) {
				synQuant3D_real();
			}
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}

	public boolean showDialog() 
	{
		// Get input parameter
		GenericDialog gd = new GenericDialog("3D Particles - Data and Parameter Setting");
		gd.addNumericField("Z-Score for Particle Detection: ", 10,  2);//2.5-3.5 good
		gd.addNumericField("Min Particle Size: ",10, 0);//2.5-3.5 good
		gd.addNumericField("Max Particle Size: ", 200, 0);//2.5-3.5 good
		gd.addNumericField("Min fill: ", 0.5, 2);//2.5-3.5 good
		gd.addNumericField("Max WH Ratio: ", 4, 0);//2.5-3.5 good
		gd.addNumericField("z distance multiplier: ", 1, 2);// 1, treat z distance the same as x-y distance
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		zscore_thres = gd.getNextNumber();
		MinSize = (int) gd.getNextNumber();
		MaxSize = (int) gd.getNextNumber();
		minFill = gd.getNextNumber();
		maxWHRatio = gd.getNextNumber();
		if(zscore_thres < 0) {//Double.isNaN(fdr) || (fdr<=0) || (fdr>=1)){
			IJ.showMessage("Invalid zscore parameter(s).\n" + "positive value needed.");
			return false;
		}
		zAxisMultiplier = gd.getNextNumber();
		numChannels = 1; // we only care two channels: post- and pre-synaptic channel
		impVec = new ImagePlus[1]; // we only care two channels: pre- and post-synaptic channe
		impVec[0] = WindowManager.getCurrentImage();

		return true;
	}

	/***
	 * The main function of synQuant.
	 * Do synapse detection channel by channel
	 * ***/
	public void synQuant3D_real() {
		//// parameter initialization
		paraQ3D q = new paraQ3D(numChannels, way2combinePostPre, 0.8);
		
		q.ExtendedDistance=ExtendedDistance;
		q.zAxisMultiplier=zAxisMultiplier;
		
		BasicMath bm = new BasicMath();
		//// data saving final results
		slideThrZ = 1000;
		ppsd3D particle3D_det = null;
        for (ImagePlus imagePlus : impVec) {
            if (imagePlus == null)
                continue;
            imp = imagePlus;
            stack = imp.getStack();
            timePts = imp.getNFrames();
            //numChannels = imp.getNChannels();
            zSlice = imp.getNSlices();
            width = imp.getWidth();
            height = imp.getHeight();
            if (q.NumChannelProcessed == 0) {
                q.synZscore = new double[timePts][zSlice][height][width];
                synIdx = new int[timePts][zSlice][height][width];
                synZscore = new double[timePts][zSlice][height][width];
            }

            double vox_x = imp.getCalibration().pixelWidth;
            if (vox_x == 1)//simulated data
                vox_x = 2.0757e-7;
            else
                vox_x = vox_x * 1e-6;//real data

            int type = imp.getType();
            long startTime1 = System.nanoTime();
            ////particle detection
            for (int i = 1; i <= timePts; i++) {
                q.curTps = i - 1;
                short[][] Arr3D = stack2array(type, stack, i); // #zstack*#pixels in one slice
                paraP3D p = new paraP3D(fdr, zscore_thres, (int) bm.matrix2DMin(Arr3D), (int) bm.matrix2DMax(Arr3D), MinSize, MaxSize, minFill, maxWHRatio);
                particle3D_det = new ppsd3D(Arr3D, width, height, vox_x, p, q);

                // for the post channel, we save the output results
                if (q._NumChannel == q.NumChannelProcessed + 1) { // last one
                    synIdx[i - 1] = particle3D_det.ppsd_main.kMap;
                    synZscore[i - 1] = particle3D_det.ppsd_main.zMap;
                    double tmpThrZ = particle3D_det.ppsd_main.thrZ;
                    if (slideThrZ > tmpThrZ)
                        slideThrZ = tmpThrZ;
                } else {// pre-channel
                    q.synZscore[q.curTps] = particle3D_det.ppsd_main.zMap;
                }
            }
            long endTime1 = System.nanoTime();
            System.out.println("Finished. Data size: " + zSlice + " * " + height + " * " + width + " * " + timePts + " Total running time: " + (endTime1 - startTime1) / 1e9);
            q.NumChannelProcessed++;
            q.var = 0; // reset variance, preparing for new channel
        }

		//// display synapse detection results: RGB or single channel
		outputImp = IJ.createHyperStack("Synapse detection results", width, height, 1, zSlice, timePts,24/*bitdepth*/);
		outputImp.show();
		for (int i=1; i <= timePts; i++){
			//display data
			SynapticDisplay(synIdx[i-1], i);
			System.out.println(i + "-th Frame SynNum: "+bm.matrix3DMax(synIdx[i-1]));
			outputImp.updateAndDraw();
		}
		//// wait to listen the change of z-score threshold
		GenericDialog gd = new NonBlockingGenericDialog("Tune zscore threshold");
		gd.addSlider("zscore threshold tuning", Math.max(0, slideThrZ), 100, Math.max(0, slideThrZ));
		// wait to listen for the changing of zscore threshold
		sliderSynMap = null;
		gd.addDialogListener(this);
		gd.showDialog();
		outputImp.close();

		boolean [][][] synMap3dBin = new boolean[synZscore[0].length][synZscore[0][0].length][synZscore[0][0][0].length];
		if (sliderSynMap == null)
		{
			for (int k = 0; k<synZscore[0].length; k++)
			{
				for (int i=0; i<synZscore[0][0].length; i++)
				{
					for (int j=0; j<synZscore[0][0][0].length;j++)
					{
                        synMap3dBin[k][i][j] = synZscore[0][k][i][j] >= zscore_thres;
					}
				}
			}
		}
		else
		{
			for (int k = 0; k<sliderSynMap.length; k++)
			{
				for (int i=0; i<sliderSynMap[0].length; i++)
				{
					for (int j=0; j<sliderSynMap[0][0].length;j++)
					{
                        synMap3dBin[k][i][j] = sliderSynMap[k][i][j] > 0;
					}
				}
			}
		}

		outputImp1 = IJ.createHyperStack("Synapse mask", width, height, 1, zSlice, 1, 8/*bitdepth*/);
		outputImp1.show();
		ImageProcessor outIP =  null;
		ImageStack curStack = outputImp1.getStack();
		for(int k=1; k<=zSlice; k++) {
			outIP = curStack.getProcessor(k);
			for (int i = 0; i<width; i++) {
				for (int j = 0; j<height; j++) {
					int x0 = synMap3dBin[k-1][j][i] ? 255 : 0;
					outIP.set(i, j, x0);
				}
			}
		}
		System.out.println("Done");
	}

	/** Listener to modifications of the input fields of the dialog.
	 *  Here the parameters should be read from the input dialog.
	 *  @param gd The GenericDialog that the input belongs to
	 *  @param e  The input event
	 *  @return whether the input is valid and the filter may be run with these parameters
	 */
	public boolean dialogItemChanged (GenericDialog gd, AWTEvent e) {
		zscore_thres = gd.getNextNumber();
		if (Math.abs(zscore_thres-slideThrZ) < 0.2)
			return true;
		else
			slideThrZ = zscore_thres;
		outputImp.show();
		for (int i=1; i <= timePts; i++)
		{
			sliderSynMap = new int [zSlice][height][width];
			for (int zz=0;zz<zSlice; zz++) 
			{
				for(int yy=0;yy<height; yy++) 
				{
					for(int xx=0;xx<width; xx++) 
					{
						if(synZscore[i-1][zz][yy][xx] >= zscore_thres) 
						{
							sliderSynMap[zz][yy][xx] = 1;//synZscore[i-1][zz][yy][xx];
						}
					}
				}
			}
			//int [][][] tmpSynIdx = imh.bwlabel3D(tmpSynZ, 26);
			SynapticDisplay(sliderSynMap, i);
			outputImp.updateAndRepaintWindow();
		}
		return true;
	}

	public void SynapticDisplay(int [][][] kSynR1, int curTimePt){
		ImageStack curStack = outputImp.getStack();
		ImageStack curInStack = imp.getStack();
		ImageProcessor inIP = null;
		ImageProcessor outIP =  null;

		double max_intensity = 0;
		if(imp.getType() == ImagePlus.GRAY16) 
		{
			for(int stackNum=1; stackNum<=curStack.getSize(); stackNum++) 
			{
				double tmp_max = curInStack.getProcessor(stackNum).getStatistics().max;
				if (max_intensity < tmp_max) 
				{
					max_intensity = tmp_max;
				}
			}
		}

		for(int stackNum=1; stackNum<=curStack.getSize(); stackNum++) 
		{
			outIP = curStack.getProcessor(stackNum);
			inIP = curInStack.getProcessor(stackNum);
			for (int i = 0; i < width; i++) 
			{
				for(int j = 0;j<height;j++)
				{
					// put channel values in an integer
					if(imp.getType() == ImagePlus.GRAY8) 
					{
						int tmp_val = (int)inIP.get(i, j) & 0xff;
						if(kSynR1[stackNum-1][j][i]>0)
							outIP.set(i, j, (((int) 200 & 0xff) << 16) + (tmp_val << 8));
						else
							outIP.set(i, j, (tmp_val << 8));
					}else 
					{
						byte tmp_float_val = (byte) ((((double)inIP.get(i, j))/max_intensity)*255);
						int tmp_val = (int)tmp_float_val & 0xff;
						if(kSynR1[stackNum-1][j][i]>0)
							outIP.set(i, j, (((int) 200 & 0xff) << 16) + (tmp_val << 8));
						else
							outIP.set(i, j, (tmp_val << 8));
					}
				}
			}
		}
	}

	public short[][] stack2array(int type, ImageStack stack, int frameNum){
		int mask;
		int nPixels=width*height;
		short[][] imArray = new short[zSlice][nPixels];
		if (type == ImagePlus.GRAY16)
		{
			mask=0xffff;
			for(int zz=1;zz<=zSlice;zz++) {
				int curSliceNum = imp.getStackIndex(1, zz, frameNum);// default all data are one-channel data
				short[] pixels = (short[])stack.getPixels(curSliceNum);
				int intP = (int)(mask&pixels[0]);
				for (int i=0; i<nPixels; i++)
				{
					intP=(int)(mask&pixels[i]);
					short p = (short)(intP/2); // for sake of range, short <= 32,767 (65535/2)
					imArray[zz-1][i]=p;
				}
			}
		}  
		else if (type == ImagePlus.GRAY8) 
		{
			mask=0xff;
			for(int zz=1;zz<=zSlice;zz++) {
				int curSliceNum = imp.getStackIndex(1, zz, frameNum); // default all data are one-channel data
				byte[] pixels = (byte[])stack.getPixels(curSliceNum);
				for (int i=0; i<nPixels; i++)
				{
					short p=(short)(mask&pixels[i]);
					imArray[zz-1][i]=p;
				}
			}
		}
		else
		{
			IJ.log("Pixel format not supported");
			return null;
		}
		return imArray;
	}

	/**
	 * Main method for debugging.
	 * <p>
	 * For debugging, it is convenient to have a method that starts ImageJ,
	 * loads an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args
	 *            unused
	 */
	public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins
		// menu
		Class<?> clazz = SynQuantSimple.class;
		String url = Objects.requireNonNull(clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class")).toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		ImagePlus image = IJ.openImage("C:\\Users\\yzwang\\proj\\Java\\data\\example synapse.tif");
		//C2-3-weak-z_Maximum intensity projection.tif");
		image.show();

		 //run the plugin
		IJ.runPlugIn(clazz.getName(), "");
	}

	void error() {
		IJ.showMessage("3D particle", "Error");
	}
}
