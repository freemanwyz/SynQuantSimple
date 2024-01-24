import ij.*;
import ij.plugin.PlugIn;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.process.ImageProcessor;

import java.util.Objects;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

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

public class SynQuantBatch implements PlugIn {
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp = null; // Current Image
	protected ImagePlus[] impVec; // Current Image
	protected int width; //image width
	protected int height; //image height
	protected int timePts; //image height
	protected int zSlice; //image height
	protected double fdr = 0.05; // fdr control threshold 
	protected double zscore_thres = 0;
	protected int MinSize,MaxSize; // synapse size range
	protected double minFill, maxWHRatio;
	protected ImagePlus outputImp1=null; // detection results output
	int [][][][] synIdx;
	double [][][][] synZscore;
	double zAxisMultiplier=1; // z axis extended distance multiplier
	
	/***Show the dialog for parameter input, then start synQuant***/
	public void run(String arg) {
		zscore_thres = 20;
		MinSize = 10;
		MaxSize = 200;
		minFill = 0.5;
		maxWHRatio = 4;
		zAxisMultiplier = 1;

		String file_name = "";
		file_name = Macro.getOptions();

		// C:\\Users\\yzwang\\proj\\Java\\data\\param.txt
		if(file_name == null) {
			GenericDialog gd = new GenericDialog("Specify configuration file");
			gd.addStringField("Filename", "", 20);
			gd.showDialog();
			file_name = gd.getNextString();
		}

		System.out.println("SynQuant");
		System.out.println(file_name);

		BufferedReader reader;
		try {
			if(!file_name.isEmpty()) {
				reader = new BufferedReader(new FileReader(file_name));
				String line = reader.readLine();
				while (line != null) {
					line = line.replaceAll("\\s+","");
					String[] arrOfStr = line.split("=", 2);
					if (Objects.equals(arrOfStr[0], "zscore_thres")) {
						zscore_thres = Float.parseFloat(arrOfStr[1]);
					}
					if (Objects.equals(arrOfStr[0], "MinSize")) {
						MinSize = Integer.parseInt(arrOfStr[1]);
					}
					if (Objects.equals(arrOfStr[0], "MaxSize")) {
						MaxSize = Integer.parseInt(arrOfStr[1]);
					}
					if (Objects.equals(arrOfStr[0], "minFill")) {
						minFill = Float.parseFloat(arrOfStr[1]);
					}
					if (Objects.equals(arrOfStr[0], "maxWHRatio")) {
						maxWHRatio = Float.parseFloat(arrOfStr[1]);
					}
					if (Objects.equals(arrOfStr[0], "zAxisMultiplier")) {
						zAxisMultiplier = Float.parseFloat(arrOfStr[1]);
					}
					line = reader.readLine();
				}
				reader.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}

		System.out.printf("%f, %d, %d, %f, %f, %f\n", zscore_thres, MinSize, MaxSize, minFill, maxWHRatio, zAxisMultiplier);

		impVec = new ImagePlus[1]; // we only care two channels: pre- and post-synaptic channe
		impVec[0] = WindowManager.getCurrentImage();

		try {
			synQuant3D_real();
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}

	/***
	 * The main function of synQuant.
	 * Do synapse detection channel by channel
	 * ***/
	public void synQuant3D_real() {
		//// parameter initialization
		paraQ3D q = new paraQ3D(1, 1, 0.8);
		
		q.ExtendedDistance=0;
		q.zAxisMultiplier=zAxisMultiplier;
		
		BasicMath bm = new BasicMath();
		//// data saving final results
		ppsd3D particle3D_det = null;
        for (ImagePlus imagePlus : impVec) {
            if (imagePlus == null)
                continue;
            imp = imagePlus;
            stack = imp.getStack();
            timePts = imp.getNFrames();
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
		boolean [][][] synMap3dBin = new boolean[synZscore[0].length][synZscore[0][0].length][synZscore[0][0][0].length];
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
		Class<?> clazz = SynQuantBatch.class;
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

}
