

/**
 * This plugin 
 * 1. segments synapse/puncta using Probability principled synapse detection
 * algorithm modified from "PPSD: Probability Principled Synapse Detection"
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @version 1.0
 * @date 2016-05-31
 *
 * Update: 
 * version 1.1. Add a fast version based on component tree; 
 * Major idea of component tree building and zscore calculation is from Dr.Petter Ranefall(petter.ranefall@it.uu.se).
 * version 1.2. Add one more noise estimation and stabilization method
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @date 2019-06-17
 *
 */

public class ppsd3D {
	protected int zSlice; //image width
	protected int width; //image width
	protected int height; //image height
	public fastppsdcore3D ppsd_main;

	public ppsd3D(short[][] imArray, int inwidth, int inheight, double vox_x, paraP3D p, paraQ3D q) {
		
		/*** STEP 1. Data preparation*******/
		width = inwidth;
		height = inheight;
		zSlice = imArray.length;
		int ntry = 1;
		q.init(imArray, vox_x, ntry, height, width);
		
		double[][][] G = new double[zSlice][height][width];
		double[][][] Gt = new double[zSlice][height][width];
		// boolean isDebug = java.lang.management.ManagementFactory.getRuntimeMXBean(). getInputArguments().toString().contains("-agentlib:jdwp");
		//// put all zstack in a single frame 
		double[][] tmpG = new double[height][width*zSlice];
		for (int zz = 0; zz<zSlice; zz++) {
			int shiftPixels = zz*width;
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					tmpG[i][j+shiftPixels] = (double) imArray[zz][i * width + j] * 255 / q.maxVal;// the data is saved by line
					G[zz][i][j] = (double) imArray[zz][i * width + j] * 255 / q.maxVal;
				}
			}
		}

		/*** STEP 2. estimate variance and do stabilization ********/
		long startTime1=System.nanoTime();     // Testing running time: start
		ImageHandling IM = new ImageHandling();
		double [][] tmpGt;
		if (q.var == 0) { // the first frame
			// the first frame, we use 3 stacks to estimate noise first
			double [][] stack4noiseEst;
			if(zSlice <= 3) // use all stacks
			{
				stack4noiseEst = new double[height][width];
				
				for (int i = 0; i < height; i++) {
					for (int j = 0; j < width; j++) {
						stack4noiseEst[i][j] = G[zSlice/2][i][j];
					}
				}
			}
			else {
				stack4noiseEst = new double[height][width*3];
				for (int i = 0; i < height; i++) {
					for (int j = 0; j < width; j++) {
						stack4noiseEst[i][j] = G[0][i][j];
						stack4noiseEst[i][j+width] = G[zSlice/2][i][j];
						stack4noiseEst[i][j+width*2] = G[zSlice-1][i][j];
					}
				}
			}
			IM.getNoiseRobust(stack4noiseEst,  q); // the noise distribution has been saved in q
			q.var = q.varRatioBased; // directly use a estimated noise variance
		}
		/***choice 1: do stabilization***/
		//tmpGt = IM.getNoiseVSTwithAlpha(tmpG,  q);
		/***choice 2: no stabilization***/
		tmpGt = tmpG;
		
		long endTime1=System.nanoTime();
		System.out.println("Running time: "+(endTime1-startTime1)/1e9); 
		//change the tmpGt to Gt and imArray
		for (int zz = 0; zz<zSlice; zz++) {
			int shiftPixels = zz*width;
			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					Gt[zz][i][j] = tmpGt[i][j+shiftPixels];// the data is saved by line, what if we use orginal data?
					imArray[zz][i * width + j] = (short) Math.round(Gt[zz][i][j]);
				}
			}
		}
		/*** STEP 3. start ppsd for synapse detection *
		 * key functions: ppsd_core3D and fastppsdcore3D
		 * *******/
//		boolean [][][] kMask = new boolean[zSlice][height][width];
		for (int i = 0; i < zSlice; i++)
//			for (int j = 0; j < height; j++)
//				for (int k = 0; k < width; k++)
//					kMask[i][j][k] = true;

			/*Fast version of PPSD(No zscore updating)*/
			startTime1=System.nanoTime();
			ppsd_main = new fastppsdcore3D(imArray, width,height,zSlice, p, q);
			endTime1=System.nanoTime();
			System.out.println("Running time: "+(endTime1-startTime1)/1e9); 
	}
}
