# SynQuant, with some extra functions

This repository provides a fork of the synapse quantification tool SynQuant that provides a simplified command, "SynQuantSimple".
It only works with single-channel images.
However, this makes it easier to integrate into other software tools in the Fiji ecosystem, like [SynBot](https://github.com/Eroglu-Lab/Syn_Bot).
For more details, please visit the [main repo of SynQuant](https://github.com/yu-lab-vt/SynQuant).

## SynQuant for single channel
This is the simplified version that processes a single channel. We support 3D images with multiple z-slices.
SynQuant will automatically read the current image as input.

## SynQuant for batch processing
This version allows using a configuration file as input to run SynQuant in batch.

For example, in ImageJ scripts, for each current image, use
```
run("SynQuantBatch", "C:\\param.txt");
```

The configuration file looks like
```
zscore_thres=2
MinSize=10
MaxSize=200  
minFill=0.5  
maxWHRatio=4
zAxisMultiplier=1
```

The meaning of these parameters can be found in the documentation of SynQuant.

If no configuration file is specified, a dialogue will be prompt for you to specify the file.

## Other choices for batch processing
For processing large amounts of data, another choice is to call SynQuant Java classes directly from MATLAB. An example is given [here](https://github.com/freemanwyz/SynQuant_MATLAB_Java). Note that this only contains a subset of the features of the Fiji plug-in, and does not provide a GUI. For a smaller amount of images, it is better to use the Fiji plug-in.

You may also try to call SynQuant using the Python-ImageJ interface [PyImageJ](https://github.com/imagej/pyimagej), but we have not tested that yet.

## Reference
[1] Yizhi Wang*, Congchao Wang*, Petter Ranefall, Gerard Joey Broussard, Yinxue Wang, Guilai Shi, Boyu Lyu, Chiung-Ting Wu, Yue Wang, Lin Tian, Guoqiang Yu. (2020). SynQuant: An Automatic Tool to Quantify Synapses from Microscopy Images, Bioinformatics, 36(5), 1599–1606
