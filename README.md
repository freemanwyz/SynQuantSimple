# SynQuant, with some extra functions

This repository provides a fork of the synapse quantification tool SynQuant that provides a simplified command, "SynQuantSimple".
It only works with single-channel images.
However, this makes it easier to integrate into other software tools in the Fiji ecosystem, like [SynBot](https://github.com/Eroglu-Lab/Syn_Bot).

For more details, please visit the main [repo](https://github.com/yu-lab-vt/SynQuant) of SynQuant.

## Batch processing
For processing large amounts of data, you may wish to write some scripts in ImageJ. Another choice is to call SynQuant Java classes directly from MATLAB. An example is given [here](https://github.com/freemanwyz/SynQuant_MATLAB_Java). Note that this only contains a subset of the features of the Fiji plug-in, and does not provide a GUI. For a smaller amount of images, it is better to use the Fiji plug-in.

You may also try to call SynQuant using the Python-ImageJ interface [PyImageJ](https://github.com/imagej/pyimagej), but we have not tested that yet.

## Reference
[1] Yizhi Wang*, Congchao Wang*, Petter Ranefall, Gerard Joey Broussard, Yinxue Wang, Guilai Shi, Boyu Lyu, Chiung-Ting Wu, Yue Wang, Lin Tian, Guoqiang Yu. (2020). SynQuant: An Automatic Tool to Quantify Synapses from Microscopy Images, Bioinformatics, 36(5), 1599–1606
