# 9.4T brain tissue segmentation using continous max-flow
## Maintainer: YingLi Lu yingliu@gmail.com

* matlab_nifti:/ matlab scripts for reading/writing nifti.
* max_flow_94t.m: main matlab script.
* fcm_image_seg_3d.m: fuzzy c-mean segmentation called by max_flow_94t.m
* gmm_image_seg_3d.m: GMM segmentation called by max_flow_94t.m

### Getting started:
* download max-flow toolbox: https://github.com/ASETS/asetsMatlabMaxFlow, add it to matlab path
* add matlab_nifti to matlab path
* edit max_flow_94t.m, modify the 'input_filename' and 'output_filename', then run it.


 
