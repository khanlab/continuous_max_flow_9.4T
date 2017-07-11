%% 9.4T data segmentation: Use continuous Max-Flow method.

% Script to do tissue segmentation on 9.4T images with continous Max-Flow

% input: masked image.
% output: segmented image.(bg:0,wm:1,gm:2,csf:3)

% Author: YingLi Lu, yinglilu@gmail.com
% Date: 2017-07-06
% Note: need install, and add to path
%    1. matlab_toolbox_matlab_nifti
%    2. https://github.com/ASETS/asetsMatlabMaxFlow
%    3. support .nii file only.

% Acknowledgement: John Baxter
   
clear all; close all;

input_filename = 'D:\shared_folder\9.4T_projects\reports\case_1_which_method_and_tune_parameters\asetsHMF3D_MRI_PreExReg_masked.nii';
output_filename = 'D:\shared_folder\9.4T_projects\reports\case_1_which_method_and_tune_parameters\asetsHMF3D_MRI_PreExReg_masked_fcm_seg_alpha_1.5.nii';

% 1. Load a volume image 
nii =  MRIread(input_filename);
img = nii.vol;

[r,c,s] = size(img);

% 2. Normalize the image intensity to [0,1]:
img = single(img);
img_n = (img - min(img(:))) / (max(img(:)) - min(img(:)));

mask = img_n>0;

% 3. build models 
%
% data term: intensity distance: too rough.
%
% method 1. intensity based.
% C_WM = abs(img_n-0.58);
% C_GM = abs(img_n-0.48);
% C_CSF = abs(img_n-0.20);
% C_K = abs(img_n-0.12);
 
% method 2. gaussian mixture model, need segment the object first. too
% slow and sensitive to intensity changes cross datasets(since based on
% histogram).
% [seg,PPP] =gmm_image_seg_3d(uint8(img_n*255),mask,3);
% %cost function.
% C_WM = -log(PPP(:,:,:,3));
% C_GM = -log(PPP(:,:,:,2));
% C_CSF = -log(PPP(:,:,:,1));
% C_K = mask;
 
% method 3. fuzzy c means
[seg,PPP]=fcm_image_seg(img_n,mask,3);
%cost function.
C_WM = -log(PPP(:,:,:,3));
C_GM = -log(PPP(:,:,:,2));
C_CSF = -log(PPP(:,:,:,1));
C_K = mask;

%smooth term:
[Gmag,Gazimuth,Gelevation] = imgradient3(img_n);
Gmag = Gmag + 0.01;

% alphaWM = 0.05./Gmag;
% alphaGM = 1.2./Gmag; %%with smaller alpha weight: shows scattered CSF on the brain boundary.  %2.2
% alphaCSF = 3./Gmag; %with smaller alpha weight: shows scattered CSF on the brain boundary.
% alphaBrain = 0.05./Gmag;
% alphaK = 0.05./Gmag;

%Ali want keep same weights(simple)
alpha=1.5;
alphaWM = alpha./Gmag;
alphaGM = alpha./Gmag; %%with smaller alpha weight: shows scattered CSF on the brain boundary.  %2.2
alphaCSF = alpha./Gmag; %with smaller alpha weight: shows scattered CSF on the brain boundary.
alphaBrain = alpha./Gmag;
alphaK = alpha./Gmag;

whitematter = asetsHMF3D({},alphaWM,C_WM);
graymatter = asetsHMF3D({},alphaGM,C_GM);
CSF = asetsHMF3D({},alphaCSF,C_CSF);
background = asetsHMF3D({},alphaK,C_K);

% if(gpuDeviceCount)
%     alphaWM = gpuArray(single(alphaWM));
%     alphaGM = gpuArray(single(alphaGM));
%     alphaCSF = gpuArray(single(alphaCSF));
%     alphaBrain = gpuArray(single(alphaBrain));
%     alphaK = gpuArray(single(alphaK));
%     
%     C_WM = gpuArray(single(C_WM));
%     C_GM = gpuArray(single(C_GM));
%     C_CSF = gpuArray(single(C_CSF));
%     C_K = gpuArray(single(C_K));
%   
% else
%     error('No CUDA devices detected. Skipping computation.');
% end

brain = asetsHMF3D({whitematter,graymatter,CSF},alphaBrain);
source = asetsHMF3D({brain,background},0);

%keep the second parameter to 0.1 (John said)
%third parameter: related to onvergence process(pretty much like gradient
%descent step), if the data intensity is 0-1, 0.25 works very well).
source.MaxFullFlow(100,0.1,0.25);

wm_mask = whitematter.u > 0.5;
gm_mask = graymatter.u > 0.5;
csf_mask = CSF.u>0.5;
background_mask = background.u>0.5;

all_mask = wm_mask + gm_mask *2 + csf_mask*3;% + background_mask*4;

%write nii
nii.vol=all_mask;
MRIwrite(nii,output_filename,'uchar');
