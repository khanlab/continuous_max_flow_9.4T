function  [seg,UU] = fcm_image_seg(im,mask,c_num)
%image segmentation using fuzzy c-means clustering fcm
%imput: im, 2D image,intensity between [0,1] or [0,255] or others. note: double class needed
%       mask, 0 background, 1 forground
%       c_num, cluster number
%output: seg, segmentation results, 0, background, 1-c_num different clusters corresponding to the intensitis. cluter c_num is corresponds to the brightest
%        UU, membership function matrix for each cluster
%  fcm.m(distfcm.m,initfcm.m,stepfcm.m) copyed from matlab's fuzzy toolbox

%author:Yingli Lu yinglilu@gmail.com
%last modified date: 2017/07/11

%demo codes
% im = double(mat2gray(imread('a.png')));
% mask = im>0.1;
% [seg,UU]=fcm_image_seg(im,mask,3);

mask_index = find(mask(:)==1);
im_1d=im(mask_index); %reshape to 1d according to the mask, note, this should be n x 1
% options = [2;	% exponent for the partition matrix U
% 		100;	% max. number of iteration
% 		1e-5;	% min. amount of improvement
% 		0];	% info display during iteration
% [center,U,obj_fcn] = fcm(im_1d, c_num,options); %note: center is c_num x 1, U is c_num x n
[center,U,obj_fcn] = fcm(im_1d, c_num); %note: center is c_num x 1, U is c_num x n
maxU = max(U);

seg=zeros(size(im));
U_i=zeros(size(im));
UU=zeros([size(im),c_num]);

%sort center, then to setup the segmentation results as increasing order:
[s_center, s_index]=sort(center);

for i=1:c_num
    seg(mask_index(U(s_index(i),:)==maxU))=i;
    U_i(mask_index)=U(s_index(i),:);
    UU(:,:,:,i)=U_i;
end




