function [seg,UUU,mu,v,p]=gmm_image_seg_3d(im,mask,c_num)
%image segmentation using gaussian mixture model
%imput: im, 3D image,intensity between [0,255] note: use uint8(ime*255), if im is between [0,1]
%       mask, 0 background, 1 forground
%       c_num, cluster number
%output: seg, segmentation results, 0, background, 1-c_num different clusters corresponding to the intensitis. cluter c_num is corresponds to the brightest
%        UU, proabability matrix for each cluster
%        mu: vector of class means 
%        v: vector of class variances
%        p: vector of class proportions   

%Author:Yingli Lu yinglilu@gmail.com
%last modified date: 2017/07/11

%original author
% Author: Prof. Jose Vicente Manjon Herrera
% Email: jmanjon@fis.upv.es
% Date: 02-05-2006

%demo codes
%  im = double(mat2gray(imread('a.png')));
%  mask = im>0.1;
%  [seg,UU,mu,v,p]=gmm_image_seg(uint8(im*255),mask,3);

% check image
im=double(im);
copy=im;           % make a copy
im=im(mask(:)==1);         % vectorize im
mi=min(im);        % deal with negative 
im=im-mi+1;       % and zero values
m=max(im);
s=length(im);

h=histogram(im);
x=find(h);
h=h(x);
x=x(:);h=h(:);

% initiate parameters

mu=(1:c_num)*m/(c_num+1);
v=ones(1,c_num)*m;
p=ones(1,c_num)*1/c_num;

% start process
%figure;
sml = mean(diff(x))/1000;
while(1)
        % Expectation
        prb = distribution(mu,v,p,x);
        scal = sum(prb,2)+eps;
        loglik=sum(h.*log(scal));
        
        %Maximizarion
        for j=1:c_num
                pp=h.*prb(:,j)./scal;
                p(j) = sum(pp);
                mu(j) = sum(x.*pp)/p(j);
                vr = (x-mu(j));
                v(j)=sum(vr.*vr.*pp)/p(j)+sml;
        end
        p = p + 1e-3;
        p = p/sum(p);

        % Exit condition
        prb = distribution(mu,v,p,x);
        scal = sum(prb,2)+eps;
        nloglik=sum(h.*log(scal));                
        if((nloglik-loglik)<0.0001) break; end;        
            if 1  %show plots
               figure(100)
               clf
               %plot(x,h);
               bar(x,h)
               hold on
               plot(x,prb,'g--')
               plot(x,sum(prb,2),'r')
               drawnow
            end
end

% calculate seg
mu=mu+mi-1;   % recover real range
s=size(copy);
seg=zeros(s);

for i=1:s(1),
for j=1:s(2),
for k=1:s(3),
  c=zeros(c_num,1);
  for n=1:c_num
    c(n)=distribution(mu(n),v(n),p(n),copy(i,j,k)); 
  end
  %ying's codes 2011/03/21
  seg(i,j,k)=find(c==max(c));
  UUU(i,j,k,:)=c(:)/sum(c);
end
end
end

%sort mu, then to setup the segmentation results as increasing order:
seg = seg.*double(mask);
[s_mu,s_index]=sort(mu);
for i=1:c_num
    seg(seg==s_index(i))=i;
end

for i=1:c_num
    UUU(:,:,:,i) = UUU(:,:,:,i).*double(mask);
end




function y=distribution(m,v,g,x)
x=x(:);
m=m(:);
v=v(:);
g=g(:);
for i=1:size(m,1)
   d = x-m(i);
   amp = g(i)/sqrt(2*pi*v(i));
   y(:,i) = amp*exp(-0.5 * (d.*d)/v(i));
end

function[h]=histogram(datos)
datos=datos(:);
ind=find(isnan(datos)==1);
datos(ind)=0;
ind=find(isinf(datos)==1);
datos(ind)=0;
tam=length(datos);
m=ceil(max(datos))+1;
h=zeros(1,m);
for i=1:tam,
    f=floor(datos(i));    
    if(f>0 & f<(m-1))        
        a2=datos(i)-f;
        a1=1-a2;
        h(f)  =h(f)  + a1;      
        h(f+1)=h(f+1)+ a2;                          
    end;
end;
h=conv(h,[1,2,3,2,1]);
h=h(3:(length(h)-2));
h=h/sum(h);

