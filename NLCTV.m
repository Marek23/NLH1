function NLCTV()
clc; clear;
mex main.c;

f0=imread('Banan.bmp');

figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);
[m,n,c]=size(f0);

sw =1;%%mój parametr(gdy =1 ==> oryginalny algorytm)

p_r=2;%%promieñ patch'a  jak poprzednio
s_r=4;%%promieñ window'a jak poprzednio


p_s =p_r*2+1;
p_sw=p_r*2*sw+1;
s_s=s_r*2+1;
t_r=p_r*sw+s_r;

BrokenAreaColor=255;

lamda=.01;sigma=5;h=7;%%parametry jak poprzednio

kernel= fspecial('gaussian',p_s,sigma);
kernelk=fspecial('gaussian',p_sw,sigma);

phi=double(1-((f0(:,:,1)==0) & ...
              (f0(:,:,2)==BrokenAreaColor) & ...
              (f0(:,:,3)==0)));

phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u0=f0;
   
u0r = main(m,n,c,u0(:),...
    m+2*t_r,n+2*t_r,h,p_s,kernel(:),...
    p_sw,kernelk(:),t_r,s_r,p_r,sw,phi(:),...
    PHI(:),s_s,lamda,f0(:));

u0 = reshape(u0r,[m,n,c]);
figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;

end