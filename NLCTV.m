function NLCTV()
clc; clear;
mex main.c;

images = dir('C:\MAREK\MAGISTERKA\Obrazy\imgmask\*.png');

for image=1:length(images)

sw =1; %%mój parametr(gdy =1 ==> oryginalny algorytm)

h=7;
while h<8

s_r = 1;
while s_r <22

p_r = 1;
while p_r <= s_r && p_r <16

clearvars -except image images h s_r p_r sw
tic

f0=imread(['C:\MAREK\MAGISTERKA\Obrazy\imgmask\' images(image).name]);

%figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);
[m,n,c]=size(f0);


p_s =p_r*2+1;
p_sw=p_r*2*sw+1;
s_s =s_r*2+1;
t_r =p_r*sw+s_r;

BrokenAreaColor=240;


lamda=.01;sigma=5; %%parametry jak poprzednio

kernel= fspecial('gaussian',p_s,sigma);
kernelk=fspecial('gaussian',p_sw,sigma);

phi=double(1-((f0(:,:,1) < 10) & ...
              (f0(:,:,2) >BrokenAreaColor) & ...
              (f0(:,:,3) < 10)));

phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u0=f0;
   
u0r = main(m,n,c,u0(:),...
    m+2*t_r,n+2*t_r,h,p_s,kernel(:),...
    p_sw,kernelk(:),t_r,s_r,p_r,sw,phi(:),...
    PHI(:),s_s,lamda,f0(:));

u0 = reshape(u0r,[m,n,c]);
%figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;

imwrite(uint8(u0),['C:\MAREK\MAGISTERKA\Obrazy\nlctvtest\' 'h' num2str(h) 'ps' num2str(p_r) 'sw' num2str(s_r) images(image).name]);

images(image).name
toc

DiaryName = strcat('C:\MAREK\MAGISTERKA\Obrazy\nlctvtest\','h',num2str(h),'ps',num2str(p_r),'sw',num2str(s_r),images(image).name,'.log');
diary(DiaryName)
diary on

clc

p_r=p_r+2;
end

s_r=s_r+2;
end

h=h+1;
end

end