function NLCTV()
clc
clear
close all

%%
%I = imread('C:\MAREK\MAGISTERKA\Obrazy\Obr6m3.png');
%imwrite(I,'Obr6m3.bmp') ;

%% 
f0=imread('Obr4m4.bmp');

%%
% f0=imnoise(f0,'salt & pepper',0.3);
% imwrite(uint8(f0),[ '..\ren1\q' num2str(1) '.bmp']);
% f0=f(:,:,1);

figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);

[m,n,ny]=size(f0);
  
p_r=2;
s_r=5;
p_s=p_r*2+1;
s_s=s_r*2+1;
t_r=p_r+s_r;

BrokenAreaColor=255;
lamda=.01;sigma=5;h=2;

kernel=fspecial('gaussian',p_s,sigma);

maska = double((f0(:,:,1)==BrokenAreaColor) & ...
            (f0(:,:,2)==BrokenAreaColor) & ...
            (f0(:,:,3)==BrokenAreaColor));

R = f0(:,:,1);
%% pobieram od u¿ytkownika salient structures
Ps = getPencils();
for i=1:length(Ps)
    SPLinked(:,:,i) = zeros(m,n);
    P = Ps{i};
    for p=1:length(P)
        SPLinked(int32(P(p,2)),int32(P(p,1)),i)=1;
    end
    se = strel('disk',10);
    SPLinked(:,:,i) = imclose(SPLinked(:,:,i),se);
    SPLinked(:,:,i) = imclose(SPLinked(:,:,i),se);
    SPLinked(:,:,i) = imclose(SPLinked(:,:,i),se);
    
    figure
    imshow(SPLinked(:,:,i));
end

%% Poczatek algorytm Criminisi
minC=0.05;
[m,n,nz] = size(f0);
mask_fill = double(~single(maska)); %% zmieniam maskê, algorytm tam gdzie ma 0 widzi maskê
%% maska musi byæ o dwa piksele szersza ni¿ maska na zniszczonym obrazie, potrzebne do ró¿niczek
mask = imerode(imerode(mask_fill,ones(3,3)),ones(3,3));
%% inicjalizacja Confidence Term
C = initializeC(mask, minC);
%% parametr rozmmiaru patch'a
patch_size = 8; %jako promien
%% patch nie mo¿e wychodziæ poza granice obrazu
ograniczenie = ones(m,n);ograniczenie(1:patch_size,:) = 0;ograniczenie(m-patch_size:m ,:) = 0;
ograniczenie(:,1:patch_size) = 0;ograniczenie(:,n-patch_size:n) = 0;
%% wzorzec mówi sk¹d mogê czerpaæ dane
wzorzec = mask_fill;
for i = 1: 1+patch_size
    wzorzec = imerode(wzorzec,ones(3,3)); %%maska poszerzona, ¿eby nie trafiaæ szukaj¹c wzorca w pole objête mask¹
end
maska = wzorzec;
for k=1:40
    wzorzec = imerode(wzorzec,ones(3,3));
end
wzorzec = (maska - wzorzec).*ograniczenie;
figure
imshow(wzorzec);
title('Punkty do czerpania Fi_q')

%% ograniczenie wyznaczam do drugiej czêœci algorytmu po SalientStructure. 
% Nie chcê korzystaæ w drugiej czêœci z wzorców g³ównych struktur
wzorzec_bez_Salient_Structure = zeros(m,n);
%% wype³nianie Salient Structures
for k = 1: length(Ps)
    %% poszerzam Salient Structures do czerpania ograniczonego wzorca
    wzorzec_linii = SPLinked(:,:,k).*wzorzec;
    for dil=1:5
        wzorzec_linii = imdilate(wzorzec_linii,ones(3,3));
    end
    wzorzec_linii = wzorzec.*wzorzec_linii;
    wzorzec_bez_Salient_Structure = wzorzec_bez_Salient_Structure + wzorzec_linii;
    %% wype³niam Salient Structure
    maska = SPLinked(:,:,k).*~mask_fill;%maska = zeros(nx,ny);
    while(any(maska(:)))
        [px,py] = find(maska > 0);
        px = px(1);
        py = py(1);
        [fixx, fixy] = ustawWspolrzednePatcha(px,py,patch_size,m,n);
        FI_p.I = f0(px+fixx-patch_size:px+fixx+patch_size,py+fixy-patch_size:py+fixy+patch_size,:);
        FI_p.m = mask_fill(px+fixx-patch_size:px+fixx+patch_size,py+fixy-patch_size:py+fixy+patch_size);
        [qx, qy] = SSD3(f0,FI_p,wzorzec_linii,patch_size);
        FI_q = f0(qx-patch_size:qx+patch_size,qy-patch_size:qy+patch_size,:);
        FR_q = R(qx-patch_size:qx+patch_size,qy-patch_size:qy+patch_size,:);
        for x= 1:2*patch_size+1
            for y= 1:2*patch_size+1
                if(mask_fill(px+fixx-1+x-patch_size,py+fixy-1+y-patch_size) == 0)
                    f0(px+fixx-1+x-patch_size,py+fixy-1+y-patch_size,:) = FI_q(x,y,:);
                    R(px-1+x-patch_size,py-1+y-patch_size,:) = FR_q(x,y,:);
                    mask_fill(px+fixx-1+x-patch_size,py+fixy-1+y-patch_size) = 1;
                end
                maska(px+fixx-1+x-patch_size,py+fixy-1+y-patch_size) = 0;
            end
        end
        figure; imagesc(uint8(f0)); colormap(gray); axis off; axis equal;
    end
    maska = mask_fill;
end

phi=maska;
phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

f01 = f0(:,:,1);
f02 = f0(:,:,2);
f03 = f0(:,:,3);

w1=zeros(s_s,s_s,m,n);
w2=zeros(s_s,s_s,m,n);
w3=zeros(s_s,s_s,m,n);

u01 = f0(:,:,1);
u02 = f0(:,:,2);
u03 = f0(:,:,3);

tic
for step=1:350
    step
    
    u1=padarray(u01,[t_r t_r],'symmetric');
    u2=padarray(u02,[t_r t_r],'symmetric');
    u3=padarray(u03,[t_r t_r],'symmetric');
    
    if step==1
        
        w1=updateWeight(u01,u1,h,kernel,t_r,s_r,p_r,phi,w1);
        w2=updateWeight(u02,u2,h,kernel,t_r,s_r,p_r,phi,w2);
        w3=updateWeight(u03,u3,h,kernel,t_r,s_r,p_r,phi,w3);
        
    end
    
    if mod(step,10)==0
        
        phi=1-(u1==BrokenAreaColor);
        w1=updateWeight2(u01,u1,h,kernel,t_r,s_r,p_r,phi,PHI,w1);
        w2=updateWeight2(u02,u2,h,kernel,t_r,s_r,p_r,phi,PHI,w2);
        w3=updateWeight2(u03,u3,h,kernel,t_r,s_r,p_r,phi,PHI,w3);
        
    end
    
    for i=1:m
        for j=1:n
                
                i0=i+t_r;
                j0=j+t_r;
                
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    
                    sum_yw1=sum(sum(u1(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w1(:,:,i,j)));
                    sum_w1=sum(sum(w1(:,:,i,j)));
                    u01(i,j)=(lamda*PHI(i0,j0)*f01(i,j)+sum_yw1)/(lamda*PHI(i0,j0)+sum_w1);
                    
                    sum_yw2=sum(sum(u2(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w2(:,:,i,j)));
                    sum_w2=sum(sum(w2(:,:,i,j)));
                    u02(i,j)=(lamda*PHI(i0,j0)*f02(i,j)+sum_yw2)/(lamda*PHI(i0,j0)+sum_w2);
                    
                    sum_yw3=sum(sum(u3(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w3(:,:,i,j)));
                    sum_w3=sum(sum(w3(:,:,i,j)));
                    u03(i,j)=(lamda*PHI(i0,j0)*f03(i,j)+sum_yw3)/(lamda*PHI(i0,j0)+sum_w3);
                    
                else
                    
                    u01(i,j)=f01(i,j);
                    u02(i,j)=f02(i,j);
                    u03(i,j)=f03(i,j);
                    
                end
                
                
        end
    end
    
    
    if mod(step,10)==0
        
        uwyn(:,:,1) = u01;
        uwyn(:,:,2) = u02;
        uwyn(:,:,3) = u03;
        
        imwrite(uint8(uwyn),[ 'testy\q' num2str(step) '.bmp']);

        
        figure; imagesc(uint8(uwyn)); colormap(gray); axis off; axis equal;
        pause(1)
        
    end
    toc
    
end

uwyn(:,:,1) = u01;
uwyn(:,:,2) = u02;
uwyn(:,:,3) = u03;

figure; imagesc(uint8(uwyn)); colormap(gray); axis off; axis equal;

function weight=updateWeight(u0,u,h,kernel,t_r,s_r,p_r,phi,w)

weight=w;

for i=1:size(u0,1)
    for j=1:size(u0,2)
            
            i0=i+t_r;
            j0=j+t_r;
            
            W1=u(i0-p_r:i0+p_r,j0-p_r:j0+p_r);
            
            if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                
                ii=1;
                for r=i0-s_r:i0+s_r
                    jj=1;
                    for s=j0-s_r:j0+s_r
                        
                        if phi(r,s)~=0
                            
                            W2=u(r-p_r:r+p_r,s-p_r:s+p_r);
                            d=sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel.*(W1-W2).^2));
                            weight(ii,jj,i,j)=exp(-d/(h*h));
                            
                        end
                        
                        jj=jj+1;
                    end
                    ii=ii+1;
                end
                
            end
    end
end

function weight=updateWeight2(u0,u,h,kernel,t_r,s_r,p_r,phi,PHI,w)

weight=w;

for i=1:size(u0,1)
    for j=1:size(u0,2)
        for k=1:size(u0,3)
            
            i0=i+t_r;
            j0=j+t_r;
            
            if PHI(i0,j0)==0
                
                W1=u(i0-p_r:i0+p_r,j0-p_r:j0+p_r);
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    
                    ii=1;
                    for r=i0-s_r:i0+s_r
                        jj=1;
                        for s=j0-s_r:j0+s_r
                            
                            if phi(r,s)~=0
                                
                                W2=u(r-p_r:r+p_r,s-p_r:s+p_r);
                                d=sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel.*(W1-W2).^2));
                                weight(ii,jj,i,j)=exp(-d/(h*h));
                                
                            end
                            
                            jj=jj+1;
                        end
                        ii=ii+1;
                    end
                    
                end
                
            end
            
        end
    end
end
