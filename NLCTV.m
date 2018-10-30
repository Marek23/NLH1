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

f01 = f0(:,:,1);
f02 = f0(:,:,2);
f03 = f0(:,:,3);

[m,n]=size(f01);
  
p_r=2;
s_r=5;
p_s=p_r*2+1;
s_s=s_r*2+1;
t_r=p_r+s_r;

BrokenAreaColor=255;
lamda=.01;sigma=5;h=2;

kernel=fspecial('gaussian',p_s,sigma); %¸ßË¹ºËº¯Êý

phi=double(1-((f0(:,:,1)==BrokenAreaColor) & ...
            (f0(:,:,2)==BrokenAreaColor) & ...
            (f0(:,:,3)==BrokenAreaColor)));

phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u01=f01;
u02=f02;
u03=f03;

w1=zeros(s_s,s_s,m,n);
w2=zeros(s_s,s_s,m,n);
w3=zeros(s_s,s_s,m,n);

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
        
%         imwrite(uint8(u0),[ '..\ren1\q' num2str(step) '.bmp']);
        uwyn(:,:,1) = u01;
        uwyn(:,:,2) = u02;
        uwyn(:,:,3) = u03;
        
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
        for k=1:size(u0,3)
            
            i0=i+t_r;
            j0=j+t_r;
            
            W1=u(i0-p_r:i0+p_r,j0-p_r:j0+p_r,k);
            
            if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                
                ii=1;
                for r=i0-s_r:i0+s_r
                    jj=1;
                    for s=j0-s_r:j0+s_r
                        
                        if phi(r,s)~=0
                            
                            W2=u(r-p_r:r+p_r,s-p_r:s+p_r,k);
                            d=sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel.*(W1-W2).^2));
                            weight(ii,jj,i,j,k)=exp(-d/(h*h));
                            
                        end
                        
                        jj=jj+1;
                    end
                    ii=ii+1;
                end
                
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
                
                W1=u(i0-p_r:i0+p_r,j0-p_r:j0+p_r,k);
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    
                    ii=1;
                    for r=i0-s_r:i0+s_r
                        jj=1;
                        for s=j0-s_r:j0+s_r
                            
                            if phi(r,s)~=0
                                
                                W2=u(r-p_r:r+p_r,s-p_r:s+p_r,k);
                                d=sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel.*(W1-W2).^2));
                                weight(ii,jj,i,j,k)=exp(-d/(h*h));
                                
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
