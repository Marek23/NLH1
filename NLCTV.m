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

kernel=fspecial('gaussian',p_s,sigma); %��˹�˺���

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

v1=zeros(s_s,s_s,m,n);
v2=zeros(s_s,s_s,m,n);
v3=zeros(s_s,s_s,m,n);

b1=zeros(s_s,s_s,m,n);
b2=zeros(s_s,s_s,m,n);
b3=zeros(s_s,s_s,m,n);

a1=zeros(s_s,s_s,m,n);
a2=zeros(s_s,s_s,m,n);
a3=zeros(s_s,s_s,m,n);

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
            
            nlgu1=zeros(s_s,s_s);
            nlgu2=zeros(s_s,s_s);
            nlgu3=zeros(s_s,s_s);
            for ii=1:s_s
                for jj=1:s_s
                    
                    nlgu1(ii,jj)=(u1(i0-(s_r+1)+ii,j0-(s_r+1)+jj)-u01(i,j))*sqrt(w1(ii,jj,i,j));
                    nlgu2(ii,jj)=(u2(i0-(s_r+1)+ii,j0-(s_r+1)+jj)-u02(i,j))*sqrt(w2(ii,jj,i,j));
                    nlgu3(ii,jj)=(u3(i0-(s_r+1)+ii,j0-(s_r+1)+jj)-u03(i,j))*sqrt(w3(ii,jj,i,j));
                    
                end %end for jj
            end %end for ii
                
            b1(:,:,i,j)=b1(:,:,i,j)+sqrt(w1(:,:,i,j)).*nlgu1-v1(:,:,i,j);
            b2(:,:,i,j)=b2(:,:,i,j)+sqrt(w2(:,:,i,j)).*nlgu2-v2(:,:,i,j);
            b3(:,:,i,j)=b3(:,:,i,j)+sqrt(w3(:,:,i,j)).*nlgu3-v3(:,:,i,j);
            
        end
    end
    
    for uiter=1:2
        for i=1:m
            for j=1:n
            
                i0=i+t_r;
                j0=j+t_r;
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    
                    divVB1=0;
                    divVB2=0;
                    divVB3=0;
                    for ii=1:s_s
                        for jj=1:s_s
                            lii=s_s-ii+1; hi=i-(s_r+1)+ii;
                            ljj=s_s-jj+1; hj=j-(s_r+1)+jj;
                        
                            if(hi<1)
                                hi = 1;
                            end
                            if(hi>m)
                                hi = m;
                            end
                            if(hj<1)
                                hj = 1;
                            end
                            if(hj>n)
                                hj = n;
                            end
                        
                            divVB1 = divVB1 + sqrt(w1(ii,jj,i,j))*(v1(ii,jj,i,j)-v1(lii,ljj,hi,hj));
                            divVB1 = divVB1 - sqrt(w1(ii,jj,i,j))*(b1(ii,jj,i,j)-b1(lii,ljj,hi,hj));
                            divVB2 = divVB2 + sqrt(w2(ii,jj,i,j))*(v2(ii,jj,i,j)-v2(lii,ljj,hi,hj));
                            divVB2 = divVB2 - sqrt(w2(ii,jj,i,j))*(b2(ii,jj,i,j)-b2(lii,ljj,hi,hj));
                            divVB3 = divVB3 + sqrt(w3(ii,jj,i,j))*(v3(ii,jj,i,j)-v3(lii,ljj,hi,hj));
                            divVB3 = divVB3 - sqrt(w3(ii,jj,i,j))*(b3(ii,jj,i,j)-b3(lii,ljj,hi,hj));
                           
                        end
                    end
                    
                    sum_uw1=sum(sum(u1(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w1(:,:,i,j)));
                    sum_uw2=sum(sum(u2(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w2(:,:,i,j)));
                    sum_uw3=sum(sum(u3(i0-s_r:i0+s_r,j0-s_r:j0+s_r).*w3(:,:,i,j)));
                    
                    sum_w1=sum(sum(w1(:,:,i,j)));
                    sum_w2=sum(sum(w2(:,:,i,j)));
                    sum_w3=sum(sum(w3(:,:,i,j)));
                    
                    u01(i,j)=(lamda*PHI(i0,j0)*f01(i,j)+sum_uw1-0.5*divVB1)/(lamda*PHI(i0,j0)+sum_w1);
                    u02(i,j)=(lamda*PHI(i0,j0)*f02(i,j)+sum_uw2-0.5*divVB2)/(lamda*PHI(i0,j0)+sum_w2);
                    u03(i,j)=(lamda*PHI(i0,j0)*f03(i,j)+sum_uw3-0.5*divVB3)/(lamda*PHI(i0,j0)+sum_w3);
                    
                else
                    
                    u01(i,j)=f01(i,j);
                    u02(i,j)=f02(i,j);
                    u03(i,j)=f03(i,j);
                    
                end
                
                
            end %end for j
        end %end for i
    end %end for uiter
    
    beta1=0;
    beta2=0;
    beta3=0;
    for i=1:m
        for j=1:n
            
            beta1=beta1+sqrt(sumsqr(v1(:,:,i,j)));
            beta2=beta2+sqrt(sumsqr(v2(:,:,i,j)));
            beta3=beta3+sqrt(sumsqr(v3(:,:,i,j)));
            
        end %end for j
    end %end for i
    
    mbeta=sqrt(beta1^2+beta2^2+beta3^2)+eps;
    
    for i=1:m
        for j=1:n
            
            i0=i+t_r;
            j0=j+t_r;
            
            for ii=1:s_s
                for jj=1:s_s
                    
                    a1(ii,jj,i,j)=(u1(i0-(s_r+1)+ii,j0-(s_r+1)+jj)...
                        -u01(i,j))*sqrt(w1(ii,jj,i,j))+b1(ii,jj,i,j);
                    a2(ii,jj,i,j)=(u2(i0-(s_r+1)+ii,j0-(s_r+1)+jj)...
                        -u02(i,j))*sqrt(w2(ii,jj,i,j))+b2(ii,jj,i,j);
                    a3(ii,jj,i,j)=(u3(i0-(s_r+1)+ii,j0-(s_r+1)+jj)...
                        -u03(i,j))*sqrt(w3(ii,jj,i,j))+b3(ii,jj,i,j);
                    
                end %end for jj
            end %end for ii
            
        end %end for j
    end %end for i
    
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
