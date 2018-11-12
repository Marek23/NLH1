function NLCTV()
clc
clear
close all

%%
%I = imread('C:\MAREK\MAGISTERKA\Obrazy\Obr6m3.png');
%imwrite(I,'Obr6m3.bmp') ;

%% 
f0=imread('031.bmp');

%%
% f0=imnoise(f0,'salt & pepper',0.3);
% imwrite(uint8(f0),[ '..\ren1\q' num2str(1) '.bmp']);
% f0=f(:,:,1);

figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);

[m,n,c]=size(f0);
  
p_r=2;
s_r=5;
p_s=p_r*2+1;
s_s=s_r*2+1;
t_r=p_r+s_r;

BrokenAreaColor=255;
lamda=.01;sigma=5;h=2;

kernel=fspecial('gaussian',p_s,sigma); %¸ßË¹ºËº¯Êý

phi=double(1-((f0(:,:,1)==0) & ...
            (f0(:,:,2)==BrokenAreaColor) & ...
            (f0(:,:,3)==0)));

phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u0=f0;

w=zeros(s_s,s_s,m,n,c);

v=zeros(s_s,s_s,m,n,c);

b=zeros(s_s,s_s,m,n,c);

a=zeros(s_s,s_s,m,n,c);

tic
for step=1:350
    step
    
    u=padarray(u0,[t_r t_r],'symmetric');
    
    if step==1
        
        w=updateWeight(u0,u,h,kernel,t_r,s_r,p_r,phi,w);
        
    end
    
    if mod(step,10)==0
        
        phi=1-((u(:,:,1)==0) & ...
            (u(:,:,2)==BrokenAreaColor) & ...
            (u(:,:,3)==0));
        w=updateWeight2(u0,u,h,kernel,t_r,s_r,p_r,phi,PHI,w);
        
    end
    
    for i=1:m
        for j=1:n
            for k=1:c
            
                i0=i+t_r;
                j0=j+t_r;
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    nlgu=zeros(s_s,s_s);
                    for ii=1:s_s
                        for jj=1:s_s
                    
                            nlgu(ii,jj)=(u(i0-(s_r+1)+ii,j0-(s_r+1)+jj,k)-u0(i,j,k))*sqrt(w(ii,jj,i,j,k));
                    
                        end %end for jj
                    end %end for ii
            
                    b(:,:,i,j,k)=b(:,:,i,j,k)+sqrt(w(:,:,i,j,k)).*nlgu-v(:,:,i,j,k);
                end
            end
        end
    end
    
    for uiter=1:2
        for i=1:m
            for j=1:n
                for k=1:c
            
                    i0=i+t_r;
                    j0=j+t_r;
                
                    if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0 && PHI(i0,j0)<1
                    
                        divVB=0;
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
                        
                                divVB = divVB + sqrt(w(ii,jj,i,j,k))*(v(ii,jj,i,j,k)-v(lii,ljj,hi,hj,k));
                                divVB = divVB - sqrt(w(ii,jj,i,j,k))*(b(ii,jj,i,j,k)-b(lii,ljj,hi,hj,k));
                                
                            end
                        end
                    
                        sum_uw=sum(sum(u(i0-s_r:i0+s_r,j0-s_r:j0+s_r,k).*w(:,:,i,j,k)));
                    
                        sum_w=sum(sum(w(:,:,i,j,k)));
                    
                        u0(i,j,k)=(lamda*PHI(i0,j0)*f0(i,j,k)+sum_uw-0.5*divVB)/(lamda*PHI(i0,j0)+sum_w);
                    
                    else
                    
                        u0(i,j,k)=f0(i,j,k);
                    
                    end
                
                end %end for k
            end %end for j
        end %end for i
    end %end for uiter
    
    beta=zeros(c);
    for i=1:m
        for j=1:n
            for k=1:c
                i0=i+t_r;
                j0=j+t_r;
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    beta(k)=beta(k)+sqrt(sumsqr(v(:,:,i,j,k)));
                end
            end %end for k
        end %end for j
    end %end for i
    
    mbeta=0;
    for k=1:c
        mbeta=mbeta+beta(k)^2;
    end
    mbeta = sqrt(mbeta)+eps;
    
    for i=1:m
        for j=1:n
            for k=1:c
            
                i0=i+t_r;
                j0=j+t_r;
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    for ii=1:s_s
                        for jj=1:s_s
                    
                            a(ii,jj,i,j,k)=(u(i0-(s_r+1)+ii,j0-(s_r+1)+jj,k)...
                                -u0(i,j,k))*sqrt(w(ii,jj,i,j,k))+b(ii,jj,i,j,k);
                    
                        end %end for jj
                    end %end for ii
                end
            end %end for k
        end %end for j
    end %end for i
    
    for i=1:m
        for j=1:n
            for k=1:c
                i0=i+t_r;
                j0=j+t_r;
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    for ii=1:s_s
                        for jj=1:s_s
                    
                            moda=sqrt(sumsqr(a(:,:,i,j,k)))+eps;
                    
                            v(ii,jj,i,j,k)=max(moda-...
                                (beta(k)/mbeta),0)*a(ii,jj,i,j,k)/moda;
                    
                        end %end for jj
                    end %end for ii
                end
            end %end for k
        end %end for j
    end %end for i
    
    if mod(step,10)==0
        
%         imwrite(uint8(u0),[ '..\ren1\q' num2str(step) '.bmp']);
        figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;
        pause(1)
        
    end
    toc
    
end

figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;

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
