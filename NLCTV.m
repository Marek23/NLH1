function NLCTV()
clc

%%
%I = imread('C:\MAREK\MAGISTERKA\Obrazy\Obr6m3.png');
%imwrite(I,'Obr6m3.bmp') ;

%% 
f0=imread('Obr6m3.bmp');

%%
% f0=imnoise(f0,'salt & pepper',0.3);
% imwrite(uint8(f0),[ '..\ren1\q' num2str(1) '.bmp']);
% f0=f(:,:,1);

figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);
[m,n,c]=size(f0);

sw =4;
p_r=2;
s_r=5;
p_s =p_r*2+1;
p_sw=p_r*2*sw+1;
s_s=s_r*2+1;
t_r=p_r*sw+s_r;

BrokenAreaColor=255;
lamda=.01;sigma=5;h=4;

kernel=fspecial('gaussian',p_s,sigma); %¸ßË¹ºËº¯Êý
kernelk=fspecial('gaussian',p_sw,sigma);

phi=double(1-(f0(:,:,1)==BrokenAreaColor));
% phi=double(1-((f0(:,:,1)==0 | f0(:,:,1)==255) |(f0(:,:,2)==0 | f0(:,:,2)==255)|(f0(:,:,3)==0 | f0(:,:,3)==255)));
phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u0=f0;
w=zeros(s_s,s_s,m,n);

tic
for step=1:350
    step
    
    u=padarray(u0,[t_r t_r],'replicate');
    
    if step==1
        
        w=updateWeight(u0,u,h,kernel,kernelk,t_r,s_r,p_r,sw,phi,w);
        
    end
    
    if mod(step,10)==0
        
        phi=1-(u(:,:,1)==BrokenAreaColor);
        w=updateWeight2(u0,u,h,kernel,kernelk,t_r,s_r,p_r,sw,phi,PHI,w);
        
    end
    
    for i=1:m
        for j=1:n
            for k=1:c
                
                i0=i+t_r;
                j0=j+t_r;
                
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))>10
                    
                    sum_yw=sum(sum(u(i0-s_r:i0+s_r,j0-s_r:j0+s_r,k).*w(:,:,i,j,k)));
                    sum_w=sum(sum(w(:,:,i,j,k)));
                    u0(i,j,k)=(lamda*PHI(i0,j0)*f0(i,j,k)+sum_yw)/(lamda*PHI(i0,j0)+sum_w);
                    
                else
                    
                    u0(i,j,k)=f0(i,j,k);
                    
                end
                
                
            end
        end
    end
    
    
    if mod(step,10)==0
        
%         imwrite(uint8(u0),[ '..\ren1\q' num2str(step) '.bmp']);
        
        figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;
        pause(1)
        
    end
    toc
    
end
figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;

function weight=updateWeight(u0,u,h,kernel,kernelk,t_r,s_r,p_r,sw,phi,w)

weight=w;
k_r=sw*p_r;

for i=1:size(u0,1)
    for j=1:size(u0,2)
        for k=1:size(u0,3)
            
            i0=i+t_r;
            j0=j+t_r;
            
            W1 =u(i0-p_r:i0+p_r,j0-p_r:j0+p_r,k);
            W1k=u(i0-k_r:i0+k_r,j0-k_r:j0+k_r,k);
            
            if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                
                ii=1;
                for r=i0-s_r:i0+s_r
                    jj=1;
                    for s=j0-s_r:j0+s_r
                        
                        if phi(r,s)~=0
                            
                            W2 =u(r-p_r:r+p_r,s-p_r:s+p_r,k);
                            W2k=u(r-k_r:r+k_r,s-k_r:s+k_r,k);
                            d =sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel .*(W1 -W2 ).^2));
                            dk=sum(sum(phi(i0-k_r:i0+k_r,j0-k_r:j0+k_r).*kernelk.*(W1k-W2k).^2));
                            weight(ii,jj,i,j,k)=exp(-d/(h*h))*exp(-dk/(h*h));
                            
                        end
                        
                        jj=jj+1;
                    end
                    ii=ii+1;
                end
                
            end
            
        end
    end
end

function weight=updateWeight2(u0,u,h,kernel,kernelk,t_r,s_r,p_r,sw,phi,PHI,w)

weight=w;
k_r=sw*p_r;

for i=1:size(u0,1)
    for j=1:size(u0,2)
        for k=1:size(u0,3)
            
            i0=i+t_r;
            j0=j+t_r;
            
            if PHI(i0,j0)==0
                
                W1 =u(i0-p_r:i0+p_r,j0-p_r:j0+p_r,k);
                W1k=u(i0-k_r:i0+k_r,j0-k_r:j0+k_r,k);
                
                if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))~=0
                    
                    ii=1;
                    for r=i0-s_r:i0+s_r
                        jj=1;
                        for s=j0-s_r:j0+s_r
                            
                            if phi(r,s)~=0
                                
                                W2 =u(r-p_r:r+p_r,s-p_r:s+p_r,k);
                                W2k=u(r-k_r:r+k_r,s-k_r:s+k_r,k);
                                d =sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r).*kernel .*(W1 -W2).^2));
                                dk=sum(sum(phi(i0-k_r:i0+k_r,j0-k_r:j0+k_r).*kernelk.*(W1k-W2k).^2));
                                
                                weight(ii,jj,i,j,k)=exp(-d/(h*h))*exp(-dk/(h*h));
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
