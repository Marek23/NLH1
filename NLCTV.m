function NLCTV()

images = dir('C:\MAREK\MAGISTERKA\Obrazy\msk\*.bmp');

% for image=1:length(images)

clc
clearvars -except images image
close all

%% 
%f0=imread(['C:\MAREK\MAGISTERKA\Obrazy\msk\' images(image).name]);
f0=imread('032g.bmp');
%%
% f0=imnoise(f0,'salt & pepper',0.3);
% imwrite(uint8(f0),[ '..\ren1\q' num2str(1) '.bmp']);
% f0=f(:,:,1);

figure; imagesc(f0); colormap(gray); axis off; axis equal;
f0=double(f0);

[m,n,c]=size(f0);
  
p_r=1;
s_r=1;
p_s=p_r*2+1;
s_s=s_r*2+1;
t_r=p_r+s_r;

BrokenAreaColor=255;
sigma=5;h=20;

kernel=fspecial('gaussian',p_s,sigma); %¸ßË¹ºËº¯Êý

phi=double(1-((f0(:,:,1)==0) & ...
            (f0(:,:,2)==BrokenAreaColor) & ...
            (f0(:,:,3)==0)));

phi=padarray(phi,[t_r t_r],'symmetric');
PHI=phi;

u0=f0;

w=zeros(s_s,s_s,m,n,c);

tic
for step=1:5000
    step
    
    u=padarray(u0,[t_r t_r],'symmetric');
    
    if step==1
        
        w=updateWeight(u0,u,h,kernel,t_r,s_r,p_r,phi,w);
        
    end
    
    if step>1
        
        phi=1-((u(:,:,1)==0) & ...
            (u(:,:,2)==BrokenAreaColor) & ...
            (u(:,:,3)==0));
        w=updateWeight2(u0,u,h,kernel,t_r,s_r,p_r,phi,PHI,w);
        
    end
    
	if(sum(phi(:)) == size(phi,1)*size(phi,2))
        break
    end
    
	for i=1:m
        for j=1:n
            
            i0=i+t_r;
            j0=j+t_r;
                
            if sum(sum(phi(i0-p_r:i0+p_r,j0-p_r:j0+p_r)))>(p_r*p_s-1) && PHI(i0,j0)<1
                    
                max_ii=1;
                max_jj=1;
                for ii=1:s_s
                    for jj=1:s_s
                        
                        x1=w(ii,jj,i,j,1);
                        x2=w(ii,jj,i,j,2);
                        x3=w(ii,jj,i,j,3);
                        
                        m1=w(max_ii,max_jj,i,j,1);
                        m2=w(max_ii,max_jj,i,j,2);
                        m3=w(max_ii,max_jj,i,j,3);
                        
                        if(sumsqr([x1,x2,x3])...
                          >sumsqr([m1,m2,m3]))
                            max_ii=ii;
                            max_jj=jj;
                        end
                    end
                end
                if(step==1)
                    u0(i,j,:)=f0(i-(s_r+1)+max_ii,j-(s_r+1)+max_jj,:);
                else
                    u0(i,j,:)=u0(i-(s_r+1)+max_ii,j-(s_r+1)+max_jj,:);
                end
            end
              
        end %end for j
	end %end for i

    
%     if mod(step,10)==0
        
%         imwrite(uint8(u0),['C:\MAREK\MAGISTERKA\Obrazy\test\' 'step' num2str(step) images(image).name]);
        figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;
        pause(1)
        
%     end
    toc
    
end

%figure; imagesc(uint8(u0)); colormap(gray); axis off; axis equal;
% end

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
