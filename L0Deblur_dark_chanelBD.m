function S = L0Deblur_dark_chanelBD(Im, kernel, lambda, wei_grad, kappa)
%%
% Image restoration with L0 regularized intensity and gradient prior
% The objective function:
% S = argmin ||I*k - B||^2 + \lambda |D(I)|_0 + wei_grad |\nabla I|_0
%% Input:
% @Im: Blurred image
% @kernel: blur kernel
% @lambda: weight for the L0 intensity prior
% @wei_grad: weight for the L0 gradient prior
% @kappa: Update ratio in the ADM
%% Output:
% @S: Latent image
%
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
%        Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 

if ~exist('kappa','var')
    kappa = 2.0;
end
%% pad image
% H = size(Im,1);    W = size(Im,2);
% Im = wrap_boundary_liu(Im, opt_fft_size([H W]+size(kernel)-1));
%%
S = Im;
betamax = 1e5;
fx = [1, -1];
fy = [1; -1];
[N,M,D] = size(Im);
sizeI2D = [N,M];
otfFx = psf2otf(fx,sizeI2D);
otfFy = psf2otf(fy,sizeI2D);
%%
KER = psf2otf(kernel,sizeI2D); %第s层的核函数
Den_KER = abs(KER).^2;
%%
Denormin2 = abs(otfFx).^2 + abs(otfFy ).^2;
if D>1
    Denormin2 = repmat(Denormin2,[1,1,D]);%复制成三层
    KER = repmat(KER,[1,1,D]);
    Den_KER = repmat(Den_KER,[1,1,D]);
end
Normin1 = conj(KER).*fft2(S);%原始图像S*模糊核K
%% pixel sub-problem
%%
dark_r = 45; %% Fixed size!
%mybeta_pixel = 2*lambda;
%[J, J_idx] = dark_channel(S, dark_r);
mybeta_pixel = lambda/(graythresh((S).^2));%graythresh计算全局阈值，使阈值化的黑白像素的类内方差最小化
maxbeta_pixel = 8;%2^3;
while mybeta_pixel< maxbeta_pixel
    %% 
    [J, J_idx] = dark_channel(S, dark_r);
    u = J;
    if D==1
        t = u.^2<lambda/mybeta_pixel;%u<thre
    else
        t = sum(u.^2,3)<lambda/mybeta_pixel;
        t = repmat(t,[1,1,D]);
    end
    u(t) = 0;%灰度通道中小于阈值的置为0
    
    clear t;
    u = assign_dark_channel_to_pixel(S, u, J_idx, dark_r);
    %返回的是根据u的值进行修改的u
    BS=1-S;
    [BJ, BJ_idx] = dark_channel(BS, dark_r);%亮通道
    bu = BJ;
    if D==1
        t = bu.^2<lambda/mybeta_pixel;
    else
        t = sum(bu.^2,3)<lambda/mybeta_pixel;
        t = repmat(t,[1,1,D]);
    end
    bu(t) = 0;
    %
    clear t;
    bu = assign_dark_channel_to_pixel(BS, bu, BJ_idx, dark_r);
    %% Gradient sub-problem
    beta = 2*wei_grad;
    %beta = 0.01;
    while beta < betamax
        Denormin   = Den_KER + beta*Denormin2 + 2*mybeta_pixel;
        %
        h = [diff(S,1,2), S(:,1,:) - S(:,end,:)];%diff(S,n,dim) :沿 dim 指定的维计算的第 n 个差分
        %前n-1列为列的差分，最后一列为第一列与最后一列的差值，相当于全差值
        
        v = [diff(S,1,1); S(1,:,:) - S(end,:,:)];
        if D==1
            t = (h.^2+v.^2)<wei_grad/beta;
        else
            t = sum((h.^2+v.^2),3)<wei_grad/beta;
            t = repmat(t,[1,1,D]);
        end
        h(t)=0; v(t)=0;
        clear t;
        %
        Normin2 = [h(:,end,:) - h(:, 1,:), -diff(h,1,2)];%二阶差分？
        Normin2 = Normin2 + [v(end,:,:) - v(1, :,:); -diff(v,1,1)];
        %
        FS = (Normin1 + beta*fft2(Normin2) + mybeta_pixel*fft2(u)+mybeta_pixel*fft2(1-bu))./Denormin;
        S = real(ifft2(FS));
        %%
        beta = beta*kappa;
        if wei_grad==0
            break;
        end
    end
    mybeta_pixel = mybeta_pixel*kappa;
end
%
end
