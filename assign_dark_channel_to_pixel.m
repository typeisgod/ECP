function [outImg] = assign_dark_channel_to_pixel(S, dark_channel_refine, dark_channel_index, patch_size)
%% assign dark channel value to image pixel
% The Code is created based on the method described in the following paper 
%   [1] Jinshan Pan, Deqing Sun, Hanspteter Pfister, and Ming-Hsuan Yang,
%        Blind Image Deblurring Using Dark Channel Prior, CVPR, 2016. 
% 
[M N C] = size(S);
%outImge = zeros(M, N); % 

% pad original image
padsize = floor(patch_size./2);
S_padd = padarray(S, [padsize padsize], 'replicate');%重复

% assign dark channel to pixel
m = 1:M;
        for n = 1:N%dark_channel_refine：小于阈值置为0的暗通道
            patch = S_padd(m:(m+patch_size-1), n:(n+patch_size-1),:);
            if ~isequal(min(patch(:)), dark_channel_refine(m,n))
                patch(dark_channel_index(m,n)) = dark_channel_refine(m,n);
            end
            for cc = 1:C
                S_padd(m:(m+patch_size-1), n:(n+patch_size-1),cc) = patch(:,:,cc);
            end
        end%根据dark_channel_refine回修正S_padd



outImg = S_padd(padsize + 1: end - padsize, padsize + 1: end - padsize,:);%裁剪
%% boundary processing
outImg(1:padsize,:,:) = S(1:padsize,:,:);  outImg(end-padsize+1:end,:,:) = S(end-padsize+1:end,:,:);
outImg(:,1:padsize,:) = S(:,1:padsize,:);  outImg(:,end-padsize+1:end,:) = S(:,end-padsize+1:end,:);
%边界值替换
%figure(2); imshow([S, outImg],[]);
end


