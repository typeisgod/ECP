function [px, py, threshold]= threshold_pxpy_v1(latent,psf_size,threshold)

% 

if ~exist('threshold', 'var')%工作区中是否存在变量threshold %第三个参数有无输入，代码为没有输入
    threshold = 0;
    b_estimate_threshold = true;
else
    b_estimate_threshold = false;
end
%{
denoised = bilateral_filter(latent, 2.0, param.sigma, 'replicate', 2);
if opts.verbose
    figure(1992);imshow(denoised);title('denoised');
end
%}
denoised=latent;

%%
% derivative filters
dx = [-1 1; 0 0];
dy = [-1 0; 1 0];
%%
% px = imfilter(denoised, [0 -1 1], 'same', 'replicate');
% py = imfilter(denoised, [0;-1;1], 'same', 'replicate');
px = conv2(denoised, dx, 'valid');
py = conv2(denoised, dy, 'valid');
pm = px.^2 + py.^2;%梯度


% if this is the first prediction, then we need to find an appropriate
% threshold value by building a histogram of gradient magnitudes
if b_estimate_threshold
    pd = atan(py./px);%梯度角度
    pm_steps = 0:0.00006:2;
    H1 = cumsum(flipud(histc(pm(pd >= 0 & pd < pi/4), pm_steps)));%c=histc(a,b)统计a再b(端点)范围内的值的数量c
    %histc([0.5,1,2,3,3,4],[1,4,5])
    % 4 1 0   
    %划分为三个段 [1,4) [4,5) [5]  每个段个数分别为4 1 0
    %flipud：上下翻转
    H2 = cumsum(flipud(histc(pm(pd >= pi/4 & pd < pi/2), pm_steps)));
    H3 = cumsum(flipud(histc(pm(pd >= -pi/4 & pd < 0), pm_steps)));
    H4 = cumsum(flipud(histc(pm(pd >= -pi/2 & pd < -pi/4), pm_steps)));

    th = max([max(psf_size)*20, 10]);
    
    for t=1:numel(pm_steps)%numel 数组元素的数目 一维度等于length
        min_h = min([H1(t) H2(t) H3(t) H4(t)]);
        if min_h >= th
            threshold = pm_steps(end-t+1);
            break
        end
    end
end

% thresholding
m = pm < threshold;
while all(m(:)==1) %%当pm全部小于thresold(m全为1)
    threshold = threshold * 0.81;%阈值变小
    m = pm < threshold;
end
px(m) = 0;%pm小于threshold的x方向和y方向梯度为0
py(m) = 0;


% % update prediction parameters
% threshold = threshold * 0.9;
%% my modification
if b_estimate_threshold
    threshold = threshold;
else
    threshold = threshold./1.1;
end

end

