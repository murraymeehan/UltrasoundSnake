function mask_blur
% clear all; 
% close all; 
format compact; beep off;

% file=[    'US001.dcm'; 'US003.dcm'; 'US005.dcm'; 'US007.dcm'; 'US009.dcm'; 'US011.dcm';    
%     'US002.dcm';     'US004.dcm'; 'US006.dcm'; 'US008.dcm'; 'US010.dcm'];
% maskfile='US1_overlay_mask.png';

% x=dicomread(file(1,:));
% subplot(4,4,1); imshow(x); title(int2str(1));

mask=imread('US1_overlay_mask.png');
figure(1); 
imshow(mask)

[overlay_x overlay_y] = find(mask==0);

I=dicomread(file(5,:));
O=I;
% for k=1:(length(overlay_x)*length(overlay_y))
%     O(overlay_points(k)) = 

% I=mask;
% % High pass filter
% HighKernel = [ -1 -1 -1; -1 4 -1; -1 -1 -1 ];
% Conv_high = conv2(HighKernel, I);
% Scaled_high = imagesc(Conv_high);
% figure, imshow(Scaled_high);


end
