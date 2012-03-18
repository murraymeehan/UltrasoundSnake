clear all; close all; format compact; beep off;

file=[    'US001.dcm'; 'US003.dcm'; 'US005.dcm'; 'US007.dcm'; 'US009.dcm'; 'US011.dcm';    
    'US002.dcm';     'US004.dcm'; 'US006.dcm'; 'US008.dcm'; 'US010.dcm'];
maskfile='US1_overlay_mask.png';

x=dicomread(file(1,:));
% subplot(4,4,1); imshow(x); title(int2str(1));
mask=imread('US2_overlay_mask.png');

% for i=2:11
%     y=dicomread(file(i,:));
%     mask(find(x~=y))=0;
% %     subplot(4,4,i); imshow(y); title(int2str(i));
% end
% % subplot(4,4,14); imshow(mask); title('final mask');
% imwrite(mask,'US1_overlay_mask.png');
