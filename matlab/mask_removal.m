function output = mask_removal(input, mask)
% mask=imread('US1_overlay_mask.png');
% input = dicomread('US005.dcm');

% tr=57; 
% mask = mask(trim:end-trim,trim:end-trim);
% input = input(trim:end-trim,trim:end-trim);
% mask=trim(mask,tr);
% input=trim(input,tr);
% output = trim(output,tr);

[sx sy] = size(input);
[X Y] = find(mask==255);

output = input;
for t=1:length(X)
	x=X(t); y=Y(t);
	if (x==0)||(x==sx)||(y==0)||(y==sx)
		output(x,y)=input(x,y);
	else
		output(x,y)=input(x+1,y);%+input(x-1,y); % other more complex interpolation schemes can go here, but it's optional.
		output(x,y)=output(x,y)/2;
	end
end

% figure; 
% hold on;
% subplot(2,2,1);
% imshow(input);
% subplot(2,2,2);
% imshow(mask);
% subplot(2,2,3);
% imshow(output);
imwrite(output,'US1_mask_removed.png');

% function O = trim(I,t)
% 	O = I(t:end-t,t:end-t);
	