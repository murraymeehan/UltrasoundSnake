function wav3 = wavelet_decompose(source)
	clear all; close all; format compact; beep off;
	% Read in images
% 	source=imread('US1_mask_removed.png');
	[row1,col1] = size(source);
	N=3;						% Levels of discrete 2-D wavelet transforms to perform
	SOURCE = cell(4,N);			% each column will hold the 4 matrices from each level of 2D DWT

	for i=1:N
			[SOURCE{1,i},SOURCE{2,i},SOURCE{3,i},SOURCE{4,i}] = dwt2(source,'haar');
			source=SOURCE{1,i};
	end

	for i=1:N
		for j=1:4
			SOURCE_plot{j,i}=(1/512)*(2^(1-i))*abs(SOURCE{j,i}); 		%scale results of DWT for plotting
		end
% 		imshow(SOURCE{1,i})
% 		print -dpng wavet_decomposed_US1_level_3.png;
	end

	imshow(SOURCE_plot{1,3}),% daspect([1 1 1]);
% 	print -dpng wavet_decomposed_US1_level_3.png;
	imwrite(SOURCE_plot{1,3},'wavet_decomposed_US1_level_3.png');
	wav3 = SOURCE_plot{1,3};
% 	screen_size = get(0, 'ScreenSize');			
% 	set(figure(1), 'Position', [0 0 screen_size(3) (screen_size(4)*.85)] );
end