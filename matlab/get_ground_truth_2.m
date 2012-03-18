% function get_ground_truth() %% A snake (active contour model) segmention algorithm for TRUS images, 
	close all; clear all; format compact;
	disp('-------- START ----------');
	list = file_list();
	N=length(list);
	N=3;
	figure; hold on;
	for i=1:N

		filename = list{i};
		clear vertex; 
		I=dicomread(filename);
		[x,y,BW{i},xi,yi]=roipoly(I,'Closed','false'); 

		xi=xi(1:end-1); % remove the redundant last point which ROIPOLY returns
		yi=yi(1:end-1);
		vertex(:,2,1) = xi;
		vertex(:,1,1) = yi;
		vertex = round(vertex);
		disp(['image: ' filename])
		disp(vertex');
		contour{i} = vertex;
		centroid{i} = [mean(xi) mean(yi)];
		plot(contour{i}(:,1),contour{i}(:,2));
		
		
	end
% 	close all;

	
% end

% function compare_truths(BW)
	prime = primes(20);
	figure; 
	IMG = uint8(zeros(size(BW{1})));
	for i=[1 2 3]
		IMG = IMG + prime(i)*uint8(BW{i});
	end
	
	IMG2 = (255/max(max(max(IMG))))*IMG;
	imshow(IMG2);
	disp('-------- END ----------');
% end


% print -dpng get_active_contour_output.png

% TODO: 
% Enter a bunch of ground truth contours manually
% Save the results for comparison in our plots
% compare them to measured contours.
