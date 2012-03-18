function get_ground_truth() %% A snake (active contour model) segmention algorithm for TRUS images, 
close all; clear all; format compact;

list = file_list();
N=length(list);

for i=2:N
	clear vertex; 
	I=dicomread(list{i});
	[x,y,mask,xi,yi]=roipoly(I); 
	
	xi=xi(1:end-1); % remove the redundant last point which ROIPOLY returns
	yi=yi(1:end-1);
	vertex(:,2,1) = xi;
	vertex(:,1,1) = yi;
	vertex = round(vertex);
	disp(['image: ' list{i}])
	disp(vertex');
	
end


% print -dpng get_active_contour_output.png
