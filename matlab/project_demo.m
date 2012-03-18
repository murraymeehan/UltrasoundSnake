function demo_to_submit_as_code_demo
	disp('--- start of demo ---');
	tic
	clear all
	close all
	format compact
	beep off
	
	dataset		=	file_list(1);				% Using dataset #1, datasets 1 to 6 available
	MASK		=	ac_fetch_mask(dataset); 
	CONTOURS	=	[];							% Initialize the contours array

	disp('Contour the median slice, which is usually clearet and easiest to identify.');
	MEDIAN_POS	= floor(length(dataset)/2);
	
	[AC_MEDIAN] = ac_contour(MEDIAN_POS,dataset,MASK,CONTOURS);
	CONTOURS{MEDIAN_POS} = AC_MEDIAN;	
	title(['image #' num2str(MEDIAN_POS) ' (first image contoured)'])
	
	disp('Contour the first half of the images in the data set.');
	AC=AC_MEDIAN; % Start with the median image's mask and contour outwards to start of dataset
	for k=MEDIAN_POS-1:-1:1
		[AC] = ac_contour(k,dataset,MASK,CONTOURS,AC);
		CONTOURS{k} = AC;	
	end

	disp('Contour the last half of the images in the data set.');
	AC=AC_MEDIAN; % Start with the median image's mask and contour outwards to end of dataset
	for k=MEDIAN_POS+1:1:length(dataset)
		[AC] = ac_contour(k,dataset,MASK,CONTOURS,AC);
		CONTOURS{k} = AC;	
	end

	disp('Constructuing 3D model from contours.');
	figure; ac_plot_3d(CONTOURS);
	figure(1); % Raise the 2D contours window.
	toc
	disp('--- end of demo ---');
end

function AC = ac_initalize(I,MODE)
%% define the initial contour to use when contouring the first slice
	if nargin < 2
		MODE = 1;
	end
	switch MODE
		case {1}
		%% Mode 1: Use contour values previously measured using ROIPOLY.
			AC=[177   231   296   398   460   460   459   427   378   310   245   199   184;
			   288   184   137   169   220   290   385   428   451   448   441   394   339]';
		case {2}
		%% Mode 2: Using ROIPOLY to get a new starting contour by
		%% asking the user to click to select 18 points, click on the
		%% first point to close the contour, then double click on the
		%% middle of the shape to finish.
			[x,y,mask,xi,yi]=roipoly(I); 
			AC(:,2) = xi(1:end-1);
			AC(:,1) = yi(1:end-1);
			AC = round(AC);
			
		%% Interpolation to 30 points (Optional)
			len=size(AC,1);
			factor = 60/len;
			AC=round(interp1([1:len],AC,[1:1/factor:len],'spline'));
	end
end

function ac_plot_3d(CONTOURS)
	%% Plot a 3D mesh object based on the contours of all images in the
	%% dataset
	NumImages = length(CONTOURS);
	NumPoints = size(CONTOURS{1},1);
	t = 0:pi/((NumImages-1)/2):2*pi;
	[X,Y,Z] = cylinder(10*sin(1+t/4),NumPoints);
	for i = 1:NumImages
		for p=1:NumPoints
			X(i,p) = CONTOURS{i}(p,1);
			Y(i,p) = CONTOURS{i}(p,2);
		end
		for p=NumPoints+1 % Handled seperately so that mesh object remains reasonable
			X(i,p) = CONTOURS{i}(1,1);
			Y(i,p) = CONTOURS{i}(1,2);
		end
		Z(i,:) = 5*i;
	end
	surf(X,Y,Z)
	axis square
	xlabel('dataset image pixel location (x)');
	ylabel('dataset image pixel location (y)');
	zlabel('depth image was taken at (mm)');
end

function O = ac_fetch_mask(dataset)
%% ac_fetch_mask accepts a list of files in the dataset, and returns a mask
%% image showing which pixels are the same in all datasets. This is based
%% on the assumption that all pixels which are the same in the entire
%% dataset of images compose the undesirable annotation grid.
	I=dicomread(dataset{1});
	O=I;
	O(:)=255;
	for i=2:length(dataset)
		I2=dicomread(dataset{i});
		O(I~=I2)=0; 
	end
	imwrite(O,'mask_US1.png');
end

function ac_plot(ac,COLOR)
	if nargin < 2
		COLOR = 0.5+0.5*rand(1,3); % constraining colour to [0.5, 1] gives only bright colors which are clearly visible.
	end
	ac = [[ac(1:end,1); ac(1,1) ],[ac(1:end,2); ac(1,2)]];
	plot( ac(:,2), ac(:,1),'Color',COLOR );
end

function O = ac_preprocessing_filter(I)
	%% The mean filter removes some speckle noise, but neccesitates a previous
	%% mask removal step. 
% 	h = fspecial('average',5);
% 	O = imfilter(I,h); 
	
	%% The median filter effectively removes mask and decreases speckle noise in one step.
	O = medfilt2(I,[10 10]);
end

function O = ac_remove_mask(I,MASK)
	%% Remove the mask pixels from an image by replacing them with the
	%% average of their neighrbouring pixels.
	[x y] = find(MASK);
	
	% bounds checking 
	[xmax ymax] = size(I);
	x(x==1)=[];
	y(y==1)=[];
	x(x==xmax)=[];
	y(y==ymax)=[];
	
	O=I;
	for i=1:length(x)
		xi = x(i); 
		yi = y(i);
		% Because of the particular mask pattern found in the dataset, all
		% mask pixels have non-mask pixels abobe, below, to the left, and
		% to the right of them.
		O(xi,yi) = mean([ O(xi,yi+1) O(xi,yi-1) O(xi-1,yi) O(xi+1,yi) ]);
	end
end

function list = file_list(set_num)
	%% Choose which of the supplied 6 datasets to operate on, defaulting to
	%% dataset #1
	if(~exist('set_num','var')) 
		set_num = 1; 
	end
	switch set_num
		case 1
			list = {'./data/US1/US001.dcm',	'./data/US1/US002.dcm',	'./data/US1/US003.dcm',	'./data/US1/US004.dcm',	'./data/US1/US005.dcm',	'./data/US1/US006.dcm',	'./data/US1/US007.dcm',	'./data/US1/US008.dcm',	'./data/US1/US009.dcm',	'./data/US1/US010.dcm',	'./data/US1/US011.dcm',	'./data/US2/US001.dcm',};
		case 2
			list = {'./data/US2/US002.dcm',	'./data/US2/US003.dcm',	'./data/US2/US004.dcm',	'./data/US2/US005.dcm',	'./data/US2/US006.dcm',	'./data/US2/US007.dcm',	'./data/US2/US008.dcm',	'./data/US2/US009.dcm',	'./data/US2/US010.dcm',	'./data/US2/US011.dcm',	'./data/US2/US012.dcm',	'./data/US2/US013.dcm',};
		case 3
			list = {'./data/US#3/US001.dcm','./data/US#3/US002.dcm','./data/US#3/US003.dcm','./data/US#3/US004.dcm','./data/US#3/US005.dcm','./data/US#3/US006.dcm','./data/US#3/US007.dcm','./data/US#3/US008.dcm','./data/US#3/US009.dcm','./data/US#3/US010.dcm','./data/US#3/US011.dcm','./data/US#3/US012.dcm',};
		case 4
			list = {'./data/US#4/US001.dcm','./data/US#4/US002.dcm','./data/US#4/US003.dcm','./data/US#4/US004.dcm','./data/US#4/US005.dcm','./data/US#4/US006.dcm','./data/US#4/US007.dcm','./data/US#4/US008.dcm','./data/US#4/US009.dcm','./data/US#4/US010.dcm','./data/US#4/US011.dcm',};
		case 5
			list = {'./data/US#5/US001.dcm','./data/US#5/US002.dcm','./data/US#5/US003.dcm','./data/US#5/US004.dcm','./data/US#5/US005.dcm','./data/US#5/US006.dcm','./data/US#5/US007.dcm','./data/US#5/US008.dcm','./data/US#5/US009.dcm','./data/US#5/US010.dcm',};
		case 6
			list = {'./data/US#6/US001.dcm','./data/US#6/US002.dcm','./data/US#6/US003.dcm','./data/US#6/US004.dcm','./data/US#6/US005.dcm','./data/US#6/US006.dcm','./data/US#6/US007.dcm','./data/US#6/US008.dcm','./data/US#6/US009.dcm','./data/US#6/US010.dcm'};
	end
end

function [AC CONTOURS] = ac_contour(POS,dataset,MASK,CONTOURS,AC) 
	%% Contouring function, which preprocesses an input image and
	%% returns the final contour resulting from MAX_ITER # of iterations of
	%% the algorithm.
	MAX_ITER	= 3;		% Note: no convergence checking is used.
	START_IMAGE = dicomread(dataset{POS}); 
	CLEAN		= ac_remove_mask(START_IMAGE,MASK); 
	I			= ac_preprocessing_filter(CLEAN); 
	if(~exist('AC','var')) 
		% Second arguement is initialization mode. 
		% Enter 1 to usepredetermined contour 
		% or enter 2 to submit your own using ROIPOLY.
		AC	= ac_initalize(CLEAN,1);  
	end

	subplot(3,4,POS); 
	imshow(I);  
	hold on;	
	title(['image #' num2str(POS)]);
	
	for i=1:MAX_ITER
		AC2 = ac_contour_algorithm(I,AC);
		AC=AC2;
		ac_plot(AC);
	end
end 

function vertex_final = ac_contour_algorithm(I,vertex)
	%% The Active Contour algorithm, one iteration only. This function is
	%% run by ac_contour. Note: There is no convergence checking used.
	
	%% weighting coefficients for the energy functions
	w_int = 1;		w_ext = 20;  	w_tot = 1/(w_int+w_ext);  w_tot=1;
	
	%% Image Gradient, used for External Energy
	num_vertices = length(vertex(:,1,1));
	[Gxy(1,:,:),Gxy(2,:,:)] = gradient(double(I),5);
	
	MaxIterations = 1; 
	for t=1:(MaxIterations)
		%% Internal Energy (Elastic energy)
		E_int(1,:,t) = (vertex(end,:,t) + vertex(2,:,t))/2 - vertex(1,:,t);
		for n = 2:num_vertices-1 
			E_int(n,:,t) = (vertex(n-1,:,t) + vertex(n+1,:,t))/2 - vertex(n,:,t);
		end
		E_int(num_vertices,:,t) = (vertex(num_vertices-1,:,t) + vertex(1,:,t))/2 - vertex(num_vertices,:,t);

		%% External Energy (Image Gadient energy)
		
		direction = circshift(vertex(:,:,t),[-1 0]) - circshift(vertex(:,:,t),[1 0]);
		normal = 0.5*[-direction(:,2) direction(:,1) ];
		for n=1:num_vertices
			xc = vertex(n,1,t)';	yc = vertex(n,2,t)';
			
			xM = abs(normal(n,1));	xm = -abs(normal(n,1));
			yM = abs(normal(n,2));	ym = -abs(normal(n,2));
			
			xM=uint16(xM);			xm=uint16(xm);
			yM=uint16(yM);			ym=uint16(ym);

			LocalImage = I((xc+xm):(xc+xM),(yc+ym):(yc+yM));
			
			xdata = [1:size(LocalImage,1)]';
			ydata = [1:size(LocalImage,2)]';
			
			xm = 1; xM = 2*abs(normal(n,1));
			ym = 1; yM = 2*abs(normal(n,2));

			innerRegion = roipoly(xdata,ydata,LocalImage,[xm xM xM]',[yM ym yM]');
			innerSum=sum(sum(LocalImage.*uint8(innerRegion)));
			
			outerRegion = roipoly(xdata,ydata,LocalImage,[xm xM xm],[yM ym ym]);
			outerSum=sum(sum(LocalImage.*uint8(outerRegion)));

			E_ext(n,:,t) = (innerSum - outerSum)*direction(n,:);
		end
		E_ext(:,:,t) = E_ext(:,:,t) / max(max(E_ext(:,:,t)));
		
		%% Define next iteration of contour
		for n=1:num_vertices
			E_tot(n,:,t) = w_int*E_int(n,:,t) + w_ext*E_ext(n,:,t);
			vertex(n,:,t+1) = vertex(n,:,t) + round(w_tot*E_tot(n,:,t)/2); 
		end
	end
% 	ac_plot(vertex(:,:,end));
	vertex_final = vertex(:,:,end);
end

