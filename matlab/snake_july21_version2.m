%% A snake (active contour model) segmention algorithm for TRUS images,
%% version #2 which uses the active contour model described in the recent
%% elec 436 reading.
%% by Murray Meehan and Irina Morozov for ELEC 435 201105 @ UVic

%% Begin the program
tic
beep off; 
close all; 
clear all;  
disp('starting program')

%% select parameters

%% read the image and apply a preprocessing filter
I = dicomread('US005.dcm');
[xMax yMax] = size(I); 
% imshow(I); 

%% define the initial contour to work from
    %% Using ROIPOLY to get a set of starting points from the user
    %     [x,y,mask,xi,yi]=roipoly(I);  
    %     [x,y,mask,xi,yi]=roipoly(I); 
    %     remove the redundant last point which ROIPOLY returns
    %     xi=xi(1:end-1);     
    %     yi=yi(1:end-1);
    %     vertex(:,2,1) = xi;
    %     vertex(:,1,1) = yi;
    %     vertex = round(vertex); % possibly optional, but made debugging easier
    %% a set of initial vertices from ROIPOLY: (for image US005.dcm)
    vertex(:,:,1) = [    93   118;   120   111;   149   106;   181   107;   219   102;   261   102;   298   102;   325   103;   357   101;   393   106;   431   102;   454   110;   472   145;   475   174;   476   205;   477   256;   477   309;   478   351;   478   383;   476   430;   475   462;   463   493;   435   500;   377   505;   334   513;   280   517;   239   514;   179   514;   134   512;    70   501;    61   472;    58   414;    56   336;    63   277;    67   182;    70   141 ];
	contour{1} = int16(vertex);
	figure; hold on;
	N=5; % number of contours to test
	r = [0.2 1.0]; % percentile range to interpolate contours on. 0% = point in center, 100% = original curve, 200% = twice as large as original, and may go over edges of image.
	for n=2:N
		X=contour{n-1}(:,1);
		Y=contour{n-1}(:,2);
		X2=X+2*n/N*(X-mean(X));
		Y2=Y+2*n/N*(Y-mean(Y));
		contour{n} = [X2 Y2];
		plot(X2,Y2);
	end
		
		
% 	figure; hold on; plot(X,Y,X2,Y2)

% % % %% gradient of the image, used for E_ext = external force on vertices
% % % num_vertices = length(vertex(:,1,1));
% % % [Gxy(1,:,:),Gxy(2,:,:)] = gradient(double(I));
% % % range_of_Gxy = [min(min(min(Gxy))) max(max(max(Gxy)))];
% % % 
% % % %% draw initial snake boundary
% % % % figure(plot_of_image,'Visible','Off'); hold on;
% % % figure(plot_of_image); hold on;
% % % plot([vertex(1:end,2,1); vertex(1,2,1) ], [vertex(1:end,1,1); vertex(1,1,1)] ,'r*-'); 
% % % 
% % % for t=1:(MaxIterations)
% % % 	
% % %     	%% Internal Energy (currently just elastic energy) (old version)
% % %         E_int(1,:,t) = (vertex(end,:,t) + vertex(2,:,t))/2 - vertex(1,:,t);
% % %     for n = 2:num_vertices-1 
% % %         E_int(n,:,t) = (vertex(n-1,:,t) + vertex(n+1,:,t))/2 - vertex(n,:,t);
% % %     end
% % %         E_int(num_vertices,:,t) = (vertex(num_vertices-1,:,t) + vertex(1,:,t))/2 - vertex(num_vertices,:,t);
% % % 
% % %     %% External Energy (currently just gradient-based force)
% % %     for n=1:num_vertices
% % %         E_ext(n,:,t) = Gxy(:,vertex(n,1,t) , vertex(n,2,t))';
% % %     end
% % %     
% % %     %% Defining next iteration of the snake 
% % %     for n=1:num_vertices
% % %         E_tot(n,:,t) = w_int*E_int(n,:,t) + w_ext*E_ext(n,:,t);
% % %         vertex(n,:,t+1) = vertex(n,:,t) + round(w_tot*E_tot(n,:,t)); 
% % %     end
% % % 
% % %     %% plot the next contour
% % %     figure(plot_of_image); hold on; 
% % % 	plot([vertex(1:end,2,t+1); vertex(1,2,t+1)] ,[vertex(1:end,1,t+1); vertex(1,1,t+1) ], 'Color', [ (t/MaxIterations) (1-(t/MaxIterations)) 1] );
% % %     
% % % end
% % % 
% % % %% plot min and max energy levels
% % % minE_int = permute( min(min(E_int(:,:,:))) , [3 1 2]);
% % % maxE_int = permute( max(max(E_int(:,:,:))) , [3 1 2]);
% % % minE_ext = permute( min(min(E_ext(:,:,:))) , [3 1 2]);
% % % maxE_ext = permute( max(max(E_ext(:,:,:))) , [3 1 2]);

% figure(plots_misc); 
% figure
% subplot(2,2,3);  hold on;   plot(minE_int,'b');     plot(maxE_int,'r');     title('E_int range');       legend('min','max');
% subplot(2,2,4); hold on;    plot(minE_ext,'b');     plot(maxE_ext,'r');     title('E_ext range');       legend('min','max');

% figure(plot_of_image,'Visible','Off');
toc
disp('end of program')
% print -dpng get_active_contour_output.png