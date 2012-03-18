% function get_active_contour(input) 
	%% A snake (active contour model) segmention algorithm for TRUS images, 
	%% NOTE: to look at the effects of each function by itself, set w_int=1,
	%% w_(others)=0, or set w_ext=1,%% w_(others)=0.
tic
close all; clear all;
%% select parameters
MaxIterations = 3; 
w_int = 1;
w_ext = 0.5;  
w_tot = 1;   % values of w_total = (0.9-1.2) seem to work okay. larger is bad. 
plot_of_image = figure;% ('Visible','Off'); 

%% read the image and apply a preprocessing filter
% I = dicomread('US005.dcm');
% I=imread('US1_mask_removed.png');
I=imread('wavet_decomposed_US1_level_3.png');
[xMax yMax] = size(I); 
imshow(I); % title('This snake is considering only internal energy terms');
hold on;
screen_size = get(0, 'ScreenSize');			
set(figure(1), 'Position', [0 0 screen_size(3) (screen_size(4)*.85)] );

% I_pp = imfilter(I,Preprocessing_Filter); 

%% define the initial contour to work from
    %% Using ROIPOLY to get a set of starting points from the user
%         [x,y,mask,xi,yi]=roipoly(I);  
    %     [x,y,mask,xi,yi]=roipoly(I); 
    %     remove the redundant last point which ROIPOLY returns
%         xi=xi(1:end-1);     
%         yi=yi(1:end-1);
%         vertex(:,2,1) = xi;
%         vertex(:,1,1) = yi;
%         vertex = round(vertex) % possibly optional, but made debugging easier
    %% a set of initial vertices from ROIPOLY: (for image US005.dcm)
%     vertex(:,:,1) = [    93   118;   120   111;   149   106;   181   107;   219   102;   261   102;   298   102;   325   103;   357   101;   393   106;   431   102;   454   110;   472   145;   475   174;   476   205;   477   256;   477   309;   478   351;   478   383;   476   430;   475   462;   463   493;   435   500;   377   505;   334   513;   280   517;   239   514;   179   514;   134   512;    70   501;    61   472;    58   414;    56   336;    63   277;    67   182;    70   141 ];
	%% a set of initial vertices from ROIPOLY: (for image wavet_decomposed_US1_level_3.png)
	vertex = [
		21    40;		31    22;		36    19;		46    19;		55    25;		57    32;
		57    48;		52    54;		43    57;		33    56;		28    55;		24    44];

%% gradient of the image, used for E_ext = external force on vertices
num_vertices = length(vertex(:,1,1));
[Gxy(1,:,:),Gxy(2,:,:)] = gradient(double(I));
range_of_Gxy = [min(min(min(Gxy))) max(max(max(Gxy)))];

%% draw initial snake boundary
plot([vertex(1:end,2,1); vertex(1,2,1) ], [vertex(1:end,1,1); vertex(1,1,1)] ,'r*-'); 

for t=1:(MaxIterations)
    	%% Internal Energy (currently just elastic energy) (old version)
        E_int(1,:,t) = (vertex(end,:,t) + vertex(2,:,t))/2 - vertex(1,:,t);
    for n = 2:num_vertices-1 
        E_int(n,:,t) = (vertex(n-1,:,t) + vertex(n+1,:,t))/2 - vertex(n,:,t);
    end
        E_int(num_vertices,:,t) = (vertex(num_vertices-1,:,t) + vertex(1,:,t))/2 - vertex(num_vertices,:,t);
    %% External Energy (currently just gradient-based force)
    for n=1:num_vertices
        E_ext(n,:,t) = Gxy(:,vertex(n,1,t) , vertex(n,2,t))';
    end
    %% Defining next iteration of the snake 
    for n=1:num_vertices
        E_tot(n,:,t) = w_int*E_int(n,:,t) + w_ext*E_ext(n,:,t);
        vertex(n,:,t+1) = vertex(n,:,t) + round(w_tot*E_tot(n,:,t)); 
    end
    %% plot the next contour
	if (mod(t,1)==0)
		plot([vertex(1:end,2,t+1); vertex(1,2,t+1)] ,[vertex(1:end,1,t+1); vertex(1,1,t+1) ], 'Color', [ (t/MaxIterations) (1-(t/MaxIterations)) 1] );
	end
end 

toc
%% plot min and max energy levels
minE_int = permute( min(min(E_int(:,:,:))) , [3 1 2]);
maxE_int = permute( max(max(E_int(:,:,:))) , [3 1 2]);
minE_ext = permute( min(min(E_ext(:,:,:))) , [3 1 2]);
maxE_ext = permute( max(max(E_ext(:,:,:))) , [3 1 2]);

print -dpng get_active_contour_output.png
% end