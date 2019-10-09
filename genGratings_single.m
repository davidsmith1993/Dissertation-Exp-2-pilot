function [grating] = genGratings_single(window,stim,trial,data,color,where,trialWithinBlock)
%function [grating] = genGratings_single(window,stim,trial,data,color,where)
%
%generates sine-wave gratings for single display
%
%4/2/08     swe     written for relational experiment for single display (much of this based on gratingdemo.m)
%9/20/08    swe     modified to be more intuitive
%11/16/08   swe     randomize position of grating (left/right). present gray disc in empty location
%1/20/2010 swe    modified for temporal, rather than spatial, split
%3/7/2017   swe     modified for online stressor exp
%3/10/2017  swe     modified to present gratings on the fly

% *** If the grating is clipped on the sides, increase widthOfGrid.
gridEndPoints = [-1 1];
halfWidthOfGrid = abs(diff(gridEndPoints)) / 2;
widthArray = linspace(gridEndPoints(1),gridEndPoints(2),stim.widthOfGrid);% widthArray is used in creating the meshgrid.
% monitor_dist = data.monitor_dist; %monitor to edge of desk (inches)
% monitor_width = data.monitor_size; %width of screen (inches)
% npixelswide = data.res_x; %width of screen (pixels)
% npixelstall = data.res_y; %height of screen (pixels)
% pixels_per_inch = npixelswide/monitor_width;
% data.visual_angle_of_array = atan(stim.widthOfGrid/pixels_per_inch/monitor_dist)*180/pi;

% Taking the absolute value of the difference between white and gray will
% help keep the grating consistent regardless of whether the CLUT color
% code for white is less or greater than the CLUT color code for black.
absoluteDifferenceBetweenWhiteAndGray = abs(color.white - color.gray);

% ---------- Image Setup ----------
% Stores the image in a two dimensional matrix.

% Creates a two-dimensional square grid.  For each element i = i(x0, y0) of
% the grid, x = x(x0, y0) corresponds to the x-coordinate of element "i"
% and y = y(x0, y0) corresponds to the y-coordinate of element "i"
[x y] = meshgrid(widthArray, widthArray);

%transform to circular grid
circlespace = sqrt(((x).^2)+((y).^2));
%positions falling outside the circle(make =0)
[d,c] = find(circlespace>halfWidthOfGrid);
toss = [d,c];

[a,b]=find(circlespace<=halfWidthOfGrid);
keep = [a,b];

for i = 1:length(d)
    w = toss(i,1);
    e = toss(i,2);
    circlespace(w,e)=0;
end

% now the remaining values are 1;
for i = 1:length(a)
    t = keep(i,1);
    z = keep(i,2);
    circlespace(t,z)=1;
end

%generate training stim

    % *** To rotate the grating, set tiltInDegrees to a new value.
    arbitrary_or = stim.stim(trialWithinBlock,stim.orient_col);
    arbitrary_sf = stim.stim(trialWithinBlock,stim.sf_col); 

    %Helie et al. 2010 stimuli and scaling
     scaled_or = data.o_gain(1)+data.o_gain(2)*arbitrary_or; %DEGREES FROM HORIZONTAL
     scaled_sf = data.sf_gain(1)+data.sf_gain(2)*arbitrary_sf; %CYCLES/DEGREE
    %Ell et al. (2009) stimuli and scaling
%     scaled_or = arbitrary_or*180/500; %DEGREES FROM HORIZONTAL
%     scaled_sf = .25+arbitrary_sf/50 %CYCLES/DEGREE

%     %Hutchinson diss stimuli and scaling
%     scaled_or = arbitrary_or*180/800; %DEGREES FROM HORIZONTAL
%     scaled_sf = .1+arbitrary_sf/80 %CYCLES/DEGREE
 
    tiltInDegrees = scaled_or-90; % The tilt of the grating in degrees. COUNTERCLOCKWISE FROM VERTICAL
    sf=scaled_sf * data.visual_angle_of_array; % spatial freq in cycles per image

	% ---------- Image Setup ----------
    y2=sin(widthArray*pi*sf(1));
    y2=scaleif(y2, 0, 1);
    y3=ones(size(y2));
    img=(y3'*y2);
    newmax=1; newmin=-1;
    delta = (newmax-newmin);
    scaled_img = delta*img + newmin;
    
    %rotate image
    scaled_img = imrotate(scaled_img,tiltInDegrees,'bilinear','crop');

     %now have circle defined by matrix  
    scaled_img = scaled_img .* circlespace;
	
    
	% Since each entry of imageMatrix is a fraction between minus one and
	% one, multiplying imageMatrix by absoluteDifferenceBetweenWhiteAndGray
	% and adding the gray CLUT color code baseline
	% converts each entry of imageMatrix into a shade of gray:
	% if an entry of "m" is minus one, then the corresponding pixel is black;
	% if an entry of "m" is zero, then the corresponding pixel is gray;
	% if an entry of "m" is one, then the corresponding pixel is white.
	grayscaleImageMatrix1 = color.gray + absoluteDifferenceBetweenWhiteAndGray * scaled_img;
    
    disc_location = [where.xc+where.jittered_horizontal_center(trial)   where.yc+where.jittered_vertical_center(trial)];   %x and y coords of grating center
    rect = [disc_location - size(grayscaleImageMatrix1)/2     disc_location + size(grayscaleImageMatrix1)/2];

    temp = Screen('MakeTexture',window, grayscaleImageMatrix1);
    %temp = grayscaleImageMatrix1; %significant slow down if MakeTextures online
    grating.g = temp;
    grating.rect = rect;

