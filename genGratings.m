function [data,grating] = genGratings(window,stim,data,color,where,structname,grating)
%function [data,grating] = genGratings(window,stim,data,color,where,stimmat,structname,grating)
%
%generates sine-wave gratings for single display
%e.g., grating50.structname.g = image
%      grating50.structname.rect = grating location
%           for stimulus with index = 50 (last column in stim files)
%the resulting structure serves as an image lookup table, indexed by the stimulus index, to avoid stimulus generation on the fly
%
%4/2/08     swe     written for relational experiment for single display (much of this based on gratingdemo.m)
%9/20/08    swe     modified to be more intuitive
%11/16/08   swe     randomize position of grating (left/right). present gray disc in empty location
%1/20/2010 swe    modified for temporal, rather than spatial, split
%1/16/2012  swe     modified for inference experiment

stimmat = stim.stim;

% *** If the grating is clipped on the sides, increase widthOfGrid.
gridEndPoints = [-1 1];
halfWidthOfGrid = abs(diff(gridEndPoints)) / 2;
widthArray = linspace(gridEndPoints(1),gridEndPoints(2),stim.widthOfGrid);% widthArray is used in creating the meshgrid.
monitor_dist = data.monitor_dist; %monitor to edge of desk (inches)
monitor_width = data.monitor_size; %width of screen (inches)
npixelswide = data.res_x; %width of screen (pixels)
npixelstall = data.res_y; %height of screen (pixels)
pixels_per_inch = npixelswide/monitor_width;
data.visual_angle_of_array = atan(stim.widthOfGrid/pixels_per_inch/monitor_dist)*180/pi;

%define center of gratings
y_shift = data.possible_shift_range*rand(1,size(stimmat,1)) - data.possible_shift_range/2;


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
for i=1:size(stimmat,1)
    % *** To rotate the grating, set tiltInDegrees to a new value.
    tiltInDegrees = stim.orient * 180/data.o_gain;
    tiltInDegrees = tiltInDegrees-90; % The tilt of the grating in degrees. COUNTERCLOCKWISE FROM VERTICAL

    arbitrary_sf = stimmat(i,stim.sf_col); %order matters: match sf outside stim range with actual orient and vice versa (in tiltInDegrees)
    scaled_sf = data.sf_gain(1)+(arbitrary_sf/data.sf_gain(2)); %CYCLES/DEGREE
    sf=scaled_sf * data.visual_angle_of_array; % spatial freq in cycles per image

	% ---------- Image Setup ----------
    %set up L and R gratings separately for ease of understanding
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
    
    disc_location = [where.xc- data.x_shift*stim.widthOfGrid where.yc+y_shift(1,i)];   %x and y coords of grating center
    
    rect = [disc_location - size(grayscaleImageMatrix1)/2     disc_location + size(grayscaleImageMatrix1)/2];

    temp = Screen('MakeTexture',window, grayscaleImageMatrix1);
    %temp = grayscaleImageMatrix1; %significant slow down if MakeTextures online
    eval(['grating.' structname num2str(stimmat(i,stim.index_col)) '.g = temp;']);
    eval(['grating.' structname num2str(stimmat(i,stim.index_col)) '.rect = rect;']);
end
