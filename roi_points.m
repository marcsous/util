function [x y s v] = roi_points(im,radius)
%[x y s v] = roi_points(im,radius)
% Receives an stack of slices for user to drop ROI points onto.
% Returns x, y and slice locations, and mean values.
%
% -im is a stack of 2D images [nx ny ns]
% -radius is the no. pixels for ROIs [5]
%
% Controls:
% Left click to drop an ROI. Right click to quit. Mouse wheel to
% scroll slices. Esc to quit. Q to quit. X to delete. B to brighten.
% D to darken.

[nx ny ns ne] = size(im);

if ne~=1
    error('im can only be 2D or 3D array');
end
if nargin<2
    radius = 5;
elseif ~isscalar(radius)
    error('radius must be a scalar');
end

%% image for display of ROIs
imd = double(im);

figure(gcf); clf reset; subplot(1,1,1);

%% capture user events

button = 0;
slice = ceil(ns/2);
x = []; y = []; s = []; old = [];

% quit with right click, Esc, q or Q
while ~ismember(button,[3 27 81 113])

    % display current slice
    figure(gcf); imagesc(imd(:,:,slice));
    title('left click to select, right click to quit',slice);

    % capture stray keyboard events that otherwise go to the commandline
    set(gcf,'WindowKeyPressFcn',@(src, event) disp(''));
    
    % wait for user input
    drawnow;
    [myy myx button] = ginputc(1);
    drawnow;

    % handle user input
    switch button

        case 1; % left click
                x(end+1) = round(myx);
                y(end+1) = round(myy);
                s(end+1) = slice;

        case 4; % roll wheel up
                slice = min(slice+1,ns);

        case 5; % roll wheel down
                slice = max(slice-1,1);
    
        case {66,98} % brighten (b or B)
                brighten(+0.1);
                
        case {68,100} % darken (d or D)
                brighten(-0.1);

    end

    % draw or delete ROI on the image
    for j = -radius:radius
        for k = -radius:radius
            if hypot(j,k)<radius+0.5 && hypot(j,k)>radius-0.5
                switch button
                    case 1; % left click
                        old(end+1) = imd(x(end)+j,y(end)+k,slice);
                        imd(x(end)+j,y(end)+k,s(end)) = NaN;
                    case {88,120}; % x or X
                        if numel(old)>0
                            imd(x(end)-j,y(end)-k,s(end)) = old(end);
                            old(end) = [];
                        end
                end
            end
        end
    end

    % now can delete point
    if button==88 || button==120 % x or X
        x = x(1:end-1);
        y = y(1:end-1);
        s = s(1:end-1);
    end

end

%% extract rois

v = [];
vs = [];
range = -radius:radius;

mask = false(2*radius+1,2*radius+1);
for j = -radius:radius
    for k = -radius:radius
        if hypot(j,k)<=radius
            mask(1+radius+j,1+radius+k) = true;
        end
    end
end

for n = 1:numel(x)

    tmp = double(im(x(n)+range,y(n)+range,s(n)));

    v(n) = mean(tmp(mask));
    vs(n) = std(tmp(mask));

end
