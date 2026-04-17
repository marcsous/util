function [x y s v] = roi_points(im,radius)
%[x y s v] = roi_points(im,radius)
% Receives a stack of slices for user to drop ROI points onto.
% Returns x, y and (s)lice locations, and mean (v)alues.
%
% -im is a stack of 2D images [nx ny ns]
% -radius is the no. pixels for ROIs [5]
%
% Controls:
% -Left click to drop an ROI
% -Right click to quit
% -Mouse wheel to scroll slices
% -Esc to quit
% -Q to quit
% -X to delete last ROI
% -B to brighten
% -D to darken

[nx ny ns ne] = size(im);

if ne~=1 % no extra dimesions allowed
    error('im can only be 2D or 3D array');
end
if nargin<2
    radius = 5;
elseif ~isscalar(radius) || ~isfinite(radius)
    error('radius must be a scalar');
end

%% initialize

button = 0;
slice = ceil(ns/2);
range = -radius:radius;
x = []; y = []; s = []; v = [];

% create circular mask
mask = false(numel(range));
for j = range
    for k = range
        if hypot(j,k)<=radius
            mask(1+radius+j,1+radius+k) = true;
        end
    end
end

% reset figure
figure(gcf); clf reset; subplot(1,1,1);

% capture keyboard events that otherwise go to the commandline
set(gcf,'WindowKeyPressFcn',@(src, event) disp(''));

%% capture user events

% quit with right click, Esc, q or Q
while ~ismember(button,[3 27 81 113])

    % update image
    figure(gcf); imagesc(im(:,:,slice));
    title('left click to select, right click to quit',slice);

    % draw ROIs on the image
    for n = 1:numel(x)

        % mean value inside mask
        tmp = im(x(n)+range,y(n)+range,s(n));
        v(n) = mean(tmp(mask));

        % draw ROI on the image
        if slice==s(n)
            text(y(n),x(n),'◯','Fontsize',2*radius+1,'HorizontalAlignment','center','Color','r');
            if     abs(v(n))>10; fmt = '%.0f';
            elseif abs(v(n))> 1; fmt = '%.1f';
            else fmt = '%.2f'; end
            text(y(n)+radius,x(n)+1,num2str(v(n),fmt),'Color','r','FontSize',14);
        end
    end

    % wait for user input
    [myy myx button] = ginputc(1,'Hide',true);

    % handle user input
    switch button

        case 1; % left click
            if round(myx)>=1 && round(myy)>=1 && round(myx)<=nx && round(myy)<=ny
                x(end+1) = round(myx);
                y(end+1) = round(myy);
                s(end+1) = slice;
            end

        case 4; % roll wheel up
                slice = max(slice-1,1);

        case 5; % roll wheel down
                slice = min(slice+1,ns);
    
        case {66,98} % brighten (b or B)
                brighten(+0.1);
                
        case {68,100} % darken (d or D)
                brighten(-0.1);

        case {88,120} % delete point (x or X)
                if numel(x)>=1
                    x = x(1:end-1);
                    y = y(1:end-1);
                    s = s(1:end-1);
                    v = v(1:end-1);
                end
        end

end
