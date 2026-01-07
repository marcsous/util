function [x y s] = roi_points(im,radius)
%[x y s] = roi_points(im,radius)
% Receives an stack of slices for user to drop ROI points onto.
% Left click to drop an ROI. Right click to finish. Mouse wheel
% scrolls through slices. Returns x, y and slice locations.
%
% -im is a stack of 2D images [nx ny ns]
% -radius is the no. pixels to display drawn ROIs [5]

[nx,ny,ns,~] = size(im);

if nargin<2
    radius = 5;
elseif ~isscalar(radius)
    error('radius must be a scalar');
end

% anything larger than 3 dimensions is averaged
im = mean(reshape(im,nx,ny,ns,[]),4);

figure(gcf); clf reset; subplot(1,1,1);

button = 0;
slice = ceil(ns/2);
x = []; y = []; s = [];

while button~=3 && button~=27 % right click && Esc key

    imagesc(sqrt(im(:,:,slice,:))); % sqrt to level the range a bit
    title('left click to select, right click to end',slice);
    
    drawnow;
    [myy myx button] = ginputc(1);
    drawnow;

    switch button

        case 1; % left click
                x(end+1) = round(myx);
                y(end+1) = round(myy);
                s(end+1) = slice;

        case 4; % roll wheel up
                slice = min(slice+1,ns);

        case 5; % roll wheel down
                slice = max(slice-1,1);
    
    end

    % draw ROI on the image
    for j = -radius:radius
        for k = -radius:radius
            if hypot(j,k)<radius+0.5 && hypot(j,k)>radius-0.5
                im(x(end)+j,y(end)+k,slice) = NaN;
            end
        end
    end

end
