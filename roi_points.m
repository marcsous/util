function [x y s] = roi_points(im,radius)
%[x y s] = roi_points(im,radius)
% Receives an stack of slices for user to drop ROI points onto.
% Left click to drop an ROI. Right click to finish. Mouse wheel
% scrolls through slices. Esc to quit. X to delete.
%
% Returns x, y and slice locations.
%
% -im is a stack of 2D images [nx ny ns]
% -radius is the no. pixels to display drawn ROIs [5]

[nx,ny,ns,~] = size(im);

if nargin<2
    radius = 5;
elseif ~isscalar(radius)
    error('radius must be a scalar');
end

% anything above 3 dimensions is averaged. sqrt to level the dynamic range a bit
im = sqrt(mean(reshape(im,nx,ny,ns,[]),4));

figure(gcf); clf reset; subplot(1,1,1);

button = 0;
slice = ceil(ns/2);
x = []; y = []; s = []; old = [];

while button~=3 && button~=27 % right click && Esc key

    figure(gcf); imagesc(im(:,:,slice,:));
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
                switch button
                    case 1;
                        old(end+1) = im(x(end)+j,y(end)+k,slice);
                        im(x(end)+j,y(end)+k,s(end)) = NaN;
                    case {88,120};
                        if numel(old)>0
                            im(x(end)-j,y(end)-k,s(end)) = old(end);
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
