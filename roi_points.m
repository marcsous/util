function [x y s v] = roi_points(im,radius,txt,clim)
%[x y s v] = roi_points(im,radius,txt,clim)
% Receives a stack of slices for user to drop ROI points onto.
% Returns x, y and (s)lice locations, and mean (v)alues.
%
% -im is a stack of 2D images [nx ny ns]
% -radius in pixels for circular ROIs [7]
% -txt for title ('  ' flexible spacing) ['']
% -clim is passed to imagesc []
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
% -P to print (~/Desktop/image.jpg)

[nx ny ns ne] = size(im);

if ne~=1 % extra dims not allowed
    if ns==1 % except if ns==1
        [ns ne] = deal(ne,ns);
        im = reshape(im,[nx ny ns]);
    else
        error('im can only be 2D or 3D array');
    end
end
if nargin<2 || isempty(radius)
    radius = 7;
elseif ~isscalar(radius) || ~isfinite(radius)
    error('radius must be a scalar');
end
if nargin<3
    txt = '';
end
if nargin<4
    clim = [-Inf Inf];
elseif ~isnumeric(clim) || numel(clim)<1 || numel(clim)>2
    error('argument ''clim'' invalid');
end

%% initialize

button = 0;
slice = ceil(ns/2);
x = []; y = []; s = []; v = [];

% reset figure
figure(gcf); clf reset; subplot(1,1,1);

% capture keyboard events that otherwise go to the commandline
set(gcf,'WindowKeyPressFcn',@(~,~) disp(''));

%% capture user events

% quit with right click, Esc, q or Q
while ~ismember(button,[3 27 81 113])

    % update image
    figure(gcf); imagesc(im(:,:,slice),clim); colorbar; %title(txt,slice);
 
    % stretch title (treat '  ' as flexible spacing)
    n = numel(strfind(txt,'  '));
    if n==0
        title(txt,slice);
    else
        nspaces = 3; % 2 spaces => nspaces
        h = title(txt,slice,'Units','Normalized');
        while h.Extent(3)<n/(n+1)
            stxt = strrep(txt,'  ',repmat(' ',1,nspaces));
            h = title(stxt,slice);
            nspaces = nspaces+1;
        end
    end

    % draw ROIs on the image
    for n = 1:numel(s)
        if slice==s(n)
            h = drawcircle('Color','red','FaceAlpha',0,'Radius',radius,'Center',[y(n) x(n)],'LineWidth',1,'InteractionsAllowed','none');
            
            mask = createMask(h);
            tmp = im(:,:,slice);
            v(n) = mean(tmp(mask(:)));

            if     abs(v(n))>10; fmt = '%.0f';
            elseif abs(v(n))> 1; fmt = '%.1f';
            else fmt = '%.2f'; end
            text(y(n)+radius+1,x(n)+1,num2str(v(n),fmt),'Color','r','FontSize',14);
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

        case {80,112} % print (p or P)
            [~] = lpr('~/Desktop/image.jpg');

    end

    % validation - etch ROI into the image
    %if button==1
    %    for j = -radius:radius
    %        for k = -radius:radius
    %            if hypot(j,k)<radius+0.5 && hypot(j,k)>radius-0.5
    %                im(x(end)+j,y(end)+k,s(end)) = NaN;
    %            end
    %        end
    %    end
    %end

end
