function h = ims(im,CLIM,TITLE,TILING)
%h = ims(im,CLIM,TITLE,TILING)
%
% CLIM:
%   [FRAC] = scale to fraction of the range
%   [CLOW CHIGH] = same as in imagesc
%   [0 1; NaN 100; ...] = individual ranges (NaN = no range)
%
%   Trickiness: if using FRAC, use sign to toggle colorbar  (negative=off)
%
% TITLE:
%   {'plot1', 'plot2', ...}
%
% TILING:
%   [N1 N2] = tiling for multi image display
%
% h returns handles to all the subplots

% add flexibility to TITLE (e.g. can use 1:N as titles)
if exist('TITLE','var')
    if isnumeric(TITLE)
        TITLE = num2str(reshape(TITLE,[],1));
    else
        % ignore certain formatting characters
        TITLE = strrep(TITLE,'_','\_');
    end
    if iscell(TITLE)
        TITLE = reshape(TITLE,[],1);
    end
else
    TITLE = '';
end

usecolorbar= true; % add colorbar

if ~exist('CLIM','var') || isempty(CLIM)
    CLIM = 0.99; % default to 99% of range
elseif numel(CLIM)==1
    if CLIM<0
        CLIM = -CLIM;        
        usecolorbar = false;
    end
    if CLIM==0 %|| CLIM>1
        error('CLIM (scalar) must be between 0 and 1.');
    end
elseif numel(CLIM)>2
    error('CLIM must be [low high; ... ].')
else
    if CLIM(2)<CLIM(1)
        CLIM = -CLIM; % assume it's the -ve flag for colorbar
        usecolorbar = false;
    end
end
if ~isreal(CLIM)
    error('CLIM must be real valued.')
end

% for complex produce 2 maps
if ~isreal(im)
    figure; ims(angle(im),[-pi pi],'phase'); % phase
    figure; im = abs(im); % magnitude
end

% make into 3D array
im = squeeze(im);
[x y n] = size(im);
if n==0
    error('Image array cannot be empty.')
end
im = reshape(im,x,y,n);

% try and figure out a nice tiling
if exist('TILING','var')
    if numel(TILING)~=2 || ~isnumeric(TILING) || any(mod(TILING,1))
        error('TILING must be a 2 vector of integers.');
    end
    if prod(TILING)<n || prod(TILING)>100
        error('TILING probably too large or too small (%ix%i).',TILING);
    end
else
    % catch likely error (too many images!)
    nmax = 16;
    if n>nmax
        mid = round(n/2-nmax/2);
        %warning('too many images (%i)... showing slices %i-%i only.',n,mid,mid+nmax)
        im = im(:,:,mid:mid+nmax);
        n = nmax;
    end
    rows = floor(sqrt(n));
    cols = ceil(n/rows);
    TILING = [rows cols];
end

    
% clear existing plots, unless a single image
if n>1
    clf reset
end

% plot here
for i = 1:TILING(1)
    for j = 1:TILING(2)
        
        % subplot number
        k = (i-1)*TILING(2)+j;
        if k>n; break; end
        if n>1 % retile (unless 1 image)
            h(k) = subplot(TILING(1),TILING(2),k);
        else
            h = gca;
        end

        % use double to prevent any matlab weirdness
        im_k = double(gather(im(:,:,k)));
        
        % display range (CLIM)
        lo = min(im_k(isfinite(im_k)));
        hi = max(im_k(isfinite(im_k)));
        
        if ~isempty(CLIM) && ~isempty(lo) && ~isempty(hi)
            if numel(CLIM)==1 % use fraction of the range
                edges = linspace(lo,hi,256);
                C = histc(im_k(isfinite(im_k)),edges);
                C = cumsum(C) / nnz(isfinite(im_k));
                [C IA] = unique(C);
                idx = interp1(C,IA,(1-CLIM)/2,'nearest','extrap');
                lo = interp1(IA,edges(IA),idx);
                idx = interp1(C,IA,1-(1-CLIM)/2,'nearest','extrap');
                hi = interp1(IA,edges(IA),idx);
            elseif size(CLIM,1)==1 % use same setting for all
                if ~isnan(CLIM(1,1));lo = CLIM(1,1);end
                if ~isnan(CLIM(1,2));hi = CLIM(1,2);end
            elseif k<=size(CLIM,1) % individual CLIMs
                if ~isnan(CLIM(k,1));lo = CLIM(k,1);end
                if ~isnan(CLIM(k,2));hi = CLIM(k,2);end
            end
        end
        if isequal(lo,hi)
            imagesc(im_k);
        else
            imagesc(im_k,[lo hi]);
        end
        set(gca,'XTickLabel','');
        set(gca,'YTickLabel','');
        set(gca,'TickLength',[0 0]);
        if usecolorbar; colorbar; end

        if ~isempty(TITLE)
            if size(TITLE,1)==1 % use same for all
                title(TITLE(1,:),'FontSize',10)
            else
                if k<=size(TITLE,1) % individual titles
                    if ~isequal(TITLE(k,:),'NaN')
                        title(TITLE(k,:),'FontSize',10);
                    end
                end
            end
        else
            if exist('mid','var') % user sent too many images, label them
                title(mid+k-1,'FontSize',10)
            end
        end
    end
end

% if extra space on the subplot, put focus on "next" one
if k>n && k<=TILING(1)*TILING(2)
    h(k) = subplot(TILING(1),TILING(2),k);
    axis off;
end

% prevent output to screen
if nargout==0; clear h; end
