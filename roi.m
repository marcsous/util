function [ind val] = roi(im,N,mytitle)
%function [ind val] = roi(im,N,mytitle)
% Interatively draw ROI on current figure (or on supplied image im).
% Displays number, mean and 95% confidence interval.
% Returns a cell array of vectors of indices for each roi.
% Returns a cell array of vectors of values for each roi.
% If present, N is the number of ROIs to draw.

% user chooses figure
if nargin==0 || isempty(im)
    h = 1;
    waitforbuttonpress;
else
    f = figure;
    h = ims(im);
    if numel(h)==1
        h = 1;
    else
        set(f,'name','Choose Subplot');
        waitforbuttonpress;
        for k = 1:numel(h)
            if isequal(h(k).Position,get(gca,'Position'))
                break;
            end
        end
        h = k; % the subplot number
    end
end

% grab data
hImage = imhandles(gca);
cdata = get(hImage,'cdata');

% make it bigger
%if ~exist('f','var'); f = figure; end
%clf; set(f,'name','Draw ROI');
%ims(cdata);
if exist('mytitle','var'); title(mytitle); drawnow; end

% loop until user selects an empty ROI (or N is reached)
if nargin<2; N = Inf; end

ind = {};
val = {};
while numel(ind)<N

    % draw roi
    bw = roipoly;
    ix = find(bw);
    if isempty(ix); break; end
    ind{end+1} = ix;
    val{end+1} = cdata(ix);

    % do calculations
    [B,BINT,R2,RINT,STATS] = regress(cdata(ix),ones(size(ix)));
    n = numel(ix);
    m = B;
    ci95 = B-BINT(1);
    stats = [n m std(cdata(ix)) ci95 median(cdata(ix)) median(abs(cdata(ix)-median(cdata(ix))))];
    disp(['n / mean / std / ci95 / median / mad :   ' num2str(stats,'%.4f\t')])

    % convert ind to include slice offset
    ind{end} = ind{end}+(h-1)*numel(cdata);
    
    % calm CPU down
    pause(0.05);
    
end
if exist('f','var'); close(f); end
if nargout<1; clear ind; end
