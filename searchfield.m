function searchfield(str,txt,show)
%function searchfield(str,txt,show)
% search within a structure
%
% str = structure/object
% txt = text to find (char)
% show = show value (true/false)

if ~isstruct(str) && ~isobject(str); error('str must be a struct/object'); end
if ~isa(txt,'char'); error('txt must be char'); end
if ~exist('show'); show=false; end
if ~isequal(show,1) && ~isequal(show,0); error('show must be 0 or 1'); end

global TXT;
TXT = txt;

unfold(str,inputname(1),show==1);

%% 
%
% https://www.mathworks.com/matlabcentral/answers/94122-is-there-a-function-that-will-display-all-the-fields-of-a-nested-structure-in-matlab

function unfold(SC,varargin)
%UNFOLD Unfolds a structure.
%   UNFOLD(SC) displays the content of a variable. If SC is a structure it
%   recursively shows the name of SC and the fieldnames of SC and their
%   contents. If SC is a cell arraythe contents of each cell are displayed.
%   It uses the caller's workspace variable name as the name of SC. 
%   UNFOLD(SC,NAME) uses NAME as the name of SC.
%   UNFOLD(SC,SHOW) If SHOW is false only the fieldnames and their sizes
%   are  shown, if SHOW is true the contents are shown also.
%   UNFOLD(SC,NAME,SHOW)

%R.F. Tap
%15-6-2005, 7-12-2005, 5-1-2006, 3-4-2006

global TXT;

% check input
switch nargin
    case 1
        Name = inputname(1);
        show = true;
    case 2
        if islogical(varargin{1})
            Name = inputname(1);
            show = varargin{1};
        elseif ischar(varargin{1})
            Name = varargin{1};
            show = true;
        else
            error('Second input argument must be a string or a logical')
        end
    case 3
        if ischar(varargin{1})
            if islogical(varargin{2})
                Name = varargin{1};
                show = varargin{2};
            else
                error('Third input argument must be a logical')
            end
        else
            error('Second input argument must be a string')
        end
    otherwise
        error('Invalid number of input arguments')
end


if isstruct(SC) || isobject(SC)
    %number of elements to be displayed
    NS = numel(SC);
    if show
        hmax = NS;
    else
        hmax = min(1,NS);
    end

    %recursively display structure including fieldnames
    for h=1:hmax
        %if hmax~=1; warning('something weird... MB'); continue; end % MB  
        F = fieldnames(SC); % MB SC(h) fails on objects...!? so reject hmax>1
        NF = length(F);
        for i=1:NF
            if NS>1
                siz = size(SC);
                if show
                    Namei = [Name '(' ind2str(siz,h) ').' F{i}];
                else
                    Namei = [Name '(' ind2str(siz,NS) ').'  F{i}];
                end
            else
                Namei = [Name '.' F{i}];
            end
            if isstruct(SC(h).(F{i}))
                unfold(SC(h).(F{i}),Namei,show);
            else
                if iscell(SC(h).(F{i}))
                    siz = size(SC(h).(F{i}));
                    NC = numel(SC(h).(F{i}));
                    if show
                        jmax = NC;
                    else
                        jmax = 1;
                    end
                    for j=1:jmax
                        if show
                            Namej = [Namei '{' ind2str(siz,j) '}'];
                        else
                            Namej = [Namei '{' ind2str(siz,NC) '}'];
                        end
                        if ~isempty(strfind(lower(Namej),lower(TXT)))
                            fprintf('%s',Namej)
                            if show && ~isempty(SC(h).(F{i}){j})
                                disp(SC(h).(F{i}){j});
                            else
                                fprintf('\n');
                            end
                       end
                    end
                else
                    if ~isempty(strfind(lower(Namei),lower(TXT)))
                        fprintf('%s',Namei)
                        if show && ~isempty(SC(h).(F{i}))
                            disp(SC(h).(F{i}))
                        else
                            fprintf('\n');
                        end
                    end
                end
            end
        end
    end
elseif iscell(SC)
    %recursively display cell
    siz = size(SC);
    for i=1:numel(SC)
        Namei = [Name '{' ind2str(siz,i) '}'];
        unfold(SC{i},Namei,show);
    end
else
    warning('something went wrong')
    if show
        disp(SC)
    end
end


%local functions
%--------------------------------------------------------------------------
function str = ind2str(siz,ndx)

n = length(siz);
%treat vectors and scalars correctly
if n==2
    if siz(1)==1
        siz = siz(2);
        n = 1;
    elseif siz(2)==1
        siz = siz(1);
        n = 1;
    end
end
k = [1 cumprod(siz(1:end-1))];
ndx = ndx - 1;
str = '';
for i = n:-1:1,
    v = floor(ndx/k(i))+1;
    if i==n
        str = num2str(v);
    else
        str = [num2str(v) ',' str];
    end
    ndx = rem(ndx,k(i));
end
