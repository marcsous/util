function [raw param ref hdr] = mapVBVD(varargin)
%function [raw param ref hdr] = mapVBVD(filename,options)
%
% Simplified wrapper to mapVBVD to read in Siemens .dat file.
%
% Instructions:
% 1. Put mapVBVD files in ./mapVBVD/
% 2. Rename mapVBVD.m to mapVBVD_original.m
% 3. Then you can call this function.
%
% Inputs:
%  filename is a char array (e.g. 'myfile.dat')
%
%  varargin is one or more options (char array)
%  'removeOS'/'~removeOS'
%  'doAverage'/'~doAverage'
%  'ignoreSeg'/'~ignoreSeg'
%  'mergeRef'/'~mergeRef'
%  'noiseAdj'/`~noiseAdj'

%% always have a return argument
raw = [];
ref = [];
hdr = [];
param = [];

%% see if we have a filename
if nargin==0
    filename = '';
else
    if isequal(varargin{1},'removeOS')  || isequal(varargin{1},'~removeOS')  || ...
       isequal(varargin{1},'doAverage') || isequal(varargin{1},'~doAverage') || ...
       isequal(varargin{1},'ignoreSeg') || isequal(varargin{1},'~ignoreSeg') || ...
       isequal(varargin{1},'mergeRef')  || isequal(varargin{1},'~mergeRef')  || ...
       isequal(varargin{1},'noiseAdj')  || isequal(varargin{1},'~noiseAdj')
        filename = '';
    else
        filename = varargin{1};
        varargin(1) = [];
    end
end

% if filename is a directory
if isdir(filename)
    pathname = strcat(filename,filesep);
    filename = '';
elseif ~isempty(filename) && ~exist(filename,'file')
    error('File "%s" does not exist.',filename)
else
    pathname = '';
end

% if filename not supplied, pop up a dialog box
if isempty(filename)
    [filename pathname] = uigetfile(strcat(pathname,'*.dat'), 'Pick a .dat file');
    if isequal(filename,0)
        disp('No file selected.')
    else
        filename = fullfile(pathname,filename);
    end
end

%% handle local options
idx = strcmp(varargin,'~mergeRef') | strcmp(varargin,'~MergeRef');
mergeRef = ~any(idx); % default is to mergeRef
varargin = varargin(~idx);

idx = strcmp(varargin,'~noiseAdj') | strcmp(varargin,'~NoiseAdj');
noiseAdj = ~any(idx); % default is to noiseAdj
varargin = varargin(~idx);

%% handle options passed to original mapVBVD
idx = strcmp(varargin,'~removeOS') | strcmp(varargin,'~removeos');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'removeOS'; % default is to removeOS
end
removeOS = ~any(idx); % we use this locally

idx = strcmp(varargin,'~doAverage') | strcmp(varargin,'~doAverage');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'doAverage'; % default is to doAverage
end

idx = strcmp(varargin,'~ignoreSeg') | strcmp(varargin,'~ignoreSeg');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'ignoreSeg'; % default is to ignoreSeg
end

% remove any dupes
varargin = unique(varargin);

%% show user the key

fprintf('%s(''%s''\n',mfilename,filename);

if nargout < 2
    disp( ['Order of raw data dimensions:'])
    disp( '   1) Columns')
    disp( '   2) Channels/Coils')
    disp( '   3) Lines')
    disp( '   4) Partitions')
    disp( '   5) Slices')
    disp( '   6) Averages')
    disp( '   7) (Cardiac-) Phases')
    disp( '   8) Contrasts/Echoes')
    disp( '   9) Measurements')
    disp( '  10) Sets')
    disp( '  11) Segments')
    disp( '  12) Ida')
    disp( '  13) Idb')
    disp( '  14) Idc')
    disp( '  15) Idd')
    disp( '  16) Ide')
end

%% use mapVBVD function from Siemens IDEA Tools
addpath(mfilename('fullpath')); % tricky! folder must have same name as mfile

twix_obj = mapVBVD_original(filename,varargin{:});

% as per mapVBVD tutorial, VE11 has 2 twix objects (annoying).
% the first one is noise scans - possibly duplicated in the
% second one but who knows? to keep life simple, implant the
% noise scans in the second twix and get rid of the first one.
% if there are more than 2 then we are in a weird place.
if numel(twix_obj)>1
    if numel(twix_obj)>2
        warning('something weird is happening: %i twix objects. only using the last one.',numel(twix_obj));
    end
    if isfield(twix_obj{1},'noise') && ~isfield(twix_obj{2},'noise')
        twix_obj{end}.noise = twix_obj{1}.noise;
    end
    twix_obj = twix_obj{end};
end

fprintf('Reading raw data... ')
raw = twix_obj.image();
hdr = twix_obj.hdr;
if isfield(twix_obj,'refscan')
    ref = twix_obj.refscan.unsorted();
end
fprintf('done\n')

%% merge ref and raw data
if isfield(twix_obj,'refscan') && mergeRef
    fprintf('Sorting ref data... ')
    ref_sorted = zeros(twix_obj.refscan.dataSize,class(ref));
    for j = 1:twix_obj.refscan.NAcq
        Lin = twix_obj.refscan.Lin(j);
        Par = twix_obj.refscan.Par(j);
        Sli = twix_obj.refscan.Sli(j);
        Ave = twix_obj.refscan.Ave(j);
        Eco = twix_obj.refscan.Eco(j);
        Phs = twix_obj.refscan.Phs(j);
        Rep = twix_obj.refscan.Rep(j);
        % not sure about skipping in general but it is needed in my dat files - MB
        if twix_obj.refscan.flagSkipToFirstLine
            Lin = Lin-min(twix_obj.refscan.Lin)+1;
            Par = Par-min(twix_obj.refscan.Par)+1;
            Sli = Sli-min(twix_obj.refscan.Sli)+1;
            Ave = Ave-min(twix_obj.refscan.Ave)+1;
            Eco = Eco-min(twix_obj.refscan.Eco)+1;
            Phs = Phs-min(twix_obj.refscan.Phs)+1;
            Rep = Rep-min(twix_obj.refscan.Rep)+1;
        end
        ref_sorted(:,:,Lin,Par,Sli,Ave,Phs,Eco,Rep) = ref(:,:,j);
    end
    fprintf('done\n')
    
    if mergeRef
        fprintf('Merging ref data... ')
        try
            for j = 1:twix_obj.refscan.NAcq
                Lin = twix_obj.refscan.Lin(j);
                Par = twix_obj.refscan.Par(j);
                Sli = twix_obj.refscan.Sli(j);
                Ave = twix_obj.refscan.Ave(j);
                Eco = twix_obj.refscan.Eco(j);
                Phs = twix_obj.refscan.Phs(j);
                Rep = twix_obj.refscan.Rep(j);
                if nnz(raw(:,:,Lin,Par,Sli,Ave,Phs,Eco,Rep))
                    if nnz(raw(:,:,Lin,Par,Sli,Ave,Phs,Eco,Rep)-ref(:,:,j))
                        error('ref and raw data not identical (line %i).',j)
                    end
                else
                    raw(:,:,Lin,Par,Sli,Ave,Phs,Eco,Rep) = ref(:,:,j);
                end
            end
            fprintf('done\n')
        catch ME
            fprintf('failed (%s)\n',ME.message)
        end
    end
    ref = ref_sorted;
end

%% factor out noise correlation
if isfield(twix_obj,'noise') && noiseAdj
    noise = twix_obj.noise.unsorted();
    % concantenate multiple noise scans
    fprintf('Found %i noise scans.\n',size(noise,3));
    noise = reshape(permute(noise,[1 3 2]),[],size(noise,2));
    disp(['Noise std per coil: ' num2str(std(noise),' %.1e')]);
    % noise correlation matrix
    [~,S,V] = svd(noise,'econ');
    X = V * inv(S);% * V' * sqrt(size(noise,1));
    for j = 1:numel(raw)/size(raw,1)/size(raw,2)
        raw(:,:,j) = raw(:,:,j) * X;
    end
    if ~isempty(ref)
        for j = 1:numel(ref)/size(ref,1)/size(ref,2)
            ref(:,:,j) = ref(:,:,j) * X;
        end
    end
    % correct the noise to check std is normalized
    noise = noise * X;
    disp(['Noise std per coil: ' num2str(std(noise),' %.1e')]);
end

%% pad array to correct final size... work in progress
try
    nx = hdr.Meas.NImageCols;
    ny = hdr.Meas.NImageLins;
    ns = hdr.Meas.NSlc;

    % oversamping by factor of 2
    if ~removeOS; nx = 2*nx; end

    % not sure why nz is flakey
    if isfield(hdr.Meas,'NImagePar')
        nz = hdr.Meas.NImagePar;
    else
        nz = hdr.Meas.NPar;
    end
    
    padsize(1) = nx-size(raw,1);
    padsize(2) = 0; % coils
    padsize(3) = ny-size(raw,3);
    padsize(4) = nz-size(raw,4);
    
    if any(padsize<0)
        padsize = max(padsize,0); % skip the offending dimension(s)
        warning('Size mis-match padding final matrix - skipping.');
    end
    raw = padarray(raw,padsize,'pre');
    
    % for bipolar readout directions - not well tested
    if padsize(1)
        [~,nk] = size(raw);
        for k = 1:nk
            if twix_obj.image.IsReflected(k)
                raw(:,k) = circshift(raw(:,k),-padsize(1));
            end
        end
    end
    
catch ME
    warning(ME.message)
end

rmpath('~/Documents/MATLAB/library/mapVBVD');


%% return some helpful things
param.te = hdr.Meas.alTE*1e-6; % seconds
param.Tesla = hdr.Meas.flNominalB0; % Tesla

if isfield(hdr.Meas,'EchoTrainLength') && ~isempty(hdr.Meas.EchoTrainLength)
    param.te = param.te(1:hdr.Meas.EchoTrainLength);
else
    param.te = param.te(1:nnz(param.te)); % hope this works
end

