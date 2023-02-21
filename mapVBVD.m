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
%  varargin is one or more options (char)
%  'removeOS'/'~removeOS'
%  'doAverage'/'~doAverage'
%  'ignoreSeg'/'~ignoreSeg'
%  'mergeRef'/'~mergeRef'
%  'noiseAdj'/`~noiseAdj'

%% tricky! original mapVBVD function stored in subdirectory mapVBVD/
addpath(mfilename('fullpath'));
%% tricky! this will be removed at the end of the file 

%% always have a return argument
raw = [];
ref = [];
hdr = [];
param = [];

%% see if we have a filename
if nargin==0
    filename = '';
else
    if contains(lower(varargin{1}),'removeos')  || ...
       contains(lower(varargin{1}),'doaverage') || ...
       contains(lower(varargin{1}),'ignoreseg') || ...
       contains(lower(varargin{1}),'mergeref')  || ...
       contains(lower(varargin{1}),'noiseadj')
        filename = '';
    else
        filename = varargin{1};
        varargin(1) = [];
    end
end

% if filename is a directory
if isdir(filename) && ~contains(filename,'.dat')
    pathname = strcat(filename,filesep);
    filename = '';
elseif ~isempty(filename)
    if ~contains(filename,'dat'); filename = strcat(filename,'.dat'); end
    if ~exist(filename,'file'); error('File "%s" does not exist.',filename); end
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

mergeRef = true; % default
idx = strcmp(lower(varargin),'~mergeref');
if any(idx)
    mergeRef = false;
else
    idx = strcmp(lower(varargin),'mergeref');    
    if any(idx); mergeRef = true; end
end
varargin = varargin(~idx); % remove from varargin

noiseAdj = true; % default
idx = strcmp(lower(varargin),'~noiseadj');
if any(idx)
    noiseAdj = false;
else
    idx = strcmp(lower(varargin),'noiseadj');    
    if any(idx); noiseAdj = true; end
end
varargin = varargin(~idx); % remove from varargin

%% handle options passed to original mapVBVD
idx = strcmp(varargin,'~removeOS') | strcmp(varargin,'~removeos');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'removeOS'; % default
end
removeOS = ~any(idx);

idx = strcmp(varargin,'~doAverage') | strcmp(varargin,'~doAverage');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'doAverage'; % default
end
doAverage = ~any(idx);

idx = strcmp(varargin,'~ignoreSeg') | strcmp(varargin,'~ignoreSeg');
if any(idx)
    varargin = varargin(~idx);
else
    varargin{end+1} = 'ignoreSeg'; % default
end

% remove any dupes
varargin = unique(varargin);

%% show user the key

fprintf('%s(''%s'')\n',mfilename,filename);

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
if noiseAdj && isfield(twix_obj,'noise')
    noise = twix_obj.noise.unsorted();
    % concantenate multiple noise scans
    fprintf('Found %i noise scan(s).\n',size(noise,3));
    noise = reshape(permute(noise,[1 3 2]),[],size(noise,2));
    disp(['Noise std per coil (pre ): ' num2str(std(noise),' %.2e')]);
    
% sqrt(twix_obj.hdr.Meas.alDwellTime(1) * twix_obj.hdr.Meas.Averages)
% twix_obj.hdr.Meas.Averages
% twix_obj.hdr.Config.GlobalImageScaleFactor
% twix_obj.hdr.Config.NoiseScaleFactor
% twix_obj.hdr.Meas.GlobalImageScaleFactor
% twix_obj.hdr.Meas.dOverallImageScaleCorrectionFactor
% twix_obj.hdr.Meas.dOverallImageScaleFactor

    % noise correlation matrix
    [~,S,V] = svd(noise,'econ');
    X = V * pinv(S) * V';
    X = X * sqrt(mean(diag(S).^2)); % preserve scaling
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
    disp(['Noise std per coil (post): ' num2str(std(noise),' %.2e')]);
    end

%% pad array to correct final size... work in progress
ncols = size(raw,1); % actual no. readout points

try
    % how many similar-but-different variables do there have to be?!
    nx = hdr.Meas.lBaseResolution; %NImageCols;
    ny = hdr.Meas.PhaseEncodingLines; %NImageLins;
    ns = hdr.Meas.NSlc;
    nc = size(raw,2); % coils
    ne = size(raw,8); % echos
    nr = size(raw,9); % repetitions (half echos)
    
    % assume oversamping by factor of 2 in x-direction
    if ~removeOS; nx = 2*nx; end
    
    % not sure why nz is so flakey
    try
        nz = hdr.Meas.NPar;
    catch
        nz = hdr.Meas.NImagePar;
    end

    % oversampling factors
    %ox = hdr.Meas.lBaseResolution / size(raw,1);
    %oy = hdr.Meas.lPhaseEncodingLines / size(raw,3);
    %oz = hdr.Meas.lPartitions / size(raw,4);
    
    padsize(1) = nx-size(raw,1);
    padsize(2) = 0; % coils
    padsize(3) = ny-size(raw,3);
    padsize(4) = nz-size(raw,4);
    fprintf('Padsize [%s]\n',num2str(padsize));
    
    % if "interpolation" or oversampling in phase directions we need more hacking here
    if any(padsize<0)
        padsize = max(padsize,0); % skip the offending dimension(s)
        warning('Size mis-match padding final matrix - interpolation? parallel imaging? Ignore.');
    end
    
    % acceleration can leave odd dimension sizes...
    %if hdr.Meas.lAccelFactPE > 1
    %    ny = size(raw,3)+padsize(3); % final size
    %    padsize(3) = mod(ny,hdr.Meas.lAccelFactPE);
    %end
    %if hdr.Meas.lAccelFact3D > 1
    %    nz = size(raw,4)+padsize(4); % final size
    %    padsize(4) = mod(nz,hdr.Meas.lAccelFact3D);
    %end    
    
    raw = padarray(raw,padsize,'post'); % 'pre' sometimes...

    % fix miscentering of bipolar partial echo readouts
    if padsize(1)
       
        k = twix_obj.image.IsReflected;
        if doAverage
            k = k(twix_obj.image.Ave==1);
        end
        
        % seems IsReflected is in acqusition order (unsorted)?
        k = reshape(k,ne,ns,ny,nz,[]); % acquisition order
        k = permute(k,[3 4 2 1 5]); % same order as raw
        k = reshape(k,[],1); % back to linearized
        
        % need to use this twix_obj.image.centerCol?

        % this seems necessary with 3D data but not 2D...
        raw(:,:,k) = circshift(raw(:,:,k),-padsize(1));
        
    end

catch ME
    warning(ME.message)
end

%% return some helpful things
param.te = hdr.Meas.alTE*1e-6; % seconds
param.Tesla = hdr.Meas.flNominalB0; % Tesla
param.tr = hdr.Meas.alTR(1)*1e-6; % seconds
param.fa = hdr.Meas.FlipAngle; % degress
param.ave = hdr.Meas.lAverages;
param.dwell_us = double(hdr.Meas.alDwellTime(1)) / 1000; % us
param.ncols = ncols; % no. freq encodes (partial echo)

if exist('noise','var')
    param.std = std(noise(:)); % std over all coils
    if doAverage; param.std = param.std / sqrt(param.ave); end
    param.noise = noise; % raw noise if caller wants to be more sophisticated
end

if isfield(hdr.Meas,'EchoTrainLength') && ~isempty(hdr.Meas.EchoTrainLength)
    param.te = param.te(1:hdr.Meas.EchoTrainLength);
else
    param.te = param.te(1:nnz(param.te)); % hope this works
end

if isfield(hdr.Meas,'NoOfBvalues')
    param.Bvalues = hdr.Meas.alBValue(1:hdr.Meas.NoOfBvalues);
    param.NoDiffusionDirections = hdr.Phoenix.sDiffusion.lDiffDirections;
    param.NoDiffusionDirectionsOf1stBvalue = hdr.Meas.NoOfDirections4FirstBValue;   
    param.Bvalues = hdr.Meas.alBValue(1:hdr.Meas.NoOfBvalues);   
end

disp(['TR/TE ' num2str([param.tr param.te]*1e3,'%.2f ') ' ms'])

%% apply ICE scaling factor to put data in a slightly nicer range
K_ICE_AMPL_SCALE_FACTOR = 3200;

raw = raw * K_ICE_AMPL_SCALE_FACTOR;
if exist('ref','var'); ref = ref * K_ICE_AMPL_SCALE_FACTOR; end
if isfield(param,'std'); param.std = param.std * K_ICE_AMPL_SCALE_FACTOR; end
if isfield(param,'noise'); param.noise = param.noise * K_ICE_AMPL_SCALE_FACTOR; end

%% tricky (see top of file)
rmpath(mfilename('fullpath'));
%% tricky
