function [data params] = get_dicom(token,pathname,maxno)
% Usage: [data params] = get_dicom(pathname)
% Usage: [data params] = get_dicom(token,pathname,maxno)
%
% Locates all DICOM files in the path specified by pathname and filtered by token.
% Recurses into subdirectories.
%
% Inputs:
%  token is an identifier to parse the DICOM header:
%   -if numeric then we assume it's a seriesNumber
%   -if string then we assume it's SeriesDescription
%   -if struct then we assume it's returned from loaddcmdir.m
%  pathname is the path to search in (searches recursively)
%  maxno is the max. number of files to read in
%
% Outputs:
%  data is a high dimensional array sorted by parameters
%  params is a structure of parameters (e.g. params.te is te in seconds)

%% handle arguments
if ~exist('token','var')
    token = []; % empty string also acceptable
elseif exist(token,'dir')
    pathname = token;
    token = [];
end
if ~exist('pathname','var') || isempty(pathname)
    pathname = uigetdir();
    if isequal(pathname,0); return; end
end
if isa(token,'struct')
    pathname = token.Path;
end

% convert file separators for local machine
pathname(pathname=='\' | pathname=='/') = filesep;
if pathname(end) ~= filesep
    pathname(end+1) = filesep;
end
if ~exist(pathname,'dir')
    error(['pathname does not exist: ' pathname])
end
if ~exist('maxno','var') || isempty(maxno)
    maxno = Inf;
end
disp([mfilename '(): ' pathname])

% parse token argument
if isa(token,'struct')
    disp([mfilename '(): reading in files specified by loaddcmdir(). Found 0   '])
    info = cell2struct(token.Images,'name',numel(token.Images));
else
    if isempty(token)
        fprintf([mfilename '(): Reading in all files. Found 0   '])
    elseif isa(token,'char')
        fprintf([mfilename '(): Parsing SeriesDescription with ''' token '''. Found 0   '])
    elseif isa(token,'numeric')
        fprintf([mfilename '(): Parsing SeriesNumber with ''' num2str(token) '''. Found 0   '])
    end
    info = get_files(pathname); % recurse into directories
end

%% read in data

% create file counter using temp file (workaround for parfor)
tmpfile = [tempname '.delete.me'];
fid = fopen(tmpfile,'w');
fclose(fid);

for j = 1:numel(info) % much faster with parfor

    % current file
    filename = [pathname info(j).name];

    % test for maxno and whether it's a DICOM file
    if update_counter(tmpfile,maxno) & isdicom(filename)

        % read header
        head = dicominfo(filename);

        % filter by property
        if isempty(token) || isa(token,'struct')
            found = true;
        elseif isa(token,'char')
            found = ~isempty(strfind(head.SeriesDescription,token));
            %keyboard
            %found = ~isempty(strfind(head.PatientName,token));
        elseif isa(token,'numeric')
            found = head.SeriesNumber==token;
        else
            found = false;
        end

        % skip "REPORT" dicoms
        if isequal(head.Modality,'SR')
            found = false;
        end

        % store properties
        if found

            % update count
            update_counter(tmpfile);

            % read tags for sorting / checking
            if isempty(head.EchoTime)
                TE(j) = -1;
            else
                TE(j) = head.EchoTime;
            end
            SL(j) = head.SliceLocation;
            FA(j) = head.FlipAngle;
            TR(j) = head.RepetitionTime;
            IN(j) = head.InstanceNumber;
            NEX(j)= head.AcquisitionNumber;
            LF(j) = head.ImagingFrequency;

            %SN(j) = head.SeriesNumber;
            SE{j} = head.SeriesInstanceUID;
            ST{j} = head.StudyInstanceUID;

            % real/imag flags: 0=mag 1=phase 2=real 3=imag
            RI(j) = 0; % default if no tag is present
            if isfield(head,'Private_0043_102f') % GE
                RI(j) = head.Private_0043_102f(1);
            elseif isfield(head,'ImageType') % Siemens (used to be Private_0051_1016)
                % look for 'P' or 'M' (phase/mag) etc.
                if strfind(head.ImageType,'\M\')
                    RI(j) = 0;
                elseif strfind(head.ImageType,'\P\')
                    RI(j) = 1;
                elseif strfind(head.ImageType,'\R\')
                    RI(j) = 2;
                elseif strfind(head.ImageType,'\I\')
                    RI(j) = 3;
                else
                    warning('RI handling not working (value = %s).',head.ImageType);
                    RI(j) = 4;
                end
            end

            % inversion tag not always present
            if isfield(head,'InversionTime')
                TI(j) = head.InversionTime;
            end

            % b-values not stored in a standard tag
            if isequal(head.Manufacturer,'GE')
                if isfield(head,'Private_0043_1039')
                    %tensor = [head.Private_0019_10bb ...
                    %          head.Private_0019_10bc ...
                    %          head.Private_0019_10bd];

                    % weird GE thing to add 1e9 except when b=0
                    BV(j) = max(0,head.Private_0043_1039(1)-1e9);
                end
            elseif isequal(head.Manufacturer,'SIEMENS')
                if isfield(head,'Private_0019_100c')
                    BV(j) = head.Private_0019_100c;
                    % head.Private_0019_100d supposedly holds the direction
                end
            end

            % store images in cell array (less memory)
            data{j} = dicomread(filename);

            % store header
            dcm{j} = head;

        end % is found

    end % is dicom

end

% clean up
delete(tmpfile);

% check for empty
if ~exist('data')
    data = [];
    params = [];
    dcm = [];
    return;
end

% count non-empty data cells
index = find(~cellfun('isempty',data));
count = numel(index);

% counter increments in multiples of ncpu, trim excess
if isfinite(maxno)
    count = maxno;
    index = index(1:maxno);
end

% parfor counter is flaky so re-display actual count
fprintf('\b\b\b\b%-4d\n',count);

% check image dims are compatible - cannot solve this
[nx ny] = size(data{index(1)});
for j = 1:count
    sz = size(data{index(j)});
    if ~isequal(sz,[nx ny])
        warning('size mis-match image %i (expecting %ix%i found %ix%i)',j,[nx ny],sz);
        keyboard
    end
end

% convert cell arrays to matrices
dcm = dcm(index);
data = data(index);
data = cell2mat(data);

data = reshape(data,nx,ny,count); % image data
TE = TE(index); % echo time
SL = SL(index); % slice location
FA = FA(index); % flip angle
TR = TR(index); % repetition time
IN = IN(index); % instance number
%SN = SN(index); % series number
SE = SE(index); % series uid
ST = ST(index); % study number
RI = RI(index); % rawdata type (real/imag)
LF = LF(index); % larmor freq
NEX = NEX(index); % no. excitations (averages)
if exist('TI','var'); TI = TI(index); else TI = []; end % inversion time
if exist('BV','var'); BV = BV(index); else BV = []; end % bvalue

% can't handle mixture of studies
if numel(unique(ST)) > 1
    error('More than one study UID present - cannot proceed. Try dicomsorter.');
end

% unique properties
[uTE,~,kte] = unique(TE);
[uSL,~,ksl] = unique(round(SL*1e3)/1e3); % round floats (sometimes values have jitter)
[uFA,~,kfa] = unique(FA);
[uTR,~,ktr] = unique(TR);
%[uSN,~,ksn] = unique(SN);
[uSE,~,kse] = unique(SE);
[uTI,~,kti] = unique(TI);
[uBV,~,kbv] = unique(BV);
[uRI,~,kri] = unique(RI);
[uLF,~,klr] = unique(round(LF*1e3)/1e3); % round floats (sometimes values have jitter)
[uNEX,~,knex] = unique(NEX);

nTE = numel(uTE);
nSL = numel(uSL);
nFA = numel(uFA);
nTR = numel(uTR);
%nSN = numel(uSN);
nSE = numel(uSE);
nTI = numel(uTI);
nBV = numel(uBV);
nRI = numel(uRI);
nNEX = numel(uNEX);
nLF = numel(uLF);
if nTI==0; nTI = 1;	kti = ones(count,1); end
if nBV==0; nBV = 1;	kbv = ones(count,1); end
if nLF~=1; error('More than one field strength present!'); end

% try to be intelligent about multiple series
if isa(token,'numeric') && ~isempty(token)
    ignore_SERIES = false;
else
    ignore_SERIES = true;
end

% Siemens stores mag/phase in consecutive series
if ~isempty(strfind(dcm{1}.Manufacturer,'SIEMENS')) && nRI>1
    ignore_SERIES = true;
end

if ignore_SERIES && nSE>1
    disp([mfilename '(): Combining ' num2str(nSE) ' series. Specify SeriesNumber to avoid this.'])
    [~,k] = unique(kse);
    for j = 1:nSE
        if isfield(dcm{k(j)},'SeriesDescription')
            disp(['  ' num2str(j) '  ' dcm{k(j)}.SeriesDescription])
        end
    end
    nSE = 1;
    uSE = 1;
    kse(:) = 1;
end

% check some dicom tags
% tag = {'StudyInstanceUID' 'SeriesInstanceUID'};
% for j = 1:numel(dcm)
%     error_found = 0;
%     for k = 1:numel(tag)
%         if ~isequal(getfield(dcm{1},tag{k}),getfield(dcm{j},tag{k}))
%             fprintf('DICOM tag mismatch (%s)\n',tag{k});
%             fprintf('tag1\n'); disp(getfield(dcm{1},tag{k}));
%             fprintf('tag2\n'); disp(getfield(dcm{j},tag{k}));
%             error_found = 1; break;
%         end
%         if error_found; break; end
%     end
% end

% error checks
expected_count = nTE*nSL*nFA*nTR*nTI*nBV*nRI*nSE*nNEX;
if count~=expected_count
    warning('%s() wrong number of images (found %i but expecting %i)',mfilename,count,expected_count)
    disp(['CHECK THESE!! echos=' num2str(nTE) ...
        ' slices=' num2str(nSL) ...
        ' flips=' num2str(nFA) ...
        ' TRs=' num2str(nTR) ...
        ' TIs=' num2str(nTI) ...
        ' Bs=' num2str(nBV) ...
        ' series=' num2str(nSE) ' (ignore=' num2str(ignore_SERIES) ')'...
        ' real/imag=' num2str(nRI) ...
        ' averages=' num2str(nNEX)])
    disp('Type "return" to ignore (expect errors!) or dbquit to cancel.')
keyboard
end

% sort by property
temp = zeros(nx,ny,nSL,nTE,nFA,nTR,nTI,nBV,nSE,nRI,nNEX,'single');
dcmtemp = cell(nSL,nTE,nFA,nTR,nTI,nBV,nSE,nRI,nNEX);
mask = zeros(size(dcmtemp)); % for debugging - counts the no. of images in each slot
for j = 1:count
    temp(:,:,ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),kse(j),kri(j),knex(j)) = data(:,:,j);
    dcmtemp(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),kse(j),kri(j),knex(j)) = dcm(j);
    mask(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),kse(j),kri(j),knex(j)) = ...
        mask(ksl(j),kte(j),kfa(j),ktr(j),kti(j),kbv(j),kse(j),kri(j),knex(j))+1;
end
if nnz(mask~=1)
    disp('Something is wrong... duplicates or missing images - need to check')
    keyboard
end

% properties to return
data = temp;
dcm = dcmtemp;
params.te = uTE * 1e-3; % seconds
params.sl = uSL;
params.fa = uFA;
params.tr = uTR;
params.se = uSE;
params.ti = uTI;
params.B = uBV;
params.Tesla = uLF / 42.57747892;
params.nex= uNEX;
params.dcm = dcm{1};

clear temp dcmtemp mask

%% handle real/imag/mag/phase (not perfect)

% scale phase from int into float range (0 to 2pi)
if strfind(dcm{1}.Manufacturer,'SIEMENS')
    phase_scale = single(2*pi/4095);
elseif strfind(dcm{1}.Manufacturer,'GE')
    phase_scale = single(2*pi/1000);
else
    if any(uRI==1)
        warning('get-dicom() phase_scale not defined for phase images')
    end
    phase_scale = 1;
end

if isequal(uRI,1)
    % 1=phase
    data = data*phase_scale-pi;
elseif isequal(uRI,0) || nRI==1
    % 0,2,3=mag, real or imag (nothing to do)
elseif isequal(uRI,[0 1])
    % 0=mag and 1=phase
    data = data(:,:,:,:,:,:,:,:,:,1,:).*exp(i*(data(:,:,:,:,:,:,:,:,:,2,:)*phase_scale-pi));
elseif isequal(uRI,[2 3])
    % 2=real and 3=imag
    data = complex(data(:,:,:,:,:,:,:,:,:,1,:),data(:,:,:,:,:,:,:,:,:,2,:));
elseif isequal(uRI,[0 1 2 3])
    % ignore mag/phase (less reliable than real/imag)
    data = complex(data(:,:,:,:,:,:,:,:,:,3,:),data(:,:,:,:,:,:,:,:,:,4,:));
    disp([mfilename '(): real/imag/mag/phase present - returning real/imag'])
else
    disp([mfilename '(): not sure about real/imag/phase/mag - need to handle this'])
    keyboard
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% recursively fetch file info from a directory %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function info = get_files(pathname)

% contents of starting directory
info = dir(pathname);

j = 1;
while j <= numel(info)
    % recurse into subdirectories
    if info(j).isdir
        % skip subdirectories with a leading '.'
        if info(j).name(1) ~= '.'
            temp = dir([pathname info(j).name]);
            % prepend path (except for '.')
            for k = 1:numel(temp)
                if temp(k).name(1) ~= '.'
                    temp(k).name = [info(j).name filesep temp(k).name];
                end
            end
            % append contents of subdirectory
            info = [info;temp];
        end
        % delete directory from list
        info(j) = [];
    else
        % skip past files
        j = j + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% detect if file is DICOM (borrowed from old MATLAB version) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isdicom(filename)
%ISDICOM    Determine if a file is probably a DICOM file.
%    TF = ISDICOM(FILENAME) returns true if the file in FILENAME is
%    probably a DICOM file and FALSE if it is not.

%   Copyright 2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:10:36 $

% MB exclude DICOMDIR, which otherwise returns true
if strfind(filename,'DICOMDIR')
    tf = false;
    return;
end

% Open the file.
fid = fopen(filename, 'r');
if (fid < 0)
    %warning('Image:isdicom:fileOpen', ...
    %    'Couldn''t open file for reading: %s',filename)
    tf = false;
    return % MB just return, don't complain
end

% Get the possible DICOM header and inspect it for DICOM-like data.
header = fread(fid, 132, 'uint8=>uint8');
fclose(fid);

if numel(header)<132 % MB exclude easy pickings

    % It's too small
    tf = false;

elseif isequal(char(header(129:132))', 'DICM')

    % It's a proper DICOM file.
    tf = true;

else
    % MB heuristic is not good enough - exclude all.
    %
    % Use a heuristic approach, examining the first "attribute".  A
    % valid attribute will likely start with 0x0002 or 0x0008.
    %group = typecast(header(1:2), 'uint16');
    %if (isequal(group, uint16(2)) || isequal(swapbytes(group), uint16(2)) || ...
    %        isequal(group, uint16(8)) || isequal(swapbytes(group), uint16(8)))
    %
    %    tf = true;
    %else
        tf = false;
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         tricky counter for parfor loop           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flag = update_counter(tmpfile,maxno)
% Warning: "clever"
% 1 input = increment counter
% 2 inputs = display counter and test for maxno
if nargin==1

    % only proceed if file is not in use by another lab
    fid = fopen(tmpfile,'a');
    while fid==-1
        pause(rand*1e-2);
        fid = fopen(tmpfile,'a');
    end

    % append another 1 to counter file
    fwrite(fid,1);
    fclose(fid);

else

    % tmpfile is a vector (1 entry per image)
    fid = fopen(tmpfile,'r');
    count = numel(fread(fid));
    fclose(fid);

    % counter status
    flag = count<maxno;

    % display counter
    if flag
        % always use 4 characters
        fprintf('\b\b\b\b%-4d',count);
    end

end
