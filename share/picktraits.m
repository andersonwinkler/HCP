function trt = picktraits(csvfiles,trtlist,idlist,nohdr,saveas)
% Select a list of traits from the HCP data
%
% Usage:
% trt = picktraits(csvfiles,trtlist,idlist,nohdr,saveas)
%
% - csvfiles : List of restricted, unrestricted, and perhaps other
%              HCP files containing the traits. If more than one,
%              use input as a cell, e.g.:
%              {'restricted.csv','unrestricted.csv'}.
%              Leave empty, e.g. {}, to use the default files (their
%              locations can be defined inside).
% - trtlist  : List of traits to be selected, e.g.:
%              {'CardSort_AgeAdj','Flanker_AgeAdj','ListSort_AgeAdj'}
% - idlist   : Optional. List of subject IDs to be selected, e.g.:
%              [100307 100408 101107 101309 101410 102008 102311]'
% - nohdr    : Optional. A flag to indicate whether the headers for
%              rows and columns should be removed.
%              Default is false, i.e., headers are included.
% - saveas   : File name to save the table as.
%
% Note that the trait names ARE case-sensitive.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Aug/2014
% http://brainder.org

% Take input arguments
if nargin < 3 || nargin > 5,
    error('Incorrect number of input arguments');
elseif nargin == 3,
    nohdr = false;
    saveas = '';
end
%csvfiles = {'restricted_zygmerged.csv','unrestricted.csv'};
if isempty(csvfiles),
    
    % Default locations (change here as needed)
    csvfiles = {...
        '~/Dropbox/projects-dphil/hcp/hcp900/restricted.csv',...
        '~/Dropbox/projects-dphil/hcp/hcp900/unrestricted.csv'};
elseif ischar(csvfiles),
    csvfiles = {csvfiles};
end
if ischar(trtlist),
    trtlist = {trtlist};
end

% Load the files
nC = numel(csvfiles);
tab = cell(nC,1);
for c = 1:nC,
    tab{c} = strcsvread(csvfiles{c});
end

% Remove repeated fields and merge tables
for c1 = 1:nC,
    H1 = tab{c1}(1,:);
    for c2 = c1+1:nC,
        H2 = tab{c2}(1,:);
        for h1 = 1:numel(H1),
            idxrep = strcmp(H1{h1},H2); % change this line to strcmpi to ignore case differences
            H2(idxrep) = [];
            tab{c2}(:,idxrep) = [];
        end
    end
end
tab = horzcat(tab{:});

% Pick the IDs,
if ~ isempty(idlist)
    idcol = strcmp('Subject',tab(1,:));
    ids   = cell2mat(tab(2:end,idcol));
    ididx = zeros(size(idlist));
    for i = 1:numel(idlist),
        ididx(i) = find(ids == idlist(i));
    end
    %idsel = find(any(bsxfun(@eq,cell2mat(tab(2:end,idcol)),idlist(:)'),2));
    tab   = tab([1; ididx+1],:);
end

% Now pick the traits
trtlist = horzcat({'Subject'},trtlist{:});
idx = zeros(size(trtlist));
for t = 1:numel(trtlist),
    f = find(strcmp(trtlist{t},tab(1,:)));
    if isempty(f),
        fprintf('Trait "%s" not found\n',trtlist{t});
    else
        idx(t) = f;
    end
end
idx(~idx) = [];
trt = tab(:,idx);

% Remove the headers if asked
if nohdr,
    trt = cell2mat(trt(2:end,2:end));
end

% Save as a table:
if ~ isempty(saveas),
    strcsvwrite(trt,saveas);
end