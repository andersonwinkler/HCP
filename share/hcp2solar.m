function hcp2solar(restrfile,unrestrfile,pedfile,hhoption)
% Takes a "restricted" and an "unrestricted" CSV files from the HCP
% and generates a pedigree file that can be used in SOLAR.
% 
% Usage:
% hcp2solar(restrfile,unrestrfile,pedfile,hhoption)
% 
% restrfile   : CSV file downloaded from https://db.humanconnectome.org/
%               containing the "Restricted Data"
%               For some releases, eg HCP500, some gawk parsing is needed.
% unrestrfile : CSV file downloaded from https://db.humanconnectome.org/
%               containing the "Unrestricted Data"
% pedfile     : File name for the pedigree file to be created.
% hhoption    : (Optional) Specify whether a HH should be included in the
%               pedigree file. It can be:
%               - 'none'   : No household (default)
%               - 'famid'  : HH is for each family, i.e., subjects with at
%                            least one common relative.
%               - 'mother' : HH is for each mother.
%               This option affects only half-sibs that share a common
%               father.
% 
% The pedigree file generated (a CSV file) can be loaded into SOLAR
% with "ped load pedfile.csv", and it will automatically generate
% the pedindex.cde, pedindex.out, and phi2.gz.
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Mar/2014
% http://brainder.org

% Load the restricted data
tmp = strcsvread(restrfile);

% If there is no Zygosity field, create it from ZygSR and ZygGT
zygo_idx = find(strcmpi(tmp(1,:),'Zygosity'));
if isempty(zygo_idx),
    zygoSR_idx = find(strcmpi(tmp(1,:),'ZygositySR'));
    zygoGT_idx = find(strcmpi(tmp(1,:),'ZygosityGT'));
    tmp(:,end+1) = cell(size(tmp,1),1);
    tmp{1,end} = 'Zygosity';
    for s = 2:size(tmp,1),
        if  (numel(tmp{s,zygoGT_idx}) == 1 && isnan(tmp{s,zygoGT_idx})) || ...
            (ischar(tmp{s,zygoGT_idx}) && strcmpi(tmp{s,zygoGT_idx},' ')) || ...
            isempty(tmp{s,zygoGT_idx}),
            tmp{s,end} = tmp{s,zygoSR_idx};
        else
            tmp{s,end} = tmp{s,zygoGT_idx};
        end
    end
end

% Locate the columns with the relevant pieces of info
egid_idx = find(strcmpi(tmp(1,:),'Subject'));
moid_idx = find(strcmpi(tmp(1,:),'Mother ID') | strcmpi(tmp(1,:),'Mother_ID'));
faid_idx = find(strcmpi(tmp(1,:),'Father ID') | strcmpi(tmp(1,:),'Father_ID'));
zygo_idx = find(strcmpi(tmp(1,:),'Zygosity'));
tabr = tmp(2:end,[egid_idx moid_idx faid_idx zygo_idx]);

% Load the unrestricted data
tmpu = strcsvread(unrestrfile);
egid_idx = find(strcmpi(tmpu(1,:),'Subject'));
sex_idx  = find(strcmpi(tmpu(1,:),'Sex') | strcmpi(tmpu(1,:),'Gender'));
tabu = tmpu(2:end,[egid_idx sex_idx]);

% If subjects have elementary info missing, remove them
todel = cellfun(@isempty,tabr) | cellfun(@chknan,tabr);
todel = any(todel,2);
tabr(todel,:) = [];
tabu(todel,:) = [];

% Convert to ordinary array
idsr = cell2mat(tabr(:,1:3));
idsr = sortrows(idsr);
idsu = cell2mat(tabu(:,1));
idsu = sortrows(idsu);

% Remove the ids that don't have their counterparts in the other file
cmp = bsxfun(@eq,idsr(:,1),idsu');
idsr(~any(cmp,2),:) = [];
tabr(~any(cmp,2),:) = [];
tabu(~any(cmp,1))   = [];

% Vars for later
N = size(idsr,1);
z = zeros(N,1);
o = ones(N,1);

% Create HHID.
if nargin < 4,
    hhoption = [];
end
if strcmpi(hhoption,'famid'),
    % All subjects with a relative in common get the same FAMID,
    % which initially is the lowest ID of whatever parent.
    % This FAMID is used as a proxy for HHID.
    hhid = zeros(N,1);
    for s = 1:N,
        mimo = min(idsr(idsr(s,2) == idsr(:,2),3));
        mifa = min(idsr(idsr(s,3) == idsr(:,3),2));
        hhid(s) = min(mimo,mifa);
    end
    [~,~,hhid] = unique(hhid);
    pedheader = 'id,fa,mo,sex,mztwin,hhid\n';
    printfmt = '%d,%d,%d,%d,%d,%d\n';
elseif strcmpi(hhoption,'mother'),
    % Subjects receive an HHID that is the same as their mother.
    % The HHID for the fathers the same as its spouse if just one
    mo = idsr(:,2);
    [~,~,hhid] = unique(mo);
    pedheader = 'id,fa,mo,sex,mztwin,hhid\n';
    printfmt = '%d,%d,%d,%d,%d,%d\n';
else
    % No HHID is created at all.
    hhid = [];
    pedheader = 'id,fa,mo,sex,mztwin\n';
    printfmt = '%d,%d,%d,%d,%d\n';
end

% Founders. Note that the sex of the father here is 2, but when SOLAR runs,
% it reassings it to 1 (and vice versa for the mother).
mo = idsr(:,2);
mo = horzcat(mo,z,z,o,z,hhid);   % mother: sex = 1
mo = unique(mo,'rows');
fa = idsr(:,3);
fa = horzcat(fa,z,z,2*o,z,hhid); % father: sex = 2
fa = unique(fa,'rows');

% Remove duplicated fathers (more than one spouse). Arbitrarily, keep
% those with the lowest mother ID (i.e., lowers HHID). This arbitrariness
% makes no difference for the HCP data because the founders have no data
% anyway (the are in the pedigree for complenetess). 
if strcmpi(hhoption,'mother'),
    farep = unique(fa(sum(bsxfun(@eq,fa(:,1),fa(:,1)'),2) > 1,1));
    for f = 1:numel(farep),
        rowrep = find(fa(:,1) == farep(f));
        fa(rowrep(2:end),:) = [];
    end
end

% MZ subjects
mzidx = strcmpi(tabr(:,4),'MZ');
mfm = horzcat(idsr(:,2:3),mzidx);
mz = z;
for m = find(mzidx)',
    mz(all(bsxfun(@eq,mfm,mfm(m,:)),2)) = m;
end

% MZ subjects that have their pair missing are treated as non-twin
U = unique(mz(mz > 0));
for u = 1:numel(U),
    uidx = mz == U(u);
    if sum(uidx) == 1,
        mz(uidx) = 0;
    end
end

% Assign the sex of the subjects
sx = z;
sx(strcmpi(tabu(:,2),'F')) = 1;
sx(strcmpi(tabu(:,2),'M')) = 2;
su = horzcat(idsr,sx,mz,hhid);

% Pedigree, ready to save
ped = vertcat(mo,fa,su);

% Now save!
fid = fopen(pedfile,'w');
fprintf(fid,pedheader);
fprintf(fid,printfmt,ped');
fclose(fid);

% ------------------------------------------------
function y = chknan(x)
% Depending on the version of strcsvread,
% there may be many NaNs
y = any(isnan(x(:)));