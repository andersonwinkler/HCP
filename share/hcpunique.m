function hcpunique(restrfile,listfile)
% Takes a "restricted" CSV file and produce a random list of "unique"
% subjects, i.e., subjects that are not known to be related one to another.
%
% Usage:
% hcp2solar(restrfile,listfile)
%
% restrfile   : CSV file downloaded from https://db.humanconnectome.org/
%               containing the "Restricted Data"
%               For some releases, eg HCP500, some gawk parsing is needed.
% listfile    : File to be created containing the IDs.
%
% _____________________________________
% Anderson M. Winkler
% UTRGV
% Aug/2022
% http://brainder.org

% Load the restricted data
tmp = strcsvread(restrfile);

% Locate the columns with the relevant pieces of info
egid_idx = find(strcmpi(tmp(1,:),'Subject'));
moid_idx = find(strcmpi(tmp(1,:),'Mother ID') | strcmpi(tmp(1,:),'Mother_ID'));
faid_idx = find(strcmpi(tmp(1,:),'Father ID') | strcmpi(tmp(1,:),'Father_ID'));
tabr = cell2mat(tmp(2:end,[egid_idx moid_idx faid_idx]));
rnd  = rand(size(tabr,1),1);
tabr = horzcat(tabr,rnd);

% For every mother, get one offspring
Umo = unique(tabr(:,2));
selected = zeros(size(Umo,1),4);
m = 1;
for mo = Umo'
    idx = tabr(:,2) == mo;
    tmp = tabr(idx,:);
    idx = tmp(:,4) == max(tmp(:,4));
    selected(m,:) = tmp(idx,:);
    m = m + 1;
end

% Repeat for the fathers
tabr = selected;
Ufa  = unique(tabr(:,3));
selected = zeros(size(Ufa,1),4);
f = 1;
for fa = Ufa'
    idx = tabr(:,3) == fa;
    tmp = tabr(idx,:);
    idx = tmp(:,4) == max(tmp(:,4));
    selected(f,:) = tmp(idx,:);
    f = f + 1;
end
selected = sort(selected(:,1));

% Now save!
fid = fopen(listfile,'w');
fprintf(fid,'%d\n',selected');
fclose(fid);
