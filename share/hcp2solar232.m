function hcp2solar(restrfile,unrestrfile,pedfile)
% Takes a "restricted" and an "unrestricted" CSV files from the HCP
% and generates a pedigree file that can be used in SOLAR.
% 
% Usage:
% hcp2solar(restrfile,unrestrfile,pedfile)
% 
% restrfile   : CSV file downloaded from https://db.humanconnectome.org/
%               containing the "Restricted Data"
% unrestrfile : CSV file downloaded from https://db.humanconnectome.org/
%               containing the "Unrestricted Data"
% pedfile     : File name for the pedigree file to be created.
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
tmpr = strcsvread(restrfile);
tabr = tmpr(2:end,1:5);
idsr = cell2mat(tabr(:,1:3));
idsr = sortrows(idsr);

% Load the unrestricted data
tmpu = strcsvread(unrestrfile);
tabu = tmpu(2:end,1:2);
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

% Create FAMID. All subjects with a relative in common get the same FAMID,
% which initially is the lowest ID of whatever parent.
famid = zeros(N,1);
for s = 1:N,
    mimo = min(idsr(idsr(s,2) == idsr(:,2),3));
    mifa = min(idsr(idsr(s,3) == idsr(:,3),2));
    famid(s) = min(mimo,mifa);
end
[~,~,famid] = unique(famid);

% Founders. Note that the sex of the father here is 2, but when SOLAR runs,
% it reassings it to 1 (and vice versa for the mother).
mo = idsr(:,2);
mo = horzcat(mo,z,z,o,z,famid);   % mother: sex = 1
mo = unique(mo,'rows');
fa = idsr(:,3);
fa = horzcat(fa,z,z,2*o,z,famid); % father: sex = 2
fa = unique(fa,'rows');

% MZ subjects
mzidx = strcmpi(tabr(:,5),'MZ');
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
su = horzcat(idsr,sx,mz,famid);

% Pedigree, ready to save
ped = vertcat(mo,fa,su);

% Now save!
fid = fopen(pedfile,'w');
fprintf(fid,'id,fa,mo,sex,mztwin,famid\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d\n',ped');
fclose(fid);
