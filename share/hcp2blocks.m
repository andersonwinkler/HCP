function varargout = hcp2blocks(restrfile,blocksfile,dz2sib,ids)
% Takes a "restricted data" CSV file from the HCP and generates
% a block file that can be used to make permutations in PALM.
%
% Usage:
% [EB,tabout] = hcp2blocks(restrfile,blocksfile,dz2sib,ids)
%
% Inputs:
% restrfile  : CSV file downloaded from https://db.humanconnectome.org/
%              containing the "Restricted Data"
%              For some releases, eg HCP500, some gawk processing is needed.
% blocksfile : CSV file to be created, with the exchangeability blocks,
%              ready for used with PALM.
% dz2sib     : (Optional) Defines whether dizygotic twins should be
%              treated as ordinary siblings (true), or be a category
%              on its own (false). Default = false.
% ids        : (Optional) A vector of subject IDs. If supplied, only the
%              subjects with the indicated IDs will be used.
%
% Outputs (if requested):
% EB      : Block definitions, ready for use, in the original order
%           as in the CSV file.
% tabout  : (Optional) A table containing:
%           - 1st col: subject ID
%           - 2nd col: mother ID
%           - 3rd col: father ID
%           - 4th col: sib type
%           - 5th col: family ID
%           - 6th col: family type
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / University of Oxford
% Dec/2013 (first version)
% Dec/2015 (this version)
% http://brainder.org

warning off backtrace

% Load the data and select what is now needed
tmp = strcsvread(restrfile);

% Locate the columns with the relevant pieces of info, i.e.,
% egoid, moid, faid, twin status and zygozity, in this exact
% order, and take just them
egid_idx = find(strcmpi(tmp(1,:),'Subject'));
moid_idx = find(strcmpi(tmp(1,:),'Mother ID')   | strcmpi(tmp(1,:),'Mother_ID'));
faid_idx = find(strcmpi(tmp(1,:),'Father ID')   | strcmpi(tmp(1,:),'Father_ID'));
twst_idx = find(strcmpi(tmp(1,:),'Twin Status') | strcmpi(tmp(1,:),'Twin_Stat'));
zygo_idx = find(strcmpi(tmp(1,:),'Zygosity'));
agey_idx = strcmpi(tmp(1,:),'Age_in_Yrs');
tab = tmp(2:end,[egid_idx moid_idx faid_idx twst_idx zygo_idx]);
age = cell2mat(tmp(2:end,agey_idx));

% If subjects have these elementary info missing, remove them
tab0a =  cellfun(@isnan, tab(:,1:3));
tab0b = ~cellfun(@ischar,tab(:,4:5));
tab0  = any(horzcat(tab0a,tab0b),2);
idstodel = cell2mat(tab(tab0,1));
warning([ ...
    'These subjects have data missing in the restricted file and will be removed: \n' ...
    repmat('         %d\n',1,numel(idstodel))],idstodel);
if nargin == 4 && ~ isempty(ids) && ~ isempty(idstodel),
    ids(any(bsxfun(@eq,ids(:),idstodel'),2)) = [];
end
tab(tab0,:) = [];
age(tab0)   = [];
N = size(tab,1);

% Treat non-monozygotic twins as ordinary siblings.
if nargin >= 3 && dz2sib,
    for n = 1:N,
        if strcmpi(tab{n,5},'notmz'),
            tab{n,4} = 'NotTwin';
            tab{n,5} = 'NotTwin';
        end
    end
end

% Instead of strings, use sib identifiers. These are based on the
% fact that it's unlikely any family on the HCP have more than 10
% siblings overall. If this happens, it'll be necessary to make
% changes here.
sibtype = zeros(N,1);
for n = 1:N,
    if strcmpi(tab{n,4},'nottwin');
        sibtype(n) = 10;
    elseif strcmpi(tab{n,5},'notmz');
        sibtype(n) = 100;
    elseif strcmpi(tab{n,5},'mz');
        sibtype(n) = 1000;
    end
end
tab = cell2mat(tab(:,1:3));

% Subselect subjects as needed
if nargin == 4 && ~isempty(ids) && islogical(ids(1)),
    tab        = tab(ids,:);
    sibtype    = sibtype(ids,:);
elseif nargin == 4 && ~ isempty(ids),
    idx = bsxfun(@eq,tab(:,1),ids');
    idx = ~ any(idx,1);
    warning([ ...
        'These subjects don''t exist in the restricted file and will be removed: \n' ...
        repmat('         %d\n',1,sum(idx))],ids(idx));
    ids(idx)   = [];
    tabnew     = zeros(length(ids),size(tab,2));
    sibtypenew = zeros(length(ids),1);
    agenew     = zeros(length(ids),1);
    for n = 1:length(ids),
        idx = tab(:,1) == ids(n);
        tabnew(n,:)     = tab(idx,:);
        sibtypenew(n,:) = sibtype(idx);
        agenew(n,:)     = age(idx);
    end
    tab        = tabnew;
    sibtype    = sibtypenew;
    age        = agenew;
end
N = size(tab,1);

% Create family IDs
famid = zeros(N,1);
U = unique(tab(:,2:3),'rows');
for u = 1:size(U,1),
    uidx = all(bsxfun(@eq,tab(:,2:3),U(u,:)),2);
    famid(uidx) = u;
end

% For parents that belong to more than one family, merge
% their families into just one, the one with lowest famid.
par = tab(:,2:3);
for p = par(:)', % for each parent
    pidx = any(par == p,2);
    famids = unique(famid(pidx)); % families that he/she belong to
    for f = 1:numel(famids),
        famid(famid == famids(f)) = famids(1);
    end
end

% Label each family according to their type. The "type" is
% determined by the number and type of siblings.
F = unique(famid);
famtype = zeros(N,1);
for f = 1:numel(F);
    fidx = F(f) == famid;
    famtype(fidx) = sum(sibtype(fidx)) + numel(unique(tab(fidx,2:3)));
end

% Twins which pair data isn't available should be treated as
% non-twins, so fix and repeat computing the family types
idx = (sibtype == 100  & (famtype >= 100  & famtype <= 199)) ...
    | (sibtype == 1000 & (famtype >= 1000 & famtype <= 1999));
sibtype(idx) = 10;
for f = 1:numel(F);
    fidx = F(f) == famid;
    famtype(fidx) = sum(sibtype(fidx)) + numel(unique(tab(fidx,2:3)));
end

% Append the new info to the table.
tab = horzcat(tab,sibtype,famid,famtype);
if nargout == 2,
    varargout{2} = [tab age];
end

% Families of the same type can be shuffled, as well as sibs of the same
% type. To do this, the simplest is to construct the blocks within each
% family type, then replicate across the families of the same type.
% Start by sorting
[~,idx] = sortrows([famid sibtype age]);
[~,idxback] = sort(idx);
tab     = tab(idx,:);
sibtype = sibtype(idx);
famid   = famid(idx);
famtype = famtype(idx);
age     = age(idx,:);

% Now make the blocks for each family
B = cell(numel(F),1);
for f = 1:numel(F),
    fidx = F(f) == famid;
    ft = famtype(find(fidx,1));
    if any(ft == [(12:10:92) 23 202 224 2002]),
        B{f} = horzcat(famid(fidx),sibtype(fidx),tab(fidx,1));
    else
        B{f} = horzcat(-famid(fidx),sibtype(fidx),tab(fidx,1));
        
        % Some particular cases of complicated families
        if ft == 33,
            tabx = tab(fidx,2:3);
            for s = 1:size(tabx,1),
                if (sum(tabx(:,1) == tabx(s,1)) == 2 && ...
                        sum(tabx(:,2) == tabx(s,2)) == 3) || ...
                        (sum(tabx(:,1) == tabx(s,1)) == 3 && ...
                        sum(tabx(:,2) == tabx(s,2)) == 2),
                    B{f}(s,2) = B{f}(s,2) + 1;
                end
            end
        elseif ft == 34,
            tabx = tab(fidx,2:3);
            for s = 1:size(tabx,1),
                if sum(tabx(:,1) == tabx(s,1)) == 2 ...
                        && sum(tabx(:,2) == tabx(s,2)) == 2,
                    B{f}(s,2) = B{f}(s,2) + 1;
                end
            end
        elseif ft == 43,
            tabx = tab(fidx,2:3);
            for s = 1:size(tabx,1),
                if (sum(tabx(:,1) == tabx(s,1)) == 3 && ...
                        sum(tabx(:,2) == tabx(s,2)) == 4) || ...
                        (sum(tabx(:,1) == tabx(s,1)) == 4 && ...
                        sum(tabx(:,2) == tabx(s,2)) == 3),
                    B{f}(s,2) = B{f}(s,2) + 1;
                end
            end
        elseif ft == 223,
            sibx = sibtype(fidx);
            B{f}(sibx == 10,2) = -B{f}(sibx == 10,2);
        elseif ft == 302,
            famtype(fidx) = 212;
            tmpage = age(fidx);
            if tmpage(1) == tmpage(2),
                B{f}(3,2) = 10;
            elseif tmpage(1) == tmpage(3),
                B{f}(2,2) = 10;
            elseif tmpage(2) == tmpage(3),
                B{f}(1,2) = 10;
            end
        end
    end
end

% Concatenate all. Prepending the famtype ensures that the
% families of the same type can be shuffled whole-block. Also,
% add column with -1, for within-block at the outermost level
B = horzcat(-ones(N,1),famtype,cell2mat(B));

% Sort back to the original order
B = B(idxback,:);

% Drop columns that are redundant (useful when the supplied ids
% contain just a few subjects)
for c = size(B,2):-1:2,
    if numel(unique(B(:,c))) == 1,
        B(:,c) = [];
    end
end
if nargout > 0,
    varargout{1} = B;
end

% Save as CSV
if nargin >= 2 && ~isempty(blocksfile) && ischar(blocksfile),
    dlmwrite(blocksfile,B,'precision','%d');
end

warning on backtrace
