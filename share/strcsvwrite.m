function strcsvwrite(varargin)
% Write a CSV file containing non-numeric fields (strings) from a
% cell array or array.
%
% Usage:
% strcsvwrite(table,filename,delimiter,end-of-line,numformat)
%
% - table    = Table to be saved (array or cell array)
% - filename = CSV file to be created
% - delim    = Field separator.  Default = ','
% - eol      = Record separator. Default = '\n'
% - numfmt   = Number format.    Default = '%g'
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Oct/2010
% http://brainder.org

% Accept the inputs
error(nargchk(1,5,nargin));
table    = varargin{1};
filename = varargin{2};
if nargin == 5,
    delimiter = varargin{2};
    eolmark   = varargin{3};
    numformat = varargin{4};
elseif nargin == 4,
    delimiter = varargin{2};
    eolmark   = varargin{3};
    numformat = '%g';
elseif nargin == 3,
    delimiter = varargin{2};
    eolmark   = '\n';
    numformat = '%g';
elseif nargin == 2,
    delimiter = ',';
    eolmark   = '\n';
    numformat = '%g';
end

% Write the file
[nR,nC] = size(table);
fid = fopen(filename,'w');
if iscell(table),
    for r = 1:nR,
        for c = 1:nC,
            if ischar(table{r,c}),
                fprintf(fid,'%s',table{r,c});
            else
                fprintf(fid,numformat,table{r,c});
            end
            if c < nC,
                fprintf(fid,delimiter);
            else
                fprintf(fid,eolmark);
            end
        end
    end
else
    for r = 1:nR,
        for c = 1:nC,
            fprintf(fid,numformat,table(r,c));
            if c < nC,
                fprintf(fid,delimiter);
            else
                fprintf(fid,eolmark);
            end
        end
    end
end
fclose(fid);
