% 20180507
% write log file in a 'log.txt' file in fpath containing varname: varvalue
function writelog(fpath,varname,varvalue)
% varname={'a','b','c'};
% varvalue={1,'B',[1,2]};
if length(varname)~=length(varvalue)
    error('# of variable mismatch')
end
% fpath='~/Dropbox (MIT)/grad/research/codes/Time crystal/figure_data/';
fname=[fpath,'log.txt'];
fid=fopen(fname,'a');
fprintf(fid,'\n\n');
for p=1:length(varname)
    fprintf(fid,[varname{p},': ']);
    value=varvalue{p};
    switch class(value)
        case {'double','logical'}
            fprintf(fid,num2str(value));
        case 'char'
            fprintf(fid,value);
        otherwise
            error(['unidentified class of variable: ',varname{p}])
    end
    fprintf(fid,'\n');
end
fclose(fid);
end