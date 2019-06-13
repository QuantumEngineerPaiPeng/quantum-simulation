% 20180803
% function to load NMR data into Matlab
% Input: 
% numList: list containing the expt numbers to be loaded
% fdname: folder name, optional, default 'expt1'
% win: optional. True for window detection
% Created by Pai Peng

function data_analysis(numList, fdname, win)

if nargin<3
    window=false;
else
    window=win;
end
if nargin<2
    fdname='expt1';
end
FID=[];D={};L=[];NS=[];P=[];PL={};PULP={};SFO1=[];SW={};DW=[];TD=[];

for p=numList
    
    try
        filename=['~/Dropbox (MIT)/grad/research/data/',fdname,'/',num2str(p),'/'];
        try
            fid=fopen([filename,'/ser'],'r','b');
            FID=[FID,fread(fid,'int32')];
            dimension=2; % for experiment dimension=2 or 3
        catch ME
            fid=fopen([filename,'/fid'],'r','b');
            FID=[FID,fread(fid,'int32')];
            dimension=1; % for experiment dimension=1
        end
    catch
        filename=sprintf('/Volumes/nmrsu1/expt11/%d',p);
        try
            fid=fopen([filename,'/ser'],'r','b');
            FID=[FID,fread(fid,'int32')];
            dimension=2; % for experiment dimension=2 or 3
        catch ME
            fid=fopen([filename,'/fid'],'r','b');
            FID=[FID,fread(fid,'int32')];
            dimension=1; % for experiment dimension=1
        end
    end
    
    fclose(fid);
    facqu=fopen([filename,'/acqu'],'r','b');
    A=[];
    TDtemp=[];
    while ~strcmp(A,'##$D= (0..63)')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    D{end+1}=str2num(fgetl(facqu));
    while ~strcmp(A,'##$L= (0..31)')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    L=[L;str2num(fgetl(facqu))];
    while ~contains(A,'##$NS=')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    NS=[NS;str2num(A(8:end))];
    while ~strcmp(A,'##$P= (0..63)')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    P=[P;str2num(fgetl(facqu))];
    while ~strcmp(A,'##$PL= (0..63)')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    PL{end+1}=str2num(fgetl(facqu));
    while ~contains(A,'##$PULPROG=')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    PULP{end+1}=A(13:end);
    while ~contains(A,'##$SFO1= ')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    SFO1=[SFO1;str2num(A(10:end))];
    while ~contains(A,'##$SW= ')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    SW{end+1}=(A(8:end));
    DW=[DW;(1e6/(2*str2num((A(8:end)))*str2num((A(10:end)))))];
    while ~contains(A,'##$TD=')
        A=fgetl(facqu);
        if feof(facqu)
            return
        end
    end
    TDtemp=str2num(A(8:end));
    fclose(facqu);
    facqu=fopen([filename,'/acqu2'],'r','b');
    if facqu>0
        A=fgetl(facqu);
        while ~contains(A,'##$TD=')
            A=fgetl(facqu);
            if feof(facqu)
                return
            end
        end
        TDtemp=[TDtemp,str2num(A(8:end))];
        fclose(facqu);
    end
    
    facqu=fopen([filename,'/acqu3'],'r','b');
    if facqu>0
        A=fgetl(facqu);
        while ~contains(A,'##$TD=')
            A=fgetl(facqu);
            if feof(facqu)
                return
            end
        end
        TDtemp=[TDtemp,str2num(A(8:end))];
        fclose(facqu);
    end
    if isempty(TD)
        TD=TDtemp;
    else
        if (TD-TDtemp)*(TD-TDtemp)'>1e-10
            error('Dimentions of the experiments do not match')
        end
    end
end

if ~window
    %     B=FID(1:TD(1):size(FID,1),:)+1i*FID(2:TD(1):size(FID,1),:);
    if length(TD)==3
        FID=reshape(FID,TD(1),TD(2),TD(3),length(numList));
        %         B=reshape(B,TD(2),TD(3),length(numList));
        B=squeeze(FID(1,:,:,:,:)+1i*FID(2,:,:,:,:));
        err=squeeze(std(FID(end-1-19*2:2:end-1,:,:,:,:),0,1));
    else
        if length(TD)==2
            FID=reshape(FID,TD(1),TD(2),length(numList));
            B=squeeze(FID(1,:,:,:)+1i*FID(2,:,:,:));
            err=squeeze(std(FID(end-1-19*2:2:end-1,:,:,:),0,1));
        end
    end
else
    B=FID(1:2:size(FID,1),:)+1i*FID(2:2:size(FID,1),:);
    if length(TD)==3
        B=reshape(B,TD(1)/2,TD(2),TD(3),length(numList));
    else
        if length(TD)==2
            B=reshape(B,TD(1)/2,TD(2),length(numList));
        end
    end
end
% dimension
if ~(length(TD)==1)
assignin('base','B',B);
assignin('base','err',err);
end
assignin('base','FID',FID);
assignin('base','D',D);
assignin('base','D',D);
assignin('base','L',L);
assignin('base','NS',NS);
assignin('base','P',P);
assignin('base','PL',PL);
assignin('base','PULP',PULP);
assignin('base','SFO1',SFO1);
assignin('base','SW',SW);
assignin('base','DW',DW);
assignin('base','TD',TD);