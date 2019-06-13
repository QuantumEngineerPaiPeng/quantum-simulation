%20180711
% compute two-point correlator or OTOC with trasverse field for different disorder strength
% and time
% Written using OperatorClass
% Created by Pai Peng
tic
N_atom=12;

bc='o';
% cplist=[1,2,3].^-3;
% cplist=1;
cplist=(1:(N_atom-1)).^-3;
model='dipolar';
% sym='kpz';
tlist=logspace(1,2,50);
% glist=linspace(0,10,20); % strength of the disorder field
glist=logspace(-1,1,21);
%glist=[0,0.3,0.6,0.8,0.9,0.95,1.05,1.1,1.2,1.4,1.7,2];

rho_0=OperatorClass(N_atom,'z',1);
% V0=OperatorClass(N_atom,'z',1);
% rho_0=LocalPauli(N_atom,1,'z');
% rho_0=LocalPauli(N_atom,1,'x')-1i*LocalPauli(N_atom,1,'y');
% V0=LocalPauli(N_atom,6,'x');
V0=(copy(rho_0));

% rho_0.symmetrize(sym);
% V0.symmetrize(sym);

C=zeros(length(glist),1);

H_int=OperatorClass(N_atom,model,-1,bc,cplist,1);
% H_int.symmetrize(sym);

Hz=OperatorClass(N_atom,'z',1);
% Hz.symmetrize(sym);

parfor p=1:length(glist)
    g=glist(p);
    H=H_int+g*Hz;
    H.diagonalize()
    %         C(:,col,pp)=TPC(rho_0,H,V0,tlist);
%     C(:,col)=OTOC(rho_0,H,V0,tlist);
%     averho=AveOpe(rho_0, H, tlist);
%     C(p)=TwoNorm(averho*V0-V0*averho);
    temp=0;
    for q=1:length(H.matrix)
        for pp=1:size(H.eigsys{q}.V,1)
            temp=temp+(H.eigsys{q}.V(:,pp)'*rho_0.matrix{q}*H.eigsys{q}.V(:,pp))^2;
        end
    end
    
%     C(p)=TwoNorm(averho);
    C(p)=temp;
end

C=C/(N_atom*2^N_atom);

figure(1)
hold on
plot(glist,C)
xlabel('g')
ppStyle(30,2,10)

toc
