% 20180827
% Constructing l-bits using Ken's method
% Written using OperatorClass
% Created by Pai Peng
tic
N_atom=8;
N_rand=100;
bc='o';
% cplist=[1,2,3].^-3;
cplist=1;
model=[1,1,1];
DisorderDir=[0,0,1]; %direction of the disorder field
DisorderDir=DisorderDir/norm(DisorderDir);

w=20; % strength of the disorder field

H_int=OperatorClass(N_atom,model,1,bc,cplist,1);

LD=cell(N_atom,N_atom,N_rand);
Dtau=zeros(1,N_rand);

Hstrength=zeros(N_atom-1,N_rand);

rng(0)
randmat=rand(N_rand,N_atom);

parfor pp=1:N_rand
    DisorderAmp=(randmat(pp,:)*2-1);
    Disorder=DisorderDir.'*DisorderAmp; % field uniformly distributed in [-w,w]
    H_rand=randOpe(N_atom,Disorder);
    %     H_rand=OperatorClass(N_atom,'x',1);
    
    H=H_int+w*H_rand;
    H.diagonalize()
    tau=cell(1,N_atom);
    tauD=zeros(2^N_atom,N_atom);
    LDtemp=cell(N_atom,N_atom);
    for p=1:N_atom
        sp=LocalPauli(N_atom,p,'z');
        olap=zeros(1,2^N_atom);
        for q=1:2^N_atom
            olap(q)=H.eigsys{1}.V(:,q)'*sp.matrix{1}*H.eigsys{1}.V(:,q);
        end
        [olap,ind]=sort(olap);
        temp=ones(1,2^N_atom);
        temp(ind(1:end/2))=-1;
        %-sum(olap(1:end/2))+sum(olap((end/2+1):end))
        tau{p}=H.eigsys{1}.V*sparse(diag(temp))*H.eigsys{1}.V';
        amp=DecomPauli(tau{p});
        LDc=LDcell(N_atom,4,p);
        for q=1:N_atom
            LDtemp{p,q}=amp(LDc{q});
        end
        tauD(:,p)=temp;
    end
    LD(:,:,pp)=LDtemp;
%     temp1=GetOverlapMat(tau);
%     Dtau(pp)=sum(sum(abs(temp1-diag(diag(temp1))))); % sum of abs of overlaps
%     temp2=zeros(N_atom-1,1);
%     count2=temp2;
%     for p=1:N_atom-1
%         for q=p+1:N_atom
%             temp2(q-p)=temp2(q-p)+abs(sum(tauD(:,p).*tauD(:,q).*H.eigsys{1}.D));
%             count2(q-p)=count2(q-p)+1;
%         end
%     end
%     Hstrength(:,pp)=temp2./count2;
end
Hstrength_ave=mean(Hstrength,2);

if N_rand==1
    figure
    plot(0:(N_atom-1),LD_ave)
else
    LDmed=zeros(N_atom,N_atom);
    figure(1)
    for p=1:N_atom
        for q=1:max((N_atom+1-p),p)
            temp=abs(reshape(cell2mat(LD(p,q,:)),1,[]));
            temp(temp<1e-13)=[];
            LDmed(p,q)=median(temp);
        end
        plot(0:N_atom-1,LDmed(p,:),'b')
        set(gca,'yscale','log')
        hold on
    end
    xlabel('Distance')
    ylabel('Weight of \tau_z')
    ppStyle(30,2,10)
end
toc