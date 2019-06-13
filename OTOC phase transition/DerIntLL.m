% 20190320
% Get the derivative of $\kappa^{-1}$.
% Run LLanalysis first. $\kappa^{-1}$ is stored in variable 'gamma2'.
% Created by Pai Peng

% indlist=[2:2:19,20:length(wlist)];
indlist=[1:length(wlist)];
dgamma2=gamma2(indlist(2:end))-gamma2(indlist(1:end-1)); % Finite difference
errdgamma2=sqrt(errgamma2(indlist(2:end)).^2+errgamma2(indlist(1:end-1)).^2);
wlist2=(wlist(indlist(2:end))+wlist(indlist(1:end-1))).'/2; 
dwlist=(wlist(indlist(2:end))-wlist(indlist(1:end-1))).';
dgamma2dw=dgamma2./dwlist; % d(\kappa^{-1})/dw
errdgamma2dw=errdgamma2./dwlist;
figure
hold on
errorbar(wlist2.',dgamma2dw,errdgamma2dw)
xlabel('$W$','Interpreter','latex')
ylabel('$d{\kappa^{-1}}/dW$','Interpreter','latex')
ppStyle(20,2,10)