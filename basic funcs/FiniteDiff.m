% 20190325
% get dy/dx using finite difference
% Input:
% x: variable (Chinese pinyin: zi bian liang)
% y: dependent (Chinese pinyin: yin bian liang)
% err: error bar of y
% x,y,err should have the same size
% Output:
% xbar: x values where the derivative is evaluated
% dydx: dy/dx
% errdydx: error bar of dy/dx
% Created by Pai Peng
function [xbar, dydx,errdydx]=FiniteDiff(x,y,err)
dy=y(2:end)-y(1:(end-1)); % Finite difference
errdy=sqrt(err(2:end).^2+err(1:(end-1)).^2);
xbar=(x(2:end)+x(1:(end-1)))/2; 
dx=(x(2:end)-x(1:(end-1)));
dydx=dy./dx; % d(\kappa^{-1})/dw
errdydx=errdy./dx;