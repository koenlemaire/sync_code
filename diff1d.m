function [ dxdt ] = diff1d( x,fs)
%[ dxdt ] = diff1d( x,dt)
%   differentiates N-Dimensional signal using central difference. Only
%   signals with equidistant time sampling are suited for this function.
%   Note: the first and last sample are unreliable estimates of the
%   derivative. 
% input: x [n_channel x n_sample] signal, fs sample frequency.
% output: dxdt 
%
% example:
% fs=1000;
% t=0:(1/fs):1;
% x=[sin(2*pi*t); t.^2+t-1];
% dxdt=diff1d(x,fs);
% figure;subplot(211); plot(t,x)
% subplot(212);plot(t,dxdt)

% koen lemaire 9/2015

dxdt=(x(:,3:end)-x(:,1:end-2))*.5*fs;% derivative
dxdt=[dxdt(:,1) dxdt dxdt(:,end)]; % correct dimension
end

