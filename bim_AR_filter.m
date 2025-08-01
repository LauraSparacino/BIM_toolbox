%% Autoregressive low-pass, high-pass filter

%%% INPUT
% sig: input data matrix N*2 where 2 is the number of processes
% can: channel to be filtered among the 2 channels
% p: order of the AR filter (e.g., p=0.94)

%%% OUTPUT
% fia: high-pass filtered signal
% fib: low-pass filtered signal

function [fia,fib]=bim_AR_filter(Y,can,p)

s=Y(:,can);
s1=size(s,1);

usc = zeros(s1,1);
fib = zeros(s1,1);

usc(1)= s(1);
for i = 2 : s1
   usc(i)= p*usc(i-1)+(1-p)*s(i);
end 

fib(s1)=usc(s1);
for i = (s1-1):-1:1
     fib(i)= p*fib(i+1)+(1-p)*usc(i);
end
   
fia = s-fib+mean(s);
