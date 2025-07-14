%% Quantization of the input series Y with c quantization levels
% Y: input series, column data
% c: number of quantization levels

function x=bim_quantization(Y,c)

n=size(Y,1);

x=zeros(n,1);
ma=max(Y); mi=min(Y);
q=(ma-mi)/c; % amplitude of quantization level

l=zeros(c,1);
for i=1:c %quantization levels
   l(i)=mi+i*q;
end

for i=1:n
   j=1;
   while (Y(i)>=l(j))&&(j<c)
      j=j+1;
   end
   x(i)=j;
end

end

