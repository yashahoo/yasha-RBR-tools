function [Fserie] = OpSmooth(serie,minB,maxB,ratio)
[Lserie,Wserie] = size(serie);
if (maxB < minB)
  minB = maxB;
end;
Band=minB;
Fserie=zeros(Lserie,Wserie+1);
i=1;
j=1;
k=1;
while k <= Lserie-maxB + 1,
  if Band < maxB
    Band=round(minB*(ratio^(i-1)));
    i=i+1;
  else
    Band=maxB;
  end;
  kl=k-1+Band;
  if (Band > 1)
    Bmean=sum(serie(k:kl,:))./Band;
  else 
    Bmean=serie(k:kl,:);
  end;
  Fserie(j,:)= [Bmean Band];
  j = j + 1;
  k=k+Band;
  if (minB == 1) & (maxB ~= 1)
    minB = 2;
  end;
end;

Fserie=Fserie(1:j-1,:);

