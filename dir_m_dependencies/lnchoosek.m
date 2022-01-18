function output = lnchoosek(n,k);
% test with: ;
%{
  n=21; kra=0:n;
  ln_1=zeros(length(kra),1);
  ln_2=zeros(length(kra),1);
  for nk=1:length(kra);
  k=kra(nk);
  ln_1(nk) = log(nchoosek(n,k));
  ln_2(nk) = lnchoosek(n,k);
  end;%  for nk=1:length(kra);
  plot(kra,ln_1,'rx-',kra,ln_2,'bo');
  %or ;
  n=25; kra=0:1:25;
  plot(kra,lnchoosek(n,kra),'.-');  
  %}
if (n<15);
output = nfact(n) - nfact(n-k) - nfact(k);
output(find(~isfinite(output)))=-16;
 else;
output = sfact(n) - sfact(n-k) - sfact(k);
output(find(~isfinite(output)))=-16;
end;%if (n<15);

function output = nfact(d);
dij=find(d>0);
output(dij) = log(factorial(d(dij)));
output(find(d<=0)) = 0;

function output = sfact(d);
dij=find(d>0);
output(dij) = d(dij).*log(d(dij)) - d(dij) + 0.5.*log(2*pi*d(dij)) + 1/12./d(dij) - 1/360./d(dij).^3 + 1/1260./d(dij).^5 - 1/1680./d(dij).^7 ;
output(find(d<=0)) = 0;
