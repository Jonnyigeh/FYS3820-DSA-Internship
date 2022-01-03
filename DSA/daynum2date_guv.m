function [aar,maaned,dag]=daynum2date_guv(yy,dagnummer)
% Beregner dag måned ut fra dagnummer gitt i SCAN.Dagnummer
%
% returnerer vektorer av dato

dagerimnd=[31 28 31 30 31 30 31 31 30 31 30 31];
if (rem(yy,4)==0) & (rem(yy,400)~=0)
  dagerimnd(2)=29;
end

if yy==2000
	dagerimnd(2)=29;
end


cumdays=cumsum(dagerimnd);


if dagnummer <=cumdays(1)        % <=
    a=0;
    maaned=1;
    dag=dagnummer;
    aar=yy; 
else
    a=find(cumdays<dagnummer);    %siste element inneholder maaned-1   % <=
    maaned=(length(a)+1);
    dag=(dagnummer-cumdays(length(a)) );
    aar=yy;
end