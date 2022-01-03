function dagnr=dagnr(aar,mnd,dag)
%	Function dagnr=dagnr(aar,mnd,dag)
% Finner dagnummer i året, 1.jan er dagnr 1 ;
%
% se også:
%function [aar,maaned,dag]=daynum2date(yy,SCAN) og [aar,maaned,dag]=daynum2date_multidays(yy,dagnummer)
% som beregner dag og måned ut fra dagnummer

%dagnr = 000;

dagerimnd=[31 28 31 30 31 30 31 31 30 31 30 31];
if (rem(aar,4)==0) & (rem(aar,400)~=0)
  dagerimnd(2)=29;
end

if aar==2000
	dagerimnd(2)=29;
end

if mnd>1
 dagnr=dag+sum(dagerimnd(1:(mnd-1)));
else
 dagnr=dag;
end

