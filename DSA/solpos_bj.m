function [azi,zen]=solpos_bj(year,month,day,time,minutt,sekund,breddegrad,lengdegrad)
%Purpose:
%	Calculate the azimuth and zenith angle of the sun	
%Syntax:
%	[azi,zen]=solpos(year,month,date,time,min,sek); caclulates the azimuth and zenith angle 
%	of the sun on year mont date at the given time time:min:sek
%Description:
%	The algoritm gives the solar zenith angle in radians relative to zentih. This means that
%	at sunset the sentih angle is pi/2. The azimuth angle is relative to south with a negative
%	azimuth angle to the west of south.
%Example:
%	To calculate the position of the sun at june 21. 1997, 12:00;
%	[a s]=solpos(1997,6,21,12,0,0);
%	a=a*180/pi;s=s*180/pi; % converts the angle to degrees.	
%Remarks: 
%	Calculates the solar position in radians, not degrees.
%	Assumes that the time given is local winter time (sun time!).
%Algorithm:
%
%Diagnostics:
%
%See also:
%	zenang, ol752_addsunangle, guv_addsunangle
%References:
%	Muhammad Iqbal, An introduction to solar radiation. Academic Press 1983
%Author: 
%	?, ?, 1995, Oddbjørn Grandum 
%Files needed: 
%
%Group:
% 	| general |	 	
% NB! (Bjørn des. 2000) hour må adderes 1, antakelig er det glemt å ta hensyn til tidssonen
%

%tar inn array av tidsverdier

N=length(time);
bredde=breddegrad*pi/180;
lengde=lengdegrad*pi/180;

sekundnr=time*3600+minutt*60+sekund;
dagogtime=dagnr(year,month,day) + (sekundnr/(60*60*24));
dagvinkel=2*pi*(dagogtime-1)/365;

tidsjevn=round( (0.000075+0.001868*cos(dagvinkel)...
    -0.032077*sin(dagvinkel)...
    -0.014615*cos(2*dagvinkel)...
    -0.04089*sin(2*dagvinkel))*229.18*60 );

%  lokallengdekorr har MINUSFORTEGN foran seg (i forhold til "Iqbal")
%  fordi lengdegraden vår (Norge) er ØST for Greenwich!!!
lokallengdekorr= -4*(15-lengdegrad)*60;	%Skal være lokallengdekorr-Standardtid= lign. til hoyre
lokalsoltid=round(sekundnr+lokallengdekorr+tidsjevn);

timevinkel= -( (lokalsoltid-60*60*12)/(60*60*24) )*2*pi ;
deklin=   0.006918-0.399912*cos(dagvinkel)...
    +0.070257*sin(dagvinkel)...
    -0.006758*cos(2*dagvinkel)...
    +0.000907*sin(2*dagvinkel)  ;

cos_zen = sin(deklin)*sin(bredde) + cos(deklin)*cos(bredde).*cos(timevinkel);
sin_zen = sqrt(1.0-(cos_zen.*cos_zen));
tan_zen = sin_zen./cos_zen ;
zen = atan(tan_zen) ;

p=find(tan_zen<0.0);
if ~isempty(p)
    zen(p) = zen(p) + pi ;
end

% if tan_zen<0.0
%     zen = zen + pi ;
% end
cos_azi = (cos(zen)*sin(bredde)-sin(deklin))./(sin(zen)*cos(bredde));
sin_azi = sqrt(1.0-(cos_azi.*cos_azi));

q1=find(timevinkel>=0.0);
if ~isempty(q1)
    sin_azi(q1) =  abs(sin_azi(q1));
end

q2=find(timevinkel<0.0);
if ~isempty(q2)
    sin_azi(q2) =  -abs(sin_azi(q2));
end

% if timevinkel>=0.0
%     sin_azi =  abs(sin_azi);
% else
%     sin_azi = -abs(sin_azi);
% end
tan_azi = sin_azi./cos_azi ;
azi = atan(tan_azi) ;

r=find(timevinkel<0.0 & tan_azi<0.0 & azi>0.0);
if ~isempty(r)
    azi(r) = azi(r) - pi ;
end

s=find(timevinkel<0.0 & tan_azi>0.0 & azi>0.0);
if ~isempty(s)
    azi(s) = azi(s) - pi ;
end

t=find(timevinkel>0.0 & tan_azi<0.0 & azi<0.0);
if ~isempty(t)
    azi(t) = azi(t) + pi ;
end

u=find(timevinkel>0.0 & tan_azi>0.0 & azi<0.0);
if ~isempty(u)
    azi(u) = azi(u) + pi ;
end


% if (timevinkel<0.0) & (tan_azi<0.0) & (azi>0.0)
%     azi = azi - pi ;
% end
% if (timevinkel<0.0) & (tan_azi>0.0) & (azi>0.0)
%     azi = azi - pi ;
% end
% if (timevinkel>0.0) & (tan_azi<0.0) & (azi<0.0)
%     azi = azi + pi ;
% end
% if (timevinkel>0.0) & (tan_azi>0.0) & (azi<0.0)
%     azi = azi + pi ;
% end


