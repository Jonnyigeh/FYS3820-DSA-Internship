function A=porfyri_return_horisontal_weighted(A,ozone,clod,H,Transmittance)
%A=NOSNOW med flere felter lagt til.
%Transmittance er en kolonne som gir transmittans i f.eks. et solseil. To kolonner, x, y, må dekke
%minst 290-800 nm. Normalt er T=1
%returnerer horisontalkomponenter, vektet med aksjonsspektrum, og kolonneposisjon som svarer til
%ønsket ozonverdi 
T=interp1(Transmittance(:,1),Transmittance(:,2),A.Wavelength);

%Vekt alle spektrumkomponentene med porfyrispekteret
%A=SNOW;eller A=NOSNOW;
flow=find(A.Wavelength < H(1,1));
fmid=find(A.Wavelength >= H(1,1) & A.Wavelength <= H(end,1) );
fhigh=find(A.Wavelength > H(end,1));
ftop=[];

w=[zeros(length(flow),1);...
    H(:,2);...
    zeros(length(fhigh),1);...
    zeros(length(ftop),1)]';
if isfield(A,'CloudOpticalDepth')&& numel(A.CloudOpticalDepth)>1
    w=repmat(w,length(A.CloudOpticalDepth),1);
else
    w=repmat(w,length(A.O3),1);
end

%multipliser med transmittansen
w=w.*repmat(T',size(w,1),1);

% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.EGLO(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EGLO_w=Fw;
% clear Fw;

%Beregn horisontal,diffus ned og reflektert opp komponenter for en gitt ozon-verdi, alternativt
%clod-verdi:
for i=1:length(A.SZA)
    fixSZA_var2(:,:)=A.EDIR(i,:,:);    % [O3/clod x wvl]
    Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
end
A.EDIR_w=Fw;
clear Fw;

for i=1:length(A.SZA)
    fixSZA_var2(:,:)=A.EDN(i,:,:);    % [O3/clod x wvl]
    Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
end
A.EDN_w=Fw;
clear Fw;

for i=1:length(A.SZA)
    fixSZA_var2(:,:)=A.EUP(i,:,:);    % [O3/clod x wvl]
    Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
end
A.EUP_w=Fw;
clear Fw;

for i=1:length(A.SZA)
    fixSZA_var2(:,:)=A.EGLO(i,:,:);    % [O3/clod x wvl]
    Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
end
A.EGLO_w=Fw;
clear Fw;

% 
% 
% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.EDN(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EDN_w=Fw;
% clear Fw;
% 
% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.EDN(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EDN_w=Fw;
% clear Fw;
% 


%Lag globalplott med spektre vektet med HEP
% wz=[zeros(length(flow),1);...
%     H(:,2);...
%     zeros(length(fhigh),1);...
%     zeros(length(ftop),1)]';
% wz=repmat(wz,length(A.SZA),1);size(wz)


% 
% for i=1:length(A.SZA)
%     fixO3(i,:)=A.EGLO(i,pos,:);    % [O3 x wvl]
%     %Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% p=find(A.SZA>=92);
% fixO3(p,:)=0;
% HEP_w=fixO3.*wz;
% %HEP_w(p,:)=0;
% 
% figure(9999)
% subplot(2,1,1)
% plot(A.Wavelength,HEP_w,'-'), grid on
% title(sprintf('HEP-weigthed spectra as function of SZA'))
% xlabel('Wavelength')
% 
% subplot(2,1,2)%plot ratio SZA 85/SZA 40:
% p40=find(A.SZA==40);
% p80=find(A.SZA==80);
% plot(A.Wavelength,HEP_w(p80,:)./HEP_w(p40,:),'r-'), grid on
% title(sprintf('Ratio HEP-weigthed at SZA 80\\circ and 40\\circ'))
% xlabel('Wavelength')
% 
% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.EDN(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EDN_w=Fw;
% clear Fw;
% 
% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.EDIR(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EDIR_w=Fw;
% clear Fw;
% 
% for i=1:length(A.SZA)
%     %fixSZA_var2(:,:)=A.EDIR(i,:,:);    % [O3 x wvl]
%     fixSZA_var2(:,:)=A.EUP(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.EUP_w=Fw;
% clear Fw;
% 
% for i=1:length(A.SZA)
%     fixSZA_var2(:,:)=A.UAVG(i,:,:);    % [O3 x wvl]
%     Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% A.UAVG_w=Fw;
% 
% %Studer også vektete uavg spektre:
% for i=1:length(A.SZA)
%     fixO3Uavg(i,:)=A.UAVG(i,pos,:);    % [O3 x wvl]
%     %Fw(i,:)=trapz(A.Wavelength',(w(:,:).*fixSZA_var2),2 )';
% end
% p=find(A.SZA>=92);
% fixO3Uavg(p,:)=0;
% HEP_Uavg_w=fixO3Uavg.*wz;
% 
% figure(10000)
% subplot(2,1,1)
% plot(A.Wavelength,HEP_Uavg_w,'-'), grid on
% title(sprintf('HEP-weigthed spectra as function of SZA'))
% xlabel('Wavelength')
% 
% subplot(2,1,2)%plot ratio SZA 85/SZA 40:
% p40=find(A.SZA==40);
% p80=find(A.SZA==80);
% plot(A.Wavelength,HEP_Uavg_w(p80,:)./HEP_Uavg_w(p40,:),'r-'), grid on
% title(sprintf('Ratio HEP-weigthed at SZA 80\\circ and 40\\circ'))
% xlabel('Wavelength')
% 
% 
% clear Fw;

pZ=find(A.SZA>90);%sletter uendelig-verdiene for SZA>91
if ~isempty(pZ)
A.EGLO_w(pZ,:)=0;
A.EDN_w(pZ,:)=0;
A.EDIR_w(pZ,:)=0;
A.EUP_w(pZ,:)=0;
%A.UAVG_w(pZ,:)=0;
end

% 
% pw=find(A.Wavelength==wl);
% pO=find(A.O3==ozone);
% %pSZ=find(A.SZA<=90);
% %først en person som dreier i takt med sola:Dvs cos(Az-Gamma)=1
% Az=Gamma;
% cosT=cos(A.SZA*pi/180)*cos(Beta)+sin(A.SZA*pi/180)*sin(Beta)*cos(Az-Gamma);
% p=find(cosT<0);
% cosT(p)=0;
% 
% 
% Bphi=squeeze(A.EDIR(:,pO,pw)).*cosT./cos(A.SZA*pi/180);
% Dphi=squeeze(A.EDN(:,pO,pw)).*(1+cos(Beta))/2;
% Rphi=squeeze(A.EUP(:,pO,pw)).*(1-cos(Beta))/2;
% 
% pO=find(A.O3==ozone);
% 
% if  isfield(A,'CloudOpticalDepth')&& numel(A.CloudOpticalDepth)>1
%     Bphi_w=A.EDIR_w.*repmat(cosT./cos(A.SZA*pi/180),1,length(A.CloudOpticalDepth));
% else
%     Bphi_w=A.EDIR_w.*repmat(cosT./cos(A.SZA*pi/180),1,length(A.O3));
% end
% 
% Dphi_w=A.EDN_w*(1+cos(Beta))/2;
% Rphi_w=A.EUP_w*(1-cos(Beta))/2;
% Gphi_w=A.EGLO_w;
% Sphi_w=A.UAVG_w*4*pi;

%siden det er svært liten ozonavhengighet gjør vi som over og velger en
%ozonkolonne:
% Bphi_w=Bphi_w(:,pO);
% Dphi_w=Dphi_w(:,pO);
% Rphi_w=Rphi_w(:,pO);
% Gphi_w=Gphi_w(:,pO);
% Sphi_w=Sphi_w(:,pO);

