function FLATE = porfyri_returner_horisontalkomponenter(sky_type,clod,ozone,albedo_type,lat,lon,flate_type,Beta,Gamma,action_spec,Transmittance)
% Returnerer irradianskomponenter:
% Transmittance er teltdukens spektrale trasnsmittans (x,y, fra 290 til 800 nm
% Først for en horisontal flate: FLATE.EDIR365 osv.
% Så for en gitt flaterepresentasjon: FLATE.DIR_surface_w osv.

%FLATE = porfyri_returner_horisontalkomponenter(clear_sky,0,350,'beach',41+54/60/3600,12+27/60/3600,'cone',60,180,'HEP')

% %Flatevalg:
%             case 'sphere_element'
%             case 'cone_element'
%             case 'vertical_element'    %Beta =90
%             case 'horizontal_element'    %Beta =0, gamma er hvasomhelst
%             case 'tilted_element'        %gamma=180 er sør, gamma=0 er nord
%             case 'sphere'
%             case 'cone'
%             case 'vertical_cylinder'
% case 'semi_cone_S'
% case 'semi_cone_N'

% % Valg av albedo_type:
%             case 'snow_free'
%             case 'beach'
%             case 'fresh_snow'

% Valg av sky_type:
%             case 'clear_sky';
%             case 'overcast'

% % Valg av action_spec
%     case 'HEP'
%     case 'CIE'

%action_spec='HEP';
%action_spec='CIE';
FLATE.Ozone=ozone;
FLATE.CLOD=clod;
FLATE.Lat=lat;
FLATE.Lon=lon;
FLATE.Sky_type=sky_type;
FLATE.Albedo_type=albedo_type;
FLATE.Flate_type=flate_type;
FLATE.Beta=Beta;
FLATE.Gamma=Gamma;
FLATE.Action_spec=action_spec;

% ozone=350;
% clod=10;%cloud optical depth som det skal gjøres beregninger for  hvis overcast

%sky_type='clear_sky';
%sky_type='overcast'; %Kan bare kjøres et stykke. Brukes kun for å lage de første figurene
%
% %%site='Oslo';%for beregning på skråstilt flate
% %site='Landvik';%for beregning på skråstilt flate
% %site ='Trondheim';
% %site ='Tromsø';
% %site ='Düsseldorf';
% site ='Rome';
% %site ='Equator';
% switch site
%     case 'Oslo'
%         lat=59.938;
%         lon=10.719;
%     case 'Landvik'
%         lat=58.33;
%         lon=8.52;
%     case 'Trondheim'
%         lat=63.42;
%         lon=10.40;
%     case 'Tromsø'
%         lat=69+41/60;
%         lon=18+58/60;
%     case 'Düsseldorf'
%         lat=51+13/60+18/3600;
%         lon=6+46/60+34/3600;
%     case 'Rome'
%         lat=41+54/60/3600;
%         lon=12+27/60/3600;
%     case 'Equator'
%         lat=0;
%         lon=0;
%     otherwise
%         return
% end

switch action_spec
    case 'EPP'
        %if strcmp(action_spec,'HEP')
        HEP=load('C:\Users\jonny\Documents\MATLAB\DSA\HEP_optimizeActionSpectrumTable.txt');%HEP actionspektrum
        Hx=[HEP(1,1):0.5:HEP(end,1)]';
        Hy=interp1(HEP(:,1),HEP(:,2),Hx,'linear');
        H=[Hx';Hy']';
        
    case 'CIE'
        %elseif strcmp(action_spec,'CIE')%Jukser det til for å få beregninger også for UV. UV-veide spektre blir fortsatt kalt HEP
        Hx=[290:0.5:800]';
        flow=find(Hx <= 298);
        fmid=find(Hx > 298 & Hx <= 328 );
        fhigh=find(Hx > 328 & Hx <= 400 );
        ftop=find(Hx > 400);
        
        % CIE-vekter:
        Hy=[ones(length(flow),1);...
            exp(0.094*log(10)*(298*ones(length(fmid),1)-Hx(fmid)));...
            exp(0.015*log(10)*(139*ones(length(fhigh),1)-Hx(fhigh)));...
            zeros(length(ftop),1)];
        H=[Hx';Hy']';
        %
        %     figure(1)
        %     set(gca,'Fontsize',14);
        %     semilogy(H(:,1),H(:,2),'k-','LineWidth',3)
        %     xlabel('Wavelength, nm')
        %     ylabel('Relative units')
        %     title('HEP aksjonsspektrum')
        % end
end
load('uvspec_alb_05_clear_with_aerosols.mat');
DAG=UVSPEC.Dagvar;
DAGVAR=UVSPEC.Dagvar; clear UVSPEC

switch sky_type
    case 'clear_sky';
        load('uvspec_alb_05_Clearsky_porfyri_pp.mat')
        NOSNOW=UVSPEC; clear UVSPEC
        NOSNOW.Dag=DAG;
        NOSNOW.Dagvar=DAGVAR;
        
        load('uvspec_alb_90_Clearsky_porfyri_pp.mat')
        SNOW=UVSPEC; clear UVSPEC
        SNOW.Dag=DAG;
        SNOW.Dagvar=DAGVAR;
        
        %Laster inn en tredje UVSPEC for albedo som ligger mellom snø og snøfritt:
        load('uvspec_alb_30_Clearsky_porfyri_pp.mat');%tilsvarer albedo for sandstrand i den synlige delen av spekeret
        BEACH=UVSPEC; clear UVSPEC
        BEACH.Dag=DAG;
        BEACH.Dagvar=DAGVAR;
        clear DAG DAGVAR
        
    case 'overcast'
        load('uvspec_alb_05_CloudOpticalDepth_cloudcover100_porfyri.mat');
        NOSNOW=UVSPEC; clear UVSPEC
        NOSNOW.Dag=DAG;
        NOSNOW.Dagvar=DAGVAR;
        
        load('uvspec_alb_90_CloudOpticalDepth_cloudcover100_porfyri');
        SNOW=UVSPEC; clear UVSPEC
        SNOW.Dag=DAG;
        SNOW.Dagvar=DAGVAR;
        
        load('uvspec_alb_30_CloudOpticalDepth_cloudcover100_porfyri_pp.mat');%tilsvarer albedo for sandstrand i den synlige delen av spekeret
        BEACH=UVSPEC; clear UVSPEC
        BEACH.Dag=DAG;
        BEACH.Dagvar=DAGVAR;
        clear DAG DAGVAR
    otherwise
        return
end

%Skråflate med vinkel Beta:
NOSNOW=porfyri_return_horisontal_weighted(NOSNOW,ozone,clod,H,Transmittance);
BEACH=porfyri_return_horisontal_weighted(BEACH,ozone,clod,H,Transmittance);
SNOW=porfyri_return_horisontal_weighted(SNOW,ozone,clod,H,Transmittance);

% [NOSNOW,ABphi_w,ADphi_w,ARphi_w,AGphi_w,ASphi_w,ABphi,ADphi,ARphi]=porfyri_calc_irradiance_components(NOSNOW,ozone,wl,Beta*pi/180,Gamma*pi/180,H);%snøfritt skråflate Beta
% [SNOW,BBphi_w,BDphi_w,BRphi_w,BGphi_w,BSphi_w,BBphi,BDphi,BRphi]=porfyri_calc_irradiance_components(SNOW,ozone,wl,Beta*pi/180,Gamma*pi/180,H);%snø skråflate Beta
% [BEACH,BBphi_w,BDphi_w,BRphi_w,BGphi_w,BSphi_w,BBphi,BDphi,BRphi]=porfyri_calc_irradiance_components(BEACH,ozone,wl,Beta*pi/180,Gamma*pi/180,H);%snø skråflate Beta

for k=1:3
    if k==1
        UV=NOSNOW;
    elseif k==2
        UV=SNOW;
    elseif k==3
        UV=BEACH;
    end
    
    %     pw=find(NOSNOW.Wavelength==wl);
    pSZA=find(NOSNOW.SZA<90);
    
    if  isfield(UV,'CloudOpticalDepth')&& numel(UV.CloudOpticalDepth)>1
        pos=find(UV.CloudOpticalDepth==clod);%ISCCPs klassifisering av skytype etter Cloud optical thickness: Cs,Ac,Sc
        %pos=find(A.CloudOpticalDepth==1);%ISCCPs klassifisering av skytype etter Cloud optical thickness: Ci,Ac,Cu
    else
        pos=find(UV.O3==ozone);
    end
    
    
    DecH=(0:1/60:24-1/60)';%hvert 30 minutt i døgnet
    Zen=zeros(length(DecH),365);
    Azi=Zen;
    EDIR365=Zen;
    EDN365=Zen;
    EUP365=Zen;
    for m=1:365 %dagnummer
        dnum=m;
        [aar,maaned,dag]=daynum2date_guv(2014,dnum);
        [Azi(:,m),Zen(:,m)]=solpos_bj(aar,maaned,dag,DecH+1,0,0,lat,lon);
        
        EDIR365(:,m)=interp1(UV.SZA,UV.EDIR_w(:,pos),Zen(:,m)*180/pi,'linear')*UV.Dagvar(dnum);
        EDN365(:,m)=interp1(UV.SZA,UV.EDN_w(:,pos),Zen(:,m)*180/pi,'linear')*UV.Dagvar(dnum);
        EUP365(:,m)=interp1(UV.SZA,UV.EUP_w(:,pos),Zen(:,m)*180/pi,'linear')*UV.Dagvar(dnum);
    end
    
    if k==1
        NOSNOW.SZA365=Zen*180/pi;
        NOSNOW.AZI365=Azi*180/pi;
        NOSNOW.EDIR365=EDIR365;
        NOSNOW.EDN365=EDN365;
        NOSNOW.EUP365=EUP365;
    elseif k==2
        SNOW.SZA365=Zen*180/pi;
        SNOW.AZI365=Azi*180/pi;
        SNOW.EDIR365=EDIR365;
        SNOW.EDN365=EDN365;
        SNOW.EUP365=EUP365;
    elseif k==3
        BEACH.SZA365=Zen*180/pi;
        BEACH.AZI365=Azi*180/pi;
        BEACH.EDIR365=EDIR365;
        BEACH.EDN365=EDN365;
        BEACH.EUP365=EUP365;
    end
    
end %for k=1:3

%alb=5;
%alb=30;
%alb=90;
switch albedo_type
    case 'snow free'
        UV=NOSNOW;
        UV.Albedo=0.05;
    case 'beach'
        UV=BEACH;
        UV.Albedo=0.30;        
    case 'fresh snow'
        UV=SNOW;
        UV.Albedo=0.90;        
end

FLATE.Albedo=UV.Albedo;

SZA365=UV.SZA365;
AZI365=UV.AZI365;
EDIR365=UV.EDIR365;
EDN365=UV.EDN365;
EUP365=UV.EUP365;

FLATE.SZA365=SZA365;
FLATE.AZI365=AZI365;
FLATE.EDIR365=EDIR365;
FLATE.EDN365=EDN365;
FLATE.EUP365=EUP365;

FLATE.DecH=DecH;

%         Dphi365=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
%         Rphi365=EUP365.*(1-cos(Beta*pi/180))/2;

switch flate_type
    case 'sphere_element'
        cosKule=ones(size(SZA365))./cos(SZA365*pi/180);
        pKule=find(cosKule<0);cosKule(pKule)=0;
        BKule=EDIR365.*cosKule;
        DKule=EDN365.*(1+cos(SZA365*pi/180))/2;%flatelemente er hele tiden vinkelrett på sola - Beta=SZA i øyeblikket
        DKule(isnan(DKule)==1)=0;
        RKule=EUP365.*(1-cos(SZA365*pi/180))/2;%flatelemente er hele tiden vinkelrett på sola - Beta=SZA i øyeblikket
        RKule(isnan(RKule)==1)=0;
        
        B=BKule;
        D=DKule;
        R=RKule;
        
        ShapeDir=cosKule;
        ShapeDn=(1+cos(SZA365*pi/180))/2;
        ShapeUp=(1-cos(SZA365*pi/180))/2;
        
    case 'cone_element'
        cosKjegle=cos((SZA365-Beta)*pi/180)./cos(SZA365*pi/180);
        pKjegle=find(cosKjegle<0);cosKjegle(pKjegle)=0;
        BKjegle=EDIR365.*cosKjegle;
        DKjegle=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
        RKjegle=EUP365.*(1-cos(Beta*pi/180))/2;
        
        B=BKjegle;
        D=DKjegle;
        R=RKjegle;
        
        ShapeDir=cosKjegle;
        ShapeDn=ones(size(SZA365))*(1+cos(Beta*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(Beta*pi/180))/2;
        
    case 'vertical_element'    %Beta =90
        cosV=tan(SZA365*pi/180);
        pV=find(cosV<0);cosV(pV)=0;
        BV=EDIR365.*cosV;
        DV=EDN365.*(1+cos(90*pi/180))/2;
        RV=EUP365.*(1-cos(90*pi/180))/2;
        
        B=BV;
        D=DV;
        R=RV;
        
        ShapeDir=cosV;
        ShapeDn=ones(size(SZA365))*(1+cos(90*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(90*pi/180))/2;
        
    case 'horizontal_element'    %Beta =0, gamma er hvasomhelst
        cosH=ones(size(SZA365));%Horisontal 
        BH=EDIR365;
        DH=EDN365;
        RH=EUP365*0;
        
        B=BH;
        D=DH;
        R=RH;
        
        ShapeDir=cosH;
        ShapeDn=ones(size(SZA365));
        ShapeUp=ones(size(SZA365))*(1-cos(0*pi/180))/2;        
        
    case 'tilted_element'        %gamma=180 er sør, gamma=0 er nord
        cosTilt=(cos(SZA365*pi/180).*cos(Beta*pi/180)+sin(SZA365*pi/180).*sin(Beta*pi/180).*cos(AZI365*pi/180 - Gamma))./cos(SZA365*pi/180);
        pcosTilt=find(cosTilt<0);cosTilt(pcosTilt)=0;
        BTilt=EDIR365.*cosTilt;
        DTilt=EDN365.*(1+cos(Beta*pi/180))/2;
        RTilt=EUP365.*(1-cos(Beta*pi/180))/2;
        
        B=BTilt;
        D=DTilt;
        
        ShapeDir=cosTilt;
        ShapeDn=ones(size(SZA365))*(1+cos(Beta*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(Beta*pi/180))/2;    
        
    case 'sphere'
        cosKuleMean=ones(size(SZA365))./(4*sin((90-SZA365)*pi/180));
        pKuleMean=find(cosKuleMean<0);cosKuleMean(pKuleMean)=0;
        BKuleMean=EDIR365.*cosKuleMean;
%         DKuleMean=EDN365.*(1+0)/2;
%         RKuleMean=EUP365.*(1-0)/2;
        DKuleMean=EDN365.*1/2;%all stråling som utgår fra bakken (horisontalplan) går gjennom halvkula over. E(halvkule)=( A(sirkel på bakken)/A(halvkule) )*Ehorisontal =Ehorisontal/2
        RKuleMean=EUP365.*1/2;
        
        B=BKuleMean;
        D=DKuleMean;
        R=RKuleMean;
        
        ShapeDir=cosKuleMean;
%         ShapeDn=ones(size(SZA365))*(1+0)/2;
%         ShapeUp=ones(size(SZA365))*(1-0)/2;
        ShapeDn=ones(size(SZA365))/2;
        ShapeUp=ones(size(SZA365))/2;
        
    case 'cone'
        cosKjegleMean=zeros(size(SZA365));%shapefaktor
        
        p_o=find((90-SZA365)>Beta);%høyvinkelområde
        if ~isempty(p_o)
            cosKjegleMean(p_o)=cos(Beta*pi/180);
        end
        
        %lav-vinkelområdet:
        p_u=find((90-SZA365)<=Beta & (90-SZA365)>0);
        th_0=ones(size(SZA365))*pi/2;
        th_0(p_u)=real(acos( tan( (90-SZA365(p_u))*pi/180).*cot(Beta*pi/180) ) );
        th_0(p_o)=0;%når sola står over kjeglen er skyggesektoren th_0 lik 0 grader
        
        sin_th_0=0*SZA365;
        sin_th_0=sin(th_0);
        
        cosKjegleMean(p_u)=(1-th_0(p_u)/pi).*cos( Beta*pi/180) + sin(Beta*pi/180).* cot((90-SZA365(p_u))*pi/180).*sin_th_0(p_u)./pi;
        pcosKjegleMean=find(cosKjegleMean<0);cosKjegleMean(pcosKjegleMean)=0;
        BKjegleMean=EDIR365.*cosKjegleMean;
        
        DKjegleMean=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
        RKjegleMean=EUP365.*(1-cos(Beta*pi/180))/2;
        
        B=BKjegleMean;
        D=DKjegleMean;
        R=RKjegleMean;
        
        ShapeDir=cosKjegleMean;
        ShapeDn=ones(size(SZA365))*(1+cos(Beta*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(Beta*pi/180))/2;

    case 'semi_cone_S' %halv kjegle, orientert med krumme siden mot sør (som om nesa peker mot sør og får bidrag +/- 90grader
        %cos_semi_cone_S=1+tan(SZA365*pi/180)*tan(Beta*pi/180)*2/pi;%shapefaktor
        %cos_semi_cone_S=cos(Beta*pi/180)+tan(SZA365*pi/180)*sin(Beta*pi/180)*2/pi;%shapefaktor
        cos_semi_cone_S=0.5+tan(SZA365*pi/180)*tan(Beta*pi/180)/pi;%ny variant, avledet for totalfluks fordelt på en halv kjegleflate
        p=find(cos_semi_cone_S<0);
        cos_semi_cone_S(p)=0;
%         p_o=find((90-SZA365)>Beta);%høyvinkelområde
%         if ~isempty(p_o)
%             cosKjegleMean(p_o)=cos(Beta*pi/180);
%         end
%         
%         %lav-vinkelområdet:
%         p_u=find((90-SZA365)<=Beta & (90-SZA365)>0);
%         th_0=ones(size(SZA365))*pi/2;
%         th_0(p_u)=real(acos( tan( (90-SZA365(p_u))*pi/180).*cot(Beta*pi/180) ) );
%         th_0(p_o)=0;%når sola står over kjeglen er skyggesektoren th_0 lik 0 grader
%         
%         sin_th_0=0*SZA365;
%         sin_th_0=sin(th_0);
%         
%         cosKjegleMean(p_u)=(1-th_0(p_u)/pi).*cos( Beta*pi/180) + sin(Beta*pi/180).* cot((90-SZA365(p_u))*pi/180).*sin_th_0(p_u)./pi;
%         pcosKjegleMean=find(cosKjegleMean<0);cosKjegleMean(pcosKjegleMean)=0;
        B_semi_cone_S=EDIR365.*cos_semi_cone_S;
        
        D_semi_cone_S=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
        R_semi_cone_S=EUP365.*(1-cos(Beta*pi/180))/2;
        
        B=B_semi_cone_S;
        D=D_semi_cone_S;
        R=R_semi_cone_S;
        
        ShapeDir=cos_semi_cone_S;
        ShapeDn=ones(size(SZA365))*(1+cos(Beta*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(Beta*pi/180))/2;


    case 'semi_cone_N' %halv kjegle, orientert med krumme siden mot nord (som om nesa peker mot nord og får bidrag +/- 90grader

        %cos_semi_cone_N=zeros(size(SZA365)); %0 når sola ikke når over kjeglekanten        
        cos_semi_cone_N=cos(Beta*pi/180)-tan(SZA365*pi/180)*sin(Beta*pi/180)*2/pi;%shapefaktor
        p=find(cos_semi_cone_N<0);
        cos_semi_cone_N(p)=0;
                
%         p_o=find((90-SZA365)>Beta);%høyvinkelområde
%         if ~isempty(p_o)
%             %cos_semi_cone_N(p_o)=1+tan(SZA365(p_o)*pi/180)*tan(Beta*pi/180)*2/pi;
%             cos_semi_cone_N(p_o)=cos(Beta*pi/180)+tan(SZA365(p_o)*pi/180)*sin(Beta*pi/180)*2/pi;%shapefaktor
%             %cos_semi_cone_N=reshape(cos_semi_cone_N,
%         end
% %         
% %         %lav-vinkelområdet:
% %         p_u=find((90-SZA365)<=Beta & (90-SZA365)>0);
% %         th_0=ones(size(SZA365))*pi/2;
% %         th_0(p_u)=real(acos( tan( (90-SZA365(p_u))*pi/180).*cot(Beta*pi/180) ) );
% %         th_0(p_o)=0;%når sola står over kjeglen er skyggesektoren th_0 lik 0 grader
% %         
% %         sin_th_0=0*SZA365;
% %         sin_th_0=sin(th_0);
% %         
% %         cosKjegleMean(p_u)=(1-th_0(p_u)/pi).*cos( Beta*pi/180) + sin(Beta*pi/180).* cot((90-SZA365(p_u))*pi/180).*sin_th_0(p_u)./pi;
% %         pcosKjegleMean=find(cosKjegleMean<0);cosKjegleMean(pcosKjegleMean)=0;
        B_semi_cone_N=EDIR365.*cos_semi_cone_N;
        
        D_semi_cone_N=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
        R_semi_cone_N=EUP365.*(1-cos(Beta*pi/180))/2;
        
        B=B_semi_cone_N;
        D=D_semi_cone_N;
        R=R_semi_cone_N;
        
        ShapeDir=cos_semi_cone_N;
        ShapeDn=ones(size(SZA365))*(1+cos(Beta*pi/180))/2;
        ShapeUp=ones(size(SZA365))*(1-cos(Beta*pi/180))/2;
        
    case 'vertical_cylinder'
        cosVCyl=cot( (90-SZA365)*pi/180)/pi;%vertikal sirkulær sylinder
        pcosVCyl= cosVCyl<0;cosVCyl(pcosVCyl)=0;
        BVCyl=EDIR365.*cosVCyl;
        
        DVCyl=EDN365.*(1+0)/2;%flaten hele tiden vendt en fast vinkel Beta
        RVCyl=EUP365.*(1-0)/2;
        
        B=BVCyl;
        D=DVCyl;
        R=RVCyl;
        
        ShapeDir=cosVCyl;
        ShapeDn=ones(size(SZA365))*(1+0)/2;
        ShapeUp=ones(size(SZA365))*(1-0)/2;                 
end

FLATE.DIR_surface_w=B;
FLATE.DDN_surface_w=D;
FLATE.RUP_surface_w=R;

FLATE.ShapeDir=ShapeDir;
FLATE.ShapeDn=ShapeDir;
FLATE.ShapeUp=ShapeUp;

%
%         cosKule=ones(size(SZA365))./cos(SZA365*pi/180);
%         cosKjegle=cos((SZA365-Beta)*pi/180)./cos(SZA365*pi/180);
%         cosSor=(cos(SZA365*pi/180).*cos(Beta*pi/180)+sin(SZA365*pi/180).*sin(Beta*pi/180).*cos(AZI365*pi/180 - 0))./cos(SZA365*pi/180);%Fast mot sør
%         cosNor=(cos(SZA365*pi/180).*cos(Beta*pi/180)+sin(SZA365*pi/180).*sin(Beta*pi/180).*cos(AZI365*pi/180 - pi))./cos(SZA365*pi/180);%Fast mot nord
%         cosKuleMean=ones(size(cosKule))./(4*sin((90-SZA365)*pi/180));
%         cosV=tan(SZA365*pi/180);
%         cosH=ones(size(cosKule));%Horisontal
%         cosVCyl=cot( (90-SZA365)*pi/180)/pi;%vertikal sirkulær sylinder
%
%         cosKjegleMean=0*cosV;%shapefaktor
%
%         p_o=find((90-SZA365)>Beta);%høyvinkelområde
%         if ~isempty(p_o)
%             cosKjegleMean(p_o)=cos(Beta*pi/180);
%         end
%
%         %lav-vinkelområdet:
%         p_u=find((90-SZA365)<=Beta & (90-SZA365)>0);
%         th_0=ones(size(cosV))*pi/2;
%         th_0(p_u)=real(acos( tan( (90-SZA365(p_u))*pi/180).*cot(Beta*pi/180) ) );
%         th_0(p_o)=0;%når sola står over kjeglen er skyggesektoren th_0 lik 0 grader
%
%         sin_th_0=0*cosV;
%         sin_th_0=sin(th_0);
%
%         cosKjegleMean(p_u)=(1-th_0(p_u)/pi).*cos( Beta*pi/180) + sin(Beta*pi/180).* cot((90-SZA365(p_u))*pi/180).*sin_th_0(p_u)./pi;%        (cos(Beta*pi/180)+sin(Beta*pi/180*cot(SZA365(p_u)*pi/180).*sin_th_0).*(pi-th_0)/pi);
%
%
%
%
%

%
%
%
%
%         % Dphi365Kule=EDN365.*(1+cos(SZA365*pi/180))/2;%flaten justerer seg så den hele tiden er vendt mot sola: Beta=SZA
%         % Rphi365Kule=EUP365.*(1-cos(SZA365*pi/180))/2;
%
%         %
%         % %samme for snø:
%         % Bphi365S=EDIR365S.*cosT365; %./cos(SZA*pi/180);
%         % Bphi365SSouth=EDIR365S.*cosT365South; %./cos(SZA*pi/180);
%         % Bphi365SNorth=EDIR365S.*cosT365North; %./cos(SZA*pi/180);%vendt sør
%         % Bphi365SEast=EDIR365S.*cosT365East; %./cos(SZA*pi/180);%vendt sør
%         % Bphi365SWest=EDIR365S.*cosT365West; %./cos(SZA*pi/180);%avg azi
%         % Bphi365SAvg=EDIR365S.*cosT365Avg; %./cos(SZA*pi/180);%avg azi
%         %
%         % BKuleS=EDIR365S.*cosKule;
%         % BKjegleS=EDIR365S.*cosKjegle;
%         % BSorS=EDIR365S.*cosSor;
%         % BNorS=EDIR365S.*cosNor;
%         % BAzimeanS=EDIR365S.*cosAzimean;
%         % BVS=EDIR365S.*cosV;
%         %
%         % Dphi365S=EDN365S.*(1+cos(Beta*pi/180))/2;
%         % Rphi365S=EUP365S.*(1-cos(Beta*pi/180))/2;
%         %
%         % Dphi365KuleS=EDN365S.*(1+cos(SZA365*pi/180))/2;%flaten justerer seg så den hele tiden er vendt mot sola: Beta=SZA
%         % Rphi365KuleS=EUP365S.*(1-cos(SZA365*pi/180))/2;
%
%         %**
%         %EKule=BKule+Dphi365Kule+Rphi365Kule;EKule(find(isnan(EKule)==1))=0;
%         EKule=BKule+DphiKule+RphiKule;EKule(find(isnan(EKule)==1))=0;
%         EKjegle=BKjegle+Dphi365+Rphi365;EKjegle(find(isnan(EKjegle)==1))=0;
%         ESor=BSor+Dphi365+Rphi365;ESor(find(isnan(ESor)==1))=0;
%         ENor=BNor+Dphi365+Rphi365;ENor(find(isnan(ENor)==1))=0;
%         EKuleMean=BKuleMean+Dphi365_halvkule+Rphi365_halvkule;EKuleMean(find(isnan(EKuleMean)==1))=0;
%         EV=BV+EDN365/2+EUP365/2;EV(find(isnan(EV)==1))=0;
%         EKjegleMean=BKjegleMean+Dphi365+Rphi365;EKjegleMean(find(isnan(EKjegleMean)==1))=0;
%         EH=BH+EDN365+EUP365;EH(find(isnan(EH)==1))=0;
%         EVCyl=BVCyl+Dphi365_halvkule+Rphi365_halvkule;EVCyl(find(isnan(EVCyl)==1))=0;
%         % EKuleS=BKuleS+Dphi365S+Rphi365S;EKuleS(find(isnan(EKuleS)==1))=0;
%         % EKjegleS=BKjegleS+Dphi365S+Rphi365S;EKjegleS(find(isnan(EKjegleS)==1))=0;
%         % ESorS=BSorS+Dphi365S+Rphi365S;ESorS(find(isnan(ESorS)==1))=0;
%         % ENorS=BNorS+Dphi365S+Rphi365S;ENorS(find(isnan(ENorS)==1))=0;
%         % EAzimeanS=BAzimeanS+Dphi365S+Rphi365S;EAzimeanS(find(isnan(EAzimeanS)==1))=0;
%         % EVS=BVS+Dphi365S+Rphi365S;EVS(find(isnan(EVS)==1))=0;
%
%         %Plott irradians vs tid for de 5 flaterepresentasjonene:vVelg dag 175
%         noonP=zeros(size(EKule,2),1);
%         for i=1:size(EKule,2)
%             noonP(i)=find(EKule(:,i)==max(EKule(:,i)));
%         end
%
%
%         pd=172;%21 juni, sommersolhverv
%         pd_mars=dagnr(2015,3,1);
%         %newH=DecH+(DecH(noonP(pd))-12);
%         newH=DecH+(12-DecH(noonP(pd)));
%
%         figure(1) %plotter HEP for ulike flategeometrier
%         plot(newH,EKule(:,pd),'r-')
%         hold on,grid on
%         plot(newH,EKjegle(:,pd),'b-')
%         %plot(newH,ESor(:,pd),'k-')
%         plot(newH,EV(:,pd),'m-')
%         plot(newH,ENor(:,pd),'k-')
%         plot(newH,EH(:,pd),'k--')%(1:10:end)
%
%         plot(newH(1:30:end),EKuleMean(1:30:numel(DecH),pd),'ko-')
%         plot(newH(1:30:end),EKjegleMean(1:30:numel(DecH),pd),'r--','LineWidth',2)
%         plot(newH(1:30:end),EVCyl(1:30:numel(DecH),pd),'b*--')
%         hold off
%
%         xlabel('Time')
%         ylabel(sprintf('%s, W/m^2',action_spec))
%         if strcmp(action_spec,'HEP')
%             title(sprintf('Sun and surface representations in %s\nSummer solstice and Slope %3.0f\\circ ',site, Beta'),'FontSize',12)
%             hl=legend('Sphere-element following sun',sprintf('Cone-element %3.0f\\circ rotating',Beta),'Vertical rotating',sprintf('Tilted North%3.0f\\circ',Beta),sprintf('Horisontal'),sprintf('Sphere'),sprintf('Cone%3.0f\\circ',Beta),sprintf('Vertical Cylinder'),'Location','northeastoutside');
%             axis([2 22 0 ceil(max(EH(:,pd))/10)*10])
%         elseif strcmp(action_spec,'CIE')
%             title(sprintf('Sun and surface representations in %s\nSummer solstice and Slope %3.0f\\circ ',site, Beta'),'FontSize',12)
%             hl=legend('Sphere-element following sun',sprintf('Cone-element %3.0f\\circ rotating',Beta),'Vertical rotating',sprintf('Tilted North%3.0f\\circ',Beta),sprintf('Horisontal'),sprintf('Sphere'),sprintf('Cone%3.0f\\circ',Beta),sprintf('Vertical Cylinder'),'Location','northeastoutside');
%             axis([2 22 0 ceil(max(EH(:,pd))*100)/100])
%         end
%
%         % figure,plot(newH,BKjegleMean(:,pd),'ko--')
%         %
%         % figure,plot(DecH,th_0(:,pd),'-'),grid
%
%         set(gca,'Fontsize',12);
%         set(hl,'FontSize',10)
%         hFig=figure(1);
%         hgexport(hFig,'-clipboard')
%
%         figure
%         %
%         %
%         %
%         % % figure(100)
%         % % plot(DecH,L1(:,pd)./nevner(:,pd),'r-')
%         % % hold on, grid on
%         % % plot(DecH,L2(:,pd)./nevner(:,pd),'b-')
%         % % plot(DecH,L3(:,pd)./nevner(:,pd),'g-')
%         % % plot(DecH,sum_cos_kjegle(:,pd),'ko--')
%         % % plot(DecH,hoyvinkel(:,pd),'m-','LineWidth',3)
%         % %
%         % % hold off
%         % % legend('sirkel','trekant','sektor','sum skygge','helbelyst')
%         % % axis([0 24 0 5])
%         %
%         % %vekter med direkte, diffus ned og opp komponenter:
%         % BK=EDIR365.*sum_cos_kjegle;
%         % DK=EDN365.*(1+cos(Beta*pi/180))/2;%flaten hele tiden vendt en fast vinkel Beta
%         % RK=EUP365.*(1-cos(Beta*pi/180))/2;
%         % EK=BK+DK+RK;
%         %
%         % figure(101)
%         % plot(DecH,BK(:,pd),'k-')
%         % hold on
%         % plot(DecH,DK(:,pd),'r-')
%         % plot(DecH,RK(:,pd),'g-')
%         % plot(DecH,EK(:,pd),'k.-')
%         % hold off
%         % legend('direkte','diffus ned','reflektert opp','Kjeglesum')
%         %
%         % figure(1)
%         % hold on
%         % plot(newH,EK(:,pd),'r.-')
%         % hold off
%         %
%         % figure
%         % test=cot((90-SZA365)*pi/180);
%         % plot(DecH,test(:,pd))
%         % grid on
%         % axis([0 24 -5 5])