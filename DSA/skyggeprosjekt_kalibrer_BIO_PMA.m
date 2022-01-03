function skyggeprosjekt_kalibrer_BIO_PMA
% Laster opp måledata fra perioden 3 Solar light instrumenter stod på
% solplattformen på Østerås sammen med UV-nettverkets referanseinstrumentet
% GUV9280. Beregner et sett av kalibreringsfunkjsjoner som gjør at solar light
% insturmentene måler mest mulig likt med referansen.
% Kalibreringsfunksjonen: Korrigert måling [UVI]=Råverdi
% [MED/H]*kalibreringsfaktor [UVI/MED/H]*korreksjon som funksjon av zenithvinkelen på måletidspunktet

%Første del er to ulike opplastingsfunksjoner for referansedata: Første
%trinn går gjennom en fullstendig beregning av alle målinger gjennom hele året, som bygger opp en struct array
%med 11 ulike doseprodukter UVI, UVA, UVB etc, og tilsvarende klarværsmodellerte verdier for dagens ozon-tykkelse.

%Andre del, gir en forenklet beregning av kun UVI, uten å ta hensyn til
%zenitvinkelen på måletidspunktet. Gjør dette for å se hvor stor forskjell
%det gir i UVI. Her ble det liten forskjell, men velger likevel den mer
%omstendelige beregningen i første trinn.

%Går deretter over til å laste opp rådat fra PMA og UV-Biometrene, beregner
%kaliberingsfunksjoner, og lagrer til struct PMA for det håndholdte insturmentet, BIO for de to UV-Biometrene (serienummer 0616 tak og 3724 bakke), og RES, som er Østerås-referansen.

clean %clc,clear,close all

station='oesteraas';
stasjons_id=1;
year=2021;

load_mat_file=0; %Østerås-GUV data er allere lagret på mat-form
%load_mat_file=1;
if load_mat_file
    RES=database_get_one_year(year,stasjons_id);%år og stasjons_id
    lat=RES.Latitude;%59.946;%Østerås
    long=RES.Longitude; %10.598;
    
    RES=database_get_one_year_SZA(RES);%Legger til SZA, korrigerer kanaldata for temperatur og offset
    
    RES=database_get_one_year_CLOD(RES);
    
    [OMI,~]=http_upload_ozone_overpass(station,stasjons_id);
    ozone_get_global_maps(year);%Laster ned alle dager i året med OMI griddete data
    
    %slett indeksfiler fra html nedlastingen:
    listing = dir(sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Ozone\\ozone_toms_omi\\Y%04d\\index*',year));
    
    for k = 1 : size(listing,1)
        baseFileName = listing(k).name;
        fullFileName = fullfile(listing(k).folder, baseFileName);
        fprintf(1, 'Now deleting %s\n', fullFileName);
        delete(fullFileName);
    end
    
    RES=database_get_one_year_ozone(RES,OMI);
    close all
    RES=database_get_one_year_doses(RES);
    close all
    figure(1)
    plot(RES.JT,RES.CIE_FARIN*40,'r-'),grid on,grid minor
    
    tmp='RES';
    mat_file=sprintf('L:\\Optisk Lab\\Uvnet\\Prosjekter\\Skyggeprosjekt\\%04d\\Kalibreringsdata_Oesteras\\GUV_%s_%04d.mat',year,station,year);
    save(mat_file,tmp);
    
else
    mat_file=sprintf('L:\\Optisk Lab\\Uvnet\\Prosjekter\\Skyggeprosjekt\\%04d\\Kalibreringsdata_Oesteras\\GUV_%s_%04d.mat',year,station,year);
    load(mat_file);
end

DT=RES.DateHour;
Y=str2num(DT(:,1:4));
M=str2num(DT(:,5:6));
D=str2num(DT(:,7:8));
h=str2num(DT(:,10:11));
m=str2num(DT(:,13:14));
s=zeros(size(m));
d=datestr([Y M D h m s]);
DTnum=datenum(d);
RES.dvec=datetime(DTnum, 'ConvertFrom', 'datenum');

clear Y M D h m s DTnum d

%%Gjentar beregning av dvec for klarværsmodellert:
DecH=RES.MODEL.HTime;
h=floor(DecH);
m=floor((DecH-h)*60);
s=floor((DecH-(h+m/60))*60);
Y=[];M=[];D=[];
%Y=RES.MODEL.Year*ones(size(DecH));
dager=unique(RES.MODEL.Daynum);
for q=1:numel(dager)
    p=find(RES.MODEL.Daynum==dager(q));
    [aar,maaned,dag]=daynum2date_guv(year,dager(q));
    Y=[Y;aar*ones(numel(p),1)];
    M=[M;maaned*ones(numel(p),1)];
    D=[D;dag*ones(numel(p),1)];
end

dm=datestr([Y M D h m s]);
DTnum=datenum(dm);
RES.MODEL.dvec=datetime(DTnum, 'ConvertFrom', 'datenum');
tmp='RES';
mat_file=sprintf('L:\\Optisk Lab\\Uvnet\\Prosjekter\\Skyggeprosjekt\\%04d\\Kalibreringsdata_Oesteras\\GUV_%s_%04d.mat',year,station,year);
save(mat_file,tmp);
clear Y M D h m s DTnum d

lat=RES.Latitude;
long=RES.Longitude;

%Referansedata: GUV på Østerås
ref_tab = skyggeprosjekt_import_GUVDB("L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\Testsett_Oesteraas_data_august_2021.txt", [3, Inf]);
DT=char(ref_tab.Datotid);
YG=str2num(DT(:,1:4));
MG=str2num(DT(:,5:6));
DG=str2num(DT(:,7:8));
hG=str2num(DT(:,10:11));
mG=str2num(DT(:,13:14));
sG=zeros(size(mG));
dG=datestr([YG MG DG hG mG sG]);
DTnum_G=datenum(dG);
dvecG=datetime(DTnum_G, 'ConvertFrom', 'datenum');

REF.dvec=dvecG;
REF.Y=YG;
REF.M=MG;
REF.D=DG;
REF.h=hG;
REF.m=mG;
REF.s=sG;
REF.CIE=ref_tab.Cie;
REF.Dnum=day(datetime([YG';MG';DG']'),'dayofyear');

%Beregn SZA og AZI ut fra UTC-verdier:
REF.SZA=zeros(size(REF.h)); % zenith vinkel
REF.AZI=zeros(size(REF.h)); % zenith vinkel
for q=1:numel(REF.h)
    [azi,zen]=solpos_bj(REF.Y(q),REF.M(q),REF.D(q),REF.h(q)+1,REF.m(q),REF.s(q),lat,long);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
    REF.SZA(q)=zen*180/pi;
    REF.AZI(q)=azi*180/pi;
end

%Sammenlign UVI hentet rett fra databasen mot beregning i matlab:
pos=find(RES.dvec>=REF.dvec(1) & RES.dvec<=REF.dvec(end));
figure(10)
subplot(2,1,1)
plot(RES.dvec(pos),RES.CIE_FARIN(pos)*40,'r-'),hold on,grid on,grid minor
plot(REF.dvec,REF.CIE*40,'k-'),hold off
legend('matlab','database')

subplot(2,1,2)
[~,ia,ib]=intersect(RES.dvec,REF.dvec);
plot(RES.dvec(ia),RES.CIE_FARIN(ia)./REF.CIE(ib),'k.-')
ylim([0.95 1.05])

%Ser at begge samsvarer +/-1 prosent,som er ok.
figure(1)
subplot(2,1,1)
plot(RES.dvec(pos),RES.CIE_FARIN(pos)*40,'r-'),grid on
legend('Referanse Østerås-GUV')
xlabel('Tid, UTC=Lokaltid-2H')
ylabel('UVI')

clear REF %Trenger ikke REF mer
% I fortsettelsen brukes heretter RES og følgende data felter:
%RES.dvec,RES.CIE_FARIN*40, RES.SZA
%RES.MODEL.dvec,RES.MODEL.UVI

%Lagrer først feltene til to filer: målt og modellert, for perioden juli-august 2021:
save2txtfile=0;
if save2txtfile
    fname_meas='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\Refdata_Oesteraas juli-august_2021.txt';
    [filepath,name,ext] = fileparts(fname_meas);
    tmp_name=sprintf('C:\\Users\\bjohnsen\\Downloads\\%s%s',name,ext);    %lagrer midlertidig lokalt på downloads
    
    fid=fopen(tmp_name,'wt');
    headtekst=sprintf('%% Referansedata %s - %d, Tid=UTC+0',station,year);
    headtekst_2=sprintf('%% Date Time UVI SZA');
    fprintf(fid,'%s\n',headtekst);
    fprintf(fid,'%s\n',headtekst_2);
    pos_m=find(RES.dvec>=datetime(datenum([2021,7,25,0,0,0]), 'ConvertFrom', 'datenum') & ...
        RES.dvec<=datetime(datenum([2021,8,25,0,0,0]), 'ConvertFrom', 'datenum'));
    for row=1:numel(pos_m)
        fprintf(fid,'%08s %02s:%02s %6.4f %6.3f\n',RES.DateHour(pos_m(row),1:8),RES.DateHour(pos_m(row),10:11),RES.DateHour(pos_m(row),13:14),RES.CIE_FARIN(pos_m(row))*40,RES.SZA(pos_m(row)));
    end
    fclose(fid)
    
    status = movefile(tmp_name,filepath);%flytter fila over på fellesområdet:
    
    fname_modelled='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\klarvaersdata_Oesteraas juli-august_2021.txt';
    [filepath,name,ext] = fileparts(fname_modelled);
    tmp_name=sprintf('C:\\Users\\bjohnsen\\Downloads\\%s%s',name,ext);    %lagrer midlertidig lokalt på downloads
    fid=fopen(tmp_name,'wt');
    headtekst=sprintf('%% Klarværsmodellerte data %s - %d, Tid=UTC+0',station,year);
    headtekst_2=sprintf('%% Date Time UVI SZA');
    fprintf(fid,'%s\n',headtekst);
    fprintf(fid,'%s\n',headtekst_2);
    pos_modelled=find(RES.MODEL.dvec>=datetime(datenum([2021,7,25,0,0,0]), 'ConvertFrom', 'datenum') & ...
        RES.MODEL.dvec<=datetime(datenum([2021,8,25,0,0,0]), 'ConvertFrom', 'datenum'));
    
    %Regner om dvec til datestring
    ds=datetime(datenum(RES.MODEL.dvec(pos_modelled)),'ConvertFrom', 'datenum', 'Format', 'ddMMyyyy hh:mm');
    
    for row=1:numel(pos_modelled)
        fprintf(fid,'%s %6.4f %6.3f\n',ds(row),RES.MODEL.UVI(pos_modelled(row)),RES.MODEL.SZA(pos_modelled(row)));
    end
    fclose(fid)
    
    status = movefile(tmp_name,filepath);%flytter fila over på fellesområdet:
    
end

%%%%%%%%%%%%%%%%% Last inn PMA %%%%%%%%%%%%%%%%%%%%%%%%%%%
PMAUVB = Skyggeprosjekt_import_PMA_file("L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\PMA_ery_SN05855_2021.log", [2, Inf]);
PMAUVB(1,:)

[Y,M,D] = ymd(PMAUVB.loggerdato);
[h,m,s] = hms(PMAUVB.loggertid);
ds=datestr([Y M D h m s]);

DTnum_PMA=datenum(ds);
p=find(unique(DTnum_PMA));
dvec=datetime(DTnum_PMA, 'ConvertFrom', 'datenum');

PMA.dvec=dvec;
PMA.Y=Y;
PMA.M=M;
PMA.D=D;
PMA.h=h;
PMA.m=m;
PMA.s=s;
PMA.Raw=PMAUVB.mleverdi;
PMA.Dnum=day(datetime([Y';M';D']'),'dayofyear');

%Beregn SZA og AZI ut fra UTC-verdier:
PMA.SZA=zeros(size(PMA.h)); % zenith vinkel
PMA.AZI=zeros(size(PMA.h)); % zenith vinkel
for q=1:numel(PMA.h)
    [azi,zen]=solpos_bj(PMA.Y(q),PMA.M(q),PMA.D(q),PMA.h(q)+1,PMA.m(q),PMA.s(q),lat,long);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
    PMA.SZA(q)=zen*180/pi;
    PMA.AZI(q)=azi*180/pi;
end

figure(1)
subplot(2,1,2)
plot(PMA.dvec,PMA.Raw,'r.-'),grid on
legend('PMA')
xlabel('Tid, UTC=Lokaltid-2H')
ylabel('MED/H')
title('Målinger på Østerås')


%%%%%%%%%%%%%%%%%% Last opp UV-Bios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UVBios = Skyggeprosjekt_import_UVBios_file("L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\UV_Bios_Oesteraas_2021.log", [2, Inf]);
UVBios(1:10,:)
%Date         Time       AVGBio0616    AVGFin3724    Temp1    Temp2
[YB,MB,DB] = ymd(UVBios.Date);

Timestring=char(UVBios.Time);
hB=str2num(Timestring(:,1:2));
mB=str2num(Timestring(:,4:5));
sB=zeros(size(mB));
dsB=datestr([YB MB DB hB mB sB]);%datestring UVBios
DTnum_UVBios=datenum(dsB);
dvecUVBios=datetime(DTnum_UVBios, 'ConvertFrom', 'datenum');
p=find(unique(DTnum_UVBios));

BIO.dvec=dvecUVBios;
BIO.Raw0616=UVBios.AVGBio0616;
BIO.Raw3724=UVBios.AVGFin3724;
BIO.Y=YB;
BIO.M=MB;
BIO.D=DB;
BIO.h=hB;
BIO.m=mB;
BIO.s=sB;
BIO.Dnum=day(datetime([YB';MB';DB']'),'dayofyear');

%Beregn SZA og AZI ut fra UTC-verdier:
BIO.SZA=zeros(size(BIO.h)); % zenith vinkel
BIO.AZI=zeros(size(BIO.h)); % zenith vinkel
for q=1:numel(BIO.h)
    [azi,zen]=solpos_bj(BIO.Y(q),BIO.M(q),BIO.D(q),BIO.h(q)+1,BIO.m(q),BIO.s(q),lat,long);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
    BIO.SZA(q)=zen*180/pi;
    BIO.AZI(q)=azi*180/pi;
end


figure(1)
subplot(2,1,2)
hold on
plot(BIO.dvec,BIO.Raw0616,'k-'),grid on,grid minor
plot(BIO.dvec,BIO.Raw3724,'g-'), hold off
ylabel('MED/H')
legend('PMA','Tak0616','Bakke3724')
xlabel('Tid, UTC=Lokaltid-2H')


%Vis GUV, PMA og UVBios for samme tidspunkt
%[~,ia,ib]=intersect(PMA.dvec,BIO.dvec);
[~,ia,ib]=intersect(PMA.dvec,RES.dvec);
[~,ic,id]=intersect(BIO.dvec,RES.dvec);

Rat_PMA=(RES.CIE_FARIN(ib)*40)./PMA.Raw(ia);SZ_PMA=PMA.SZA(ia);
Rat_Bio0616=(RES.CIE_FARIN(id)*40)./BIO.Raw0616(ic);SZ_Bio=BIO.SZA(ic);
Rat_Bio3724=(RES.CIE_FARIN(id)*40)./BIO.Raw3724(ic);

posZ_PMA=find(SZ_PMA<=70);
posZ_Bio=find(SZ_Bio<=70);

Calfact_PMA_tmp=mean(Rat_PMA(posZ_PMA));%Faktor å multipliser PMA med for å få UVI
Calfact_Bio0616_tmp=mean(Rat_Bio0616(posZ_Bio));
Calfact_Bio3724_tmp=mean(Rat_Bio3724(posZ_Bio));

figure(2) %plot ratios vs SZA
plot(SZ_PMA,Rat_PMA,'ro','markerfacecolor','y'),grid on,hold on
plot(SZ_Bio,Rat_Bio0616,'k.')
plot(SZ_Bio,Rat_Bio3724,'g.')
plot([0 90],[Calfact_PMA_tmp Calfact_PMA_tmp],'r--')
plot([0 90],[Calfact_Bio0616_tmp Calfact_Bio0616_tmp],'k--')
plot([0 90],[Calfact_Bio3724_tmp Calfact_Bio3724_tmp],'g--')
legend('PMA','Bio0616 Tak','Bio3724 Bakke')
hold off
xlabel('SZA')
ylabel('Calfact UVI/MED/H')
xlim([40 90])
ylim([0 5])

%Lag kurvetilpasning vs SZA og beregn kalfaktorer på nytt
%last opp curvefitting toolboksen og prøv ulike funksjoner

%Bruker curvefitting toolboks til å tilpasse en eksponentiell funksjon
%til x,y data, av formen c-a*exp(b*x)
%Last inn x, y data med j=1;         xData=DATA{j}.median_tau_pyro; yData=DATA{j}.quartile_tau_UVI;
%Skriv i kommando vinduet: cftool(xData,yData). Her kan du velge
%kurvefunksjoner, og velge File, generate function.

x_trinn=40:1:90;
xData=SZ_PMA;p1=find(xData<90);xData=xData(p1);
yData=Rat_PMA;yData=yData(p1);
[fitresult, gof] = skyggeprosjekt_create_exp_Fit(xData, yData);
figure(3)
subplot(2,1,1)
plot(SZ_PMA,Rat_PMA,'ro','markerfacecolor','y'),grid on,hold on
%plot(xData,yData,'go-','markersize',2)
%[fitresult, gof] = create_exp1_Fit(xData, yData );
fit_normalisert=fitresult(x_trinn')/fitresult(40); %normaliserer til SZA 40
SZ_corr=[[0,1];[x_trinn;fit_normalisert']';[180 fit_normalisert(end)]];%gjelder nå SZA 0:180 grader
%plot(x_trinn,fit_normalisert,'b.-')
plot(SZ_PMA,Rat_PMA./interp1(SZ_corr(:,1),SZ_corr(:,2),SZ_PMA,'linear'),'go-','markersize',2)
%Calfact_PMA=median((ones(size(Rat_PMA))./Rat_PMA)./interp1(SZ_corr(:,1),SZ_corr(:,2),SZ_PMA,'linear'))
Calfact_PMA=median(yData./interp1(SZ_corr(:,1),SZ_corr(:,2),xData,'linear'));
hold off
ylabel('UVI/MED/H')
title('ratio reference GUV [UVI]/PMA [MED/H]')
ylim([0 6])

%Korrigerte PMA-målinger i enheter av UVI
subplot(2,1,2)
plot(RES.dvec,RES.CIE_FARIN*40,'r.-'),grid on,grid minor
hold on
%plot(PMA.dvec,PMA.Raw*Calfact_PMA./(interp1(SZ_corr(:,1),SZ_corr(:,2),PMA.SZA,'linear')),'k.-')
plot(PMA.dvec,PMA.Raw*Calfact_PMA.*interp1(SZ_corr(:,1),SZ_corr(:,2),PMA.SZA,'linear'),'k.-')
xlim([dateshift(PMA.dvec(1),'start','day') dateshift(PMA.dvec(1),'end','day')])
legend('Reference GUV','PMA')
ylabel('UVI')
title('Calibrated data at DSA Østerås')
PMA.UVI_calfactor=Calfact_PMA;
PMA.SZA_correction=SZ_corr;
PMA.Calib_info='UVI=PMA.Raw*PMA.UVI_calfactor*interp1(SZ_corr(:,1),SZ_corr(:,2),PMA.SZA)';

%Gjenter for UVBio0616 tak og bakke:
xData=SZ_Bio;p1=find(xData<90);xData=xData(p1);
yData=Rat_Bio0616;yData=yData(p1);
[fitresult, gof] = skyggeprosjekt_create_exp_Fit(xData, yData);
figure(4) %BIO0616
subplot(2,1,1)
plot(SZ_Bio,Rat_Bio0616,'ro','markerfacecolor','y'),grid on,hold on
fit_normalisert=fitresult(x_trinn')/fitresult(40); %normaliserer til SZA 40
SZ_corr=[[0,1];[x_trinn;fit_normalisert']';[180 fit_normalisert(end)]];%gjelder nå SZA 0:180 grader
plot(SZ_Bio,Rat_Bio0616./interp1(SZ_corr(:,1),SZ_corr(:,2),SZ_Bio,'linear'),'go-','markersize',2)
Calfact_Bio0616=median(yData./interp1(SZ_corr(:,1),SZ_corr(:,2),xData,'linear'));
hold off
ylabel('UVI/MED/H')
title('ratio reference GUV [UVI]/BIO0616 (Tak) [MED/H]')
ylim([0 6])
xlim([40 100])

%Korrigerte BIO-målinger i enheter av UVI
subplot(2,1,2)
plot(RES.dvec,RES.CIE_FARIN*40,'r.-'),grid on,grid minor
hold on
plot(BIO.dvec,BIO.Raw0616*Calfact_Bio0616.*interp1(SZ_corr(:,1),SZ_corr(:,2),BIO.SZA,'linear'),'k.-')
xlim([dateshift(BIO.dvec(1),'start','day') dateshift(BIO.dvec(end),'end','day')])
legend('Reference GUV','BIO0616')
ylabel('UVI')
title('Calibrated data at DSA Østerås')
BIO.UVI_calfactor0616=Calfact_Bio0616;
BIO.SZA_correction0616=SZ_corr;
BIO.Calib_info='UVI=BIO.Raw*BIO.UVI_calfactorXXXX*interp1(SZ_corr(:,1),SZ_corr(:,2),BIO.SZA)';


%UVBIO på bakken - 3724
xData=SZ_Bio;p1=find(xData<90);xData=xData(p1);
yData=Rat_Bio3724;yData=yData(p1);
[fitresult, gof] = skyggeprosjekt_create_exp_Fit(xData, yData);
figure(5)%BIO3724
subplot(2,1,1)
plot(SZ_Bio,Rat_Bio3724,'ro','markerfacecolor','y'),grid on,hold on
fit_normalisert=fitresult(x_trinn')/fitresult(40); %normaliserer til SZA 40
SZ_corr=[[0,1];[x_trinn;fit_normalisert']';[180 fit_normalisert(end)]];%gjelder nå SZA 0:180 grader
plot(SZ_Bio,Rat_Bio3724./interp1(SZ_corr(:,1),SZ_corr(:,2),SZ_Bio,'linear'),'go-','markersize',2)
Calfact_Bio3724=median(yData./interp1(SZ_corr(:,1),SZ_corr(:,2),xData,'linear'));
hold off
ylabel('UVI/MED/H')
title('ratio reference GUV [UVI]/BIO3724 (Bakke) [MED/H]')
ylim([0 6])
xlim([40 100])

%Korrigerte BIO-målinger i enheter av UVI
subplot(2,1,2)
plot(RES.dvec,RES.CIE_FARIN*40,'r.-'),grid on,grid minor
hold on
plot(BIO.dvec,BIO.Raw3724*Calfact_Bio3724.*interp1(SZ_corr(:,1),SZ_corr(:,2),BIO.SZA,'linear'),'k.-')
xlim([dateshift(BIO.dvec(1),'start','day') dateshift(BIO.dvec(end),'end','day')])
legend('Reference GUV','BIO3724')
ylabel('UVI')
title('Calibrated data at DSA Østerås')
BIO.UVI_calfactor3724=Calfact_Bio3724;
BIO.SZA_correction3724=SZ_corr;
BIO.Calib_info='UVI=BIO.Raw*BIO.UVI_calfactorXXXX*interp1(SZ_corr(:,1),SZ_corr(:,2),BIO.SZA)';

%Plott UVI for alle instrumenter for hele perioden, delt opp i 4 første og
%4 siste døgn. Inkludert klarværsberegnet:

figure(6)
subplot(2,1,1)
plot(RES.dvec,RES.CIE_FARIN*40,'r.-'),grid on,grid minor
hold on
plot(PMA.dvec,PMA.Raw*PMA.UVI_calfactor.*interp1(PMA.SZA_correction(:,1),PMA.SZA_correction(:,2),PMA.SZA,'linear'),'m.-')
plot(BIO.dvec,BIO.Raw0616*BIO.UVI_calfactor0616.*interp1(BIO.SZA_correction0616(:,1),BIO.SZA_correction0616(:,2),BIO.SZA,'linear'),'b-')
plot(BIO.dvec,BIO.Raw3724*BIO.UVI_calfactor3724.*interp1(BIO.SZA_correction3724(:,1),BIO.SZA_correction3724(:,2),BIO.SZA,'linear'),'g-')
plot(RES.MODEL.dvec,RES.MODEL.UVI,'r--')
hold off
legend('Reference GUV','PMA','BIO0616','BIO3724','Clearsky modelled')
start_dato=datetime(datenum([2021,7,28,0,0,0]), 'ConvertFrom', 'datenum');
slutt_dato=datetime(datenum([2021,7,31,0,0,0]), 'ConvertFrom', 'datenum');
%xlim([dateshift(cal_dato,'start','day') dateshift(cal_dato,'end','day')])
xlim([dateshift(start_dato,'start','day') dateshift(slutt_dato,'end','day')])
ax = gca;
ax.YLim(1)=0;
title('First 4 days of period at Østerås' )
ylabel('UVI')

subplot(2,1,2)
plot(RES.dvec,RES.CIE_FARIN*40,'r.-'),grid on,grid minor
hold on
plot(PMA.dvec,PMA.Raw*PMA.UVI_calfactor.*interp1(PMA.SZA_correction(:,1),PMA.SZA_correction(:,2),PMA.SZA,'linear'),'k-')
plot(BIO.dvec,BIO.Raw0616*BIO.UVI_calfactor0616.*interp1(BIO.SZA_correction0616(:,1),BIO.SZA_correction0616(:,2),BIO.SZA,'linear'),'b-')
plot(BIO.dvec,BIO.Raw3724*BIO.UVI_calfactor3724.*interp1(BIO.SZA_correction3724(:,1),BIO.SZA_correction3724(:,2),BIO.SZA,'linear'),'g-')
plot(RES.MODEL.dvec,RES.MODEL.UVI,'r--')
hold off
legend('Reference GUV','PMA','BIO0616','BIO3724','Clearsky modelled')
% cal_dato=datetime(datenum([2021,8,22,0,0,0]), 'ConvertFrom', 'datenum');
% xlim([dateshift(cal_dato,'start','day') dateshift(cal_dato,'end','day')])
start_dato=datetime(datenum([2021,8,19,0,0,0]), 'ConvertFrom', 'datenum');
slutt_dato=datetime(datenum([2021,8,22,0,0,0]), 'ConvertFrom', 'datenum');
xlim([dateshift(start_dato,'start','day') dateshift(slutt_dato,'end','day')])
ax = gca;
ax.YLim(1)=0;
title('Last 4 days of period at Østerås' )
ylabel('UVI')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Beregn også daglig ratio mellom Bios og RES for å se drift fram til
%målinger i barnehagen på dag 25.08.2021

dager=unique(BIO.Dnum);
dv_PMA=[];DM_PMA_2_Ref=[];DSTD_PMA_2_Ref=[];%rat_PMA=[];dv_rat_PMA=[];
dv_B=[]; DM_BIO0616_2_Ref=[];DM_BIO3724_2_Ref=[];
DSTD_BIO0616_2_Ref=[];DSTD_BIO3724_2_Ref=[];
for k=1:numel(dager)
    [aar,maaned,dag]=daynum2date_guv(year,dager(k));
    ds=datestr([aar maaned dag 12 0 0]);%datestring
    DTnum=datenum(ds);
    DV=datetime(DTnum, 'ConvertFrom', 'datenum');
    
    pB=find(BIO.Dnum==dager(k) & BIO.SZA<=75);
    pR=find(RES.Daynums==dager(k) & RES.SZA<=75);
    pPMA=find(PMA.Dnum==dager(k) & PMA.SZA<=75);
    if ~isempty(pB)
        [~,r,s]=intersect(BIO.dvec(pB),RES.dvec(pR));
        dv_B=[dv_B;DV];
        r0616=BIO.Raw0616(pB(r))*BIO.UVI_calfactor0616.*interp1(BIO.SZA_correction0616(:,1),BIO.SZA_correction0616(:,2),BIO.SZA(pB(r)),'linear')./(RES.CIE_FARIN(pR(s))*40);
        r3724=BIO.Raw3724(pB(r))*BIO.UVI_calfactor3724.*interp1(BIO.SZA_correction3724(:,1),BIO.SZA_correction3724(:,2),BIO.SZA(pB(r)),'linear')./(RES.CIE_FARIN(pR(s))*40);
        DM_BIO0616_2_Ref=[DM_BIO0616_2_Ref;mean(r0616)];
        DM_BIO3724_2_Ref=[DM_BIO3724_2_Ref;mean(r3724)];
        DSTD_BIO0616_2_Ref=[DSTD_BIO0616_2_Ref;std(r0616)];DSTD_BIO3724_2_Ref=[DSTD_BIO3724_2_Ref;std(r3724)];
        %daily_mean_ratios_B=[daily_calfacts_B;[1/mean(r0616) 1/mean(r3724)] ];
    end
    
    if ~isempty(pPMA)
        dv_PMA=[dv_PMA;DV];
        [~,t,u]=intersect(PMA.dvec(pPMA),RES.dvec(pR));
        %rPMA=PMA.Raw(pPMA(t))./(RES.CIE_FARIN(pR(u))*40);
        rPMA=PMA.Raw(pPMA(t))*PMA.UVI_calfactor.*interp1(PMA.SZA_correction(:,1),PMA.SZA_correction(:,2),PMA.SZA(pPMA(t)),'linear')./(RES.CIE_FARIN(pR(u))*40);
        DM_PMA_2_Ref=[DM_PMA_2_Ref;mean(rPMA)];
        DSTD_PMA_2_Ref=[DSTD_PMA_2_Ref;std(rPMA)];
        %         dv_rat_PMA=[dv_rat_PMA;PMA.dvec(pPMA(t))];
        %         rat_PMA=[rat_PMA;ones(size(rPMA))./rPMA];
    end
    
end

%plott faktorer:
figure(7)
errorbar(dv_PMA,DM_PMA_2_Ref,DSTD_PMA_2_Ref,'mo','markersize',10,'markerfacecolor','m'),grid on,grid minor,hold on
errorbar(dv_B,DM_BIO0616_2_Ref,DSTD_BIO0616_2_Ref,'bs','markerfacecolor','b','markersize',10)
errorbar(dv_B,DM_BIO3724_2_Ref,DSTD_BIO0616_2_Ref,'go','markerfacecolor','g','markersize',10)
xmin=dv_B(1);
xmax=dv_B(end);
%plot([xmin xmax],[Calfact_PMA Calfact_PMA],'y--')
% plot([xmin xmax],[Calfact_Bio0616 Calfact_Bio0616],'g--')
% plot([xmin xmax],[Calfact_Bio3724 Calfact_Bio3724],'m--')
% plot(dv_rat_PMA,rat_PMA,'y.')
hold off
legend('PMA','Bio0616 Tak','Bio3724 Bakke')
title('Daily mean ratios to reference GUV, SZA<75')
xlabel('Date')
ylabel('Rel. units')


%Lagre PMA og BIO data
PMA.Station=station;
tmp='PMA';
mat_file=sprintf('L:\\Optisk Lab\\Uvnet\\Prosjekter\\Skyggeprosjekt\\%04d\\Kalibreringsdata_Oesteras\\PMA_05855_%s_%04d.mat',year,station,year);
save(mat_file,tmp);

BIO.Station=station;
tmp='BIO';
mat_file=sprintf('L:\\Optisk Lab\\Uvnet\\Prosjekter\\Skyggeprosjekt\\%04d\\Kalibreringsdata_Oesteras\\UV_Bios_0616_3724_%s_%04d.mat',year,station,year);
save(mat_file,tmp);

