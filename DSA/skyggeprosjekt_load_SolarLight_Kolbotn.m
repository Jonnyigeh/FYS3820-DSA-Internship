function skyggeprosjekt_load_SolarLight_Kolbotn
%function skyggeprosjekt_load_PMA_UVBios
%laster opp PMA-loggerdata og UVBiodata fra målinger i barnehagen
clean
%https://www.google.no/maps/@59.8114997,10.8031707,17z?hl=no
station='Kolbotn Skolebakken';
year=2021;
lat=59.8114997;
long=10.8031707;

%Last først opp structs med kalibreringsdata for solar light instrumentene

load('L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\GUV_oesteraas_2021.mat'); %RES GUV
load('L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\PMA_05855_oesteraas_2021.mat'); %PMA
load('L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\UV_Bios_0616_3724_oesteraas_2021.mat'); %BIO

PMAUVB = Skyggeprosjekt_import_PMA_file("L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Målinger\PMA_UVB_25082021.log", [2, Inf]);
PMAUVB(1,:)

[Y,M,D] = ymd(PMAUVB.loggerdato);
[h,m,s] = hms(PMAUVB.loggertid);
ds=datestr([Y M D h m s]);

DTnum_PMA=datenum(ds);
p=find(unique(DTnum_PMA));
dvec=datetime(DTnum_PMA, 'ConvertFrom', 'datenum');

BHG.PMA.dvec=dvec;
BHG.PMA.Y=Y;
BHG.PMA.M=M;
BHG.PMA.D=D;
BHG.PMA.h=h;
BHG.PMA.m=m;
BHG.PMA.s=s;
BHG.PMA.Raw=PMAUVB.mleverdi;
BHG.PMA.Dnum=day(datetime([Y';M';D']'),'dayofyear');

clear PMAUVB

%Beregn SZA og AZI ut fra UTC-verdier:
BHG.PMA.SZA=zeros(size(BHG.PMA.h)); % zenith vinkel
BHG.PMA.AZI=zeros(size(BHG.PMA.h)); % zenith vinkel
for q=1:numel(BHG.PMA.h)
    [azi,zen]=solpos_bj(BHG.PMA.Y(q),BHG.PMA.M(q),BHG.PMA.D(q),BHG.PMA.h(q)+1,BHG.PMA.m(q),BHG.PMA.s(q),lat,long);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
    BHG.PMA.SZA(q)=zen*180/pi;
    BHG.PMA.AZI(q)=azi*180/pi;
end

BHG.PMA.UVI=BHG.PMA.Raw*PMA.UVI_calfactor.*interp1(PMA.SZA_correction(:,1),PMA.SZA_correction(:,2),BHG.PMA.SZA,'linear');

figure(1)
plot(BHG.PMA.dvec,BHG.PMA.UVI,'k.-'),grid on,grid minor
xlabel('Tid, UTC')
ylabel('UVI')
title('Kolbotn barnehage, avd. Solbakken')
ax = gca;
ax.YLim(1)=0;

% subplot(2,1,2)
% bar(BHG.PMA.dvec,BHG.PMA.UVI,'k'),grid on


%Lst opp UV-Bios
UVBios = Skyggeprosjekt_import_UVBios_file("L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Målinger\UVBios_Kolbotn_24-25_august_2021.log", [2, Inf]);
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

BHG.BIO.dvec=dvecUVBios;
BHG.BIO.Raw0616=UVBios.AVGBio0616;
BHG.BIO.Raw3724=UVBios.AVGFin3724;
BHG.BIO.Y=YB;
BHG.BIO.M=MB;
BHG.BIO.D=DB;
BHG.BIO.h=hB;
BHG.BIO.m=mB;
BHG.BIO.s=sB;
BHG.BIO.Dnum=day(datetime([YB';MB';DB']'),'dayofyear');

%Beregn SZA og AZI ut fra UTC-verdier:
BHG.BIO.SZA=zeros(size(BHG.BIO.h)); % zenith vinkel
BHG.BIO.AZI=zeros(size(BHG.BIO.h)); % zenith vinkel
for q=1:numel(BHG.BIO.h)
    [azi,zen]=solpos_bj(BHG.BIO.Y(q),BHG.BIO.M(q),BHG.BIO.D(q),BHG.BIO.h(q)+1,BHG.BIO.m(q),BHG.BIO.s(q),lat,long);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
    BHG.BIO.SZA(q)=zen*180/pi;
    BHG.BIO.AZI(q)=azi*180/pi;
end

clear UVBios

BHG.BIO.UVI0616=BHG.BIO.Raw0616*BIO.UVI_calfactor0616.*interp1(BIO.SZA_correction0616(:,1),BIO.SZA_correction0616(:,2),BHG.BIO.SZA,'linear');
BHG.BIO.UVI3724=BHG.BIO.Raw3724*BIO.UVI_calfactor3724.*interp1(BIO.SZA_correction3724(:,1),BIO.SZA_correction3724(:,2),BHG.BIO.SZA,'linear');

figure(2)
subplot(2,1,1)
plot(BHG.BIO.dvec,BHG.BIO.UVI0616,'k.-'),grid on,grid minor,hold on
plot(BHG.BIO.dvec,BHG.BIO.UVI3724,'b.-'), hold off
xlabel('Tid, UTC')
ylabel('UVI')
legend('0616 Tak','3724 Bakke')
title([{'Solar ligt metre plassert på taket og på skråtak / bakken'},{'Kolbotn barnehage'}])

subplot(2,1,2)
plot(BHG.BIO.dvec,BHG.BIO.UVI3724./BHG.BIO.UVI0616,'b.-'), hold off,grid on,grid minor
ylabel('Ratio Bakke/tak')

%Vis PMA og UVBios for samme tidspunkt
%[~,ia,ib]=intersect(BHG.PMA.dvec,BHG.BIO.dvec);
pos_1=find(BHG.BIO.dvec>=datetime(datenum([2021,8,24,1,0,0]), 'ConvertFrom', 'datenum'));
pos_2=find(BHG.BIO.dvec>=BHG.PMA.dvec(1));
pos_3=find(RES.MODEL.dvec>=datetime(datenum([2021,8,24,1,0,0]), 'ConvertFrom', 'datenum') & ...
    RES.MODEL.dvec<=datetime(datenum([2021,8,26,1,0,0]), 'ConvertFrom', 'datenum'));

figure(3)
plot(BHG.PMA.dvec,BHG.PMA.UVI,'ro-'),grid on
xlabel('Tid, UTC')
title('Kolbotn barnehage, avd. Solbakken')
hold on
plot(BHG.BIO.dvec(pos_2),BHG.BIO.UVI0616(pos_2),'k.-'),grid on,grid minor,hold on
plot(BHG.BIO.dvec(pos_2),BHG.BIO.UVI3724(pos_2),'b.-'), hold off
xlabel('Tid, UTC=Lokaltid-2H')
ylabel('UVI')
legend('PMA','Tak 0616','Bakke 3724')

figure(4) %Bare Bio 0616 og Bio 3724. Legg til klarværsmodellert
plot(BHG.BIO.dvec(pos_1),BHG.BIO.UVI0616(pos_1),'k.-'),grid on,grid minor,hold on
plot(BHG.BIO.dvec(pos_1),BHG.BIO.UVI3724(pos_1),'b.-'), 
plot(RES.MODEL.dvec(pos_3),RES.MODEL.UVI(pos_3),'r-')
hold off
xlabel('Tid, UTC')
title('Kolbotn barnehage, avd. Solbakken')
xlabel('Tid, UTC=Lokaltid-2H')
ylabel('UVI')
legend('Tak 0616','Bakke 3724','Klarværsmodellert')

%Lagre data til mat-fil, og som tekstfiler
BHG.Station=station;
BHG.Year=year;
BHG.Lat=lat;
BHG.Long=long;
tmp='BHG';
mat_file='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Målinger\Barnehage_maalinger.mat';
save(mat_file,tmp);

%lagre UVI data
%save2txtfile=0;
save2txtfile=1;
if save2txtfile
    fname_meas='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Målinger\PMA_UVI_Kolbotn_2021.txt';
    [filepath,name,ext] = fileparts(fname_meas);
    tmp_name=sprintf('C:\\Users\\bjohnsen\\Downloads\\%s%s',name,ext);    %lagrer midlertidig lokalt på downloads
    
    fid=fopen(tmp_name,'wt');
    headtekst=sprintf('%% UVI beregnet for PMA %s - %d, Tid=UTC+0',station,year);
    headtekst_2=sprintf('%% Date Time UVI SZA');
    fprintf(fid,'%s\n',headtekst);
    fprintf(fid,'%s\n',headtekst_2);
    pos_m=find(BHG.PMA.dvec>=datetime(datenum([2021,8,24,0,0,0]), 'ConvertFrom', 'datenum'));
    
    %Regner om dvec til datestring
    ds=datetime(datenum(BHG.PMA.dvec(pos_m)),'ConvertFrom', 'datenum', 'Format', 'ddMMyyyy hh:mm');
    
    for row=1:numel(pos_m)
        %fprintf(fid,'%08s %02s:%02s %6.4f %6.3f\n',RES.DateHour(pos_m(row),1:8),RES.DateHour(pos_m(row),10:11),RES.DateHour(pos_m(row),13:14),RES.CIE_FARIN(pos_m(row))*40,RES.SZA(pos_m(row)));
        fprintf(fid,'%s %6.4f %6.3f\n',ds(row),BHG.PMA.UVI(pos_m(row)),BHG.PMA.SZA(pos_m(row)));
    end
    fclose(fid)
    
    status = movefile(tmp_name,filepath);%flytter fila over på fellesområdet:
    
    %Gjenta for BIO
    fname_meas='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Målinger\UVBios_UVI_Kolbotn_2021.txt';
    [filepath,name,ext] = fileparts(fname_meas);
    tmp_name=sprintf('C:\\Users\\bjohnsen\\Downloads\\%s%s',name,ext);    %lagrer midlertidig lokalt på downloads
    
    fid=fopen(tmp_name,'wt');
    headtekst=sprintf('%% UVI beregnet for UV-Biometre %s - %d, Tid=UTC+0',station,year);
    headtekst_2=sprintf('%% Date Time UVI_0616 UVI_3724 SZA');
    fprintf(fid,'%s\n',headtekst);
    fprintf(fid,'%s\n',headtekst_2);
    pos_m=find(BHG.BIO.dvec>=datetime(datenum([2021,8,24,0,0,0]), 'ConvertFrom', 'datenum'));
    
    %Regner om dvec til datestring
    ds=datetime(datenum(BHG.BIO.dvec(pos_m)),'ConvertFrom', 'datenum', 'Format', 'ddMMyyyy hh:mm');
    
    for row=1:numel(pos_m)
        %fprintf(fid,'%08s %02s:%02s %6.4f %6.3f\n',RES.DateHour(pos_m(row),1:8),RES.DateHour(pos_m(row),10:11),RES.DateHour(pos_m(row),13:14),RES.CIE_FARIN(pos_m(row))*40,RES.SZA(pos_m(row)));
        fprintf(fid,'%s %6.4f %6.4f %6.3f\n',ds(row),BHG.BIO.UVI0616(pos_m(row)),BHG.BIO.UVI3724(pos_m(row)),BHG.BIO.SZA(pos_m(row)));
    end
    fclose(fid)
    
    status = movefile(tmp_name,filepath);%flytter fila over på fellesområdet:
    

    
%     fname_modelled='L:\Optisk Lab\Uvnet\Prosjekter\Skyggeprosjekt\2021\Kalibreringsdata_Oesteras\klarvaersdata_Oesteraas juli-august_2021.txt';
%     [filepath,name,ext] = fileparts(fname_modelled);
%     tmp_name=sprintf('C:\\Users\\bjohnsen\\Downloads\\%s%s',name,ext);    %lagrer midlertidig lokalt på downloads
%     fid=fopen(tmp_name,'wt');
%     headtekst=sprintf('%% Klarværsmodellerte data %s - %d, Tid=UTC+0',station,year);
%     headtekst_2=sprintf('%% Date Time UVI SZA');
%     fprintf(fid,'%s\n',headtekst);
%     fprintf(fid,'%s\n',headtekst_2);
%     pos_modelled=find(RES.MODEL.dvec>=datetime(datenum([2021,7,25,0,0,0]), 'ConvertFrom', 'datenum') & ...
%         RES.MODEL.dvec<=datetime(datenum([2021,8,25,0,0,0]), 'ConvertFrom', 'datenum'));
%     
%     %Regner om dvec til datestring
%     ds=datetime(datenum(RES.MODEL.dvec(pos_modelled)),'ConvertFrom', 'datenum', 'Format', 'ddMMyyyy hh:mm');
%     
%     for row=1:numel(pos_modelled)
%         fprintf(fid,'%s %6.4f %6.3f\n',ds(row),RES.MODEL.UVI(pos_modelled(row)),RES.MODEL.SZA(pos_modelled(row)));
%     end
%     fclose(fid)
%     
%     status = movefile(tmp_name,filepath);%flytter fila over på fellesområdet:
    
end


