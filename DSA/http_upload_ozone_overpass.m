function [OMI,Status]=http_upload_ozone_overpass(station,time_out_problems)
%laster ned fra http hvis det ikke er timeout problemer med serveren.
%Alternativt, leses overpass data som allerede er lastet ned, men som er
%gamle

% [OMI,Status]=http_upload_ozone_overpass('oesteraas',1)
% [OMI,Status]=http_upload_ozone_overpass('blindern',1)
% [OMI,Status]=http_upload_ozone_overpass('tromso',1)
% [OMI,Status]=http_upload_ozone_overpass('andoya',1)
% [OMI,Status]=http_upload_ozone_overpass('nyaal',1)
% [OMI,Status]=http_upload_ozone_overpass('kise',1)
% [OMI,Status]=http_upload_ozone_overpass('trondheim',1)
% [OMI,Status]=http_upload_ozone_overpass('bergen',1)
% [OMI,Status]=http_upload_ozone_overpass('finse',1)
% [OMI,Status]=http_upload_ozone_overpass('landvik',1)

path_overpass='N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\Satelitt_overpass_time';
%path_overpass='D:\work\uvnett\satelitt_overpass_time';

switch station
    case 'oesteraas'
        omi_file='aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt
    case 'Oesteraas'
        omi_file='aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt
    case 'Østerås'
        omi_file='aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt        
    case 'blindern'
        omi_file='aura_omi_l2ovp_omuvb_v03_oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_oslo.txt
    case 'kjeller'
        omi_file='aura_omi_l2ovp_omuvb_v03_oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_oslo.txt
    case 'Blindern'
        omi_file='aura_omi_l2ovp_omuvb_v03_oslo.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_oslo.txt
    case 'tromso'
        omi_file='aura_omi_l2ovp_omuvb_v03_tromso.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_tromso.txt
    case 'Tromso'
        omi_file='aura_omi_l2ovp_omuvb_v03_tromso.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_tromso.txt
    case 'Tromsø'
        omi_file='aura_omi_l2ovp_omuvb_v03_tromso.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_tromso.txt        
    case 'andoya'
        omi_file='aura_omi_l2ovp_omuvb_v03_alomar.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_alomar.txt
        %omi_file='aura_omi_l2ovp_omaeruv_v03_andenes.txt';
    case 'Andoya'
        omi_file='aura_omi_l2ovp_omuvb_v03_alomar.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_alomar.txt
        %omi_file='aura_omi_l2ovp_omaeruv_v03_andenes.txt';
    case 'Andøya'
        omi_file='aura_omi_l2ovp_omuvb_v03_alomar.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_alomar.txt
        %omi_file='aura_omi_l2ovp_omaeruv_v03_andenes.txt';        
    case 'nyaal'
        omi_file='aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt
    case 'Nyaalesund'
        omi_file='aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt
    case 'trondheim'
        omi_file='aura_omi_l2ovp_omuvb_v03_trondheim.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_trondheim.txt
    case 'Trondheim'
        omi_file='aura_omi_l2ovp_omuvb_v03_trondheim.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_trondheim.txt
    case 'kise'
        omi_file='aura_omi_l2ovp_omuvb_v03_kise.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_kise.txt
    case 'Kise'
        omi_file='aura_omi_l2ovp_omuvb_v03_kise.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_kise.txt
    case 'bergen'
        omi_file='aura_omi_l2ovp_omuvb_v03_bergen.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_bergen.txt
    case 'Bergen'
        omi_file='aura_omi_l2ovp_omuvb_v03_bergen.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_bergen.txt
    case 'finse'
        omi_file='aura_omi_l2ovp_omuvb_v03_finse.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_finse.txt
    case 'Finse'
        omi_file='aura_omi_l2ovp_omuvb_v03_finse.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_finse.txt
    case 'landvik'
        omi_file='aura_omi_l2ovp_omuvb_v03_landvik.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_landvik.txt
    case 'Landvik'
        omi_file='aura_omi_l2ovp_omuvb_v03_landvik.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_landvik.txt
    case 'ekofisk'
        omi_file='aura_omi_l2ovp_omuvb_v03_landvik.txt';
    case 'Dar-es-salaam'
        omi_file='aura_omi_l2ovp_omuvb_v03_dar.es.salaam.txt';
        %http://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_dar.es.salaam.txt
    case 'Jokioinen'
        omi_file='aura_omi_l2ovp_omuvb_v03_jokioinen.txt';
        %https://avdc.gsfc.nasa.gov/download_2.php?site=2057856112&id=79&go=download&path=&file=aura_omi_l2ovp_omuvb_v03_jokioinen.txt
    otherwise
        return
end


%untar(tarfilename,outputdir)

fname=sprintf('%s\\%s',path_overpass,omi_file);
%%%%%%%%%%%%%%%%%%%%%
% Direkte nedlasting av fil på csv-format:
%web('https://veret.gfi.uib.no/taarn/view.php?f=2021-08-01&t=Minutt&csv=1')
url=[('https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/'),(sprintf('%s',omi_file))];
%web(url);
outfilename = websave(fname,url);%laster ned en fil fra internet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Kjør manuelt disse wgetkommandoene i cmd-vinduet. filer havne i downloads. Flytt derfra til N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\Satelitt_overpass_time
%%wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_andoya.txt -P /path/to/folder C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_oslo.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_tromso.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_alomar.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_trondheim.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_kise.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_bergen.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_finse.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_landvik.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_dar.es.salaam.txt -P C:/Users/bjohnsen/Downloads
% % wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_jokioinen.txt -P C:/Users/bjohnsen/Downloads
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Last ned data fra alle stasjonene manuelt først: Kjøre linje 130-153. deretter kommenteres linjene
% %ut, slik at http_upload_ozone_overpass.m kan kjøres som et kall
% %Last inn overpassdata fra uvnettverksstasjonene på en gang:
% name_list={'aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt',...
%     'aura_omi_l2ovp_omuvb_v03_oslo.txt',...
%     'aura_omi_l2ovp_omuvb_v03_tromso.txt',...
%     'aura_omi_l2ovp_omuvb_v03_alomar.txt',...
%     'aura_omi_l2ovp_omuvb_v03_ny.alesund.2.txt',...
%     'aura_omi_l2ovp_omuvb_v03_trondheim.txt',...
%     'aura_omi_l2ovp_omuvb_v03_kise.txt',...
%     'aura_omi_l2ovp_omuvb_v03_bergen.txt',...
%     'aura_omi_l2ovp_omuvb_v03_finse.txt',...
%     'aura_omi_l2ovp_omuvb_v03_landvik.txt'};
% 
% wd=pwd;
% path_wget='c:\users\bjohnsen';
% cd(path_wget);
% %b=sprintf('uvspec <clear_aerosols_290_500_parprocessing.inp> syntspec.out');
% %b=sprintf('uvspec <%s> syntspec_%d.out',tempname,j);
% for j=1:size(name_list,2)
%     b=['wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/',sprintf('%s',name_list{j}),' -P C:/Users/bjohnsen/Downloads/ozone_overpass'];
%     %b=('wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_osteras.oslo.txt -P C:/Users/bjohnsen/Downloads/ozone_overpass');5
%     dos(b)%samme som !og resten av uvspec kommandoen
% end
% 
% status = movefile('C:/Users/bjohnsen/Downloads/ozone_overpass/*',path_overpass);
% cd(wd)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% if time_out_problems
%
% else
% if ~time_out_problems
%     nytter ikke å laste ned ozon-data automatisk - da det kjøres et skript på serveren.
%     Last derfor ned OMI overpass data manuelt på forhånd.
%
%     %https://gs614-avdc1-pz.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/
%     [f, status] = urlwrite(sprintf('http://avdc.gsfc.nasa.gov/download_2.php?site=595385375&id=79&go=download&path=&file=%s',omi_file),sprintf('%s',fname));
%     [f, status] = urlwrite(sprintf('http://avdc.gsfc.nasa.gov/index.php?site=2057856112&id=79&go=download&path=&file=%s',omi_file),sprintf('%s',fname));
%
%     https://gs614-avdc1-pz.gsfc.nasa.gov/index.php?site=1560054543&id=85&go=list&path=/txt
%     Bruk heller websave: F.eks. nedlasting av data fra github
%     url_tst='https://raw.githubusercontent.com/uvnrpa/Minute_Data/master/AND/UVI_Andoya_2000.txt';
%     fname='E:\temp_download\uvi_andoya_2001.txt';
%     outfilename = websave(fname,url_tst)
%     url_folder='https://avdc.gsfc.nasa.gov/index.php?site=2057856112&id=79';
%
%     mypath=pwd;
%     cd c:\users\bjohnsen
%     command_line=sprintf('wget -c -i --no-check-certificate %s/%s',url_folder,omi_file)
%     system(command_line)
%     cd(mypath)
%
%     url='https://avdc.gsfc.nasa.gov/index.php?site=2057856112&id=79/aura_omi_l2ovp_omuvb_v03_kise.txt';
%
% %     api = url_folder;%'http://www.ngdc.noaa.gov/stp/space-weather/';
% %     url = [api 'solar-data/solar-indices/sunspot-numbers/' ...
% %         'american/lists/list_aavso-arssn_yearly.txt'];
%     filename = 'OMI_ozon_kise.txt';
%     options = weboptions('Timeout',Inf);
%     outfilename = websave(filename,url,options)
%
%
%
%     wget -nd -r --header='Authorization: Bearer CF486EFA-B84E-11E9-9166-9BCBF9A6CF41'
% https://avdc.gsfc.nasa.gov/index.php?site=2057856112&id=79/aura_omi_l2ovp_omuvb_v03_finse.txt
%
% wget -c -i https://avdc.gsfc.nasa.gov/index.php?site=2057856112&id=79/aura_omi_l2ovp_omuvb_v03_finse.txt
% wget -c -i https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_andoya.txt
% wget -nd  https://avdc.gsfc.nasa.gov/pub/data/satellite/Aura/OMI/V03/L2OVP/OMUVB/aura_omi_l2ovp_omuvb_v03_andoya.txt
%
%
%     %%url = 'https://gs614-avdc1-pz.gsfc.nasa.gov/index.php?site=1560054543&id=85&go=list&path=/txt';
%     filename = sprintf('%s\\%s',path_overpass,omi_file);
%     options = weboptions('Timeout',Inf);
%     outfilename = websave(filename,url,options);
%
%     S = webread(url)
%
%     % else
%     %     status=1;
%end

stat_=1;
Status=stat_;
if stat_  %fila er lastet ned
    [~,~,year,daynum,sec,~,~,lat,lon,~,sza,~,~,~,~,~,o3,~,~,~,~,~,~,CldOpt,~,~,~,~,~,~,~,~,~,~,~,~,~,uvi,LER,SurfAlb,~] = textread(sprintf('%s\\%s',path_overpass,omi_file),'%s%f%d%d%d%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',51);
    %dropp rader med ekstremverdier av ozon
    ph=find(o3>600);
    if ~isempty(ph)
        year(ph)=[];
        daynum(ph)=[];
        sec(ph)=[];
        lat(ph)=[];
        lon(ph)=[];
        sza(ph)=[];
        o3(ph)=[];
        CldOpt(ph)=[];
        uvi(ph)=[];
        LER(ph)=[];
        SurfAlb(ph)=[];
    end
    OMI.Year=year;
    OMI.Daynum=daynum;
    OMI.Seconds=sec;
    OMI.Lat=lat;
    OMI.Long=lon;
    OMI.SZA=sza;
    OMI.O3=o3;
    OMI.UVI=uvi;
    OMI.CldOpt=CldOpt;
    OMI.LER=LER;
    OMI.SurfAlb=SurfAlb;
    OMI.JT=[];
    
    skuddaar=[1900:4:2100];
    allyears=unique(OMI.Year);
    
    %     daysyear=365*ones(size(allyears));
    %     decyear=[];
    %     Yyear=[];
    %     Dyear=[];   %dagnummervektor
    %     O3year=[];   %medianverdi for dagen
    %     decday=[];   %desimal dagnummervektor
    OMI.Yyear=[];
    OMI.Dyear=[];
    OMI.JTyear=[];
    OMI.O3year=[];
    OMI.CldOptyear=[];
    OMI.LERyear=[];
    OMI.SurfAlbyear=[];
    for i=1:length(allyears)
        if ismember(allyears(i),skuddaar)
            daysyear=366;
        else
            daysyear=365;
        end
        
        py=find(OMI.Year==allyears(i));
        OMI.JT=[OMI.JT;OMI.Year(py)+ (OMI.Daynum(py)+ OMI.Seconds(py)/(24*3600))/daysyear];
        punday=unique(OMI.Daynum(py));
        %beregn median ozon for hver dag
        for k=1:length(punday)
            pm=find( OMI.Daynum(py) == punday(k));
            OMI.Yyear=[OMI.Yyear;median(OMI.Year(py(pm)))];
            OMI.Dyear=[OMI.Dyear;median(OMI.Daynum(py(pm)))];
            
            p_ok=find(OMI.O3(py(pm))>200 & OMI.O3(py(pm))<600);
            if ~isempty(p_ok)
                OMI.O3year=[OMI.O3year;median(OMI.O3(py(pm(p_ok))))];
            else
                OMI.O3year=[OMI.O3year;0];
            end
            
            %OMI.O3year=[OMI.O3year;median(OMI.O3(py(pm)))];
            OMI.CldOptyear=[OMI.CldOptyear;median(OMI.CldOpt(py(pm)))];
            OMI.LERyear=[OMI.LERyear;median(OMI.LER(py(pm)))];
            OMI.SurfAlbyear=[OMI.SurfAlbyear;median(OMI.SurfAlb(py(pm)))];
            OMI.JTyear=[OMI.JTyear;median(OMI.Year(py(pm)))+ median(OMI.Daynum(py(pm)))/daysyear];
            [allyears(i) punday(k)]
        end
        
    end
    
    %     %interpoler ozon til en medianverdi for hver dag
    %     skuddaar=[1900:4:2100];
    % allyears=unique(SAT.Year);
    % daysyear=365*ones(size(allyears));
    % decyear=[];
    % Yyear=[];
    % Dyear=[];   %dagnummervektor
    % decday=[];   %desimal dagnummervektor
    % for i=1:length(allyears)
    %     if ismember(allyears(i),skuddaar)
    %         daysyear(i)=366;
    %     end
    %     %    Xyear=[Xyear;allyears(i)*ones(size(daysyear(i))) + (/daysyear(i))];
    %     Yyear=[Yyear;allyears(i)*ones(daysyear(i),1)];  %vektor med 365 og 366 årstall etter hverandre
    %     Dyear=[Dyear;(1:daysyear(i))'];
    %     decday=[decday;((0.5+(1:daysyear(i)))/daysyear(i))'];
    %
    %     p=find(SAT.Year==allyears(i));
    %     decyear=[decyear;SAT.Year(p)+(SAT.Daynum(p)+SAT.Seconds(p)/(24*3600))/daysyear(i)];
    % end
    %
    % SAT.DecYear=decyear;
    % SAT.Yyear=Yyear;
    % SAT.Dyear=Dyear;
    % SAT.decday=decday;
    % SAT.O3year=interp1(SAT.DecYear,SAT.O3,SAT.Yyear+SAT.decday,'linear');%interpoler ozon til en verdi for hver dag i hele perioden
    
    %     figure
    %     subplot(3,1,1)
    %     plot(OMI.JT,OMI.O3,'k.-'),grid on
    %     subplot(3,1,2)
    %     plot(OMI.JT,OMI.CldOpt,'k.-'),grid on
    %     subplot(3,1,3)
    %     plot(OMI.JT,OMI.SurfAlb,'k.-'),grid on
    %
    %     figure
    %     p2010=find(OMI.Year==2010);
    %     plot(OMI.Daynum(p2010)+OMI.Seconds(p2010)/(24*3600),OMI.O3(p2010),'k.-'),grid on
else
    disp('Stasjonsdata mangler')
end
