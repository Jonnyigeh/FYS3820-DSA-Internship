function RES=database_get_one_year_ozone(RES,OMI)
%Beregner daglig ozon fra daily gridded data, supplert med OMI overpass data
%Beregner RES.MODEL.O3 for modellering av klarværsdoser


pc='dolly';
stasjons_id=RES.Stasjons_id;
station=RES.Stasjon;
lat=RES.Latitude;
long=RES.Longitude;

switch stasjons_id
    case 1
        stasjon='Østerås';
        stat='OST';
    case 2
        stasjon='Kise';
        stat='KIS';
    case 3
        stasjon='Bergen';
        stat='BRG';
    case 4
        stasjon='Landvik';
        stat='LAN';
    case 5
        stasjon='Blindern';
        stat='BLI';
    case 6
        stasjon='Trondheim';
        stat='TRH';
    case 7
        stasjon='Tromsø';
        stat='TSO';
    case 8
        stasjon='Nyaal';
        stat='NYA';
    case 9
        stasjon='Andøya';
        stat='AND';
    case 10
        stasjon='Finse';
        stat='FIN';
    case 11
        stasjon='Kjeller';
        stat='KJE';
    otherwise
        return
end


% station=RES.Stasjon;
% switch stasjons_id
%     case 8
%         station='nyaal';
%     case 7
%         station='tromso';
%     case 9
%         station='andoya';
%     case 6
%         station='trondheim';
%     case 3
%         station='bergen';
%     case 10
%         station='finse';
%     case 2
%         station='kise';
%     case 5
%         station='blindern';
%     case 1
%         station='oesteraas';
%     case 11
%         station='kjeller';
%     case 4
%         station='landvik';
%     case 12
%         station='ekofisk';
%     case 100
%         station='Dar-es-salaam';
%     case 110
%         station='Jokioinen';
%     otherwise
%         return
% end
% 
% switch station
%     case 'tromso'
%         lat=69.65;
%         long=18.93;
%         stat='tso';
%         %         year=
%         %         daynum=[];
%     case 'andoya'
%         lat=69.28;
%         long=16.01;
%         stat='and';
%         %         year=
%         %         daynum=[];
%     case 'nyaal'
%         lat=78.92;
%         long=11.92;
%         stat='nya';
%         %         year=
%         %         daynum=[];
%     case 'trondheim'
%         lat=63.42;
%         long=10.40;
%         stat='trh';
%         %         year=
%         %         daynum=[];
%     case 'kise'
%         lat=60.78;
%         long=10.82;
%         stat='kis';
%         %         year=
%         %         daynum=[];
%     case 'bergen'
%         lat=60.38;
%         long=5.33;
%         stat='brg';
%         %         year=
%         %         daynum=[];
%     case 'finse'
%         lat=60.593;%koordinater fra google https://www.google.com/maps/d/viewer?mid=1upEzn-I3GlSP6zgs4pnZMuQgcNA&hl=en&ll=60.59299437534447%2C7.523682013244638&z=18
%         long=7.524;
%         stat='fin';
%         
%     case 'oesteraas'
%         lat=59.946;
%         long=10.598;
%         stat='ost';
%         %         year=
%         %         daynum=[];
%     case 'blindern'
%         lat=59.938;
%         long=10.719;
%         stat='bli';
%         %         year=
%         %         daynum=[];
%     case 'landvik'
%         lat=58.33;
%         long=8.52;
%         stat='lan';
%         %         year=
%         %         daynum=[];
%     case 'kjeller'
%         lat=59.975;
%         long=11.053;
%         stat='kje';
%     case 'ekofisk'
%         lat=-6-(32+49.991/60)/60;   %hotellplattformen heter 2/4 H på Ekofisk
%         long=3+(12+47.285/60)/60;
%         stat='EKO';
%         
%     case 'Dar-es-salaam';
%         lat=-6.870;
%         long=39.200;
%         stat='DAR';
%         
%     case 'Jokioinen'
%         lat=60.814;
%         long=23.499;
%         stat='JOK';
%         
%     otherwise
%         return
% end
% 
% station





year=RES.Year;
ozone=zeros(366,1);
days=(1:366)';

if year<=2004
    switch pc
        case 'dolly'
            path=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Ozone\\ozone_toms_earthprobe\\Y%04i\\',year);
        case 'lokal'
            path=sprintf('c:\\ozone_toms_earthprobe\\Y%04i\\',year);
        otherwise
            return
    end
else
    switch pc
        case 'dolly'
            path=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\ozone\\ozone_toms_omi\\Y%04i\\',year);
            path_overpass='N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\Satelitt_overpass_time';
        case 'lokal'
            path=sprintf('c:\\ozone\\ozone_toms_omi\\Y%04i\\',year);
        otherwise
            return
    end
    
end

[pf,nf,extf]=fileparts(path);
%N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\ozone_toms_omi\Y2009

upload_gridded_data=1;
if upload_gridded_data
    file_content=dir(pf);
    eksisterer=0;
    for i=3:size(file_content,1)    %catalog 1 og 2 er . og ..
        filename=sprintf('%s\\%s',pf,file_content(i).name);
        aar=str2num(file_content(i).name(14:17));
        maaned=str2num(file_content(i).name(18:19));
        dag=str2num(file_content(i).name(20:21));
        dnum=dagnr(aar,maaned,dag);
        wd=cd;  %matlab directoriet
        [path,name,ext]=fileparts('N:\data_fra_photon\LAB\OLAB\BENTHAM\Matlab\TOMS_mapping\toms\');
        %%        cd(path);
        if year <= 2004
            bins=288;
            res=1.25;
        else
            bins=360;
            res=1.00;
        end
        [toms,titulo]=read_toms_BJ(filename,bins);   % gamle data har 1.25 grader oppløsning (288 bins), mens nye har 1.00 grader (360)
        %%        cd(wd);
        %%        [X,Y]=meshgrid(-179.375:1.25:179.375,-89.5:1:89.5);
        [X,Y]=meshgrid(-179.5:res:179.5,-89.5:1:89.5);
        O3=interp2(X,Y,toms,long,lat)    %input: lengde,bredde
        %finn nærmeste to rader i Y (en rad for hver breddegrad)
        %finn nærmeste to kolonner i X (en kolonne for hver lengdegrad)
        %Sjekk om de fire nærmeste punktene er fri for 0-verdier
        pY=find(   abs(Y(:,1)-lat)   ==min(abs(Y(:,1)-lat)));
        tY=Y(:,1);
        tY(pY)=[];
        pY2=find(   abs(Y(:,1)-lat)   ==min(abs(tY(:,1)-lat)));
        clear tY
        
        pX=find(   abs(X(1,:)-long)   ==min(abs(X(1,:)-long)));
        tX=X(1,:);
        tX(pX)=[];
        pX2=find(   abs(X(1,:)-long)   ==min(abs(tX(1,:)-long)));
        clear tX
        O3_4points=[toms(pY,pX),toms(pY,pX2);toms(pY2,pX),toms(pY2,pX2)];
        if sum(sum(O3_4points>0))<4
            O3=0;
        end
        
        days(dnum)=dnum;    %=[days;dnum];
        ozone(dnum)=O3;  %=[ozone;O3];
        
        
    end
end %if ~isempty(nf)
%erstatt 0-verdier i ozone med data fra satelitt overpass: Må først
%lage en komplett serie av overpass ozon for alle dager mellom først og
%siste obsevasjon i omi-serien:

%lag vektor med alle dager fra start til slutt omi perioden:
TY1=OMI.Yyear(1);TY2=OMI.Yyear(length(OMI.Yyear));
TD1=OMI.Dyear(1);TD2=OMI.Dyear(length(OMI.Dyear));
ys=[TY1:TY2];
daysyear=365*ones(size(ys));
skuddaar=[1900:4:2100];
yvect=[];
dayvect=[];
jtvect=[];
ozone_vect=[];
for i=1:length(ys)
    if ismember(ys(i),skuddaar)
        daysyear(i)=366;
    end
    if ys(i)==TY1
        dager=[TD1:daysyear(i)]';
    elseif ys(i)==TY2
        dager=[1:TD2]';
    else
        dager=[1:daysyear(i)]';
    end
    
    yvect=[yvect;ys(i)*ones(length(dager),1)];
    dayvect=[dayvect;dager];
    jtvect=[jtvect;ys(i)+dager/daysyear(i)];
end

figure(year)
plot(year+days/365,ozone,'ro-','markerfacecolor','r'),
grid on,hold on %interpolert fra satelittkart
plot(OMI.JTyear,OMI.O3year,'ko-','markerfacecolor','k')
hold off
title(sprintf('%s in %d',station,year))
axis([year year+1 200 600])

if year>=2004.752  %OMI data fins herfra - fyll inn hull i ozone-vektoren fra hyppige omi overpasses
    pz=find(ozone==0);
    pomi=find(OMI.Yyear==year);
    [C,IA,IB]=intersect(days(pz),OMI.Dyear(pomi));
    ozone(pz(IA))=OMI.O3year(pomi(IB));
end
hold on
plot(year+days/365,ozone,'gs-','markersize',12),hold off
axis([year year+1 200 600])

figure(year*10)
plot(days,ozone,'ko-'),grid
%title(sprintf('ozone %s, year %04i',stat,aar))
title(sprintf('ozone %s, year %04i',stat,year))
axis([1 366 200 600])

if year<=2004
    name=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_epc_%s_%04i.txt',stat,stat,year);    %TOMS corrected, fra L3_ozone_epc_yyyymmdd.txt
else
    name=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_omi_%s_%04i.txt',stat,stat,year);
end

fid=fopen(name,'wt');
for j = 1 : length(ozone)
    fprintf(fid,'%i %6.4f\n',days(j),ozone(j));
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Beregner RES.Dailymean_O3
%Fyll først vektor med GUV_O3. Hvis daglig griddedete satelitt-data mangler, fyll inn huller
%med normalozon for Oslo, beregnet av Arne Dahlback.
%Hvis år>=1996 og år <=2004. Fyll inn EPC gridded data, https://ozonewatch.gsfc.nasa.gov/data/eptoms. Hvis EPC ikke
%finnes, men OMI gridded data finnes, fyll inn OMI (gjelder 2004).
%Hvis år>2004, bruk OMI overpass.


% Beregner ozon: To serier: En kvasi-serie RES.Dailymean_O3 (365x1) basert på guv hvor guvdata fins, og der ikke fins er
% det satt: O3 fra satelitt hvis omi/epc finnes (griddet ozon + overpass der det er hull i griddete), og ellers satt inn normalozon fra Oslo.

%Satelittbasert ozon: RES.MODEL.O3. Før satelitperioden er det brukt guv (<1996 og dag 207), ellers
%epc og omi fra kart og overpass. Innsetting av tomme satelitthull med normalozon for Oslo kan
%gjenkjennes ved at ozon-verdiene har heltallsverdier.

%year<1996: Hvis guvdata fins, bruk guvO3. Der det er hull-fyll inn normal ozon for Oslo.
% year>=2006 Last inn EPTOMS. Der TOMS dager fins, erstatt normalozon med TOMS.

%Hvis guv=O3 mangler, bruk default er normalozon Oslo. Last opp EPtoms (year <=2005). Erstatt med satelittdata der det fins.
% year>=2005: bruk default normalozon for oslo. Fyll inn med OmI data der de fins.

%Felter:RES.O3 (minuttverdier) og RES.Dailymean_O3. Legges over i RES.MODEL.O3, som brukes som ozonverdier for
%modellberegninger

%Se også: database_calc_O3_seasonal_mean_OMI.m
%Der beregnes sesongmidlet ozon fra OMI for perioden 2005 til 2010. Legges til struct OMI
% og RES.MODEL.O3_seasonal_mean. Brukes kun til beregning av normal-uv for stasjonene.
% I praksis er normalozon Oslo og sesongmidlet ozon fra OMI ganske like.

%clean

%run_uncalibrated=1;%kjører mot RAW.Raw_offset_corr som inneholder daswindata (W/m^2/nm) som ikke er korrigert for drift
run_uncalibrated=0;%Normalmodus - beregninger gjort på driftskorrigerte data RES.Raw_offset_corr
%

% switch stasjon
%     case 'Tromsø'
%         lat=69.68;
%         long=18.97;
%         stat='tso';
%         prim_sno='9276';
%         stasjon_norsk='Tromso';
%     case 'Andøya'
%         lat=69.28;
%         long=16.01;
%         stat='and';
%         prim_sno='9276';
%         stasjon_norsk='Andoya';
%     case 'Nyaal'
%         lat=78.92;
%         long=11.92;
%         stat='nya';
%         prim_sno='9275';
%         stasjon_norsk='Nyaalesund';
%     case 'Trondheim'
%         lat=63.42;
%         long=10.40;
%         stat='trh';
%         prim_sno='9274';
%         stasjon_norsk='Trondheim';
%     case 'Kise'
%         lat=60.78;
%         long=10.82;
%         stat='kis';
%         prim_sno='9272';
%         prim_sno_2='29229';
%         prim_sno_3='x';
%         stasjon_norsk='Kise';
%     case 'Bergen'
%         lat=60.38;
%         long=5.33;
%         stat='brg';
%         prim_sno='9270';
%         stasjon_norsk='Bergen';
%     case 'Østerås'
%         lat=59.946;
%         long=10.598;
%         stat='ost';
%         prim_sno='29222';
%         prim_sno_2='29229';
%         prim_sno_3='29243';
%         %stasjon_norsk='Østerås';
%         stasjon_norsk='Oesteraas';
%     case 'Blindern'
%         lat=59.938;
%         long=10.719;
%         stat='bli';
%         prim_sno='9222';
%         stasjon_norsk='Blindern';
%     case 'Landvik'
%         lat=58.33;
%         long=8.52;
%         stat='lan';
%         prim_sno='9271';
%         stasjon_norsk='Landvik';
%     case 'Finse'
%         lat=60.593;%koordinater fra google https://www.google.com/maps/d/viewer?mid=1upEzn-I3GlSP6zgs4pnZMuQgcNA&hl=en&ll=60.59299437534447%2C7.523682013244638&z=18
%         long=7.524;
%         stat='FIN';
%         prim_sno='29237';
%         stasjon_norsk='Finse';
%     case 'Kjeller'
%         lat=59.975;
%         long=11.053;
%         stat='KJE';
%         prim_sno='9222';
%         stasjon_norsk='Kjeller';
%     otherwise
%         return
% end

if size(RES.Raw,1)>0
    RES.O3=nan(size(RES.Raw,1),1);
    
    %            finn serienummer til instrument og tilhørende drifts og ANGF:
    SN=unique(RES.Inst);
    RES.Dailymean_O3=nan(size(RES.Daynums_year));
    for i=1:length(SN)  %for alle instrumenter som har stått på stasjonen i denne perioden
        id_pos=find(RES.Inst==SN(i));
        %drift_facts=load(sprintf('C:\\Bjorn\\Jobb\\uvnett\\Calib\\chdrift05_%i.txt',SN(i)));
        drift_facts=load(sprintf('N:\\uvnet\\guv\\Calib\\chdrift%02i_%i.txt',05,SN(i)));   %returnerer dd mm yyyy drift1..5 og offset1..5
        normfact=drift_facts(147,4:8);
        if SN(i)==9277
            SN(i)=29222;
        elseif SN(i)==9278
        elseif SN(i)==9279
            SN(i)=29229;
            SN(i)=29237;
        elseif SN(i)==9280
            SN(i)=29243;
        end
        %load(sprintf('C:\\Bjorn\\Jobb\\FARIN2005\\characterisation\\angle\\GUV\\%i\\ANGF_%i.mat',SN(i),SN(i)));
        load(sprintf('N:\\data_fra_photon\\LAB\\FARIN2005\\characterisation\\angle\\GUV\\%i\\ANGF_%i.mat',SN(i),SN(i)));
        
        %                 if SN(i)==9222  %sett inn en dummy 313-kanal i CIE-faktorene:
        %                     temp=CIE_fact(2:4);
        %                     CIE_fact(2)=0;
        %                     CIE_fact(3:5)=temp; clear temp
        %                 end
        [nx,ny]=size(RES.Raw_offset_corr(id_pos,:));
        
        % renormaliser rådata til 2005 - lag ozon rutine... beregn
        % for hver dag...
        days=unique(RES.Daynums(id_pos));
        Raw=100*RES.Raw_offset_corr(id_pos,:).*repmat(normfact,nx,1);   %renormalisert til Farin 2005 og regner om til uW/cm2/nm
        
        if SN(i)==9222   %mangler 313-kanal, fra databasen fås verdier 0 i 313-kolonnen
            temp=RES.Raw_offset_corr;temp(:,2:4)=temp(:,3:5);temp(:,5)=[];
            Raw=100*temp(id_pos,:).*repmat(normfact([1,3:5]),nx,1);   %renormalisert til Farin 2005 og regnet om til uW/cm2/nm
        end
        
        for j=1:length(days)
            daynum=days(j);
            p=find(RES.Daynums(id_pos)==daynum);
            
            % Landvik:            if year==2020 && SN(i)==29229 && daynum<91 %GUV9278. Dropp korreksjon, gjør bare at ozonverdiene blir svært feile
            %                         vinterkorr=[3 1.015 0.85 0.625 0.9]./[5.4 1.68 1.09 0.695 0.975];%dag 1-90
            %                         RAW.Raw=Raw(p,:)./repmat(vinterkorr,numel(p),1);       %døgndata
            %                     else
            
            RAW.Raw=Raw(p,:);       %døgndata
            %                     end
            RAW.Stasjon=stasjon;
            RAW.Stat=stat;
            RAW.Daynum=daynum;
            RAW.Year=year;
            RAW.Instrument=num2str(unique(RES.Inst(id_pos(p))));
            RAW.SZA=RES.SZA(id_pos(p));
            RAW.AZI=RES.AZI(id_pos(p));
            RAW.Time=(RES.JT(id_pos(p))-daynum)*24;
            RAW.SkyTrans=RES.SkyTrans(id_pos(p));
            %RAW.OMI_ovp_data=RES.Dailymean_O3_satelite_EPC_OMI(daynum,:);
            
            [tmp_time,tmp_sza,tmp_O3]=farin_calc_O3(ANGF,RAW,daynum,'asløkfjaf');   %for SZA<90, fyll inn resten av posisjonene:
            k=find(RAW.SZA<90);
            RAW.O3=nan(size(RAW.SZA));
            RAW.O3(k)=tmp_O3;
            pnoon=find(tmp_sza==min(tmp_sza));
            posnoon=find( (tmp_time>tmp_time(pnoon)-0.5) & (tmp_time<tmp_time(pnoon)+0.5) );
            RAW.MeanO3=mean(tmp_O3(posnoon));
            RES.O3(id_pos(p),1)=RAW.O3;
            RES.Dailymean_O3(daynum,1)=mean(tmp_O3(posnoon));
            
            clear RAW

        end   %                for j=1:length(days)
        
    end
    
    figure(year*10)
    yyaxis left
    plot(RES.Daynums_year,RES.Dailymean_O3,'ko-','MarkerSize',4,'MarkerFaceColor','k')%basert på guv
    
    norm_o3=load(sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\normalozon_oslo_79_90.txt'));
    norm_o3=[norm_o3;[366,norm_o3(365,2)]];
    sat=norm_o3;    %initaliserer, overskriver med satelitt data etterpå
    
    pnan=find(isnan(RES.Dailymean_O3)==1);
    
    if year<1996 %Bruk ozondata beregnet fra GUV på stasjonen - erstat nans med normalozon for Oslo 1979-1990.
        if ~isempty(pnan)
            RES.Dailymean_O3(pnan)=norm_o3(pnan,2);
        end
        d=sat(:,1); %normalozon for Oslo
        O=sat(:,2);
        
        %elseif year>=1996 & year<2004
    elseif year==1996% & year<2004 %huller fylles inn med EPC TOMS der slike data fins, ellers, bruk normalozon for Oslo, uansett stasjon
        %Bruk samme TOMS EPC data på Østerås som på Blindern
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            %fp=sprintf('C:\\Bjorn\\Jobb\\uvnett\\ozon\\%s\\O3_BLI_%i.txt',stat,year);
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
            %fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_BLI_%i.txt',stat,year);
        else %Bruk TOMS EPC data tilhørende hver stasjon
            %fp=sprintf('C:\\Bjorn\\Jobb\\uvnett\\ozon\\%s\\O3_%s_%i.txt',stat,stat,year);
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
            %fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat_stat,year);
        end
        epc=load(fp);
        pepc=find(epc(:,2)>200);
        day_meas=epc(pepc,1);
        sat(day_meas,2)=epc(pepc,2);
        d=sat(:,1);O=sat(:,2);clear sat
        
        if year==1996
            if ~isempty(pnan)
                p1=find(RES.Daynums_year(pnan)<=207);
                if ~isempty(p1)
                    RES.Dailymean_O3(pnan(p1))=norm_o3(pnan(p1),2);
                end
                
                p2=find(RES.Daynums_year(pnan)>207);
                if ~isempty(p2)
                    %sjekk om epc data fins på disse dagene:
                    if O(pnan(p2))>200
                        RES.Dailymean_O3(pnan(p2))=O(pnan(p2));
                    else
                        RES.Dailymean_O3(pnan(p2))=norm_o3(pnan(p2),2);
                    end
                end
                
            end
            
        end
        
    elseif year>1996 & year<2004
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
        else
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
        end
        epc=load(fp);
        pepc=find(epc(:,2)>200);
        day_meas=epc(pepc,1);
        sat(day_meas,2)=epc(pepc,2);
        d=sat(:,1);O=sat(:,2);clear sat
        if ~isempty(pnan)
            p2=find(RES.Daynums_year(pnan)<=366);
            if ~isempty(p2)
                %sjekk om epc data fins på disse dagene:
                if O(pnan(p2))>200
                    RES.Dailymean_O3(pnan(p2))=O(pnan(p2));
                else
                    RES.Dailymean_O3(pnan(p2))=norm_o3(pnan(p2),2);
                end
            end
            
        end
        
    elseif year==2004 %hent data fra epc, suppler med OMI hvis EPC har brudd:
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            fepc=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
            fomi=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_omi_BLI_%i.txt',year);
        else
            fepc=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
            fomi=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_omi_%s_%i.txt',stat,stat,year);
        end
        
        epc=load(fepc);
        omi=load(fomi);
        pepc=find(epc(:,2)>200);
        pomi=find(omi(:,2)>200);
        sat(pepc,2)=epc(pepc,2);
        sat(pomi,2)=omi(pomi,2);
        d=sat(:,1);O=sat(:,2);clear sat
        
        if ~isempty(pnan)
            p2=find(RES.Daynums_year(pnan)<=366);
            if ~isempty(p2)
                %sjekk om epc data fins på disse dagene:
                if O(pnan(p2))>200
                    RES.Dailymean_O3(pnan(p2))=O(pnan(p2));
                else
                    RES.Dailymean_O3(pnan(p2))=norm_o3(pnan(p2),2);
                end
            end
            
        end
        
    else
        %fp=sprintf('C:\\Bjorn\\Jobb\\uvnett\\ozon\\%s\\O3_omi_%s_%i.txt',stat,stat,year);
        fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_omi_%s_%i.txt',stat,stat,year);
        omi=load(fp);
        pomi=find(omi(:,2)>200);
        day_meas=omi(pomi,1);
        sat(day_meas,2)=omi(pomi,2);
        d=sat(:,1);
        O=sat(:,2);clear sat
        
        if ~isempty(pnan)
            p2=find(RES.Daynums_year(pnan)<=366);
            if ~isempty(p2)
                %sjekk om epc data fins på disse dagene:
                if O(pnan(p2))>200
                    RES.Dailymean_O3(pnan(p2))=O(pnan(p2));
                else
                    RES.Dailymean_O3(pnan(p2))=norm_o3(pnan(p2),2);
                end
            end
            
        end
        
    end
    
    hold on %fortsatt figure(year*10)
    yyaxis left
    plot(d,O,'ro-'),grid on,
    %             yyaxis right
    %             plot(RES.Daynums_year,RES.Dailymean_O3./O,'b.-')%basert på guv
    hold off %fra interpolerte satelittkart
    %axis([0 370 200 600])
    
    figure(year)
    plot(RES.JT,RES.O3,'k.-'), hold on
    plot(d+0.5,O,'rs-','markerfacecolor','r')
    hold off, grid on
    legend('GUV 1 min','OMI')
    title(sprintf('%s %d',stasjon,year))
    axis([0 370 200 600])
    ylabel('Ozone, DU')
    xlabel('Daynum')
    set(gca,'FontSize',14)
    
    RES.Days_Dailymean_O3_satelite_EPC_OMI=d;
    RES.Dailymean_O3_satelite_EPC_OMI=[d';O']';
    
    % Velg ozon-serie det skal klarværsmodelleres for:
    if year<1996
        RES.MODEL.O3=RES.Dailymean_O3(1:daysyear);
        
    elseif year==1996
        if ~isfield(RES,'MODEL')
            RES.MODEL=[];
            RES.MODEL.O3=[];
        end
        %                 if isfield(RES.MODEL,'O3')
        %                     RES.MODEL.O3=[];
        %                 end
        p=find(RES.Daynums_year<=207);
        RES.MODEL.O3(p,1)=RES.Dailymean_O3(p);
        p=find(RES.Daynums_year>207 & RES.Daynums_year <= daysyear);
        RES.MODEL.O3(p,1)=O(p);
    else
        RES.MODEL.O3=O(1:daysyear);
    end
    clear O
    
    %close(year*10)
    
    figure(year)
    hold on
    %plot(RES.Daynums_year,RES.MODEL.O3,'rs-','MarkerSize',5),grid on
    legend('GUV daily mean','OMI')
    title(sprintf('%s %d',stasjon,year))
    hold off
    
else %if size(RES.Raw,1)==0
    
    %Legg til felt RES.MODEL.O3 basert på satelitt data eller normaldata
    norm_o3=load(sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\normalozon_oslo_79_90.txt'));
    norm_o3=[norm_o3;[366,norm_o3(365,2)]];
    sat=norm_o3;    %initaliserer, overskriver med satelitt data etterpå
    
    if year < 1996 %Bruk ozondata beregnet fra GUV på stasjonen - erstat nans med normalozon for Oslo 1979-1990.
        %                 if ~isempty(pnan)
        %                     RES.Dailymean_O3(pnan)=norm_o3(pnan,2);
        %                 end
        d=sat(:,1); %normalozon for Oslo
        O=sat(:,2);
        RES.MODEL.O3(1:daysyear,1)=O(1:daysyear);
        
    elseif year==1996% & year<2004 %huller fylles inn med EPC TOMS der slike data fins, ellers, bruk normalozon for Oslo, uansett stasjon
        %Bruk samme TOMS EPC data på Østerås som på Blindern
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
        else %Bruk TOMS EPC data tilhørende hver stasjon
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
        end
        epc=load(fp);
        pepc=find(epc(:,2)>200 & ~isnan(epc(:,2)));
        day_meas=epc(pepc,1);
        sat(day_meas,2)=epc(pepc,2);
        d=sat(:,1);O=sat(:,2);clear sat
        
        %if year==1996
        RES.MODEL.O3(207:daysyear,1)=O(207:daysyear);
        RES.MODEL.O3(1:207-1,1)=norm_o3(1:207-1,2);
        pn=find(isnan(RES.MODEL.O3)==1);
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end

    elseif year>1996 & year<2004
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
        else
            fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
        end
        epc=load(fp);
        pepc=find(epc(:,2)>200 & ~isnan(epc(:,2)));
        day_meas=epc(pepc,1);
        %                 sat(day_meas,2)=epc(pepc,2);
        %                 d=sat(:,1);O=sat(:,2);clear sat
        
        %NB! skal bruke hele årsserien når year>1996
        %                 RES.MODEL.O3(207:daysyear,1)=O(207:daysyear);
        %                 RES.MODEL.O3(1:207-1,1)=norm_o3(1:207-1);
        RES.MODEL.O3=zeros(daysyear,1);
        RES.MODEL.O3(day_meas)=epc(pepc,2);
        pn=find(isnan(RES.MODEL.O3)==1)
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end
        pn=find(RES.MODEL.O3==0);
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end

        
    elseif year==2004 %hent data fra epc, suppler med OMI hvis EPC har brudd:
        if strcmp(stat,'ost')||strcmp(stat,'OST')
            fepc=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_EPC_BLI_%i.txt',year);
            fomi=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\BLI\\O3_omi_BLI_%i.txt',year);
        else
            fepc=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_EPC_%s_%i.txt',stat,stat,year);
            fomi=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_omi_%s_%i.txt',stat,stat,year);
        end
        
        epc=load(fepc);
        omi=load(fomi);
        pepc=find(epc(:,2)>200);
        pomi=find(omi(:,2)>200);
        day_meas=epc(pepc,1);
        %sat=norm_o3;
        sat(pepc,2)=epc(pepc,2);
        sat(pomi,2)=omi(pomi,2);
        d=sat(:,1);O=sat(:,2);clear sat
        
        RES.MODEL.O3=zeros(daysyear,1);
        RES.MODEL.O3(day_meas)=epc(pepc,2);
        pn=find(isnan(RES.MODEL.O3)==1)
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end
        pn=find(RES.MODEL.O3==0);
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end
        
    else
        fp=sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\ozon\\%s\\O3_omi_%s_%i.txt',stat,stat,year);
        omi=load(fp);
        pomi=find(omi(:,2)>200);
        day_meas=omi(pomi,1);
        sat(pomi,2)=omi(pomi,2);
        d=sat(:,1);O=sat(:,2);clear sat
        O=O(1:daysyear);
        
        RES.MODEL.O3=zeros(daysyear,1);
        RES.MODEL.O3=O;
        pn=find(isnan(RES.MODEL.O3)==1)
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end
        pn=find(RES.MODEL.O3==0);
        if ~isempty(pn)
            RES.MODEL.O3(pn)=norm_o3(pn,2);
        end

    end
    
    figure(year)
    hold on
    plot(RES.MODEL.O3,'k.-'),grid on,hold off
    axis([0 370 200 600])
    ylabel('Ozone, DU')
    xlabel('Daynum')
    set(gca,'FontSize',14)

    clear O
    
end  %if size(RES.Raw,1)>0

if(any(findall(0,'Type','Figure')==1))
    close(1);
end

%figure(1), plot([1:numel(RES.MODEL.O3)],RES.MODEL.O3,'r.-')

% 
% if run_uncalibrated==1
%     RAW=RES;
%     clear RES
%     tmp = 'RAW';
%     save(fpath,tmp)
%     clear RAW
% else
%     tmp = 'RES';
%     save(fpath,tmp)
%     clear RES;
% end

figure(year*10)
set(gcf, 'Position', get(0, 'Screensize'));%fullscreen size - for best visning av graf
pause(2)
%proot=sprintf('E:\\Data\\work\\uvnett\\UV_uncorrected_GUV_data\\');
%gpath = sprintf('%sfigures\\Ozone\\%sozone_%s_%04i.fig',proot,prename,stasjon_norsk,year);
%savefig(gpath)
%
%         sprintf('%s\\figures\\Ozone\\',proot,year
%         E:\Data\work\uvnett\UV_uncorrected_GUV_data\figures\Ozone
%         E:\Data\work\uvnett\UV\Figures\Ozone

%end %for year=[1996:1997]
%
%     close all
%
% end %for stasjons_id=1:2



