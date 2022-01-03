function RES=database_get_one_year_SZA(RES)
%leser drifts-korrigerte rådata fra mat-fil (eksportert fra databasen) år
%for år, beregner SZA og tid.
%Beregner residual offset (sza>=100), og trekker fra rådata
% Beregner offsetkorrigerte kanalverdier RES.Raw_offset_corr
%Gjør samtidig korreksjon for temperatur som avviker fra 40C. RES.Raw_offset_corr er her
%utgangspunkt for alle beregninger av doseprodukter
%clean

%Les inn temperaturfunksjon:
angf_path=sprintf('N:\\data_fra_photon\\LAB\\FARIN2005\\characterisation\\angle\\GUV\\9273\\ANGF_9273.mat');
load(angf_path);
temp_respons=ANGF.Temp_respons_qth_nettosignal;%Antar at alle guv'en har lik temperaturerespons, og at den er lik i sol som mot en QTH lampe
clear ANGF

T_lim_L=0; %korrigerer ikke for temperatur hvis T ligger under 0 grader
T_lim_U=50;%korrigerer ikke for temperatur hvis T er over 50 grader
%temp_respons er kun målt opp i intervallet 1C til 46C. Ukjent hva responsen er utenfor dette
%temperaturintervallet, så legg til en sikker lavere og høyere skranke:

lavT=[-50,temp_respons(1,2:6)];
hoyT=[100,temp_respons(end,2:6)];
temp_respons=[lavT;temp_respons;hoyT];
year=RES.Year;
stasjons_id=RES.Stasjons_id;
stasjon=RES.Stasjon;

% %%for stasjons_id=1:11
%     %for stasjons_id=9:11
%     %for stasjons_id=8:11
% %    for stasjons_id=8:8%Nya
%     %for stasjons_id=9:9
%     %for stasjons_id=10:10%Finse
%     %for stasjons_id=11:11 %Kjeller
%     %for stasjons_id=7:7
%     %for stasjons_id=1:5
%     for stasjons_id=1:1
%     %for stasjons_id=2:2
%     %for stasjons_id=5:5%Blindern
%     %for stasjons_id=3:3%bergen
%     %for stasjons_id=4:4%Landvik
%     %for stasjons_id=6:6%trondheim
%     %for stasjons_id=3:3
%     %%for stasjons_id=9:9
%     %%for stasjons_id=10:10
%     switch stasjons_id
%         case 1
%             stasjon='Østerås';
%             stasjon_norsk='Oesteraas';
%             stat='OST';
%         case 2
%             stasjon='Kise';
%             stasjon_norsk='Kise';
%             stat='KIS';
%         case 3
%             stasjon='Bergen';
%             stasjon_norsk='Bergen';
%             stat='BRG';
%         case 4
%             stasjon='Landvik';
%             stasjon_norsk='Landvik';
%             stat='LAN';
%         case 5
%             stasjon='Blindern';
%             stasjon_norsk='Blindern';
%             stat='BLI';
%         case 6
%             stasjon='Trondheim';
%             stasjon_norsk='Trondheim';
%             stat='TRH';
%         case 7
%             stasjon='Tromsø';
%             stasjon_norsk='Tromso';
%             stat='TSO';
%         case 8
%             stasjon='Nyaalesund';
%             stasjon_norsk='Nyaalesund';
%             stat='NYA';
%         case 9
%             stasjon='Andøya';
%             stasjon_norsk='Andoya';
%             stat='AND';
%         case 10
%             stasjon='Finse';
%             stasjon_norsk='Finse';
%             stat='FIN';
%         case 11
%             stasjon='Kjeller';
%             stasjon_norsk='Kjeller';
%             stat='KJE';
%             
%         otherwise
%             return
%     end
%     
    %tic
    %for year=1995:1996
    %for year=1994:2013
    %    for year=2001:2001
    %for year=[2006,2008]
    %    for year=1995:1995
    %for year=2010:2011
    %for year=2012:2013
    %for year=2013:2013
    %for year=2016:2016
    %for year=2017:2017
    %for year=2018:2018
    %for year=2019:2019
    %for year=2020:2021
    %for year=2020:2020
    %%for year=2021:2021
        %         RES.Station=[];
        %         RES.DateHour=[];
        %         RES.Inst=[];
        %         RES.CIEDB=[];
        %         RES.Raw=[];
        %         RES.Daynums=[];
        %     RES.Latitude=rlat;
        %     RES.Longitude=rlong;
        
        %            path = sprintf('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\UV\\%03s\\UV_DB_%s_%04i.mat',stat,stasjon,year);
        %	path = sprintf('C:\\Bjorn\\Jobb\\uvnett\\UV\\%03s\\UV_DB_%s_%04i.mat',stat,stasjon,year);
        
        if (rem(year,4)==0) & (rem(year,400)~=0)
            daysyear=366;
        else
            daysyear=365;
        end
        if year==2000
            daysyear=366;
        end
        
        
%         path = sprintf('%s%03s\\UV_DB_%s_%04i.mat',proot,stat,stasjon_norsk,year);
%         load(path)
        
        RES.Daynums_year=(1:daysyear)';
        
        %sprintf('%s\t%i',stasjon,year)
        %        disp([stasjon,year])
        if size(RES.Raw,1)>0
            RES.JT=RES.Daynums+(str2num(RES.DateHour(:,10:11))+str2num(RES.DateHour(:,13:14))/60)/24;
            
            %sorter etter økende datotid:
            A=RES.JT;
            [~,index] = sortrows(A);
            figure
            plot(diff(index),'.-')
            %sorter alle felter i RES
            RES.JT=RES.JT(index);
            RES.DateHour=RES.DateHour(index,:);
            RES.InstRES=RES.Inst(index);
            RES.CIEDB=RES.CIEDB(index);
            RES.Daynums=RES.Daynums(index);
            RES.Driftsfactors=RES.Driftsfactors(index);
            RES.Offsetfactors=RES.Offsetfactors(index);
            RES.Temperatur=RES.Temperatur(index);
            
            %Beregn SZA og AZI
            Days=unique(RES.Daynums);   %unike dager med målinger
            RES.JD_unique=Days;
            
            %fjern første minuttverdi i hvert døgn, da den potensielt har
            %avvikende verdi (tid 00:00), ses tydeligst for Nyaalesund:
            tic
            q=[];
            for j=1:length(Days)
                p=find(floor(RES.JT)==Days(j));
                if ~isempty(p)
                    p1=p(1);
                    tid_start=(RES.JT(p1)-Days(j))*24;
                    if tid_start < 1/60
                        q=[q;p1];
                    end
                end
            end
            if ~isempty(q)
                RES.DateHour(q,:)=[];
                RES.Inst(q)=[];
                RES.CIEDB(q)=[];
                RES.Raw(q,:)=[];
                RES.Daynums(q)=[];
                RES.Driftsfactors(q)=[];
                RES.Offsetfactors(q)=[];
                RES.Temperatur(q)=[];
                RES.JT(q)=[];
            end
            toc
            
            Days=unique(RES.Daynums);   %unike dager med målinger
            RES.JD_unique=Days;
            RES.Offset=zeros(366,size(RES.Raw,2)); %Residual offset, lik signal over sza 100 grader
            raw_corrected=RES.Raw;  %erstatter verdier etterpå med offsetkorrigerte
            SZA=zeros(size(RES.Raw,1),1);
            AZI=SZA;
            
            for j=1:length(Days)
                daynum=Days(j);
                %mnd og dag fra dagnr:
                [aar,maaned,dag]=daynum2date_guv(year,daynum);
                % hh,mm,ss:
                p=find(RES.Daynums==daynum);
                Dechour=(RES.JT(p)-daynum)*24;
                time=floor(Dechour);
                minutt=round( (Dechour-time)*60);   %korrigert for 30 sek for tidlig minuttangivelse
                sekund=round(( (Dechour - time)*60-minutt)*60);
                aar=year.*ones(size(p));
                maaned=maaned*ones(size(p));
                dag=dag*ones(size(p));
                
                azi=zeros(size(p)); % azimuth vinkel
                zen=zeros(size(p)); % zenith vinkel
                %                 for i=1:length(p) %NB! i kallet maa hour økes med 1 pga feil i solpos rutina
                %                     [azi(i),zen(i)]=solpos(aar(i),maaned(i),dag(i),time(i)+1,minutt(i),sekund(i),RES.Latitude,RES.Longitude);
                %                 end
                [azi,zen]=solpos_bj(aar,maaned,dag,time+1,minutt,sekund,RES.Latitude,RES.Longitude);    %tar inn en hel tidsvektor, mens rutina over virker på skalar tid
                
                %                RES.SZA=[RES.SZA;zen.*180/pi];
                %                RES.AZI=[RES.AZI;azi.*180/pi];
                SZA(p,1)=zen*180/pi;
                AZI(p,1)=azi*180/pi;
                
                sz_max=max(zen.*180/pi);
                if sz_max>=100
                    q=find(zen.*180/pi>=100);
                    %mn=mean(RES.Raw(p(q),:));
                    %md=median(RES.Raw(p(q),:));
                    Offset=median(RES.Raw(p(q),:));
                else
                    if j==1
                        Offset=zeros(1,size(RES.Raw,2));
                    else
                        Offset=RES.Offset(j-1,:);   %settes lik verdiene dagen før
                    end
                end
                %                RES.Offset(Days(j),:)=Offset;
                RES.Offset(Days(j),:)=Offset;
            end
            RES.SZA=SZA;
            RES.AZI=AZI;
            
            %             %             %Beregn residual offset før CIE: finn først SZA_max for støy til hver kanal
            %             figure(year)
            %             for i=1:size(RES.Raw,2);
            %                 subplot(2,3,i)
            %                 semilogy(RES.SZA,RES.Raw(:,i),'k.'),grid on
            %             end
            
            
            
            %raw_corrected=[];
            for j=1:length(RES.JD_unique)
                daynum=RES.JD_unique(j);
                p=find(RES.Daynums==daynum);
                %                raw_corrected=[raw_corrected;RES.Raw(p,:)-repmat(RES.Offset(daynum,:),length(p),1)];
                %raw_corrected=[raw_corrected;RES.Raw(p,:)-repmat(RES.Offset(j,:),length(p),1)];
                raw_corrected(p,:)=RES.Raw(p,:)-repmat(RES.Offset(j,:),length(p),1);
            end
            RES.Raw_offset_corr=raw_corrected;
            
            
            if RES.Inst~=9222
                p=find(RES.Temperatur>T_lim_L & RES.Temperatur<T_lim_U);
                TCorr=ones(size(RES.Raw));
                TCorr(p,:)=interp1(temp_respons(:,1),temp_respons(:,2:end),RES.Temperatur(p));
                RES.Raw_offset_corr=RES.Raw_offset_corr./TCorr;
            end
            
%             tmp = 'RES';
%             save(path,tmp)
            
            %             for i=1:size(RES.Raw,2);
            %                 figure(year+i)
            %                 hold on
            %                 semilogy(RES.SZA,RES.Raw_offset_corr(:,i),'r.'),grid on
            %                 hold off
            %             end
            %
            %             close all
        end  %if size(RES.Raw,1)>0
        %        close all
        
%         if year>=2020 && ~ismember(stasjons_id,[5,7])
%             figure(stasjons_id)
%             subplot(2,1,1)
%             plot(RES.JT,40*RES.CIEDB,'k-'),grid on
%             title(sprintf('%s %d',stasjon,year))
%             ylabel('UVI')
%             plot_name=sprintf('%s\\tmp_graphs\\%s_%d.png',proot,stasjon,year);
%             %export_fig(sprintf('%s',plot_name),'-m1');
%             
%             subplot(2,1,2)
%             plot(RES.JT,RES.Temperatur,'k-'),grid on
%             title(sprintf('%s %d',stasjon,year))
%             ylabel('Temperatur')
%             ylim([0 50])
%         end
        
%         clear RES;
%     end %for year=[1996:1997]
%     %close all
    
%     toc
% end %for stasjons_id=1:2




