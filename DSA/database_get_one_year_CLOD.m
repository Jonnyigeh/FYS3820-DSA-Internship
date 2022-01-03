function RES=database_get_one_year_CLOD(RES)
%Beregner skytransmittans
pc='Dell790';  %DEL Win7 stasjonær

switch pc
    case 'dolly'
        proot=('N:\\data_fra_photon\\LAB\\OLAB\\UV_NETT\\Totalstraaling\\UV\\');
    case 'bb'
        proot=sprintf('C:\\Bjorn\\Jobb\\uvnett\\UV\\');
    case 'Dell'
        proot=sprintf('D:\\work\\uvnett\\UV\\');
        
    case 'Dell790'
        proot=sprintf('D:\\Data\\work\\uvnett\\UV\\');
        
    otherwise
        return
end

stasjons_id=RES.Stasjons_id;
station=RES.Stasjon;
lat=RES.Latitude;
long=RES.Longitude;
year=RES.Year;

load('uvspec_clear_albedo.mat');
UVSPEC_Albedo=UVSPEC; %SZA, Albedo, bølgelengde 9 * 6 * 241
clear UVSPEC
[X,Y]=meshgrid(UVSPEC_Albedo.Albedo,UVSPEC_Albedo.SZA);%X=alb 0:0.20:1.00, repetert SZA rader, Y= kolonner sza 5:85, repetert albedo kolonner
Z=squeeze(UVSPEC_Albedo.Uvspec(:,:,UVSPEC_Albedo.Wavelength==380));%uvspec ved 380nm som funksjon av sza og albedo, 9*6
UVSPEC_Albedo.X=X;
UVSPEC_Albedo.Y=Y;
UVSPEC_Albedo.Z=Z;
clear X Y Z

toleranse=1.01; %tillater at N-verdi (signal på kanal 380) inntil 3% større enn teoretisk klarværs N
%clear_OD=1.5;   %skille mellom clear/cloudy - må være lik med farin_clod_cosinefit
%%clear_OD=3.0;
clear_OD=0;    %betyr at det benyttes en fast diffus-korreksjon på begge datasettene

%for stasjons_id=1:1
%%for stasjons_id=1:11
%for stasjons_id=2:2
%for stasjons_id=4:4
%for stasjons_id=6:6
%%for stasjons_id=2:3
%for stasjons_id=5:11
%for stasjons_id=8:10
%for stasjons_id=5:5
%for stasjons_id=10:10
%for stasjons_id=11:11
%for stasjons_id=9:9
%%for stasjons_id=8:8
%for stasjons_id=1:11
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
        %stasjon='Nyaalesund';
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

switch stasjon
    case 'Tromsø'
        lat=69.68;
        long=18.97;
        stat='tso';
        prim_sno='9276';
        stasjon_norsk='Tromso';
    case 'Andøya'
        lat=69.28;
        long=16.01;
        stat='and';
        prim_sno='9276';
        stasjon_norsk='Andoya';
    case 'Nyaal'
        lat=78.92;
        long=11.92;
        stat='nya';
        prim_sno='9275';
        stasjon_norsk='Nyaalesund';
    case 'Trondheim'
        lat=63.42;
        long=10.40;
        stat='trh';
        prim_sno='9274';
        stasjon_norsk='Trondheim';
    case 'Kise'
        lat=60.78;
        long=10.82;
        stat='kis';
        prim_sno='9272';
        prim_sno_2='29229';
        prim_sno_3='x';
        stasjon_norsk='Kise';
    case 'Bergen'
        lat=60.38;
        long=5.33;
        stat='brg';
        prim_sno='9270';
        stasjon_norsk='Bergen';
    case 'Østerås'
        lat=59.946;
        long=10.598;
        stat='ost';
        prim_sno='29222';
        prim_sno_2='29229';
        prim_sno_3='29243';
        %stasjon_norsk='Østerås';
        stasjon_norsk='Oesteraas';
    case 'Blindern'
        lat=59.938;
        long=10.719;
        stat='bli';
        prim_sno='9222';
        stasjon_norsk='Blindern';
    case 'Landvik'
        lat=58.33;
        long=8.52;
        stat='lan';
        prim_sno='9271';
        stasjon_norsk='Landvik';
    case 'Finse'
        lat=60.60;
        long=7.52;
        stat='FIN';
        prim_sno='29237';
        stasjon_norsk='Finse';
    case 'Kjeller'
        lat=59.975;
        long=11.053;
        stat='KJE';
        prim_sno='9222';
        stasjon_norsk='Kjeller';
    otherwise
        return
end

tic
%%lastyear=2013;
%for year=2012:2016
%for year=1995:2018
%for year=1996:2018
%for year=2016:2016
%for year=2014:2014
%for year=2018:2018
%for year=2019:2019
%for year=2019:2019
%for year=2020:2020
%for year=2009:2009
%for year=2020:2020
%for year=1995:1996
%%for year=2021:2021

%         %for year=lastyear:lastyear
%         yy=num2str(year); yy=yy(3:4);yy=str2num(yy);
%         %path = sprintf('C:\\Bjorn\\Jobb\\uvnett\\UV\\%03s\\UV_DB_%s_%04i.mat',stat,stasjon,year);
%         path = sprintf('%s%03s\\UV_DB_%s_%04i.mat',proot,stat,stasjon_norsk,year);
%         load(path)
%
%         sprintf('%s\t%i',stasjon,year)

if size(RES.Raw,1)>0
    
    %            finn serienummer til instrument og tilhørende drifts og ANGF:
    SN=unique(RES.Inst);
    RES.Dailymean_CLOD=nan(size(RES.Daynums_year));
    RES.Dailymean_SkyTrans=nan(size(RES.Daynums_year));
    for i=1:length(SN)  %for alle instrumenter som har stått på stasjonen i denne perioden
        id_pos=find(RES.Inst==SN(i));
        %drift_facts=load(sprintf('C:\\Bjorn\\Jobb\\uvnett\\Calib\\chdrift05_%i.txt',SN(i)));
        drift_facts=load(sprintf('N:\\uvnet\\guv\\Calib\\chdrift%02i_%i.txt',05,SN(i)));   %returnerer dd mm yyyy drift1..5 og offset1..5
        normfact=drift_facts(147,4:8);
        if SN(i)==9277
            SN(i)=29222;
        elseif SN(i)==9278
            SN(i)=29229;
        elseif SN(i)==9279
            SN(i)=29237;
        elseif SN(i)==9280
            SN(i)=29243;
        end
        
        load(sprintf('N:\\data_fra_photon\\LAB\\FARIN2005\\characterisation\\angle\\GUV\\%i\\ANGF_%i.mat',SN(i),SN(i)));
        
        [nx,ny]=size(RES.Raw_offset_corr(id_pos,:));
        
        % renormaliser rådata til 2005 - lag ozon rutine... beregn
        % for hver dag...
        days=unique(RES.Daynums(id_pos));
        Raw=100*RES.Raw_offset_corr(id_pos,:).*repmat(normfact,nx,1);   %renormalisert til Farin 2005 og regnet om til uW/cm2/nm
        if SN(i)==9222   %mangler 313-kanal, fra databasen fås verdier 0 i 313-kolonnen
            temp=RES.Raw_offset_corr;temp(:,2:4)=temp(:,3:5);temp(:,5)=[];
            Raw=100*temp(id_pos,:).*repmat(normfact([1,3:5]),nx,1);   %renormalisert til Farin 2005 og regnet om til uW/cm2/nm
        end
        
        for j=1:length(days)
            daynum=days(j)
            p=find(RES.Daynums(id_pos)==daynum);
            RAW.Raw=Raw(p,:);       %døgndata
            RAW.SZA=RES.SZA(id_pos(p));
            RAW.Time=(RES.JT(id_pos(p))-daynum)*24;
            
            % beregn CLOD
            albedo=0.05;
            % alb_corr=ALB.alb_corr;
            alb_corr=1.036; %faktor som Uvspec_0.05 må mulitiplseres med for at klarværstransmittans målt med 380-kanalen skal matche albedo 5 porsent under FARIN 2005
            %%tmp_CLOD=farin_calc_OD_old(RAW.SZA,RAW.Time,ANGF,RAW.Raw(:,ANGF.chan_OD),daynum,toleranse,'asløkfjaf',0);
            [tmp_CLOD,tmp_SkyTrans]=farin_calc_OD_old(albedo,alb_corr,UVSPEC_Albedo,RAW.SZA,RAW.Time,ANGF,RAW.Raw(:,ANGF.chan_OD),daynum,toleranse,'asløkfjaf',0);
            
            %albedo,alb_corr,UVSPEC_Albedo,
            %[tmp_CLOD,n_clear]=farin_calc_OD(albedo,alb_corr,UVSPEC_Albedo,SZA,Time,ANGF,n_meas,dnum,toleranse,mat_name,do_plot)
            
            
            %                    [tmp_time,tmp_sza,tmp_O3]=farin_calc_O3(ANGF,RAW,daynum,'asløkfjaf');   %for SZA<90, fyll inn resten av posisjonene:
            k=find(RAW.SZA<90);
            %RAW.CLOD=nan(size(RAW.SZA));
            %RAW.CLOD(k)=tmp_CLOD;
            RAW.CLOD=tmp_CLOD;
            tmp_sza=RAW.SZA;
            pnoon=find(tmp_sza==min(tmp_sza));
            tmp_time=RAW.Time;
            posnoon=find( (tmp_time>tmp_time(pnoon)-0.5) & (tmp_time<tmp_time(pnoon)+0.5) );
            RAW.MeanCLOD=mean(tmp_CLOD(posnoon));
            RES.CLOD(id_pos(p),1)=RAW.CLOD;
            RES.Dailymean_CLOD(daynum,1)=mean(tmp_CLOD(posnoon));
            
            RAW.SkyTrans=tmp_SkyTrans;
            RAW.MeanSkyTrans=mean(tmp_SkyTrans(posnoon));
            RES.SkyTrans(id_pos(p),1)=RAW.SkyTrans;
            RES.Dailymean_SkyTrans(daynum,1)=mean(tmp_SkyTrans(posnoon));
            
            clear RAW
            
        end   %                for j=1:length(days)
        
    end
    
    figure(year*10)
    yyaxis left
    plot(RES.Daynums_year,RES.Dailymean_CLOD,'ko-')
    yyaxis right
    plot(RES.Daynums_year,RES.Dailymean_SkyTrans,'ro-')
    legend('CLOD','Transmittans')
    grid on
    title(sprintf('%s %d',stat,year))
    
    
    figure(year)
    plot(RES.JT,RES.CLOD,'r-')
    hold on
    plot(RES.Daynums_year,RES.Dailymean_CLOD,'ko-')
    hold off
    g=gca;
    g.YLim(1)=-1;
    
    figure
    plot(RES.JT,RES.SkyTrans,'k-')
    hold on
    plot(RES.Daynums_year,RES.Dailymean_SkyTrans,'ro-')
    hold off, grid on
    title(sprintf('%s %d',stat,year))
    g=gca;
    g.YLim(1)=0;
    set(gcf, 'Position', get(0, 'Screensize'));%fullscreen size - for best visning av graf
    pause(2)
    
    %             tmp = 'RES';
    %             save(path,tmp)
    
end  %if size(RES.Raw,1)>0
%
%         clear RES;
%     end %for year=[1996:1997]
%     close all

%end %for stasjons_id=1:2
