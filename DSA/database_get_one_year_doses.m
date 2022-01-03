function RES=database_get_one_year_doses(RES)
%Beregner 11 sett av produktdata
%%renormaliserer data til dag 147 2005, beregner CIE, uvb og uva med faktoerer fra
%%Farin

%run_uncalibrated=1;%kjører mot RAW.Raw_offset_corr som inneholder daswindata (W/m^2/nm) som ikke er korrigert for drift
run_uncalibrated=0;%Normalmodus - beregninger gjort på driftskorrigerte data RES.Raw_offset_corr

stasjons_id=RES.Stasjons_id;
station=RES.Stasjon;
lat=RES.Latitude;
long=RES.Longitude;
year=RES.Year;

%pc='bb';
%pc='dolly';
%pc='Dell';  %Dell XP arbeidsstasjon
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
        prename='';
        if run_uncalibrated==1
            proot=sprintf('E:\\Data\\work\\uvnett\\UV_uncorrected_GUV_data\\');
            prename='RAW_';
        end
        
    otherwise
        return
end

switch stasjons_id
    case 1
        stasjon='Østerås';
        stasjon_norsk='Oesteraas';
        stat='OST';
        prim_sno='29222';
    case 2
        stasjon='Kise';
        stasjon_norsk='Kise';
        stat='KIS';
        prim_sno='9272';
    case 3
        stasjon='Bergen';
        stasjon_norsk='Bergen';
        stat='BRG';
        prim_sno='9270';
    case 4
        stasjon='Landvik';
        stasjon_norsk='Landvik';
        stat='LAN';
        prim_sno='9271';
    case 5
        stasjon='Blindern';
        stasjon_norsk='Blindern';
        stat='BLI';
        prim_sno='9222';
    case 6
        stasjon='Trondheim';
        stasjon_norsk='Trondheim';
        stat='TRH';
        prim_sno='9274';
    case 7
        stasjon='Tromsø';
        stasjon_norsk='Tromso';
        stat='TSO';
        prim_sno='9276';
    case 8
        stasjon='Nyaalesund';
        stasjon_norsk='Nyaalesund';
        stat='NYA';
        prim_sno='9275';
    case 9
        stasjon='Andøya';
        stasjon_norsk='Andoya';
        stat='AND';
        prim_sno='9276';
    case 10
        stasjon='Finse';
        stasjon_norsk='Finse';
        stat='FIN';
        prim_sno='29237';
    case 11
        stasjon='Kjeller';
        stasjon_norsk='Kjeller';
        stat='KJE';
        prim_sno='9222';
    otherwise
        return
end
annual=[];
numC=5;

if size(RES.Raw,1)>0
    RES.CIE_FARIN=0*RES.JT;%CIE-vektet irradians
    RES.UVB_FARIN=0*RES.JT;%uvektet uvb
    RES.UVA_FARIN=0*RES.JT;%uvektet uva
    RES.DNA_FARIN=0*RES.JT;%DNA-vektet
    RES.VITD_FARIN=0*RES.JT;%
    RES.ANCHOVY_FARIN=0*RES.JT;%
    RES.FLINTPLANT_FARIN=0*RES.JT;
    RES.HEP_FARIN=0*RES.JT;%
    RES.PAR_FARIN=0*RES.JT;%
    RES.CIE1998_FARIN=0*RES.JT;%
    RES.NMSC_FARIN=0*RES.JT;%
    
    %            finn serienummer til instrument og tilhørende drifts og cie faktorer:
    SN=unique(RES.Inst);
    for i=1:length(SN)  %for alle instrumenter som har stått på stasjonen i denne perioden
        id_pos=find(RES.Inst==SN(i));
        
        sno_long=SN(i);
        if SN(i)==9277
            %SN(i)=29222;
            sno_long=29222;
        elseif SN(i)==9278
            %SN(i)=29229;
            sno_long=29229;
        elseif SN(i)==9279
            %SN(i)=29237;
            sno_long=29237;
        elseif SN(i)==9280
            %SN(i)=29243;
            sno_long=29243;
        end
        
        d=load(sprintf('N:\\uvnet\\guv\\Calib\\chdrift%02i_%i.txt',05,RES.Inst(id_pos(1))));   %returnerer dd mm yyyy drift1..5 og offset1..5
        Norm=d(147,[1:numC]+3);
        
        %opplasting av ANGF:
        if sno_long==29243
            load(sprintf('N:\\data_fra_photon\\LAB\\FARIN2005\\characterisation\\angle\\GUV\\%i\\ANGF_%i_HMG.mat',sno_long,sno_long));
        else
            load(sprintf('N:\\data_fra_photon\\LAB\\FARIN2005\\characterisation\\angle\\GUV\\%i\\ANGF_%i.mat',sno_long,sno_long));
        end
        
        w_CIE=ANGF.CIE*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
        w_uvb=ANGF.uvb_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
        w_uva=ANGF.uva_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
        w_DNA=ANGF.DNA_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
        w_vitD=ANGF.vitD_coeffs*100;
        w_anchovy=ANGF.anchovy_coeffs*100;
        w_flintplant=ANGF.flintplant_coeffs*100;
        w_HEP=ANGF.HEP_coeffs*100;
        w_PAR=ANGF.PAR_coeffs*100;
        w_CIE1998=ANGF.CIE1998_coeffs*100;
        w_NMSC=ANGF.NMSC_coeffs*100;
        %end
        
        %if strcmp(guv_name,'9222')  %sett inn en dummy 313-kanal i CIE-faktorene:
        if sno_long==9222  %sett inn en dummy 313-kanal i CIE-faktorene:
            temp=w_CIE(2:4);
            w_CIE(2)=0;
            w_CIE(3:5)=temp; clear temp
            
            temp=w_uvb(2:4);
            w_uvb(2)=0;
            w_uvb(3:5)=temp; clear temp
            
            temp=w_uva(2:4);
            w_uva(2)=0;
            w_uva(3:5)=temp; clear temp
            
            temp=w_DNA(2:4);
            w_DNA(2)=0;
            w_DNA(3:5)=temp; clear temp
            
            temp=w_vitD(2:4);
            w_vitD(2)=0;
            w_vitD(3:5)=temp; clear temp
            
            temp=w_anchovy(2:4);
            w_anchovy(2)=0;
            w_anchovy(3:5)=temp; clear temp
            
            temp=w_flintplant(2:4);
            w_flintplant(2)=0;
            w_flintplant(3:5)=temp; clear temp
            
            temp=w_HEP(2:4);
            w_HEP(2)=0;
            w_HEP(3:5)=temp; clear temp
            
            temp=w_PAR(2:4);
            w_PAR(2)=0;
            w_PAR(3:5)=temp; clear temp
            
            temp=w_CIE1998(2:4);
            w_CIE1998(2)=0;
            w_CIE1998(3:5)=temp; clear temp
            
            temp=w_NMSC(2:4);
            w_NMSC(2)=0;
            w_NMSC(3:5)=temp; clear temp
            
        end
        
        w_CIE=w_CIE.*Norm; %CIE fkatorer normalisert til FARIN dag 147 2005
        w_uvb=w_uvb.*Norm;
        w_uva=w_uva.*Norm;
        w_DNA=w_DNA.*Norm;
        w_vitD=w_vitD.*Norm;
        w_anchovy=w_anchovy.*Norm;
        w_flintplant=w_flintplant.*Norm;
        w_HEP=w_HEP.*Norm;
        w_PAR=w_PAR.*Norm;
        w_CIE1998=w_CIE1998.*Norm;
        w_NMSC=w_NMSC.*Norm;

        epsilon_cie=ANGF.UVI_SZA_corr;  %Faktorer å dividere med. Legger først til 180 grader
        epsilon_uvb=ANGF.UVB_SZA_corr;
        epsilon_uva=ANGF.UVA_SZA_corr;
        epsilon_DNA=ANGF.DNA_SZA_corr;
        epsilon_vitD=ANGF.vitD_SZA_corr;
        epsilon_anchovy=ANGF.anchovy_SZA_corr;
        epsilon_flintplant=ANGF.flintplant_SZA_corr;
        epsilon_HEP=ANGF.HEP_SZA_corr;
        epsilon_PAR=ANGF.PAR_SZA_corr;
        epsilon_cie1998=ANGF.CIE1998_SZA_corr;
        epsilon_nmsc=ANGF.NMSC_SZA_corr;
        %end
        
        %ekstrapoler til SZA = 180:
        L=size(epsilon_cie,1);
        epsilon_cie(L+1,:)=[180 epsilon_cie(L,2)];
        epsilon_uvb(L+1,:)=[180 epsilon_uvb(L,2)];
        epsilon_uva(L+1,:)=[180 epsilon_uva(L,2)];
        epsilon_DNA(L+1,:)=[180 epsilon_DNA(L,2)];
        epsilon_vitD(L+1,:)=[180 epsilon_vitD(L,2)];
        epsilon_anchovy(L+1,:)=[180 epsilon_anchovy(L,2)];
        epsilon_flintplant(L+1,:)=[180 epsilon_flintplant(L,2)];
        epsilon_HEP(L+1,:)=[180 epsilon_HEP(L,2)];
        epsilon_PAR(L+1,:)=[180 epsilon_PAR(L,2)];
        epsilon_cie1998(L+1,:)=[180 epsilon_cie1998(L,2)];
        epsilon_nmsc(L+1,:)=[180 epsilon_nmsc(L,2)];
        
        [nx,ny]=size(RES.Raw_offset_corr(id_pos,:));
        epsilon_cie_corr=interp1(epsilon_cie(:,1),epsilon_cie(:,2),RES.SZA(id_pos),'linear');
        epsilon_uvb_corr=interp1(epsilon_uvb(:,1),epsilon_uvb(:,2),RES.SZA(id_pos),'linear');
        epsilon_uva_corr=interp1(epsilon_uva(:,1),epsilon_uva(:,2),RES.SZA(id_pos),'linear');
        
        epsilon_DNA_corr=interp1(epsilon_DNA(:,1),epsilon_DNA(:,2),RES.SZA(id_pos),'linear');
        epsilon_vitD_corr=interp1(epsilon_vitD(:,1),epsilon_vitD(:,2),RES.SZA(id_pos),'linear');
        epsilon_anchovy_corr=interp1(epsilon_anchovy(:,1),epsilon_anchovy(:,2),RES.SZA(id_pos),'linear');
        epsilon_flintplant_corr=interp1(epsilon_flintplant(:,1),epsilon_flintplant(:,2),RES.SZA(id_pos),'linear');
        epsilon_HEP_corr=interp1(epsilon_HEP(:,1),epsilon_HEP(:,2),RES.SZA(id_pos),'linear');
        epsilon_PAR_corr=interp1(epsilon_PAR(:,1),epsilon_PAR(:,2),RES.SZA(id_pos),'linear');
        epsilon_cie1998_corr=interp1(epsilon_cie1998(:,1),epsilon_cie1998(:,2),RES.SZA(id_pos),'linear');
        epsilon_nmsc_corr=interp1(epsilon_nmsc(:,1),epsilon_nmsc(:,2),RES.SZA(id_pos),'linear');

        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_cie_corr)*w_CIE(k);
        end
        RES.CIE_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_uvb_corr)*w_uvb(k);
        end
        RES.UVB_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_uva_corr)*w_uva(k);
        end
        RES.UVA_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_DNA_corr)*w_DNA(k);
        end
        RES.DNA_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_vitD_corr)*w_vitD(k);
        end
        RES.VITD_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_anchovy_corr)*w_anchovy(k);
        end
        RES.ANCHOVY_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_flintplant_corr)*w_flintplant(k);
        end
        RES.FLINTPLANT_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_HEP_corr)*w_HEP(k);
        end
        RES.HEP_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_PAR_corr)*w_PAR(k);
        end
        RES.PAR_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_cie1998_corr)*w_CIE1998(k);
        end
        RES.CIE1998_FARIN(id_pos,1)=sum(A,2);clear A
        
        A=zeros(nx,ny);
        for k=1:ny
            A(:,k)=(RES.Raw_offset_corr(id_pos,k)./epsilon_nmsc_corr)*w_NMSC(k);
        end
        RES.NMSC_FARIN(id_pos,1)=sum(A,2);clear A
        
        
        %Gjenta for perioden hvor andre koeffisienter skulle vært benyttet:
        %switch sno_long
        run_again=0;
        if sno_long==9276 && year==2016
            pd=find(RES.Daynums(id_pos)<=76); %%dag 1-76 2016 hadde 9276 defekt 320 - kanal, bruk alternative koeffiseinter uten C3
            run_again=1;
        elseif sno_long==9276 && year==2015
            pd=find(RES.Daynums(id_pos)>312); %%dag 313 2015 hadde 9276 defekt 320 - kanal, bruk alternative koeffiseinter uten C3
            run_again=1;
        end
        
        if run_again==1
            w_CIE=ANGF.C3_DEFECT.CIE*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
            w_uvb=ANGF.C3_DEFECT.uvb_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
            w_uva=ANGF.C3_DEFECT.uva_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
            w_DNA=ANGF.C3_DEFECT.DNA_coeffs*100;   %*100 fordi ANGF er beregnet for rådata microW/cm^2/nm, mens databasen har W/m^2/nm
            w_vitD=ANGF.C3_DEFECT.vitD_coeffs*100;
            w_anchovy=ANGF.C3_DEFECT.anchovy_coeffs*100;
            w_flintplant=ANGF.C3_DEFECT.flintplant_coeffs*100;
            w_HEP=ANGF.C3_DEFECT.HEP_coeffs*100;
            w_PAR=ANGF.C3_DEFECT.PAR_coeffs*100;
            w_CIE1998=ANGF.C3_DEFECT.CIE1998_coeffs*100;
            w_NMSC=ANGF.C3_DEFECT.NMSC_coeffs*100;
            
            epsilon_cie=ANGF.C3_DEFECT.UVI_SZA_corr;  %Faktorer å dividere med. Legger først til 180 grader
            epsilon_uvb=ANGF.C3_DEFECT.UVB_SZA_corr;
            epsilon_uva=ANGF.C3_DEFECT.UVA_SZA_corr;
            epsilon_DNA=ANGF.C3_DEFECT.DNA_SZA_corr;
            epsilon_vitD=ANGF.C3_DEFECT.vitD_SZA_corr;
            epsilon_anchovy=ANGF.C3_DEFECT.anchovy_SZA_corr;
            epsilon_flintplant=ANGF.C3_DEFECT.flintplant_SZA_corr;
            epsilon_HEP=ANGF.C3_DEFECT.HEP_SZA_corr;
            epsilon_PAR=ANGF.C3_DEFECT.PAR_SZA_corr;
            epsilon_cie1998=ANGF.C3_DEFECT.CIE1998_SZA_corr;
            epsilon_nmsc=ANGF.C3_DEFECT.NMSC_SZA_corr;
            
            %ekstrapoler til SZA = 180:
            L=size(epsilon_cie,1);
            epsilon_cie(L+1,:)=[180 epsilon_cie(L,2)];
            epsilon_uvb(L+1,:)=[180 epsilon_uvb(L,2)];
            epsilon_uva(L+1,:)=[180 epsilon_uva(L,2)];
            epsilon_DNA(L+1,:)=[180 epsilon_DNA(L,2)];
            epsilon_vitD(L+1,:)=[180 epsilon_vitD(L,2)];
            epsilon_anchovy(L+1,:)=[180 epsilon_anchovy(L,2)];
            epsilon_flintplant(L+1,:)=[180 epsilon_flintplant(L,2)];
            epsilon_HEP(L+1,:)=[180 epsilon_HEP(L,2)];
            epsilon_PAR(L+1,:)=[180 epsilon_PAR(L,2)];
            epsilon_cie1998(L+1,:)=[180 epsilon_cie1998(L,2)];
            epsilon_nmsc(L+1,:)=[180 epsilon_nmsc(L,2)];
            
            [nx,ny]=size(RES.Raw_offset_corr(id_pos(pd),:));
            epsilon_cie_corr=interp1(epsilon_cie(:,1),epsilon_cie(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_uvb_corr=interp1(epsilon_uvb(:,1),epsilon_uvb(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_uva_corr=interp1(epsilon_uva(:,1),epsilon_uva(:,2),RES.SZA(id_pos(pd)),'linear');
            
            epsilon_DNA_corr=interp1(epsilon_DNA(:,1),epsilon_DNA(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_vitD_corr=interp1(epsilon_vitD(:,1),epsilon_vitD(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_anchovy_corr=interp1(epsilon_anchovy(:,1),epsilon_anchovy(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_flintplant_corr=interp1(epsilon_flintplant(:,1),epsilon_flintplant(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_HEP_corr=interp1(epsilon_HEP(:,1),epsilon_HEP(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_PAR_corr=interp1(epsilon_PAR(:,1),epsilon_PAR(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_cie1998_corr=interp1(epsilon_cie1998(:,1),epsilon_cie1998(:,2),RES.SZA(id_pos(pd)),'linear');
            epsilon_nmsc_corr=interp1(epsilon_nmsc(:,1),epsilon_nmsc(:,2),RES.SZA(id_pos(pd)),'linear');
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_cie_corr)*w_CIE(k);
            end
            RES.CIE_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_uvb_corr)*w_uvb(k);
            end
            RES.UVB_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_uva_corr)*w_uva(k);
            end
            RES.UVA_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_DNA_corr)*w_DNA(k);
            end
            RES.DNA_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_vitD_corr)*w_vitD(k);
            end
            RES.VITD_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_anchovy_corr)*w_anchovy(k);
            end
            RES.ANCHOVY_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_flintplant_corr)*w_flintplant(k);
            end
            RES.FLINTPLANT_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_HEP_corr)*w_HEP(k);
            end
            RES.HEP_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_PAR_corr)*w_PAR(k);
            end
            RES.PAR_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_cie1998_corr)*w_CIE1998(k);
            end
            RES.CIE1998_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
            A=zeros(nx,ny);
            for k=1:ny
                A(:,k)=(RES.Raw_offset_corr(id_pos(pd),k)./epsilon_nmsc_corr)*w_NMSC(k);
            end
            RES.NMSC_FARIN(id_pos(pd),1)=sum(A,2);clear A
            
        end %if run_again==1        
        
        clear ANGF
        %%Metode for å teste Arne Dahlbacks egne beregninger med
        %%Farin skalen:
        %                 if strcmp(stasjon, 'Nyaalesund')
        %                     %beregn UVI med Arne Dahlbacks koeffisienter
        %                     AD_CIE=zeros(nx,1);
        %                     AD_CIEfact=[0.2485*10,0,0.2596,-0.488e-1,0.3934e-1];
        %                     B=zeros(nx,ny);
        %                     for k=1:ny
        %                         B(:,k)=RES.Raw_offset_corr(id_pos,k)*AD_CIEfact(k);
        %                     end
        %                     AD_CIE(id_pos,1)=sum(B,2);clear B
        %                     AD_aarsdose=trapz(RES.JT*24*3600,AD_CIE);
        %                 end
        %
        %                 if strcmp(stasjon,'Andøya')
        %                     %beregn UVI med Arne Dahlbacks koeffisienter
        %                     AD_CIE=zeros(nx,1);
        %                     AD_CIEfact=[0.2580*10,0,0.2330,-0.3736e-1,0.3570e-1];
        %                     B=zeros(nx,ny);
        %                     for k=1:ny
        %                         B(:,k)=RES.Raw_offset_corr(id_pos,k)*AD_CIEfact(k);
        %                     end
        %                     AD_CIE(id_pos,1)=sum(B,2);clear B
        %                     AD_aarsdose=trapz(RES.JT*24*3600,AD_CIE);
        %                 end
        %
        %                 if strcmp(stasjon,'Blindern')
        %                     %beregn UVI med Arne Dahlbacks koeffisienter
        %                     AD_CIE=zeros(nx,1);
        %                     AD_CIEfact=[0.6981*10,0,0.2463,-0.3206e-1,0.4090e-1];
        %                     B=zeros(nx,ny);
        %                     for k=1:ny
        %                         B(:,k)=RES.Raw_offset_corr(id_pos,k)*AD_CIEfact(k);
        %                     end
        %                     AD_CIE(id_pos,1)=sum(B,2);clear B
        %                     AD_aarsdose=trapz(RES.JT*24*3600,AD_CIE);
        %                 end
        
    end     %for i=1:length(SN)  %for alle instrumenter som har stått på stasjonen i denne perioden
    
    %les inn normalverdier for UVI til aktuell stasjon
    %clsky_data=load(sprintf('C:\\Bjorn\\Jobb\\uvnett\\UV\\uvclsky_%03s.txt',stat));
    if strcmp(stasjon,'Nyaalesund')
        %clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt','Nyaal'));
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt',stasjon));
    elseif strcmp(stasjon,'Andøya')
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt','Andoya'));
    elseif strcmp(stasjon,'Østerås')
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt','Oesteraas'));
    elseif strcmp(stasjon,'Kjeller')
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt','Blindern'));
    elseif strcmp(stasjon,'Tromsø')
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt','Tromso'));
    else
        clsky_data=load(sprintf('N:\\uvnet\\guv\\normal_uvi\\normal_uvi_%s.txt',stasjon));
    end
    
    clsky_CIE=[];
    clsky_dailydose=[];
    clsky_maxUVI=[];
    for i=1:366
        clsky_CIE=[clsky_CIE;[i+[0:10/60:24-10/60]/24;clsky_data(i,:)/40]' ]; %regnet om til CIE irradians
        clsky_dailydose=[clsky_dailydose;trapz(  24*3600*(0:10/60:24-10/60)/24,clsky_data(i,:)/40)   ];
        clsky_maxUVI=[clsky_maxUVI;max(clsky_data(i,:))];
    end
    
    RES.CIE_Norm_ClSky_Station=clsky_CIE;
    RES.maxUVI_Norm_ClSky_Station=clsky_maxUVI;
    RES.DailyDose_Norm_ClSky_Station=clsky_dailydose;
    RES.YearlySum_Norm_ClSky_Station=sum(RES.DailyDose_Norm_ClSky_Station);
    %            RES.YearlySum_Norm_ClSky_Station=trapz(RES.CIE_Norm_ClSky_Station(:,1)*24*3600,RES.CIE_Norm_ClSky_Station(:,2));
    
    %Beregning av daily max UVI, dose og yearlysum for GUV målinger:
    
    if (find(RES.JD_unique==366) & daysyear==365)
        RES.JD_unique(RES.JD_unique==366)=[];
    end
    
    %fjern fra RES.JD_unique dager med færre samplinger enn 2:
    q=[];
    for i=1:length(RES.JD_unique)
        daynum=RES.JD_unique(i);
        p=find(RES.Daynums==daynum);
        if numel(p)==1
            q=[q;i];
        end
    end
    if ~isempty(q)
        RES.JD_unique(q)=[];
    end
    
    daydose_cie=[];
    daydose_cie_DB=[];
    daydose_uvb=[];
    daydose_uva=[];
    daydose_DNA=[];
    daydose_vitD=[];
    daydose_anchovy=[];
    daydose_flintplant=[];
    daydose_HEP=[];
    daydose_PAR=[];
    daydose_CIE1998=[];
    daydose_NMSC=[];
    
    max_UVI_DB=[];
    max_UVI_FARIN=[];
    max_UVB_FARIN=[];
    max_UVA_FARIN=[];
    max_DNA_FARIN=[];
    max_VITD_FARIN=[];
    max_ANCHOVY_FARIN=[];
    max_FLINTPLANT_FARIN=[];
    max_HEP_FARIN=[];
    max_PAR_FARIN=[];
    max_UVI1998_FARIN=[];
    max_NMSC_FARIN=[];
    
    for i=1:length(RES.JD_unique)
        daynum=RES.JD_unique(i)
        p=find(RES.Daynums==daynum);
        daydose_cie=[daydose_cie;trapz((RES.JT(p)-daynum)*24*3600,RES.CIE_FARIN(p) )];
        daydose_cie_DB=[ daydose_cie_DB;trapz((RES.JT(p)-daynum)*24*3600,RES.CIEDB(p) )];
        daydose_uvb=[daydose_uvb;trapz((RES.JT(p)-daynum)*24*3600,RES.UVB_FARIN(p) )];
        daydose_uva=[daydose_uva;trapz((RES.JT(p)-daynum)*24*3600,RES.UVA_FARIN(p) )];
        daydose_DNA=[daydose_DNA;trapz((RES.JT(p)-daynum)*24*3600,RES.DNA_FARIN(p) )];
        daydose_vitD=[daydose_vitD;trapz((RES.JT(p)-daynum)*24*3600,RES.VITD_FARIN(p) )];
        daydose_anchovy=[daydose_anchovy;trapz((RES.JT(p)-daynum)*24*3600,RES.ANCHOVY_FARIN(p) )];
        daydose_flintplant=[daydose_flintplant;trapz((RES.JT(p)-daynum)*24*3600,RES.FLINTPLANT_FARIN(p) )];
        daydose_HEP=[daydose_HEP;trapz((RES.JT(p)-daynum)*24*3600,RES.HEP_FARIN(p) )];
        daydose_PAR=[daydose_PAR;trapz((RES.JT(p)-daynum)*24*3600,RES.PAR_FARIN(p) )];
        daydose_CIE1998=[daydose_CIE1998;trapz((RES.JT(p)-daynum)*24*3600,RES.CIE1998_FARIN(p) )];
        daydose_NMSC=[daydose_NMSC;trapz((RES.JT(p)-daynum)*24*3600,RES.NMSC_FARIN(p) )];
        
        max_UVI_FARIN=[max_UVI_FARIN;max(RES.CIE_FARIN(p))*40];
        max_UVI_DB=[max_UVI_DB;max(RES.CIEDB(p))*40];
        max_UVB_FARIN=[max_UVB_FARIN;max(RES.UVB_FARIN(p))];
        max_UVA_FARIN=[max_UVA_FARIN;max(RES.UVA_FARIN(p))];
        max_DNA_FARIN=[max_DNA_FARIN;max(RES.DNA_FARIN(p))];
        max_VITD_FARIN=[max_VITD_FARIN;max(RES.VITD_FARIN(p))];
        max_ANCHOVY_FARIN=[max_ANCHOVY_FARIN;max(RES.ANCHOVY_FARIN(p))];
        max_FLINTPLANT_FARIN=[max_FLINTPLANT_FARIN;max(RES.FLINTPLANT_FARIN(p))];
        max_HEP_FARIN=[max_HEP_FARIN;max(RES.HEP_FARIN(p))];
        max_PAR_FARIN=[max_PAR_FARIN;max(RES.PAR_FARIN(p))];
        max_UVI1998_FARIN=[max_UVI1998_FARIN;max(RES.CIE1998_FARIN(p))*40];
        max_NMSC_FARIN=[max_NMSC_FARIN;max(RES.NMSC_FARIN(p))];
    end
    
    %RES.Daynums_year=(1:daysyear)';
    
    RES.MaxUVI_DB=zeros(size(RES.Daynums_year));%nullstiller alle felter
    RES.MaxUVI_FARIN=RES.MaxUVI_DB;
    RES.MaxUVB_FARIN=RES.MaxUVI_DB;
    RES.MaxUVA_FARIN=RES.MaxUVI_DB;
    RES.MaxDNA_FARIN=RES.MaxUVI_DB;
    RES.MaxVITD_FARIN=RES.MaxUVI_DB;
    RES.MaxANCHOVY_FARIN=RES.MaxUVI_DB;
    RES.MaxFLINTPLANT_FARIN=RES.MaxUVI_DB;
    RES.MaxHEP_FARIN=RES.MaxUVI_DB;
    RES.MaxPAR_FARIN=RES.MaxUVI_DB;
    RES.MaxUVI1998_FARIN=RES.MaxUVI_DB;
    RES.MaxNMSC_FARIN=RES.MaxUVI_DB;
    
    RES.MaxUVI_DB(RES.JD_unique)=max_UVI_DB;
    RES.MaxUVI_FARIN(RES.JD_unique)=max_UVI_FARIN;
    RES.MaxUVB_FARIN(RES.JD_unique)=max_UVB_FARIN;
    RES.MaxUVA_FARIN(RES.JD_unique)=max_UVA_FARIN;
    RES.MaxDNA_FARIN(RES.JD_unique)=max_DNA_FARIN;
    RES.MaxVITD_FARIN(RES.JD_unique)=max_VITD_FARIN;
    RES.MaxANCHOVY_FARIN(RES.JD_unique)=max_ANCHOVY_FARIN;
    RES.MaxFLINTPLANT_FARIN(RES.JD_unique)=max_FLINTPLANT_FARIN;
    RES.MaxHEP_FARIN(RES.JD_unique)=max_HEP_FARIN;
    RES.MaxPAR_FARIN(RES.JD_unique)=max_PAR_FARIN;
    RES.MaxUVI1998_FARIN(RES.JD_unique)=max_UVI1998_FARIN;
    RES.MaxNMSC_FARIN(RES.JD_unique)=max_NMSC_FARIN;
    
    
    RES.DailyDose_DB=zeros(size(RES.Daynums_year));%NULLSTILLER
    RES.DailyDose_FARIN=RES.DailyDose_DB;
    RES.DailyDose_UVB_FARIN=RES.DailyDose_DB;
    RES.DailyDose_UVA_FARIN=RES.DailyDose_DB;
    RES.DailyDose_DNA_FARIN=RES.DailyDose_DB;
    RES.DailyDose_VITD_FARIN=RES.DailyDose_DB;
    RES.DailyDose_ANCHOVY_FARIN=RES.DailyDose_DB;
    RES.DailyDose_FLINTPLANT_FARIN=RES.DailyDose_DB;
    RES.DailyDose_HEP_FARIN=RES.DailyDose_DB;
    RES.DailyDose_PAR_FARIN=RES.DailyDose_DB;
    RES.DailyDose_CIE1998_FARIN=RES.DailyDose_DB;
    RES.DailyDose_NMSC_FARIN=RES.DailyDose_DB;
    
    RES.DailyDose_DB(RES.JD_unique)= daydose_cie_DB;
    RES.DailyDose_FARIN(RES.JD_unique)=daydose_cie;
    RES.DailyDose_UVB_FARIN(RES.JD_unique)= daydose_uvb;
    RES.DailyDose_UVA_FARIN(RES.JD_unique)=daydose_uva;
    RES.DailyDose_DNA_FARIN(RES.JD_unique)=daydose_DNA;
    RES.DailyDose_VITD_FARIN(RES.JD_unique)=daydose_vitD;
    RES.DailyDose_ANCHOVY_FARIN(RES.JD_unique)=daydose_anchovy;
    RES.DailyDose_FLINTPLANT_FARIN(RES.JD_unique)=daydose_flintplant;
    RES.DailyDose_HEP_FARIN(RES.JD_unique)=daydose_HEP;
    RES.DailyDose_PAR_FARIN(RES.JD_unique)=daydose_PAR;
    RES.DailyDose_CIE1998_FARIN(RES.JD_unique)=daydose_CIE1998;
    RES.DailyDose_NMSC_FARIN(RES.JD_unique)=daydose_NMSC;
    
    RES.YearlySum_CIEDB=sum(RES.DailyDose_DB);
    RES.YearlySum_CIE_FARIN=sum(RES.DailyDose_FARIN);
    RES.YearlySum_UVB_FARIN=sum(RES.DailyDose_UVB_FARIN);
    RES.YearlySum_UVA_FARIN=sum(RES.DailyDose_UVA_FARIN);
    RES.YearlySum_DNA_FARIN=sum(RES.DailyDose_DNA_FARIN);
    RES.YearlySum_VITD_FARIN=sum(RES.DailyDose_VITD_FARIN);
    RES.YearlySum_ANCHOVY_FARIN=sum(RES.DailyDose_ANCHOVY_FARIN);
    RES.YearlySum_FLINTPLANT_FARIN=sum(RES.DailyDose_FLINTPLANT_FARIN);
    RES.YearlySum_HEP_FARIN=sum(RES.DailyDose_HEP_FARIN);
    RES.YearlySum_PAR_FARIN=sum(RES.DailyDose_PAR_FARIN);
    RES.YearlySum_CIE1998_FARIN=sum(RES.DailyDose_CIE1998_FARIN);
    RES.YearlySum_NMSC_FARIN=sum(RES.DailyDose_NMSC_FARIN);
    
    if run_uncalibrated==1
        RAW=RES;
        clear RES
%         tmp = 'RAW';
%         save(fpath,tmp)
        RES=RAW;clear RAW
        
        %clear RAW
    else
%         tmp = 'RES';
%         save(fpath,tmp)
        %clear RES;
    end
    
end  %if size(RES.Raw,1)>0
%end %        if fidsjekk>0
%Beregn clearsky UVI for dagens totalozon, tilpasset sesongmessig
%avvik mellom modellert og målt UVI på hver stasjon

%Må først hente inn UVI-korreksjon til modellberegninger, ut fra
%ratio MODEL UVI / GUV UVI:

measurement_name=sprintf('N:\\uvnet\\guv\\calc_model_2_meas_uvi\\model_2_meas_%s_%s.mat',stasjon_norsk,prim_sno);

if stasjons_id==1
    %measurement_name=sprintf('N:\\uvnet\\guv\\calc_model_2_meas_uvi\\model_2_meas_%s_combined.mat',stasjon_norsk);
    measurement_name=sprintf('N:\\uvnet\\guv\\calc_model_2_meas_uvi\\model_2_meas_%s_%s.mat',stasjon_norsk,prim_sno);%bruker heller 9277!
elseif stasjons_id==8
    measurement_name=sprintf('N:\\uvnet\\guv\\calc_model_2_meas_uvi\\model_2_meas_%s_%s.mat','Nyaal',prim_sno);
elseif stasjons_id==11 %Kjeller
    measurement_name=sprintf('N:\\uvnet\\guv\\calc_model_2_meas_uvi\\model_2_meas_%s_%s.mat','Blindern',prim_sno);
end
GUV=load(measurement_name); %Gir struct GUV.MODEL med innholdet i struct MODEL
Climat_mean_UVI_model2GUV=GUV.RES.Climat_mean_UVI_model2GUV;
Climat_mean_Dose_model2GUV=GUV.RES.Climat_mean_Dose_model2GUV;

RES.MODEL.Climat_mean_UVI_model2GUV=GUV.RES.Climat_mean_UVI_model2GUV;
RES.MODEL.Climat_mean_Dose_model2GUV=GUV.RES.Climat_mean_Dose_model2GUV;

clear GUV

%Beregn clearsky irradianser - for statisk sommeralbedo 5%, variabel SZA og ozon -
%korriger deretter UVI etter climat_mean UVI:
%         if year>=1995
%             norm_ozone=0;   %skal bruke dagens ozone, fra satelitt data
%         else
%             norm_ozone=1;   %bruker seasonal mean ozon for Oslo
%         end

if year>=1900
    norm_ozone=0;   %skal alltid bruke allerede beregnet ozon som ligger i RES.OMI.O3
else
    norm_ozone=1;   %bruker seasonal mean ozon for Oslo
end
[stasjon year]
MODEL=uvspec_calc_yeardata_new(RES.MODEL.O3,year,stasjon_norsk,stat,RES.Latitude,RES.Longitude,RES.MODEL.Climat_mean_UVI_model2GUV,norm_ozone);%best resultat
%MODEL=uvspec_calc_yeardata_new(RES.MODEL.O3,year,stasjon,stat,RES.Latitude,RES.Longitude,RES.MODEL.Climat_mean_Dose_model2GUV,norm_ozone);
%(O3_vector,ozone_guv,year,stasjon,stat,lat,long,norm_ozone)

MODEL.Climat_mean_UVI_model2GUV=RES.MODEL.Climat_mean_UVI_model2GUV;
MODEL.Climat_mean_Dose_model2GUV=RES.MODEL.Climat_mean_Dose_model2GUV;

MODEL.JT=MODEL.Daynum+MODEL.HTime/24;

RES.MODEL=MODEL;

%         figure(year)
%         hold on
%         plot(RES.MODEL.JT,RES.MODEL.UVI,'g-')
%         hold off
%         legend('UVI-DB','UVI-Farin','Clearsky model')

%RES.YearlySum_ModelClSky_Station=trapz(RES.MODEL.JT*24*3600,RES.MODEL.UVI/40);
RES.YearlySum_ModelClSky_Station=sum(RES.MODEL.DailyDose);

%         tmp = 'RES';
%         save(fpath,tmp)

if run_uncalibrated==1
    RAW=RES;
%     clear RES
%     tmp = 'RAW';
%     save(fpath,tmp)
    RES=RAW;clear RAW
    
    %clear RAW
else
%     tmp = 'RES';
%     save(fpath,tmp)
    %clear RES;
end

if isfield(RES,'YearlySum_CIEDB')
    annual=[annual;year,RES.YearlySum_ModelClSky_Station,RES.YearlySum_Norm_ClSky_Station,RES.YearlySum_CIEDB, RES.YearlySum_CIE_FARIN];
else
    annual=[annual;year,0,RES.YearlySum_ModelClSky_Station, 0,0];
end

%         figure(year*10)
%         subplot(2,1,1)
%         plot((1:length(RES.MODEL.Max_UVI)),RES.MODEL.Max_UVI,'m-'),grid on
%         hold on
%         if isfield(RES,'MaxUVI_DB')
%             plot((1:length(RES.maxUVI_Norm_ClSky_Station)),RES.maxUVI_Norm_ClSky_Station,'r-')
%             plot(RES.Daynums_year,RES.MaxUVI_DB,'b.-')
%             plot(RES.Daynums_year,RES.MaxUVI_FARIN,'k.-')
%         end
%         hold off
%         ylabel(sprintf('Daily UVI %s %i',stasjon,year))
%         if isfield(RES,'MaxUVI_DB')
%             legend('Model clearsky','Norm clearsky','DB','FARIN')
%             y2=max([ max(RES.maxUVI_Norm_ClSky_Station) max(RES.MODEL.Max_UVI) max(RES.MaxUVI_DB) max(RES.MaxUVI_FARIN)]);
%         else
%             legend('Model clearsky')
%             y2=max([max(RES.MODEL.Max_UVI)]);
%         end
%         y1=-0.1;
%         axis([1 366 y1 y2*1.1 ])
%
%         subplot(2,1,2)
%         plot((1:length(RES.MODEL.DailyDose)),RES.MODEL.DailyDose,'m-'), grid on
%         hold on
%         if isfield(RES,'DailyDose_DB')
%             plot((1:length(RES.DailyDose_Norm_ClSky_Station)),RES.DailyDose_Norm_ClSky_Station,'r-')
%             plot(RES.Daynums_year,RES.DailyDose_DB,'b.-')
%             plot(RES.Daynums_year,RES.DailyDose_FARIN,'k.-')
%         end
%         hold off
%         ylabel(sprintf('Daily CIE dose %s %i',stasjon,year))
%         if isfield(RES,'MaxUVI_DB')
%             legend('Model clearsky','Norm clearsky','DB','FARIN')
%             y2=max([max(RES.MODEL.DailyDose),max(RES.DailyDose_Norm_ClSky_Station),max(RES.DailyDose_DB),max(RES.DailyDose_FARIN)]);
%         else
%             legend('Model clearsky')
%             y2=max(max(RES.MODEL.DailyDose));
%         end
%         y1=-100;
%         axis([1 366 y1 y2*1.1 ])

if ~isempty(RES.Daynums)
    
    figure(1)
    plot(RES.Daynums_year,RES.MaxUVI_FARIN,'k.-'),grid on,title(sprintf('Max UVI %s %d',stat,year))
    figure(2)
    plot(RES.Daynums_year,RES.MaxUVB_FARIN,'k.-'),grid on,title(sprintf('Max UVB %s %d',stat,year))
    figure(3)
    plot(RES.Daynums_year,RES.MaxUVA_FARIN,'k.-'),grid on,title(sprintf('Max UVA %s %d',stat,year))
    figure(4)
    plot(RES.Daynums_year,RES.MaxDNA_FARIN,'k.-'),grid on,title(sprintf('Max DNA %s %d',stat,year))
    figure(5)
    plot(RES.Daynums_year,RES.MaxVITD_FARIN,'k.-'),grid on,title(sprintf('Max VITD %s %d',stat,year))
    figure(6)
    plot(RES.Daynums_year,RES.MaxANCHOVY_FARIN,'k.-'),grid on,title(sprintf('Max ANCHOVY EGG & LARVAE %s %d',stat,year))
    figure(7)
    plot(RES.Daynums_year,RES.MaxFLINTPLANT_FARIN,'k.-'),grid on,title(sprintf('Max FLINT & CALDWELL PLANT %s %d',stat,year))
    figure(8)
    plot(RES.Daynums_year,RES.MaxHEP_FARIN,'k.-'),grid on,title(sprintf('Max HEP %s %d',stat,year))
    figure(9)
    plot(RES.Daynums_year,RES.MaxPAR_FARIN,'k.-'),grid on,title(sprintf('Max PAR %s %d',stat,year))
    
    figure(10)
    subplot(3,4,1)
    plot(RES.Daynums_year,RES.MaxUVI_FARIN,'k.-'),grid on,title(sprintf('Max UVI %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,2)
    plot(RES.Daynums_year,RES.MaxUVB_FARIN,'k.-'),grid on,title(sprintf('Max UVB %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,3)
    plot(RES.Daynums_year,RES.MaxUVA_FARIN,'k.-'),grid on,title(sprintf('Max UVA %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,4)
    plot(RES.Daynums_year,RES.MaxDNA_FARIN,'k.-'),grid on,title(sprintf('Max DNA %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,5)
    plot(RES.Daynums_year,RES.MaxVITD_FARIN,'k.-'),grid on,title(sprintf('Max VITD %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,6)
    plot(RES.Daynums_year,RES.MaxANCHOVY_FARIN,'k.-'),grid on,title(sprintf('Max ANCHOVY EGG & LARVAE %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,7)
    plot(RES.Daynums_year,RES.MaxFLINTPLANT_FARIN,'k.-'),grid on,title(sprintf('Max FLINT & CALDWELL PLANT %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,8)
    plot(RES.Daynums_year,RES.MaxHEP_FARIN,'k.-'),grid on,title(sprintf('Max HEP %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,9)
    plot(RES.Daynums_year,RES.MaxPAR_FARIN,'k.-'),grid on,title(sprintf('Max PAR %s %d',stat,year))
    ylabel('mole/s/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,10)
    plot(RES.Daynums_year,RES.MaxUVI1998_FARIN,'k.-'),grid on,title(sprintf('Max UVI1998 %s %d',stat,year))
    ylabel('mole/s/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,11)
    plot(RES.Daynums_year,RES.MaxNMSC_FARIN,'k.-'),grid on,title(sprintf('Max NMSC %s %d',stat,year))
    ylabel('mole/s/m^2')
    set(gca,'FontSize',8)
    %export_fig('Z:\Seksjon laboratorier\Optisk Lab\Uvnet\Prosjekter\Porfyri_prosjektet\figurer\Data_products_max','-png');
    
    
    figure(11)%dagsdoser
    subplot(3,4,1)
    plot(RES.Daynums_year,RES.DailyDose_FARIN,'k.-'),grid on,title(sprintf('Daily CIE %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,2)
    plot(RES.Daynums_year,RES.DailyDose_UVB_FARIN,'k.-'),grid on,title(sprintf('Daily UVB %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,3)
    plot(RES.Daynums_year,RES.DailyDose_UVA_FARIN,'k.-'),grid on,title(sprintf('Daily UVA %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,4)
    plot(RES.Daynums_year,RES.DailyDose_DNA_FARIN,'k.-'),grid on,title(sprintf('Daily DNA %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,5)
    plot(RES.Daynums_year,RES.DailyDose_VITD_FARIN,'k.-'),grid on,title(sprintf('Daily VITD %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,6)
    plot(RES.Daynums_year,RES.DailyDose_ANCHOVY_FARIN,'k.-'),grid on,title(sprintf('Daily ANCHOVY EGG & LARVAE %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,7)
    plot(RES.Daynums_year,RES.DailyDose_FLINTPLANT_FARIN,'k.-'),grid on,title(sprintf('Daily FLINT & CALDWELL PLANT %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,8)
    plot(RES.Daynums_year,RES.DailyDose_HEP_FARIN,'k.-'),grid on,title(sprintf('Daily HEP %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,9)
    plot(RES.Daynums_year,RES.DailyDose_PAR_FARIN,'k.-'),grid on,title(sprintf('Daily PAR %s %d',stat,year))
    ylabel('mole/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,10)
    plot(RES.Daynums_year,RES.DailyDose_CIE1998_FARIN,'k.-'),grid on,title(sprintf('Daily CIE1998 %s %d',stat,year))
    ylabel('mole/m^2')
    set(gca,'FontSize',8)
    subplot(3,4,11)
    plot(RES.Daynums_year,RES.DailyDose_NMSC_FARIN,'k.-'),grid on,title(sprintf('Daily NMSC %s %d',stat,year))
    ylabel('mole/m^2')
    set(gca,'FontSize',8)
    %export_fig('Z:\Seksjon laboratorier\Optisk Lab\Uvnet\Prosjekter\Porfyri_prosjektet\figurer\Data_products_doses','-png');
    
    %modellert klarvær mot målt: Max-verdier hver dag:
    
    if stasjons_id==5%Blindern - sammenlign PAR-kanal med proxy PAR
        %Om PAR-data ikke allerede er gjort klar, kjør: porfyri_GUV_load_PAR_data.m
        PAR=load(sprintf('D:\\Data\\work\\uvnett\\PAR\\%s\\PAR_%s_%04i.mat',stat,stat,year));
        undn=unique(floor(PAR.RES.JD));
        if ~isempty(undn)
            V=zeros(numel(undn),3);
            for pp=1:numel(undn)
                kk=find(floor(PAR.RES.JD)==pp);
                if~isempty(kk)
                    pp
                    V(pp,1)=pp;
                    V(pp,2)=max(PAR.RES.PAR(kk))/100;%enheter mol/m^2/s
                    V(pp,3)=trapz((PAR.RES.JD(kk) -pp)*24*3600,PAR.RES.PAR(kk))/100;
                end
            end
        end
        qq=find(V(:,1)==0);
        if ~isempty(qq)
            V(qq,:)=[];
        end
    end
    
    figure(12)
    subplot(3,4,1)
    plot(RES.Daynums_year,RES.MaxUVI_FARIN,'k.-'),grid on,title(sprintf('UVI %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_UVI,'r-')
    hold off
    ylabel('UVI')
    set(gca,'FontSize',8)
    
    subplot(3,4,2)
    plot(RES.Daynums_year,RES.MaxUVB_FARIN,'k.-'),grid on,title(sprintf('UVB %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_UVB,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,3)
    plot(RES.Daynums_year,RES.MaxUVA_FARIN,'k.-'),grid on,title(sprintf('UVA %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_UVA,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,4)
    plot(RES.Daynums_year,RES.MaxDNA_FARIN,'k.-'),grid on,title(sprintf('DNA %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_DNA,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,5)%reservert ACGIH
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_ACGIH,'r-'),grid on,title(sprintf('ACGIH %s %d',stat,year))
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,6)
    plot(RES.Daynums_year,RES.MaxVITD_FARIN,'k.-'),grid on,title(sprintf('VITD %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_VITD,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,7)
    plot(RES.Daynums_year,RES.MaxANCHOVY_FARIN,'k.-'),grid on,title(sprintf('ANCHOVY EGG & LARVAE %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_ANCHOVY,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,8)
    plot(RES.Daynums_year,RES.MaxFLINTPLANT_FARIN,'k.-'),grid on,title(sprintf('FLINT & CALDWELL PLANT %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_FLINTPLANT,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,9)
    plot(RES.Daynums_year,RES.MaxHEP_FARIN,'k.-'),grid on,title(sprintf('HEP %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_HEP,'r-')
    hold off
    ylabel('W/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,10)
    plot(RES.Daynums_year,RES.MaxPAR_FARIN,'k.-'),grid on,title(sprintf('PAR %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_PAR,'r-')
    if stasjons_id==5
        plot(V(:,1),V(:,2),'b-')
    end
    hold off
    ylabel('mole/s/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,11)
    plot(RES.Daynums_year,RES.MaxUVI1998_FARIN,'k.-'),grid on,title(sprintf('UVI1998 %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_UVI1998,'r-')
    if stasjons_id==5
        plot(V(:,1),V(:,2),'b-')
    end
    hold off
    ylabel('UVI')
    set(gca,'FontSize',8)
    
    subplot(3,4,12)
    plot(RES.Daynums_year,RES.MaxNMSC_FARIN,'k.-'),grid on,title(sprintf('NMSC %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.Max_NMSC,'r-')
    if stasjons_id==5
        plot(V(:,1),V(:,2),'b-')
    end
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    figure(12)
    set(gcf, 'Position', get(0, 'Screensize'));%fullscreen size - for best visning av graf
    pause(2)%må legge inn en ventepause for at figurvinduet skal rekke å maksimaliseres før lagring
    gpath = sprintf('%sfigures\\DailyMax\\%sDailyMax_%s_%04i.fig',proot,prename,stasjon_norsk,year);
    savefig(gpath)
    
    %sjekk modellert klarvær mot faktisk målt: dagsdoser
    figure(13)
    subplot(3,4,1)
    plot(RES.Daynums_year,RES.DailyDose_FARIN,'k.-'),grid on,title(sprintf('CIE %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,2)
    plot(RES.Daynums_year,RES.DailyDose_UVB_FARIN,'k.-'),grid on,title(sprintf('UVB %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_UVB,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,3)
    plot(RES.Daynums_year,RES.DailyDose_UVA_FARIN,'k.-'),grid on,title(sprintf('UVA %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_UVA,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,4)
    plot(RES.Daynums_year,RES.DailyDose_DNA_FARIN,'k.-'),grid on,title(sprintf('DNA %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_DNA,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,5)%reservert ACGIH
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_ACGIH,'r-'),grid on,title(sprintf('ACGIH %s %d',stat,year))
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    
    subplot(3,4,6)
    plot(RES.Daynums_year,RES.DailyDose_VITD_FARIN,'k.-'),grid on,title(sprintf('VITD %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_VITD,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,7)
    plot(RES.Daynums_year,RES.DailyDose_ANCHOVY_FARIN,'k.-'),grid on,title(sprintf('ANCHOVY EGG & LARVAE %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_ANCHOVY,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,8)
    plot(RES.Daynums_year,RES.DailyDose_FLINTPLANT_FARIN,'k.-'),grid on,title(sprintf('FLINT & CALDWELL PLANT %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_FLINTPLANT,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,9)
    plot(RES.Daynums_year,RES.DailyDose_HEP_FARIN,'k.-'),grid on,title(sprintf('HEP %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_HEP,'r-')
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,10)
    plot(RES.Daynums_year,RES.DailyDose_PAR_FARIN,'k.-'),grid on,title(sprintf('PAR %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_PAR,'r-')
    
    if stasjons_id==5
        plot(V(:,1),V(:,3),'b-')
    end
    
    hold off
    ylabel('mole/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,11)
    plot(RES.Daynums_year,RES.DailyDose_CIE1998_FARIN,'k.-'),grid on,title(sprintf('CIE1998 %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_CIE1998,'r-')
    
    if stasjons_id==5
        plot(V(:,1),V(:,3),'b-')
    end
    
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    subplot(3,4,12)
    plot(RES.Daynums_year,RES.DailyDose_NMSC_FARIN,'k.-'),grid on,title(sprintf('NMSC %s %d',stat,year))
    hold on
    plot(RES.MODEL.Day_vector,RES.MODEL.DailyDose_NMSC,'r-')
    
    if stasjons_id==5
        plot(V(:,1),V(:,3),'b-')
    end
    
    hold off
    ylabel('J/m^2')
    set(gca,'FontSize',8)
    
    figure(13)
    set(gcf, 'Position', get(0, 'Screensize'));%fullscreen size - for best visning av graf
    pause(2)
    %hFig=figure(13);
%     gpath = sprintf('%sfigures\\DailyDoses\\%sDailyDoses_%s_%04i.fig',proot,prename,stasjon_norsk,year);
%     savefig(gpath)
%     %hgexport(hFig,'-clipboard')
%     
%     h=gca;
%     saveas(gca,'test.svg');%vector graphics    
    
    figure(2000)
    subplot(2,1,1)
    plot(RES.Daynums_year,RES.MODEL.Max_UVI./RES.MaxUVI_FARIN,'k.-'),grid on,title(sprintf('UVI %s %d',stat,year))
    ylabel('UVI model/meas')
    ylim([0.5 1.5])
    
    subplot(2,1,2)
    plot(RES.Daynums_year,RES.MODEL.DailyDose./RES.DailyDose_FARIN,'k.-'),grid on,title(sprintf('Dose %s %d',stat,year))
    ylabel('dose model/meas')
    ylim([0.5 1.5])
    
end %if isfield(RES,'Daynums_year')


