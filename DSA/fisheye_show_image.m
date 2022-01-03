function fisheye_show_image
%Beregner skyview faktor for fisheye fotos.
%kamera: circular fisheye, equidistant (polar angle) projection: punkter med samme
%polarvinkel avtegnes på en sirkel i bildeplanet
%Tidsangivelse er norsk sommertid (UTC+2)

%{
From line 138 - 249 code has been changed to introduce separate bw-filtration.
Line 794 - 810 also changed to make meshgrids so the plots can be surfaced over year and day
Line 926 - 1032 removed indexation from
porfyri_returner_horisontalkomponent to properly extract data for the
entire year. 
Some extra code added to surfaceplot the various results. Also contruction
of a struct so the plots can easily be remade without the need of running
entire script. - Jonny (2021)
%}

 
system('echo off | clip') %tømmer windows clipboardet

%kalibreringsbilde DSCN*.JPG: nikon coolpix 5400 - beregner bildesentrum, radius ut til sza 90 grader og
%radius ut til bilde kanten

outpath='C:\Users\jonny\Pictures\DSA\Bearbeidet';%output figs lagres her
im_cal='C:\Users\jonny\Pictures\DSA\DSCN2333.jpg';
rgb = imread(im_cal);
[im_cntrx,im_cntry,im_rad,im_radH]=fisheye_find_center_radius(rgb);
close(1),close(2),close(6)

%imfile='Z:\Seksjon Overvåkning og miljødata\imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2336.jpg';83054 overvåkning UV og soling\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2335.jpg';
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2334.jpg';
% imfile='C:\Users\jonny\Pictures\DSA\unprocessed\DSCN2430.JPG';limpix=150; % Kolbotn bhage - første bilde under solseil
% imfile='C:\Users\jonny\Pictures\DSA\unprocessed\DSCN2431.JPG';limpix=150;
% imfile='C:\Users\jonny\Pictures\DSA\unprocessed\DSCN2432.JPG';limpix=150;
imfile='C:\Users\jonny\Pictures\DSA\unprocessed\DSCN2430.JPG';

%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2336.jpg';limpix=130;%uthever lys himmel. limpix=35 uthever både himmel og solseilene, mens husvegger og trær blir mørke;%limpix=130;%limpix=110;%under solseil i Helset barnehag
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2338.jpg';limpix=85;
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2342.jpg';
%%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2349.jpg';limpix=80;%Bryn skole-bra eksempel
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2355.jpg';limpix=107;
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2356.jpg';limpix=90;
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2357.jpg';limpix=90;%Høvik barnehage

%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2350.jpg';limpix=160;%Bryn skole  skur
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2351.jpg';
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2352.jpg';
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2364.jpg';limpix=80;%Sandvikastorsenterplass
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2365.jpg';limpix=150;%Sandvikastorsenter McDonaldsgata
%%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2363.jpg';limpix=90;%Under trekrone-bra bilde

%imfile='Z:\Seksjon Overvåkning og miljødata\83054 overvåkning UV og soling\BILDERogILLUSTRASJONER\Bilder\WWW_photos\Stasjoner\Oesteraas\plattform2004\DSCN0021 - Oesteraas.jpg';%fisheye bilde klarvær Østerås
%imfile='Z:\Seksjon Overvåkning og miljødata\83054 overvåkning UV og soling\BILDERogILLUSTRASJONER\Bilder\WWW_photos\Stasjoner\Landvik\2005\DSCN0324.JPG';
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2394_Parasoll.JPG';limpix=100;%Parasoll
%imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2395_Parasoll.JPG';limpix=100;%Parasoll
%parasollbildene viser at en stor parasoll dekker omtrent en polarvinkel på 60 grader, regnet fra hodehøyde. Skyview blir
%da pi*(1-sin^2(60))= 0.25. Andre studier har vist at en middels parasoll dekker rundt

% %bilde hentet fra nettet-tatt med et annet kamera
% imfile='J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\fishey_citysquare_1024px-28_-_New_York_-_Octobre_2008.jpg';limpix=253;

[im_date,filestring,filetype]=file_creation_date(imfile); %2016           4          13          11          38          5

rgb = imread(imfile);
org_rgb=rgb;

if ~strcmp(filestring(1:4),'DSCN') %ikke nikon coolpix strålevernskamera. %anta at bildet allerede er sentrert og croppet til et sirkulært bilde
    im_cntrx = 1+size(rgb,2)/2;
    im_cntry = 1+size(rgb,1)/2;
    im_rad=im_cntrx-1;
    im_radH=im_cntrx-1;
end

%utsnitt, kvadratisk som omkranser sirkulært område:
rect=[im_cntrx-im_rad im_cntry-im_rad 2*im_rad 2*im_rad];
fish_rgb_RH = imcrop(rgb,rect);
rgb=fish_rgb_RH;
fish_bw_RH=rgb2gray(fish_rgb_RH);
M = size(fish_rgb_RH,1);%samme
N = size(fish_rgb_RH,2);

im_cntrx=round(M/2);
im_cntry=im_cntrx;

r=im_rad; % radius ut til bildesirkelen
r90=im_radH;%radius ut til 90-graders horisonten
%Finn piksler innenfor sirkelen:Beregn først avstand fra sentrum ut til hvert bildepunkt:
[Y,X]=meshgrid(1:N,1:M);
R=sqrt((Y-im_cntrx).^2+(X-im_cntry).^2);
pr=find(R<=r90);

% figure(100)
% subplot(3,1,1),ir=fish_rgb_RH(:,:,1);ih=imhist(ir(pr));bar(ih,'r'),axis([10 260 0 max(ih(2:end-1))])
% subplot(3,1,2),ig=fish_rgb_RH(:,:,2);ih=imhist(ig(pr));bar(ih,'g'),axis([10 260 0 max(ih(2:end-1))])
% subplot(3,1,3),ib=fish_rgb_RH(:,:,3);ih=imhist(ib(pr));bar(ih,'b'),axis([10 260 0 max(ih(2:end-1))])


% figure(101)%anbefales å velge terskel 80 percentil
% tr=cumsum(imhist(ir(pr)));
% tg=cumsum(imhist(ig(pr)));
% tb=cumsum(imhist(ib(pr)));
% hold on
% plot(tr/tr(end),'r.-')
% plot(tg/tg(end),'g.-')
% plot(tb/tb(end),'b.-')
% hold off
% grid on


%figure(1)
%imshow(org_rgb);
%imshow(fish_rgb_RH);
% Enlarge figure to full screen.
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

%r90=im_radH;
center=[im_cntrx,im_cntry]; % center cordinates
t=-pi:2*pi/3600:pi;
x=round(r*cos(t)+center(1));
y=round(r*sin(t)+center(2));

xH=round(r90*cos(t)+center(1));
yH=round(r90*sin(t)+center(2));

% figure(2)
% imshow(fish_bw_RH)
% %axis off          % Remove axis ticks and numbers
% axis image
% hold on
% plot(im_cntrx,im_cntry,'wo','MarkerSize',10)
% plot([im_cntrx im_cntrx],[1 M],'Color','w','LineWidth',1)
% plot([1 N],[im_cntry im_cntry],'Color','w','LineWidth',1)
% plot(x,y,'w--','LineWidth',0.5),
% plot(xH,yH,'b--','LineWidth',0.1),
% hold off
%
% figure(3);
% p=find(fish_bw_RH>1);
% imhist(fish_bw_RH(p));

% figure(1)
% hold on
% plot(x,y,'w--','LineWidth',0.5),
% plot(xH,yH,'r--','LineWidth',0.1)
% hold off

%test_im=fish_bw_RH;
test_im=0*fish_bw_RH;
test2_im = test_im;      % This is for the split-image
test3_im = test_im;
test4_im = test_im;
%limpix=160;
%limpix=154;
%limpix=110;
%limpix=107;

sav = fish_rgb_RH;
u = fish_rgb_RH; uu = u; ul = u;        %
l = fish_rgb_RH; lu = l; ll = l;        %
u(:,:,1) = triu(fish_rgb_RH(:,:,1));    % This part will "cut" the image in half, finding the upper 
u(:,:,2) = triu(fish_rgb_RH(:,:,2));    % triangular part of the matrix - to better try to filter 
u(:,:,3) = triu(fish_rgb_RH(:,:,3));    % the image. 

ul(:,:,1) = tril(flipdim(u(:,:,1),2));
ul(:,:,2) = tril(flipdim(u(:,:,2),2));
ul(:,:,3) = tril(flipdim(u(:,:,3),2));

uu(:,:,1) = triu(flipdim(u(:,:,1),2));
uu(:,:,2) = triu(flipdim(u(:,:,2),2));
uu(:,:,3) = triu(flipdim(u(:,:,3),2));

l(:,:,1) = tril(fish_rgb_RH(:,:,1));    % Same, but for the lower triangular part of the matrix
l(:,:,2) = tril(fish_rgb_RH(:,:,2));
l(:,:,3) = tril(fish_rgb_RH(:,:,3));

lu(:,:,1) = triu(flipdim(l(:,:,1),2));
lu(:,:,2) = triu(flipdim(l(:,:,2),2));
lu(:,:,3) = triu(flipdim(l(:,:,3),2));

ll(:,:,1) = tril(flipdim(l(:,:,1),2));
ll(:,:,2) = tril(flipdim(l(:,:,2),2));
ll(:,:,3) = tril(flipdim(l(:,:,3),2));

% % PD = upper piece
% % PD2 = left piece
% % PD3 = right piece
% % PD4 = bottom piece
upper_quadrant = ["> 0", "> 0", "> 245"]; % These parameters are individual for each fisheye photo
left_quadrant = ["< 200", "< 200", "> 200"]; % and they are saved in the IMG.struct for each of the already processed pictures
right_quadrant = ["> 200", "> 200", "> 0"]; % To create new ones, remove the % from line 1010 after tweaking parameters
lower_quadrant = ["> 200", "> 0", "> 0"]; % to best bw-filter the desired picture.

pD = find(eval(join("uu(:,:,3)" + upper_quadrant(3))) & eval(join("uu(:,:,2)" + upper_quadrant(2))) & eval(join("uu(:,:,1)" + upper_quadrant(1))));
pD2 = find(eval(join("ul(:,:,3)" + left_quadrant(3))) & eval(join("ul(:,:,2)" + left_quadrant(2))) & eval(join("ul(:,:,1)" + left_quadrant(1))));
pD3 = find(eval(join("lu(:,:,3)" + right_quadrant(3))) & eval(join("lu(:,:,2)" + right_quadrant(2))) & eval(join("lu(:,:,1)" + right_quadrant(1))));
pD4 = find(eval(join("ll(:,:,3)" + lower_quadrant(3))) & eval(join("ll(:,:,2)" + lower_quadrant(2))) & eval(join("ll(:,:,1)" + lower_quadrant(1))));

% pD = find(uu(:,:,3) > 200 & uu(:,:,2) > 0 & uu(:,:,1) > 0);     % Choose parameteres  for the 4 pieces
% pD2 = find(ul(:,:,3) > 200 & ul(:,:,2) < 210 & ul(:,:,1) < 200);          %  of cake that is the image, making the user
% pD3 = find(lu(:,:,3) > 200 & lu(:,:,2) > 0 & lu(:,:,1) > 0);    % able to handpick more efficiently what rgb-values to
% pD4 = find(ll(:,:,3) > 140 & ll(:,:,2) > 0 & ll(:,:,1) > 0);          % differentiate between sky and non-sky.

if ~strcmp(filestring(1:4),'DSCN') %dette er new-york bildet fra nettet
    %pD=find(fish_rgb_RH(:,:,3)>=252 & fish_rgb_RH(:,:,1)>164);
    %pD=find(fish_rgb_RH(:,:,3)>=252);%blå kanal alene
    %pD=find(fish_rgb_RH(:,:,1)>=170);%rød kanal alene
    %pD=find(fish_rgb_RH(:,:,2)>=216);%grønn kanal alene
    pD=find(fish_rgb_RH(:,:,3)>=252 & fish_rgb_RH(:,:,2)>=216 & fish_rgb_RH(:,:,1)>=170);
end
%pD=find(fish_rgb_RH(:,:,3)>fish_rgb_RH(:,:,1) & fish_rgb_RH(:,:,3)>limpix );
%test_im(pD)=0;%blå kanal: blå himmmel har høy blåkomponent, mens bygninger har lav komponent
test_im(pD)=255;
test2_im(pD2)=255;           % This is for the split-image
test3_im(pD3)=255;
test4_im(pD4)=255;


%pL=find(fish_rgb_RH(:,:,3)>=limpix);
%test_im(pL)=255;

% % test_im(find(test_im<limpix))=0;
% % test_im(find(test_im>=limpix))=255; 
% figure(4)
% imshow(test_im);
% hold on
% plot(xH,yH,'y--','LineWidth',0.1),
% hold off

%kopi av bilde innenfor r90-sirkelen:
test_im2=0*test_im;
test_im2(pr)=1;

test2_im2=0*test2_im;
test2_im2(pr)=1;

test3_im2=0*test3_im;
test3_im2(pr)=1;

test4_im2=0*test4_im;
test4_im2(pr)=1;


% figure(5)
% imshow(test_im2.*test_im)
% hold on
% plot(xH,yH,'y--','LineWidth',0.1),
% hold off
%
% close(2),close(3),close(4),close(5)

%Speilvend bildet - dvs flip kolonner. Dermed blir det som å se himmelen projisert ovenfra og ned,
%med riktig NESW orientering
I1= test_im2.*test_im; %et bilde med verdier 0 utenfor r90 sirkelen, og innenfor er blåverdier >limpix=255, resten 0
I2 = test2_im2.*test2_im;
I3 = test3_im2.*test3_im;
I4 = test4_im2.*test4_im;
clear test_im test_im2 test2_im test2_im2 test3_im test3_im2 test4_im test4_im2

I = I1 + I2 + I3 + I4;

fish_bw_RH_mirror= I; %flipdim(I,2);%speilvendt binærbilde
clear I
fish_rgb_RH_mirror=flipdim(fish_rgb_RH,2);%speilvendt fish_rgb_RH bilde

figure(1)
imshow(fish_rgb_RH_mirror) % Viser speilvendt normalt rgb bilde

figure(2)
imshow(fish_bw_RH_mirror)
% hold on
% plot(xH,yH,'y--','LineWidth',0.1)
% plot(im_cntrx,im_cntry,'yo','MarkerSize',10)
% plot([1+xH(1) N-xH(1)],[im_cntry im_cntry],'Color','y','LineWidth',1)
% plot([im_cntrx im_cntrx],[1+xH(1) M-xH(1)],'Color','y','LineWidth',1)

pW=find(t==-pi);
pN=find(t==-pi/2);
pE=find(t==0);
pS=find(t==pi/2);
dpix=round((xH(pW)-x(pW))/2);

plot(xH(pW),yH(pW),'yo','MarkerFaceColor','y','MarkerSize',4)
txt=text(x(pW),y(pW),'W','Color','y','FontSize',15);
d_cx=txt.Extent(3)/2;%pixellengde til bokstaven W
d_cy=txt.Extent(3)/2;%pixelhøyde til W
delete(txt)
text(x(pW)+dpix-d_cx,y(pW),'W','Color','y','FontSize',15);

plot(xH(pE),yH(pE),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pE)-dpix-d_cx,y(pE),'E','Color','y','FontSize',15)

plot(xH(pN),yH(pN),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pN)-d_cx,y(pN)+dpix,'N','Color','y','FontSize',15)

plot(xH(pS),yH(pS),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pS)-d_cx,y(pS)-dpix,'S','Color','y','FontSize',15)
hold off

%lag også et bilde med firkanten utenfor sirkelen markert hvitt
fish_bw_RH_grayspace=fish_bw_RH_mirror;
p_out=find(R>r90);
fish_bw_RH_grayspace(p_out)=155;

figure(3)
imshow(fish_bw_RH_grayspace)
hold on
plot(xH,yH,'y--','LineWidth',0.1)
plot(im_cntrx,im_cntry,'yo','MarkerSize',10)
plot([1+xH(1) N-xH(1)],[im_cntry im_cntry],'Color','y','LineWidth',1)
plot([im_cntrx im_cntrx],[1+xH(1) M-xH(1)],'Color','y','LineWidth',1)

pW=find(t==-pi);
pN=find(t==-pi/2);
pE=find(t==0);
pS=find(t==pi/2);
dpix=round((xH(pW)-x(pW))/2);

plot(xH(pW),yH(pW),'yo','MarkerFaceColor','y','MarkerSize',4)
txt=text(x(pW),y(pW),'W','Color','y','FontSize',15);
d_cx=txt.Extent(3)/2;%pixellengde til bokstaven W
d_cy=txt.Extent(3)/2;%pixelhøyde til W
delete(txt)
text(x(pW)+dpix-d_cx,y(pW),'W','Color','y','FontSize',15);

plot(xH(pE),yH(pE),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pE)-dpix-d_cx,y(pE),'E','Color','y','FontSize',15)

plot(xH(pN),yH(pN),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pN)-d_cx,y(pN)+dpix,'N','Color','y','FontSize',15)

plot(xH(pS),yH(pS),'yo','MarkerFaceColor','y','MarkerSize',4)
text(x(pS)-d_cx,y(pS)-dpix,'S','Color','y','FontSize',15)
hold off


%monter hvert bilde i en figur
%test medianfiltrering og differanse fra orgbilde

%K=medfilt2(fish_rgb_RH(:,:,3),[5 5]);
%K=medfilt2(fish_bw_RH_mirror,[5 5]);
%imshow(fish_rgb_RH,K,'montage')

%prøv å hent ut tabeller med rg og b verdier
%velg logisk indeksering:r(p_out)=uint8(100)
%sett ¨så sammen igjen en tre-dim bildefil
r_mirr=fish_rgb_RH_mirror(:,:,1);r_mirr(R>r90)=255;
g_mirr=fish_rgb_RH_mirror(:,:,2);g_mirr(R>r90)=255;
b_mirr=fish_rgb_RH_mirror(:,:,3);b_mirr(R>r90)=255;
fish_rgb_RH_grayspace(:,:,1)=r_mirr;
fish_rgb_RH_grayspace(:,:,2)=g_mirr;
fish_rgb_RH_grayspace(:,:,3)=b_mirr;
% figure(103)
% imshowpair(fish_rgb_RH_grayspace,medfilt2(fish_bw_RH_grayspace,[5 5]),'montage')

% %differanse mellom 2 biler:
% figure
% imshowpair(fish_bw_RH_mirror,K,'diff')

%Hvert bildepunkt svarer til vinkelen som objektet på himmelen danner i forhold til innfallsloddet
%Radius til hvert bildepunkt:
[X,Y]=meshgrid(1:M,1:N);
[Az_pix,R_pix] = cart2pol(X-im_cntrx,Y-im_cntry);%theta og azi vinkel til hvert bildepunkt. theta =0 til pi/2 og phi=0 +/-pi
%skyview=sum(sin*cos*dr*dphi*Skymatrise 0 og 1)/sum(sin*cos*dr*dphi*himmelhalvkulen (1))
Msky=fish_bw_RH_mirror;Msky(Msky>0)=1;
theta=R_pix*(pi/2)/r90;%theta er polarvinkel - (theta,azi) danner polarkoordinatene til hvert bildepunkt
th_x=theta.*cos(Az_pix);%skalar x og y koordinatene til bildepunktene
th_y=theta.*sin(Az_pix);

%beregne areal dtheta*dphi rundt hvert bildepunkt:=dx*dy
dth_x=abs(([zeros(size(th_x,1),1),th_x(:,1:end-1)]-[th_x(:,2:end),zeros(size(th_x,1),1)]))/2;
dth_y=abs([zeros(1,size(th_y,2));th_y(1:end-1,:)]-[th_y(1:end-1,:);zeros(1,size(th_y,2))])/2;

sv_mat_t=abs(sin(theta).*cos(theta).*dth_x.*dth_y.*double(Msky));
sv_mat_n=abs(sin(theta).*cos(theta).*dth_x.*dth_y.*ones(size(Msky)));
svf=sum(sum(sv_mat_t(pr)))/sum(sum(sv_mat_n(pr)));%fish_bw_RH_mirror-view faktor hvor innenfor 90-sirkelen

sa_mat_t=abs(sin(theta).*dth_x.*dth_y.*double(Msky));
sa_mat_n=abs(sin(theta).*dth_x.*dth_y.*ones(size(Msky)));
solid_ang_fraction=sum(sum(sa_mat_t(pr)))/sum(sum(sa_mat_n(pr)));%hvor stor romvinkelen til fri himmel er i forhold til en åpen horisont

figure(3)
hold on
title([sprintf('Skye-view factor %4.1f %%. ',svf*100),sprintf('Solid angle fraction %4.1f %%',solid_ang_fraction*100)],'fontSize',16)
% saveas(gcf, fullfile(outpath, "Skyview",join([filestring, "skyviewfactor.png"],"_"))) % lagrer skyviewfactor bildet
hold off
%test om det går an å få skilt himmelen ved å subtrahere blå-rød kanal

figure(102)
imshowpair(fish_rgb_RH_mirror,fish_bw_RH_grayspace,'montage')

figure(4)
mesh(R_pix)
axis square

figure(5)
mesh(Az_pix*180/pi)
axis square
xlabel('x')

%Beregn solbane for en gitt dato i Oslo og tegn inn:
site='Oslo';%for beregning på skråstilt flate
%site ='Trondheim';
%site ='Tromsø';
%site ='Düsseldorf';
%site ='Equator';
switch site
    case 'Oslo'
        lat=59.938;
        lon=10.719;
    case 'Trondheim'
        lat=63.42;
        lon=10.40;
    case 'Tromsø'
        lat=69+41/60;
        lon=18+58/60;
    case 'Düsseldorf'
        lat=51+13/60+18/3600;
        lon=6+46/60+34/3600;
    case 'Equator'
        lat=0;
        lon=0;
    otherwise
        return
end

DecH=(0:1/60:24-1/60)';%hvert 30 minutt i døgnet
Zen=zeros(length(DecH),365);
Azi=Zen;
EDIR365=Zen;
EDN365=Zen;
EUP365=Zen;

%her må det gjøres endringer...
actioncurve='EPP';
%ozone=400;
ozone=350;
wl=405;%ymax=7W/m^2/nm

for m=1:365 %dagnummer
    dnum=m;
    [aar,maaned,dag]=daynum2date_guv(2021,dnum);
    [Azi(:,m),Zen(:,m)]=solpos_bj(aar,maaned,dag,DecH+1,0,0,lat,lon);
    SZA365=Zen*180/pi;
    AZI365=Azi*180/pi;
    %     EDIR365(:,m)=interp1(UV.SZA,UV.EDIR_w(:,pO),SZA365(:,m),'linear')*UV.Dagvar(dnum);
    %     EDN365(:,m)=interp1(UV.SZA,UV.EDN_w(:,pO),SZA365(:,m),'linear')*UV.Dagvar(dnum);
    %     EUP365(:,m)=interp1(UV.SZA,UV.EUP_w(:,pO),SZA365(:,m),'linear')*UV.Dagvar(dnum);
    %
    %     %Gjør det samme for snø - brukes senere til å beregne faktor økning i
    %     %EPP som funksjon SZA, relativt til snøfritt
    %     EDIR365S(:,m)=interp1(B.SZA,B.EDIR_w(:,pO),SZA365(:,m),'linear')*B.Dagvar(dnum);%S for snow
    %     EDN365S(:,m)=interp1(B.SZA,B.EDN_w(:,pO),SZA365(:,m),'linear')*B.Dagvar(dnum);
    %     EUP365S(:,m)=interp1(B.SZA,B.EUP_w(:,pO),SZA365(:,m),'linear')*B.Dagvar(dnum);
    
end

%et punkt langs solbanen er gitt som (rcosphi,rsinphi), der r=sza*skaleringsfaktor
r_SZA=r90*(SZA365/90);%avstandsvektor fra sentrum ut til SZA til solen
pUH=find(SZA365>90);%sola under horisonten
r_SZA(pUH)=0;
xS=im_cntrx+r_SZA.*sin(AZI365*pi/180);
yS=im_cntry+r_SZA.*cos(AZI365*pi/180);


%Beregn solpos da bildet ble tatt- klokken var her stilt etter norsk sommertid. UTC-må trekke fra 2 timer (ikke tilfellet for 2021 bilder, 
% derfor er +2 blitt fjernet i: ,imdate(4)+1-0,imdate(5) etc. Ble satt
% tilbake da dette ikke ser ut til å stemme 
dayIm=dagnr(im_date(1),im_date(2),im_date(3));%dagen fisheyebildene ble tatt i barnehager og skole
[AzIm,ZenIm]=solpos_bj(im_date(1),im_date(2),im_date(3),im_date(4)+1-0,im_date(5),im_date(6),lat,lon);
xIm=im_cntrx+r90*(ZenIm/(pi/2))*sin(AzIm);
yIm=im_cntry+r90*(ZenIm/(pi/2))*cos(AzIm);
%
% %tester om radius til solposisjon er lineært økende med SZA
% xImT=im_cntrx+.*sin(AZI365*pi/180);
% yImT=im_cntry+r90*(SZA365/90).*cos(AZI365*pi/180);

%,im_date(4),im_date(5),im_date(6)
% %solbane i piksel koordinater
% xS=im_cntrx+r_SZA(:,day).*sin(AZI365(:,day)*pi/180);
% yS=im_cntry+r_SZA(:,day).*cos(AZI365(:,day)*pi/180);

%solskiven dekker en vinkel på 0.5 grader sett fra jorden
std_sun=round((r90/90)*0.5);
%std_sun=round(1);
hfilter=ones(std_sun,std_sun)/(std_sun*std_sun);
%sunfiltered=imgaussfilt(fish_bw_RH_grayspace,std_sun);
sunfiltered=imfilter(fish_bw_RH_grayspace,hfilter);
% figure
% imshowpair(fish_bw_RH_grayspace,sunfiltered,'montage')

[~,pT,~]=intersect(floor(DecH*1000),[6:2:16]*1000);%vil markere punkter langs solbanen for disse UTC tidene
% 
% figure(3)
% hold on
% plot(xS(:,dayIm),yS(:,dayIm),'y.','MarkerSize',5,'LineWidth',0.1)
% plot(xIm,yIm,'or','MarkerSize',18,'MarkerFaceColor','r')
% hold off

% figure(30)
% imshow(fish_bw_RH_grayspace)
% hold on
% plot(xS(:,dayIm),yS(:,dayIm),'ro','MarkerSize',8,'MarkerFaceColor','r')
% plot(xIm,yIm,'pr','MarkerSize',30,'MarkerFaceColor','r')
% plot(xS(pT,dayIm),yS(pT,dayIm),'ko','MarkerSize',5,'MarkerFaceColor','y')
% hold off

% figure(31)%solbane i fisheye farge
% imshow(fish_rgb_RH_mirror)
% hold on
% plot(xS(:,dayIm),yS(:,dayIm),'r.','MarkerSize',5,'LineWidth',0.1)
% plot(xIm,yIm,'pr','MarkerSize',15)
% plot(xS(pT,dayIm),yS(pT,dayIm),'ko','MarkerSize',5,'MarkerFaceColor','y')
% % saveas(gcf, fullfile(outpath, "Solbaner", join([filestring, "solbane.png"],"_"))) % Lagrer solbane bildet
% hold off

%tegn inn solbane i panormabilde figure(104)
p=find(AZI365(:,dayIm)>=0);
%u=[AZI365(p,dayIm)';SZA365(p,dayIm)']';
u=[360-AZI365(p,dayIm)';SZA365(p,dayIm)']';
%v=sortrows(u,1);
%v=sortrows(u,-1);

p=find(AZI365(:,dayIm)<0);
%w=flipud([360+AZI365(p,dayIm)';SZA365(p,dayIm)']');
w=([AZI365(p,dayIm)';SZA365(p,dayIm)']');
u_comb=[abs(w);u];
xtemp=floor(u_comb(:,1)*3600/360);%*size(imP,2)/size(u_comb,1)
ytemp=round(900-(90-u_comb(:,2))*900/90);

%%Panoramabilder:
% % figure(44)%prøver å strekke bildet ut til et panoramabilde
% % %plot3(Az_pix(:)+pi,R_pix(:)*(pi/2)/r90,double(fish_bw_RH_mirror(:)))
% % C=fish_rgb_RH;
% % test1=double(C(:,:,1))/255;
% % test2=double(C(:,:,2))/255;
% % test3=double(C(:,:,3))/255;
% %
% % test=[test1(:)';test2(:)';test3(:)']';
% % scatter(Az_pix(:)+pi,R_pix(:)*(pi/2)/r90,20,test,'filled')
% % %interpolere til felles grid?
% % scatter(1:50,cos((1:50)*pi/50),30,colors,'filled')
% %
% % imshow(Az_pix*pi/2,R_pix)
%
%

panorama=0;
%panorama=1;
if panorama
    %Metode som fungerer til å lage panorama:
    %tettere beskåret utsnitt, kvadratisk som omkranser sirkulært område:
    rect_ny=[im_cntrx-r90 im_cntry-r90 2*r90 2*r90];
    fish_bw_R90_grayspace = imcrop(fish_bw_RH_grayspace,rect_ny);
    
    %im = double(fish_bw_RH_grayspace)/255.0;
    %im = double(fish_bw_R90_grayspace)/255.0;
    im = double(fish_bw_R90_grayspace)/255.0;
    im=flipdim(im,2);%flipper tilbake til original - blir et rettvendt bilde, som må til for å vise panorama riktig
    figure(1); imshow(im);
    tic
    imP = ImToPolar(im, 0, 1, 900, 3600);
    toc
    figure(104); imshow(imP);
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    %text(size(imP,2)/2, 100, 'Flattended fisheye', 'FontSize', 20);
    set(gca,'Position',[0.02 0.02 0.96 0.96])
    hFig=figure(104)
    hgexport(hFig,'-clipboard')
    
    %Det samme for croppet fish_rgb_RH-bilde:
    %%fish_rgb_R90_grayspace = imcrop(fish_rgb_RH_grayspace,rect_ny);
    fish_rgb_R90_grayspace = imcrop(fish_rgb_RH_grayspace,rect_ny);
    %im_max=255;
    im_max=1;
    im_r=double(fish_rgb_R90_grayspace(:,:,1))/im_max; im_r=flipdim(im_r,2);%flipper til rettvendt bilde
    im_g=double(fish_rgb_R90_grayspace(:,:,2))/im_max; im_g=flipdim(im_g,2);%flipper til rettvendt bilde
    im_b=double(fish_rgb_R90_grayspace(:,:,3))/im_max; im_b=flipdim(im_b,2);%flipper til rettvendt bilde
    
    tic
    imP_col(:,:,1) = ImToPolar(im_r, 0, 1, 900, 3600);
    imP_col(:,:,2) = ImToPolar(im_g, 0, 1, 900, 3600);
    imP_col(:,:,3) = ImToPolar(im_b, 0, 1, 900, 3600);
    toc
    
    figure(105); imshow(uint8(imP_col));
    hold on
    rectangle('Position', [1, 1, size(imP_col,2)-1, size(imP_col,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    hFig=figure(105)
    hgexport(hFig,'-clipboard')
    
    
    figure(106)
    imshowpair(uint8(imP_col),imP,'montage')
    
    hfig=figure(107)
    subplot(2,1,1)
    imshow(uint8(imP_col));
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    
    subplot(2,1,2)
    imshow(imP);
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    
    %function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    %ha = tight_subplot(2,1,[.01 .03],[.1 .01],[.01 .01])
    ha = tight_subplot(2,1,[.001 .003],[.01 .01],[.01 .01])
    axes(ha(1));imshow(uint8(imP_col));
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    
    axes(ha(2));imshow(imP);
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 2)
    hold off
    
    hFig=figure(107);
    hgexport(hFig,'-clipboard')
    savefig(sprintf('%s\\panorama_%s.fig',outpath,filestring));
    % Enlarge figure to full screen.
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    export_fig(sprintf('%s\\panorama_%s',outpath,filestring),'-png');
    
    
    
    
    hFig=figure(3)
    hgexport(hFig,'-clipboard')
    savefig(sprintf('%s\\circular_%s.fig',outpath,filestring));
    % Enlarge figure to full screen.
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    export_fig(sprintf('%s\\circular_%s',outpath,filestring),'-png');
    
    % %Fungerer ikke: alt dette er forsøk på å brette ut fisheye til panorama, uten å lykkes
    % imC = Polar2Im(fish_rgb_RH_grayspace,3600,'linear');
    % imP = FISHCOLOR2(fish_bw_RH_grayspace);
    %
    %
    % r90
    % %Az_pix(pr)
    % A=uint8(0*R_pix);
    % A(pr)=uint8(255);
    % figure
    % imshow(A)
    % w=[];
    % e=[];
    % for i=1:M
    %    pos=find(A(:,i)==255);
    %    if numel(pos)>0
    %        if numel(pos)>1
    %           w=[w;Az_pix(pos(1),i)];
    %           e=[e;Az_pix(pos(end),i)];
    %        end
    %    end
    % end
    %
    % az_grid=[w;sortrows(e)];
    % % figure
    % % plot(az_grid)
    %
    % %theta spenn fra 0 til 90 grader:
    % th_range=R_pix((im_cntrx:im_cntrx+r90),im_cntry)*90/r90;
    % [pan_phi,pan_th]=meshgrid(az_grid,th_range);
    % % pan_th(:,1)
    % % pan_phi(1,:)
    %
    % % %lag panoramamatrise phixtheta og interpoler gråverdier
    % % pan=0*pan_phi;
    % % for i=1:size(pan,2)
    % %     pan(:,i)=interp1(Az_pix(i,:),double(fish_bw_RH_mirror(i,:)),az_grid,'linear');
    % % end
    %
    % xi=Az_pix'; xi=xi(:);
    % yi=R_pix';yi=yi(:);
    % zi_r=double(fish_rgb_RH(:,:,1))';zi_g=double(fish_rgb_RH(:,:,2))';zi_b=double(fish_rgb_RH(:,:,3))';
    % zi_r=zi_r(:);zi_g=zi_g(:);zi_b=zi_b(:);
    % tic
    % F_r = scatteredInterpolant(xi,yi,zi_r);
    % F_g = scatteredInterpolant(xi,yi,zi_g);
    % F_b = scatteredInterpolant(xi,yi,zi_b);
    % toc
    %
    %
    % %double(fish_rgb_RH(:));
    % %[xq,yq]=meshgrid(-pi:(2*pi/360):pi,(0:1:90));
    % [xq,yq]=meshgrid(-pi:(2*pi/360)/10:pi,(90:-1/10:0));
    % xq_ny=xq(:); yq_ny=yq(:);
    % [ii,jj]=size(xq);
    % tic
    % V_r=F_r(xq_ny,yq_ny);
    % V_g=F_r(xq_ny,yq_ny);
    % V_b=F_r(xq_ny,yq_ny);
    % toc
    %
    % pan_mat(:,:,1)=reshape(V_r,[ii,jj]);pan_mat(:,:,2)=reshape(V_g,[ii,jj]);pan_mat(:,:,3)=reshape(V_b,[ii,jj]);
    % pan_mat=uint8(pan_mat);
    %
    % figure
    % imshow(pan_mat)
    %
    %
    % tic
    % vq_ny_r=griddata(xi,yi,zi_r,xq_ny,yq_ny);vq_ny_g=griddata(xi,yi,zi_g,xq_ny,yq_ny);vq_ny_b=griddata(xi,yi,zi_b,xq_ny,yq_ny);
    % toc
    %
    % col_mat(:,:,1)=reshape(vq_ny_r,[ii,jj]);col_mat(:,:,2)=reshape(vq_ny_g,[ii,jj]);col_mat(:,:,3)=reshape(vq_ny_b,[ii,jj]);
    % col_mat=uint8(col_mat);
    % %col_mat(20,:)
    % figure
    % %pcolor(xq,yq,col_mat),colorbar
    % imshow(col_mat)
    %
    % %test_r=griddata(Az_pix',R_pix',double(fish_rgb_RH(:,:,1))',xq_ny,yq_ny);
    % %test_r=griddata(Az_pix',R_pix',double(fish_rgb_RH(:,:,1))',xq,yq);
    % test_r=griddata(Az_pix',R_pix',double(fish_rgb_RH(:,:,1))',yq,xq);
    %
    % figure
    % imshow(reshape(test_r,ii,jj))
    %
    % figure,plot(xi,'.'),grid
    % figure,plot(yi,'.'),grid
    % figure,plot(zi_r,'.'),grid
    %
    % figure,plot(xq_ny,'.'),grid
    % figure,plot(yq_ny,'.'),grid
    % figure,plot(vq_ny,'.'),grid
    %
    %
    % figure
    % pcolor(reshape(xq_ny,91,361),reshape(yq_ny,91,361),reshape(vq_ny,91,361))
    %
    % %n_mat=
    %
    %
    % pan_grid=meshgrid(linspace(0),   linspace(-pi,(2*pi/360)/10,pi));%hvor mange piksler rundt r90-sirkelen?
    % % R_grid=linspace(min(min(R_pix)),max(max(R_pix)),numel(R_pix));
    % %
    % % r_chan=griddata(Az_pix,R_pix,double(fish_rgb_RH(:,:,1)),az_grid,R_grid,'nearest');
    % % g_chan=griddata(Az_pix,R_pix,double(fish_rgb_RH(:,:,2)),az_grid,R_grid,'nearest');
    % % b_chan=griddata(Az_pix,R_pix,double(fish_rgb_RH(:,:,3)),az_grid,R_grid,'nearest');
    % % res(:,:,1)=reshape(r_chan,[1343,1343]);
    % % res(:,:,2)=reshape(g_chan,[1343,1343]);
    % % res(:,:,3)=reshape(b_chan,[1343,1343]);
    % % figure
    % % imshow(uint8(res))
    % % %mesh(xpan,ypan,fish_bw_RH_mirror,'.')
    % % plot3(xpan,ypan,fish_bw_RH_mirror,'.')
    
    
    %%*****
    
    figure(32)%tegner inn solbanen i panoramabilde
    imshow(imP)
    hold on
    plot(xtemp(1:60:end),ytemp(1:60:end),'py-','MarkerSize',15,'MarkerFaceColor','y','LineWidth',4)
    hold off
    
    
    figure(107)%tegn inn solbaner
    axes(ha(1));
    %subplot(2,1,1)
    hold on
    plot(xtemp(1:60:end),ytemp(1:60:end),'py-','MarkerSize',15,'MarkerFaceColor','y','LineWidth',4)
    if AzIm>0
        xtemp_sun=round((360-AzIm*180/pi)*3600/360);
    else
        xtemp_sun=round(abs(AzIm*180/pi)*3600/360);
    end
    plot(xtemp_sun,round(900-(90-ZenIm*180/pi)*900/90),'pr','MarkerSize',15,'MarkerFaceColor','r')
    hold off
    
    axes(ha(2));
    %subplot(2,1,2)
    hold on
    plot(xtemp(1:60:end),ytemp(1:60:end),'py-','MarkerSize',15,'MarkerFaceColor','y','LineWidth',4)
    plot(xtemp_sun,round(900-(90-ZenIm*180/pi)*900/90),'pr','MarkerSize',15,'MarkerFaceColor','r')
    hold off
    
    hFig=figure(107);
    hgexport(hFig,'-clipboard')
    
end %if panorama


%Beregn bildepikselkoordinater langs solbanen og plott hvitt/svart langs banen
%svart hvitt bilde: fish_bw_RH_grayspace
[X,Y] = meshgrid(1:365,DecH);
pix_bane=0*meshgrid(1:365,DecH);
pix_bane_avg=pix_bane;
for dayIm = (1:365)
    pOH=find(SZA365(:,dayIm)<=90);
    % pix_bane=[];
    % pix_bane_avg=[];
    %pix_bane=0*DecH;
    %pix_bane_rev=0*DecH;
    for i=1:numel(pOH)
        %     pix_bane=[pix_bane;fish_bw_RH_grayspace(round(xS(i,dayIm)),round(yS(i,dayIm)),1)];
        %     pix_bane_avg=[pix_bane_avg;sunfiltered(round(xS(i,dayIm)),round(yS(i,dayIm)),1)];
        pix_bane(pOH(i),dayIm)=fish_bw_RH_grayspace(round(yS(pOH(i),dayIm)),round(xS(pOH(i),dayIm)));
        pix_bane_avg(pOH(i), dayIm)=sunfiltered(round(yS(pOH(i),dayIm)),round(xS(pOH(i),dayIm)));
        %pix_bane_avg(i)=[pix_bane_avg;sunfiltered(round(xS(i,dayIm)),round(yS(i,dayIm)),1)];
    end
end
pix_bane=double(pix_bane)/255;
pix_bane_avg=double(pix_bane_avg)/255;



%figure(51) %pikselverdier som funksjon av tid på dagen
%surf(X,Y,pix_bane) % Overflate plot av pikselverdier som funksjon av tid på dagen, gjennom alle dager i året
% plot(DecH,pix_bane,'ko-')
% hold on
% plot(DecH,pix_bane_avg,'ro-')
% hold off
%axis([floor(DecH(pOH(1))) ceil(DecH(pOH(end))) 0 1.05])

%figure(52) %pikselverdier som funksjon av tid på dagen
%surf(X,Y,pix_bane_avg) % Overflate plot av avg_pix verdier som funk av tpd, hele året.
% plot(DecH,pix_bane_avg,'ro-')
% axis([floor(DecH(pOH(1))) ceil(DecH(pOH(end))) 0 1.05])

%Finn pikselplassering av DecH i panormabilde figure(104)
p=find(AZI365(:,dayIm)>=0);
clear u
u=[360-AZI365(p,dayIm)';pix_bane_avg(p)']';

p=find(AZI365(:,dayIm)<0);
ww=([AZI365(p,dayIm)';pix_bane_avg(p)']');
DechH_comb=[abs(ww);u];
xH=floor(DechH_comb(:,1)*3600/360);%*size(imP,2)/size(u_comb,1)
yH=round(50+(1-DechH_comb(:,2))*300);

if panorama
    
    hfig=figure(108)
    
    subplot(2,1,1)
    %axes(ha(1));
    imshow(uint8(imP_col));
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 1)
    hold off
    
    subplot(2,1,2)
    %axes(ha(2));
    imshow(imP);
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 1)
    hold off
    
    %function ha = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
    %ha = tight_subplot(2,1,[.01 .03],[.1 .01],[.01 .01])
    %hb = tight_subplot(2,1,[.005 .005],[.1 .01],[.05 .05])
    hb = tight_subplot(2,1,[.001 .003],[.01 .01],[.01 .01])
    axes(hb(1));imshow(uint8(imP_col));
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 1)
    hold off
    
    axes(hb(2));imshow(imP);
    hold on
    rectangle('Position', [1, 1, size(imP,2)-1, size(imP,1)-1],...
        'EdgeColor','k', 'LineWidth', 1)
    hold off
    
    %figure(108)
    axes(hb(1));
    %subplot(2,1,1)
    hold on
    plot(xtemp(1:60:end),ytemp(1:60:end),'py--','MarkerSize',15,'MarkerFaceColor','y','LineWidth',1.5)
    plot(xtemp,ytemp,'y--','LineWidth',2)
    if AzIm>0
        xtemp_sun=round((360-AzIm*180/pi)*3600/360);
    else
        xtemp_sun=round(abs(AzIm*180/pi)*3600/360);
    end
    plot(xtemp_sun,round(900-(90-ZenIm*180/pi)*900/90),'pr','MarkerSize',15,'MarkerFaceColor','r')
    hold off
    
    %subplot(2,1,2)
    axes(hb(2));
    hold on
    plot(xtemp(1:60:end),ytemp(1:60:end),'py--','MarkerSize',15,'MarkerFaceColor','y','LineWidth',1.5)
    plot(xtemp,ytemp,'y--','LineWidth',2)
    plot(xtemp_sun,round(900-(90-ZenIm*180/pi)*900/90),'pr','MarkerSize',15,'MarkerFaceColor','r')
    plot(xH,yH,'y-','LineWidth',1.5)
    hold off
    
    
    % Enlarge figure to full screen.
    %set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    hFig=figure(108);
    hgexport(hFig,'-clipboard')
    savefig(sprintf('%s\\panorama_with_sun_%s.fig',outpath,filestring));
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    export_fig(sprintf('%s\\panorama_with_sun_%s',outpath,filestring),'-png');
    
    system('echo off | clip')
end %if panorama

%Nå kan direkte og diffusbidragt på en flate beregnes vs DecH, vektet med skyviewfkator og
%shadowfaktor:

%For klar himmel er transmittans = 1:
T_clearsky(:,1)=(290:0.5:800)';
T_clearsky(:,2)=1;

%clearsky=1;
clearsky=0;
if clearsky
    clod=0;
else
    clod=10;
end
FLATE_cone = porfyri_returner_horisontalkomponenter('clear_sky',clod,350,'snow free',lat,lon,'cone',60,0,'CIE',T_clearsky);
E_cone_free=FLATE_cone.DIR_surface_w+FLATE_cone.DDN_surface_w+FLATE_cone.RUP_surface_w;%åpen, fri himmel
E_cone_free(isnan(E_cone_free))=0;

% figure(110)
% plot(DecH,E_cone_free(:,dayIm),'r-','LineWidth',3)
% hold on

FLATE_horizontal = porfyri_returner_horisontalkomponenter('clear_sky',clod,350,'snow free',lat,lon,'horizontal_element',60,0,'CIE',T_clearsky); 
E_horizontal_free=FLATE_horizontal.DIR_surface_w+FLATE_horizontal.DDN_surface_w+FLATE_horizontal.RUP_surface_w;%åpen, fri himmel
E_horizontal_free(isnan(E_horizontal_free))=0;

% plot(DecH,E_horizontal_free(:,dayIm),'k-','LineWidth',3)

%Beregn irradians for solbane når skyline blokkerer direktekomponenten og horisonten har begrenset skyview:
%pix_bane_avg %piksler med verdi 1 betyr åpen sol, 0 er blokkert sol

% Dir_cone_sky=FLATE_cone.DIR_surface_w(:,dayIm).*pix_bane_avg;
% Dn_cone_sky=FLATE_cone.DDN_surface_w(:,dayIm)*svf;
% R_cone_opp=((FLATE_cone.EDIR365(:,dayIm).*pix_bane_avg + FLATE_cone.EDN365(:,dayIm)*svf)*FLATE_cone.Albedo).*FLATE_cone.ShapeUp(:,dayIm);
% R_cone_opp(isnan(R_cone_opp))=0;
% E_cone_sv=Dir_cone_sky+Dn_cone_sky+R_cone_opp;E_cone_sv(isnan(E_cone_sv))=0;


% FLATE_horizontal.DIR_surface_w(:,dayIm).*pix_bane_avg+FLATE_horizontal.DDN_surface_w(:,dayIm)+FLATE_horizontal.RUP_surface_w(:,dayIm);
% Dir_horizontal_sky=FLATE_horizontal.DIR_surface_w(:,dayIm).*pix_bane_avg;
% Dn_horizontal_sky=FLATE_horizontal.DDN_surface_w(:,dayIm)*svf;
% R_horizontal_opp=((FLATE_horizontal.EDIR365(:,dayIm).*pix_bane_avg + FLATE_horizontal.EDN365(:,dayIm)*svf)*FLATE_horizontal.Albedo).*FLATE_horizontal.ShapeUp(:,dayIm);
% R_horizontal_opp(isnan(R_horizontal_opp))=0;
% E_horizontal_sv=Dir_horizontal_sky+Dn_horizontal_sky+R_horizontal_opp;E_horizontal_sv(isnan(E_horizontal_sv))=0;

Dir_cone_sky=FLATE_cone.DIR_surface_w.*pix_bane_avg;
Dn_cone_sky=FLATE_cone.DDN_surface_w*svf;
R_cone_opp=((FLATE_cone.EDIR365.*pix_bane_avg + FLATE_cone.EDN365*svf)*FLATE_cone.Albedo).*FLATE_cone.ShapeUp;
R_cone_opp(isnan(R_cone_opp))=0;
E_cone_sv=Dir_cone_sky+Dn_cone_sky+R_cone_opp;E_cone_sv(isnan(E_cone_sv))=0;

FLATE_horizontal.DIR_surface_w.*pix_bane_avg+FLATE_horizontal.DDN_surface_w+FLATE_horizontal.RUP_surface_w;
Dir_horizontal_sky=FLATE_horizontal.DIR_surface_w.*pix_bane_avg;
Dn_horizontal_sky=FLATE_horizontal.DDN_surface_w*svf;
R_horizontal_opp=((FLATE_horizontal.EDIR365.*pix_bane_avg + FLATE_horizontal.EDN365*svf)*FLATE_horizontal.Albedo).*FLATE_horizontal.ShapeUp;
R_horizontal_opp(isnan(R_horizontal_opp))=0;
E_horizontal_sv=Dir_horizontal_sky+Dn_horizontal_sky+R_horizontal_opp;E_horizontal_sv(isnan(E_horizontal_sv))=0;

FLATE_vertical = porfyri_returner_horisontalkomponenter('clear_sky',clod,350,'snow free',lat,lon,'vertical_element',60,0,'CIE',T_clearsky);
E_vertical_free=FLATE_vertical.DIR_surface_w+FLATE_vertical.DDN_surface_w+FLATE_vertical.RUP_surface_w;%åpen, fri himmel
E_vertical_free(isnan(E_vertical_free))=0;
Dir_vertical_sky=FLATE_vertical.DIR_surface_w.*pix_bane_avg;
Dn_vertical_sky=FLATE_vertical.DDN_surface_w*svf;
R_vertical_opp=((FLATE_vertical.EDIR365.*pix_bane_avg + FLATE_vertical.EDN365*svf)*FLATE_vertical.Albedo).*FLATE_vertical.ShapeUp;
R_vertical_opp(isnan(R_vertical_opp))=0;
E_vertical_sv=Dir_vertical_sky+Dn_vertical_sky+R_vertical_opp;E_vertical_sv(isnan(E_vertical_sv))=0;




dose_fraksjon_cone=trapz(DecH,E_cone_sv)/trapz(DecH,E_cone_free(:,dayIm));
dose_fraksjon_horizontal=trapz(DecH,E_horizontal_sv)/trapz(DecH,E_horizontal_free(:,dayIm));

% figure(198)
% plot(DecH,E_cone_sv,'r-','LineWidth',1.5)
% figure(199)
% plot(DecH,E_horizontal_sv,'k-','LineWidth',1.5)
% 
% hold off
% legend('Cone, open horizon','Horizontal open horizon','Cone skyview','Horizontal skyview')
% xlabel('Time')
% ylabel('EPP-irradiance,W/m^2')
% title(sprintf('Irradiance on a cone and horizontal plane\nDosefraction: Cone: %2.0f %%, Horizontal: %2.0f %%',dose_fraksjon_cone*100,dose_fraksjon_horizontal*100))

%{
The following two plots will give us a surfaceplot of the UVI on both a
coned surface, and a horizontal for each day, 24 hours, throughout the
whole year.
%}

figure(201);
surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),E_cone_sv(1:10:end,1:5:end));hold on;
title({"UVI on a coned surface angled at 60 degrees", "(the bridge of the nose of a child running in circles)"})
xlabel("Days in a year [days]")
ylabel("Hours in a day [t]")

figure(202);
surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),E_horizontal_sv(1:10:end,1:5:end));hold on;
% contourf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),E_horizontal_sv(1:10:end,1:5:end));hold on;
title("UVI on a horizontal surface")
xlabel("Days in a year [days]")
ylabel("Hours in a day [t]")


IMG.fisheye_bildenavn_jpg = join([filestring, ".JPG"],""); % Filename of the photograph in question
%IMG.Bilde_tekst = ; % Some information regarding the photograph
IMG.Opptaksdatotid = im_date; % Date and time of photo
IMG.Blacknwhite = fish_bw_RH_grayspace; % Black n white of the mirrored fisheye photograph
IMG.sunfiltered = sunfiltered; % Filtered version of the BW photo
IMG.fish_rgb_RH_mirror = fish_rgb_RH_mirror; % The mirror version of the fisheye photograph
IMG.M_matrise = pix_bane_avg; % The matrix showing suns route through the sky throughout the day, for each day in a year
IMG.svf = svf; % Skyview factor for the particular photograph
IMG.UVI_cone = E_cone_sv * 40; % UVI for a coned surface, given the suns route in M_natrise
IMG.UVI_horizontal = E_horizontal_sv * 40; % UVI for a horizontal surface, given the suns route in M_natrise
IMG.UVI_vertical = E_vertical_sv * 40;
IMG.UVI_cone_free = E_cone_free * 40; % UVI for a free horizon on a coned surface
IMG.UVI_horizontal_free = E_horizontal_free * 40; % UVI for a free horizon on a horizontal surface
IMG.UVI_vertical_free = E_vertical_free * 40;
IMG.RGB_filter_values = {[upper_quadrant],[left_quadrant],[right_quadrant],[lower_quadrant]}; % Parameters to BW-filter the image and creating the bw-matrix

% save(join(filestring+"_data.mat"), "IMG");

hFig=figure(110);
hgexport(hFig,'-clipboard')

%
% %Herfra er det bare tull
% rgb_flip=flipdim(fish_rgb_RH,2);
%
% pix_val_r=[];
% pix_val_g=[];
% pix_val_b=[];
% for i=1:numel(DecH)
%     pix_val_r=[pix_val_r;rgb_flip(round(xS(i,dayIm)),round(yS(i,dayIm)),1)];
%     pix_val_g=[pix_val_g;rgb_flip(round(xS(i,dayIm)),round(yS(i,dayIm)),2)];
%     pix_val_b=[pix_val_b;rgb_flip(round(xS(i,dayIm)),round(yS(i,dayIm)),3)];
% end
%
%
% figure(10)
% plot(DecH,pix_val_r,'r-')
% hold on
% plot(DecH,pix_val_g,'g-')
% plot(DecH,pix_val_b,'b-')
% hold off
%
%
% figure(11)
% plot(AZI365(:,dayIm),pix_val_r,'r-')
% hold on
% plot(AZI365(:,dayIm),pix_val_g,'g-')
% plot(AZI365(:,dayIm),pix_val_b,'b-')
% hold off
%
% %
% %
% % % metode for å finne radius og azivinkler til bildepunkter
% % %M=fspecial('gaussian',256,32); % generate fake image. Bredde 256 piksler, sigma 32
% % M=fish_bw_RH_mirror; % generate fake image. Bredde 256 piksler, sigma 32
% % X0=size(M,1)/2; Y0=size(M,2)/2;%pikselverdier i origo, midt i bildet
% % % [Y X z]=find(M);%rad og kolonneindekser for alle bildepunkter >0
% % [Y X z]=find(fish_bw_RH_mirror);%rad og kolonneindekser for alle bildepunkter >0
% % X=X-X0; Y=Y-Y0;
% % theta = atan2(Y,X);
% % rho = sqrt(X.^2+Y.^2);
% %
% % % Determine the minimum and the maximum x and y values:
% % rmin = min(rho); tmin = min(theta);
% % rmax = max(rho); tmax = max(theta);
% %
% % % Define the resolution of the grid:
% % rres=611; % # of grid points for R coordinate. (change to needed binning)
% % tres=611; % # of grid points for theta coordinate (change to needed binning)
% %
% % F = TriScatteredInterp(rho,theta,double(z),'natural');
% %
% % %Evaluate the interpolant at the locations (rhoi, thetai).
% % %The corresponding value at these locations is Zinterp:
% %
% % [rhoi,thetai] = meshgrid(linspace(rmin,rmax,rres),linspace(tmin,tmax,tres));
% % Zinterp = F(rhoi,thetai);
% %
% % figure
% % subplot(1,2,1); imagesc(M) ; axis square
% % subplot(1,2,2); imagesc(Zinterp) ; axis square
% %
% %
% % %nytt eksempel fra internet
% % [rows, columns, numberOfColorBands] = size(fish_rgb_RH);
% % if numberOfColorBands > 1
% %     % It's not really gray scale like we expected - it's color.
% %     % Convert it to gray scale by taking only the green channel.
% %     grayImage = fish_rgb_RH(:, :, 2); % Take green channel.
% % end
% %
% % % Specify the distortion correction factor.
% % % Positive to correct pincushion distortion, and bring coreners inwards.
% % % Negative to correct barrel (fisheye) distortion, and push corners outward.
% % figure
% % df = -0.5;
% % xCenter = columns / 2;  % Assume optical axis in center of image.
% % yCenter = rows / 2;
% % samplingRate = 50;
% % [xBad, yBad] = meshgrid([1 : samplingRate : columns], [1: samplingRate   :rows]);
% % hold on;
% % plot(xBad, yBad, 'r+', 'MarkerSize', 5, 'LineWidth', 2);
% % for r = 1 : size(xBad, 1)
% %     for c = 1 : size(xBad, 2)
% %         x = xBad(r, c);
% %         y = yBad(r, c);
% %         % Get the actual distance from the optic axis.
% %         rBad = sqrt((x - xCenter)^2 + (y - yCenter)^2);
% %         % Get the delta R that is moved.
% %         deltaRadius = df * rBad / (1 + df);
% %         rTrue = rBad - deltaRadius;
% %         angle = atan2((y - yCenter), (x - xCenter));
% %         xTrue(r, c) = xCenter + rTrue * cos(angle);
% %         yTrue(r, c) = yCenter + rTrue * sin(angle);
% %         D = interp2(double(grayImage), x-xTrue, y-yTrue);
% %         % display the result
% %         imshow(D, []);
% %
% %     end
% % end
% % figure
% % plot(xTrue, yTrue, 'bo', 'MarkerSize', 5, 'LineWidth', 2);
% % legend({'+ = Distorted'; 'o = Corrected'});
%
%
%

%Beregner irradianser for Helset barnehage under solseil, hvor transmittans i teltduk er beregnet,
%og viewfaktorer er delt opp for fri sky, solseil sky og wall view:
if strcmp(imfile,'J:\UV, sol og helseeffekter\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2336.jpg') %limpix=130;%limpix=110;%under solseil i Helset barnehag
    %ved å sette pix_limit for blåkanal til 130 finnes skyview himmel. Ved å sette pix_limit til 35
    %blir alt hvitt unntatt trær og bygninger. Kan derfor separere skyview, tentview og wallview
    svf=0.339;%himmel
    vfsunsails=0.509;%solseilbidrag
    vfwalls=0.152;%Resten er view til vegger og trær
    
    %fra morgen til kveld krysser solen seilet to ganger:
    % fra indeks 462 til 564, og indeks 631 til 742
    
    psail=[462:564,631:742];
    pix_bane_sky=pix_bane_avg; pix_bane_sky(psail)=0;
    pix_bane_sunsail=0*pix_bane_avg;
    pix_bane_sunsail(psail)=1; %Når kombinerer pix_bane_sky og pix_bane_sunsail kan bidrag fra sky og transmitter i duken adderes
    
    actionspec='EPP';
    surface_type='cone';
    %surface_type='semi_cone_S';
    %surface_type='horizontal_element';
    
    beta=60; %for vertikal flate, horisontal flate settes beta automatisk til 0 og 90 grader
    
    ground_albedo='snow free';
    %ground_albedo='beach';
    %ground_albedo='fresh snow';
    
    %skytype='clear_sky';
    skytype='overcast';
    switch skytype
        case 'clear_sky'
            clod=0;
        case 'overcast'
            %clod=1;
            clod=10; %helt diffusert her, også for UV
            %clod=60
        otherwise
            return
    end
    
    
    %Transmittance of sunsail above red roof, compared with cloudy sky above and below: [RGB]=[108 87 54]./[164 169 166]= 0.69, 0.51, 0.33
    %senterbølgelengder til RGB:
    %Z:\Fagdokumenter\Optisk\UV_artikler\Master_og_doktor_avhandlinger\Thesis_skyview_Univ_Manchester_2015.pdf
    %Roberto Carrasco-Hernandez, 2015: s. 191, effective wavelengths of RGB channels: 460, 530, 600 nm
    %I publikasjon Utrillas 2010, om beach umbrella fant de en UVI transmittans på 34 %. Velger
    %30 prosent for bølgelengder 290-400 nm
    
    x=(290:0.5:800)';
    %yT=interp1([290 400 460 530 600 800],[0.3 0.3 0.33 0.51 0.69 0.75],x,'linear','extrap');%transmittans i duken
    yT=interp1([305 400 460 530 600 800],[0.04 0.2 0.33 0.51 0.69 0.75],x,'linear','extrap');%transmittans i duken
    
    figure(10)
    plot([305 460 530 600],[0.04 0.33 0.51 0.69],'ko','MarkerSize',10)
    hold on
    plot(x,yT,'k--','LineWidth',2)
    hold off
    legend('RGB channels','extrapolated')
    ylabel('Transmittance')
    xlabel('Wavelength, nm')
    %title('Transmittance of sun sail')
    axis([290 800 0 1])
    set(gca,'Fontsize',[12]);
    hFig=figure(10);
    hgexport(hFig,'-clipboard')
    
    Tsail(:,1)=x;
    Tsail(:,2)=yT;
    %Beregn EPP irrdianskomponenter uten solseil og med solseilduk, og for en antatt reflektans fra
    %trær og vegger på 10 prosent (blåkomponent):
    
    %vekt deretter bidrag fra hver, avhengig av hvor solbanen går
    
    %fri horisont og clod=0. Husk at shapefaktorer allerede er innberegnet i .DIR_surface_w,
    %.DDN_surface_w og .RUP_surface_w, mens .EDIR365 osv er horisontalkomponenten
    FLATE_cone_free_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,T_clearsky);
    E_cone_free_0=FLATE_cone_free_0.DIR_surface_w+...
        FLATE_cone_free_0.DDN_surface_w+...
        FLATE_cone_free_0.RUP_surface_w;%åpen, fri himmel
    E_cone_free_0(isnan(E_cone_free_0))=0;
    
    %fri horisont, for gitt clod
    FLATE_cone_free = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,T_clearsky);
    E_cone_free=FLATE_cone_free.DIR_surface_w +...
        FLATE_cone_free.DDN_surface_w+...
        FLATE_cone_free.RUP_surface_w;%åpen, fri himmel
    E_cone_free(isnan(E_cone_free))=0;
    
    %ingen fri horisont, clod=0, Vektet etter transmittans,
    FLATE_cone_sunsail_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Tsail);
    E_cone_sunsail_0=FLATE_cone_sunsail_0.DIR_surface_w+...
        FLATE_cone_sunsail_0.DDN_surface_w+...
        FLATE_cone_sunsail_0.RUP_surface_w;%solseildekket himmel,høy transmittans
    E_cone_sunsail_0(isnan(E_cone_sunsail_0))=0;
    
    
    %flat reflektans på 10 % for trær og vegger
    Rwalls(:,1)=x;
    Rwalls(:,2)=0.1;
    
    %ingen fri horisont, clod=0, solseil totalblokkerer nesten, transmittans lik reflektans til
    %bygninger (10 %)
    FLATE_cone_lowT_sunsail_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_lowT_sunsail_0=FLATE_cone_lowT_sunsail_0.DIR_surface_w+...
        FLATE_cone_lowT_sunsail_0.DDN_surface_w+...
        FLATE_cone_lowT_sunsail_0.RUP_surface_w;%solseildekket himmel,lav transmittans
    E_cone_lowT_sunsail_0(isnan(E_cone_lowT_sunsail_0))=0;
    
    
    %Vektet etter transmittans, clod
    FLATE_cone_sunsail = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Tsail);
    E_cone_sunsail=FLATE_cone_sunsail.DIR_surface_w+...
        FLATE_cone_sunsail.DDN_surface_w+...
        FLATE_cone_sunsail.RUP_surface_w;%åpen, fri himmel
    E_cone_sunsail(isnan(E_cone_sunsail))=0;
    
    %Vektet etter bygningers reflektans på 10 prosent:
    
    %clod=0
    FLATE_cone_walls_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_walls_0=FLATE_cone_walls_0.DIR_surface_w.*FLATE_cone_walls_0.ShapeDir+...
        FLATE_cone_walls_0.DDN_surface_w.*FLATE_cone_walls_0.ShapeDn+...
        FLATE_cone_walls_0.RUP_surface_w.*FLATE_cone_walls_0.ShapeUp;%trær og bygninger
    E_cone_walls_0(isnan(E_cone_walls_0))=0;
    
    %aktuell clod
    FLATE_cone_walls = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_walls=FLATE_cone_walls.DIR_surface_w+...
        FLATE_cone_walls.DDN_surface_w+...
        FLATE_cone_walls.RUP_surface_w;%trær og bygninger
    E_cone_walls(isnan(E_cone_walls))=0;
    
    figure(111)
    plot(DecH,E_cone_free_0(:,dayIm),'k--','LineWidth',1.5)
    hold on
    plot(DecH,E_cone_sunsail_0(:,dayIm),'k:','LineWidth',1.5)
    %%plot(DecH(1:30:end),E_cone_lowT_sunsail_0(1:30:end,dayIm),'ko-','LineWidth',1.5)
    hold off
    
    %Nå: beregn hver av tre irradianskomponenter, vekt med skyview,sailview og wallview og legg sammen
    
    %     figure,
    %     plot(DecH,pix_bane_avg,'ro-')
    %     hold on
    %     plot(DecH(psail),pix_bane_avg(psail),'k.-')
    %     hold off
    
    %Beregn irradians for solbane når skyline blokkerer direktekomponenten og horisonten har begrenset skyview:
    %pix_bane_avg %piksler med verdi 1 betyr åpen sol, 0 er blokkert sol
    
    %Først for klarvær, clod 0:
    %Husk at når himmelen er oppstykket, må RUP beregnes på grunnlag av horisontalkomponentene til
    %direkte og diffust ned, multiplisert med albedo og shapefaktor for opp-komponenten
    Dir_cone_obstructed_0=(FLATE_cone_free_0.DIR_surface_w(:,dayIm).*pix_bane_sky + ...
        FLATE_cone_sunsail_0.DIR_surface_w(:,dayIm).*pix_bane_sunsail); %bidrag fra vegger og trær er innbakt i pix_bane_sky
    
    Dn_cone_obstructed_0=(FLATE_cone_free_0.DDN_surface_w(:,dayIm)*svf + ...
        FLATE_cone_sunsail_0.DDN_surface_w(:,dayIm)*vfsunsails + ...
        FLATE_cone_walls_0.DDN_surface_w(:,dayIm)*vfwalls);
    
    R_cone_obstructed_0 =((FLATE_cone_free_0.EDIR365(:,dayIm).*pix_bane_sky + FLATE_cone_sunsail_0.EDIR365(:,dayIm).*pix_bane_sunsail + ...
        FLATE_cone_free_0.EDN365(:,dayIm)*svf + FLATE_cone_sunsail_0.EDN365(:,dayIm)*vfsunsails + FLATE_cone_walls_0.EDN365(:,dayIm)*vfwalls )*...
        FLATE_cone_free_0.Albedo).*FLATE_cone_free_0.ShapeUp(:,dayIm);
    R_cone_obstructed_0(isnan(R_cone_obstructed_0))=0;
    
    E_cone_obstructed_0=Dir_cone_obstructed_0+Dn_cone_obstructed_0+R_cone_obstructed_0;
    E_cone_obstructed_0(isnan(E_cone_obstructed_0))=0;
    
    %så for aktuell clod
    Dir_cone_obstructed=(FLATE_cone_free.DIR_surface_w(:,dayIm).*pix_bane_sky + ...
        FLATE_cone_sunsail.DIR_surface_w(:,dayIm).*pix_bane_sunsail); %bidrag fra vegger og trær er innbakt i pix_bane_sky
    
    Dn_cone_obstructed=(FLATE_cone_free.DDN_surface_w(:,dayIm)*svf + ...
        FLATE_cone_sunsail.DDN_surface_w(:,dayIm)*vfsunsails + ...
        FLATE_cone_walls.DDN_surface_w(:,dayIm)*vfwalls);
    
    R_cone_obstructed = (FLATE_cone_free.EDIR365(:,dayIm).*pix_bane_sky + FLATE_cone_sunsail.EDIR365(:,dayIm).*pix_bane_sunsail + ...
        FLATE_cone_free.EDN365(:,dayIm)*svf + FLATE_cone_sunsail.EDN365(:,dayIm)*vfsunsails + FLATE_cone_walls.EDN365(:,dayIm)*vfwalls )*...
        FLATE_cone_free.Albedo.*FLATE_cone_free.ShapeUp(:,dayIm);
    R_cone_obstructed(isnan(R_cone_obstructed))=0;
    
    E_cone_obstructed=Dir_cone_obstructed+Dn_cone_obstructed+R_cone_obstructed;
    E_cone_obstructed(isnan(E_cone_obstructed))=0;
    
    
    plot_title=0;
    figure(111)
    hold on
    plot(DecH,E_cone_obstructed_0,'k-','LineWidth',1.5)
    plot(DecH(1:30:end),E_cone_obstructed(1:30:end),'ko-','LineWidth',1)
    hold off
    %legend('Non-obstructed','Sunsail-covered, high-T','Sunsail-covered, T 10 %','Obstructed, clearsky',sprintf('obstructed, overcast'))
    legend('Non-obstructed','Sunsail-covered, high-T','Obstructed, clearsky',sprintf('obstructed, overcast'))
    legend('boxoff')
    xlabel('Time')
    
    %%ylabel(sprintf('%s-irradiance,W/m^2',actionspec))
    ylabel(sprintf('PE,W/m^2',actionspec))
    %title(sprintf('Irradiance on a cone and horizontal plane\nDosefraction: Cone: %2.0f %%, Horizontal: %2.0f %%',dose_fraksjon_cone*100,dose_fraksjon_horizontal*100))
    if plot_title
        title(sprintf('Irradiance on a cone in a free and obstructed sky environment'))
    end
    xlim([4 19])
    ylim([0 max(ceil(E_cone_free_0(:,dayIm)))])
    
    
    trapz(DecH,E_cone_obstructed_0)/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio med alle obstruksjoner/åpen horisont
    trapz(DecH,E_cone_sunsail_0(:,dayIm))/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio under et heldekkende solseil/åpen horisont
    trapz(DecH,E_cone_sunsail(:,dayIm))/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio under et heldekkende solseil med overskyet/åpen klarværs horisont
    
    
    pD=find(DecH>9.45 & DecH<=10.45);%direkte sol
    pS=find(DecH>10.57 & DecH<=11.57);%skygget av solseil
    trapz(DecH(pD),E_cone_obstructed_0(pD))/trapz(DecH(pD),E_cone_free_0(pD,dayIm))
    trapz(DecH(pS),E_cone_obstructed_0(pS))/trapz(DecH(pS),E_cone_free_0(pS,dayIm))
    
    hFig=figure(111);
    hgexport(hFig,'-clipboard')
    
    %Gjenta for CIE-vektet UV:
    %actionspec='EPP';
    actionspec='CIE';
    if strcmp(actionspec,'CIE')
        scale_fact=40;
    else
        scale_fact=1;
    end
    
    clear FLATE_cone_free FLATE_cone_sunsail FLATE_cone_walls
    
    FLATE_cone_free_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,T_clearsky);
    E_cone_free_0=FLATE_cone_free_0.DIR_surface_w+...
        FLATE_cone_free_0.DDN_surface_w+...
        FLATE_cone_free_0.RUP_surface_w;%åpen, fri himmel
    E_cone_free_0(isnan(E_cone_free_0))=0;
    
    %fri horisont, for gitt clod
    FLATE_cone_free = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,T_clearsky);
    E_cone_free=FLATE_cone_free.DIR_surface_w +...
        FLATE_cone_free.DDN_surface_w+...
        FLATE_cone_free.RUP_surface_w;%åpen, fri himmel
    E_cone_free(isnan(E_cone_free))=0;
    
    %ingen fri horisont, clod=0, Vektet etter transmittans,
    FLATE_cone_sunsail_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Tsail);
    E_cone_sunsail_0=FLATE_cone_sunsail_0.DIR_surface_w+...
        FLATE_cone_sunsail_0.DDN_surface_w+...
        FLATE_cone_sunsail_0.RUP_surface_w;%solseildekket himmel,høy transmittans
    E_cone_sunsail_0(isnan(E_cone_sunsail_0))=0;
    
    
    %ingen fri horisont, clod=0, solseil totalblokkerer nesten, transmittans lik reflektans til
    %bygninger (10 %)
    FLATE_cone_lowT_sunsail_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_lowT_sunsail_0=FLATE_cone_lowT_sunsail_0.DIR_surface_w+...
        FLATE_cone_lowT_sunsail_0.DDN_surface_w+...
        FLATE_cone_lowT_sunsail_0.RUP_surface_w;%solseildekket himmel,lav transmittans
    E_cone_lowT_sunsail_0(isnan(E_cone_lowT_sunsail_0))=0;
    
    
    %Vektet etter transmittans, clod
    FLATE_cone_sunsail = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Tsail);
    E_cone_sunsail=FLATE_cone_sunsail.DIR_surface_w+...
        FLATE_cone_sunsail.DDN_surface_w+...
        FLATE_cone_sunsail.RUP_surface_w;%åpen, fri himmel
    E_cone_sunsail(isnan(E_cone_sunsail))=0;
    
    
    %clod=0
    FLATE_cone_walls_0 = porfyri_returner_horisontalkomponenter('clear_sky',0,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_walls_0=FLATE_cone_walls_0.DIR_surface_w.*FLATE_cone_walls_0.ShapeDir+...
        FLATE_cone_walls_0.DDN_surface_w.*FLATE_cone_walls_0.ShapeDn+...
        FLATE_cone_walls_0.RUP_surface_w.*FLATE_cone_walls_0.ShapeUp;%trær og bygninger
    E_cone_walls_0(isnan(E_cone_walls_0))=0;
    
    %aktuell clod
    FLATE_cone_walls = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_walls=FLATE_cone_walls.DIR_surface_w+...
        FLATE_cone_walls.DDN_surface_w+...
        FLATE_cone_walls.RUP_surface_w;%trær og bygninger
    E_cone_walls(isnan(E_cone_walls))=0;
    
    figure(112)
    plot(DecH,scale_fact*E_cone_free_0(:,dayIm),'k--','LineWidth',1.5)
    hold on
    plot(DecH,scale_fact*E_cone_sunsail_0(:,dayIm),'k:','LineWidth',1.5)
    %%plot(DecH(1:30:end),scale_fact*E_cone_lowT_sunsail_0(1:30:end,dayIm),'ko-','LineWidth',1.5)
    hold off
    
    %Nå: beregn hver av tre irradianskomponenter, vekt med skyview,sailview og wallview og legg sammen
    
    %     figure,
    %     plot(DecH,pix_bane_avg,'ro-')
    %     hold on
    %     plot(DecH(psail),pix_bane_avg(psail),'k.-')
    %     hold off
    
    %Beregn irradians for solbane når skyline blokkerer direktekomponenten og horisonten har begrenset skyview:
    %pix_bane_avg %piksler med verdi 1 betyr åpen sol, 0 er blokkert sol
    
    %Først for klarvær, clod 0:
    %Husk at når himmelen er oppstykket, må RUP beregnes på grunnlag av horisontalkomponentene til
    %direkte og diffust ned, multiplisert med albedo og shapefaktor for opp-komponenten
    Dir_cone_obstructed_0=(FLATE_cone_free_0.DIR_surface_w(:,dayIm).*pix_bane_sky + ...
        FLATE_cone_sunsail_0.DIR_surface_w(:,dayIm).*pix_bane_sunsail); %bidrag fra vegger og trær er innbakt i pix_bane_sky
    
    Dn_cone_obstructed_0=(FLATE_cone_free_0.DDN_surface_w(:,dayIm)*svf + ...
        FLATE_cone_sunsail_0.DDN_surface_w(:,dayIm)*vfsunsails + ...
        FLATE_cone_walls_0.DDN_surface_w(:,dayIm)*vfwalls);
    
    R_cone_obstructed_0 =((FLATE_cone_free_0.EDIR365(:,dayIm).*pix_bane_sky + FLATE_cone_sunsail_0.EDIR365(:,dayIm).*pix_bane_sunsail + ...
        FLATE_cone_free_0.EDN365(:,dayIm)*svf + FLATE_cone_sunsail_0.EDN365(:,dayIm)*vfsunsails + FLATE_cone_walls_0.EDN365(:,dayIm)*vfwalls )*...
        FLATE_cone_free_0.Albedo).*FLATE_cone_free_0.ShapeUp(:,dayIm);
    R_cone_obstructed_0(isnan(R_cone_obstructed_0))=0;
    
    E_cone_obstructed_0=Dir_cone_obstructed_0+Dn_cone_obstructed_0+R_cone_obstructed_0;
    E_cone_obstructed_0(isnan(E_cone_obstructed_0))=0;
    
    %så for aktuell clod
    Dir_cone_obstructed=(FLATE_cone_free.DIR_surface_w(:,dayIm).*pix_bane_sky + ...
        FLATE_cone_sunsail.DIR_surface_w(:,dayIm).*pix_bane_sunsail); %bidrag fra vegger og trær er innbakt i pix_bane_sky
    
    Dn_cone_obstructed=(FLATE_cone_free.DDN_surface_w(:,dayIm)*svf + ...
        FLATE_cone_sunsail.DDN_surface_w(:,dayIm)*vfsunsails + ...
        FLATE_cone_walls.DDN_surface_w(:,dayIm)*vfwalls);
    
    R_cone_obstructed = (FLATE_cone_free.EDIR365(:,dayIm).*pix_bane_sky + FLATE_cone_sunsail.EDIR365(:,dayIm).*pix_bane_sunsail + ...
        FLATE_cone_free.EDN365(:,dayIm)*svf + FLATE_cone_sunsail.EDN365(:,dayIm)*vfsunsails + FLATE_cone_walls.EDN365(:,dayIm)*vfwalls )*...
        FLATE_cone_free.Albedo.*FLATE_cone_free.ShapeUp(:,dayIm);
    R_cone_obstructed(isnan(R_cone_obstructed))=0;
    
    E_cone_obstructed=Dir_cone_obstructed+Dn_cone_obstructed+R_cone_obstructed;
    E_cone_obstructed(isnan(E_cone_obstructed))=0;
    
    
    figure(112)
    hold on
    plot(DecH,scale_fact*E_cone_obstructed_0,'k-','LineWidth',1.5)    
    plot(DecH(1:30:end),scale_fact*E_cone_obstructed(1:30:end),'ko-','LineWidth',1)
    hold off
    %legend('Free horizon','Sunsail covered','obstructed clear sky',sprintf('obstructed CLOD %d',clod))
    legend('Non-obstructed','Sunsail-covered, high-T','Sunsail-covered, T 10 %','Obstructed, clearsky',sprintf('obstructed, overcast'))
    legend('boxoff')
    xlabel('Time')
    %%ylabel(sprintf('%s x %d,W/m^2',actionspec,scale_fact))
    ylabel(sprintf('UVery x %d',scale_fact))
    %title(sprintf('Irradiance on a cone and horizontal plane\nDosefraction: Cone: %2.0f %%, Horizontal: %2.0f %%',dose_fraksjon_cone*100,dose_fraksjon_horizontal*100))
    if plot_title
        title(sprintf('Irradiance on a cone in a free and obstructed sky environment'))
    end
    xlim([4 19])
    %ylim([0 max(ceil(scale_fact*E_cone_free_0(:,dayIm)))])
    ylim([0 2.5])
    
    hFig=figure(112);
    hgexport(hFig,'-clipboard')
    
    trapz(DecH,E_cone_obstructed_0)/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio med alle obstruksjoner/åpen horisont
    trapz(DecH,E_cone_sunsail_0(:,dayIm))/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio under et heldekkende solseil/åpen horisont
    trapz(DecH,E_cone_sunsail(:,dayIm))/trapz(DecH,E_cone_free_0(:,dayIm))%dagsratio under et heldekkende solseil med overskyet/åpen klarværs horisont
    
    
    pD=find(DecH>9.45 & DecH<=10.45);%direkte sol
    pS=find(DecH>10.36 & DecH<=11.36);%skygget av solseil
    trapz(DecH(pD),E_cone_obstructed_0(pD))/trapz(DecH(pD),E_cone_free_0(pD,dayIm))
    trapz(DecH(pS),E_cone_obstructed_0(pS))/trapz(DecH(pS),E_cone_free_0(pS,dayIm))
    
    %Gammelt    ......
    %vekt deretter bidrag fra hver, avhengig av hvor solbanen går
    FLATE_cone_free = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,T_clearsky);
    E_cone_free=FLATE_cone_free.DIR_surface_w+FLATE_cone_free.DDN_surface_w+FLATE_cone_free.RUP_surface_w;%åpen, fri himmel
    E_cone_free(isnan(E_cone_free))=0;
    
    %Vektet etter transmittans
    FLATE_cone_sunsail = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Tsail);
    E_cone_sunsail=FLATE_cone_sunsail.DIR_surface_w+FLATE_cone_sunsail.DDN_surface_w+FLATE_cone_sunsail.RUP_surface_w;%åpen, fri himmel
    E_cone_sunsail(isnan(E_cone_sunsail))=0;
    
    %     %Vektet etter bygningers reflektans på 10 prosent:
    %     Rwalls(:,1)=x;
    %     Rwalls(:,2)=0.1;
    
    FLATE_cone_walls = porfyri_returner_horisontalkomponenter(skytype,clod,350,ground_albedo,lat,lon,surface_type,beta,0,actionspec,Rwalls);
    E_cone_walls=FLATE_cone_walls.DIR_surface_w+FLATE_cone_walls.DDN_surface_w+FLATE_cone_walls.RUP_surface_w;%trær og bygninger
    E_cone_walls(isnan(E_cone_walls))=0;
    
    figure(112)
    plot(DecH,scale_fact*E_cone_free(:,dayIm),'r-','LineWidth',1.5)
    hold on
    plot(DecH,scale_fact*E_cone_sunsail(:,dayIm),'b-','LineWidth',1.5)
    hold off
    
    %Nå: beregn hver av tre irradianskomponenter, vekt med skyview,sailview og wallview og legg sammen
    
    %     figure,
    %     plot(DecH,pix_bane_avg,'ro-')
    %     hold on
    %     plot(DecH(psail),pix_bane_avg(psail),'k.-')
    %     hold off
    
    %Beregn irradians for solbane når skyline blokkerer direktekomponenten og horisonten har begrenset skyview:
    %pix_bane_avg %piksler med verdi 1 betyr åpen sol, 0 er blokkert sol
    
    Dir_cone_obstructed=FLATE_cone_free.DIR_surface_w(:,dayIm).*pix_bane_sky + ...
        FLATE_cone_sunsail.DIR_surface_w(:,dayIm).*pix_bane_sunsail; %bidrag fra vegger og trær er innbakt i pix_bane_sky
    
    Dn_cone_obstructed=FLATE_cone_free.DDN_surface_w(:,dayIm)*svf + ...
        FLATE_cone_sunsail.DDN_surface_w(:,dayIm)*vfsunsails + ...
        FLATE_cone_walls.DDN_surface_w(:,dayIm)*vfwalls;
    
    R_cone_obstructed =((FLATE_cone_free.EDIR365(:,dayIm).*pix_bane_sky + FLATE_cone_sunsail.EDIR365(:,dayIm).*pix_bane_sunsail + ...
        FLATE_cone_free.EDN365(:,dayIm)*svf + FLATE_cone_sunsail.EDN365(:,dayIm)*vfsunsails + FLATE_cone_walls.EDN365(:,dayIm)*vfwalls )*...
        FLATE_cone_free.Albedo).*FLATE_cone_free.ShapeUp(:,dayIm);
    R_cone_obstructed(isnan(R_cone_obstructed))=0;
    
    E_cone_obstructed=Dir_cone_obstructed+Dn_cone_obstructed+R_cone_obstructed;
    E_cone_obstructed(isnan(E_cone_obstructed))=0;
    
    figure(112)
    hold on
    plot(DecH,scale_fact*E_cone_obstructed,'k-','LineWidth',1.5)
    hold off
    legend('Free horizon','Sunsail','obstructed sky')
    xlabel('Time')
    %ylabel('CIE-irradiance,W/m^2')
    ylabel(sprintf('%s x %d,W/m^2',actionspec,scale_fact))
    %title(sprintf('Irradiance on a cone and horizontal plane\nDosefraction: Cone: %2.0f %%, Horizontal: %2.0f %%',dose_fraksjon_cone*100,dose_fraksjon_horizontal*100))
    title(sprintf('Irradiance on a cone in a free and obstructed sky environment'))
    xlim([4 19])
    ylim([0 3])
    
    hFig=figure(112);
    hgexport(hFig,'-clipboard')
    
    
    %Gammelt hit
    
    
    
    
end