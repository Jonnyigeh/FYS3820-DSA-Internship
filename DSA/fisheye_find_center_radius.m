function [im_cntrx,im_cntry,im_rad,im_radH]=fisheye_find_center_radius(rgb)


%les inn fisheye bilde - bestem sentrum og radius til bildekant og radius til 90 graders horisont

%imfile='Z:\Seksjon Overvåkning og miljødata\83054 overvåkning UV og soling\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2342.jpg';
%im_cal='Z:\Seksjon Overvåkning og miljødata\83054 overvåkning UV og soling\BILDERogILLUSTRASJONER\Bilder\2016\Skyggeprosjekt_bilder\Fisheye_originalbilder\DSCN2333.jpg';
%rgb = imread(im_cal);
im_gray=rgb2gray(rgb);
[ny,nx]=size(rgb(:,:,1));
o=round([nx,ny]./[2 2]);

M = size(rgb,1);%samme
N = size(rgb,2);

% figure(20)
% ax(1) = subplot(2,2,1);
% image(rgb); title('RGB image')
% axis off          % Remove axis ticks and numbers
% axis image
% ax(2) = subplot(2,2,2);
% %im = mean(rgb,3);
% im = rgb(:,:,1);
% %im = rgb(:,:,2);
% %im = rgb(:,:,3);
% image(im); title('Intensity Heat Map')
% %colormap(Gray)
% axis off          % Remove axis ticks and numbers
% axis image        % Set aspect ratio to obtain square pixels
% %axis square
% set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
% % ax(2) = subplot(122);
% % im = mean(rgb,3);
% % image(im); title('Intensity Heat Map')
% % colormap(hot(256))
% % linkaxes(ax,'xy')
% % axis(ax,'image')
% im_gray=rgb2gray(rgb);
% ax(3) = subplot(2,2,3);
% %image(im_gray); title('gray image')
% imshow(im_gray); title('gray image')
% axis off          % Remove axis ticks and numbers
% axis image

figure(1)
imshow(rgb);

figure(2)
semilogy([1:nx],im_gray(o(2),:),'k.-'),grid on
axis([1 nx 0 255])

im_wx=find(im_gray(o(2),:)>5);
im_cntrx=round([im_wx(1)+im_wx(end)]/2);

% figure(1)
% hold on
% % line([im_cntrx im_cntrx],[1 ny],'Color','w','LineWidth',4)
% % line([1500 1500],[1 1000],'Color','w')
%
% for k = 1:25:N
%     x = [k k];
%     y = [1 M];
%     plot(x,y,'Color','w','LineStyle','-');
%     plot(x,y,'Color','k','LineStyle',':');
% end
%
% hold off
%
% figure(3)
% semilogy([1:ny],im_gray(:,im_cntrx),'k.-'),grid on
% axis([1 ny 0 255])
%
% figure(4)
% semilogy([1:ny],squeeze(rgb(:,im_cntrx,1)),'k.-'),grid on
% axis([1 ny 0 255])
% %ser at lasertopp på R-kanal ligger på kolonne 1610=90-graders horisont
%
% figure(5)
% semilogy([1:ny],squeeze(rgb(:,im_cntrx,2)),'k.-'),grid on
% axis([1 ny 0 255])
%
figure(6)
semilogy([1:ny],squeeze(rgb(:,im_cntrx,3)),'k.-'),grid on
axis([1 ny 0 255])

im_wy=find(im_gray(:,im_cntrx)>10);
im_cntry=round([im_wy(1)+im_wy(end)]/2);

im_radx=round([im_wx(end)-im_cntrx]);
im_rady=round([im_wy(end)-im_cntry]);
im_rad_90=round([1610-im_cntry]);

im_rad=im_radx;
im_radH=im_rad_90;
r=im_radx; % radius
r90=im_rad_90;
center=[im_cntrx,im_cntry]; % center cordinates
t=-pi:0.001:pi;
x=round(r*cos(t)+center(1));
y=round(r*sin(t)+center(2));

xH=round(r90*cos(t)+center(1));
yH=round(r90*sin(t)+center(2));

figure(1)
hold on
plot(im_cntrx,im_cntry,'wo','MarkerSize',10)
plot([im_cntrx im_cntrx],[1 ny],'Color','w','LineWidth',1)
plot([1 nx],[im_cntry im_cntry],'Color','w','LineWidth',1)

plot(x,y,'w--','LineWidth',0.5),
plot(xH,yH,'b--','LineWidth',0.1),

%plot([1 nx],[o(2) o(2)],'w-')
%plot([o(1) o(1)],[1 ny],'b--')

%plot(o(2),o(1),'ro','MarkerSize',20)
% ax = copyobj(gca, gcf);
% set(ax,'color','none','xgrid','off', 'xcolor','w', 'ygrid','off', 'ycolor','w')%NB: denne
% forhindrer plotting
hold off

figure(1)
hold on

for k = 1:100:N
    x = [k k];
    y = [1 M];
    plot(x,y,'Color','w','LineStyle','-');
    plot(x,y,'Color','k','LineStyle',':');
end

hold off
