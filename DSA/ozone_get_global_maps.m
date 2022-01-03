function ozone_get_global_maps(year)
%Laster ned hele årskatalogen med global, griddete ozondata

%clean
%year=2019;
%year=2020;
%year=2021;

ws=pwd; %'E:\Data\Matlab-programs\m-filer'

wget_folder='c:\users\bjohnsen' %wget er installert her

cd (wget_folder)

wget_cmd=sprintf('wget -r --no-parent --cut-dirs=2 https://ozonewatch.gsfc.nasa.gov/data/omi/Y%d/',year);

tic
dos(wget_cmd)
toc

%alle dalig global-maps lastes ned til C:\Users\bjohnsen\ozonewatch.gsfc.nasa.gov\Y2019

%flytt mappen herfra og over til N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\ozone_toms_omi

from_folder=sprintf('%s\\ozonewatch.gsfc.nasa.gov\\Y%d',wget_folder,year);
to_folder='N:\data_fra_photon\LAB\OLAB\UV_NETT\Ozone\ozone_toms_omi';

movefile(from_folder,to_folder)

cd(ws)



% wget -r --no-parent --cut-dirs=2 https://ozonewatch.gsfc.nasa.gov/data/omi/Y2019/
% filer havner i C:\Users\bjohnsen\ozonewatch.gsfc.nasa.gov\Y2019,
% samtidig som vi unngår /data/omi undermapper.

