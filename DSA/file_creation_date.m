function [file_date,filestring,filetype]=file_creation_date(file)
%returnerer filecreation date. Bjørn Johnsen

[dummy1,filestring,filetype]= fileparts(file);
 allfiles = dir(dummy1)
 filenames = {allfiles(:).name};
 [~,idx] = ismember(upper(sprintf('%s%s',filestring,filetype)),filenames)%må ha store bokstaver
 if idx==0
      [~,idx] = ismember(sprintf('%s%s',filestring,filetype),filenames)%har små bokstaver
      if idx==0
          printf('bildefila mangler')
      end
 end
 %file_date=datevec(allfiles(idx).datenum);
 file_date = datevec(imfinfo(char(fullfile(dummy1,filenames(idx)))).DateTime, 'yyyy:mm:dd HH:MM:SS');
 %clear dummy1 filestring filetype allfiles filenames