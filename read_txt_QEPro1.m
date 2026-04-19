% Read the QEPro data from a txt file
% Input: 
%   filename        - filename
%
% Output: 
%   table           - data matrix extracted from the input file
% 
% 03/06/2020 (Nghi Truong)

function table = read_txt_QEPro1(filename)

fid = fopen(filename,'r');
str = textscan(fid,'%s','Delimiter','\r');
str = str{1};
fclose(fid);

for i=19:length(str)
    table(i-18,:) = str2double(strsplit(str{i}));
end