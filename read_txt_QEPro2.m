% Read the QEPro data from a txt file
% Input: 
%   filename        - filename
%
% Output: 
%   table           - data matrix extracted from the input file
% 
% 03/06/2020 (Nghi Truong)

function table = read_txt_QEPro2(filename)

fid = fopen(filename,'r');
str = textscan(fid,'%s','Delimiter','\r');
str = str{1};
fclose(fid);

for i=14:length(str)
    table(i-13,:) = str2double(strsplit(str{i}));
end