function out3(m,textstr,fp)
% Usage: out3(m,textstr,fp)
% 
% Purpose: Print formated 3x3 matrix with header text
%
% Input: m       -- 3x3 matrix
%        textstr -- text string
%        fp      -- file pointer
%
% Output: Text written to fp

    sep(1:63) = '*';
    % Print header
    fprintf(fp,'\n%s\n%s\n%s\n',sep,textstr,sep);
    % Print matrix
    fprintf(fp,'% 4.3e\t % 4.3e\t % 4.3e\n',m'); 
end