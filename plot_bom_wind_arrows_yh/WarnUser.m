function [fullFileName]=WarnUser(ME)
  % save error file when catch statement is executed
%   e.g. try; expression... catch ME  ;[fullFileName]=WarnUser(ME); end
%   yasha march 2019

%  fullFileName = ['ErrorLog_' datestr(now,'yyyymmddHHMM') '.txt'];
  fullFileName = ['ErrorLog.txt'];
  fid = fopen(fullFileName, 'at');

fprintf(fid,'There was an error! on: %s\n',datestr(now));
fprintf(fid,'The identifier was: %s\n ',ME.identifier);
        fprintf(fid,'The message was: %s\n',ME.message);
        fprintf(fid,'file: %s\n',ME.stack.file);
        fprintf(fid,'name: %s \n',ME.stack.name);
        fprintf(fid,'line: %s \n',ME.stack.line);

  
  % Open the Error Log file for appending.
 
  fclose(fid);
  return;