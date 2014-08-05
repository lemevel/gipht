function printps(psfilename,ghandle)
% write current graphics window to a PS file

if nargin < 1
    psfilename = mfilename;
end

if numel(strfind(psfilename,'.ps')) <= 0
    psfilename = sprintf('%s.ps',psfilename);
end

% print current figure to ps file name with 1200 DPI and TIFF
t0=pwd;
%t1=strrep(psfilename,'_','\_');
t1=psfilename;
%t2=date;
t2 = datestr(now,31); %31             'yyyy-mm-dd HH:MM:SS'    2000-03-01 15:45:17 
tu=getenv('USER');
t3=sprintf('%s %s %s %s',t1,t0,t2,tu);
t4= strrep(t3,'_','\_');

% label with file name, date and time
subplot('Position',[0 0 10 0.05],'Units','Centimeters');
axis off
text(0,0,t4...
    ,'Units','Centimeters'...
    ,'VerticalAlignment','Bottom'...
    ,'HorizontalAlignment','Left'...
    ,'Clipping','off'...
    ,'FontName','Courier','FontSize',9 ...
    ,'Rotation',0);

% Label the figure
% Does not work on Hengill
% mycomputer = computer;
% if strcmp(mycomputer, 'GLNXA64')==0
% http://www.mathworks.com/support/solutions/en/data/1-703XB9/?solution=1-703XB9
% Subject:
% 
% Is it possible to programmatically check whether MATLAB has been started
% with the "-nodisplay" option? Problem Description:
% 
% I have a program that needs to behave differently depending on whether
% MATLAB has a display or not. However, I cannot determine a way to
% programmatically check this.
% 
% The reason I need to do this is that I sometimes start MATLAB in batch
% mode from a shell script, for testing my programs.
% 
% Solution:
% 
% The ability to programmatically check whether MATLAB has a display is not
% available in MATLAB 7.6 (R2008a).
% 
% As a workaround, you can do one of the following:
% 
% 1. Manually set or unset an environment variable in your shell script
% that launches MATLAB, so you can query it from inside MATLAB. For
% instance, setting the ISDISPLAY environment variable in your shell script
% before launching MATLAB:
% 
% 
% setenv ISDISPLAY no matlab -r "foo; quit"
% 
% (here shown with C Shell syntax), means that the MATLAB command
% getenv('ISDISPLAY') will return the string 'no' in any MATLAB processes
% that are spawned from this shell.
% 
% 2. Query the "ScreenSize" property of the root object inside MATLAB:
% get(0, 'ScreenSize')
% When there is no display, this returns [1 1 1 1] instead of an actual screen size. However, this relies on behavior that isn't actually specified (by the doc, for instance) to work in any particular way, so may be subject to change in the future. If you were going to use this many times, it might be wise to wrap it in a function (e.g. create an "isdisplay.m" function file), so you can easily change the implementation in the future, if needed. (This method worked as of MATLAB R2008a.) 
% ss4 = get(0, 'ScreenSize');


% if ss4 == [1 1 1 1 ]
%    print(gcf,strrep(psfilename,'ps','ps'),'-dpsc2','-r1200'); % print PS if no display  
% else
%    print(gcf,psfilename,'-dps','-r1200'); % otherwise, print PS
% end

if exist('ghandle','var') == 1
    print(ghandle,psfilename,'-dpsc','-r1200'); % otherwise, print PS
else
    print(gcf,psfilename,'-dpsc','-r1200'); % otherwise, print PS
end

return


