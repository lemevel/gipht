function giphtpath
% set path for GIPHT

mycomputer = computer;
% % Intel Mac
if strcmp(mycomputer, 'MACI64')==1
    setenv('GIPHT_HOME','/Users/feigl/gipht');
elseif strcmp(mycomputer, 'GLNXA64')==1
    if exist('/usr1/gipht','dir') == 7
        % standard version
        %setenv('GIPHT_HOME','/usr1/gipht');
        % development version
        %setenv('GIPHT_HOME','/usr1/feigl/GIPhT2.4');
        % development version
        setenv('GIPHT_HOME','/usr1/feigl/gipht');
    else
        % default
        home = getenv('HOME');
        %setenv('GIPHT_HOME',strcat(home,'/GIPhT2.4'));
        setenv('GIPHT_HOME','./');
    end
else
    warning(sprintf('Unknown computer %s\n. Path left unchanged.\n',mycomputer));
    setenv('GIPHT_HOME','./');
end
fprintf(1,'Environment variable GIPHT_HOME is set to %s\n',getenv('GIPHT_HOME'));
p = [strcat(getenv('GIPHT_HOME'),'/src:'),strcat(getenv('GIPHT_HOME'),'/extern')];
addpath(p,'-BEGIN');
% fprintf(1,'Matlab command search path is now: %s\n',path);
return
