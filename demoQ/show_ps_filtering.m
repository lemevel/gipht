close all;
clear all;
giphtpath
printfun = 'printjpg'
format compact

% number of columns in each interferogram
ncols =         1420
% number of lines in each interferogram
nrows =         1230


pha  =double(read_pha('pha_11176_21540_ort.pha',ncols))/256.;
psp  =double(read_pha('psp_11176_21540_ort.pha',ncols))/256.;

pha2=pha(400:600,600:800);
psp2=psp(400:600,600:800);

figure;imagesc(pha2);colorbar;cmapblackzero;
title('before filtering: wrapped phase [cycles]');
xlabel('column index');
ylabel('row index');
feval(printfun,sprintf('pha.jpg',mfilename));

figure;imagesc(psp2);colorbar;cmapblackzero;
title('after filtering: wrapped phase [cycles]');
xlabel('column index');
ylabel('row index');
feval(printfun,sprintf('psp.jpg',mfilename));

