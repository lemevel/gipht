%test_quadphase

% demostrate quadtree compression on wrapped phase
% lap11jul09 removed coherence references
% lap11jul09 added code to skip over signature NX & NY in .qls file
% lap14jul09 added code to test quadesahp reconstructor
close all;
clear all;
giphtpath
nf = 0;

% number of columns in each interferogram
 ncols =         1420
% number of lines in each interferogram
 nrows =         1230

% Now check the C version of the reconstructor
%[ssx,srx] = unix('../src/qls2pha.maci64 psp_11176_21540_ort.qls -o trx.pha -d 1');
[ssx,srx] = unix('../src/qls2pha.a64    psp_11176_21540_ort.qls -o trx.pha -d 1');
if ( ssx ~= 0 )
	txt=sprintf('FAILURE of quadphase reconstruction program\n====Reason====\n%s==============\n',srx);
    disp(txt);
end

%rrx=zeros(nrows,ncols,'int8');
%rrx=double(read_pha('trx.pha',ncols));
rrx=read_pha('trx.pha',ncols);
nf=nf+1;h(nf)=figure;imagesc(rrx);colorbar;cmapblackzero;
title('wrapped phase after quad-tree reconstruction by qls2pha (256 DN per cycle)')
xlabel('column index');ylabel('row index');

