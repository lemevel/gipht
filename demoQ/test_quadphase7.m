%test_quadphase

% demostrate quadtree compression on wrapped phase
% lap11jul09 removed coherence references
% lap11jul09 added code to skip over signature NX & NY in .qls file
% lap14jul09 added code to test quadesahp reconstructor
clear all; close all;
giphtpath
nf = 0;

% number of columns in each interferogram
 ncols =         1420
% number of lines in each interferogram 
 nrows =         1230

cmd1='../src/pha2qls.maci64 psp_11176_21540_ort.pha  1420 1230 qsp_11176_21540_ort.pha grx_11176_21540_ort.i2 gry_11176_21540_ort.i2 qsp_11176_21540_ort.qls'
% usage: pha2qls     input.pha                nx   ny   output.pha              output.qls              ithresh[8] minpix[16]
                     
%unix(sprintf('%s',cmd1))  % default values
%unix(sprintf('%s 25',cmd1)) % set max circular mean deviation to 25/256 of a cycle
unix(sprintf('%s  8 9',cmd1)) % set min number of OK pixels in patch to 9
%unix(sprintf('%s  16 25',cmd1)) % set min number of OK pixels in patch to 25
%unix(sprintf('%s psp_11176_21540_ort.qls 8 16 ../demoQ/coh_11176_21540_ort.oct 1',cmd1)) % use only pixels with coh > 1st percentile

p  =2.0*pi*double(read_pha('psp_11176_21540_ort.pha',ncols))/256.;
q  =2.0*pi*double(read_pha('qsp_11176_21540_ort.pha',ncols))/256.;
grx=2.0*pi*double(read_i2 ('grx_11176_21540_ort.i2', ncols))/256./256.;
%lap c=read_oct('coh_11176_21540_ort.oct',ncols);

%lap c=read_oct('coh_11176_21540_ort.oct',ncols);

nf=nf+1;h(nf)=figure;imagesc(p);colorbar;cmapblackzero;
title('wrapped phase (radians)');
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));
%printjpg('psp_11176_21540_ort.jpg');
nf=nf+1;h(nf)=figure;imagesc(q);colorbar;cmapblackzero;
title('wrapped phase after quad-tree partitioning (radians per pixel)')
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;imagesc(grx);colorbar;cmapblackzero;
title('x gradient after quad-tree partitioning (radians per pixel)')
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%printjpg('qsp_11176_21540_ort.jpg');
%lap nf=nf+1;h(nf)=figure;imagesc(c);colorbar;cmapblackzero;
%lap colormap('gray');
%lap title('coherence (0 to 255)')
%lap xlabel('column index');ylabel('row index');
%printjpg('coh_11176_21540_ort.jpg');

%lap ibad=find(c>0);c=c(ibad);medcoh=median(c)
%lap figure
%lap hist(double(reshape(c,numel(c),1)),256);
%lap xlabel('coherence low (0) to high (255)');
%lap ylabel('Number of pixels');
%lap title('Histogram');


%[im,jm,qp,nok,nnull,i1,i2,j1,j2,tr]=textread('qsp_11176_21540_ort.qls','%d%d%d%d%d%d%d%d%d%d','headerlines',1);
ilist = read_i2('qsp_11176_21540_ort.qls',6);

% write a program in C to expand .phalist into .pha file, as done in the
% following lines
whos ilist
nx=ilist(1,3); % Number of Cols
ny=ilist(1,4); % Number of Rows
if ( ncols ~= nx || nrows ~= ny )
    disp('Warning nx & nx dont match qls file header');
end

i1=ilist(2:end,1);  % Index to row of first pixel in patch
j1=ilist(2:end,2);  % Index to col of first pixel in patch
kw=ilist(2:end,3);  % width (and height) of square patch
% qp=double(ilist(2:end,4));  % ohase value coded -128 to +127 such that 256 DN = 1 cycle
% gx=double(ilist(2:end,5));  % X-ward gradient of phase value coded -128 to +127 such that 256^2 DN = 1 cycle
% gy=double(ilist(2:end,6));  % Y-ward gradient of phase value coded -128 to +127 such that 256^2 DN = 1 cycle
% % 2011-MAR-24 - GREAT BIG BUG - NOW fixed in pha2qls3.c
% qqp=2*pi*double(ilist(2:end,4))/256/256;% phase value
% qgx=2*pi*double(ilist(2:end,5))/256/256;  % X-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
% qgy=2*pi*double(ilist(2:end,6))/256/256;  % Y-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
qp=2*pi*double(ilist(2:end,4))/256./256.;  % phase value
gx=2*pi*double(ilist(2:end,5))/256./256.;  % X-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
gy=2*pi*double(ilist(2:end,6))/256./256.;  % Y-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel

%printjpg('rsp_11176_21540_ort.jpg');

% i2=i1+kw; % patches are square
% j2=j1+kw;
i2=i1+kw-1; % patches are square
j2=j1+kw-1;
%ibad = find(abs(qp) > 0);
%im=im(ibad);
%jm=jm(ibad);
%qp=qp(ibad);
%i1=i1(ibad);
%i2=i2(ibad);
%j1=j1(ibad);
%j2=j2(ibad);
min(qp)
max(qp)
mean(qp)
std(qp)
npatch = numel(qp)
figure;
nf=nf+1;h(nf)=figure;hist(double(colvec(qp)),64);axis tight;
xlabel('phase value (radians)');
ylabel('Number of phase values');
title('After quadtree partitioning');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
hist(double(colvec(gx)),64);axis tight;
xlabel('x gradient value (radians per pixel)');
ylabel('Number of phase values');
title('X-gradient');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

%nf=nf+1;h(nf)=figure;qqplotvonmises(double(colvec(gx))/256.,'X-gradient');
%nf=nf+1;h(nf)=figure;qqplot(double(colvec(gx))/2./pi);title('X-gradient');


nf=nf+1;h(nf)=figure;
hist(double(colvec(gy)),64);axis tight;
xlabel('y gradient value (radians per pixel)');
ylabel('Number of phase values');
title('y-gradient');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));
% nf=nf+1;h(nf)=figure;qqplotvonmises(double(colvec(gy))/256.0,'Y-gradient');
% nf=nf+1;h(nf)=figure;qqplot(double(colvec(gy))/2./pi);title('Y-gradient');

% convert indices from C to Fortran convention
i1 = i1+1;
i2 = i2+1;
j1 = j1+1;
j2 = j2+1;



kbad=find(i1<1)    ;i1(kbad)=1;
kbad=find(i1>ncols);i1(kbad)=ncols;
kbad=find(j1<1)    ;j1(kbad)=1;
kbad=find(j1>nrows);j1(kbad)=nrows;
fprintf(1,'Extrema of i1 %d %d\n',min(i1),max(i1));
fprintf(1,'Extrema of j1 %d %d\n',min(j1),max(j1));

kbad=find(i2<1)    ;i2(kbad)=1;
kbad=find(i2>ncols);i2(kbad)=ncols;
kbad=find(j2<1)    ;j2(kbad)=1;
kbad=find(j2>nrows);j2(kbad)=nrows;
fprintf(1,'Extrema of i2 %d %d\n',min(i2),max(i2));
fprintf(1,'Extrema of j2 %d %d\n',min(j2),max(j2));

%r=zeros(nrows,ncols,'int8');
r=zeros(nrows,ncols);
grxpatch=zeros(nrows,ncols);
disp('dimensions of r');size(r)
disp('dimensions of q');size(q)
%r=zeros(nrows,ncols);
for k=1:npatch
%     i0 = double(i1(k)+i2(k))/2.0;
%     j0 = double(j1(k)+j2(k))/2.0;
    i0 = double(i1(k));
    j0 = double(j1(k));
    for i=i1(k):i2(k)
        for j=j1(k):j2(k)          
            r(j,i)=qp(k) + (i-i0)*gx(k) + (j-j0)*gy(k);
            grxpatch(j,i) = gx(k);
        end
    end
end

r = rwrapm(r);
disp('dimensions of r');size(r)
disp('dimensions of q'); size(q)
% for k=1:npatch
%     r(j1(k):j2(k),i1(k):i2(k))=int8(qp(k));
%     if i2-i1 ~= j2-j1
%       fprintf(1,'%#4d %#4d %#4d %#4d\n',i1,i2,j1,j2)
%     end
% %   for j=j1(k)+1:j2(k)+1
% %      for i=i1(k)+1:i2(k)+1
% %         fprintf(1,'%#4d %#4d %#4d %#4d\n',i,j,k,qp(k));
% %         r(j,i)=qp(k);
% %      end
% %   end
% end



nf=nf+1;h(nf)=figure;
imagesc(r);colorbar;cmapblackzero;
title('wrapped phase after reconstructing quad-tree partitioning (radians)')
xlabel('column index');ylabel('row index');
%printjpg('rsp_11176_21540_ort.jpg');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
imagesc(grxpatch);colorbar;cmapblackzero;
title('X phase gradient from QLS list (radians per pixel)')
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

% reconstruction error in phase
e=r-p;
%e=rarcm(r,p);
iempty=find(r==0);
e(iempty)=0;
iempty=find(p==0);
e(iempty)=0;

nf=nf+1;h(nf)=figure;
hist(colvec(e));
xlabel('phase value (radians)');
ylabel('Number of pixels');
title('Reconstruction errors');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
imagesc(e);colorbar;%cmapblackzero;
title('reconstruction errors (radians)')
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

figure
ibad=find(abs(e)>2.0*pi/256.0);e=e(ibad);
hist(colvec(e),256);
xlabel('reconstruction error in phase (radians)');
ylabel('Number of pixels');
title('After quadtree partitioning');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));
fprintf(1,'number   of errors %#10d\n',numel(e));
fprintf(1,'fraction of errors %10.2g\n',numel(e)/nrows/ncols);

% reconstruction error in phase gradient
eg = grxpatch - grx;
min(min(eg))
max(max(eg))


nf=nf+1;h(nf)=figure;
hist(colvec(eg));
xlabel('phase phase gradient (radians per pixel)');
ylabel('Number of pixels');
title('Reconstruction errors');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
imagesc(eg);colorbar;%cmapblackzero;
title('reconstruction errors (radians per pixel)')
xlabel('column index');ylabel('row index');
printpdf(sprintf('%s_%02d.pdf',mfilename,nf));


% % Now check the C version of the reconstructor
% [ssx,srx] = unix('../src/qls2pha qsp_11176_21540_ort.qls trx');
% if ( ssx ~= 0 )
% 	txt=sprintf('FAILURE of quadphase reconstruction program\n====Reason====\n%s==============\n',srx);
%     disp(txt);
% end
% 
% rrx=zeros(nrows,ncols,'int8');
% rrx=read_pha('trx',ncols);
% nf=nf+1;h(nf)=figure;imagesc(rrx);colorbar;cmapblackzero;
% title('wrapped phase after qls2pha quad-tree reconstruction (256 DN per cycle)')
% xlabel('column index');ylabel('row index');
% 
% % Now check the C version of the reconstructor with aux_int.pha & ncol
% [ssy,sry] = unix('../src/qls2pha qsp_11176_21540_ort.qls trz psp_11176_21540_ort.pha ncol');
% if ( ssy ~= 0 )
%     disp('FAILURE of quadphase aux_input reconstruction program');
% 	txt=sprintf('FAILURE of quadphase aux_input reconstruction\n====Reason====\n%s==============\n',sry);
%     disp(txt);
% end

% rry=zeros(nrows,ncols,'int8');
% rry=read_pha('trz',ncols);
% nf=nf+1;h(nf)=figure;imagesc(rry);colorbar;cmapblackzero;
% title('wrapped phase after qls2pha aux_input reconstruction (256 DN per cycle)')
% xlabel('column index');ylabel('row index');
