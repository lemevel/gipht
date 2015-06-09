%test_quadphase

% demostrate quadtree compression on wrapped phase
% lap11jul09 removed coherence references
% lap11jul09 added code to skip over signature NX & NY in .qls file
% lap14jul09 added code to test quadesahp reconstructor
clear all; close all;
giphtpath

% number of columns in each interferogram
 ncols =         1420
% number of lines in each interferogram 
 nrows =         1230

cmd1='../src/pha2qls3 psp_11176_21540_ort.pha  1420 1230 qsp_11176_21540_ort.pha grx_11176_21540_ort.pha gry_11176_21540_ort.pha qsp_11176_21540_ort.qls'
% usage: pha2qls     input.pha                nx   ny   output.pha              output.qls              ithresh[8] minpix[16]
                     
%unix(sprintf('%s',cmd1))  % default values
%unix(sprintf('%s 25',cmd1)) % set max circular mean deviation to 25/256 of a cycle
unix(sprintf('%s  8 9',cmd1)) % set min number of OK pixels in patch to 9
%unix(sprintf('%s  16 25',cmd1)) % set min number of OK pixels in patch to 25
%unix(sprintf('%s psp_11176_21540_ort.qls 8 16 ../demoQ/coh_11176_21540_ort.oct 1',cmd1)) % use only pixels with coh > 1st percentile

p=read_pha('psp_11176_21540_ort.pha',ncols);
q=read_pha('qsp_11176_21540_ort.pha',ncols);
grx=read_pha('grx_11176_21540_ort.pha',ncols);
%lap c=read_oct('coh_11176_21540_ort.oct',ncols);

%lap c=read_oct('coh_11176_21540_ort.oct',ncols);

figure;imagesc(p);colorbar;cmapblackzero;
title('wrapped phase (256 DN per cycle)');
xlabel('column index');ylabel('row index');
%printjpg('psp_11176_21540_ort.jpg');
figure;imagesc(q);colorbar;%cmapblackzero;
title('wrapped phase after quad-tree partitioning (256 DN per cycle)')
xlabel('column index');ylabel('row index');
figure;imagesc(10*grx);colorbar;%cmapblackzero;
title('x gradient after quad-tree partitioning (256 DN per cycle)')
xlabel('column index');ylabel('row index');

%printjpg('qsp_11176_21540_ort.jpg');
%lap figure;imagesc(c);colorbar;cmapblackzero;
%lap colormap('gray');
%lap title('coherence (0 to 255)')
%lap xlabel('column index');ylabel('row index');
%printjpg('coh_11176_21540_ort.jpg');

%lap iok=find(c>0);c=c(iok);medcoh=median(c)
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
qp=ilist(2:end,4);  % ohase value coded -128 to +127 such that 256 DN = 1 cycle
% gx=ilist(2:end,5);  % X-ward gradient of phase value coded -128 to +127 such that 256 DN = 1 cycle
% gy=ilist(2:end,6);  % Y-ward gradient of phase value coded -128 to +127 such that 256 DN = 1 cycle
gx=double(ilist(2:end,5))/1000;  % X-ward gradient of phase value coded -128 to +127 such that 256 DN = 1 cycle
gy=double(ilist(2:end,6))/1000;  % Y-ward gradient of phase value coded -128 to +127 such that 256 DN = 1 cycle
i2=i1+kw; % patches are square
j2=j1+kw;
%iok = find(abs(qp) > 0);
%im=im(iok);
%jm=jm(iok);
%qp=qp(iok);
%i1=i1(iok);
%i2=i2(iok);
%j1=j1(iok);
%j2=j2(iok);
min(qp)
max(qp)
mean(qp)
std(single(qp))
npatch = numel(qp)
figure
hist(double(colvec(qp)),64);axis tight;
xlabel('phase value (256 DN per cycle)');
ylabel('Number of phase values');
title('After quadtree partitioning');

figure
hist(double(colvec(gx)),64);axis tight;
xlabel('x gradient value (256 DN per cycle)');
ylabel('Number of phase values');
title('X-gradient');
figure;qqplotvonmises(double(colvec(gx))/256.,'X-gradient');
figure;qqplot(double(colvec(gx))/256.);title('X-gradient');

figure
hist(double(colvec(gy)),64);axis tight;
xlabel('y gradient value (256 DN per cycle)');
ylabel('Number of phase values');
title('y-gradient');
figure
qqplotvonmises(double(colvec(gy))/256.0,'Y-gradient');
figure;qqplot(double(colvec(gy))/256.);title('Y-gradient');

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

r=zeros(nrows,ncols,'int8');
disp('dimensions of r');size(r)
disp('dimensions of q');size(q)
%r=zeros(nrows,ncols);
for k=1:npatch
    for i=i1(k):i2(k)
        for j=j1(k):j2(k)
            i0 = double(i1(k)+i2(k))/2.0;
            j0 = double(j1(k)+j2(k))/2.0;
            r(j,i)=int8(qp(k) + (i-i0)*gx(k) + (j-j0)*gy(k));
            %r(j,i)=int8(qp(k)                 + (j-j0)*gy(k));
        end
    end
end
disp('dimensions of r');size(r)
disp('dimensions of q');size(q)
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

figure;imagesc(r);colorbar;cmapblackzero;
title('wrapped phase after reconstructing quad-tree partitioning (256 DN per cycle)')
xlabel('column index');ylabel('row index');
%printjpg('rsp_11176_21540_ort.jpg');

figure
hist((reshape(double(r-q),numel(r),1)));
xlabel('phase value (256 DN per cycle)');
ylabel('Number of pixels');
title('Reconstruction errors');

e=r-q;
figure;imagesc(e);colorbar;cmapblackzero;
title('reconstruction errors (256 DN per cycle)')
xlabel('column index');ylabel('row index');

figure
iok=find(abs(double(e))>1e-6);e=e(iok);
hist(double(reshape(e,numel(e),1)),256);
xlabel('reconstruction error (256 DN per cycle)');
ylabel('Number of pixels');
title('After quadtree partitioning');
fprintf(1,'number   of errors %#10d\n',numel(e));
fprintf(1,'fraction of errors %10.2g\n',numel(e)/nrows/ncols);

% % Now check the C version of the reconstructor
% [ssx,srx] = unix('../src/qls2pha qsp_11176_21540_ort.qls trx');
% if ( ssx ~= 0 )
% 	txt=sprintf('FAILURE of quadphase reconstruction program\n====Reason====\n%s==============\n',srx);
%     disp(txt);
% end
% 
% rrx=zeros(nrows,ncols,'int8');
% rrx=read_pha('trx',ncols);
% figure;imagesc(rrx);colorbar;cmapblackzero;
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
% figure;imagesc(rry);colorbar;cmapblackzero;
% title('wrapped phase after qls2pha aux_input reconstruction (256 DN per cycle)')
% xlabel('column index');ylabel('row index');
