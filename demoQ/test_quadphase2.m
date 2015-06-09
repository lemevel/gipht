%test_quadphase

% demostrate quadtree compression on wrapped phase
clear all; close all;
% number of columns in each interferogram
 ncols =         1420
% number of lines in each interferogram 
 nrows =         1230

cmd1='../src/quadphase psp_11176_21540_ort.pha  1420 1230 qsp_11176_21540_ort.pha qsp_11176_21540_ort.phalst'
unix(sprintf('%s',cmd1))  % default values
%unix(sprintf('%s 25',cmd1)) % set max circular mean deviation to 25/256 of a cycle
%unix(sprintf('%s  8 20',cmd1)) % set min number of OK pixels in patch to 20
%unix(sprintf('%s psp_11176_21540_ort.phalst 8 16 ../demoQ/coh_11176_21540_ort.oct 1',cmd)) % use only pixels with coh > 1st percentile

p=read_pha('psp_11176_21540_ort.pha',ncols);
q=read_pha('qsp_11176_21540_ort.pha',ncols);
c=read_oct('coh_11176_21540_ort.oct',ncols);

figure;imagesc(p);colorbar;cmapblackzero;
title('wrapped phase (256 DN per cycle)');
xlabel('column index');ylabel('row index');
%printjpg('psp_11176_21540_ort.jpg');
figure;imagesc(q);colorbar;cmapblackzero;
title('wrapped phase after quad-tree partitioning (256 DN per cycle)')
xlabel('column index');ylabel('row index');
%printjpg('qsp_11176_21540_ort.jpg');
figure;imagesc(c);colorbar;cmapblackzero;
colormap('gray');
title('coherence (0 to 255)')
xlabel('column index');ylabel('row index');
%printjpg('coh_11176_21540_ort.jpg');

iok=find(c>0);c=c(iok);medcoh=median(c)
figure
hist(double(reshape(c,numel(c),1)),256);
xlabel('coherence low (0) to high (255)');
ylabel('Number of pixels');
title('Histogram');


%[im,jm,qp,nok,nnull,i1,i2,j1,j2,tr]=textread('qsp_11176_21540_ort.phalst','%d%d%d%d%d%d%d%d%d%d','headerlines',1);
ilist = read_i2('qsp_11176_21540_ort.phalst',10);
whos ilist
im=ilist(:,1);
jm=ilist(:,2);
qp=ilist(:,3);
i1=ilist(:,4);
i2=ilist(:,5);
j1=ilist(:,6);
j2=ilist(:,7);
tr=ilist(:,8);
nok  =ilist(:,9);
nnull=ilist(:,10);
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
hist(double(reshape(qp,numel(qp),1)),64);axis tight;
xlabel('phase value (256 DN per cycle)');
ylabel('Number of phase values');
title('After quadtree partitioning');

% convert indices from C to Fortran convention
i1 = i1+1;
i2 = i2+1;
j1 = j1+1;
j2 = j2+1;
ibad=find(i1<1)    ;i1(ibad)=1;
ibad=find(i1>ncols);i1(ibad)=ncols;
ibad=find(j1<1)    ;j1(ibad)=1;
ibad=find(j1>nrows);j1(ibad)=nrows;
r=zeros(nrows,ncols,'int8');
%r=zeros(nrows,ncols);
for k=1:npatch
    r(j1(k):j2(k),i1(k):i2(k))=int8(qp(k));
%   for j=j1(k)+1:j2(k)+1
%      for i=i1(k)+1:i2(k)+1
%         fprintf(1,'%#4d %#4d %#4d %#4d\n',i,j,k,qp(k));
%         r(j,i)=qp(k);
%      end
%   end
end


figure;imagesc(r);colorbar;cmapblackzero;
title('wrapped phase after reconstructing quad-tree partitioning (256 DN per cycle)')
xlabel('column index');ylabel('row index');
%printjpg('rsp_11176_21540_ort.jpg');

figure
hist(double(reshape(r-q,numel(r),1)));
xlabel('phase value (256 DN per cycle)');
ylabel('Number of pixels');
title('Reconstruction errors');

e=r-q;
figure;imagesc(e);colorbar;cmapblackzero;
title('reconstruction errors (256 DN per cycle)')
xlabel('column index');ylabel('row index');
%printjpg('esp_11176_21540_ort.jpg');

figure
iok=find(r~=0);r=r(iok);
hist(double(reshape(r,numel(r),1)),256);
xlabel('phase value (256 DN per cycle)');
ylabel('Number of pixels');
title('After quadtree partitioning');

   


