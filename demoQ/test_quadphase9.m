%test_quadphase

% demostrate quadtree compression on wrapped phase
% lap11jul09 removed coherence references
% lap11jul09 added code to skip over signature NX & NY in .qls file
% lap14jul09 added code to test quadesahp reconstructor
close all;
clear all;
giphtpath
printfun = 'printpdf'
format compact
nf = 0;

% number of columns in each interferogram
ncols =         1420
% number of lines in each interferogram
nrows =         1230

% run the resampling - good for estimating gradients on a few patches
%cmd1='../src/pha2qls.maci64 psp_11176_21540_ort.pha  1420 1230 -P qsp_11176_21540_ort.pha -G -L 127 -N 9 -M 32 >&! pha2qls.out'
% run the resampling - good for image compression
%cmd1='../src/pha2qls.maci64 psp_11176_21540_ort.pha  1420 1230 -P qsp_11176_21540_ort.pha -G -L 4 -N 9 -M 127 >&! pha2qls.out'
% run the resampling - compromise
%cmd1='../src/pha2qls.maci64 psp_11176_21540_ort.pha  1420 1230 -P qsp_11176_21540_ort.pha -G -L 32 -N 9 -M 16 >&! pha2qls.out'
%cmd1='../src/pha2qls.a64 psp_11176_21540_ort.pha  1420 1230 -P qsp_11176_21540_ort.pha -X grx_11176_21540_ort.i2 -Y gry_11176_21540_ort.i2 -L 32 -N 9 -M 16 >&! pha2qls.out'
%cmd1='../src/pha2qls.a64 pha_11176_21540_ort.pha  1420 1230 -P qha_11176_21540_ort.pha -X grx_11176_21540_ort.i2 -Y gry_11176_21540_ort.i2 -L 32 -N 16 -M 16 >&! pha2qls.out'
 cmd1='../src/pha2qls.maci64 psp_11176_21540_ort.pha  1420 1230 -P qha_11176_21540_ort.pha -X grx_11176_21540_ort.i2 -Y gry_11176_21540_ort.i2 -L 32 -N 9 -M 16 >&! pha2qls.out'

tstart=tic;
%unix(sprintf('%s',cmd1))  % default values
%unix(sprintf('%s 25',cmd1)) % set max circular mean deviation to 25/256 of a cycle
unix(sprintf('%s  8 9',cmd1)) % set min number of OK pixels in patch to 9
%unix(sprintf('%s  16 25',cmd1)) % set min number of OK pixels in patch to 25
%unix(sprintf('%s psp_11176_21540_ort.qls 8 16 ../demoQ/coh_11176_21540_ort.oct 1',cmd1)) % use only pixels with coh > 1st percentile
elapsed_time_in_seconds = toc(tstart)

p  =2.0*pi*double(read_pha('psp_11176_21540_ort.pha',ncols))/256.;
q  =2.0*pi*double(read_pha('qha_11176_21540_ort.pha',ncols))/256.;
gry=2.0*pi*double(read_i2 ('gry_11176_21540_ort.i2', ncols))/256./256.;
grx=2.0*pi*double(read_i2 ('grx_11176_21540_ort.i2', ncols))/256./256.;

nf=nf+1;h(nf)=figure;imagesc(p);colorbar;cmapblackzero;
title('wrapped phase (radians)');
xlabel('column index');ylabel('row index');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));
%printjpg('psp_11176_21540_ort.jpg');
nf=nf+1;h(nf)=figure;imagesc(q);colorbar;cmapblackzero;
title('wrapped phase after quad-tree partitioning (radians per pixel)')
xlabel('column index');ylabel('row index');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;imagesc(grx);colorbar;cmapblackzero;
title('row-ward gradient after quad-tree partitioning (radians per pixel)')
xlabel('column index');ylabel('row index');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;imagesc(gry);colorbar;cmapblackzero;
title('col-ward gradient after quad-tree partitioning (radians per pixel)')
xlabel('column index');ylabel('row index');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));


%[im,jm,qpv,nok,nnull,i1,i2,j1,j2,tr]=textread('qha_11176_21540_ort.qls','%d%d%d%d%d%d%d%d%d%d','headerlines',1);
qlist = read_i2('psp_11176_21540_ort.qls',6);


% write a program in C to expand .phalist into .pha file, as done in the
% following lines
whos qlist
% nx=qlist(1,3); % Number of Cols
% ny=qlist(1,4); % Number of Rows
% if ( ncols ~= nx || nrows ~= ny )
%     disp('Warning nx & nx dont match qls file header');
% end
if typecast(qlist(1,3:4),'int32') ~= ncols  || typecast(qlist(1,5:6),'int32') ~= nrows
    error(sprintf('Number of columns (%d %d) or rows (%d %d) incorrect.\n'...
        ,typecast(qlist(1,3:4),'int32'),ncols...
        ,typecast(qlist(1,5:6),'int32'),nrows));
end


i1=double(qlist(2:end,1));  % Index to col of first pixel in patch
j1=double(qlist(2:end,2));  % Index to row of first pixel in patch
kw=double(qlist(2:end,3));  % width (and height) of square patch
% qpv=double(qlist(2:end,4));  % ohase value coded -128 to +127 such that 256 DN = 1 cycle
% gx=double(qlist(2:end,5));  % X-ward gradient of phase value coded -128 to +127 such that 256^2 DN = 1 cycle
% gy=double(qlist(2:end,6));  % Y-ward gradient of phase value coded -128 to +127 such that 256^2 DN = 1 cycle
% % 2011-MAR-24 - GREAT BIG BUG - NOW fixed in pha2qls3.c
% qqp=2*pi*double(qlist(2:end,4))/256/256;% phase value
% qgx=2*pi*double(qlist(2:end,5))/256/256;  % X-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
% qgy=2*pi*double(qlist(2:end,6))/256/256;  % Y-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
qpv=2*pi*double(qlist(2:end,4))/256./256.;  % phase value
grx=2*pi*double(qlist(2:end,5))/256./256.;  % X-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel
gry=2*pi*double(qlist(2:end,6))/256./256.;  % Y-ward gradient of phase value coded as 256^2 DN = 1 cycle per pixel

%printjpg('rsp_11176_21540_ort.jpg');

% i2=i1+kw; % patches are square
% j2=j1+kw;
i2=i1+kw-1; % patches are square
j2=j1+kw-1;
%ibad = find(abs(qpv) > 0);
%im=im(ibad);
%jm=jm(ibad);
%qpv=qpv(ibad);
%i1=i1(ibad);
%i2=i2(ibad);
%j1=j1(ibad);
%j2=j2(ibad);
min(qpv)
max(qpv)
mean(qpv)
std(qpv)
npatch = numel(qpv)

figure;
nf=nf+1;h(nf)=figure;
hist(double(qpv(:)),64);axis tight;
xlabel('phase value (radians)');
ylabel('Number of phase values');
title('After quadtree partitioning');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
hist(double(grx(:)),64);axis tight;
xlabel('x gradient value (radians per pixel)');
ylabel('Number of phase values');
title('X-gradient');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));

nf=nf+1;h(nf)=figure;
hist(double(gry(:)),64);axis tight;
xlabel('y gradient value (radians per pixel)');
ylabel('Number of phase values');
title('Y-gradient');
feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));

% convert indices from C to Fortran convention
i1 = i1+1;
i2 = i2+1;
j1 = j1+1;
j2 = j2+1;
imid = double(i1+i2)/2.0;
jmid = double(j1+j2)/2.0;

% kbad=find(i1<1)    ;i1(kbad)=1;
% kbad=find(i1>ncols);i1(kbad)=ncols;
% kbad=find(j1<1)    ;j1(kbad)=1;
% kbad=find(j1>nrows);j1(kbad)=nrows;
fprintf(1,'Extrema of i1 %d %d\n',min(i1),max(i1));
fprintf(1,'Extrema of j1 %d %d\n',min(j1),max(j1));

% kbad=find(i2<1)    ;i2(kbad)=1;
% kbad=find(i2>ncols);i2(kbad)=ncols;
% kbad=find(j2<1)    ;j2(kbad)=1;
% kbad=find(j2>nrows);j2(kbad)=nrows;
fprintf(1,'Extrema of i2 %d %d\n',min(i2),max(i2));
fprintf(1,'Extrema of j2 %d %d\n',min(j2),max(j2));

% Test reconstruction two ways, then compare
for kk = 1:3    
    if kk == 1
        r1=zeros(nrows,ncols);
        grxpatch = zeros(nrows,ncols);
        %%% Recoded to avoid triple loop - goes 10x faster PES
        for k=1:npatch;
            % i2buf[0]=(short) i1;		        	             /* Column Index (X coordinate) of Upper Left Corner   */
            % i2buf[1]=(short) j1;                               /* Row    Index (Y coordinate) of Upper Left Corner   */
            % i2buf[2]=(short) i2-i1+1;		                     /* Number of pixels, width and height of square patch */
            % i2buf[3]=(short) rint(256.0f*(double)cmdr);        /* phase value in 256^2 DN */
            % i2buf[4]=(short) rint(256.0f*slopevector.slopex);  /* value of phase gradient in +X direction in 256^2 DN == 1 cycle/pixel */
            % i2buf[5]=(short) rint(256.0f*slopevector.slopey);  /* value of phase gradient in -Y direction in 256^2 DN == 1 cycle/pixel */
            
            
            %             jvec=double(j1(k):j2(k)); % row index - over Y
            %             ivec=double(i1(k):i2(k)); % column index - over X
            %             r(jvec,ivec)=qpv(k); % add in DC offset
            %   distance from center of patch
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),gry(k)*(jvec-(jvec(1)+jvec(end))/2.0));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),grx(k)*(ivec-(ivec(1)+ivec(end))/2.0));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),gry(k)*(jvec-double(j1(k)+j2(k))/2.0));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),grx(k)*(ivec-double(i1(k)+i2(k))/2.0));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),gry(k)*(jvec-jmid(k)));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),grx(k)*(ivec-imid(k)));
            %dr = grx(k)*(ivec-double(i1(k)+i2(k))/2.0) + gry(k)*(jvec-double(j1(k)+j2(k))/2.0);
            %             dr = grx(k)*(ivec-imid(k)) + gry(k)*(jvec-jmid(k));
            %             r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),dr);
            %            r(jvec,ivec)=bsxfun(@plus,r(jvec,ivec),grx(k)*(ivec-imid(k))+gry(k)*(jvec-jmid(k)));
            %grxpatch(jvec,ivec)=grx(k); % gradient
            for i=i1(k):i2(k)               
                for j=j1(k):j2(k)
                    r1(j,i) = qpv(k) + (i-imid(k))*grx(k) + (j-jmid(k))*gry(k);
                    grxpatch(j,i) = grx(k);
                end
            end
            %if k < 5
            if abs(grx(k)) > 0.05 && abs(gry(k)) > 0.05
%                 disp('jvec'); jvec
%                 disp('ivec'); ivec
%                 disp('r');    r(jvec,ivec)
                 fprintf(1,'%3d %3d %#10.4f %#10.4f %#10.4f\n',jmid(k),imid(k),qpv(k),grx(k),gry(k));
%                  disp('jvec'); j1(k):j2(k)
%                  disp('ivec'); i1(k):i2(k)
%                disp('r1');    
                 r1(j1(k):j2(k),i1(k):i2(k));
            end
        end
        r1 = rwrapm(r1);
        r  = r1;
    elseif kk == 2
        % Now check the C version of the reconstructor
       %[ssx,srx] = unix('../src/qls2pha.maci64 psp_11176_21540_ort.qls -o trx.pha -d 1');
       %[ssx,srx] = unix('../src/qls2pha.a64 psp_11176_21540_ort.qls -o rsp_11176_21540_ort.pha');
       [ssx,srx] = unix('../src/qls2pha.maci64 psp_11176_21540_ort.qls -o rsp_11176_21540_ort.pha');
        if ( ssx ~= 0 )
            txt=sprintf('FAILURE of quadphase reconstruction program\n====Reason====\n%s==============\n',srx);
            disp(txt);
        end
        
        r2=2.0*pi*double(read_pha('rsp_11176_21540_ort.pha',ncols))/256.0;
        r =r2;
    end
        
    disp('dimensions of r');size(r)
    disp('dimensions of q'); size(q)
    
    if kk < 3
        nf=nf+1;h(nf)=figure;
        imagesc(r);colorbar;colormap('jet');cmapblackzero;
        switch kk
            case 1
                title('wrapped phase after reconstructing quad-tree partitioning from QLS list (radians) by MATLAB');
            case 2
                title('wrapped phase after reconstructing quad-tree partitioning from QLS list (radians) by QLS2PHA');
        end
        xlabel('column index');ylabel('row index');
        feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));      
    end
    
    
    if kk == 1
        nf=nf+1;h(nf)=figure;
        imagesc(grxpatch);colormap('jet');colorbar;cmapblackzero;
        title('X phase gradient from QLS list (radians per pixel)')
        xlabel('column index');ylabel('row index');
        feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));
    end
    
    % reconstruction deviations in phase
    if kk == 1 || kk == 2
        %e=r-p;
        e=rarcm(r,p);
        %e=rwrapm(r-p);
        iempty=find(abs(r)<2.0*pi/256.0);
        e(iempty)=NaN;
        iempty=find(abs(p)<2.0*pi/256.0);
        e(iempty)=NaN;
    elseif kk == 3
        e = rarcm(r1,r2);
        %e=rwrapm(r1-r2);
        iempty=find(abs(r1)<2.0*pi/256.0);
        e(iempty)=NaN;
        iempty=find(abs(r2)<2.0*pi/256.0);
        e(iempty)=NaN;
    end
    
    nf=nf+1;h(nf)=figure;
    hist(colvec(e),256);
    xlabel('phase value (radians)');
    ylabel('Number of pixels');
    if kk == 1
        title('deviations after reconstruction by MATLAB  (radians)');
    elseif kk == 2
        title('deviations after reconstruction by QLS2PHA (radians)');
    else
        title('deviations between two reconstructions (radians)');
    end
    feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));
    
    nf=nf+1;h(nf)=figure;
    imagesc(e);colormap('jet');colorbar;cmapblackzero;
    if kk == 1
        title('nonzero deviations after reconstruction by MATLAB  (radians)');
    elseif kk == 2
        title('nonzero deviations after reconstruction by QLS2PHA (radians)');
    else
        title('nonzero deviations between two reconstructions (radians)');
    end
    xlabel('column index');ylabel('row index');
    feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));
    
    if kk == 1 || kk == 2
        nf=nf+1;h(nf)=figure;
        ibad=find(abs(e)>2.0*pi/256.0);
        e=e(ibad);
        ibad2=find(abs(e)>4.0*pi/256.0); % above threshold
        hist(colvec(e),256);
        xlabel('reconstruction deviations in phase (radians)');
        ylabel('Number of pixels');
        if kk == 1
            title('nonzero deviations after reconstruction by MATLAB  (radians)');
        elseif kk == 2
            title('nonzero deviations after reconstruction by QLS2PHA (radians)');
        end
        feval(printfun,sprintf('%s_%02d.pdf',mfilename,nf));
    end
    
    fprintf(1,'number   of deviations                 %#10d\n' ,numel(e));
    fprintf(1,'number   of deviations above threshold %#10d\n' ,numel(ibad2));
    fprintf(1,'fraction of deviations                 %12.4f\n',numel(e)/nrows/ncols);
    fprintf(1,'fraction of deviations above threshold %12.4f\n',numel(ibad2)/nrows/ncols);
    
end


