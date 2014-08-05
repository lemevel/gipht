function [rng0,TSTP] = generate_partials(DST,PST,TST)
% generate partials

fprintf(1,'\n\n----------------   %s begins at %s ----------\n',upper(mfilename),datestr(now,31));

fitfun = char(rowvec(PST.fitfun));
mparam = PST.mparam;
ndata  = numel(DST.phaobs);
pnames = PST.names;

%copy initial values to final values
PST.p1 = PST.p0;

% evaluate exact fitting function at initial estimate of parameters
rng0 = feval(fitfun,DST,PST,TST);
DST.phamod = rng0;

% count free parameters
ifree = find(abs(PST.ub-PST.lb) > 0.);
% initialize structure for speed
TSTP=struct('partial_wrt_1param',zeros(ndata,mparam));
%TSTP=struct('partial_wrt_1param',zeros(ndata,numel(ifree)));
%fprintf(1,'ID, name, min, max, mean of partial derivative, delta (radians)\n');

ptlcol = zeros(ndata,1);

mfree = 0;
for j=1:mparam
    %for j = ifree
    % separation between lower and upper bounds
    db = abs(PST.ub(j)-PST.lb(j));
    if (db > 0.)
        fprintf(1,'Calculating partial derivatives for parameter %5d %s %12.4e\n',j,char(pnames{j}),PST.scale(j));
        mfree = mfree+1;
        % half a step down (left) in parameter
        PSTP1 = PST;
        PSTP1.p1(j) = PSTP1.p1(j) - 1.0d0 * PST.scale(j)/2.0;
        PSTP1.flag{j} = 'F#';
        %PSTP1.p1(j) = PSTP1.p1(j) - 1.0d0 * db/4.0;
        %PSTP1.p1(j) = -1.0d0 * PST0.scale(j);
        %ierr=check_struct(PSTP1);
        %rng1 = feval(fitfun,DST,PSTP1,TST);
        rng1 = feval(fitfun,DST,PSTP1,TST);
        %disp rng1; size(rng1)
        %fprintf(1,'Extrema for rng1 %g %g\n',min(rng1),max(rng1));
        % half a step up (right) in parameter
        PSTP2 = PST;
        PSTP2.p1(j) = PSTP2.p1(j) + 1.0d0 * PST.scale(j)/2.0;
        PSTP2.flag{j} = 'F#';
        %PSTP2.p1(j) = PSTP2.p1(j) + 1.0d0 * db/4.0;
        %PSTP2.p1(j) = +1.0d0 * PST0.scale(j);
        %ierr=check_struct(PSTP1);
        rng2 = feval(fitfun,DST,PSTP2,TST);
        %disp rng2; size(rng2)
        %fprintf(1,'%3d %s %g %g %g %g %g %g\n',j,char(pnames{j}),nanmin(rng1),nanmax(rng1),nanmin(rng2),nanmax(rng2),nanmin(rng2-rng1),nanmax(rng2-rng1));
        % partial derivative is difference (right minus left)
        der1 = colvec((rng2 - rng1) / PST.scale(j));
        %der1 = colvec((rng2 - rng1) / db / 2.0);
        
        %disp der1; size(der1)
        
        %     scl1 = PST.scale(j)*ones(ndata,1);
        %     der2 = zeros(ndata,1);
        %     deltap = ones(ndata,1);
        nbad = numel(find(isfinite(rng1)==0)) + numel(find(isfinite(rng2)==0));
     else
        nbad = ndata;
    end
    
    
    if nbad == 0
        for i=1:ndata
            if abs(der1(i)) > 0.0
                ptlcol(i) = der1(i);
            else
                ptlcol(i) = 0.0;
            end
        end
    else
        if nbad ~= ndata
            warning(sprintf('nbad is %d. Replacing partial derivatives with zero.\n',nbad));
        end
        ptlcol = zeros(ndata,1);
    end
    
    
    %     fprintf(1,' %03d %s %12.5E %12.5E %12.5E %12.5E\n'...
    %         ,j,char(pnames{j})...
    %         ,nanmin(ptlcol),nanmax(ptlcol),nanmean(ptlcol));
    %
    
    % return partial derivative wrt 1 parameter as 1 column in TSTP structure
    % include all columns
    TSTP.partial_wrt_1param(1:ndata,j) = ptlcol;
    % include columns only for free parameters
    %TSTP.partial_wrt_1param(1:ndata,mfree) = ptlcol;
    %end
end

fprintf(1,'\nFinished generating partial derivatives for %d free parameters\n',mfree);

% % display non-zero elements of the Hessian matrix of partial derivatives
% figure;
% spy(TSTP.partial_wrt_1param);
% xlabel(sprintf('mparam = %d columns',mparam));
% ylabel(sprintf('ndata = %d rows',ndata));

TSTP.partial_wrt_1param = sparse(TSTP.partial_wrt_1param);

% make a 3-D plot
%h = plot_obscal3(DST);


fprintf(1,'\n\n----------------   %s ends   at %s ----------\n',upper(mfilename),datestr(now,31));

return
end


