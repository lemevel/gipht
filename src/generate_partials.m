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


% Jacobian matrix is partial derivative of observable with respect to parameter
dDdP = zeros(ndata,mparams);
%fprintf(1,'ID, name, min, max, mean of partial derivative, delta (radians)\n');

%ptlcol = zeros(ndata,1);

% count free parameters
ifree = find(abs(PST.ub-PST.lb) > 0.);

mfree = 0;
parfor j=ifree
    fprintf(1,'--- Calculating partial derivatives for parameter %5d %s %12.4e\n',j,char(pnames{j}),PST.scale(j));
    mfree = mfree+1;
    
    % half a step down (left) in parameter
    PSTP1 = PST;
    PSTP1.p1(j) = PSTP1.p1(j) - 1.0d0 * PST.scale(j)/2.0;
    PSTP1.flag{j} = 'F#';
    rng1 = feval(fitfun,DST,PSTP1,TST);

    % half a step up (right) in parameter
    PSTP2 = PST;
    PSTP2.p1(j) = PSTP2.p1(j) + 1.0d0 * PST.scale(j)/2.0;
    PSTP2.flag{j} = 'F#';
    rng2 = feval(fitfun,DST,PSTP2,TST);

    % partial derivative is difference (right minus left)
    der1 = colvec((rng2 - rng1) / PST.scale(j));
    %der1 = colvec((rng2 - rng1) / db / 2.0);
        
    % find valid elements
    iok1 = find(isfinite(der1)==1);  % finite value
    iok2 = find(abs(der1)>0.0);      % non-zero
    iok  = intersect(iok1,iok2);     % both of above
    
    % set all elements to zero
    ptlcol = zeros(ndata,1);
    
    % overwrite with valid elements
    ptlcol(iok) = der1(iok);
    
    % count bad elements
    nbad = ndata - numel(iok1);
    if nbad ~= 0
        warning(sprintf('nbad is %d. Replacing partial derivatives with zero.\n',nbad));
    end
    
    
    %     fprintf(1,' %03d %s %12.5E %12.5E %12.5E %12.5E\n'...
    %         ,j,char(pnames{j})...
    %         ,nanmin(ptlcol),nanmax(ptlcol),nanmean(ptlcol));
    %
    
    % return partial derivative wrt 1 parameter as 1 column in TSTP structure
    dDdP(:,j) = ptlcol;
end

fprintf(1,'\nFinished generating partial derivatives for %d free parameters\n',mfree);

% % display non-zero elements of the Hessian matrix of partial derivatives
% figure;
% spy(TSTP.partial_wrt_1param);
% xlabel(sprintf('mparam = %d columns',mparam));
% ylabel(sprintf('ndata = %d rows',ndata));

TSTP.partial_wrt_1param = sparse(dDdP);


fprintf(1,'\n\n----------------   %s ends   at %s ----------\n',upper(mfilename),datestr(now,31));

return
end


