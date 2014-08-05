function h1 = plotpairs3(tm,ts,Qdiff,Qdsig,imast,islav,Qdmod ...
, tfit,Qdfit,Qfsig ...
, xlab,ylab,titlestring,mparam)
% function h1 = plotpairs2(tm,ts,Qdiff,Qfsig,imast,islav,Qdmod ...
% , tfit,Qfit,Qfsig ...
% , xlab,ylab,titlestring,mparam)%
% 
% inputs:
% all these are column vectors with length equal to the number of pairs
%   tm           master epochs in decimal years
%   ts           slave  epochs in decimal years
%   Qdiff        observed values of differential quantity (NOT rate)
%   Qdsig        uncertainty of observed volume values
%   Qdmod        modeled values of differential quantity
%   imast        indices to master epochs
%   islave       indices to slave epochs
% Column vectors 
%  tfit time epochs
%  Qfit differential quantity as modeled by fit
%  Qfsig uncertainty on model fit
%
% 2014-06-27 Kurt Feigl



error(nargchk(14,14,nargin));

ndat = numel(Qdiff);
tmid = (tm+ts)/2;
tu = unique([tm ts]);
th = (ts-tm)/2.; % half interval
t0 = mean(tmid); % reference epoch

h1 = figure('color','w'); hold on;     

% if exist('tfit','var') == 1 && exist('ppredm','var') == 1 && exist('pest','var') == 1
%     fprintf(1,'Using incoming model\n');
% end

%plot modeled values  in black
disp('tfit');size(tfit)
if numel(Qdfit) > 0
    %fprintf(1,'i,  tfit,  Qfit\n');
%     for i=1:numel(tfit)
%         fprintf(1,'%3d %#10.4f %#12.4g\n',i,tfit(i),Qdfit(i));
%     end
    plot(tfit,Qdfit,      'k-', 'LineWidth',2); % modeled value
    plot(tfit,Qdfit-Qfsig,'k--','LineWidth',1); % lower envelope 
    plot(tfit,Qdfit+Qfsig,'k--','LineWidth',1); % upper envelope 
end


% plot observed values of differential change in red
if numel(Qdfit) > 0
    %fprintf(1,'i,tm(i),ts(i),ppredm(i),ppreds(i),Qdsig(i)\n');
    for i=1:ndat
%       plot observed differential quantity so that master epoch falls on modeled (fit) line
        ppredm(i) = interp1(tfit,Qdfit,tm(i),'linear');
        ppreds(i) = ppredm(i) +  Qdiff(i);
%        fprintf(1,'%3d %#10.4f %#10.4f %#12.4g %#12.4g %#12.4g\n',i,tm(i),ts(i),ppredm(i),ppreds(i),Qdsig(i));
        plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]         ,'ro-','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',4,'LineWidth',2); % segment between master and slave        
        plot([ts(i),ts(i)],[ppreds(i)-Qdsig(i)   ,ppreds(i)+Qdsig(i)],'b-', 'LineWidth',1);                                                                % error bar
        plot(ts(i)                               ,ppreds(i)-Qdsig(i) ,'bv', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % bottom of error bar
        plot(ts(i)                               ,ppreds(i)+Qdsig(i) ,'b^', 'MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1);                       % top of error bar
    end
end
axis([floor(min(tu)) ceil(max(tu)) -Inf +Inf]);
%axis fill;
axis auto

set(gca,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','Xcolor','k','Ycolor','k');
h=title (titlestring); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','HorizontalAlignment','Center');
h=xlabel(xlab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');
h=ylabel(ylab); set(h,'FontName','Helvetica','Fontsize',12,'FontWeight','bold','color','k');


return


