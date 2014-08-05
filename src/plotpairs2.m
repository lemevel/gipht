function h1 = plotpairs2(tm,ts,pa,ps,imast,islav,dmod ...
, tfit,pfit,sfit ...
, xlab,ylab,titlestring,mparam)
% function h1 = plotpairs2(tm,ts,pa,ps,imast,islav,dmod ...
% , tfit,pfit,sfit ...
% , xlab,ylab,titlestring,mparam)%
% 
% inputs:
% all these are column vectors with length equal to the number of pairs
%   tm           master epochs in decimal years
%   ts           slave  epochs in decimal years
%   pa           observed volume values
%   ps           uncertainty of observed volume values
%   imast        indices to master epochs
%   islave       indices to slave epochs
% Column vectors with length equal to the number of unique epochs
%  tfit, pfit, sfit
%
% 2012-JUN-27 Kurt Feigl



error(nargchk(14,14,nargin));

ndat = numel(pa);
tmid = (tm+ts)/2;
tu = unique([tm ts]);
th = (ts-tm)/2.; % half interval
t0 = mean(tmid); % reference epoch

h1 = figure('color','w'); hold on;     

if exist('tfit','var') == 1 && exist('ppredm','var') == 1 && exist('pest','var') == 1
    fprintf(1,'Using incoming model\n');
end

%plot modeled values of Volume in black
disp('tfit');size(tfit)
if numel(tfit) > 0
    fprintf(1,'i,  tfit,  pfit\n');
    for i=1:numel(tfit)
        fprintf(1,'%3d %#10.4f %#12.4g\n',i,tfit(i),pfit(i));
    end
    plot(tfit,pfit,'k-','LineWidth',2);
end

if numel(pfit) > 0
    % plot observed values of Volume change in red
    fprintf(1,'i,tm(i),ts(i),ppredm(i),ppreds(i),ps(i)\n');
    for i=1:ndat
%         % Volume value at second epoch of each pair shows CHANGE in Volume
%         ppreds(i) = ppredm(i) + pa(i);
%         fprintf(1,'%3d %#10.4f %#10.4f %#12.4g %#12.4g %#12.4g\n',i,tm(i),ts(i),ppredm(i),ppreds(i),ps(i));
%         plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]      ,'ro-','MarkerSize',3);
%        %plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]      ,'r+-','MarkerFaceColor',[1 0 0],'MarkerSize',3);
%         plot([ts(i),ts(i)],[ppreds(i)-ps(i)      ,ppreds(i)+ps(i)],'b-'); % error bar
%         plot(ts(i)                               ,ppreds(i)-ps(i) ,'bv','MarkerFaceColor',[0 0 1],'MarkerSize',3); % bottom of error bar
%         plot(ts(i)                               ,ppreds(i)+ps(i) ,'b^','MarkerFaceColor',[0 0 1],'MarkerSize',3); % top of error bar
        % Volume value at second epoch of each pair shows CHANGE in Volume
        ppredm(i) = interp1(tfit,pfit,tm(i),'linear');
        ppreds(i) = ppredm(i) + pa(i);
        fprintf(1,'%3d %#10.4f %#10.4f %#12.4g %#12.4g %#12.4g\n',i,tm(i),ts(i),ppredm(i),ppreds(i),ps(i));
        plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]      ,'ro-','MarkerFaceColor',[1 0 0],'MarkerEdgeColor','k','MarkerSize',3,'LineWidth',1);
        %plot([tm(i),ts(i)],[ppredm(i)            ,ppreds(i)]      ,'r+-','MarkerFaceColor',[1 0 0],'MarkerSize',4,'LineWidth',1);
        plot([ts(i),ts(i)],[ppreds(i)-ps(i)      ,ppreds(i)+ps(i)],'b-','LineWidth',1); % error bar
        plot(ts(i)                               ,ppreds(i)-ps(i) ,'bv','MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1); % bottom of error bar
        plot(ts(i)                               ,ppreds(i)+ps(i) ,'b^','MarkerFaceColor',[0 0 1],'MarkerSize',1,'LineWidth',1); % top of error bar
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


