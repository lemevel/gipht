function ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, plotts,ylab)
%function ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, plotts,ylab)
%
% plot pseudo-absolute Bperp as a function of time 
% and return ktours the mininum-lenthg traveling salesman path
%
% Automatically handle 3 cases
% 1  X = t      Y = Bperp
% 2  X = t      Y = Ddop
% 3  X = Bperp  Y = Ddop
%
%
% ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, plotts)
%
%
% ktours = plotbp(tepochs, bpest, DD, trees)
%  
% ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates)
%
% ktours = plotbp(tepochs, bpest, DD, trees, iuniqorbs, uniqdates, plotts)
%     plotts = 0 do not connect points with traveling salesman net 
%     plotts = 1 connect points with traveling salesman net 
%     plotts = 2 connect points with traveling salesman net AND existence
%     plotts = 3 connect points with existence only
%
% Kurt Feigl CNRS 
% 2005 January
% 2005 JUL 14 add orbit numbers and dates as option
% 2006 FEB 23 fix mod bug
% 2007 NOV 17 connect points with traveling salesman trajectory
% 2008-MAR-29 correct annoying bug
% 2014-JUL-06 include ylab as input

fidtxtout = fopen(sprintf('%sout.txt',mfilename),'a+t');
for ifile = [1 fidtxtout]
   fprintf(ifile,'%s begins at %s\n',mfilename,datestr(now,31));
end


nargchk(4,8,nargin);
if nargin <= 6
   plotts = 0;
end
nargoutchk (1,1,nargout);

if (max(bpest)-min(bpest) < 1.0)
   if max(tepochs) < 1990
      icase = 3;
      xlab = 'Bperp (m)';
      ylab = 'Azimuthal Doppler/PRF';
      Tend = 10;
      yrbp = 100/0.1; % 100 m of Bperp is like 0.1 PRF of Doppler separatio
   else
      icase = 2;
      xlab = 'year';
      ylab = 'Azimuthal Doppler/PRF';
      Tend  = 30;
      yrbp = 1/0.001; % 1 year of time separation is like 0.001 PRF of Doppler separation
   end
else
    if exist(ylab,'var') == 1
        icase = 1;
        xlab = 'year';
        Tend = 25;
        yrbp = 1; % 1 year of time separation is like 1000 m of orbits separation
    else
        icase = 1;
        xlab = 'year';
        ylab = 'Bperp (m)';
        Tend = 25;
        yrbp = 1/1000; % 1 year of time separation is like 1000 m of orbits separation
    end
end


% define symbols to use

%            y     yellow        .     point              -     solid
%            m     magenta       o     circle             :     dotted
%            c     cyan          x     x-mark             -.    dashdot 
%            r     red           +     plus               --    dashed   
%            g     green         *     star
%            b     blue          s     square
%            w     white         d     diamond
%            k     black         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram

mysyms = {'gx-' 'ro-' 'b*-' 'ks-' 'md-' 'cv-'};
mysols = {'g-'  'r-'  'b-'  'k-'  'm-'  'c-'};
mylins = {'g:'  'r:'  'b:'  'k:'  'm:'  'c:'};
mydash = {'g--' 'r--' 'b--' 'k--' 'm--' 'c--'};
mysym0 = {'gx'  'ro'  'b*'  'ks'  'md'  'cv'};

% graphics handle to return
h=figure;hold on; 

% number of trees
[ntrees,ndum] = size(trees);

% number of pairs and number of epochs
[np,me] = size(DD); 

if plotts == 1
   titl1 = 'Optimal set of pairs';
else
   titl1 = 'Selected pairs'
end
titl = sprintf('%s np = %d me = %d ntrees = %d\n',titl1,np,me,ntrees);

%plot origin to make legend come out right
if ntrees < 10
    for j=1:ntrees
        plot(min(tepochs),0,mysyms{1+mod(j,length(mysyms))});
        tree = trees(j,:);
        k=isfinite(tree);
        tree=tree(k);
        me = length(find(k == 1));
        if nargin >= 6
            tree_name{j} = strcat(sprintf('trees %s orbits:',char(j+64)),sprintf(' %7d',iuniqorbs(tree(1:me))));
        else
            tree_name{j} = strcat(sprintf('trees %s ID:',char(j+64)),sprintf('%3d',char(j+64),tree(1:me)));
        end
    end
    ktours = zeros(size(trees));
    
    legend(tree_name,'Location','NorthOutside');
    % over plot origin with white
    plot(min(tepochs),0,'sw');
    plot(min(tepochs),0,'ow');
    plot(min(tepochs),0,'xw');
    plot(min(tepochs),0,'*w');
    plot(min(tepochs),0,'dw');
    plot(min(tepochs),0,'vw');
end


% draw end points of available pairs
for i=1:np      
   ddcol = DD(i,:); 
   j=find(abs(ddcol)>0);
   id0(i) = min(j);  
   id1(i) = max(j);
end
id2 = unique([id0 id1]);

for ifile = [1 fidtxtout]
   fprintf(ifile,'Pair trees Member0 Member1 orbn0 orbn1 year0 year1 %s %s\n',xlab,ylab);
end
i=0;

% label available epochs
for j=1:ntrees
    tree = trees(j,:);
    k=isfinite(tree);
    me = length(find(k == 1));
    tree_name{j} = sprintf('%3d',tree(1:me));
    for i=1:me
       %           plot(tree(i),tepochs(id2(tree(i))),mysyms{mod(j,length(mysy
       %           ms))}); hold on;
       px = tepochs(id2(tree(i)));
       py = bpest(id2(tree(i)));
       
       if py >  mean(bpest)
          if py >  mean(bpest) + std(bpest)
             rotang = 30;
             aline = 'left';
             dpy = 0.05 * std(bpest);
          else
             rotang = 15;
             aline = 'left';
             dpy = 0.05 * std(bpest);
          end
       else
          if py <  mean(bpest) - std(bpest)
             rotang = -30;
             aline = 'left';
             dpy = -0.05 * std(bpest);
          else
             rotang = -15;
             aline = 'left';
             dpy = -0.05 * std(bpest);
          end
       end

       % plot solid line for TSP
       %plot(px,py,mysols{1+mod(j,length(mysols))}); hold on;
       %plot(px,py,mysyms{1+mod(j,length(mysyms))}); hold on;
       if nargin >= 6
          if plotts == 0
             if mod(i,3) == 1
                hh=text (px,1.1*max(py),sprintf('%7d',iuniqorbs(id2(tree(i)))));
                set(hh,'HorizontalAlignment','left','rotation',45);
             elseif mod(i,3) == 2
                hh=text (px,1.2*max(py),sprintf('%7d',iuniqorbs(id2(tree(i)))));
                set(hh,'HorizontalAlignment','left','rotation',45);
             else
                hh=text (px,1.3*max(py),sprintf('%7d',iuniqorbs(id2(tree(i)))));
                set(hh,'HorizontalAlignment','left','rotation',45);
             end
          else
             hh=text (px,py+dpy,sprintf('%7d',iuniqorbs(id2(tree(i)))));
             set(hh,'HorizontalAlignment',aline,'rotation',rotang,'BackgroundColor',[1 1 1],'Margin',0.1);
             set(hh,'FontName','Helvetica','Fontsize',10);
          end
       end
    end
end

if plotts > 0
   for j = 1:ntrees
      % find the members of this trees
      kkeep = find(isnan(trees(j,:)) == 0);
      tspxy = zeros(numel(kkeep),2);
      % traveling salesman coordinates are time and Bperp
      tspxy(:,1) = tepochs(trees(j,kkeep));
      tspxy(:,2) = bpest(trees(j,kkeep));

      % rescale
      %yrbp = (max(tspxy(:,2))-min(tspxy(:,2)))/(max(tspxy(:,1))-min(tspxy(:,1)));
      %yrbp = 1;
      tspxy(:,1) = tspxy(:,1)*yrbp;
      nmem = numel(kkeep);
      
      %fprintf (1,'Traveling Salesman on trees %d with %d members and scale %.3f and Tend = %f\n',j,nmem,yrbp,Tend);
      
      if nmem > 3
         % traveling salesman problem
         ktour = tspsiman(tspxy,Tend);
         ktour = ktour(1:length(ktour)-1); % remove -1
      else
         if nmem == 3
            dista = sqrt( (tspxy(2,1))-(tspxy(1,1))^2 + (tspxy(2,2))-(tspxy(1,2))^2);
            distb = sqrt( (tspxy(3,1))-(tspxy(2,1))^2 + (tspxy(3,2))-(tspxy(2,2))^2);
            distc = sqrt( (tspxy(3,1))-(tspxy(1,1))^2 + (tspxy(3,2))-(tspxy(1,2))^2);
            if dista + distb < dista + distc
               ktour = [1 2 3];
            else
               ktour = [1 3 2];
            end
         else
            ktour = [1 2];
         end
      end
      ktours(j,1:numel(ktour)) = ktour;
 
      % overwrite values in ktour order
      %tree_name{j} = strcat(sprintf('trees %s orbits:',char(j+64)),sprintf(' %7d',iuniqorbs(tree(1:me))));
      %tree_name{j} = strcat(sprintf('trees %s orbits:',char(j+64)),sprintf(' %7d',iuniqorbs(trees(j,kkeep(ktour)))));
 
      for k=1:numel(ktour)-1
         i=i+1;  % count pairs
         i0=trees(j,kkeep(ktour((k))));
         i1=trees(j,kkeep(ktour((k+1))));
         % Traveling Salesman Pairs
         for ifile = [1 fidtxtout]
            fprintf(ifile,'%3d %3d %3d %3d %5d %5d %12.4f %12.4f %12.4f %12.4f\n',i,j,k,k+1 ...
               ,iuniqorbs(i0),iuniqorbs(i1),tepochs(i0),  tepochs(i1)...
               ,tepochs(i1)-tepochs(i0)...
               ,bpest(i1)-  bpest(i0));
         end
         % connect TSP with dotted line
         if plotts < 3
            plot([tepochs(i0) tepochs(i1)],[bpest(i0) bpest(i1)]...
               ,mylins{1+mod(j,length(mylins))},'Linewidth',2);
         end

%          plot([tepochs(ktour(id2(k)))  tepochs(ktour(id2(k+1)))]...
%             ,  [ bpest(ktour(id2(k)))   bpest(ktour(id2(k+1)))  ]...
%             ,    mydash{1+mod(j,length(mydash))},'Linewidth',2);        
       end

      % connect existing pairs as dashed lines
      %plot(tspxy(:,1)/yrbp,tspxy(:,2),mysols{1+mod(j,length(mysols))},'Linewidth',2); hold on;
   end
end


for i=1:np
   for j = 1:ntrees
      if sum(ismember(trees(j,:),id0(i))) == 1 & sum(ismember(trees(j,:),id1(i))) == 1
   
         %                plot([id0(i) id1(i)],[tepochs(id0(i)) tepochs(id1(i))],mysyms{1+mod(j,length(mysyms))}); hold on;
         % draw symbol
         plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mysym0{1+mod(j,length(mysym0))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
         
         if plotts == 0  || plotts == 2 || plotts == 3 % draw dashed line for possible pairs
              plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mydash{1+mod(j,length(mydash))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
              %plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mylins{1+mod(j,length(mylins))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
         end
         
%          if plotts == 2  % draw dashed line
%             plot([tepochs(id0(i)) tepochs(id1(i))],[bpest(id0(i)) bpest(id1(i))],mydash{1+mod(j,length(mydash))},'Linewidth',2,'MarkerFaceColor','k'); hold on;
%          end
      end
   end
end



plot([min(tepochs) min(tepochs)],[min(bpest)-0.1*(max(bpest)-min(bpest)) max(bpest)+0.1*(max(bpest)-min(bpest))],'w.'); % draw a white dot to stretch scales

%legend(tree_name,'Location','NorthOutside');

h2=title (titl); set(h2,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
h2=xlabel(xlab);       set(h2,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
h2=ylabel(ylab);         set(h2,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');
set(gca,'FontName','Helvetica','Fontsize',14,'FontWeight','bold');

if nargout == 0
   hold off
end

for ifile = [1 fidtxtout]
   fprintf(ifile,'%s ended at %s\n',mfilename,datestr(now,31));
end
fclose(fidtxtout);


return;

