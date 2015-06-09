function cmap=cmapblackzero(ival)
% make a color map with zero set to black
% updated to work with any clim PES
if nargin > 0
    tmap= colormap(jet(256));
    [nlevels,ndum] = size(tmap);
    cmap=zeros(nlevels+1,3);
    for i=1:nlevels
        cmap(i+1,:)=tmap(i,:);
    end
    cmap(1,:) = [0 0 0]; % black
else
    cmap=colormap(jet(256));
    clim=sort(get(gca,'clim'));
    vec=abs(linspace(clim(1),clim(2),length(cmap)));
    
    [smallest zerolevel]=min(vec);
    zerolevel=find(vec==smallest);  % this captures case when their are two bins equidistant from zer0
   
    cmap(zerolevel,:)=zeros(size(zerolevel,2),3);
end
colormap(gca,cmap);

