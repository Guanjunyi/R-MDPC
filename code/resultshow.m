function resultshow(data,CL,centers,rho,delta,icl,rho_min,realcenter)

if nargin<3
    isCenter = 0;
else
    isCenter = 1;
end
if nargin<4
    isDecisionGraph = 0;
else
    isDecisionGraph = 1;
end
if nargin<5
    isRealcenters = 0;
else
    isRealcenters = 1;
end

PtSize = 2;

NC = length(unique(CL));
label_set = unique(CL);


[N,M] = size(data);

figure('Position',[400 400 350 300]);
cmap = UiGetColormap(NC);


for i=1:NC
    l=label_set(i);
    if M~=3
        if l~=0
            scatter(data((CL==l),1),data((CL==l),2),PtSize+5,'o','filled','MarkerFaceColor',cmap(l,:),'MarkerEdgeColor',cmap(l,:));
        else
            scatter(data((CL==l),1),data((CL==l),2),PtSize+50,'x','MarkerEdgeColor','k');
        end
    else
        if l~=0
            scatter3(data((CL==l),1),data((CL==l),2),data((CL==l),3),PtSize+5,'o','filled','MarkerFaceColor',cmap(l,:),'MarkerEdgeColor',cmap(l,:));
        else
            scatter3(data((CL==l),1),data((CL==l),2),data((CL==l),3),PtSize+5,'x','filled','MarkerEdgeColor','k');
        end
    end
    hold on
end

if isCenter
    for i=1:NC
        if M~=3
            %                 plot(data(centers(i),1),data(centers(i),2),'o','MarkerSize',PtSize+1,'MarkerFaceColor','r');
            scatter(data(centers(i),1),data(centers(i),2),PtSize+200,'p','filled','r','MarkerEdgeColor','w');
            hold on
        else
            %             scatter(data(centers(i),1),data(centers(i),2),data(centers(i),3),PtSize+80,'p','filled','r','MarkerEdgeColor','w');
            hold on
        end
    end
end

set(gca,'XTickLabel','');
set(gca,'YTickLabel','');
set(gca,'ZTickLabel','');
if M~=3
    axis off
end

if isDecisionGraph
    figure('Position',[0 400 350 300]); %% »æÖÆ¾ö²ßÍ¼
    plot(rho, delta,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k')
    hold on
    grid on;
    axis([rho_min max(rho) 0 max(delta)]);
    for i=1:NC
        if ismember(icl(i),realcenter)
        plot(rho(icl(i)), delta(icl(i)),'s','MarkerSize',10,'MarkerEdgeColor','g')
        else
        plot(rho(icl(i)), delta(icl(i)),'s','MarkerSize',10,'MarkerEdgeColor','b')
        end
    end
    
    if isRealcenters
        N_re_cts = length(unique(realcenter));
        for i=1:N_re_cts
            plot(rho(realcenter(i)), delta(realcenter(i)),'o','MarkerSize',4,'MarkerFaceColor','r','MarkerEdgeColor','r')
        end
    end
        
end


function [cmap]=UiGetColormap(NC)
%Get colormap for plotting below
 colormap hsv
cmap=colormap;
cmap=cmap(round(linspace(1,length(cmap),NC+1)),:);
cmap=cmap(1:end-1,:);

