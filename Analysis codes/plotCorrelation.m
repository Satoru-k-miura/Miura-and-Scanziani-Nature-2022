function plotCorrelation(arocX,arocY,sigX,sigY,depth,fastreg,R)

if isempty(R)
    R=true(size(sigX));
end

if isempty(depth) && isempty(fastreg)
    arocA = arocX(sigX & sigY & R);
    arocB = arocY(sigX & sigY & R);
    arocA2 = arocX(~(sigX & sigY) & R);
    arocB2 = arocY(~(sigX & sigY) & R);

    figure;hold on
    line([0 0],[-1 1])
    line([-1 1],[0 0])
    plot(arocA,arocB,'or','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
    plot(arocA2,arocB2,'or','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);

    mdl = fitlm(arocA,arocB,'linear');
    xnew = linspace(-1,1,1000)';
    [ypred,yci] = predict(mdl,xnew);
    plot(xnew,ypred,'-r');
    plot(xnew,yci(:,1),':k');
    plot(xnew,yci(:,2),':k');

    xlabel('AROC')
    ylabel('AROC')
    xlim([-1 1])
    ylim([-1 1])
    pbaspect([2 2 1])

    if length(arocA)>=2
        [R,P] = corrcoef(arocA,arocB);
        R=R(1,2);
        P=P(1,2);
        text(-0.8,0.8,['R = ',num2str(R),', p = ',num2str(P)])
    end

elseif isempty(depth) && ~isempty(fastreg)
    arocA = arocX(sigX & sigY & R);
    arocB = arocY(sigX & sigY & R);
    frSig = fastreg(sigX & sigY & R);
    arocA2 = arocX(~(sigX & sigY) & R);
    arocB2 = arocY(~(sigX & sigY) & R);
    frNonsig = fastreg(~(sigX & sigY) & R);

    figure;hold on
    line([0 0],[-1 1])
    line([-1 1],[0 0])
    for i=1:length(arocA)
        fr = frSig(i);
        if strcmp(fr,'fast')
            plot(arocA(i),arocB(i),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
        else
            plot(arocA(i),arocB(i),'^','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
        end
    end
    for i=1:length(arocA2)
        fr = frNonsig(i);
        if strcmp(fr,'fast')
            plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
        else
            plot(arocA2(i),arocB2(i),'^','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
        end
    end
    
    mdl = fitlm(arocA,arocB,'linear');
    xnew = linspace(-1,1,1000)';
    [ypred,yci] = predict(mdl,xnew);
    plot(xnew,ypred,'-r');
    plot(xnew,yci(:,1),':k');
    plot(xnew,yci(:,2),':k');

    xlabel('AROC')
    ylabel('AROC')
    xlim([-1 1])
    ylim([-1 1])
    pbaspect([2 2 1])

    if length(arocA)>2
        [R,P] = corrcoef(arocA,arocB);
        R=R(1,2);
        P=P(1,2);
        text(-0.8,0.8,['R = ',num2str(R),', p = ',num2str(P)])
    end
    
elseif ~isempty(depth) && isempty(fastreg)
    arocA = arocX(sigX & sigY & R);
    arocB = arocY(sigX & sigY & R);
    dSig = depth(sigX & sigY & R);
    arocA2 = arocX(~(sigX & sigY) & R);
    arocB2 = arocY(~(sigX & sigY) & R);
    dNonsig = depth(~(sigX & sigY) & R);
    
    figure;hold on
    line([0 0],[-1 1])
    line([-1 1],[0 0])
    for i=1:length(arocA)
        d = dSig(i);
        if d<-175 %L2/3
            plot(arocA(i),arocB(i),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);
        elseif d<0 && d>=-175 %L4
            plot(arocA(i),arocB(i),'o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',8);
        elseif d>=0 && d<=200 % L5
            plot(arocA(i),arocB(i),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
        else
            plot(arocA(i),arocB(i),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8);
        end
    end
    for i=1:length(arocA2)
        d = dNonsig(i);
        if d<-175 %L2/3
            plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','r','MarkerSize',8);
        elseif d<0 && d>=-175 %L4
            plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','g','MarkerSize',8);
        elseif d>=0 && d<=200 % L5
            plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
        else
            plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',8);
        end
    end
    
    mdl = fitlm(arocA,arocB,'linear');
    xnew = linspace(-1,1,1000)';
    [ypred,yci] = predict(mdl,xnew);
    plot(xnew,ypred,'-r');
    plot(xnew,yci(:,1),':k');
    plot(xnew,yci(:,2),':k');

    xlabel('AROC')
    ylabel('AROC')
    xlim([-1 1])
    ylim([-1 1])
    pbaspect([2 2 1])

    if length(arocA)>2
        [R,P] = corrcoef(arocA,arocB);
        R=R(1,2);
        P=P(1,2);
        text(-0.8,0.8,['R = ',num2str(R),', p = ',num2str(P)])
    end
    
else
    arocA = arocX(sigX & sigY & R);
    arocB = arocY(sigX & sigY & R);
    dSig = depth(sigX & sigY & R);
    frSig = fastreg(sigX & sigY & R);
    arocA2 = arocX(~(sigX & sigY) & R);
    arocB2 = arocY(~(sigX & sigY) & R);
    dNonsig = depth(~(sigX & sigY) & R);
    frNonsig = fastreg(~(sigX & sigY) & R);
    
    figure;hold on
    line([0 0],[-1 1])
    line([-1 1],[0 0])
    for i=1:length(arocA)
        d = dSig(i);
        fr = frSig(i);
        if strcmp(fr,'fast')
            if d<-175 %L2/3
                plot(arocA(i),arocB(i),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);
            elseif d<0 && d>=-175 %L4
                plot(arocA(i),arocB(i),'o','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',8);
            elseif d>=0 && d<=200 % L5
                plot(arocA(i),arocB(i),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
            else
                plot(arocA(i),arocB(i),'o','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8);
            end
        else
            if d<-175 %L2/3
                plot(arocA(i),arocB(i),'^','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);
            elseif d<0 && d>=-175 %L4
                plot(arocA(i),arocB(i),'^','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',8);
            elseif d>=0 && d<=200 % L5
                plot(arocA(i),arocB(i),'^','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
            else
                plot(arocA(i),arocB(i),'^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',8);
            end
        end
    end
    for i=1:length(arocA2)
        d = dNonsig(i);
        fr = frNonsig(i);
        if strcmp(fr,'fast')
            if d<-175 %L2/3
                plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','r','MarkerSize',8);
            elseif d<0 && d>=-175 %L4
                plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','g','MarkerSize',8);
            elseif d>=0 && d<=200 % L5
                plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
            else
                plot(arocA2(i),arocB2(i),'o','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',8);
            end
        else
            if d<-175 %L2/3
                plot(arocA2(i),arocB2(i),'^','MarkerFaceColor','none','MarkerEdgeColor','r','MarkerSize',8);
            elseif d<0 && d>=-175 %L4
                plot(arocA2(i),arocB2(i),'^','MarkerFaceColor','none','MarkerEdgeColor','g','MarkerSize',8);
            elseif d>=0 && d<=200 % L5
                plot(arocA2(i),arocB2(i),'^','MarkerFaceColor','none','MarkerEdgeColor','b','MarkerSize',8);
            else
                plot(arocA2(i),arocB2(i),'^','MarkerFaceColor','none','MarkerEdgeColor','k','MarkerSize',8);
            end
        end
    end
    
    mdl = fitlm(arocA,arocB,'linear');
    xnew = linspace(-1,1,1000)';
    [ypred,yci] = predict(mdl,xnew);
    plot(xnew,ypred,'-r');
    plot(xnew,yci(:,1),':k');
    plot(xnew,yci(:,2),':k');

    xlabel('AROC')
    ylabel('AROC')
    xlim([-1 1])
    ylim([-1 1])
    pbaspect([2 2 1])

    if length(arocA)>2
        [R,P] = corrcoef(arocA,arocB);
        R=R(1,2);
        P=P(1,2);
        text(-0.8,0.8,['R = ',num2str(R),', p = ',num2str(P)])
    end
end