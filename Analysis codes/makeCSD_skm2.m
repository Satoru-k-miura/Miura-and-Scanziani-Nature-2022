function makeCSD_skm2(avLFPbySweep)

numch = length(avLFPbySweep);
chsp = 25;
timescale = [-2000 2000];
depth = [chsp*2 (numch-1)*chsp-chsp*2] - 14*25;
ylen = depth(2)-depth(1);
ystep = ylen/(numch-4); % 28 channels to plot
yinit = ystep/2;
Fs = 1000;

% Display average LFP
figure();
subplot(2,1,1);
cumOffset=0;
for i=1:length(avLFPbySweep)
    o=max(avLFPbySweep{i});
    cumOffset=cumOffset+o;
    plot(0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs),avLFPbySweep{i}-cumOffset);
    hold on;
end

% Calculate CSD
csd=zeros(length(avLFPbySweep),size(avLFPbySweep{1},2));
csd(1,:)=(2*avLFPbySweep{1}+avLFPbySweep{2})/3;
csd(end,:)=(2*avLFPbySweep{end}+avLFPbySweep{end-1})/3;
for i=2:size(csd,1)-1
    csd(i,:)=(avLFPbySweep{i-1}+2*avLFPbySweep{i}+avLFPbySweep{i+1})/4;
end
for i=1:size(csd,1)
    csd(i,:)=smooth(csd(i,:),17);
end
data1=csd(1:2:size(csd,1),:);
data2=csd(2:2:size(csd,1),:);
% data3 = csd;
% 
% data3=diff(data3,2,1)/0.050^2/1000;
data1=diff(data1,2,1)/0.050^2/1000;
data2=diff(data2,2,1)/0.050^2/1000; % mV/mm^2

data=zeros(length(data1),2*size(data1,1));
DataCounter=1;
for i=1:size(data1,1)
    data(:,DataCounter)=data1(i,:);
    DataCounter=DataCounter+1;
    data(:,DataCounter)=data2(i,:);
    DataCounter=DataCounter+1;
end
data=-data';
% data=-data3;


% 2D linear interpolation
% [x y]=size(data);
% y=1:y;
% x=1:x;
% [xi yi]=meshgrid(1:1:max(x),1:0.05:max(y));
% dataInt=interp2(x,y,data',xi,yi);
[x,y]=size(data);
y=1:y;
x=1:x;
[xi,yi]=meshgrid(1:0.25:max(x),1:0.25:max(y));
dataInt=interp2(x,y,data',xi,yi);

subplot(2,1,2);
dataInt=dataInt';
imagesc(timescale,depth,dataInt);
drawnow;
grid on;
ylim(depth);xlim([-500 500])
[map] = diverging_map(0:0.01:1,[0 0 1],[1 0 0]);
colormap(map);

figure;
yyaxis left
imagesc(timescale,depth,dataInt);
drawnow;
colormap(map);
ylim(depth);xlim([-500 500])
% ylim([-25 825])
hold on
yyaxis right
cumOffset=chsp;
scl = 1;
avLFPbySweep = flipud(avLFPbySweep);
for i=3:length(avLFPbySweep)-2
    o=chsp;
    cumOffset=cumOffset+o;
    plot(((0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs))-2)*1e3,avLFPbySweep{i}*scl+cumOffset,'-k');
end
ylim(depth)
ylim auto
yl = ylim;
yyaxis left
ylim([depth(1)-(yl(2)-depth(2)) depth(2)+(depth(1)-yl(1))])

figure;
imagesc(timescale,depth,dataInt);
drawnow;
colormap(map);
ylim(depth);xlim([-500 500])

figure; hold on
cumOffset=chsp;
scl = 1;
avLFPbySweep = flipud(avLFPbySweep);
for i=3:length(avLFPbySweep)-2
    o=chsp;
    cumOffset=cumOffset+o;
    plot(((0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs))-2)*1e3,avLFPbySweep{i}*scl+cumOffset,'-k');
end

depth = [3 30];
figure;
imagesc(timescale,depth,data);
drawnow;
colormap(map);
ylim(depth);xlim([-500 500])






% Calculate CSD
% % csd=zeros(length(avLFPbySweep)-2,size(avLFPbySweep{1},2));
% % for i=2:length(avLFPbySweep)-1
% %     csd(i-1,:)=avLFPbySweep{i+1}-avLFPbySweep{i}+avLFPbySweep{i-1};
% % end
% % csd=zeros(length(avLFPbySweep)-4,size(avLFPbySweep{1},2));
% for i=3:length(avLFPbySweep)-2
%     csd(i-2,:)=avLFPbySweep{i+2}-avLFPbySweep{i}+avLFPbySweep{i-2};
% end
% 
% % Smoothing
% % for j=1:50
% %     for i = 1:size(csd,1)
% %         csd(i,:) = smooth(csd(i,:),17);
% %     end
% % end
%     
% % Plot CSD
% subplot(2,1,2);
% % heatmap(csd,'Standardize','column');
% heatmap(csd);

% subplot(3,1,3);
