%% generate spreadsheet with files and frames

labl = 'finger_drag_sandpaper';

txt_files = getfilenamese(pwd,'*.txt')';

temp = strfind(txt_files,labl);
sel = cellfun(@isempty,temp);
sel = ~sel;

txt_sel = txt_files(sel);

frames = [];
matfiles = {};

for ff = 1:length(txt_sel)
    clear temp take
    temp = readcell(txt_sel{ff});
%     take = temp(1:2:end);
%     take(:,2) = temp(2:2:end);
    matfiles = [matfiles;temp(:,1)];
    frames = [frames;cell2mat(temp(:,3))];
end

%% load data from spreadsheet

% t_wind = [-40 40]; % step cycle 
% t_wind = [-60 60]; % digging 
% t_wind = [-100 100]; % wet dog shake
% t_wind = [-200 300]; % lift
% t_wind = [-40 40]; % finger tap
t_wind = [-100 150]; % drag
% 

u = unique(matfiles);
window = {};
labels = {};

for ff = 1:length(u)
    load(strcat(u{ff}(1:15),'.h5_summary.mat'),'delay','Fs','spikes');
    times.load = frames(ismember(matfiles,u{ff}),:);

    for ii = 1:size(times.load,1)
    
    f = times.load(ii,:);
    f = [f f];
    f = (f(1)+t_wind(1)):(f(2)+t_wind(2));

    temp = zeros(Fs*17,1);
    temp = temp(delay:end);
    temp(round(spikes.times)) = 1;

    t = f*0.005;
    t = [t(1),max(t)];
    t = round(t.*Fs);

    if t(1)>0 && t(2)<length(temp)
    window = [window; temp(t(1):(t(end)))];
    labels = [labels; strcat(u{ff},num2str(f(1)),':',num2str(f(2)))];
    else
        strcat('drop @',num2str(ff),'-',num2str(ii),': ',num2str(times.load(ii)))
    end

    end

clearvars -except u frames matfiles window labels txt_sel t_wind shift backup labl

end


b = cellfun(@length,window);
b = max(b);

for ii = 1:length(window)
    if length(window{ii})<b
        window{ii} = [window{ii};nan(b-length(window{ii}),1)];
    end
end

data = cell2mat(window');


% generate rasters and histogram
Fs = 50000;
raster.x = [];
raster.y = [];
% figure
% for f = 1:size(data,2)
%     t.spikes{f} = find((diff(data(:,f)>0.5)>0) == 1)./length(nonzeros(~isnan(data(:,f))));
%     raster.x = [raster.x;t.spikes{f}];
%     raster.y = [raster.y;repmat(f,length(t.spikes{f}),1)];
% end
% plot(raster.x,raster.y,'|black')
% xlim([0 1])
% ylim([0 max(raster.y)+1])


raster.x = [];
raster.y = [];

figure 
for f = 1:size(data,2)
    t.spikes{f} = find((diff(data(:,f)>0.5)>0) == 1);
    raster.x = [raster.x;t.spikes{f}];
    raster.y = [raster.y;repmat(f,length(t.spikes{f}),1)];
end
plot(raster.x./Fs./0.005 + t_wind(1),raster.y,'|black'),hold on
plot([0 0],[0 max(raster.y)],'--r')
set(gca,'ydir','reverse')
% xlim([0 1])
ylim([0 max(raster.y)+1])
xlabel('Frames')

figure,hold on
temp = hist(raster.x./Fs/0.005,(t_wind(1):t_wind(2))-t_wind(1))./max(raster.y);

% temp = hist(raster.x./Fs,0:0.005:1.5)./max(raster.y);
times.hist = (1:length(temp))./0.005;
temp = temp./0.005;
bar((t_wind(1):t_wind(2)),temp,'facecolor','black')
plot([0 0],[0 max(temp)],'--r')


%%  generate figures for export for step cycle

% exclusions

sel = [];


kill = ones(length(raster.y),1);
if ~isempty(sel)
    for ff = 1:(length(sel))
        temp = raster.y ~= sel(ff);
        kill = kill.*temp;
    end
end

raster.xSel = raster.x(~~kill);
raster.ySel = raster.y(~~kill);

if ~isempty(sel)
    raster.ySel(raster.ySel>max(sel)) = raster.ySel(raster.ySel>max(sel))-(max(sel)-min(sel));
end

figure
set(gcf,'position',[440 395 524 488]);
subplot(2,1,1)
plot(raster.xSel./Fs+t_wind(1)*0.005,raster.ySel,'|black'),hold on
plot([0 0],[0 max(raster.y)],'--r')
set(gca,'ydir','reverse')
set(gca,'fontsize',14)
ylim([0 max(raster.ySel)+1])
ylabel('Event #')
box off

xlim([-0.2 0.2])


temp = hist(raster.xSel./Fs,(t_wind(1):t_wind(2))*0.005-t_wind(1)*0.005)./max(raster.ySel);

subplot(2,1,2), hold on
times.hist = (1:length(temp))./0.005;
temp = temp./0.005;
bar((t_wind(1):t_wind(2))*0.005,temp,'facecolor','black')
plot([0 0],[0 200],'--r')
set(gca,'fontsize',14)
ylabel('Hz')

xlim([-0.2 0.2])

xlabel('Time (s)')

% saveas(gcf,strcat(txt_sel{1}(1:11),labl,'.fig'))
% saveas(gcf,strcat(txt_sel{1}(1:11),labl,'.jpg'),'jpg')


%%  generate figures for export for drag

% exclusions

sel = [];


kill = ones(length(raster.y),1);
if ~isempty(sel)
    for ff = 1:(length(sel))
        temp = raster.y ~= sel(ff);
        kill = kill.*temp;
    end
end

raster.xSel = raster.x(~~kill);
raster.ySel = raster.y(~~kill);

if ~isempty(sel)
    raster.ySel(raster.ySel>max(sel)) = raster.ySel(raster.ySel>max(sel))-(max(sel)-min(sel));
end

figure
set(gcf,'position',[440 395 524 488]);
subplot(2,1,1)
plot(raster.xSel./Fs+t_wind(1)*0.005,raster.ySel,'|black'),hold on
plot([0 0],[0 max(raster.y)],'--r')
set(gca,'ydir','reverse')
set(gca,'fontsize',14)
ylim([0 max(raster.ySel)+1])
ylabel('Event #')
box off

xlim([-0.2 0.3])

temp = hist(raster.xSel./Fs,(t_wind(1):t_wind(2))*0.005-t_wind(1)*0.005)./max(raster.ySel);

subplot(2,1,2), hold on
times.hist = (1:length(temp))./0.005;
temp = temp./0.005;
bar((t_wind(1):t_wind(2))*0.005,temp,'facecolor','black')
plot([0 0],[0 200],'--r')
set(gca,'fontsize',14)
ylabel('Hz')

xlim([-0.2 0.3])

xlabel('Time (s)')

saveas(gcf,strcat(txt_sel{1}(1:11),labl,'_PROP.fig'))
saveas(gcf,strcat(txt_sel{1}(1:11),labl,'_PROP.jpg'),'jpg')



%% align events manually
figure
plot(raster.x./Fs./0.005 + t_wind(1),raster.y,'|black'),hold on
plot([0 0],[0 max(raster.y)],'--r')
set(gca,'ydir','reverse')

for ii = 1:max(raster.y)
    [x,y] = ginput(1);
    shift_frames(ii,:) = [x y];
end

backup.frames = frames;
frames = frames+shift_frames(:,1);
frames = round(frames);

