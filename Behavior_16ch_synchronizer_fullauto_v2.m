
% clear
close all

avi_files = getfilenamese(pwd,'*.avi')';
h5_files = getfilenamese(pwd,'*.h5')';

master.C1 = 'on'; % create adjusted video for C1
master.C2 = 'on'; % create adjusted video for C2
master.combined = 'on';  % create adjusted video for C1 & C2 combined

for ff = [1,3:22]
%     try
    clearvars -except master avi_files h5_files ff 

rec_no = h5_files{ff}(12:15);

temp = strfind(avi_files,rec_no);


sel = cellfun(@isempty,temp);
sel = ~sel;

temp = avi_files(sel);
test = [];
if length(temp)>2
    for ii = 1:length(temp)
       test(ii) = ~isempty(strfind(temp{ii}(1:15),rec_no));
    end
temp = temp(~~test);
end

vid{1} = temp{1};
vid{2} = temp{2};

temp = strfind(master.AVI,vid{1});
temp = ~cellfun(@isempty,temp);
masterSel(1) = find(temp == 1);

temp = strfind(master.AVI,vid{2});
temp = ~cellfun(@isempty,temp);
masterSel(2) = find(temp == 1);

% r = zeros(16,2);

% r(1,1) = 58; % set frames for detecing LED timer, cam 1
% r(1,2) = 56;

dr = [0,0]; % x, y

Orient = 'horizontal';
% Orient = 'vertical';



Fs = 50000; % sampling rate
FR = 200; % frame rate

coords = [];

fileno = vid{1}(13:15);
% fileno(2) = '0';
% fileno = vid{1}(17:18);

% a = num2str(fileno);
tags = '';
f = strcat(vid{1}(1:11),tags,vid{1}(12:15),'.h5');
fsweep = strcat('/sweep_0',fileno,'/analogScans');
data = h5read(f,fsweep);
data = double(data);

filena = f;

%
% measure LED intensities
for mm = 1:2

    clearvars -except  master masterSel data delay h5_files fsweep fileno r vid tags mm filena coords summary dr FR Fs Orient rect rectSpike  rec_no avi_files

    if mm == 1
        vidFile = vid{1};
    else
        vidFile = vid{2};
    end

    rect = master.ROI(masterSel(mm),:);
    rectSpike = master.ROIspike(masterSel(mm),:);
    coords = master.coords(:,:,masterSel(mm));

    coords(:,1) = coords(:,1)-rect(1);
    coords(:,2) = coords(:,2)-rect(2);

    vidReader = VideoReader(vidFile);
    vidPlayer = vision.DeployableVideoPlayer;
    n.Frames = vidReader.FrameRate.*vidReader.Duration;

    temp = [];

    mult = [-1 0 1]';
    mult = repmat(mult,3,1);
    mult = [mult,[-1 -1 -1 0 0 0 1 1 1]'];
    I{1} = read(vidReader,1);
    x = size(I{1},1);
    y = size(I{1},2);

%     clear temp roi
%     roi = coords(:,:,mm);

%     for ff = 1:5
%         sel.row(:,ff) = (roi(ff,1)-2):(roi(ff,1)+2);
%         sel.col(:,ff) = (roi(ff,2)-2):(roi(ff,2)+2);
%     end

    hw = waitbar(0,strcat('Collecting LED timing from cam',num2str(mm)));

    rect = round(rect);
    rectSpike = round(rectSpike);

    for ff = 1:n.Frames
        I = read(vidReader,ff);
            if (rect(2)+rect(4))>size(I,1)
                while (rect(2)+rect(4))>size(I,1)
                rect(4) = rect(4)-1;
                end
            end
            ROI = I(rect(:,2):(rect(:,2)+rect(:,4)),rect(:,1):(rect(:,1)+rect(:,3)),:);
            ROI = imgaussfilt(ROI,0.5);
            ROI = ROI(:,:,2);
            if (rectSpike(2)+rectSpike(4))>size(I,1)
%                 while (rectSpike(2)+rectSpike(4))>size(I,1)
                rectSpike(4) = rectSpike(4)-1;
%                 end
            end
            ROIspike = I(rectSpike(:,2):(rectSpike(:,2)+rectSpike(:,4)),rectSpike(:,1):(rectSpike(:,1)+rectSpike(:,3)),:);
            ROIspike = imgaussfilt(ROIspike,0.5);
            ROIspike = ROIspike(:,:,2);
        fr = ff./n.Frames;
%         test = test(:,:,2);
        LED(ff,:) = ROI(sub2ind(size(ROI),round(coords(:,2)),round(coords(:,1))));
        LEDspike(ff,:) = mean(ROIspike(:));
%         LED(ff,1) = mean(mean(test(sel.col(:,1),sel.row(:,1))));
%         LED(ff,2) = mean(mean(test(sel.col(:,2),sel.row(:,2))));
%         LED(ff,3) = mean(mean(test(sel.col(:,3),sel.row(:,3))));
%         LED(ff,4) = mean(mean(test(sel.col(:,4),sel.row(:,4))));
%         test = read(vidReader,ff);
%         test = test(:,:,3);
%         LED(ff,5) = mean(mean(test(sel.col(:,5),sel.row(:,5))));
        %     LED(ff,2) = nanmean(test(ind(:,2)));
        %     LED(ff,3) = nanmean(test(ind(:,3)));
        %     LED(ff,4) = nanmean(test(ind(:,4)));
        if mod(ff,floor(n.Frames/10))<1e-2
            waitbar(ff/n.Frames,hw);
        end

    end
    close(hw)
    summary.LED{mm} = LED;
    summary.LEDspike{mm} = LEDspike;
    close all
end
% generate thresholded LED values
for mm = 1:2

    if mm == 1
        vidFile = vid{1};
    else
        vidFile = vid{2};
    end

    vidReader = VideoReader(vidFile);
    vidPlayer = vision.DeployableVideoPlayer;
    n.Frames = vidReader.FrameRate.*vidReader.Duration;
    
    LED = summary.LED{mm};
    LEDspike = summary.LEDspike{mm};

    % for ff = 1:n.Frames
    figure,
    ax1 = subplot(2,1,1);plot(LED)
    ylim([0 255])
    ax2 = subplot(2,1,2);plot(LEDspike);
    linkaxes([ax1 ax2],'x');
%     xlim([400 450])
    
    pause(0.5)
    % threshold light data

    th.timer = 200;%input('timer threshold?');
    th.trig = mean(LEDspike)+10;% input('trigger threshold?');
%     clc
    clear thresh
    thresh = LED>th.timer;
    threshSpike = LEDspike>th.trig;
%     thresh(:,5) = LED(:,5)>th.trig;
    thresh = double(thresh);
    threshSpike = double(threshSpike);
    close all

    %
    %
    % ax5 = subplot(3,1,3);plot(thresh);
    % % thresh = thresh;
    % % thresh(:,[2,4]) = -thresh(:,[2,4]);
    % % plot(thresh)
    %
    % ax3 = subplot(3,1,1);plot(cumsum(thresh(:,1)+thresh(:,2)))
    % ax4 = subplot(3,1,2);plot(cumsum(thresh(:,3)+thresh(:,4)))
    %
    %
    % linkaxes([ax3 ax4 ax5],'x')
% for loop model comparison

    clear rec
    % thresh = LED>205;
    if n.Frames <3000
        L = 2000;
    elseif n.Frames >2000
        L = 2992;
    end
%     thresh = thresh(:,1:5);
    
    standard = zeros([L,16]);
%         standard = zeros([2000,16]);

%     temp = [1 0 0 0]';
temp = zeros(16,1);
temp(1) = 1;
    
    standard(:,1) = repmat(temp,floor(L./16),1);

    for ii = 2:16
        standard(:,ii) = circshift(standard(:,ii-1),1);
    end
    standard(:,17) = 0;

    rec = [thresh,threshSpike];
    times.skip = [];
    for ff = 1:L
        %     cla
        %     plot(rec(:,4)./2,'r')
        %     hold on
        %     plot(standard(:,4),'b')
        %     xlim([-10 10]+ff)
        %     plot([ff,ff],[0 1],'--g')
        %     pause(0.2)
        if sum(rec(ff,1:16)) == 0
            if ff == 1
                rec(ff,1) = 1;
            else
                temp = find(rec(ff-1,1:16) == 1) ;
            end
            if temp == 16
                rec(ff,1) = 1;
            else
                rec(ff,temp+1) = 1;
            end
        end

        if standard(ff,2) ~= rec(ff,2)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,3) ~= rec(ff,3)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,4) ~= rec(ff,4)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,5) ~= rec(ff,5)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,6) ~= rec(ff,6)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,7) ~= rec(ff,7)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,8) ~= rec(ff,8)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,9) ~= rec(ff,9)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,10) ~= rec(ff,10)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,11) ~= rec(ff,11)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,12) ~= rec(ff,12)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,13) ~= rec(ff,13)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,14) ~= rec(ff,14)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,15) ~= rec(ff,15)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
        elseif standard(ff,16) ~= rec(ff,16)
            rec = [rec(1:(ff-1),:);standard(ff,:);rec(ff:end,:)];
            times.skip = [times.skip;ff];
            %         disp(length(times.skip));
        else
        end
    end

    temp = [0;diff(times.skip)];
    temp(temp ~= 1) = 0;
    mult = 0:(length(times.skip)-1);
    
    times.skip = times.skip-mult'+temp;
    times.raw = zeros(round(n.Frames),1);
    times.raw(times.skip) = 1;
    
    summary.skipped{mm} = times.skip;
    
    disp(strcat(num2str(length(rec)),' total frames for cam',num2str(mm)));
    
    save(strcat(filena(1:18),'_rec',num2str(mm),'.mat'));
    close all
    summary.rec{mm} = rec;
end


% synchronize frames to traces

clearvars -except master masterSel data delay h5_files fsweep fileno r vid tags mm filena coords summary Fs Orient FR  rec_no avi_files


times.cam = data(:,2)>6000;
times.cam = diff(times.cam)>0;
times.cam = find(times.cam == 1);
times.cam = times.cam./Fs;



times.trig = data(:,3)>3000;
times.trig = diff(times.trig)>0;
times.trig = find(times.trig == 1);
times.trig = times.trig./Fs;

times.camAdj = times.cam-times.cam(1);
times.trigAdj = times.trig;
% times.trigAdj = times.trig-times.cam(1);

% times.trigAdj = times.trigAdj(1:end);

% plot(times.trigAdj)

trace.trig = zeros(length(data),1);
% trace.trig(round((times.trig(2:end)-0.0548).*Fs)) = 1;
% trace.trig(round((times.trigAdj(1:end)).*Fs)) = 1;% -0.0548

times.trace = (1:length(trace.trig))./Fs;
times.trace = times.trace;
% plot(times.trace,trace.trig)

rec0001 = summary.rec{1};
rec0002 = summary.rec{2};

times.C1Spike = find(rec0001(:,end) == 1)./FR;
times.C2Spike = find(rec0002(:,end) == 1)./FR;

% delayC1 = (times.trig-0.1-times.C1Spike);
% delayC2 = (times.trig-0.1-times.C2Spike);

clear rec;

figure

ax1 = subplot(3,1,1); plot(times.trace,data(:,1));
grid on
grid minor

ax2 = subplot(3,1,2); plot((1:length(rec0001)).*0.005,rec0001(:,end));
hold on
% plot(times.trace,trace.trig./2,'r')
plot(times.trace,double(data(:,3))./6000,'r');
grid on
grid minor

ax3 = subplot(3,1,3);plot((1:length(rec0002)).*0.005,rec0002(:,end))
hold on
plot(times.trace,double(data(:,3))./6000,'r');
% plot(times.trace,trace.trig./2,'r')
grid on
grid minor

linkaxes([ax1 ax2 ax3],'x')
saveas(gcf,strcat(vid{1}(1:15),'_trig&cams'));

C1 = zeros(length(data),1);
D = zeros(length(data),1);

mult = times.C1Spike*Fs;
mult = repmat(mult,1,25);
mult = mult + repmat(-12:12,size(mult,1),1);
mult = mult';
mult = mult(:);
mult = round(mult);

C1(mult) = 1;

mult = times.trig*Fs;
mult = repmat(mult,1,25);
mult = mult + repmat(-12:12,size(mult,1),1);
mult = mult';
mult = mult(:);
mult = round(mult);

D(mult) = 1;

[crossCorr, lags] = xcorr(C1,D,Fs*0.3);

[~, idx] = max(abs(crossCorr));
delay = lags(idx) / Fs;

delay = abs(delay.*Fs);

temp = double(data(:,1))./5000;
temp = temp - mean(temp);
temp = temp(delay:end);

figure
plot(times.trace,C1)
hold on
plot(times.trace(1:(length(D)-delay+1))',D(delay:end),'r')
xlim([0,15])
title(strcat('Delay = ',num2str(delay/Fs), 's'))
saveas(gcf,strcat(vid{1}(1:15),'_delay_CC_calc.jpg'));
saveas(gcf,strcat(vid{1}(1:15),'_delay_CC_calc.fig'));

summary.delay = delay;



% frame drop replacement
clearvars -except master masterSel data delay h5_files fsweep fileno r vid tags mm filena coords summary Fs Orient  rec_no avi_files wr

wr.FR = 50;

if contains(master.C1,'on') && contains(master.C2,'on')
    k = 1:2;
end

if contains(master.C1,'on') && contains(master.C2,'off')
    k = 1;
end
if contains(master.C1,'off') && contains(master.C2,'on')
    k = 2;
end


for mm = k

load(strcat(filena(1:18),'_rec',num2str(mm),'.mat'));

filename = strcat(vid{mm}(1:(length(vid{mm})-4)),'_SYNCD_',num2str(wr.FR),'FPS.avi');

video = VideoWriter(filename,'Motion JPEG AVI');
video.FrameRate = wr.FR;%vidReader.FrameRate;
v.Quality = 95;

open(video)

% n.Frames = vidReader.Duration.*vidReader.FrameRate;
clear img
% img = zeros([vidReader.Height,vidReader.Width,3,round(n.Frames)]);

% for ff = 1:n.Frames
% img(:,:,:,ff) = read(vidReader,ff);
% end

% oldFrames = 1;
% newFrames = 1;

kk = 1;
clear img
times.skip(end+1) = NaN;
addedFrames = 0;
frameNo = 0;
hw = waitbar(0,strcat('Writing syncd video ',num2str(mm)));
clear img
% figure
for ff = 1:n.Frames
    clear img
    if kk <= length(times.skip)
        if times.skip(kk) == ff
            while times.skip(kk) == ff
                if times.skip(kk) == 1
                    img(:,:,:,1) = read(vidReader,1);
                    writeVideo(video,img);
                    kk = kk+1;
%                     addedFrames = addedFrames + 1;
%                     frameNo = frameNo + 1;
%                     imshow(img)
%                     disp(strcat('ADDED FRAME, ',num2str(addedFrames)))
%                     pause(0.1)
                else
                    img = read(vidReader,ff-1);
                    writeVideo(video,img);
                    kk = kk+1;
%                     addedFrames = addedFrames + 1;
%                     frameNo = frameNo + 1;
%                     imshow(img)
%                     disp(strcat('ADDED FRAME, ',num2str(addedFrames)))
%                     pause(0.1)
                    %
                end
            end
        end
            img = read(vidReader,ff);
            writeVideo(video,img);
%             imshow(img)
%             frameNo = frameNo + 1;
%             disp(strcat('ORIGINAL FRAME, ',num2str(frameNo)))
%             pause(0.1)

    end
    

    
    if mod(ff,floor(n.Frames/10))<1e-2
          waitbar(ff/n.Frames,hw);
     end
    
end
    
    
mult = 0:(length(times.skip)-1);
times.added = times.skip + mult';
    
close(hw)
close(video);


end

if contains(master.combined,'on') 

% combine videos
clearvars -except master masterSel delay h5_files data fsweep fileno r vid tags mm filena coords summary Fs Orient  rec_no avi_files wr


vidSync{1} = strcat(vid{1}(1:(length(vid{1})-4)),'_SYNCD_',num2str(wr.FR),'FPS.avi');
vidSync{2} = strcat(vid{2}(1:(length(vid{2})-4)),'_SYNCD_',num2str(wr.FR),'FPS.avi');

    vidReader1 = VideoReader(vidSync{1});
    vidPlayer1 = vision.DeployableVideoPlayer;
    n.Frames1 = vidReader1.FrameRate.*vidReader1.Duration;

    vidReader2 = VideoReader(vidSync{2});
    vidPlayer2 = vision.DeployableVideoPlayer;
    n.Frames2 = vidReader2.FrameRate.*vidReader2.Duration;

heightMax = max([vidReader1.Height,vidReader2.Height]);
widthMax = max([vidReader1.Width,vidReader2.Width]);
n.framesMax = max([n.Frames1,n.Frames2]);


filename = strcat(vid{1}(1:15),'_COMBINED','_SYNCD_',num2str(wr.FR),'FPS.avi');

video = VideoWriter(filename,'Motion JPEG AVI');
v.Quality = 95;
video.FrameRate = wr.FR;%vidReader1.FrameRate;
open(video)

hw = waitbar(0,'Writing combined video...');

temp = heightMax - vidReader1.Height;

for ff = 1:n.framesMax

    if ff>n.Frames1
        img1(:,:,:) = 0;
    else
        clear img1
        img1(:,:,:) = read(vidReader1,ff);
        img1 = [img1;zeros([temp,vidReader1.Width,3])];
    end
    
    clear img2
    if ff>n.Frames2
         img2(:,:,:) = zeros(max([vidReader2.Height,vidReader1.Height]),vidReader2.Width,3);
    else
        img2(:,:,:) = read(vidReader2,ff);
        if size(img2,1)<vidReader1.Height
            img2(vidReader2.Height:vidReader1.Height,:,:) = 0;
        end
    end

    
    
if sum(Orient == 'horizontal') >9
img = [img2, img1];
elseif  sum(Orient == 'vertical') >9
    img = [img1; img2];
end
    writeVideo(video,img);
    
    if mod(ff,floor(n.framesMax/10))<1e-2
          waitbar(ff/n.framesMax,hw);
     end
end
close(video)
close(hw)

end

% template search for spiking
template = master.template;
% template = template-mean(template(1:5));

% corr = conv(template,temp)
% thresh =  3;% 7/18
% thresh = 0.23; % for 7/16. 
% thresh = 0.8; % for 6/8
% thresh = 1;
% thresh = 0.55; 
% thresh = 0.55 % for 5/15
% thresh = 0.6; % for 4/23. 
% thresh = 0.9; % for 8/2. 
% thresh = 3; % for 8/9
thresh = 0.2; % for 8/24
% thresh = 3; %for 4/11
% thresh = 0.6; %for 3/23
% thresh = 1.2; %for 3/28
% thresh = 0.8; %for 2/27
% thresh = 1; % for 2/14
% thresh = 2.5;

spikes.times = [];

temp = double(data(:,1))./5000;
temp = temp - mean(temp);
temp = temp(delay:end);

% Perform the template search
for i = 1:length(temp) - length(template) + 1
    % Calculate the baseline of the input vector in the current window
    windowBaseline = mean(temp(i:i+length(template)-1));
    
    % Subtract the baseline from the input vector in the current window
    windowData = temp(i:i+length(template)-1) - windowBaseline;
    
    % Calculate the correlation between the template and the window data
    correlation = sum(template .* windowData);
    
    % Check if the correlation exceeds the threshold
    if correlation > thresh
        % Add the index of the detected event to the event indices array
        spikes.times = [spikes.times, i];
    end
end

sel = diff(spikes.times)>1;
sel = [sel,0];
spikes.times = spikes.times(~~sel);

times.temp = (1:length(temp))./Fs;
x = zeros(length(spikes.times),1);



figure
plot(times.temp,temp,'black')
hold on
plot(spikes.times./Fs,x,'*g')


saveas(gcf,strcat(vid{1}(1:15),'_detected_spikes.fig'));

clear sel correlation heightMax i img img1 img2 tags x
    
csvwrite(strcat(vid{1}(1:15),'.txt'),temp);
% audiowrite(strcat(vid{1}(1:15),'_',num2str(wr.FR),'FPS.wav'),temp,wr.FR.*250); % should be *250 for 50 kHz, 100 for 20 kHz

temp = zeros(length(data),1);
temp = temp(delay:end);
temp(round(spikes.times)) = 1;

temp = movmean(temp,4).*4;

audiowrite(strcat(vid{1}(1:15),'_',num2str(wr.FR),'FPS_thresh.wav'),temp,wr.FR.*250);
    
clear ii L mm mult standard temp x y
save(strcat(filena(1:18),'_summary','.mat'));

%     catch
%     csvwrite('ERROR.txt',[]);
%     end
    
end


% 
%     %% highpass, or notch filter
% 
% temp = double(data(:,1))./5000;
% temp = temp - nanmean(temp);
% % temp = temp(delay:end);
%     
% % high pass filter
% d=fdesign.highpass('N,Fc',150,200,50000); % N,6dB,sampling rate
% designmethods(d);
% Hd = design(d);
% 
% temp = filter(Hd,temp(:,1));
% temp = circshift(temp,-75);
% 
% % notch filter, 60 Hz
% w0 = 60/(Fs/2); 
% bw = w0/30;
% [b a] = iirnotch(w0,bw,30);
% temp = filter(b,a,temp);
% 
% csvwrite(strcat(vid{1}(1:15),'.txt'),temp);
% % audiowrite(strcat(vid{1}(1:15),num2str(wr.FR)*250,'FPS.wav'),temp,wr.FR.*250);
% audiowrite(strcat(vid{1}(1:15),'_',num2str(wr.FR),'FPS.wav'),temp,wr.FR.*250);
% 
% 
% % fvtool(a,b)
% % plot(double(raw(:,1))./5000)
% % hold on
% % plot(temp./5000,'r')
% 
% % data(:,1) = temp;
% 
% %% write thresholded data from detected spikes
% 
% spikes.times = spikes.times./1000;
% spikes.times = spikes.times.*Fs;
% % spikes.times = spikes.times - delay;
% 
% temp = zeros(length(data),1);
% temp = temp(delay:end);
% temp(round(spikes.times)) = 1;
% 
% audiowrite(strcat(vid{1}(1:15),'_',num2str(wr.FR),'FPS_thresh.wav'),temp,wr.FR.*250);
% 
% clear ii L mm mult standard temp x y
% save(strcat(filena(1:18),'_summary','.mat'));
%     