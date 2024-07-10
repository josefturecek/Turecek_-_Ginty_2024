

avi_files = getfilenamese(pwd,'*.avi')';

master.ROI = zeros(1,4);
master.ROIspike = zeros(1,4);

master.AVI = {};


% temp = strfind(avi_files,'C1');
% sel = cellfun(@isempty,temp);
% sel = ~sel;
% 
avi_sel = avi_files;
% avi_sel = avi_files(sel);

% select only for videos that have not been analyzed

temp = strfind(avi_sel,'SYNCD');
sel = cellfun(@isempty,temp);
sel = ~sel;

% remove SYNCD videos
rem = avi_sel(sel);

avi_sel = avi_sel(~sel);

mirrors = 'off';

%%

for ii = 1:42

    clearvars -except ii master avi_files avi_sel mirrors manual_skip rotato

    vidFile = avi_sel{ii};

    vidReader = VideoReader(vidFile);
    vidPlayer = vision.DeployableVideoPlayer;
    n.Frames = vidReader.FrameRate.*vidReader.Duration;

    for ff = 1:300
        img(:,:,:,ff) = read(vidReader,ff);
    end

    temp = strfind(vidFile,'C1');

    if ~isempty(temp)
        subplot(2,2,1)
        I = max(img,[],4);
        imshow(I);
        title(vidFile)
    else
        subplot(2,2,2)
        I = max(img,[],4);
        imshow(I);
        title(vidFile)
    end

    if ~exist('manual_skip')
        % % check if previous ROI matches current; if so, skip
        if exist('master')
            if length(nonzeros(master.ROI(:,1)))>1
                mm = 1;

                rect = master.ROI(ii-2,:);
                try
                    ROI1_max = master.ROI_image{ii-2};
                    ROI2_max = I(rect(mm,2):(rect(mm,2)+rect(mm,4)),rect(mm,1):(rect(mm,1)+rect(mm,3)),:);
                catch
                    ROI2_max = ROI1_max;
                    ROI2_max(:) = 0;
                end


                c1 = double(ROI1_max(:));
                c2 = double(ROI2_max(:));

                if isempty(c1)
                    rect = master.ROI(ii-4,:);
                    try
                        ROI1_max = master.ROI_image{ii-4};
                        ROI2_max = I(rect(mm,2):(rect(mm,2)+rect(mm,4)),rect(mm,1):(rect(mm,1)+rect(mm,3)),:);
                    catch
                        ROI2_max = ROI1_max;
                        ROI2_max(:) = 0;
                    end

                    c1 = double(ROI1_max(:));
                    c2 = double(ROI2_max(:));
                end



                [X R2] = fit(c1,c2,'poly1');

                if R2.rsquare >0.8

                    master.ROIspike(ii,:) = master.ROIspike(ii-2,:);
                    master.coords(:,:,ii) = master.coords(:,:,ii-2);
                    master.ROI_image{ii} = master.ROI_image{ii-2};
                    disp('ROI the same as previous, took parameters from last');
                    skipp = 1;
                else
                    skipp = 0;
                end
            else
                skipp = 0;
            end
        else
            skipp = 0;
        end
    else

        clear manual_skip
        skipp = 0;
        
    end

            pause(0.1)
            if skipp ~= 1
            pas = input('mark LED? 1 = y, 9 = same as last');
            
            else 
                pas = 9;
            end

            if pas == 1
                rotato = input('rotate LED array in deg:');
                % user identify array region

                mm = 1;
                temp = drawrectangle;
                rect = temp.Position;

                temp = drawrectangle;
                rectSpike = temp.Position;

                ROI_max = I(rect(mm,2):(rect(mm,2)+rect(mm,4)),rect(mm,1):(rect(mm,1)+rect(mm,3)),:);
                ROI_all = img(rect(mm,2):(rect(mm,2)+rect(mm,4)),rect(mm,1):(rect(mm,1)+rect(mm,3)),:,:);

                master.ROI_image{ii} = ROI_max;
                % figure
                % imshow(ROI_max)

                % isolate LED pixels

                clear LED_max LED_ind F stat coords

                for ff = 1:size(img,4)
                    %     imshow(ROI_max-ROI_all(:,:,:,ff),[]);
                    F(:,:,:,ff) = (ROI_max-ROI_all(:,:,:,ff));
                end

                LED_max = max(F,[],4);
                % imshow(LED_max)

                se = strel('cube',4);
                % figureR
                for ff = 1:size(img,4)
                    LED_ind(:,:,:,ff) = LED_max-(ROI_max-ROI_all(:,:,:,ff));
                    %     Ibw = im2bw(LED_ind(:,:,:,ff));
                    Ibw = LED_ind(:,:,:,ff)>150;
                    Ibw = sum(Ibw,3);
                    %     Ibw = im2bw(LED_ind(:,:,:,ff));
                    Ibw = uint8(Ibw);
                    Ibw = imerode(Ibw,se);
                    Ilabel = bwlabel(Ibw);
                    if sum(Ilabel(:))>0
                        stat = regionprops(Ilabel,'centroid');
                        coords(ff,:) = stat.Centroid;
                    else
                        coords(ff,:) = [NaN,NaN];
                    end
                    %
%                         imshow(LED_ind(:,:,:,ff),[]);
%                         hold on
%                         plot(coords(ff,1),coords(ff,2),'or');
%                         pause(1)

                end

                clear ff Ibw Ilabel A B C D stats

                % cluster points and sort

                [k,ctr,~,~] = kmeans(coords,16);
                
%                 plot(ctr(:,1),ctr(:,2),'ob')
%                 hold on

                t = rotato.*pi./180;
                R = [cos(t), -sin(t); sin(t), cos(t)];

                ctr = ctr*R;
%                 plot(ctr(:,1),ctr(:,2),'or')
            

                Ix = discretize(ctr(:,1),4);
                Iy = discretize(ctr(:,2),4);

                % plot(ctr(:,1),ctr(:,2),'-oblack');
                % hold on
                % for ff = 1:4
                % plot(ctr(Ix ==ff,1),ctr(Ix == ff,2),'-xr');
                % pause(1)
                % plot(ctr(Iy ==ff,1),ctr(Iy == ff,2),'-gr');
                % pause(1)
                % end
                order = [];
                % for xx = 1:4
                %     for yy = 4:-1:1



                if strcmp(mirrors,'on')

                    for xx = 4:-1:1
                        for yy = 4:-1:1
                            temp = find(and(Ix == xx, Iy == yy) == 1);
                            order = [order;temp];
                            %         pause(1);
                        end
                    end

                else

                


                    for xx = 1:4
                        for yy = 4:-1:1
                            temp = find(and(Ix == xx, Iy == yy) == 1);
                            order = [order;temp];
                            %         pause(1);
                        end
                    end
                end


                ctr = ctr(order,:);

                t = rotato.*-pi./180;
                R = [cos(t), -sin(t); sin(t), cos(t)];

                ctr = ctr*R;


                % imshow(LED_max),hold on
                % plot(ctr(:,1),ctr(:,2),'-r')

                ctr(:,1) = ctr(:,1)+rect(1);
                ctr(:,2) = ctr(:,2)+rect(2);

                temp = strfind(vidFile,'C1');

                if ~isempty(temp)
                    subplot(2,2,3)
                    hold off
                    cla
                    imshow(I),hold on
                    plot(ctr(:,1),ctr(:,2),'-r')
                    title(vidFile)
                else
                    subplot(2,2,4)
                    hold off
                    cla
                    imshow(I),hold on
                    plot(ctr(:,1),ctr(:,2),'-r')
                    title(vidFile)
                end



            pause(0.2)

            master.ROI(ii,:) = rect;
            master.ROIspike(ii,:) = rectSpike;
            master.coords(:,:,ii) = ctr;

        else

            master.ROI(ii,:) = master.ROI(ii-2,:);
            master.ROIspike(ii,:) = master.ROIspike(ii-2,:);
            master.coords(:,:,ii) = master.coords(:,:,ii-2);
            master.ROI_image{ii} = master.ROI_image{ii-2};

        end

        master.AVI{ii,1} = avi_sel{ii};
        % close all
end
    % master.

    sel = cellfun(@isempty,master.AVI);
    for ff = 1:length(master.AVI)
        if sel(ff) == 1
            master.AVI{ff} = {'empty'};
        end
    end

    clearvars -except master
    temp = string(datetime("now"));
    temp = strrep(temp,':','-');
    save(strcat('avi_files_summary',temp,'.mat'));