%% collect filenames and time stamps
t.h5 = {};
t.h5n = [ 0 0 0];
h5_files = getfilenamese(pwd,'*.h5')';

for ff = 1:length(h5_files)
    t.h5{ff} = dir(h5_files{ff});
    t.h5{ff} = t.h5{ff}.date;
    t.h5{ff} = t.h5{ff}((end-7):end);
    temp = (strrep(t.h5{ff},':',''));
    t.h5n(ff,:) = [str2num(temp(1:2)),str2num(temp(3:4)),str2num(temp(5:6))];
end

avi_files = getfilenamese(pwd,'*.avi')';

for ff = 1:length(avi_files)
    t.avi{ff} = dir(avi_files{ff});
    t.avi{ff} = t.avi{ff}.date;
    t.avi{ff} = t.avi{ff}((end-7):end);
    temp = (strrep(t.avi{ff},':',''));
    t.avin(ff,:) = [str2num(temp(1:2)),str2num(temp(3:4)),str2num(temp(5:6))];
end

%% rename
for ff = 1:length(avi_files)
    h = ismember(t.h5n(:,1),t.avin(ff,1));
    m = ismember(t.h5n(:,2),t.avin(ff,2));
    s = abs(t.h5n(:,3)-t.avin(ff,3))<4;
    sel = and(and(h,m),s);

    if sum(sel) == 1
        New_name = h5_files(sel);
        New_name = New_name{1}(1:15);
        New_name = [New_name,avi_files{ff}(16:end)];

        status = movefile(avi_files{ff},New_name);
    end
end


