clear
close all
clc

fs = 1/1.5;
%% load data for HbO & CCO
%fd_name = 'C:\Users\sadra\OneDrive - University of Texas at Arlington\Research\PhD\Aim 2\Data\Grouped\';
addpath('C:\Users\fizas\OneDrive - UT Arlington\Desktop\Research\Water Calculations_ Young Adults\Grouped\Raw_water\R808\post\')
fd_name = 'C:\Users\fizas\OneDrive - UT Arlington\Desktop\Research\Water Calculations_ Young Adults\Grouped\Raw_water\R808\post\';
input = 'Grouped';
output1 = 'Filtered_water';
output2 = 'Power';
freq1 = '0.005-0.2';
filter_type = 'filtfilt';
overlap = 180;

stim='R852';
time='post';

sub = [1:6,9,12:16,18,21:25,29:31,33:35];
for mm=1
%% loading data

%sub_fd = ['sub_',num2str(sub(mm))];
%sub_fd = ['sub_',num2str(1)];

%data = load([fd_name,input,'/',stim,'/',time,'/',sub_fd,'.mat']);
data = load(['sub_1.mat']);

%sub1_tls = data.Conc_all(:,1:289,:);
sub1_tls = data.Conc_all(:,:,:);  %280:560
x = (1:size(sub1_tls,2))/(fs*60);
%end

%%
for jj = 1:2 %1 is left and 2 is right
    for ii = 1:4
        xw = squeeze(sub1_tls(ii,:,jj));
        xw = xw(:);

        win = min(210, length(xw));
        if win < 10
            warning('Too short signal for Welch: ii=%d jj=%d', ii, jj);
            pxx = nan(129,1); f = linspace(0, fs/2, 129)';
        else
            noverlap = min(180, win-1);
            nfft = max(256, 2^nextpow2(win));

            [pxx,f] = pwelch(xw, win, noverlap, nfft, fs);
        end

        pxx_all(:,ii,jj) = pxx;
        f_all(:,ii,jj) = f;
        pxx_full(:,ii,jj,mm) = pxx;
    end
end

% for jj = 1:2 % right and left channel
% for ii = 1:4 % HbO, Hb, CCO and water
% 
%     [pxx,f] = pwelch(sub1_tls(ii,:,jj),210,overlap,round(1024*fs),fs);
% % [pxx2,f2] = cwt(sub1_tls(ii,1:290,jj),'morse',fs);
% % pxx2 = abs(pxx2).^2;
% % pxx2 = sum(pxx2,2);
% pxx_all(:,ii,jj) = pxx;
% f_all(:,ii,jj) = f;
% pxx_full(:,ii,jj,mm) = pxx_all(:,ii,jj);
% 
% % pxx2_all(:,ii,jj) = pxx2;
% % f2_all(:,ii,jj) = f2;
% % pxx2_full(:,ii,jj,mm) = pxx2_all(:,ii,jj);
% end
% end


%% FIR filter (ISO)
Endogenic = [0.004,0.006,0.019,0.021];
Neurogenic = [0.018,0.022,0.038,0.042];
Myogenic = [0.038,0.042,0.195,0.205];

a=[0 0 1 1 0 0];

[b_wide,d_wide] = butter(10,[0.005 0.1]/(fs/2),'bandpass');
[h,q]=freqz(b_wide,d_wide);

% Endogenic
f2_Endo=[0 Endogenic(1) Endogenic(2) Endogenic(3) Endogenic(4) fs/2]/(fs/2);
% b_Endo=firpm(n,f2_Endo,a);
% [b_Endo_high,d_Endo_high] = butter(8,0.005/(fs/2),'high');
% [b_Endo_low,d_Endo_low] = butter(6,0.02/(fs/2),'low');
[b_Endo,d_Endo] = butter(3,[0.005 0.02]/(fs/2),'bandpass');



% Neurogenic
f2_Neuro=[0 Neurogenic(1) Neurogenic(2) Neurogenic(3) Neurogenic(4) fs/2]/(fs/2);
% b_Neuro=firpm(n,f2_Neuro,a);
[b_Neuro,d_Neuro] = butter(4,[0.0205 0.039]/(fs/2),'bandpass');
[h,q]=freqz(b_Neuro,d_Neuro);
figure
freqz(b_Neuro);
% figure
% plot(f2_Neuro,a,q/pi,abs(h))
% legend('ideal','firpm Design')
% xlabel 'Radian Frequency (\omega/\pi)', ylabel 'Magnitude'

% Myogenic
f2_Myo=[0 Myogenic(1) Myogenic(2) Myogenic(3) Myogenic(4) fs/2]/(fs/2);
% b_Myo=firpm(n,f2_Myo,a);
[b_Myo,d_Myo] = butter(9,[0.041 0.15]/(fs/2),'bandpass');
[h,q]=freqz(b_Myo,d_Myo);
% figure
% freqz(b_Myo);
% figure
% plot(f2_Myo,a,q/pi,abs(h))
% legend('ideal','firpm Design')
% xlabel 'Radian Frequency (\omega/\pi)', ylabel 'Magnitude'

%% Apply filter

for ii = 1:2
for jj = 1:4
    if strcmp (filter_type,'filter') ==1
h_wide_tls(jj,:,ii) = filter(b_wide,d_wide,sub1_tls(jj,:,ii));    
h_Endo_tls(jj,:,ii) = filter(b_Endo,d_Endo,sub1_tls(jj,:,ii));
h_Neuro_tls(jj,:,ii) = filter(b_Neuro,d_Neuro,sub1_tls(jj,:,ii));
h_Myo_tls(jj,:,ii) = filter(b_Myo,d_Myo,sub1_tls(jj,:,ii));
    else 
h_wide_tls(jj,:,ii) = filtfilt(b_wide,d_wide,sub1_tls(jj,:,ii));    
h_Endo_tls(jj,:,ii) = filtfilt(b_Endo,d_Endo,sub1_tls(jj,:,ii));
h_Neuro_tls(jj,:,ii) = filtfilt(b_Neuro,d_Neuro,sub1_tls(jj,:,ii));
h_Myo_tls(jj,:,ii) = filtfilt(b_Myo,d_Myo,sub1_tls(jj,:,ii));        
    end
end
end

%% interpolate outliers
 h_Endo_tls2 = h_Endo_tls;
% figure()
for jj = 1:2
    for ii = 1:4
 constant_std = 2*std(h_Endo_tls(ii,:,jj),[],2);
%  constant_std = [1.1,1.1;0.11,0.11];
 constant_mean = mean(h_Endo_tls(ii,:,jj),2);
tick =0;
std_index = 1;
 while isempty(std_index) == 0 && tick < 10000 
 x1= find(abs(h_Endo_tls2(ii,:,jj))>constant_std);
       std_index = x1;
 for m = 1:length(x1)
 if x1(m)<10 || x1(m)>280    
      std_index(m) = 0;
 end
 std_index2 = nonzeros(std_index');
 end
if isempty(std_index)==0
 for n = -2:2
 h_Endo_tls2(ii,std_index2+n,jj) = 0.25*(h_Endo_tls(ii,std_index2+n-2,jj)+h_Endo_tls(ii,std_index2+n-1,jj)+h_Endo_tls(ii,std_index2+n+1,jj)+h_Endo_tls(ii,std_index2+n+2,jj));
 end
 tick = tick+1;
end
 end
% subplot(2,2,((ii-1)*ii)+jj)
% plot(1:size(h_Endo_tls2,2),h_Endo_tls2(ii,:,jj))
% title(mm)

    end
end

h_Neuro_tls2 = h_Neuro_tls;
% figure()
for jj = 1:2
    for ii = 1:4
 constant_std = 2*std(h_Neuro_tls(ii,:,jj),[],2);
%  constant_std = [1.1,1.1;0.11,0.11];
 constant_mean = mean(h_Neuro_tls(ii,:,jj),2);
tick =0;
std_index = 1;
 while isempty(std_index) == 0 && tick < 10000 
 x1= find(abs(h_Neuro_tls2(ii,:,jj))>constant_std);
       std_index = x1;
 for m = 1:length(x1)
 if x1(m)<10 || x1(m)>280    
      std_index(m) = 0;
 end
 std_index2 = nonzeros(std_index);
 end

 for n = -2:2
 h_Neuro_tls2(ii,std_index2'+n,jj) = 0.25*(h_Neuro_tls(ii,std_index2+n-2,jj)+h_Neuro_tls(ii,std_index2+n-1,jj)+h_Neuro_tls(ii,std_index2+n+1,jj)+h_Neuro_tls(ii,std_index2+n+2,jj));
 end
 tick = tick+1;
 end
% subplot(2,2,((ii-1)*ii)+jj)
% plot(1:size(h_Endo_tls2,2),h_Endo_tls2(ii,:,jj))
% title(mm)

h_Myo_tls2 = h_Myo_tls;
% figure()
for jj = 1:2
    for ii = 1:4
 constant_std = 2*std(h_Myo_tls(ii,:,jj),[],2);
%  constant_std = [1.1,1.1;0.11,0.11];
 constant_mean = mean(h_Myo_tls(ii,:,jj),2);
tick =0;
std_index = 1;
 while isempty(std_index) == 0 && tick < 10000 
 x1= find(abs(h_Myo_tls2(ii,:,jj))>constant_std);
       std_index = x1;
 for m = 1:length(x1)
 if x1(m)<10 || x1(m)>280    
      std_index(m) = 0;
 end
 std_index2 = nonzeros(std_index');
 end

 for n = -2:2
 %h_Myo_tls2(ii,std_index2+n,jj) = 0.25*(h_Myo_tls(ii,std_index2+n-2,jj)+h_Myo_tls(ii,std_index2+n-1,jj)+h_Myo_tls(ii,std_index2+n+1,jj)+h_Myo_tls(ii,std_index2+n+2,jj));
 end
 tick = tick+1;
 end
% subplot(2,2,((ii-1)*ii)+jj)
% plot(1:size(h_Endo_tls2,2),h_Endo_tls2(ii,:,jj))
% title(mm)

    end
end
h_wide_tls2 = h_wide_tls;
h_ISO(:,:,:,1) = h_Endo_tls2;
h_ISO(:,:,:,2) = h_Neuro_tls2;
h_ISO(:,:,:,3) = h_Myo_tls2;
% save([fd_name,output1,'/',freq1,'/',stim,'/',time,'/',sub_fd,'_',filter_type,'.mat'],'h_ISO');

    end
end

% % % %%
% % % %%%%%%%%%%%%%%%%%%% Power Analysis %%%%%%%%%%%%%%%%%%%%%%
% % % %%%% Endo
% % % 
% % % for jj=1:2
% % %     for ii = 1:4
% % %         if strcmp(time,'Post')==1
% % % pow(ii,jj,1) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.007 0.02],'psd'); 
% % %         else 
% % % pow(ii,jj,1) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.007 0.02],'psd');             
% % %         end
% % %     end
% % % end
% % % %%%% Neuro
% % % for jj=1:2
% % %     for ii = 1:4
% % %         if strcmp(time,'Post')==1
% % % pow(ii,jj,2) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.02 0.04],'psd'); 
% % %         else 
% % % pow(ii,jj,2) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.02 0.04],'psd');             
% % %         end
% % %     end
% % % end
% % % 
% % % %%%% Myo
% % % for jj=1:2
% % %     for ii = 1:4
% % %         if strcmp(time,'Post')==1
% % % pow(ii,jj,3) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.04 0.2],'psd'); 
% % %         else 
% % % pow(ii,jj,3) = bandpower(pxx_all(:,ii,jj),f_all(:,ii,jj),[0.04 0.2],'psd');             
% % %         end
% % %     end
% % % end
% % % % save([fd_name,output2,'/',freq1,'/',stim,'/',time,'/',sub_fd,'.mat'],'pow');
end
%% Normalize values

for ch = 1:4
    for channel = 1:2
h_wide_Norm(ch,:,channel,mm) = h_wide_tls2(ch,:,channel)./(max(h_wide_tls2(ch,:,channel))-min(h_wide_tls2(ch,:,channel)));
h_Endo_Norm(ch,:,channel,mm) = h_Endo_tls2(ch,:,channel)./(max(h_Endo_tls2(ch,:,channel))-min(h_Endo_tls2(ch,:,channel)));
h_Neuro_Norm(ch,:,channel,mm) = h_Neuro_tls2(ch,:,channel)./(max(h_Neuro_tls2(ch,:,channel))-min(h_Neuro_tls2(ch,:,channel)));
h_Myo_Norm(ch,:,channel,mm) = h_Myo_tls2(ch,:,channel)./(max(h_Myo_tls2(ch,:,channel))-min(h_Myo_tls2(ch,:,channel)));
    end
end
    for channel = 1:2
HbT_wide(:,channel,mm) = h_wide_tls2(1,:,channel)+h_wide_tls2(2,:,channel)
HbT_Endo(:,channel,mm) = h_Endo_tls2(1,:,channel)+h_Endo_tls2(2,:,channel);
HbT_Neuro(:,channel,mm) = h_Neuro_tls2(1,:,channel)+h_Neuro_tls2(2,:,channel);
HbT_Myo(:,channel,mm) = h_Myo_tls2(1,:,channel)+h_Myo_tls2(2,:,channel);
    end
    for channel = 1:2
h_wide_Norm(5,:,channel,mm) = (HbT_wide(:,channel)./(max(HbT_wide(:,channel))-min(HbT_wide(:,channel))))';
h_Endo_Norm(5,:,channel,mm) = (HbT_Endo(:,channel)./(max(HbT_Endo(:,channel))-min(HbT_Endo(:,channel))))';
h_Neuro_Norm(5,:,channel,mm) = (HbT_Neuro(:,channel)./(max(HbT_Neuro(:,channel))-min(HbT_Neuro(:,channel))))';
h_Myo_Norm(5,:,channel,mm) = (HbT_Myo(:,channel)./(max(HbT_Myo(:,channel))-min(HbT_Myo(:,channel))))';
    end
    for channel = 1:2
h_wide_Norm(6,:,channel,mm) = h_wide_Norm(4,:,channel)-h_wide_Norm(5,:,channel);
h_Endo_Norm(6,:,channel,mm) = h_Endo_Norm(4,:,channel)-h_Endo_Norm(5,:,channel);
h_Neuro_Norm(6,:,channel,mm) = h_Neuro_Norm(4,:,channel)-h_Neuro_Norm(5,:,channel);
h_Myo_Norm(6,:,channel,mm) = h_Myo_Norm(4,:,channel)-h_Myo_Norm(5,:,channel);
    end

    
    
h_wide_Norm_all = reshape(h_wide_Norm,6,[]);    
h_Endo_Norm_all = reshape(h_Endo_Norm,6,[]);
h_Neuro_Norm_all = reshape(h_Neuro_Norm,6,[]);
h_Myo_Norm_all = reshape(h_Myo_Norm,6,[]);
figure 
for ch = 1:6
    subplot(3,2,ch)
plot((0:1/fs:420)/60,h_Neuro_Norm(ch,1:281,1),'b.-','linewidth',2)
xlim([0 7])
set(gca,'fontsize',30,'fontweight','bold','linewidth',2)
% ylim([-0.35 0.35])
if ch == 1
    title('Neurogenic')
    ylabel('\Delta[HbO] ')
%     xlabel('Time (min)')
end
if ch == 2
%     title('Endogenic')
    ylabel('\Delta[Hb]')
    title('Neurogenic')
end
if ch == 3
    ylabel('\Delta[CCO]')

end
if ch == 4
    ylabel('\DeltaH2O')
end
if ch == 5
    ylabel('\Delta[HbT]')
    xlabel('Time (min)')    
end
if ch == 6
    ylabel('\DeltaH2O (Free)')
    xlabel('Time (min)') 
end
end
hold on
subplot(3,2,6)
plot((0:1/fs:420)/60,h_Endo_Norm(5,1:281,1),'r.-','linewidth',2)
legend('H2O (Free)','\Delta[HbT]')
hold off
%test
%% Correlation Calculation

%% Correlation Calculation

%% Correlation Calculation

% Wide Plot
figure
[R_wide(:,:), P_wide(:,:)] = corrplot([h_wide_Norm_all(3,:)', h_wide_Norm_all(6,:)']);
Z_wide(:,:) = atanh(R_wide(:,:));
% Remove the default text objects (r values)
ax = gca; % Get current axis
hText = findobj(ax, 'Type', 'Text'); % Find all text objects
delete(hText); % Delete text objects

% Endo Plot
figure
[R_Endo(:,:), P_Endo(:,:)] = corrplot([h_Endo_Norm_all(5,:)', h_Endo_Norm_all(6,:)']);
Z_Endo(:,:) = atanh(R_Endo(:,:));
% Remove text objects
ax = gca; % Get current axis
hText = findobj(ax, 'Type', 'Text'); % Find all text objects
delete(hText); % Delete text objects

% Neuro Plot
figure
[R_Neuro(:,:), P_Neuro(:,:)] =  corrplot([h_Neuro_Norm_all(3,:)', h_Neuro_Norm_all(6,:)']);
Z_Neuro(:,:) = atanh(R_Neuro(:,:));
% Remove text objects
ax = gca; % Get current axis
hText = findobj(ax, 'Type', 'Text'); % Find all text objects
delete(hText); % Delete text objects

% Myo Plot
figure
[R_Myo(:,:), P_Myo(:,:)] =  corrplot([h_Myo_Norm_all(3,:)', h_Myo_Norm_all(6,:)']);
Z_Myo(:,:) = atanh(R_Myo(:,:));
% Remove text objects
ax = gca; % Get current axis
hText = findobj(ax, 'Type', 'Text'); % Find all text objects
delete(hText); % Delete text objects


%%
%% Correlation Calculation

  figure
    [R_wide(:,:),P_wide(:,:)] = corrplot([h_wide_Norm_all(3,:)',h_wide_Norm_all(6,:)']);
    Z_wide(:,:) = atanh(R_wide(:,:));

    figure
    [R_Endo(:,:),P_Endo(:,:)] = corrplot([h_Endo_Norm_all(3,:)',h_Endo_Norm_all(6,:)']);
    Z_Endo(:,:) = atanh(R_Endo(:,:));
%     [r_Endo,Lag_Endo] = xcorr(h_Endo_Norm(5,1:281,side),h_Endo_Norm(6,1:281,side));
%     figure
%     stem(Lag_Endo,r_Endo)



    figure
    [R_Neuro(:,:),P_Neuro(:,:)] = corrplot([h_Neuro_Norm_all(3,:)',h_Neuro_Norm_all(6,:)']);
    Z_Neuro(:,:) = atanh(R_Neuro(:,:));
    % [r_Neuro,Lag_Neuro] = xcorr(h_Neuro_Norm(5,1:281,side),h_Neuro_Norm(6,1:281,side));
    % figure
    % stem(Lag_Neuro,r_Neuro)



    figure
    [R_Myo(:,:),P_Myo(:,:)] = corrplot([h_Myo_Norm_all(3,:)',h_Myo_Norm_all(6,:)']);
    Z_Myo(:,:) = atanh(R_Myo(:,:));
%     [r_Myo,Lag_Myo] = xcorr(h_Myo_Norm(5,1:281,side),h_Myo_Norm(6,1:281,side));
%     figure
%     stem(Lag_Myo,r_Myo)    

%% slope calculation
figure
scatter(h_Endo_Norm_all(5,:)',h_Endo_Norm_all(6,:)')
hold on;
coefficients = polyfit(h_Endo_Norm_all(5,:)',h_Endo_Norm_all(6,:)', 1);
fittedY = polyval(coefficients, h_Endo_Norm_all(5,:)');
plot(h_Endo_Norm_all(5,:)', fittedY, '-r', 'LineWidth', 2);
slope = coefficients(1);
intercept = coefficients(2);
slopeText = sprintf('y = %.2fx + %.2f', slope, intercept);
text(max(h_Endo_Norm_all(5,:)') * 0.6, max(h_Endo_Norm_all(1,:)') * 0.9, slopeText, 'FontSize', 12, 'Color', 'blue');
hold off;
figure
scatter(h_Neuro_Norm_all(5,:)',h_Neuro_Norm_all(6,:)')
hold on;
coefficients = polyfit(h_Neuro_Norm_all(5,:)',h_Neuro_Norm_all(6,:)', 1);
fittedY = polyval(coefficients, h_Neuro_Norm_all(5,:)');
plot(h_Neuro_Norm_all(5,:)', fittedY, '-r', 'LineWidth', 2);
slope = coefficients(1);
intercept = coefficients(2);
slopeText = sprintf('y = %.2fx + %.2f', slope, intercept);
text(max(h_Neuro_Norm_all(5,:)') * 0.6, max(h_Neuro_Norm_all(1,:)') * 0.9, slopeText, 'FontSize', 12, 'Color', 'blue');
hold off;

figure
scatter(h_Myo_Norm_all(5,:)',h_Myo_Norm_all(6,:)')
hold on;
coefficients = polyfit(h_Myo_Norm_all(5,:)',h_Myo_Norm_all(6,:)', 1);
fittedY = polyval(coefficients, h_Myo_Norm_all(3,:)');
plot(h_Myo_Norm_all(5,:)', fittedY, '-r', 'LineWidth', 2);
slope = coefficients(1);
intercept = coefficients(2);
slopeText = sprintf('y = %.2fx + %.2f', slope, intercept);
text(max(h_Myo_Norm_all(5,:)') * 0.6, max(h_Myo_Norm_all(1,:)') * 0.9, slopeText, 'FontSize', 12, 'Color', 'blue');
hold off;
% h_Norm_All(:,:,1)=h_Endo_Norm_all;
% h_Norm_All(:,:,2)=h_Neuro_Norm_all;
% h_Norm_All(:,:,3)=h_Myo_Norm_all;
% save([fd_name,'Results\',time, 'sub_1_R852.mat'],'h_Norm_All');  %5 and 6 for normalized water and water free, last 1 2 3 for endo neuro and myo