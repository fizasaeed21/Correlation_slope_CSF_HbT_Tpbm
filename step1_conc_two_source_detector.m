clear
clc
%step 1 in bbnirs calculation
%The first code to get concentrations of HbO, Hb, CCO and water
% also keep all the finction files in the same folder as the code
%This code is for two sources and two detectors

%% save raw spectrum
nb_row_cam = 402;

startw = 750;
endw = 1000;

r = 3; % distance between light source and detector (SD) (cm)
name_list = {'left','right'};
legend_list = {'Left','Right'};
cl_list = {'r.-','b.-'};
%% read camera or QEPro data

% fd_name = '/Users/nghitruong/Documents/NghiTruong_data/downloads/OneDrive_6_06-12-2020';
fd_name = 'C:\Users\fizas\OneDrive - UT Arlington\Desktop\Research\AD Study\';
sub_fd1 = 'AD-PNTS\PNT8\671\'; %should be QEP016571_
name1 = '671';
sub_fd2 = 'AD-PNTS\PNT8\972\'; %should be QEP009772_
name2 = '972';

file_list1 = dir([fd_name,sub_fd1,'*.csv']);
bbNIRS_type = 1; % camera
if isempty(file_list1)
    file_list1 = dir([fd_name,sub_fd1,'*.txt']);
    file_list11 = {file_list1.name}';
    file_list2 = dir([fd_name,sub_fd2,'*.txt']);
    file_list22 = {file_list2.name}';
    bbNIRS_type = 2; % QEPro
end


[~,reindex1] = sort(str2double(regexp({file_list1.name},'\d+','match','once')));
file_list1 = file_list1(reindex1);

[~,reindex2] = sort(str2double(regexp({file_list2.name},'\d+','match','once')));
file_list2 = file_list2(reindex2);

if exist([fd_name,sub_fd1,name1,'_raw_data.mat'],'file')==2
    load([fd_name,sub_fd1,name1,'_raw_data.mat']);
    

else  

 

            for i = 1:length(file_list1)
                table1 = read_txt_QEPro1([fd_name,sub_fd1,file_list1(i).name]);
                spectrum_1(:,i,1) = table1(find(fix(table1(:,1))==(startw),1):...
                    find(fix(table1(:,1))==(endw),1,'last'),2);
                if i==1
                    wavelength = table1(find(fix(table1(:,1))==(startw),1):...
                        find(fix(table1(:,1))==(endw),1,'last'),1);
                end
            end
            

    
            for i = 1:length(file_list2)
                table2 = read_txt_QEPro2([fd_name,sub_fd2,file_list2(i).name]);
                spectrum_2(:,i) = table2(find(fix(table2(:,1))==(startw),1):...
                    find(fix(table2(:,1))==(endw),1,'last'),2);
                
                if i==1
                    wavelength2 = table2(find(fix(table2(:,1))==(startw),1):...
                        find(fix(table2(:,1))==(endw),1,'last'),1);
                end
            end
% end

%%%%%% find and interpolate missing points

%%% find starting number
first = 0;
kk = 1;  %change first data point number based on starting number
while kk~=first
   if  strcmpi ({['QEP016571_',num2str(kk,'%05d'),'.txt']},file_list11(1,1)) == 1
first = kk;
   else
       kk = kk+1;
   end
end

for index = first:first+567 % 290 time points (change based on your data time points) % change based on your total number of points
count1_temp = {['QEP016571_',num2str(index,'%05d'),'.txt']};
count1 (index-first+1,:) = count1_temp;  

count2_temp = {['QEP009772_',num2str(index,'%05d'),'.txt']};
count2 (index-first+1,:) = count2_temp;
end
for count = 1:length (file_list11)
missing1(:,count)  = strcmpi ((count1)',file_list11(count));
end

for count = 1:length (file_list22)
missing2(:,count)  = strcmpi ((count2)',file_list22(count));
end

spectrum_data = zeros(length(spectrum_1),length(count1));
spectrum_data2 = zeros(length(spectrum_2),length(count2));

    mm = 0;
for index = 1:561  % change based on your total number of points
    if missing1 (index,:) == 0
       missing11(mm+1) = index; 
       spectrum_data (:,index) = 0.5*(spectrum_1(:,index-mm-1)+spectrum_1(:,index-mm));
       mm = mm+1;
    else
    spectrum_data (:,index) = spectrum_1(:,index-mm);    
    end
end

    mm = 0;
for index = 1:561 % change based on your total number of points
    if missing2 (index,:) == 0
       missing22(mm+1) = index; 
       spectrum_data2 (:,index) = 0.5*(spectrum_2(:,index-mm-1)+spectrum_2(:,index-mm));
       mm = mm+1;
    else
    spectrum_data2 (:,index) = spectrum_2(:,index-mm);    
    end
end


    spectrum_data(:,:,2) = spectrum_data2(1:length(spectrum_data),:);
    
    save([fd_name,sub_fd1,name1,'_raw_data.mat'],'spectrum_data','wavelength');   
end

%% deciding DPF
ua750=0.129; % obtain by ISS for forehead (change for arm)
ua830=0.129;   % check if its different based on tissue and wavelength
us750=12.60;
us830=11.40;
startw = 780;  % change wavelength interest region
endw = 950; % changed for water
HbO=(780.21*(ua750)-1545.39*(ua830))/(log(10)*(780.21*549.34-1050.43*1545.39));  %Baseline subtraction
Hb=(549.34*(ua830)-1050.43*(ua750))/(log(10)*(780.21*549.34-1050.43*1545.39));
a=exp((log(750)*log(us830/us750))/log(750/830)+log(us750));
b=log(us830/us750)/log(750/830);
w=startw*(10^(-3)):0.001:endw*(10^(-3));% wavelength range
us=a*((w*1000).^(-b));
us=us';
ext=xlsread('Extinction_coefficients_water'); % pre-estimated absorption values of HbO-Hb-CCO-fat-water at different wavelengths
ext=(ext(startw-649:endw-649,:)); % get pre-estimated data for wavelengths from 780 to 900
% uaa(:,:)=HbO*ext(:,2)+Hb*ext(:,3);
% uaa = uaa*log(10);
uaa(:,:)=HbO*ext(:,2)+Hb*ext(:,3)+ext(:,5); %changed for water
% uaa=uaa*log(10);
num_of_w=endw-startw+1;% number of wavelength that we covered
ua=uaa;
DPF=((sqrt(3*us))./(2*(sqrt(ua))));
% DPF = 6.26*ones(size(DPF)); % ask Xinlong for more details

%% extract data
extc=xlsread('Extinction_coefficients_water');
extc_1 = extc(startw-649:endw-649, 2:5); %changed for water
ext_inv = pinv(extc_1);
% ext_inv = (inv(transpose(extc_1)*extc_1))*transpose(extc_1);
% HbO=extc(startw-649:endw-649,2);
% Hb=extc(startw-649:endw-649,3);
% CCO=extc(startw-649:endw-649,4);
% H2O=extc(startw-649:endw-649,5);


% figure
% plot(ext(:,1), extc_1)

%% loop fitting

l = [1/5,1/5,1/5,1/5,1/5];

start_point = 2;

for i_channel = 1:size(spectrum_data,3)
    
    current_temp_data = spectrum_data(:,:,i_channel);
    % baseline=(current_temp_data(:,length(filename_id{1})).*...
    %     current_temp_data(:,length(filename_id{1}))).^(1/2);
    baseline=(mean(current_temp_data(:,start_point-1:start_point),2).*...
        mean(current_temp_data(:,start_point-1:start_point),2)).^(1/2);
    
    baseline = conv(baseline,l,'same');
    sigma1 = std(current_temp_data(:,start_point-1:start_point),0,2);
    
    
    for jj=start_point:size(spectrum_data,2) %for loop for all data sets apart from baselines
        
        post = current_temp_data(:,jj); % real data starts from the fourth column, first: wavelength, second and third: baseline
        post = conv(post,l,'same');
        
%         sigma = std(current_temp_data(:,jj:jj+2),0,2);
%         post = mean(current_temp_data(:,start_point+1:end,:,:),2);
        
        dOD(:,2)=log10(baseline)-log10(post);% change of signal
        
%         dOD_err = (((sigma./(post.*log(10))).^2)+(sigma1./(baseline.*log(10))).^2).^(1/2);
%         dOD_err = 0.0029;
%         dOD(:,3) = dOD_err;
        
        dOD(:,1)=wavelength(:,1);
        wavelength_modify=fix(dOD(:,1));% find which row is the starting wavelength
        [a,b]=find(wavelength_modify(:,1)==startw);
        [c,d]=find(wavelength_modify(:,1)==endw);
        x=min(a);% find the earliest row that reach starting wavelength
        y=max(c);% find the largest row that reach starting wavelength
        dOD_selected=dOD(x:c,:);% select the range we want, now just selected the range, we still need to select integer wavelengths in the next step
        dOD_selected(:,1)=fix(dOD_selected(:,1));% fix and select integer wavelegnths
        
        
        %% select only one data from multiple wavelengths reading like select one as 780nm from 780.14nm 780.75nm  780.98nm.
        for i=1:num_of_w
            [a,b]=find(dOD_selected(:,1)==(i+startw-1));
            x=min(a);
            y=min(b);
            final_dOD(i,2)=dOD_selected(x,2);
%             final_dOD_err(i,1)=dOD_selected(x,3);
        end
        final_dOD(:,1)=startw:1:endw;
        fit_dOD=final_dOD(:,2);% the data for fitting have to be without waveleng colum, so make a separate fit_dOD for fitting
        all_dOD(:,jj)=fit_dOD;
        
      Conc(:,jj)=((1./r).*ext_inv)*(fit_dOD./DPF);  
%         fit_dOD_err(:,jj)=final_dOD_err;

% for short and long separation (check the QEPro number and channel number) 
% r1 = 1;
% r2 = 3;
% if i_channel == 1
%         Conc(:,jj)=((1./r1).*ext_inv)*(fit_dOD./DPF);
% else
%     Conc(:,jj)=((1./r2).*ext_inv)*(fit_dOD./DPF);
% end
%         Conc_err(:,jj) = (1./r).*sqrt(((ext_inv./DPF').^2)*(fit_dOD_err(:,jj)).^2);
%         figure
%         plot(1:jj,10e6.*(1./r).*sqrt(((ext_inv./DPF').^2)*(fit_dOD_err(:,1:jj)).^2))
%         legend('std of \DeltaHbO (\muMolar)','std of \DeltaHHb (\muMolar)','std of \DeltaCCO (\muMolar)')
%         figure
%         plot(1:size(ext_inv,2),fit_dOD_err(:,1:jj))
%         hold on
%         figure
%         plot(1:size(ext_inv,2),ext_inv(:,:)./DPF')
%         hold on
%         plot(1:size(ext_inv,2),ext_inv(3,:)./DPF')
%         
       
%         Xinlong's way
%         Y0=[1e-6; 1e-6; 1e-6]; % fit for the least-square
%         [Conc(:,jj) , fmin4(:,jj)]=fminsearch(@oxyy,Y0,optimset('TolX',1e-20),r(i_channel),DPF,HbO,Hb,CCO,fit_dOD);%optimset('TolX',1e-8)
%         
%         ccco=Conc(3,jj);
%         chbo=Conc(1,jj);
%         chb=Conc(2,jj);
%         fitted_line=(ccco*CCO+chbo*HbO+chb*Hb).*r(i_channel).*DPF;
        % figure() % plot fitted line with original data
        % plot(final_dOD(:,1),fit_dOD,'.-','color','k');
        % hold on
        % plot(final_dOD(:,1),fitted_line,'.-','color','r');
        % legend data fit
        % xlabel( 'wavelength (nm)' );
        % ylabel( 'change of optical density' ); % Y4 and fmin4 is the fitted result of concentration of CCO HBO and HB and least square value
    end
    
    %% plotting
    
    recovery=(1:size(spectrum_data,2)-start_point+2)/40;

        chbo=Conc(1,jj);
        chb=Conc(2,jj);
        ccco=Conc(3,jj);
        ch2o=Conc(4,jj);
        
        HbT(1,jj)=Conc(1,jj)+Conc(2,jj); % HbT
        H2O_free(6,jj)=Conc(4,jj)-HbT(1,jj); % H2O (free)
        
  Conc(1:3,:)=Conc(1:3,:).*(10^6);  


    
    figure(1)
    subplot(411)
    plot(recovery,Conc(1,start_point-1:end)./max(abs(Conc(1,start_point-1:end))),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('\DeltaHbO (\muMolar)','FontSize',16 );
%     ylim([-1 1]);
    title('HbO','FontSize',18);
    hold on
    legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
    subplot(412)
    plot(recovery,Conc(2,start_point-1:end)./max(abs(Conc(2,start_point-1:end))),cl_list{i_channel});
    hold on
%      xlabel( 'Time (min)','FontSize',12 );
    ylabel( '\DeltaHb (\muMolar)','FontSize',16 );
%     ylim([-1 1]);
    title('Hb','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
    set(gca,'FontSize',16 );
    subplot(413)
    plot(recovery,Conc(3,start_point-1:end)./max(abs(Conc(3,start_point-1:end))),cl_list{i_channel});
    hold on
    xlabel( 'Time (min)','FontSize',16 );
    ylabel('\DeltaCCO (\muMolar)','FontSize',16 );
%     ylim([-1 1]);
    title('CCO','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
set(gca,'FontSize',16 );
 subplot(414)
    plot(recovery,Conc(4,start_point-1:end)./max(abs(Conc(4,start_point-1:end))),cl_list{i_channel});
    hold on
    xlabel( 'Time (min)','FontSize',16 );
    ylabel('\DeltaH2O (AU)','FontSize',16 );
%     ylim([-1 1]);
    title('Water','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
set(gca,'FontSize',16 );
    Conc_all(:,:,i_channel) = Conc(:,start_point-1:end);
    save([fd_name,sub_fd1,name1,'_HbO_CCO_water.mat'],'Conc_all');
end


%% Filtering
fs = 1/1.5;
for ii = 1:4
    for m = 1:2
[pxx,f] = pwelch(Conc_all(ii,:,m),[],[],round(1024*fs),fs);

pxx_all(:,ii,m) = pxx;
f_all(:,ii,m) = f;
    end
end

f1 = [0.004,0.006];

% n=3000;
a=[0 0 1 1];
% w=[1 1];

% Endogenic
% f2=[0 f1(1) f1(2) fs/2]/(fs/2);
% b_Endo=firpm(n,f2_Endo,a);
[b,d] = butter(3,0.005/(fs/2),'high');
[h,q]=freqz(b,d);

for jj = 1:4
    for m = 1:2
h1(jj,:,m) = filter(b,d,Conc_all(jj,:,m));
[pxx_1,f_1] = pwelch(h1(jj,:,m),[],[],round(1024*fs),fs);

pxx_1_all(:,jj,m) = pxx_1;
f_1_all(:,jj,m) = f_1;
    end
end
%     save([fd_name,sub_fd1,name1,'_HbO_CCO.mat'],'Conc_all');
% figure
% for i_channel = 1:size(spectrum_data,3)
% 
% subplot(411)
%     plot(recovery,h1(1,start_point-1:end,i_channel)./max(abs(h1(1,start_point-1:end,i_channel))),cl_list{i_channel});
% %     xlabel( 'Time (min)','FontSize',12 );
%     ylabel('\DeltaHbO (\muMolar)','FontSize',16 );
% %     ylim([-1 1]);
%     title('HbO','FontSize',18);
%     hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
%     set(gca,'FontSize',16 );
%     subplot(412)
%     plot(recovery,h1(2,start_point-1:end,i_channel)./max(abs(h1(2,start_point-1:end,i_channel))),cl_list{i_channel});
%     hold on
% %      xlabel( 'Time (min)','FontSize',12 );
%     ylabel( '\DeltaHb (\muMolar)','FontSize',16 );
% %     ylim([-1 1]);
%     title('Hb','FontSize',18);
% %     legend(legend_list,'Location','northwest','FontSize',12)
%     set(gca,'FontSize',16 );
%     subplot(413)
%     plot(recovery,h1(3,start_point-1:end,i_channel)./max(abs(h1(3,start_point-1:end,i_channel))),cl_list{i_channel});
%     hold on
%     xlabel( 'Time (min)','FontSize',16 );
%     ylabel('\DeltaCCO (\muMolar)','FontSize',16 );
% %     ylim([-1 1]);
%     title('CCO','FontSize',18);
% %     legend(legend_list,'Location','northwest','FontSize',12)
% set(gca,'FontSize',16 );
% subplot(414)
%     plot(recovery,h1(4,start_point-1:end,i_channel)./max(abs(h1(4,start_point-1:end,i_channel))),cl_list{i_channel});
% %     xlabel( 'Time (min)','FontSize',12 );
%     ylabel('\DeltaH2O (AU)','FontSize',16 );
% %     ylim([-1 1]);
%     title('Water','FontSize',18);
%     hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
%     set(gca,'FontSize',16 );
% end

%% smoothed plots

h2(:,:,:) = movmean (h1(:,:,:),5,2);
h2(5,:,:)=h2(1,:,:)+h2(2,:,:); % HbT
h2(6,:,:)=(h2(4,:,:)./max(h2(4,:,:)))-h2(5,:,:)./max(h2(5,:,:)); % H2O (free)

figure
for i_channel = 1:size(spectrum_data,3)

subplot(321)
    plot(recovery,h2(1,start_point-1:end,i_channel),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('\Delta[HbO] (\muM)','FontSize',16 );
%     ylim([-1 1]);
    title('HbO','FontSize',18);
    hold on
    legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
    subplot(323)
    plot(recovery,h2(2,start_point-1:end,i_channel),cl_list{i_channel});
    hold on
%      xlabel( 'Time (min)','FontSize',12 );
    ylabel( '\Delta[Hb] (\muM)','FontSize',16 );
%     ylim([-1 1]);
    title('Hb','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
    set(gca,'FontSize',16 );
    subplot(325)
    plot(recovery,h2(3,start_point-1:end,i_channel),cl_list{i_channel});
    hold on
    xlabel( 'Time (min)','FontSize',16 );
    ylabel('\Delta[CCO] (\muM)','FontSize',16 );
%     ylim([-1 1]);
    title('CCO','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
set(gca,'FontSize',16 );
subplot(322)
    plot(recovery,h2(5,start_point-1:end,i_channel),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('\Delta[HbT] (\muM) ','FontSize',16 );
%     ylim([-1 1]);
    title('HbT','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
subplot(324)
    plot(recovery,100*h2(4,start_point-1:end,i_channel),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('\DeltaH2O (%)','FontSize',16 );
%     ylim([-1 1]);
    title('Water','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
    subplot(326)
    plot(recovery,h2(6,start_point-1:end,i_channel),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('Norm. \DeltaH2O (Free)','FontSize',16 );
%     ylim([-1 1]);
    title('Water (Free)','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
end

figure
for i_channel = 1:size(spectrum_data,3)

subplot(321)
    plot(recovery,h2(1,start_point-1:end,i_channel)./max(abs(h2(1,start_point-1:end,i_channel))),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('Normalized \DeltaHbO','FontSize',16 );
%     ylim([-1 1]);
    title('HbO','FontSize',18);
    hold on
    legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
    subplot(323)
    plot(recovery,h2(2,start_point-1:end,i_channel)./max(abs(h2(2,start_point-1:end,i_channel))),cl_list{i_channel});
    hold on
%      xlabel( 'Time (min)','FontSize',12 );
    ylabel( 'Normalized \DeltaHb','FontSize',16 );
%     ylim([-1 1]);
    title('Hb','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
    set(gca,'FontSize',16 );
    subplot(325)
    plot(recovery,h2(3,start_point-1:end,i_channel)./max(abs(h2(3,start_point-1:end,i_channel))),cl_list{i_channel});
    hold on
    xlabel( 'Time (min)','FontSize',16 );
    ylabel('Normalized \DeltaCCO','FontSize',16 );
%     ylim([-1 1]);
    title('CCO','FontSize',18);
%     legend(legend_list,'Location','northwest','FontSize',12)
set(gca,'FontSize',16 );
subplot(322)
    plot(recovery,h2(5,start_point-1:end,i_channel)./max(abs(h2(5,start_point-1:end,i_channel))),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('Normalized \DeltaHbT ','FontSize',16 );
%     ylim([-1 1]);
    title('HbT','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
subplot(324)
    plot(recovery,h2(4,start_point-1:end,i_channel)./max(abs(h2(4,start_point-1:end,i_channel))),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('Normalized \DeltaH2O','FontSize',16 );
%     ylim([-1 1]);
    title('Water','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
    subplot(326)
    plot(recovery,h2(6,start_point-1:end,i_channel)./max(abs(h2(6,start_point-1:end,i_channel))),cl_list{i_channel});
%     xlabel( 'Time (min)','FontSize',12 );
    ylabel('Normalized \DeltaH2O (Free)','FontSize',16 );
%     ylim([-1 1]);
    title('Water (Free)','FontSize',18);
    hold on
%     legend(legend_list,'Location','northwest','FontSize',16)
    set(gca,'FontSize',16 );
end