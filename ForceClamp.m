%#######################################################################
%
%                       * Force Clamp Program *
%
%          M-File which reads data files (either new .dat or old .fcp
%     files) from muscle force clamp experiments with three steps in
%     muscle lengths.
%
%          Users pick the end of the muscle clamp and the data 20 ms
%     before the end to the picked end is averaged.  The muscle lengths
%     are converted to fractions of the total muscle length.  Muscle
%     velocities are calculated from the muscle length time histories.
%
%          The muscle forces are corrected for the passive muscle
%     forces using two different methods.  A Hill's equation model is
%     used to fit the muscle velocities and forces.
%
%          Plots of the muscle velocity and muscle power with muscle
%     are plotted for review before saving the data and fit parameters
%     to a MS-Excel spreadsheet.
%
%          The user selects whether to create a new spreadsheet file or
%     append to an user selected existing spreadsheet file.
%
%     NOTES:  1.  If the user selects an existing output MS-Excel
%             spreadsheet, the spreadsheet cannot be open in another
%             program (e.g. MS-Excel, text editor, etc.) while using
%             this program.
%
%             2.  M-files get_tens.m, rd_dat.m, rd_fcp.m, and val2idx.m
%             must be in the current path or directory.
%
%     29-Mar-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Clear Workspace
%
clc;
clear variables;
close all;
fclose all;
%
% Additional Prompts When Digitizing?
%
% lrn_steps = true;       % Prompt for each step when digitizing
lrn_steps = false;      % Just prompt for the first step when digitizing
%
% Filter Parameters
%
filterOrder = 2;        % Desired filter order/2 to account for filtfilt processing
cutoff = 100;           % 100 Hz cutoff
%
% Get Input Data .DAT or .FCP File(s)
%
[fnam,pnam,fidx] = uigetfile({'*.dat','New Force Clamp files'; ...
    '*.fcp','Old Force Clamp files';'*.*','All files (*.*)'}, ...
    'Please Select Force Clamp Files for Analysis','MultiSelect', 'on');
%
if fidx==0              % User hit "Cancel" button
  return;
end
%
% Read Files Based on the File Name Extension
%
if iscell(fnam)
%
  fnam = fnam(:);
%
% Read Multiple New Force Clamp Files (*.DAT)
%
  if all(endsWith(fnam,'.dat','IgnoreCase',true))
%
    new = true;         % New force clamp files
    nv = size(fnam,1);  % Number of valid force clamp runs
%
    srate = zeros(nv,1);     % Sampling rate in Hz
    flen = zeros(nv,1); % Muscle fiber length in mm
    slen = zeros(nv,1); % Muscle sarcomere length in um
    fdia = zeros(nv,1); % Muscle fiber diameter in mm
    t = cell(1,nv);     % Time in s
    mlen = cell(1,nv);  % Muscle length in mm
    mfrc = cell(1,nv);  % Muscle force in mN
    pass = cell(1,nv);  % Passive muscle force
    act = cell(1,nv);   % Active muscle force
    tstep = zeros(3,nv);     % Initial time of steps in s
    vstep = zeros(3,nv);     % Option value of steps
    tlens = zeros(nv,1);     % Initial time of length sample in s
    tlenh = zeros(nv,1);     % Hold (end) time of length sample in s
%
    for k = 1:nv
       [srate(k),flen(k),fdia(k),t{k},mlen{k},~,mfrc{k},~,pass{k}, ...
                     act{k},tstep(:,k),vstep(:,k),tlens(k),tlenh(k), ...
                     slen(k)] = rd_dat(fullfile(pnam,fnam{k}));
    end
%
% Sort Files by First Option Value
%
    [~,ids] = sort(vstep(1,:),'descend');
    fnam = fnam(ids);
    srate = srate(ids); 
    flen = flen(ids)*1000;             % Convert mm to um
    slen = slen(ids);
    fdia = fdia(ids); 
    fdia1 = fdia*1000;
    fdia2 = fdia*1000;
    fcsa = pi*fdia1.*fdia2/4;
    fvol = fcsa.*flen/1000.^3;         % Muscle fiber volume in mm^3
    t = t(:,ids);
    mlen = mlen(:,ids);
    mfrc = mfrc(:,ids);
    pass = cellfun(@mean,pass(:,ids))';
    act = cellfun(@mean,act(:,ids))';
    tstep = tstep(:,ids);
    vstep = 100*vstep(:,ids);          % Convert to percent
    tlens = tlens(ids);
    tlenh = tlenh(ids);
%
  else
    error([' *** ERROR in ForceClamp:  Multiple input files only', ...
           ' supported for new force clamp files ( *.dat)!']);
  end
%
else
%
% Read an Old Force Clamp File (*.FCP)
%
  if endsWith(fnam,'.fcp','IgnoreCase',true)
    new = false;        % Old force clamp file
    [t,mlen,mfrc,mstrain,mstress,n,rnums,amp,flen,fdia1,fdia2, ...
                                    fcsa] = rd_fcp(fullfile(pnam,fnam));
%
    fvol = fcsa*flen/1000;             % Volume in mm^3
    fcsa = fcsa*1000*1000;             % Convert from mm^2 to um^2
    slen = NaN;         % Sarcomere length not recorded in file

%
    vstep(3,:) = amp'+100;             % Get last (3rd) step amplitudes
    vstep(1,:) = 2*amp'+270;           % Get first (1st) step amplitudes
    vstep(2,:) = vstep(1,:)/2-10;      % Get second (2nd) step amplitudes
    idx = vstep(1,:)==90;              % Trap for first step = 90%
    vstep(2,idx) = 40;                 % Set corresponding second step to 40%
%
    [srate,pass,act,idv,nv] = get_tens(n,rnums,t);
%
% Valid Runs Only
%
    t = t(idv);
    mlen = mlen(idv);
    mfrc = mfrc(idv);
    mstrain = mstrain(idv);
    mstress = mstress(idv);
    vstep = vstep(:,idv);
%
    mlen = cellfun(@(x) rdivide(x,1000),mlen,'UniformOutput',false); % Convert um to mm
%
% Read a New Force Clamp File (*.DAT)
%
  elseif endsWith(fnam,'.dat','IgnoreCase',true)
    new = true;         % New force clamp file
    [srate,flen,fdia,t,mlen,~,mfrc,~,pass,act,tstep,vstep,tlens, ...
                              tlenh,slen] = rd_dat(fullfile(pnam,fnam));
    nv = 1;             % Number of valid force clamp runs
    t = {t};
    mlen = {mlen};
    mfrc = {mfrc};
    pass = {pass};
    act = {act};
    flen = flen*1000;   % Convert mm to um
    fdia1 = fdia*1000;
    fdia2 = fdia*1000;
    fcsa = pi*fdia1*fdia2/4;
    fvol = fcsa*flen/1000.^3;
    vstep = 100*vstep;  % Convert to percent
  else
    error([' *** ERROR in ForceClamp:  Input file name extension', ...
           ' not recognized (*.dat or *.fcp)!']);
  end
%
end
%
% Filter Data
%
if all(srate==srate(1))
%
% Setup Filter
%
  nyquist = srate(1)/2; 
  Wn = cutoff/nyquist;
  [b,a] = butter(filterOrder,Wn);
%
% Apply Filter
%
  mlen_filt = cellfun(@(x) filtfilt(b,a,x),mlen,'UniformOutput',false);
  mfrc_filt = cellfun(@(x) filtfilt(b,a,x),mfrc,'UniformOutput',false);
%
else
%
  mlen_filt = cell(1,nv);
  mfrc_filt = cell(1,nv);
%
  for k = 1:nv
%
% Setup Filter
%
     nyquist = srate(k)/2;
     Wn = cutoff/nyquist;
     [b,a] = butter(filterOrder,Wn);
%
% Apply Filter
%
     mlen_filt{k} = filtfilt(b,a,mlen{k});
     mfrc_filt{k} = filtfilt(b,a,mfrc{k});
%
  end
%
end
%
% Mean Passive and Maximum Active Muscle Forces
%
pass_avg = mean(pass);
%
act_max = max(act);
%
% Plot Muscle Force Time History and Get the Last 20 ms of Each Step 
%
mfrc_base = zeros(1,nv);               % Mean base muscle forces
mlen_base = zeros(1,nv);               % Mean base muscle lengths
idbase = zeros(2,nv);   % Indices to 25 ms base muscle force window
%
mvel = cell(1,nv);      % Normalized muscle velocities
%
mfrc_step = zeros(3,nv);               % Mean step muscle forces
mvel_step = zeros(3,nv);               % Mean step muscle velocities
idsteps = zeros(3,nv);  % Index to the start of step data
idstepe = zeros(3,nv);  % Index to the end of step data
%
hf = figure('Name','Force Time History','WindowState','maximized', ...
            'NumberTitle','off');
orient landscape;
%
for k = 1:nv
%
   plot(t{k},mfrc_filt{k},'b.-','LineWidth',1.0,'MarkerSize',7);
   xlabel('Time (s)','FontSize',12,'FontWeight','bold');
   ylabel('Muscle Force (mN)','FontSize',12,'FontWeight','bold');
   ttxt = {'Muscle Force Clamp Steps'; ['Step 1 = ', ...
           int2str(vstep(1,k)) '%, Step 2 = ' int2str(vstep(2,k)), ...
           '%, Step 3 = ' int2str(vstep(3,k)) '%']};
   title(ttxt,'FontSize',16,'FontWeight','bold');
   hold on;
   axlim = axis;
%
   if new
%
     idbase(2,k) = val2idx(t{k},tstep(1,k)-0.025);    % End of base muscle force window
     idbase(1,k) = val2idx(t{k},t{k}(idbase(2))-0.025);    % Start of base muscle force window
%
     ha = gca;
     ha.XLim = [tstep(1,k)-0.1 tlenh(k)];
     drawnow;
     axlim = axis;
%
     for ks = 1:3
        plot(repmat(tstep(ks,k),1,2),axlim(3:4),'g:', ...
             'Color',[0 0.7 0],'LineWidth',1);   % Start of each step
     end
%
     plot(repmat(tlens,1,2),axlim(3:4),'r:','LineWidth',1);     % Length sample
     plot(repmat(tlenh,1,2),axlim(3:4),'r:','LineWidth',1);     % Length hold
     axis(axlim);
%
   else
%
     idbase(1,k) = val2idx(t{k},0.005);     % Start of base muscle force window
     idbase(2,k) = val2idx(t{k},t{k}(idbase(1))+0.025);     % End of base muscle force window
%
   end
%
% Plot Window for Mean Base Muscle Force, and Mean Base Muscle Force
%
   idx = idbase(1,k):idbase(2,k);
   mfrc_base(k) = mean(mfrc_filt{k}(idx));
   mlen_base(k) = mean(mlen_filt{k}(idx));
%
   for kf = 1:2
      plot(repmat(t{k}(idbase(kf,k)),1,2),axlim(3:4),'k:', ...
           'LineWidth',1);
   end
   plot(axlim(1:2),repmat(mfrc_base(k),1,2),'k-','LineWidth',1);     % Mean base force
   axis(axlim);
%
% Calculate Normalized Muscle Velocities using Central Difference
%
   mvel{k} = mlen_filt{k}./mlen_base(k);    % Normalize muscle lengths
   mvel{k} = diff(mvel{k})./diff(t{k});     % Muscle velocities (ML/s)
   mvel{k} = [mvel{k}(1); (mvel{k}(1:end-1)+mvel{k}(2:end))/2; ...
              mvel{k}(end)];
%
% Get End of Each step
%
   for ks = 1:3         % Loop through steps
%
      hs = gobjects(0);
      he = gobjects(0);
%
      stxt = int2str(ks);
%
      kpck = true;
      while kpck
%
           delete([hs,he]);
%
           if ks==1
             msgtxt = ['Pick an End Time for Step ' stxt, ...
                       ' of the clamp.'];
             uiwait(msgbox({msgtxt},'non-modal'));
           elseif ks>1&&lrn_steps
             msgtxt = ['Pick an End Time for Step ' stxt, ...
                       ' of the clamp.'];
             uiwait(msgbox({msgtxt},'non-modal'));
           end
           figure(hf);
           [itim,~] = ginput(1);
%
           idstepe(ks,k) = val2idx(t{k},itim);   % From end time to index
           idsteps(ks,k) = val2idx(t{k},t{k}(idstepe(ks,k))-0.02);   % From start time to index
%
% Plot Times
%
           hs = plot([t{k}(idsteps(ks,k)) t{k}(idsteps(ks,k))], ...
                     axlim(3:4),'g-','Color',[0.2 1 0.2],'LineWidth',1);
           he = plot([t{k}(idstepe(ks,k)) t{k}(idstepe(ks,k))], ...
                     axlim(3:4),'g-','Color',[0 0.6 0],'LineWidth',1);
           axis(axlim);
           figure(hf);
           drawnow;
%
% Check End Time
%
           kpck = logical(2-menu({['End time for Step ' stxt, ...
                          ' OK?']; '(dark green vertical line)'}, ...
                          'No','Yes'));
%
      end
%
% Get Mean Muscle Forces and Velocities
%
      idx = idsteps(ks,k):idstepe(ks,k);
      mfrc_step(ks,k) = mean(mfrc_filt{k}(idx));
      mvel_step(ks,k) = mean(mvel{k}(idx));
      
   end
%
% Plot Muscle Force and Velocity
%
   hf2 = figure('Name','Step Mean Values','WindowState', ...
                'maximized','NumberTitle','off');
   orient landscape;
%
   ha1 = subplot(1,2,1);
   plot(t{k},mfrc_filt{k},'g.-','Color',[0 0.6 0],'LineWidth',1.0, ...
        'MarkerSize',7);
   xlabel('Time (s)','FontSize',12,'FontWeight','bold');
   ylabel('Muscle Force (mN)','FontSize',12,'FontWeight','bold');
   ha1.XLim = axlim(1:2);
   hold on;
%
   plot(t{k}([idbase(1,k),idbase(2,k)]),repmat(mfrc_base(k),1,2), ...
        'b-','LineWidth',3);
%
   for ks = 1:3
      plot(t{k}([idsteps(ks,k),idstepe(ks,k)]), ...
           repmat(mfrc_step(ks,k),1,2),'r-','LineWidth',3);
   end
%
   ha2 = subplot(1,2,2);
   plot(t{k},mvel{k},'r.-','LineWidth',1.0,'MarkerSize',7);
   xlabel('Time (s)','FontSize',12,'FontWeight','bold');
   ylabel('Muscle Velocity (ML/s)','FontSize',12, ...
          'FontWeight','bold');
   ha2.XLim = axlim(1:2);
   ha2.YLim(2) = 0.1;
   hold on;
%
   for ks = 1:3
      plot(t{k}([idsteps(ks,k),idstepe(ks,k)]), ...
           repmat(mvel_step(ks,k),1,2),'b-','Color',[0 0.6 0], ...
           'LineWidth',3);
   end
%
   sgtitle(ttxt,'FontSize',16,'FontWeight','bold');
%   
   pause;
   close(hf2);
%
   clf(hf);
%
end
%
close(hf);
%
% Normalize Muscle Forces and Correct for Passive Muscle Forces
%
Fo = max(mfrc_base);
mfrc_rel = mfrc_step(:)./Fo;
mfrc_num = mfrc_rel*act_max;
mfrc_den = act_max-pass_avg;
mfrc_pass1 = mfrc_num./mfrc_den;
mfrc_pass2 = (mfrc_num-pass_avg)./mfrc_den; % Per Mark's spreadsheet
%
% Absolute of Normalized Muscle Velocities
%
mvel_abs = abs(mvel_step(:));
Vo = max(mvel_abs);
% mvel_rel = mvel_norm./Vo;
%
% Plot Muscle Velocities versus Muscle Forces
%
figure('Name','Step Mean Values','WindowState','maximized', ...
       'NumberTitle','off');
orient landscape;
%
plot(mfrc_pass1,mvel_abs,'go','LineWidth',1,'MarkerSize',7);
hold on;
plot(mfrc_pass2,mvel_abs,'rx','LineWidth',1,'MarkerSize',7);
%
legend({'Passive force corrected','Per Mark''s spreadsheet'});
%
xlabel('Normalized Muscle Force','FontSize',12,'FontWeight','bold');
ylabel('Muscle Velocity (ML/s)','FontSize',12,'FontWeight','bold');
title('Muscle Force-Velocity Curve','FontSize',16,'FontWeight','bold');
%
pause;
close all;
%
% Fit Passive Muscle Force Corrected Data (pass1)
%
param = [1; 0.1];
%
hillfun = @(hill_pars,xdata)((hill_pars(2)*(1+hill_pars(1)))./ ...
                            (xdata+hill_pars(1)))-hill_pars(2);
[hill_pars1,~,res1] = lsqcurvefit(hillfun,param,mfrc_pass1,mvel_abs);
%
yfit1 = hillfun(hill_pars1,mfrc_pass1);
R2_1 = 1-sum((mvel_abs-yfit1).^2)/sum((mvel_abs-mean(mvel_abs)).^2);
%
Hill_a1 = hill_pars1(1);
Hill_b1 = hill_pars1(2);
%
vmax1 = Hill_b1/Hill_a1;
pow1 = mfrc_pass1.*mvel_abs;
%
% Fit Passive Muscle Force Corrected Data Per Mark's Spreadsheet (pass2)
%
[hill_pars2,~,res2] = lsqcurvefit(hillfun,param,mfrc_pass2,mvel_abs);
%
yfit2 = hillfun(hill_pars2,mfrc_pass2);
R2_2 = 1-sum((mvel_abs-yfit2).^2)/sum((mvel_abs-mean(mvel_abs)).^2);
%
Hill_a2 = hill_pars2(1);
Hill_b2 = hill_pars2(2);
%
vmax2 = Hill_b2/Hill_a2;
pow2 = mfrc_pass2.*mvel_abs;
%
% Calculate Muscle Velocities and Power Using Model Fits
%
velfit1 = hillfun(hill_pars1,mfrc_pass1);
velfit2 = hillfun(hill_pars2,mfrc_pass2);
%
powfit1 = velfit1.*mfrc_pass1;
powfit2 = velfit2.*mfrc_pass2;
%
mfrc_pass1_rng = sort([-mfrc_pass1; mfrc_pass1]);
mfrc_pass1_rng = mfrc_pass1_rng(mfrc_pass1_rng>-1);   % Plot range
%
mfrc_pass2_rng = sort([-mfrc_pass2; mfrc_pass2]);
mfrc_pass2_rng = mfrc_pass2_rng(mfrc_pass2_rng>-1);   % Plot range
%
velfitplt1 = hillfun(hill_pars1,mfrc_pass1_rng); % Velocities for plot
velfitplt2 = hillfun(hill_pars2,mfrc_pass2_rng); % Velocities for plot
%
% Plot Muscle Velocities and Power versus Muscle Forces for Negative
% Range
%
figure('Name','Muscle Velocities and Power','WindowState', ...
       'maximized','NumberTitle','off');
orient landscape;
%
ha = axes;
yyaxis left;
ha.YColor = 'k';        % Left Y-axis color is black
plot(mfrc_pass1,mvel_abs,'go','LineWidth',1,'MarkerSize',7);
hold on;
plot(mfrc_pass2,mvel_abs,'rx','LineWidth',1,'MarkerSize',7);
plot(0,vmax1,'gs','MarkerFaceColor','g','LineWidth',1,'MarkerSize',7);
plot(0,vmax2,'rs','MarkerFaceColor','r','LineWidth',1,'MarkerSize',7);
plot(mfrc_pass1_rng,velfitplt1,'g-','LineWidth',1);
plot(mfrc_pass2_rng,velfitplt2,'r-','LineWidth',1);
%
xlabel('Normalized Muscle Force','FontSize',12,'FontWeight','bold');
ylabel('Muscle Velocity (ML/s)','FontSize',12,'FontWeight','bold');
%
ha = gca;
ha.YLim = [-0.1 1.07*max([mvel_abs; vmax1; vmx2])];
%
[mfrc_pass1s,ids1] = sort(mfrc_pass1); % Sort muscle forces for plot
[mfrc_pass2s,ids2] = sort(mfrc_pass2); % Sort muscle forces for plot
%
yyaxis right;
ha.YColor = 'k';        % Right Y-axis color is black
plot(mfrc_pass1s,pow1(ids1),'bv','LineWidth',1,'MarkerSize',7);
hold on;
plot(mfrc_pass2s,pow2(ids2),'m^','LineWidth',1,'MarkerSize',7);
plot(mfrc_pass1s,powfit1(ids1),'b-','LineWidth',1);
plot(mfrc_pass2s,powfit2(ids2),'m-','LineWidth',1);
%
legend({'Passive force corrected';'Per Mark''s spreadsheet'; ...
        'Maximum velocity passive'; 'Maximum velocity per Mark''s'; ...
        'Passive force corrected fit';'Per Mark''s spreadsheet fit'; ...
        'Power passive force corrected'; ...
        'Power per Mark''s spreadsheet'; ...
        'Power passive force corrected fit'; ...
        'Power per Mark''s spreadsheet fit'});
%
ylabel('Normalized Power (ML/s)','FontSize',12,'FontWeight','bold');
%
pause;
close all;
%
% Calculate Additional Output Parameters
%
Vo_norm1 = Vo/vmax1;
Vo_norm2 = Vo/vmax2;
%
Vopt1 = Hill_b1*sqrt((1/Hill_a1)+1)-Hill_b1;
Vopt2 = Hill_b2*sqrt((1/Hill_a2)+1)-Hill_b2;
%
Topt1 = sqrt(Hill_a1.^2+Hill_a1)-Hill_a1;
Topt2 = sqrt(Hill_a2.^2+Hill_a2)-Hill_a2;
%
pmax1 = Vopt1*Topt1;
pmax2 = Vopt2*Topt2;
%
% Combine Parameters and Results
%
% pow_all = [pow1 pow2];
Hill_a_all = [Hill_a1; Hill_a2];
Hill_b_all = [Hill_b1; Hill_b2];
R2_all = [R2_1; R2_2];
vmax_all = [vmax1; vmax2];
Vo_norm_all = [Vo_norm1; Vo_norm2];
pmax_all = [pmax1; pmax2];
Vopt_all = [Vopt1; Vopt2];
Topt_all = [Topt1; Topt2];
%
% Create Array of Data and Results
%
dat1 = [vstep(:) mfrc_rel mfrc_pass2 mvel_abs pow1 velfit1 powfit1];
nr = size(dat1,1);
%
mdata1 = NaN(nr-nv,1);  % NaNs (missing data) to fill arrays
mdata2 = NaN(nr-2,1);   % NaNs (missing data) to fill arrays
%
dat2 = [[pass; mdata1] [pass_avg; pass_avg; mdata2] [act; mdata1]];
%
mdata3 = NaN(nr-2,5);   % NaNs (missing data) to fill arrays
%
dat3 = [[act_max; act_max] Hill_b_all Hill_a_all R2_all vmax_all];
dat3 = [dat3; mdata3];
%
dat4 = [Vo_norm_all pmax_all Vopt_all Topt_all [Vo; Vo]];
dat4 = [dat4; mdata3];
%
dat = [dat1 dat2 dat3 dat4];           % Array of data and results
%
% Include Muscle Fiber Geometry
%
dat0 = [flen fdia1 fdia2 fcsa fvol slen];
nr0 = size(dat0,1);
mdata0 = NaN(nr-nr0,6);
dat0 = [dat0; mdata0];
%
dat = [dat0 dat];       % Array of data and results with fiber geometry
%
% Include Filename
%
nf = size(fnam,1);
%
if iscell(fnam)
  datf = [fnam; repmat({' '},nr-nf,1)];
  xlsnam = fnam{1};     % New output spreadsheet file name
else
  datf = [cellstr(fnam); repmat({' '},nr-nf,1)];
  xlsnam = fnam;        % New output spreadsheet file name
end
%
% Create Table of Data and Results
%
hdr0 = {'Fiber length (um)','Top width (um)','Side width (um)', ...
        'Cross-Sectional Area (um2)','Volume (mm3)', ...
        'Sarcomere Length (um)'};

hdr1 = {'Target Stress (percent Tmax)', ...
        'Actual Stress (fraction Tmax)', ...
        'Actual Stress - Passive (fraction Tmax)', ...
        'Velocity (ML/s)','Power (ML/s x fraction Tmax)'};
hdr2 = {'Fit Velocity (ML/s)','Fit Power (ML/s x fraction Tmax)', ...
        'Passive Tension (mM/mm2)', ...
        'Average Passive Tension (mM/mm2)'};
hdr3 = {'Force Clamp Tension (mM/mm2)','Force Clamp Tmax (mM/mm2)', ...
        'b (ML/s)','a','Rsquare','Vmax (ML/s)','Vo/Vmax', ...
        'Pmax (ML/s x fraction Tmax)','Vopt (ML/s)','Topt/Tmax'};

hdrs = [hdr0 hdr1 hdr2 hdr3 'Vo (ML/s)'];
%
t1 = array2table(datf,'VariableNames',{'Filename'});
t2 = array2table(dat,'VariableNames',hdrs);
%
t = [t1 t2];
%
% Write Table to MS-Excel Spreadsheet
%
newfile = logical(2-menu({'Create a new spreadsheet file or'; ...
                  'append to an existing file?'},'New file', ...
                  'Existing file'));
%
if newfile
%
  idot = strfind(xlsnam,'.');          % Find "dots" in file name
  xlsnam = xlsnam(1:idot(end));        % Remove file extension
  xlsnam = [xlsnam 'xlsx'];
  xlsnam = fullfile(pnam,xlsnam);      % New output spreadsheet file name
%
  repfile = true;
  if exist(xlsnam,'file')
    repfile = logical(2-menu('New spreadsheet file already exists!', ...
                      'Replace file','Append to the existing file'));
  end
%
  if repfile
    writetable(t,xlsnam,'Sheet','Summary','WriteMode','replacefile');
  else
    writetable(t,xlsnam,'Sheet','Summary','WriteVariableNames', ...
               false,'WriteMode','append');
  end
%
else                    % Append to an existing file
%
% Get Output MS-Excel Spreadsheet File
%
  fidx = 0;
  while fidx==0
       [fnam,pnam,fidx] = uigetfile({'*.xls;*.xlsx', ...
          'MS-Excel files'},'Please Select an Output Spreadsheet File');
  end
%
  xlsnam = fullfile(pnam,fnam);
  writetable(t,xlsnam,'Sheet','Summary','WriteVariableNames', ...
             false,'WriteMode','append');
%
end
%
return