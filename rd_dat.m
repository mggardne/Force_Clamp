function [srate,flen,fdia,t,mlen_in,mlen_out,mfrc_in,mfrc_out, ...
                   pass,act,tstep,vstep,tlens,tlenh,slen] = rd_dat(fnam)
%RD_DAT    Opens a muscle clamp .DAT output file and reads the sampling
%          rate, initial muscle fiber length, muscle fiber diameter,
%          time, muscle lengths, and muscle forces in the file.
%
%          [SRATE,FLEN,FDIA,T,MLEN_IN,MLEN_OUT,MFRC_IN,MFRC_OUT] =
%          RD_FCP(FNAM) opens a muscle clamp output file (*.DAT), FNAM,
%          and reads the sampling rate, SRATE, in Hz, initial muscle
%          fiber length, FLEN, in mm, muscle fiber diameter, FDIA, in
%          mm, time, T, in s, muscle length, MLEN_IN, in mm, muscle
%          length out, MLEN_OUT, in mm, muscle force in, MFRC_IN, in mN,
%          and muscle force out, MFRC_OUT, in mN.
%
%          [SRATE,FLEN,FDIA,T,MLEN_IN,MLEN_OUT,MFRC_IN,MFRC_OUT,PASS,
%          ACT,TSTEP,VSTEP,TLENS,TLENH] = RD_FCP(FNAM) returns the
%          passive, PASS, and active, ACT, muscle forces, initial time
%          of the three steps, TSTEP, the option value of the steps,
%          VSTEP, the initial time, TLENS, and hold time, TLENH, of the
%          length sample.
%
%          [SRATE,FLEN,FDIA,T,MLEN_IN,MLEN_OUT,MFRC_IN,MFRC_OUT,PASS,
%          ACT,TSTEP,VSTEP,TLENS,TLENH,SLEN] = RD_FCP(FNAM) returns the
%          initial sarcomere length, SLEN, in um.
%
%          NOTES:  1.  The muscle clamp output file (*.DAT) must be in
%                  a standard format.
%
%                  2.  Only one single run is in an output file.
%
%         06-Mar-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Input
%
if (nargin<1)
  error(' *** ERROR in rd_dat:  No input file name!');
end
%
% Check Input File Name
%
idot = strfind(fnam,'.');
idot = idot(end);
if ~strcmpi(fnam(idot:end),'.dat')
  error(' *** ERROR in rd_dat:  Input file extension must be .DAT!');
end
%
% Open Input File
%
fid = fopen(fnam,'r');
%
% Read the Data from the File
%
umm = 'mm';             % Units of mm
%
hdr0 = 'Time (ms)	Control Function	Options';
hdr1 = 'Time (ms)  Lin (mm)  Lout (mm)  Fin (mN)  Fout (mN)';
%
frmt = '%f %f %f %f %f %f %f %f %f %f %f';       % Format of data
%
while ~feof(fid)
     lin = fgets(fid);
%
     if ~feof(fid)
%
       if startsWith(lin,'*** ')
%
% Get Setup Parameters
%
         if startsWith(lin,'*** Setup Parameters ***')
           lin = fgets(fid);
           while ~startsWith(lin,'*** ')
%
                if startsWith(lin,'Fiber Length:')
                  idu = strfind(lin,umm);
                  flen = str2double(lin(14:idu-1));
                end
%
                if startsWith(lin,'Initial Sarcomere Length:')
                  idu = strfind(lin,'um');
                  slen = str2double(lin(26:idu-1));
                end
%
                if startsWith(lin,'Diameter:')
                  idu = strfind(lin,umm);
                  fdia = str2double(lin(10:idu-1));
                end
%
                lin = fgets(fid);
%
           end
%
         end
%
% Get Test Protocol Parameters
%
         if startsWith(lin,'*** Test Protocol Parameters ***')
           lin = fgets(fid);
           n = 0;       % Counter of control function timing
%
           while ~startsWith(lin,'*** ')
%
                if startsWith(lin,'A/D Sampling Rate:')
                  idu = strfind(lin,'Hz');
                  srate = str2double(lin(19:idu-1));
                end
%
                if startsWith(lin,hdr0)
                  lins = cell(200,1);
%
                  while ~startsWith(lin,'*** ')
                       n = n+1;
                       lin = fgets(fid);
                       lins{n} = lin;
                  end
%
                  n = n-1;
                  lins = lins(1:n);
                  lins = strtrim(lins);
                  tstep = zeros(3,1);  % Initial time of step
                  vstep = zeros(3,1);  % Option value of step
                  ks = 0;
%
                  for k = 1:n
                     astr = strsplit(lins{k});
%
                     if startsWith(astr{2},'Force-Sample')
                       if str2double(astr{3})==1
                         tpass = str2double(astr{1});
                         dtpass = str2double(astr{4});
                       else
                         tact = str2double(astr{1});
                         dtact = str2double(astr{4});
                       end
                     end
%
                     if startsWith(astr{2},'Force-Step')
                       ks = ks+1;
                       tstep(ks) = str2double(astr{1});
                       vstep(ks) = str2double(astr{3});
                     end
%
                     if startsWith(astr{2},'Length-Sample')
                       tlens = str2double(astr{1});
                     end
%
                     if startsWith(astr{2},'Length-Hold')
                       tlenh = str2double(astr{1});
                     end
%
                  end   % End of for k loop
%                
                end     % End of if line starts with hdr0
%
             if n<1
               lin = fgets(fid);
             end
%
           end          % End of while loop
%
         end            % End of if Test Protocol Parameters
%
% Get Time, Force, and Length Signals
%
         if startsWith(lin,'*** Force and Length Signals vs Time ***')
           lin = fgets(fid);
           if startsWith(lin,hdr1)
             data = textscan(fid,frmt,'Delimiter','space');
             t = data{:,1};
             mlen_in = data{:,2};
             mlen_out = data{:,3};
             mfrc_in = data{:,4};
             mfrc_out = data{:,5};
           end
%
         end
%         
       end
%
     end
end
%
% Get the Passive and Active Muscle Forces
%
dt = t-tpass;
[~,idp1] = min(abs(dt));
dt = t-(tpass+dtpass);
[~,idp2] = min(abs(dt));
%
pass = mfrc_in(idp1:idp2);
%
dt = t-tact;
[~,idp1] = min(abs(dt));
dt = t-(tact+dtact);
[~,idp2] = min(abs(dt));
%
act = mfrc_in(idp1:idp2);
%
% Convert Times from ms to s
%
t = t/1000;
tstep = tstep/1000;
tlens = tlens/1000;
tlenh = tlenh/1000;
%
fclose(fid);
%
return