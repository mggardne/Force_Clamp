function [t,mlen,mforc,mstrain,mstress,n,rnums,amp,flen,fdia_maj, ...
                                           fdia_min,fcsa] = rd_fcp(fnam)
%RD_FCP    Opens a muscle clamp .FCP output file and reads the time,
%          muscle lengths, muscle forces, muscle strain, and muscle
%          stress for the series of runs in the file.
%
%          [T,MLEN,MFORC,MSTRAIN,MSTRESS,N] = RD_FCP(FNAM) opens a
%          muscle clamp output file (*.FCP), FNAM, and reads the time,
%          T, in seconds, muscle length, MLEN, in um, muscle force,
%          MFORC, in mN, muscle strain, MSTRAIN, muscle stress,
%          MSTRESS, in mN/mm^2, and the number of runs in the file, N.
%          T, MLEN, MFORC, MSTRAIN, and MSTRESS are cell arrays of
%          length N.
%
%          [T,MLEN,MFORC,MSTRAIN,MSTRESS,N,RNUMS] = RD_FCP(FNAM)
%          returns the run numbers, RNUMS, as a cell array of texts.
%
%          [T,MLEN,MFORC,MSTRAIN,MSTRESS,N,RNUMS,AMP] = RD_FCP(FNAM)
%          returns the amplitude, AMP, of the run numbers in an array.
%
%          [T,MLEN,MFORC,MSTRAIN,MSTRESS,N,RNUMS,AMP,FLEN,FDIA_MAJ,
%          FDIA_MIN,FCSA] = RD_FCP(FNAM) returns the muscle length,
%          FLEN, in um, muscle major diameter, FDIA_MAJ, in um, muscle
%          minor diameter, FDIA_MIN, in um, and muscle cross-sectional
%          area, FCSA, in mm^2.
%
%          NOTES:  1.  The muscle clamp output file (*.FCP) must be in
%                  a standard format.
%
%                  2.  Total number of runs in an output file must not
%                  exceed ten (10).
%
%         20-Feb-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Input
%
if (nargin<1)
  error(' *** ERROR in rd_fcp:  No input file name!');
end
%
% Check Input File Name
%
idot = strfind(fnam,'.');
idot = idot(end);
if ~strcmpi(fnam(idot:end),'.fcp')
  error(' *** ERROR in rd_fcp:  Input file extension must be .FCP!');
end
%
% Open Input File
%
fid = fopen(fnam,'r');
%
% Find the Number of Runs within the File
%
ifirstl = true;         % Get first muscle length in the file
ifirstd = true;         % Get first muscle diameters in the file
ifirstc = true;         % Get first muscle cross-sectional area in the file
%
n = 0;                  % Number of runs in the file
%
rnums = cell(1,10);
%
while ~feof(fid)
     lin = fgets(fid);
     if ~feof(fid)
%
       if startsWith(lin,'Saved Run Number')
         n = n+1;       % Number of runs
         idr = strfind(lin,':')+2;
         rnums{n} = strip(lin(idr:end));    % Run number as text
       end
%
       if ifirstl
         if startsWith(lin,'Muscle Length')
           idc = strfind(lin,':')+2;
           idu = strfind(lin,'um')-1;
           flen = str2double(lin(idc(1):idu(1)));     % Muscle length in um
           ifirstl = false;
         end
       end
%
       if ifirstd
         if startsWith(lin,'Muscle Diam Major')
           idc = strfind(lin,':')+2;
           idu = strfind(lin,'um')-1;
           fdia_maj = str2double(lin(idc(1):idu(1))); % Muscle major diameter in um
           fdia_min = str2double(lin(idc(2):idu(2))); % Muscle minor diameter in um
           ifirstd = false;
         end
       end
%
       if ifirstc
         if startsWith(lin,'Muscle Cross Sectional Area')
           idc = strfind(lin,':')+2;
           idu = strfind(lin,'mm')-1;
           fcsa = str2double(lin(idc(1):idu(1)));     % Muscle cross-sectional area in mm^2
           ifirstc = false;
         end
       end
%
     end
end
%
if n==0
  error(' *** ERROR in rd_fcp:  No run numbers found in the file!');
end
%
rnums = rnums(1,1:n);
%
frewind(fid);
%
% Read Data from Each Run
%
amp = zeros(n,1);
%
hdr = ['Time (s)	Muscle Length (um)	Force (mN)	Strain	', ...
       'Stress (mN/mm^2)'];
%
frmt = '%f %f %f %f %f';               % Format of data
%
t = cell(1,n);          % Time (s)
mlen = cell(1,n);       % Muscle lengths (um)
mforc = cell(1,n);      % Muscle forces (mN)
mstrain = cell(1,n);    % Muscle strain
mstress = cell(1,n);    % Muscle stress (mN/mm^2)
%
k = 0;
l = 0;
%
while ~feof(fid)
     lin = fgets(fid);
%
     if ~feof(fid)
%
       if startsWith(lin,hdr)
         k = k+1;
         data = textscan(fid,frmt,'Delimiter','tab');
         t(1,k) = data(:,1);
         mlen(1,k) = data(:,2);
         mforc(1,k) = data(:,3);
         mstrain(1,k) = data(:,4);
         mstress(1,k) = data(:,5);
       end
%
       if startsWith(lin,'Amplitude:')
         l = l+1;
         ida = strfind(lin,':')+2;
         idp = strfind(lin,'%')-1;
         amp(l) = str2double(lin(ida:idp));
       end
%
     end
%
end
%
fclose(fid);
%
return