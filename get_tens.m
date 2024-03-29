function [srate,pass,act,idv,nv] = get_tens(n,rnums,t)
%GET_TENS  Prompts user for passive and active muscle forces, and to
%          verify sampling rate and valid runs.
%
%          [SRATE,PASS,ACT,IDV] = GET_TENS(N,RNUMS,T) Given the number
%          of runs, N, the run names in a character cell array, RNUMS,
%          and the time vector in a cell array, T, prompts the user for
%          passive and active muscle forces, and to verify sampling
%          rate and valid runs.  The function returns the sampling
%          rates in Hz, SRATE, the passive muscle forces, PASS, the
%          active muscle forces, ACT, and a logical index to the valid
%          runs, IDV.
%
%          [SRATE,PASS,ACT,IDV,NV] = GET_TENS(N,RNUMS,T) returns the
%          number of valid runs, NV.
%
%          NOTES:  1.  The input cell arrays should be of length N.
%                  The length of the output arrays will equal the
%                  number of valid runs.
%
%         27-Feb-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Inputs
%
if (nargin<3)
  error(' *** ERROR in get_tens:  Three inputs are required!');
end
%
% Check Lengths of Input Cell Arrays
%
lt = length(t);
if (length(rnums)~=lt)||(lt~=n)
  error([' *** ERROR in get_tens:  Length of input cell arrays', ...
         ' must equal the number of runs!']);
end
%
% Setup Figure and Panel for Table
%
hf = figure('Position',[200 250 900 300],'Name', ...
            'Passive/Active Tensions','NumberTitle','off');
%
ptit = ['Please enter passive and active tensions.  Only checked ', ...
        'runs will be used in the analysis.'];
hp = uipanel('Parent',hf,'Title',ptit,'FontSize',24, ...
     'TitlePosition','centertop','FontSize',10,'FontWeight','bold', ...
     'BackgroundColor','white','Position',[0.01 0.01 0.98 0.98]);
%
uicontrol('Parent',hp,'Style','pushbutton','String','Proceed', ...
          'Position',[550 4 72 36],'CallBack','uiresume(gcbf)');
%
% Table Data
%
dat = cell(n,5);
dat(:,1) = rnums(:);
srate = cellfun(@diff,t(:),'UniformOutput',false);
srate = cellfun(@(x) rdivide(1,x),srate,'UniformOutput',false);
dat(:,2) = cellfun(@mean,srate,'UniformOutput',false);
vdat = true(n,1);
dat(:,5) = mat2cell(vdat,ones(n,1));
%
% Table Column and Row Names, Table Format, and Editable Columns
%
cnams =   {'Run Name','Sampling Rate','Passive Tensions', ...
           'Active Tensions','Valid Runs?'};
rnams = [repmat('Run ',n,1) int2str((1:n)')];
tfrmt = {'char','short','long','long','logical'};     % Format
coledit =  [false true(1,4)];          % Editable columns
%
% Use Table to Get Passive and Active Muscle Forces, and Verify
% Sampling Rates and Valid Runs
%
htab = uitable('Parent',hp,'Units','normalized','Position', ...
               [0.075 0.15 0.85 0.80],'Data', dat,'ColumnWidth', ...
               {'auto','auto',216,216,'auto'}, 'ColumnName', cnams, ...
               'ColumnFormat',tfrmt,'ColumnEditable',coledit, ...
               'RowName',rnams);
%
uiwait(hf);
%
% Get User Data Input
%
tabdat = get(htab,'Data');
%
close(hf);              % Close input table
%
% Get Variables from the Data Cell Array
%
idv = [tabdat{:,5}]';
srate = [tabdat{idv,2}]';
pass = [tabdat{idv,3}]';
act = [tabdat{idv,4}]';
nv = sum(idv);
%
return