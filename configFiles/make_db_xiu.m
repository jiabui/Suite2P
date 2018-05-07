% UPDATE Christmas 2016: number of clusters determined automatically, but
% do specify the "diameter" of an average cell for best results. You can do this with either
% db(iexp).diameter, or ops0.diameter.
%%
function db = make_db_xiu(mouse_str,date_str,range_session,isconcat, diameter)
if ~exist('diameter','var')
    diameter = 20;
end

if ~exist('isconcat','var')
    isconcat = 0;
end

if ~isconcat
    for i = 1:length(range_session)
        db(i).mouse_name    = mouse_str;
        db(i).date          = date_str;
        db(i).expts         = range_session(i);
        db(i).diameter      = diameter;
    end
    
else
    i = 1;
    c1 = cell(1,length(range_session));
    c2 = cell(1,length(range_session));
    c3 = cell(1,length(range_session));
    for j = 1:length(range_session)
        c1{j} = mouse_str;
        c2{j} = date_str;
        c3{j} = range_session(j);
    end
    db(i).mouse_name    = c1;
    db(i).date          = c2;
    db(i).expts         = c3;
    db(i).diameter      = diameter;
end
%     i = 0;
    
%     i = i+1;
%     db(i).mouse_name    = {'890C','890C','890C'};
%     db(i).date          = {'122117','122117','122117'};%'121917';
%     db(i).expts         = {3, 4, 5};
%     db(i).diameter      = 20;
%     
    % db(i).mouse_name    = {'MK020', 'M150416_MK020'};
    % db(i).date          = {'2015-07-30', '2015-07-30'};
    % db(i).expts         = {[2010 2107], [1 2 3]};
    % db(i).diameter      = 12;
end

%%

% i = 0;
%
% i = i+1;
% db(i).mouse_name    = 'm011';
% db(i).date          = '120717';
% db(i).expts         = [4];
% db(i).diameter      = 20;
%



%%
% i = 0;
%
% i = i+1;
% db(i).mouse_name    = 'M326013-B';
% db(i).date          = '2017-06-29';
% db(i).expts         = [1];
% db(i).diameter      = 20;

%%

% %
% i = i+1;
% db(i).mouse_name    = 'M150329_MP009';
% db(i).date          = '2015-04-10';
% db(i).expts         = [5 6 7 8 9 10 11];
% db(i).diameter      = 12;
%
% i = i+1;
% db(i).mouse_name    = 'M150824_MP019';
% db(i).date          = '2015-12-19';
% db(i).expts         = [4];
% db(i).diameter      = 6;
%
% % example of datasets, which consist of several sessions - use cell arrays
% % will be treated as subsets of experiment with the same FOV, with
% % different names/dates (for one reason or another), analyzed together
% i = i+1;
% db(i).mouse_name    = {'MK020', 'M150416_MK020'};
% db(i).date          = {'2015-07-30', '2015-07-30'};
% db(i).expts         = {[2010 2107], [1 2 3]};
% db(i).diameter      = 12;
%
% % example for datasets without folder structure
% db(i).mouse_name    = 'notImportant';
% db(i).date          = '2016';
% db(i).expts         = []; % leave empty, or specify subolders as numbers
% db(i).diameter      = 12;
% db(i).RootDir       = 'F:\DATA\neurofinder\neurofinder.01.00\images'; % specify full path to tiffs here
%
% % example extra entries
% % db(i).AlignToRedChannel= 1;
% % db(i).BiDiPhase        = 0; % adjust the relative phase of consecutive lines
% % db(i).nSVD             = 1000; % will overwrite the default, only for this dataset
% % db(i).comments      = 'this was an adaptation experiment';
% % db(i).expred        = [4]; % say one block which had a red channel
% % db(i).nchannels_red = 2; % how many channels did the red block have in total (assumes red is last)
