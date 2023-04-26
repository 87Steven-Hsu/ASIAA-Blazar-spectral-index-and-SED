%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A code to cross match ALAM Calibrator Catalog and ROMA-BZCAT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;

%%% Read file of ACC (# = 47115)
DataStartLine = 1;
NumVariables = 10;
VariableNames  = {'name','class','z','RA','DEC','flux','errflux','band','freq','date'};
VariableWidths = [11,4,8,10,10,10,10,3,10,22] ;                                                 
DataType       = {'string','double','double','double','double','double','double','double','double','char'};
opts = fixedWidthImportOptions('NumVariables',NumVariables,...
                               'DataLines',DataStartLine,...
                               'VariableNames',VariableNames,...
                               'VariableWidths',VariableWidths,...
                               'VariableTypes',DataType);
ACC = readtable(['/Users/87steven/Documents/ASIAA/Blazar SED code and data/' ...
    'ALMA calibrator catalog (ACC)/ACC.txt'],opts);

nameACC = ACC.name;
classACC = ACC.class;
zACC = ACC.z;
RAACC = ACC.RA;     % [deg]
DECACC = ACC.DEC;   % [deg]
flux = ACC.flux;
errflux = ACC.errflux;
band = ACC.band;
freq = ACC.freq;
date = ACC.date;

Date = string(date);

%%% Read file of ACC list (# = 3361)
DataStartLine = 1;
NumVariables = 6;
VariableNames  = {'name','class','z','RA','DEC','Num_obs'};
VariableWidths = [10,2,7,9,9,5] ;                                                 
DataType       = {'string','double','double','double','double','double'};
opts = fixedWidthImportOptions('NumVariables',NumVariables,...
                               'DataLines',DataStartLine,...
                               'VariableNames',VariableNames,...
                               'VariableWidths',VariableWidths,...
                               'VariableTypes',DataType);
ACClist = readtable(['/Users/87steven/Documents/ASIAA/Blazar SED code and data/' ...
    'ALMA calibrator catalog (ACC)/ACC_list.txt'],opts);

nameACClist = ACClist.name;
classACClist = ACClist.class;
zACClist = ACClist.z;
RAACClist = ACClist.RA;     % [deg]
DECACClist = ACClist.DEC;   % [deg]
Numobslist = ACClist.Num_obs;

%%% Read file of BZCAT (# = 3561)
DataStartLine = 1;
NumVariables = 20;
VariableNames  = {'Seq','type','name','RA_h','RA_m','RA_s','DEC_d','DEC_m','DEC_s','l','b','z',...
    'z_uncer','R_mag','class','flux_ratio','flux_143','flux_xray','flux_gamma','Spec_index'};
VariableWidths = [8,7,12,9,8,8,8,8,8,8,8,8,7,9,29,9,8,8,12,8] ;                                                 
DataType       = {'double','string','string','double','double','double','double','double','double',...
    'double','double','double','string','double','string','double','double','double','double','double'};
opts = fixedWidthImportOptions('NumVariables',NumVariables,...
                               'DataLines',DataStartLine,...
                               'VariableNames',VariableNames,...
                               'VariableWidths',VariableWidths,...
                               'VariableTypes',DataType);
BZCAT = readtable('/Users/87steven/Documents/ASIAA/Blazar SED code and data/ROMA-BZCAT/BZCAT.txt',opts);

Seq = BZCAT.Seq;
type = BZCAT.type;
nameCAT = BZCAT.name;
RA_h = BZCAT.RA_h;
RA_m = BZCAT.RA_m;
RA_s = BZCAT.RA_s;
DEC_d = BZCAT.DEC_d;
DEC_m = BZCAT.DEC_m;
DEC_s = BZCAT.DEC_s;
l = BZCAT.l;
b = BZCAT.b;
zCAT = BZCAT.z;
z_uncer = BZCAT.z_uncer;
R_mag = BZCAT.R_mag;
classCAT = BZCAT.class;
flux_ratio = BZCAT.flux_ratio;
flux_143 = BZCAT.flux_143;
flux_xray = BZCAT.flux_xray;
flux_gamma = BZCAT.flux_gamma;
Spec_index = BZCAT.Spec_index;

RA = RA_h + RA_m/60 + RA_s/3600;    % [hr]
RACAT = RA*15;                      % [deg]

for i = 1:max(size(DEC_d))
    
    if DEC_d(i) < 0
    
        DECCAT(i) = DEC_d(i) - DEC_m(i)/60 - DEC_s(i)/3600;
    
    else
   
        DECCAT(i) = DEC_d(i) + DEC_m(i)/60 + DEC_s(i)/3600;  % [deg]
    
    end
    
end

DECCAT = DECCAT';

%% Plot ACC and Roma-BZCAT
figure(1);
    set(gcf,'pos',[100,100,1200,800]);
    
    plot(RACAT,DECCAT,'r.','MarkerSize',6,'MarkerFaceColor','k') 
    hold on
    plot(RAACClist,DECACClist,'g.','MarkerSize',6,'MarkerFaceColor','k') 
    hold off
    
    grid on
    xlim([0 360])   
    ylim([-90 100]) 
    xlabel('RA (J2000) [deg]','FontSize', 15)
    ylabel('DEC (J2000) [deg]','FontSize', 15)   
    
    leg1 = 'Roma-BZCAT, N = '+ string(max(size(RACAT)));
    leg2 = 'ACC, N = '+ string(max(size(nameACClist)));
    legend(leg1, leg2, 'location','northeast')

    set(gca,'FontSize',15,'LineWidth',2); 

    sgt = sgtitle('Roma-BZCAT and ALMA Calibrator Catalog (ACC) Distribution'); 
    sgt.FontSize = 20;
    sgt.LineWidth = 2;
    
%% Check distribution of ACC and ACC list
figure(2);
    set(gcf,'pos',[100,100,1200,800]);
    
    plot(RAACC,DECACC,'ro','MarkerSize',6) 
    hold on
    plot(RAACClist,DECACClist,'g.','MarkerSize',6,'MarkerFaceColor','k') 
    hold off
    
    grid on
    xlim([0 360])   
    ylim([-90 100]) 
    xlabel('RA (J2000) [deg]','FontSize', 15)
    ylabel('DEC (J2000) [deg]','FontSize', 15)   
    
    leg1 = 'ACC, N = '+ string(max(size(RAACC)));
    leg2 = 'ACC List, N = '+ string(max(size(RAACClist)));
    legend(leg1, leg2, 'location','northeast')

    set(gca,'FontSize',15,'LineWidth',2); 

    sgt = sgtitle('ACC and ACC List Distribution'); 
    sgt.FontSize = 20;
    sgt.LineWidth = 2;
%% distance calculation between Roma-BZCAT and ACC List coordinates
d = zeros(max(size(nameCAT)),max(size(nameACClist)));

for i = 1:max(size(nameCAT))
    
    for j = 1:max(size(nameACClist))
        
        d(i,j) = sqrt( (RACAT(i)-RAACClist(j))^2 + (DECCAT(i)-DECACClist(j))^2 );
        
    end
    
end

%% Find the cloest distance and distance lower than 0.5 deg
Xresult = zeros(max(size(nameCAT)),2);
matrad = 10/3600;   % [deg]

for i = 1:max(size(nameCAT))
    
    dismin = find(d(i,:) == min(d(i,:))); % find minimum distance of blazar sources and WHAM point
    
    if d(i,dismin) <= matrad            % distance is less than "matrad" deg 
        
%         for j = 1:max(size(dismin))
%             
%             k = 2*j-1;
%             
%             Xresult(i,k) = RAACC(dismin(j));   %column 1 will assign ACC RA
%             Xresult(i,k+1) = DECACC(dismin(j));  %column 2 will assign ACC DEC
%             
%         end

        Xresult(i,1) = RAACClist(dismin);            
        Xresult(i,2) = DECACClist(dismin);            
        
    else                                % if distance is larger than 10"
        Xresult(i,1) = NaN;             %column 1 will assign NaN
        Xresult(i,2) = NaN;               %column 2 will assign NaN
    end
    
end
%% Xmatch result check
matchsucc = find(~isnan(Xresult(:,1))); % find match success sources' index
 
figure(3);
    set(gcf,'pos',[100,100,1200,800]);
    
    plot(RACAT(matchsucc),DECCAT(matchsucc),'ro','MarkerSize',6)    % plot Roma-BZCAT
    hold on 
    plot(Xresult(:,1),Xresult(:,2),'g.','MarkerSize',6)     % plot ACC list
    hold off
    
    xlabel('RA (J2000) [deg]','FontSize', 15)
    ylabel('DEC (J2000) [deg]','FontSize', 15)   
    xlim([0 360])   
    ylim([-90 90])   
    grid on
    legend('Roma-BZCAT','ACC','location','northeast')
    
    set(gca,'FontSize',15,'LineWidth',2); 
    
    sgt = sgtitle('Cross Match of Roma-BZCAT with ACC List, r = 10", N = '+string(max(size(matchsucc)))); 
    sgt.FontSize = 20;
    sgt.LineWidth = 2;
%% find same ACC list coordinate with ACC. Than, find out corresponding ACC source name
fs = zeros(max(size(matchsucc)),1000);

for i = 1:max(size(matchsucc))
        
    % find same coordinate of cross match result with ACC list, than, extract its index
    aa = find(RAACClist == Xresult(matchsucc(i),1) & DECACClist == Xresult(matchsucc(i),2));    
    
    findsame = strcmpi(nameACC,nameACClist(aa));    % use the index to find source name in ACC
    
    nameACCind = find(findsame ~= 0);
    
    for j = 1:max(size(nameACCind))

        fs(i,j) = nameACCind(j);    % save ACC source index into an array
        
    end

end
%% save cross match data to a file
filename = strcat('Xmatch_result_20220216.txt');
fileID = fopen(filename,'w');

for i = 1:1367 % max(size(matchsucc))

% 	for k = 1:max(size(find((~isnan(Xresult(matchsucc(i),:)) & Xresult(matchsucc(i),:) ~=0))))/2
%         
%             j = find(RAACC == Xresult(matchsucc(i),1) & DECACC == Xresult(matchsucc(i),2));
%         
%             fprintf(fileID,'%s\t %10.6f\t %10.6f\t %10.8f\t %10.8f\t %s\t\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %s\n',...
%                 nameCAT(matchsucc(i)), RACAT(matchsucc(i)), DECCAT(matchsucc(i)), zCAT(matchsucc(i)),...
%                 R_mag(matchsucc(i)), classCAT(matchsucc(i)), flux_ratio(matchsucc(i)),flux_143(matchsucc(i)),...
%                 flux_xray(matchsucc(i)), flux_gamma(matchsucc(i)), Spec_index(matchsucc(i)), flux(j(k)),...
%                 errflux(j(k)), band(j(k)), freq(j(k)), Date(j(k)));
%             
%     end
    
    for k = 1:max(size(find(fs(i,:) ~= 0)))
       
        fprintf(fileID,['%s\t %10.6f\t %10.6f\t %10.6f\t %10.6f\t %10.6f\t %10.6f\t %10.6f\t %10.6f\t %10.8f\t %10.8f\t %s\t\t ' ...
            '%10.8f\t %10.8f\t %10.8f\t %12.10f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %10.8f\t %s\n'],...
                nameCAT(matchsucc(i)), RA_h(matchsucc(i)), RA_m(matchsucc(i)), RA_s(matchsucc(i)), DEC_d(matchsucc(i)),...
                DEC_m(matchsucc(i)), DEC_s(matchsucc(i)), RACAT(matchsucc(i)), DECCAT(matchsucc(i)), zCAT(matchsucc(i)),...
                R_mag(matchsucc(i)), classCAT(matchsucc(i)), flux_ratio(matchsucc(i)),flux_143(matchsucc(i)),...
                flux_xray(matchsucc(i)), flux_gamma(matchsucc(i)), Spec_index(matchsucc(i)), flux(fs(i,k)),...
                errflux(fs(i,k)), band(fs(i,k)), freq(fs(i,k)), Date(fs(i,k)));
        
    end
    
end

fclose(fileID);

%% Cross match Roma-BZCAT with ACC List by Source Name (FULL LIST)
samename = zeros(max(size(nameCAT)),max(size(nameACClist)));

for i = 1:max(size(nameCAT))
        
    for j = 1:max(size(nameACClist))
        
        samename(i,j) = strcmpi(nameCAT(i),nameACClist(j));

    end
    
end

[row,col] = find(samename == 1);

%% Plot corss match by source name
figure(4);
    set(gcf,'pos',[100,100,1200,800]);
    
    plot(RACAT(row),DECCAT(row),'bo','MarkerSize',10)                   % plot Roma-BZCAT searce by name
    hold on
    plot(RAACClist(col),DECACClist(col),'b.','MarkerSize',10)           % plot ACC list searce by name
    %plot(RACAT(matchsucc),DECCAT(matchsucc),'go','MarkerSize',10)       % plot Roma-BZCAT search by r = 10"
    %plot(Xresult(:,1),Xresult(:,2),'g.','MarkerSize',10)                % plot ACC list search by r = 10"
    hold off
    
    grid on
    xlim([0 360])   
    ylim([-90 100]) 
    xlabel('RA (J2000) [deg]','FontSize', 15)
    ylabel('DEC (J2000) [deg]','FontSize', 15)   
    
    legend('Roma-BZCAT matched by name', 'ACC matched by name',...
        'Roma-BZCAT matched by radius','ACC matched by radius', 'location','northeast')

    set(gca,'FontSize',15,'LineWidth',2); 

    sgt = sgtitle('BZCAT and ACC Xmatch by Source Name, N = '+string(max(size(row)))); 
    sgt.FontSize = 20;
    sgt.LineWidth = 2;

%% Find sources matched by source name 
samename = zeros(max(size(nameCAT(matchsucc))), max(size(nameCAT(row))));

for i = 1:max(size(nameCAT(matchsucc)))
        
    for j = 1:max(size(nameCAT(row)))
        
        samename(i,j) = strcmpi(nameCAT(matchsucc(i)),nameCAT(row(j))); % comparing A: (Xmatch success source name) with B: (name matched with ACC and BZCAT)
        
    end
    
end

[row1,col1] = find(samename == 1); % extract index

aa = zeros(max(size(samename)), 1);
for i = 1:max(size(samename))
    
    if find(samename(:,i) == 1)    % find 1 (source matched with A and B)

        aa(i) = 1;
        
    end
  
end

cc = find(aa == 0);     % find no matched source 
%% Distance calculation of sources matched with name but not with distance
for i = 1:max(size(row(cc)))
        
    dis(i) = sqrt( (RACAT(row(cc(i)))-RAACClist(col(cc(i))))^2 + (DECCAT(row(cc(i)))-DECACClist(col(cc(i))))^2 );      
        
end

%%% find source name
    %nameCAT(row(cc))
    %nameACClist(col(cc))
 
%%
figure(5);
    set(gcf,'pos',[100,100,1200,800]);
    
    plot(RACAT(row(cc)),DECCAT(row(cc)),'bo','MarkerSize',10) 
    hold on
    plot(RAACClist(col(cc)),DECACClist(col(cc)),'b.','MarkerSize',10) 
    
    plot(RACAT(matchsucc),DECCAT(matchsucc),'go','MarkerSize',10)       % plot Roma-BZCAT search by r = 10"
    plot(Xresult(:,1),Xresult(:,2),'g.','MarkerSize',10)                % plot ACC list search by r = 10"
    
    hold off
    
    grid on
    xlim([0 360])   
    ylim([-90 100]) 
    xlabel('RA (J2000) [deg]','FontSize', 15)
    ylabel('DEC (J2000) [deg]','FontSize', 15)   
    
    legend('BZCAT matched by name, N = 1390', 'ACC matched by name, N = 1390',...
        'BZCAT matched by radius, N = 1367','ACC matched by radius, N = 1367', 'location','northeast')

    set(gca,'FontSize',15,'LineWidth',2); 

    sgt = sgtitle('Sources Xmatch Succfully by Source Name BUT Not With r = 10", N = '+string(max(size(row(cc))))); 
    sgt.FontSize = 20;
    sgt.LineWidth = 2;
%% Find source matched with r = 10", but not matched with source name
s = zeros(max(size(nameCAT(matchsucc))), max(size(nameCAT(matchsucc(row1)))));

for i = 1:max(size(nameCAT(matchsucc)))
        
    for j = 1:max(size(nameCAT(matchsucc(row1))))
        
        s(i,j) = strcmpi(nameCAT(matchsucc(i)),nameCAT(matchsucc(row1(j)))); % find BZCAT source matched by r = 10" with source matched by source name
        
    end
    
end

aa = zeros(max(size(s)), 1);
for i = 1:max(size(s))
    
    if find(s(i,:) == 1)    % find 1 (source matched by r = 10" and matched by source name)

        aa(i) = 1;
        
    end
  
end

cc = find(aa == 0);     % find no matched source 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Search RESULT %%%%
% cc = 1174
% => matchsucc(1174) == 3157 (BZCAT)
% => nameCAT(3157) == "J2117+0502" (BZCAT)
% Using code line 223-240 to find out corresponding ACC name
% => Corresponding ACC index: 42800-42802
% => nameACC(42800:42802) == "J2117+0503"
% distance bwtween "J2117+0502" and "J2117+0503" is given by,
% => d = sqrt( (RACAT(3157)-RAACC(42800))^2 + (DECCAT(3157)-DECACC(42800))^2 )
% => d == 0.0018 [deg]
% 10" = 10/3600 deg = 0.0028 deg
% => These two sources ("J2117+0502" and "J2117+0503") can be regarded as
% same source 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dis = zeros(max(size(row(cc))),1);
for i = 1:max(size(row(cc)))
        
    dis(i) = sqrt( (RACAT(row(cc(i)))-RAACClist(col(cc(i))))^2 + (DECCAT(row(cc(i)))-DECACClist(col(cc(i))))^2 );      
        
end

























