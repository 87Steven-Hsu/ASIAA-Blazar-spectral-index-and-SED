%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a code to run "ACCgetlightcurve_1000.m" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all;

fileFolder = fullfile('/Users/87steven/Documents/ASIAA/ACC variability new/ACC lightcurve flagged gt 100 points');    
dirOutput = dir(fullfile(fileFolder,'*'));
fileName = {dirOutput.name}; 
%%
for k = 1:length(fileName)-3
    
    disp(['k = ', num2str(k)])
    sourcename = char(fileName(k+3));

    dtlag = 4;
    makeplots = 1;
         
    [sourcename, tau_3, tau_6, tau_7, ipl_3, ipl_6, ipl_7, D_1000, D1000_err_up,...
    D1000_err_down, fowb, tauchar_error_up,tauchar_error_down] = ACCgetlightcurve_1000(sourcename,dtlag,makeplots);

    D_1000_3 = D_1000(1);
    D_1000_6 = D_1000(2);
    D_1000_7 = D_1000(3);

    D1000_uperr_3 = D1000_err_up(1);
    D1000_uperr_6 = D1000_err_up(2);
    D1000_uperr_7 = D1000_err_up(3);
    
    D1000_lowerr_3 = D1000_err_down(1);
    D1000_lowerr_6 = D1000_err_down(2);
    D1000_lowerr_7 = D1000_err_down(3);

    fowb_3 = fowb(1);
    fowb_6 = fowb(2);
    fowb_7 = fowb(3);

    if ~isempty( ipl_3 ) && ~isnan(fowb_3) 
        if fowb_3 > max(tau_3(ipl_3)) || isnan(fowb_3) 
            fowb_3 = max(tau_3(ipl_3));
        end
    end
    if ~isempty( ipl_6 ) && ~isnan(fowb_6) 
        if fowb_6 > max(tau_6(ipl_6)) || isnan(fowb_6) 
            fowb_6 = max(tau_6(ipl_6));
        end
    end
    if ~isempty( ipl_7 ) && ~isnan(fowb_7) 
        if fowb_7 > max(tau_7(ipl_7)) || isnan(fowb_7) 
            fowb_7 = max(tau_7(ipl_7));
        end
    end

    tauchar_uperror_3 = tauchar_error_up(1);
    tauchar_uperror_6 = tauchar_error_up(2);
    tauchar_uperror_7 = tauchar_error_up(3);

    tauchar_lowerror_3 = tauchar_error_down(1);
    tauchar_lowerror_6 = tauchar_error_down(2);
    tauchar_lowerror_7 = tauchar_error_down(3);

    sourcename_m(k,:) = table(string(sourcename(1:10)), D_1000_3, D1000_uperr_3, D1000_lowerr_3, ...
        D_1000_6, D1000_uperr_6, D1000_lowerr_6, D_1000_7, D1000_uperr_7, D1000_lowerr_7, ...
        fowb_3, tauchar_uperror_3, tauchar_lowerror_3, fowb_6, tauchar_uperror_6, tauchar_lowerror_6, ...
        fowb_7, tauchar_uperror_7, tauchar_lowerror_7);
   
end
filename = strcat('/Users/87steven/Documents/ASIAA/ACC variability new/result_0817.csv');
fileID = fopen(filename,'w');

writetable(sourcename_m, filename);

fclose(fileID);