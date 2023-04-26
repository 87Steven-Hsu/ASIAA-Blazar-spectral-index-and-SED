%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to do structure function analysis for single ALMA Calibrator Catalog (ACC) source
% Only analysis light curve in band 3, 6, and 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sourcename, tau_3, tau_6, tau_7, ipl_3, ipl_6, ipl_7, D_1000, D1000_err_up,...
    D1000_err_down, fowb, tauchar_error_up,tauchar_error_down] = ACCgetlightcurve_1000(sourcename, dtlag, makeplots)
%clc;clear;close all;

% sf = sf ampitude
% seesfamp = error of sf ampitude
% Serrmean = error of lightcurve
% dtlag = time lag
dtlag = 4;      % set in days (e.g. 4) to fix value for minimum timescale
makeplots = 1;  % set to 1 to make plots, 0 for no plotting
 
%%% test source
%sourcename = 'J0006-0623.csv';
%sourcename = 'J0348-2749.csv';

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get lightcurve data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name = sourcename(1:10)
    filename = strcat('ACC lightcurve flagged gt 100 points/',sourcename); 
    fid = fopen(filename);  % openfile

    for i = 1:1
        fgetl(fid);  % skip header of first line
    end

    Source = textscan(fid,'%f %f %f %f %f %f %f %f %f','delimiter',','); 
    fclose(fid);           
    
    figure(1);
    set(gcf,'pos',[100,100,800,950]); 

    for i = 1:3 % 1 => band 3, 2 => band 6, 3 => band 7
        disp(['i = ', num2str(i) ])

        mjd = Source{3*i-2}';    % mjd = modified Julian date   % F{1} = first row of F
        S = Source{3*i-1}';      % S = flux
        Serr = Source{3*i}';   % Serr = flux error
        
        notnanind = find(~isnan(S));
        mjd = mjd(notnanind);
        S = S(notnanind);
        Serr = Serr(notnanind);
    
        soldays = mjd - mjd(1); 
        days = soldays.*0.99726957;   % convert to units of sidereal day 

        data_length = 100;
        if length(S) > data_length
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % extract SF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            nlag = round( (max(mjd) - min(mjd))./dtlag ); % number of timelag bins
    
            %%% find/selecting the low and high states
            % calculate mean flux
            threshhilo =  mean(S); % S = flux
    
            [sf, stdsf, ~, tau, nsf] = strfnerrOVRO(days, S, Serr, dtlag, nlag);  % get SF for full lightcurve   
            % nlag = number of timelag bins % ~ = errwalt (not use in the code)
    
            nprth = 30; % ####(ADJUSTABLE PARAMETER)#### 
    
            % only select bins in which the number of pairs are above a certain threshold
            % threshold as a percentage of total number of points in lightcurve
            ipl = find(nsf > nprth & tau ~= 0);  % nsf & tau size = 926
            % insert tau into certain band 
            if i == 1
                days_3 = days;
                tau_3 = tau;
                ipl_3 = ipl;
            elseif i == 2
                days_6 = days;
                tau_6 = tau;
                ipl_6 = ipl;
            else
                days_7 = days;
                tau_7 = tau;
                ipl_7 = ipl;
            end
        
            % calculate SF error bars and errors of bin size
            SFerrorbar = 3;
            if (SFerrorbar == 1)   % Option 1 (std((Sj - Sk)^2)/sqrt(N-1)
                errsf = stdsf./(sqrt(nsf-1)); 
            elseif (SFerrorbar == 2)    % Option 2 (You et al., 2007)
                errsf = mean(sf(ipl))*(sqrt(tau/max(days)));
            elseif (SFerrorbar == 3)    % Option 3 (combination of option 1 an option 2!
                errsf = stdsf./(sqrt(nsf-1)).*sqrt(tau/max(days)); 
            end
    
            errsfbin = ones(numel(sf))*(dtlag./2);  % numel = returns the number of elements
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % return variables
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %errsfamp = errsf(ipl);  % SF amplitude error
            %sftau = tau(ipl);   % SF tau
            %dtlagused = dtlag;  % time lag 
            Serrmean = median(Serr)./mean(S);  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SF fitting
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            horline = 1:1:5000;
    
            Dnoise = 2*Serrmean^2; 
            g1000 = fittype('a*(1+b/1000)/(1+b/x) + Dnoise','coefficients',{'a','b'},'independent',{'x'},'problem',{'Dnoise'}); 
            % a = D(1000d) b = tau_char x = tau
        
            errsfuse = errsf(ipl);  
            meansf = mean(sf(ipl)); 
            weight = (meansf./errsfuse);
    
            x = tau(ipl);
            y = sf(ipl);
            xx = 1:max(x);
    
            if sum(isnan(x)) == 0 && sum(isnan(y)) == 0 && sum(isnan( weight)) == 0 && length(x) > 2 && length(y) > 2
                %%% SF fitting
                startPoints = [min(x), min(y)];
                %fo1000 = fit(x',aa',g1000,'problem',Dnoise, 'Start', startPoints); % fitting without weight
                %fo1000 = fit(x', y',g1000,'problem',Dnoise, 'Start', startPoints); % fitting without weight
                %foa = fo1000.a; fob = fo1000.b;
                fow1000 = fit(x',y',g1000, 'weight', weight, 'problem', Dnoise, 'Start', startPoints); % fitting with weight
                fowa = fow1000.a; %fowa = SF amplitude (D(1000d))
                % fowb = characteristic timescale
                fowb(1,i) = fow1000.b; 
                fityw = fowa*(1+fowb(1,i)./1000)./(1+fowb(1,i)./xx) + Dnoise;

                % charactristic timescale uncertainty calculattion
                ci = confint(fow1000,0.95);
                tauchar_error_up(1,i) = ci(4);
                tauchar_error_down(1,i) = ci(3);
  
                % get D(1000) value
                if length(fityw) > 1000
                    D_1000(1,i) = fityw(1000);
                else
                    D_1000(1,i) = nan;
                end
            
                % get D(1000) error using "predint"
                p22 = predint(fow1000,xx,0.95,'functional','on');
                if length(p22) > 1000
                    D1000_err_up(1,i) = p22(1000,2);
                    D1000_err_down(1,i) = p22(1000,1);
                else
                 D1000_err_up(1,i) = nan;
                    D1000_err_down(1,i) = nan;
                end

                disp('SF fitting proceed sucessfully')
                disp('data length = '+string(length(S)))
            else
                fityw = nan(length(xx));
                p22 = nan(length(xx), 2);
                x = nan(20, 1);
                y = nan(20, 1);

                D_1000(1,i) = nan;
                D1000_err_up(1,i) = nan;
                D1000_err_down(1,i) = nan;

                fowb(1,i) = nan;
                tauchar_error_up(1,i) = nan;
                tauchar_error_down(1,i) = nan;

                disp('nan exist in x, y, or weight, no SF fitting proceed')
            end
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot lightcurves and SF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(6,1, 2*i-1); 
            
                plot(days,S,'ko','MarkerSize',3,'MarkerFaceColor','k') 
                hold on
                errorbar(days,S,Serr,Serr,'k','MarkerSize',3,'MarkerFaceColor','k')
                plot(horline, ones(size(horline))*threshhilo,'r--','LineWidth',1.5)
                hold off
                %xlim([0 max(days_3)]) 
                %ylim([0.5*min(S) 1.2*max(S)])
                ylabel('S (Jy)','FontSize', 15) 
                set(gca,'FontSize',20,'LineWidth',2); 

                if ~isempty(ipl) || ~isnan(Dnoise) || ~isempty( sf(ipl) ) 
                    if ~isnan(max(sf(ipl)))
                        xlim([0 max(tau(ipl))])
                        ylim([-0.5*Dnoise 1.2*max(S)])
                    else
                        xlim([0 1500])
                        ylim([0 1])
                    end
                else
                    xlim([0 1500])
                    ylim([0 1])
                end

            subplot(6,1, 2*i);  
        
                plot(x, y,'k.'); % origin SF data
                hold on
                %%% SF fitting line
                
                noise_arr = [];
        
                %%% plot fitting line without weight
                %fity = foa*(1+fob./1000)./(1+fob./xx) + Dnoise;
                %plot(xx,fity,'r','LineWidth',2, 'color', 'green'); 
                %%% plot fitting line with weight
                plot(xx,fityw,'r','LineWidth',2, 'color', 'blue'); 

                %%% "predint" method uncertainty
                plot(xx,p22(:,2),'m--','LineWidth',2); % upper limit
                plot(xx,p22(:,1),'m--','LineWidth',2); % lower limit
                noise_arr(1:length(x)) = Dnoise;
                plot(x,noise_arr,'b--','LineWidth',2);    % Dnoise
                e1 = errorbar(tau(ipl),sf(ipl),errsf(ipl),'ko','MarkerSize',3,'MarkerFaceColor','k');  
                e1.Annotation.LegendInformation.IconDisplayStyle = 'off'; 
                e2 = errorbar(tau(ipl),sf(ipl), errsfbin(ipl),'horiz','ko','MarkerSize',3,'MarkerFaceColor','k'); 
                e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                hold off
                %legend('D(\tau)','Dnoise','D(\tau) fitting','D(\tau) fitting with weight') 
                
                if ~isempty(ipl) || ~isnan(Dnoise) || ~isempty( sf(ipl) ) 
                    if ~isnan(max(sf(ipl)))
                        xlim([0 max(tau(ipl))])
                        ylim([-0.5*Dnoise 1.2*max(sf(ipl))])
                    else
                        xlim([0 1500])
                        ylim([0 1])
                    end
                else
                    xlim([0 1500])
                    ylim([0 1])
                end
                
                ylabel('D(\tau)','FontSize',15)
                set(gca,'FontSize',20,'LineWidth',2);
                
                if i == 3
                    xlabel('\tau (Sidereal Days)','FontSize',15)
                    sgt = sgtitle([name,' light curve, SF & SF fitting']);
                    sgt.FontSize = 20;
                    sgt.LineWidth = 2; 
                    set(gca,'FontSize',20,'LineWidth',2);
                end

                %%% save file 
                %saveas(gcf,['/Users/87steven/Documents/ASIAA/ACC variability/light curve and SF gt 100 matlab_0522/',name,' lightcurve+SF.png']);

        else
            
            if i == 1
                disp('Band 3 data length is lower than '+string(data_length))
                days_3 = nan;
                tau_3 = nan;
                ipl_3 = nan;
            elseif i == 2
                disp('Band 6 data length is lower than '+string(data_length))
                days_6 = nan;
                tau_6 = nan;
                ipl_6 = nan;
            else
                disp('Band 7 data length is lower than '+string(data_length))
                days_7 = nan;
                tau_7 = nan;
                ipl_7 = nan;
            end

            D_1000(1,i) = nan;
            D1000_err_up(1,i) = nan;
            D1000_err_down(1,i) = nan;

            fowb(1,i) = nan;
            tauchar_error_up(1,i) = nan;
            tauchar_error_down(1,i) = nan;

        end
        
    end
    
end