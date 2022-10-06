% Generation of Figs. 2-9 in the paper
% Y. Sun, K. C. Ho, G. Wang. J. Chen, Y. Yang, L. Chen, and Q. Wan, 
% "Computationally attractive and location robust estimator for IoT device positioning," 
% IEEE Internet Things J., Nov. 2021.
%
% Yimao Sun and K. C. Ho   04-08-2022
%
%       Copyright (C) 2022
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clear all;
% close all;
% clc

rng('default');
% warning off

senPos = [
        % minimum number of sensors
        10.23    38.38    16.29
        46.64   -87.12    62.94
        124.02  -7.98     81.16
        105.02  -51.72    26.47
       -81.56    104.48  -80.49
        ]';

% target direction
theta = 22.13*pi/180;
phi = 14.41*pi/180;

[N,M] = size(senPos);
mon = 1000;

aveNse = 0;
for l=1:mon
    aveNse = aveNse + randn(M,1); 
end
aveNse = aveNse/mon/sqrt(2);
PP = aveNse(2:end) - aveNse(1);

disp('Simulation is running ...');
models = ['nse';'rag'];

for im = 1:2
    model = models(im,:);
    switch model
        case 'nse'
            % ******* vs. noise power *******
            nsePwr = -40:10:50;
            souRange = 15*1e2;
        case 'rag'
            % ******* vs. range *******
            nsePwr = 0;  % 10log(rad^2)
            souRange = [1,3,5:5:80]*1e2;
    end
    
    for ir = 1:length(souRange)
        disp(['Range: ',num2str(souRange(ir)),', ',num2str(ir),'/',num2str(length(souRange)),' ...']);
        %******** Generate Data ********
        % source location
        souLoc = souRange(ir) * [cos(theta)*cos(phi); sin(theta)*cos(phi); sin(phi)] + senPos(:,1);
        % true range
        r = sqrt(sum((souLoc-senPos).^2,1))';
        % true TDOAs
        rd = r(2:end) - r(1);
        g = 1/r(1);

        [uTh,uPh,ug,ep,er] = deal(zeros(mon,1));
        for in = 1:length(nsePwr)
            disp(['10log(NoisePower): ',num2str(nsePwr(in)),', ',num2str(in),'/',num2str(length(nsePwr)),' ...']);
            Q = 10^(nsePwr(in)/10) * (ones(M-1, M-1)+eye(M-1))/2;

            % Calculate CRLB
            CRB = ConsCRLB( senPos, souLoc, Q );
            CRLB_a(in,ir) = CRB(1,1)+CRB(2,2);
            CRLB_g(in,ir) = CRB(3,3);

            rng('default');

            for i = 1:mon
                % measured TDOAs
                tmp = randn(M,1);
                rdNse = sqrt(10^(nsePwr(in)/10)) * ((tmp(2:M)-tmp(1))/sqrt(2)-PP);
                rd_m = rd + rdNse;

                nAg = 0;
                % SCO-MPR Method
                nAg = nAg + 1;
                [mprSol, ~] = TDOA_SCO_MPR( senPos, rd_m, Q );
                uTh(i,nAg) = mprSol(1);
                uPh(i,nAg) = mprSol(2);
                ug(i,nAg) = mprSol(3);

                % SUM-MPR Method
                nAg = nAg + 1;
                [Th2, Ph2, g2, ~] = TDOA_SUM_MPR( senPos, rd_m, Q );
                uTh(i,nAg) = Th2;
                uPh(i,nAg) = Ph2;
                ug(i,nAg) = g2;

                % GTRS-MPR Method
                nAg = nAg + 1;
                [Th3, Ph3, g3, ~] = TDOA_GTRS_MPR( senPos, rd_m, Q );
                uTh(i,nAg) = Th3;
                uPh(i,nAg) = Ph3;
                ug(i,nAg) = g3;

            end

            % calculate MSE and bias
            for ia = 1:nAg
                mse_a(in,ir,ia) = mean((uTh(:,ia)-theta).^2+(uPh(:,ia)-phi).^2);
                mse_g(in,ir,ia) = mean((ug(:,ia)-g).^2);
                
                avBia_a(in,ir,ia) = sqrt(abs(mean(uTh(:,ia))-theta).^2+abs(mean(uPh(:,ia))-phi).^2);
                avBia_g(in,ir,ia) = abs(mean(ug(:,ia))-g);
            end
        end
    end

    symbs = ['o','v','s','*','^','+','x'];
    name = {'SCO-MPR','SUM-MPR','GTRS-MPR'};

    switch model
        case 'nse'
            xlabtext = '10log(\sigma^2(m^2))';
            xdata = nsePwr;
            yl_mse = [-80, 20;
                      -120, 20];
            yl_bias = [-155,3;
                       -250,0];
        case 'rag'
            xlabtext = 'Range(m)';
            xdata = souRange;
            yl_mse = [-40, 10;
                      -80,-45];
            yl_bias = [-80,0;
                       -150,-50];
    end

    %% plot results
    figure;
    for ia = 1:nAg
        plot(xdata, 10*log10(mse_a(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
    end
    plot(xdata, 10*log10(CRLB_a), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    ylim(yl_mse(1,:));
    xlabel(xlabtext, 'FontSize', 13);
    ylabel('10log(MSE(\theta,\phi)(rad^2))', 'FontSize', 13);
    lgd11 = legend('Show');
    set(lgd11, 'FontSize',11, 'Location', 'Northwest');
    set(gcf,'Position',[404 310 560 300]);

    figure;
    for ia = 1:nAg
        plot(xdata, 10*log10(mse_g(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
    end
    plot(xdata, 10*log10(CRLB_g), '-', 'LineWidth', 1.5, 'DisplayName', 'CRLB');
    ylim(yl_mse(2,:));
    xlabel(xlabtext, 'FontSize', 13);
    ylabel('10log(MSE(g)(1/m^2))', 'FontSize', 13);
    lgd2 = legend('Show');
    set(lgd2, 'FontSize',11, 'Location', 'Northwest');
    set(gcf,'Position',[404 310 560 300]);

    % Bias
    figure;
    for ia = 1:nAg
        plot(xdata, 20*log10(avBia_a(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
    end
    ylim(yl_bias(1,:));
    xlabel(xlabtext, 'FontSize', 13);
    ylabel('20log(Bias(\theta,\phi)(rad))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Northwest');
    xlim([min(xdata),max(xdata)]); %ylim([-170 10]);
    set(gcf,'Position',[404 310 560 300]);

    figure;
    for ia = 1:nAg
        plot(xdata, 20*log10(avBia_g(:,:,ia)), symbs(ia), 'LineWidth', 1.5, 'DisplayName', name{ia});hold on;grid on;
    end
    ylim(yl_bias(2,:));
    xlabel(xlabtext, 'FontSize', 13);
    ylabel('20log(Bias(g)(m^{-1}))', 'FontSize', 13);
    h3 = legend('Show');
    set(h3, 'FontSize',11, 'Location', 'Southeast');
    set(gcf,'Position',[404 310 560 300]);

    if im == 1
        clear mse_a mse_g CRLB_a CRLB_g avBia_a avBia_g;
    end
end