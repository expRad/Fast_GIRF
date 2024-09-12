classdef GIRF_plotter < handle
% This class defines a GIRF_plotter-object, which can create nice plots of the calculated GIRFs.

% Copyright (c) 2022 Hannah Scholten

    methods
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = GIRF_plotter()
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % constructor
        end
        
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plot_GSTFs(obj, H_1, name1, H_2, name2, H_3, name3, H_4, name4, numGSTFs, term2plot, ax_name)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot up to 4 GSTFs in the frequency and time domain
            figure('Units','normalized','Position',[0.05 0.02 0.5 0.85]);
            % First subplot: magnitude of the GSTFs in frequency domain
            subplot(3,1,1);
            hold on;
            grid on;
            plot(H_1.f_axis/1000, abs(H_1.gstf(:,term2plot)), 'DisplayName',[name1,' , df=',num2str(H_1.df),'Hz']);
            if numGSTFs>1
                plot(H_2.f_axis/1000, abs(H_2.gstf(:,term2plot)), 'DisplayName',[name2,' , df=',num2str(H_2.df),'Hz']);
                if numGSTFs>2
                    plot(H_3.f_axis/1000, abs(H_3.gstf(:,term2plot)), 'DisplayName',[name3,' , df=',num2str(H_3.df),'Hz']);
                    if numGSTFs>3
                        plot(H_4.f_axis/1000, abs(H_4.gstf(:,term2plot)), 'DisplayName',[name4,' , df=',num2str(H_4.df),'Hz']);
                    end
                end
            end
            xlabel('Frequency (kHz)');
            ylabel('Magnitude');
            if term2plot==1
                xlim([-4.9 4.9]);
                title(['Magnitude of GSTF term ', ax_name, '->B_0']);
            elseif (term2plot==2 && size(H_1.gstf,2)==2)
                axis([-50 50 0 1.5]);
                title(['Magnitude of GSTF term ', ax_name, '->', ax_name]);
            else
                xlim([-4.9 4.9]);
                title(['Magnitude of GSTF ', ax_name, ' (channel ', num2str(term2plot), ')']);
            end
            legend;

            % Second subplot: phase of the GSTFs in frequency domain
            subplot(3,1,2);
            hold on;
            grid on;
            plot(H_1.f_axis/1000, angle(H_1.gstf(:,term2plot)), 'DisplayName',[name1,' , df=',num2str(H_1.df),'Hz']);
            if numGSTFs>1
                plot(H_2.f_axis/1000, angle(H_2.gstf(:,term2plot)), 'DisplayName',[name2,' , df=',num2str(H_2.df),'Hz']);
                if numGSTFs>2
                    plot(H_3.f_axis/1000, angle(H_3.gstf(:,term2plot)), 'DisplayName',[name3,' , df=',num2str(H_3.df),'Hz']);
                    if numGSTFs>3
                        plot(H_4.f_axis/1000, angle(H_4.gstf(:,term2plot)), 'DisplayName',[name4,' , df=',num2str(H_4.df),'Hz']);
                    end
                end
            end
            xlabel('Frequency (kHz)');
            ylabel('Phase (rad)');
            if term2plot==1
                xlim([-4.9 4.9]);
                title(['Phase of GSTF term ', ax_name, '->B_0']);
            elseif (term2plot==2 && size(H_1.gstf,2)==2)
                axis([-50 50 -4 4]);
                title(['Phase of GSTF term ', ax_name, '->', ax_name]);
            else
                xlim([-4.9 4.9]);
                title(['Phase of GSTF ', ax_name, ' (channel ', num2str(term2plot), ')']);
            end
            legend;
            
            % Third subplot: time domain
            subplot(3,1,3);
            hold on;
            grid on;
            plot(H_1.t_axis*1000, H_1.girf(:,term2plot), 'DisplayName',[name1,' , dt=',num2str(H_1.dt),'s'],'visible','off');
            if numGSTFs>1
                plot(H_2.t_axis*1000, H_2.girf(:,term2plot), 'DisplayName',[name2,' , dt=',num2str(H_2.dt),'s']);
                if numGSTFs>2
                    plot(H_3.t_axis*1000, H_3.girf(:,term2plot), 'DisplayName',[name3,' , dt=',num2str(H_3.dt),'s']);
                    if numGSTFs>3
                        plot(H_4.t_axis*1000, H_4.girf(:,term2plot), 'DisplayName',[name4,' , dt=',num2str(H_4.dt),'s']);
                    end
                end
            end
            xlabel('Time (ms)');
            ylabel('GIRF');
            if term2plot==1
                title(['GIRF term ', ax_name, '->B_0']);
            elseif (term2plot==2 && size(H_1.gstf,2)==2)
                title(['GIRF term ', ax_name, '->', ax_name]);
            else
                title(['GIRF ', ax_name, ' (channel ', num2str(term2plot), ')']);
            end
            legend;
        end % function plot_GSTFs
        
                
    end % methods
    
end % classdef









