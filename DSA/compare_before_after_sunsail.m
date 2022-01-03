function compare_before_after_sunsail(sface) % Sface defines what surface to be plotted, Cone - Horizontal - Vertical
data1 = load("DSCN2430_data.mat", "IMG").IMG;
data2 = load("DSCN2440_data.mat", "IMG").IMG;
[X,Y] = meshgrid(1:365,(0:1/60:24-1/60));

switch sface
    case "cone"
        rel_numb = sum(data1.UVI_cone,"All") / sum(data2.UVI_cone,"All");
        figure(10)
        subplot(1,2,2)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data1.UVI_cone(1:10:end,1:5:end));hold on;
            %title({"UVI on a coned surface angled at 60 degrees", "Measures in place"})
            title("Sunsail installed")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")
            zlim([0.1, 2.5])

        subplot(1,2,1)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data2.UVI_cone(1:10:end,1:5:end));hold on;
            %title({"UVI on a coned surface angled at 60 degrees", "No measures made"})
            title("No form of sunscreening")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")
            zlim([0.1, 2.5])
        sgt = sgtitle(sprintf("Yearly reduction from UVI-reducing measures: %.2f %%", (1-rel_numb)*100),"Color", "red");
        sgt.FontSize = 15;
        
        figure(11);hold on;
        plot(Y(:,175), data1.UVI_cone(:,175))
        plot(Y(:,175), data2.UVI_cone(:,175))
        plot(Y(:,175), data1.UVI_cone_free(:,175))
        legend("With sunsail installed", "No sunscreen measures", "Theoretical clear sky")
        title({"UVI during the day at the 24th of June in Oslo (Kolbotn) on coned surface", sprintf("Yearly reduction from UVI-reducing measures: %.2f %%" , (1-rel_numb)*100)})
        ylabel("UV Index")
        xlabel("Hours in the day")
        % saveas(gcf, "Utsnitt_sthans.png")
        hold off;
        figure(12)
        imshowpair(data1.fish_rgb_RH_mirror, data2.fish_rgb_RH_mirror, "montage")
        
    case "horizontal"
        rel_numb = sum(data1.UVI_horizontal,"All") / sum(data2.UVI_horizontal,"All");
        figure(1)
        subplot(1,2,1)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data1.UVI_horizontal(1:10:end,1:5:end));hold on;
            title("UVI on a horizontal surface")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")

        subplot(1,2,2)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data2.UVI_horizontal(1:10:end,1:5:end));hold on;
            title("UVI on a horizontal surface")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")    
        sgt = sgtitle(sprintf("Total reduction (in %) from UVI-reducing measures: %.2f", (1-rel_numb)*100), "red");
        sgt.FontSize = 15;

        figure(11)
        imshowpair(data1.fish_rgb_RH_mirror, data2.fish_rgb_RH_mirror, "montage")
    case "vertical"
        rel_numb = sum(data1.UVI_vertical,"All") / sum(data2.UVI_vertical,"All");
        figure(1)
        subplot(1,2,1)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data1.UVI_vertical(1:10:end,1:5:end));hold on;
            title("UVI on a vertical surface")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")

        subplot(1,2,2)
        surf(X(1:10:end,1:5:end),Y(1:10:end,1:5:end),data2.UVI_vertical(1:10:end,1:5:end));hold on;
            title("UVI on a vertical surface")
            xlabel("Days in a year [days]")
            ylabel("Hours in a day [t]")
            zlabel("UVI")    
        sgt = sgtitle(sprintf("Total reduction (in %) from UVI-reducing measures: %.2f", (1-rel_numb)*100), "red");
        sgt.FontSize = 15;

        figure(11)
        imshowpair(data1.fish_rgb_RH_mirror, data2.fish_rgb_RH_mirror, "montage")
        
      end
    
end