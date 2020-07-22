function plot_electrodes(type)
if type == 1
    for dx = -90:60:90
        plot([dx dx],[-90 90],'k-')
        plot([-90 90],[dx dx],'k-')
    end
elseif type == 2
    for dx = -92.5:10:92.5
        for dy = -92.5:10:92.5
            rectangle('position',[dx,dy,5,5],...
                'FaceColor',[1 1 0 0.3],...
                'EdgeColor',[1 1 0 0.3])
        end
    end
end
end