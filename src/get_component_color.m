function clr = get_component_color(K)
    clr0 = orderedcolors("gem12");
    nClr = size(clr0,1);
    if K>nClr
        nsubclr = ceil(K/(nClr-1));
        clr = [];
        for c = 1:nClr-1
            red = linspace(clr0(c,1),clr0(c+1,1),nsubclr);
            green = linspace(clr0(c,2),clr0(c+1,2),nsubclr);
            blue = linspace(clr0(c,3),clr0(c+1,3),nsubclr);
            clrsc = [red',green',blue'];
            for sc = 1:nsubclr
                clr = cat(1,clr,clrsc(sc,:));
                if size(clr,1)==K
                    break
                end
            end
            if size(clr,1)==K
                break
            end
        end
    else
        clr = clr0(1:K,:);
    end
end