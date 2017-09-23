function phase = NonlocalUnwrap(handles)

I = handles.phase;

THRESH = 1*pi;
[Ny, Nx] = size(I);

groups = -1*ones(size(I));

H =  circshift(I, [-1, 0]) - 2*I - circshift(I, [1, 0]);
V =  circshift(I, [0, -1]) - 2*I - circshift(I, [0, 1]);
D1 = circshift(I, [-1, -1]) - 2*I - circshift(I, [1, 1]);
D2 = circshift(I, [-1, 1]) - 2*I - circshift(I, [1, -1]);
D = H.^2 + V.^2 + D1.^2 + D2.^2;

R = 1./D;

edges = -1*ones(length(I(:)),4);

for i = 2:Ny-1
    for j = 2:Nx-1
        n = sub2ind(size(I), i, j);
        edges(n, 1) = R(i,j) + R(i+1,j);
        edges(n, 2) = R(i,j) + R(i-1,j);
        edges(n, 3) = R(i,j) + R(i,j+1);
        edges(n, 4) = R(i,j) + R(i,j-1);
    end
end

% eind(1) contains most reliable edge
[~, eind] = sort(-edges(:));

Numgoodedges = length(edges(edges > 0));

wbar = waitbar(0,'Traversing edges...');

% Some edges we ignore, i.e. the corners of the image
for n = 1:Numgoodedges
    ind = eind(n);
    [ie, je] = ind2sub(size(edges), ind);
    [isource, jsource] = ind2sub(size(I), ie);
    if je == 1
        isink = isource+1; jsink = jsource;
    end
    if je == 2
        isink = isource-1; jsink = jsource;
    end
    if je == 3
        isink = isource; jsink = jsource+1;
    end
    if je == 4
        isink = isource; jsink = jsource-1;
    end
    
    try
        
        g1 = groups(isource, jsource);
        g2 = groups(isink, jsink);
        
        I1 = I(isource, jsource);
        I2 = I(isink, jsink);
        
        
        % Not unwrapped yet
        if g1 == -1 && g2 == -1
            while abs(I1 - I2) > THRESH
                I2 = I2 + 2*pi*sign(I1-I2);
            end
            
            I(isource, jsource) = I1;
            I(isink, jsink) = I2;
            
            gmax = max(groups(:));
            if gmax < 0
                gmax = 0;
            end
            groups(isource, jsource) = gmax+1;
            groups(isink, jsink) = gmax+1;
        end
        
        % g2 not in group
        if g1 > 0 && g2 == -1
            while abs(I1 - I2) > THRESH
                I2 = I2 + 2*pi*sign(I1-I2);
            end
            
            I(isink, jsink) = I2;
            groups(isink, jsink) = g1;
        end
        
        % g1 not in group
        if g1 == -1 && g2 > 0
            while abs(I1 - I2) > THRESH
                I1 = I1 + 2*pi*sign(I2-I1);
            end
            
            I(isource, jsource) = I1;
            groups(isource, jsource) = g2;
        end
        
        % Pixels in different groups
        if g1 > 0 && g2 > 0 && g1 ~= g2
            g1size = length(groups(groups == g1));
            g2size = length(groups(groups == g2));
            
            % Merge g2 into g1
            if g1size > g2size
                while abs(I1 - I2) > THRESH
                    I(groups == g2) = I(groups == g2) + 2*pi*sign(I1-I2);
                    I2 = I2 + 2*pi*sign(I1-I2);
                end
                groups(groups == g2) = g1;
            end
            
            % Merge g1 into g2
            if g1size < g2size
                while abs(I1 - I2) > THRESH
                    I(groups == g1) = I(groups == g1) + 2*pi*sign(I2-I1);
                    I1 = I1 + 2*pi*sign(I2-I1);
                end
                groups(groups == g1) = g2;
            end
            
        end
        
    catch
        % Couldn't access something? TODO: Check this doesn't happen
    end
    
    if ~mod(n,100000)
        waitbar(n/Numgoodedges)
        n
    end
end

close(wbar);

% Add corners manually
I(1,1) = 0.5*(I(1,2)+I(2,1));
I(1,end) = 0.5*(I(1,end-1)+I(2,end));
I(end,1) = 0.5*(I(end,2)+I(end-1,1));
I(end,end) = 0.5*(I(end,end-1)+I(end-1,end));

phase = I;

end