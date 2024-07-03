clear

files = dir;
nfiles = size(files,1);

for j = 1:nfiles
    display(['file ',num2str(j),' of ',num2str(nfiles)])
    if strfind(files(j).name,".mat")
        load(files(j).name,'iPTO');
        dataFilenames(iPTO,:) = {files(j).name};
    end
end

nPTO = size(dataFilenames,1);

%% Comparison to reference
D_w_YuJenne2017 = 0.5466; % [m^3/rad]
S_ro_YuJenne2017 = 2162; % [m^2]
q_permTotal_YuJenne2017 = (1476)/24/3600; % [m^3/day -> m^3/s]

D_w_PFF = 0.23; % [m^3/rad]
S_ro_PFF = 3700; % [m^2]
q_permTotal_PFF = 0.0175673200663364;%(500)/24/3600; % [m^3/day -> m^3/s]

switch 1
    case 1
        D_w_ref = D_w_YuJenne2017;
        S_ro_ref = S_ro_YuJenne2017;
        q_permTotal_ref = q_permTotal_YuJenne2017;
    case 2        
        D_w_ref = D_w_PFF;
        S_ro_ref = S_ro_PFF;
        q_permTotal_ref = q_permTotal_PFF;
end

% D_w_case1 = zeros(nPTO,1);
% S_ro_case2 = zeros(nPTO,1);
% q_permTotal_case3 = zeros(nPTO,1);

for iPTO = iiPTO
    display(['PTO ',num2str(iPTO),' of ',num2str(nPTO)])
    
    % load(char(dataFilenames(iPTO,:)),'data','D_wArray','S_roArray')
    
    % Case 1: Pump displacement based on equivalent membrane area and permeate
    % production
    try
        % find indices and value for 1 above target membrane area 
        S_ro_above = S_roArray(find(S_roArray > S_ro_ref,1,"first"));
        I_above = find(data(iPTO).S_ro == S_ro_above);
    
        % build array for permeate production 1 above target membrane area
        q_permTotal_above = data(iPTO).q_permTotal(I_above);
    
        % find indices and value for 1 below target membrane area
        S_ro_below = S_roArray(find(S_roArray < S_ro_ref,1,"last"));
        I_below = find(data(iPTO).S_ro == S_ro_below);
    
        % build array for permeate production 1 below target membrane area
        q_permTotal_below = data(iPTO).q_permTotal(I_below);
    
        % interpolate between curves above and below 
        x = (S_ro_ref - S_ro_below)/(S_ro_above - S_ro_below);
        q_permTotal_interp = q_permTotal_below ...
                             + x*(q_permTotal_above - q_permTotal_below); 
    
        % find index for lowest pump displacement 1 above target 
        % permeate production
        J_above = find(q_permTotal_interp > q_permTotal_ref,1,"first");
    
        % find index for lowest pump displacement 1 below target 
        % permeate production
        J_below = J_above-1;
    
        % interpolate between points 1 above and 1 below
        x = (q_permTotal_ref - q_permTotal_interp(J_below))...
            /(q_permTotal_interp(J_above) - q_permTotal_interp(J_below));
        D_w_case1(iPTO) = D_wArray(J_below) + ...
                          + x*(D_wArray(J_above) - D_wArray(J_below));
    
        figure
        plot(D_wArray,24*3600*q_permTotal_above)
        hold on
        plot(D_wArray,24*3600*q_permTotal_below)
        plot(D_wArray,24*3600*q_permTotal_interp)
        scatter(D_wArray(J_above),24*3600*q_permTotal_interp(J_above),'kx')
        scatter(D_wArray(J_below),24*3600*q_permTotal_interp(J_below),'kx')
        scatter(D_w_case1(iPTO),24*3600*q_permTotal_ref,'rx')
        title(['Membrane Area = ',...
            num2str(S_ro_ref) 'm^2: PTO ',num2str(iPTO)])
        xlabel('Displacement (m^3/rad)')
        ylabel('Permeate Production (m^3/day)')
        legend('data above','data below','interpolated',...
            'data above','data below','interpolated')

    catch
        D_w_case1(iPTO) = nan;
    end

    % Case 2: Membrane area based on equivalent pump displacement and permeate
    % production
    try
        % find indices and value for 1 above target pump displacement
        D_w_above = D_wArray(find(D_wArray > D_w_ref,1,"first"));
        I_above = find(data(iPTO).D_w == D_w_above);
    
        % build array for permeate production 1 above pump displacement
        q_permTotal_above = data(iPTO).q_permTotal(I_above);
    
        % find indices and value for 1 below target pump displacement
        D_w_below = D_wArray(find(D_wArray < D_w_ref,1,"last"));
        I_below = find(data(iPTO).D_w == D_w_below);
    
        % build array for permeate production 1 below pump displacement
        q_permTotal_below = data(iPTO).q_permTotal(I_below);
    
        % interpolate between curves above and below 
        x = (D_w_ref - D_w_below)/(D_w_above - D_w_below);
        q_permTotal_interp = q_permTotal_below ...
                             + x*(q_permTotal_above - q_permTotal_below); 
    
        % find index for lowest membrane area 1 above target permeate 
        % production
        J_above = find(q_permTotal_interp > q_permTotal_ref,1,"first");
    
        % find index for lowest membrane area 1 below target permeate 
        % production
        J_below = J_above-1;
    
        % interpolate between points 1 above and 1 below
        x = (q_permTotal_ref - q_permTotal_interp(J_below))...
            /(q_permTotal_interp(J_above) - q_permTotal_interp(J_below));
        S_ro_case2(iPTO) = S_roArray(J_below) + ...
                          + x*(S_roArray(J_above) - S_roArray(J_below));
    
        figure
        plot(1e-3*S_roArray,24*3600*q_permTotal_above)
        hold on
        plot(1e-3*S_roArray,24*3600*q_permTotal_below)
        plot(1e-3*S_roArray,24*3600*q_permTotal_interp)
        scatter(1e-3*S_roArray(J_above), ...
            24*3600*q_permTotal_interp(J_above),'kx')
        scatter(1e-3*S_roArray(J_below), ...
            24*3600*q_permTotal_interp(J_below),'kx')
        scatter(1e-3*S_ro_case2(iPTO),24*3600*q_permTotal_ref,'rx')
        title(['Displacement = ',...
            num2str(D_w_ref) 'm^3/rad: PTO ',num2str(iPTO)])
        xlabel('Membrane Area (1000 m^2)')
        ylabel('Permeate Production (m^3/day)')
        legend('data above','data below','interpolated',...
            'data above','data below','interpolated')
    catch
        S_ro_case2(iPTO) = nan;
    end

    % Case 3: Permeate production based on equivalent pump displacement and 
    % membrane area 
    x = data(iPTO).D_w(:);
    xq = D_w_ref;
    y = data(iPTO).S_ro(:);
    yq = S_ro_ref;
    v = data(iPTO).q_permTotal(:);
    q_permTotal_case3(iPTO) = griddata(x,y,v,xq,yq);

end