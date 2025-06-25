clearvars
close all

load grandaverage_HCP

% annealing variables
H = 1e04; hfrac = 10; Texp = 1-hfrac/H; T0 = 1; 
flag = 0;  % display cost
Tall = T0.*Texp.^[1:1:H];

% number of trials
T = 5;

% bag sizes (loop)
klvl = [3:30];  % bag size
K = length(klvl);

% minimize or maximize
minmax = 1;     % 1 for min, -1 for max

tic

COV = FC;

% initialize
minbagall = zeros(K,N,T);
inibagall = zeros(K,N,T);
mincostall = zeros(K,T);
inicostall = zeros(K,T);
%OIDhtall = zeros(K,3,H,T);

% loop over bag size
parfor kk=1:K
    
    k = klvl(kk);
    disp(k)
    
    mincostt = zeros(1,T);
    inicostt = zeros(1,T);
    minbagt = zeros(N,T);
    inibagt = zeros(N,T);
    %OIDht = zeros(3,H,T);
    
    % loop over trials
    for t=1:T
        
        h = 0; hcnt = 0;
        
        % initial bag of nodes
        rp = randperm(N);
        bag1 = rp(1:k);
        minbag = bag1;
        minbagglobal = minbag;
        [O1,I1,D1] = calcO_logdet(COV(minbag,minbag));
        %C1 = calcC_alt(COV(minbag,minbag));
        iniOID = [O1,I1,D1]; 
        OIDmin = iniOID; 
        OIDlow = iniOID;
        inicost = O1*minmax;
        mincost = inicost;
        lowcost = mincost;
        
        OIDh = zeros(3,H);
        
        while h<H
            
            h = h+1; hcnt = hcnt+1;
            
            % current temperature
            Tc = T0*Texp^h;
            
            if ((mod(h,H/10)==0)&&(flag==1))
                disp(['at step ',num2str(h),' - elapsed time ',num2str(toc),' - lowest cost = ',num2str(mincost)]);
            end
            
            % tweak minbag
            newbag = minbag;
            % swap how many?
            swap = min(k,ceil(abs(randn))); % no more than k
            %for s=1:swap
            candidates = setdiff(1:N,newbag);
            newbag(randperm(k,swap)) = candidates(randperm(length(candidates),swap));
            %newbag(ceil(rand*k)) = candidates(randi(length(candidates),1));
            %end
            
            [Onew,Inew,Dnew] = calcO_logdet(COV(newbag,newbag));
            OIDnew = [Onew,Inew,Dnew];
            costnew = Onew*minmax;
            
            % annealing
            randcon = rand < exp(-(costnew-lowcost)/Tc);
            %disp([num2str(costnew),' ',num2str(lowcost),' ',num2str(mincost),' ',num2str(costnew<lowcost),' ',num2str(randcon)]);
            if (costnew < lowcost) || randcon
                minbag = newbag;
                lowcost = costnew;
                OIDlow = OIDnew;
                % is this a new absolute best?
                if (lowcost < mincost)
                    minbagglobal = minbag;
                    mincost = lowcost;
                    OIDmin = OIDlow;
                    hcnt = 0;
                end
            end
            
            % time course of optimization
            %OIDh(1:3,h) = OIDmin;
            
        end
        
        mincostt(t) = mincost;
        inicostt(t) = inicost;
        temp = zeros(1,N); temp(minbagglobal) = 1; minbagt(:,t) = temp;
        temp = zeros(1,N); temp(bag1) = 1; inibagt(:,t) = temp;
        %OIDht(:,:,t) = OIDh;
        %disp([num2str(t),' ',num2str(min(mincostt(1:t)))])
        
    end
    
    % package
    minbagall(kk,:,:) = minbagt;
    inibagall(kk,:,:) = inibagt;
    mincostall(kk,:) = mincostt;
    inicostall(kk,:) = inicostt;
    %OIDhtall(kk,:,:,:) = OIDht;
    
end

toc

save anneal_grand_HCP_minO_sweep_5000 ...
    minbagall inibagall mincostall inicostall T H hfrac Texp T0 K klvl minmax FC N S lts T
