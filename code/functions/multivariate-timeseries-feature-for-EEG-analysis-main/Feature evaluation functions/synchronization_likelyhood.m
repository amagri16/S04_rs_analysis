function [RESULT, S] = synchronization_likelyhood(X)
% X should be a M x N array where M is number of samples and N is the number
% of channels 
[num_samples, num_chans] = size(X);

% initialize variables
m1 = 1;           % first sample
m2 = num_samples; % last sample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize parameters
[~,eLag,eDim] = phaseSpaceReconstruction(X,'MaxDim',num_chans*10,'MaxLag',round(0.8*num_samples/10));
lag  = eLag;    % 10;   % lag
m    = eDim;    % 10;   % embedding dimension
w1   = lag;     % 100;  % window (Theiler correction for autocorrelation)
pref = 0.05;    % 0.01;
w2   = round(10/pref+w1-1);   % 410;  % window (used to sharpen the time resolution of synchronization measure)

speed = 16;

% set active channels
usechan = zeros(num_chans,1); %FALSE = 0, TRUE = 1
usechan(1:num_chans) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% trim input data
num_usechan = sum(usechan(1:num_chans) ); % determine number of used channels
usechan_index = zeros(num_usechan,1); % create an index that relates # in-use channel to # actual channel

n = 1; % counter
for k = 1 : num_chans
    if usechan(k)
        usechan_index(n) = k;
        n = n + 1;
    end %if
end %for
X = X(:,usechan_index);

% calculate number of iterations
num_it = floor( (m2 - lag*(m-1))/speed ) - ceil( m1/speed ) + 1;


% calculate the synchronization likelihood matrix
[ S_matrix, hit_matrix ] = synchronization(X,lag,m,w1,w2,pref,speed);

%--------------------------------------------------
% calculate outputs
%--------------------------------------------------
% calculate S_ki for each channel & time, averaged over all other channels ("first file")
% i.e. an average synchronization value for driver system k, with all other response systems l
S_ki_matrix = zeros(num_usechan, num_it); % initialize
S_ki_temp = sum(S_matrix, 2); % sum across response systems (l)
S_ki_temp = (S_ki_temp - 1) / (num_usechan - 1); % average the sum ( -1 occurs to eliminate current channel count )
S_ki_matrix(:) = S_ki_temp(:); % store in matrix

%--------------------------------------------------
% calculate S_i for each time; is an averaged S_ki across k ("second file")
S_i_matrix = sum(S_ki_matrix,1)/num_usechan; % sum across k and average, size (1, num_it)

%--------------------------------------------------
% calculate pairwise time-averaged synchronization likelihood ("third file")

S_kl_temp = sum(hit_matrix, 3); % sum the hit matrix across time i, size (num_chan, num_chan)
hit_diag = diag( S_kl_temp );
% at a (k,l) position, s_kl_temp contains the number of hits occuring at both channels k & l, over all i & j
% at a (k,k) position, s_kl_temp contains the number of hits at channel k, over all i & j

S_kl_matrix = hit_diag * ones(1,num_usechan) + ones(num_usechan,1) * hit_diag';
% at a (k,l) position, s_kl_matrix contains the number of hits occuring at k and at l (hits at both k&l are counted twice)
% at a (k,k) position, s_kl_matrix contains the number of hits occuring at k, times 2

S_kl_matrix = S_kl_matrix + (S_kl_matrix == 0); % if S(k,k) == 0 & S(l,l) == 0 then S(k,l) must also be 0.
% this calculation protects against division by 0

% 2 * ( #k & l are both hit ) / ( #k is hit + #l is hit )
% = harmonic average of ( #k hit / # k & l are hit ) and ( #l hit / # k & l are hit )
S_kl_matrix = 2 * S_kl_temp ./ S_kl_matrix;

%--------------------------------------------------

% % overall synchronization
S = sum(sum(sum(S_matrix,1)-1))/(num_usechan-1)/num_usechan/num_it;
RESULT = S_kl_matrix;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYNCHRONIZATION SUBPROCEDURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a modified version of synchronization.m (display changes)

function [RESULT1, RESULT2] = synchronization(S_in,lag,m,w1,w2,pref,speed);

% read in data, determine size of data
% S_in must be 2-D, size #samples x #chans
[num_samples, num_chans] = size(S_in);

% initialize variables
m1 = 1;           % first sample
m2 = num_samples; % last sample

% calculate number of iterations
num_it = floor( (m2 - lag*(m-1))/speed ) - ceil( m1/speed ) + 1;

% initialize variables - inner loop
epsilon = ones(num_chans,1); % epsilon(of channel k) at time i

% initialize variables - outer loop
S_matrix = zeros(num_chans,num_chans,num_it); %S(k,l,i) matrix
hit_matrix = zeros(num_chans,num_chans,num_it);
i_count = 0; % iteration count, used to store matrix entries

% on_display percentage meter
pmdots = 20; % number of dots to display

for i = m1 : (m2 - lag*(m-1))
    if mod(i, speed) == 0
        
        i_count = i_count + 1;
        
        pm = floor(i_count/num_it*pmdots);
        
        % determine the valid j times, w1<|i-j|<w2
        j = m1 : (m2 - lag*(m-1));
        valid_range = abs(i-j)>w1 & abs(i-j)<w2; % vector of valid range positions
        num_validj = sum(valid_range); % number of valid range positions

        % construct compressed table of euclidean distances
        euclid4_table = zeros(num_chans,num_validj);
        n = 0; % counter
        for j = m1 : (m2 - lag*(m-1))
            if valid_range(j-m1+1)
                n = n + 1;
                for k = 1 : num_chans
                    % euclid4 not explicity called; saves about 25% of time
                    %euclid4_table(k,n) = euclid4(S_in,lag,m,k,i,j);
                    euclid4_table(k,n) = sqrt( sum(   (S_in(i+lag*(0:(m-1)),k) ...
                                                     - S_in(j+lag*(0:(m-1)),k)).^2   ) );
                end %for
            end %if
        end %for
        
        % construct table of epsilons (formerly used the "crlocal" subroutine)
        % epsilon(k) is epsilon_{k,i}: the actual threshold distance such that the fraction
        % of all distances |X_{k,i} - X_{k,j}| less than epsilon_{k,i} is Pref
        for k = 1 : num_chans
            sorted_table = sort( euclid4_table(k,:) ); % size (1,validj)
            epsilon(k) = sorted_table( ceil( pref * num_validj ) );
        end %for
        
        % construct 'hit' table, i.e. determine if |X_{k,i} - X_{k,j}| <= epsilon_x for each k & j
        % size (num_chans, num_validj)
        hit_table = ( euclid4_table <= ( epsilon(1:num_chans) * ones(1,num_validj) ) );
        hit_table = double(hit_table); %Matlab 6.5
        
        % construct alternate hit table:
        % at position (k,l), determine the number of hits occuring at both channels k & l (across all j)
        % size (num_chans, num_chans)
        hit_table2 = hit_table * hit_table';
        
        % determine number of hits for each channel, across all j
        % NOTE: this is equivalent to diag(hit_table2)
        %num_hitsperchan = sum( hit_table, 2 ); % size (num_chans,1)
        num_hitsperchan = diag(hit_table2);
        
        % store hit_table2 in a 3-D array
        hit_matrix(:,:,i_count) = hit_table2;
        
        % perform conditional probability calculation
        % divide k^th row by number of hits for channel k
        S_matrix(:,:,i_count) = hit_table2 ./ ( num_hitsperchan * ones(1,num_chans) );
        
    end %if mod(i, speed) == 0
end %for

RESULT1 = S_matrix;
RESULT2 = hit_matrix;
  
% % overall synchronization likelihood
% S = sum(sum(sum(S_matrix,1)-1))/(num_chans-1)/num_chans/num_it;