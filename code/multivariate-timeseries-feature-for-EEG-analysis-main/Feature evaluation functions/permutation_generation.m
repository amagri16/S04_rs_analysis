n = 4; % Number of possible values for each element (1 through 4)
k = 5; % Length of each permutation
totalPerms = n^k; % Total permutations with repetition, which is 4^6

% Preallocate matrix for efficiency
permsWithRep = zeros(totalPerms, k);

% Generate permutations with repetition
count = 1;
for a = 1:n
    for b = 1:n
        for c = 1:n
            for d = 1:n
                for e = 1:n
                    permsWithRep(count, :) = [a, b, c, d, e];
                    count = count + 1;
                end
            end
        end
    end
end

% Optionally shuffle rows to get random permutations
randomPermsWithRep = permsWithRep(randperm(size(permsWithRep, 1)), :);

% Display some of the permutations
disp(randomPermsWithRep(1:min(end, 10), :)); % Display the first 10 for brevity
