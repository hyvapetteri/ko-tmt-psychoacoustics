classdef KnockoutTournament < TournamentInterface
    %KNOCKOUTTOURNAMENT A class for running a knockout-type tournament of 
    % n_teams, until K best teams have been found.
    %
    %   Winner of each "game" is decided by choosing the team that has
    %   won wins number of games. For example, in a best-of-three, the
    %   number of wins required is two, and for a single game, number
    %   of wins is one.
    %
    %   Note, that the resulting K teams are an unordered set; the
    %   first team on the list may very well not be the best.
    %
    %   The tournament is played in two phases: first as a regular
    %   knock-out tournament, and then after the first winner is found,
    %   the same tournament is replayed by replacing the winner by the
    %   next team in teams_leave_out. The logic here is, that the
    %   second-best team in the first phase of the tournament is one of the
    %   teams that lost to the winner, and by replaying the tournament
    %   without the winner, the second-best should win —- no matter whether
    %   they lost to the winner in the first game or the final. However, if
    %   the team replacing the original winner is better than the original
    %   second-best, then that team will win the replayed tournament. There
    %   will always be top_K - 2 teams in teams_leave_out. After the last
    %   team from teams_leave_out has played in the tournament, the last
    %   tournament will run by replacing the last winner by an infinitely
    %   poor team, i.e. that the first opponent of the last winner will win
    %   the first game always.

    properties
        tmt_algo % either 1 for treesort or 2 for replacement selection
        N % total number of teams in the tournament
        n_rounds % number of rounds in the tournament tree
        max_games % maximum number of games in total
        games_played % games played so far in total. Updated after each full round
        top_K % how many top teams to pick
        k % the current k out of top_K
        wins % n of wins required (e.g. best of 3 -> 2 wins)
        idxs % the list of teams in random order

        % team_pos gives the starting location of a team in the tournament.
        % For example team_pos(2) will give the round and game numbers
        % where team 2 is seeded. The first row of team_pos contains the
        % starting round of each team, and second row contains the index
        % of the team in tmt{round}
        team_pos

        % tmt is the tournament tree as a cell array, where each cell of
        % the array corresponds to a tournament round. Each tournament
        % round cell contains an array of team numbers, and adjacent teams
        % in the array are paired for comparison. There are always an even
        % number of teams in each round.
        tmt

        % res has the same structure as tmt, but instead of team numbers,
        % it contains the scores for each round.
        res
        teams_leave_out % teams left out of the first phase

        round % current round in tmt (-> tmt{round})
        game % current game in round

        replay % replay is true if we are not anymore in the first phase of the tournament
        preround % is there a pre-round for getting the number of games to 2^n?
        finished % has the tournament finished?
        n_round_games % number of games in the current round
        k_winners % the final top_K winners

        name = 'Knockout' % name of the round

    end

    methods
        function obj = KnockoutTournament(n_teams, K, wins, algorithm)
            obj.N = n_teams;
            obj.top_K = K;
            obj.wins = wins;
            obj.tmt_algo = algorithm;
            
            if algorithm == 1
                obj.n_rounds = ceil(log2(obj.N));
                obj.max_games = n_teams - K + sum(arrayfun(@(x) ceil(log2(x)), (n_teams + 1 - K + 1):n_teams));
            elseif algorithm == 2
                % top_K - 2 teams are left out from the first phase of the
                % tournament
                obj.n_rounds = ceil(log2(obj.N - obj.top_K + 2));
    
                obj.max_games = n_teams - K + (K-1)*ceil(log2(n_teams + 2 - K));
            else
                error('Value of property algorithm must be either 1 (treesort) or 2 (replacement selection)');
            end

            obj.games_played = 0;
        end

        function obj = randomInit(obj,randomSeed)
            %RANDOMINIT Initialize the tournament tree with random
            %allocation of teams
            %   randomSeed can be provided as integer, which will be used
            %   for the random number generator, if repeatable results are
            %   desired.
            if nargin > 1
                rng(randomSeed);
            else
                rng('shuffle');
            end

            % randomize the team pairings
            obj.idxs = randperm(obj.N);
            [~, I] = sort(obj.idxs);
            obj.team_pos = repmat(I, 2, 1);
            obj.team_pos(1,:) = 0;

            n_leave_out = obj.top_K - 2;
            if (obj.top_K == 1) || (obj.tmt_algo == 1)
                n_leave_out = 0;
            end
            % we want a tournament size of N = 2^t, but if N is not a power
            % of 2, such that N = 2^t + k, 0 <= k < 2^t, then we first pick
            % 2k teams to play a preround of k games, and the N - 2k other
            % teams don't play in the preround. After the preround, we have
            % 2^t players left for a balanced knockout tournament. We use
            % n_extra variable name instead of k, since k is already used
            % in another meaning.
            balanced_n = 2^floor(log2(obj.N - n_leave_out));
            n_extra = obj.N - n_leave_out - balanced_n;

            obj.tmt = {};
            obj.res = {};

            if (n_extra > 0)
                % If k in N = 2^t + k is > 0

                obj.preround = true;

                % only 2k teams compete in the preround
                obj.tmt{1} = obj.idxs(1:2*n_extra);
                obj.res{1} = zeros(size(obj.tmt{1}));

                % number of games is always number of teams / 2
                n_games = ceil(length(obj.tmt{1})/2);
                % the N - 2k teams not playing in preround start at the
                % second round
                obj.tmt{2} = [zeros(1, n_games) obj.idxs((2*n_games + 1):(end-n_leave_out))];
                obj.res{2} = zeros(size(obj.tmt{2}));
                obj.n_round_games = n_games;
                % initialize also third round data structures
                obj.tmt{3} = zeros(1,ceil(length(obj.tmt{2})/2));
                obj.res{3} = zeros(size(obj.tmt{3}));

                % teams in preround are seeded in round 1
                obj.team_pos(1, obj.idxs(1:2*n_extra)) = 1;
                % teams not in preround are seeded in round 2
                obj.team_pos(1, obj.idxs((2*n_extra + 1):(end-n_leave_out))) = 2;
                % the position of 2nd round seeded teams
                obj.team_pos(2, obj.idxs((2*n_extra + 1):(end-n_leave_out))) = ...
                    obj.team_pos(2, obj.idxs((2*n_extra + 1):(end-n_leave_out))) - n_games;
            else
                % the tournament size is already a power of 2, so no need
                % for a preround; every team starts at round 1.

                obj.preround = false;

                obj.tmt{1} = obj.idxs(1:(end-n_leave_out));
                n_games = ceil(length(obj.tmt{1})/2);
                % initialize also second round data structures
                obj.tmt{2} = zeros(1,n_games);
                obj.n_round_games = n_games;
                obj.team_pos(1,:) = 1;
                obj.res{1} = zeros(size(obj.tmt{1}));
                obj.res{2} = zeros(size(obj.tmt{2}));
            end
            
            % teams left out of the first phase
            obj.teams_leave_out = obj.idxs((end-n_leave_out+1):end);

            obj.round = 1;
            obj.game = 1;
            obj.k = 1;
            obj.k_winners = [];
            obj.replay = false;
            obj.finished = false;
        end

        function pair = getCurrentGame(obj)
            %GETCURRENTGAME Returns a 2x1 array containing the team numbers
            %competing in the current game.
            if ~obj.finished
                teams = obj.tmt{obj.round};
                pair = [teams(2*obj.game - 1), teams(2*obj.game)];
            else
                pair = [NaN NaN];
            end
        end

        function n_games = getNCompletedGames(obj)
            %GETNCOMPLETEDGAMES returns the number of games that have been
            %decided in the tournament so far. A game is considered as
            %completed, when either team has won at least WINS number of
            %matches

            % get the number of games played up to the current round
            n_games = obj.games_played;
            % add to that the number of decided games on the current round
            scores = obj.res{obj.round};
            n_decided = sum(scores == obj.wins);
            n_games = n_games + n_decided;
        end

        function obj = update(obj, result)
            %UPDATE Updates the tournament tree according to the result of
            % the current game. The input to UPDATE has to be either 1 or 2 
            % (which team won the current game, corresponding to idx of 
            % getCurrentGame pair)

            keep_looking = true;
                
            while keep_looking
                scores = obj.res{obj.round};
                % teams competing in game K have indices 2*K - 1 (team 1) and
                % 2*K (team 2). Result is either 1 or 2, so if team 2 won,
                % winner_idx will be 2*K - 2 + 2 = 2*K
                winner_idx = 2*obj.game - 2 + result;
                scores(winner_idx) = scores(winner_idx) + 1;
                obj.res{obj.round} = scores;
    
                if any(scores > obj.wins)
                    error('Too many wins, something went wrong.');
                end
    
                teams = obj.tmt{obj.round};
                winners = zeros(obj.n_round_games,1);
    
                if ~obj.replay
                    for i=1:obj.n_round_games
                        s1 = scores(2*i - 1);
                        s2 = scores(2*i);
                        if (s1 == obj.wins)
                            winners(i) = teams(2*i - 1);
                        elseif (s2 == obj.wins)
                            winners(i) = teams(2*i);
                        end
                    end
                    % find games where there are less than obj.wins wins
                    unsettled_games = find(winners == 0);
                    % if all games have a winner, move to next round
                    round_finished = isempty(unsettled_games);
    
                else
                    % in the second phase, we only replay the previous winner's
                    % path, and thus there is only one game per round. We play
                    % the game until either team has obj.wins number of wins 
                    round_finished = false;
                    unsettled_games = obj.game;
    
                    s1 = scores(2*obj.game - 1);
                    s2 = scores(2*obj.game);
                    if (s1 == obj.wins)
                        winners(1) = teams(2*obj.game - 1);
                        round_finished = true;
                    elseif (s2 == obj.wins)
                        winners(1) = teams(2*obj.game);
                        round_finished = true;
                    end
                end
    
    
                if round_finished
                    if (teams(2*obj.game) > -1) && (teams(2*obj.game - 1) > -1)
                        obj.games_played = obj.games_played + sum(obj.res{obj.round}, 'omitnan');
                    end
                    % clear the scores for previous round
                    obj.res{obj.round} = nan(size(obj.res{obj.round}));
                    
                    if ~obj.replay && (obj.round < obj.n_rounds)
                        % we are in the first phase of the tournament, and the
                        % round that just finished, was not the final. For the
                        % next round, start at first game.
                        obj.game = 1;
                        obj.round = obj.round + 1;
    
                        if obj.preround && (obj.round == 2)
                            % if the round that just finished was the preround,
                            % update the next round with the winners of the
                            % preround, but keep other teams intact.
                            teams = obj.tmt{2};
                            teams(1:obj.n_round_games) = winners;
                            obj.tmt{2} = teams;
    
                            n_games = ceil(length(teams)/2);
                            obj.n_round_games = n_games;
                        else
                            % if the previous round was not the preround, then
                            % the teams competing on the next round are only
                            % the winners of the round that just finished
                            obj.tmt{obj.round} = winners;
                            n_games = ceil(length(winners)/2);
                            obj.n_round_games = n_games;
                            obj.tmt{obj.round + 1} = zeros(1,n_games);
                            obj.res{obj.round + 1} = zeros(1,n_games);
                        end
                        
    
                    elseif obj.replay && (obj.round < obj.n_rounds)
                        % if we are in the second phase, and not in the final,
                        % just move the winner of the current round to the next
                        % round, following the earlier winner's path
                        teams = obj.tmt{obj.round + 1};
                        teams(obj.game) = winners(1);
                        obj.tmt{obj.round + 1} = teams;
    
                        obj.game = ceil(obj.game/2);
                        obj.round = obj.round + 1;
    
                        % reset the score
                        scores = obj.res{obj.round};
                        scores(2*obj.game) = 0;
                        scores(2*obj.game - 1) = 0;
                        obj.res{obj.round} = scores;
    
                    elseif obj.round == obj.n_rounds
                        % the round that just finished was the final game.
                        % Thus, the winner of the round is one of the top_K
                        % teams
                        w = winners(1);
                        obj.k_winners(obj.k) = w;
    
                        if (obj.tmt_algo == 2) && (obj.k < obj.top_K - 1)
                            % there are still teams that where left out of the
                            % first phase, so replace the winner with the next
                            % one in line.
                            pos = obj.team_pos(:,w); % starting position of the current winner
                            obj.round = pos(1);
                            obj.game = ceil(pos(2)/2);
    
                            teams = obj.tmt{pos(1)};
                            % replace the winner by the next one in line
                            teams(pos(2)) = obj.teams_leave_out(obj.k);
                            % update team_pos accordingly
                            obj.team_pos(:, obj.teams_leave_out(obj.k)) = pos;
                            obj.tmt{pos(1)} = teams;
    
                            % reset score
                            scores = obj.res{pos(1)};
                            scores(2*obj.game) = 0;
                            scores(2*obj.game - 1) = 0;
                            obj.res{pos(1)} = scores;
    
                            % just play one game per round, since
                            % we are replaying the winner's path
                            obj.n_round_games = 1;
                            obj.replay = true;
    
                        elseif ((obj.tmt_algo == 1) && (obj.k < obj.top_K - 1)) || (obj.k == obj.top_K - 1)
                            % there are no more teams waiting. The winner is
                            % now replaced by an "infinitely bad" team, meaning
                            % that the opponent of the first game of the
                            % current winner will automatically win and move to
                            % the next round.
                            pos = obj.team_pos(:,w);
                            teams = obj.tmt{pos(1)};
                            teams(pos(2)) = -1; % "infinitely bad team"
                            obj.tmt{pos(1)} = teams;
                            
                            % restart the tournament, playing the winner's
                            % path
                            obj.round = pos(1);
                            obj.game = ceil(pos(2)/2);
                            
                            % reset the score
                            scores = obj.res{obj.round};
                            scores(2*obj.game) = 0;
                            scores(2*obj.game - 1) = 0;
                            obj.res{obj.round} = scores;
    
                            obj.n_round_games = 1;
                            obj.replay = true;
                        else
                            % the tournament has finished, no more games to
                            % play
                            obj.game = -1;
                            obj.round = -1;
                            obj.n_round_games = -1;
                            obj.finished = true;
                        end
    
                        obj.k = obj.k + 1;
                    else
                        % should not happen
                        error('unexpected control logic situation.');
                    end
                else
                    % this is the else to: if round_finished
                    % so either just pick a random game among the unfinished
                    % games in the current round, or if there is just one 
                    % unsettled game left, take that.
                    if length(unsettled_games) > 1
                        obj.game = randsample(unsettled_games, 1);
                    else
                        obj.game = unsettled_games;
                    end
                end

                % check if the opponent of the next match is "infinitely 
                % bad", i.e. if we should just move directly to the next 
                % game
                p = obj.getCurrentGame();
                if any(p == -1)
                    wi = 0;
                    scores = obj.res{obj.round};

                    % one of the opponents is "infinitely bad", so we make
                    % sure the other team advances immediately. We will
                    % re-run the loop inside the update function until we
                    % end up with a game where both teams are still
                    % playing. So, variable keep_looking will stay true
                    if (p(1) == -1)
                        wi = 2;
                        scores(2*obj.game - 1) = 0;
                        scores(2*obj.game) = obj.wins - 1;
                    elseif (p(2) == -1)
                        wi = 1;
                        scores(2*obj.game - 1) = obj.wins - 1;
                        scores(2*obj.game) = 0;
                    end
                    obj.res{obj.round} = scores;
                    result = wi;
                else
                    % neither team has been replaced by the "infinitely
                    % bad" team, so stop looking
                    keep_looking = false;
                end
            end

        end % <- end of function update

        function printStatus(obj)
            fprintf('Top-%d Tournament with %d teams, %d rounds\n', obj.top_K, obj.N, obj.n_rounds);
            fprintf('Winners: ');
            for i=1:length(obj.k_winners)
                fprintf('%3d  ', obj.k_winners(i));
            end
            fprintf('\n');
            fprintf('Current tournament tree:\n');
            fprintf('Teams waiting: ');
            for i=obj.k:(obj.top_K - 2)
                fprintf('%d ', obj.teams_leave_out(i));
            end
            fprintf('\n');
            for i=1:min(length(obj.tmt), obj.n_rounds)
                fprintf('round %d:\t', i);
                teams = obj.tmt{i};
                scores = obj.res{i};
                whitespace = repmat(' ', 1, 2*2^i - 3);
                for j=1:length(teams)
                    fprintf('%3d', teams(j));
                    fprintf(whitespace);
                end
                fprintf('\n');

                fprintf('scores :\t');
                for j=1:length(scores)
                    fprintf('%3d', scores(j));
                    fprintf(whitespace);
                end
                fprintf('\n\n');
            end
            fprintf('Current round: %d, game: %d\n', obj.round, obj.game);
        end

    end
end