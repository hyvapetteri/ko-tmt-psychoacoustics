classdef RoundRobinTournament < TournamentInterface
    %ROUNDROBINTOURNAMENT A class for running a knockout-type tournament of 
    % n_teams, until K best teams have been found.
    %
    %   Winner of each "game" is decided by choosing the team that has
    %   won wins number of games. For example, in a best-of-three, the
    %   number of wins required is two, and for a single game, number
    %   of wins is one.


    properties
        N % total number of teams in the tournament
        max_games % number of games in total
        reps % n of times each pair of team meets

        % tmt is the tournament tree as a cell array, where each cell of
        % the array corresponds to a tournament round. Each tournament
        % round cell contains an array of team numbers, and adjacent teams
        % in the array are paired for comparison. There are always an even
        % number of teams in each round.
        tmt

        % res has the same structure as tmt, but instead of team numbers,
        % it contains the scores for each game.
        res

        game % current game in round

        finished % has the tournament finished?

        name = 'RoundRobin' % name of the round

    end

    methods
        function obj = RoundRobinTournament(n_teams, reps)
            obj.N = n_teams;
            obj.reps = reps;
        end

        function obj = randomInit(obj,randomSeed, leaveOut)
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

            if nargin < 3
                leaveOut = [];
            end

            % randomize the game order
            obj.tmt = repmat(nchoosek(1:obj.N, 2), obj.reps, 1);

            if ~isempty(leaveOut)
                % leave out certain games
                %[~,idx] = ismember(sort(leaveOut,2), sort(obj.tmt,2), 'rows');
                
                [leaveOutBool, ~] = ismember(sort(obj.tmt,2), sort(leaveOut,2),  'rows');
                
                obj.tmt(leaveOutBool,:) = [];
            end

            obj.max_games = size(obj.tmt,1);
            obj.tmt = obj.tmt(randperm(obj.max_games),:);

            obj.res = zeros(size(obj.tmt));
      
            obj.game = 1;
            obj.finished = false;
        end

        function pair = getCurrentGame(obj)
            %GETCURRENTGAME Returns a 2x1 array containing the team numbers
            %competing in the current game.
            if ~obj.finished
                pair = obj.tmt(obj.game,:);
            else
                pair = [NaN NaN];
            end
        end

        function n_games = getNCompletedGames(obj)
            %GETNCOMPLETEDGAMES returns the number of games that have been
            %decided in the tournament so far. A game is considered as
            %completed, when either team has won at least WINS number of
            %matches

            n_games = obj.game - 1;
        end

        function obj = update(obj, result)
            %UPDATE Updates the tournament tree according to the result of
            % the current game. The input to UPDATE has to be either 1 or 2 
            % (which team won the current game, corresponding to idx of 
            % getCurrentGame pair)
            
            obj.res(obj.game, result) = 1;

            if obj.game == obj.max_games
                obj.finished = true;
            end
            
            obj.game = obj.game + 1;

        end % <- end of function update

        function printStatus(obj)

        end

    end
end