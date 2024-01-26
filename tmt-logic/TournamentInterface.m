classdef TournamentInterface < handle
    % Abstract class for creating tournaments of paired comparisons

    methods (Abstract)
        randomInit(obj, randomSeed)
        % Initialize the tournament with a random placement of teams.
        % RANDOMSEED can be provided as integer, which will be used
        % for the random number generator, if repeatable results are
        % desired.

        pair = getCurrentGame(obj)
        %GETCURRENTGAME Returns a 2x1 array containing the team numbers
        %competing in the current game.
        
        n_games = getNCompletedGames(obj)
        %GETNCOMPLETEDGAMES returns the number of games that have been
        %decided in the tournament so far.

        update(obj, result)
        %UPDATE Updates the tournament tree according to the result of
        % the current game. The input to UPDATE has to be either 1 or 2 
        % (which team won the current game, corresponding to idx of 
        % getCurrentGame pair)

        printStatus(obj)
    end
end