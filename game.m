function [  ] = game( )

decks = 4;
players = 3;
cardsPlayed = [];

suit = 1:13;
cardsToPlay = ShuffleArray(repmat(suit,1,4*decks));


    function PlayHand()
        hands = Deal();
        for playerNo = 1:players
            v = input('hit = 1, stick = 2');
            while v < 2
                hands = Hit(playerNo, hands);
                points = HandPoints(hands{playerNo});
                if prod(points > 21) == 1
                    v = 2;
                end
            end
        end
    end

    function hands = Deal()
        cardsPlayed = [cardsPlayed, cardsToPlay(1:2*(players+1))];
        hands = mat2cell(cardsToPlay(1:2*(players+1)),1,ones(1,players+1)*2);
        cardsToPlay = cardsToPlay(2*(players+1)+1:end);
    end

    function hands = Hit(playerNo, hands)
        cardsPlayed = [cardsPlayed, cardsToPlay(1))];
        hands{playerNo} = [hands{playerNo} cardsToPlay(1)];
        cardsToPlay = cardsToPlay(2:end);
    end

    function oddsOfDealerPoints = FindOddsDealerPoints(cardsInHand)
        % returns [17 18 19 20 21]
        
        % run all possibilities to find most likely outcome
        % dealer stops on 17+
        
        % dealer may be holding cards or maybe not
        
    end
end

function shuf = ShuffleArray(arr)
    idx = randperm(length(arr));
    shuf = arr(idx);
end

function points = HandPoints(hand)
    points = sum(hand + (hand > 10) .* (10 - hand));
    points = [points ismember(1,hand).*(points + 10)];
end



