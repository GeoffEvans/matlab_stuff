classdef AolConfig
    enumeration
        Orig_4
        Alt_4
        InterTrip
        CyclicPairs
    end
    methods
        function pairs = get_pairings(config)
            if config == AolConfig.Orig_4
                pairs = [1, 3; 2, 4];
            elseif config == AolConfig.Alt_4
                pairs = [1, 2; 3, 4];
            elseif config == AolConfig.InterTrip
                pairs = [1, 2; 3, 4; 5, 6];
            elseif config == AolConfig.CyclicPairs
                pairs = [1, 2; 3, 4; 5, 6];
            else
                error('unknown AOL config - is it 4 or 6 and in what order?')
            end
        end
        
        function aod_dirs = get_aod_dirs(config)
            if config == AolConfig.Orig_4
                aod_dirs = [1, 0, -1, 0; 0, 1, 0, -1];
            elseif config == AolConfig.Alt_4
                aod_dirs = [1, -1, 0, 0; 0, 0, 1, -1];
            else
                error('not implemented')
            end
        end
        
        function aod_dirs = get_num_aods(config)
            if config == AolConfig.Orig_4 || config == AolConfig.Alt_4
                aod_dirs = 4;
            else
                error('not implemented')
            end
        end
    end
end

