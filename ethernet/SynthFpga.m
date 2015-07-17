classdef SynthFpga < handle

    properties
        dest_addr = hex2dec(['aa'; 'bb'; 'cc'; 'dd'; 'ee'; 'ff']);
        send_addr
    end
    
    methods
        function this = SynthFpga()
            [~,result] = dos('getmac');
            this.send_addr = result(160:176);
        end
        
        function load(this)
            [a, b, ticksPerRamp] = integration_with_driver();
            drive_coeffs = struct('a',a,'b',b,'t',ticksPerRamp);
            this.load_plane_size();
            this.load_plane_records_queue(drive_coeffs);
        end
        
        function load_plane_size(this)
            length = 6;
            data = [split_2_bytes(length) SynthComs.load_plane_size.v 0 100 0 100 0];
            send_packet(this.dest_addr, this.send_addr, data);
        end
        
        function load_plane_records_queue(this, drive_coeffs)
            start_index = 0;
            num_recs = length(drive_coeffs);
            while start_index < num_recs
                end_index = start_index + 13;
                if end_index > num_recs
                    end_index = num_recs; % catch the end
                end
                this.load_plane_records(drive_coeffs(start_index:end_index), start_index, end_index);
                start_index = start_index + 14;
            end
        end
        
        function load_plane_records(this, drive_coeffs, start_index, end_index)
            recs = make_xy_records(drive_coeffs.a, drive_coeffs.b);
            data = [0 0 SynthComs.load_plane_records.v 0 0 100 0 100 0 split_2_bytes(start_index) split_2_bytes(end_index) reshape(recs',1,[])];
            data_length = length(data);
            data(1:2) = split_2_bytes(data_length);
            send_packet(this.dest_addr, this.send_addr, data);
        end
    end    
end

