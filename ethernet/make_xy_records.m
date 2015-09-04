function recs = make_xy_records( as, bs )
    x1 = make_records(as(:,1), bs(:,1));
    y1 = make_records(as(:,2), bs(:,2));
    x2 = make_records(as(:,3), bs(:,3));
    y2 = make_records(as(:,4), bs(:,4));

    recs_raw = [x1, y1, x2, y2];
    recs = reshape(recs_raw',1,[]);
end

function r = make_records(a, b)
    number_of_recs = size(a,1);
    r = zeros(number_of_recs,28);
    r(:,3:4) = split_2_bytes(a);
    r(:,5:6) = split_2_bytes(b);
end