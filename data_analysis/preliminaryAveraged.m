eff = [0.020243902	0.095121951	0.074796748	0.056639566	0.060704607	0.151219512	0.894308943	0.951219512	0.845528455	0.470867209	0.451219512	0.161788618	0.086585366	0.043089431	0.053794038	0.165853659	; ... 
0.018943089	0.042411587	0.139837398	0.125203252	0.070817997	0.087804878	0.766666667	0.959349593	0.93495935	0.61076056	0.593495935	0.056097561	0.144715447	0.113821138	0.072346937	0.055005461	; ... 
0.01141791	0.031923614	0.021343284	0.071486875	0.142537313	0.130597015	0.478358209	0.669988472	0.880597015	0.955223881	0.835820896	0.094029851	0.09547385	0.164179104	0.164179104	0.041044776	...  
];

freqs = [30 35 40];
angles = [0	0.9	1	1.1	1.3	1.4	1.8	1.9	2	2.1	2.2	2.6	2.9	3	3	3.4];
[anglesGrid,freqsGrid] = meshgrid(angles,freqs);

figure()
s = surf(anglesGrid,freqsGrid,eff);
colorbar;
