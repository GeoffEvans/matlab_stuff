eff = [0.019152542	0.010862069	0.068448276	0.232758621	0.322413793	0.355172414	0.267241379	0.21637931	0.186206897	0.168103448	0.028965517	0.006758621	0.032413793	; ... 
0.024830508	0.060948276	0.087068966	0.390517241	0.472413793	0.485344828	0.521551724	0.489655172	0.326724138	0.163793103	0.018965517	0.074396552	0.014827586	; ... 
0.017118644	0.062586207	0.062241379	0.279310345	0.63362069	0.913793103	0.512931034	0.235344828	0.143965517	0.131034483	0.067241379	0.134482759	0.017155172	; ... 
0.039491525	0.013793103	0.146551724	0.094827586	0.380172414	0.879310345	0.797413793	0.507758621	0.162068966	0.047155172	0.062586207	0.035517241	0.06362069	; ... 
0.01220339	0.053706897	0.05387931	0.10862069	0.148275862	0.623275862	0.922413793	0.719827586	0.254310345	0.073017241	0.156896552	0.00462069	0.005284483	; ... 
0.012627119	0.014482759	0.03	0.121551724	0.131034483	0.247413793	0.779310345	0.853448276	0.457758621	0.101724138	0.139655172	0.004663793	0.00525	; ... 
0.006059322	0.005336207	0.021034483	0.011293103	0.49137931	0.072844828	0.250862069	0.478448276	0.490517241	0.202586207	0.064827586	0.005681034	0.00312931	...  
];

freqs = [20 25 30 35 40 45 50];
angles = [0	0.6	1	1.4	1.6	1.8	2	2.2	2.4	2.6	3	3.4	4];
[anglesGrid,freqsGrid] = meshgrid(angles,freqs);

figure()
s = surf(anglesGrid,freqsGrid,eff);
colorbar;
