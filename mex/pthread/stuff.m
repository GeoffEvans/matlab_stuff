data = rand(5e7,1);
myPosixThread('test.data',data);  % start running in parallel
%data2 = fft(data);  % post-processing (pthread I/O runs in parallel)
