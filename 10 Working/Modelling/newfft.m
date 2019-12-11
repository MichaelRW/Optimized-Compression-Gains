function [ f, output ] = newfft( input, FS, varargin )

% Still have to work out the kinks~%

if (size(input,2) == 1)
  input = input(:)';
end

for i = 1:size(input,1)
  input(i,:) = input(i,:) - mean(input(i,:))/2;
end

if length(input) <= 2*FS
    N = 2*FS;
else
    N = FS*2^ceil(log2(length(input)/FS));
end

X = fft(input,N,2)/N;

% Remove half, and make 2*X(2:end) %
Pts = ceil((N+1)/2);
output=abs(4*X(:,1:Pts));
f=(0:Pts-1)*FS/N;

if (nargout == 0) | strncmpi(varargin,'y',1)
    figure;
    plot(gca,f,smooth(output));
    drawnow;
end