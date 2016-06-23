one = randi(all,1,1);

t = wav{one,S1S2}{2,1};
t = cell2mat(t);
d = diff(t);

plot(d);