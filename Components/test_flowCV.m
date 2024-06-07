deltap = linspace(0,1e6,1e3);
p_c = 1e5;
kv = 1;
dp = 1e5;

q = flowCV(deltap,kv,p_c,dp);

figure
plot(deltap,q)