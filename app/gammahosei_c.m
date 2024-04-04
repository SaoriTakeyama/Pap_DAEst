function g=gammahosei_c(g_norm)

% g_norm = g_norm.^(1/2.2);

g1 = g_norm;
g2 = g_norm;
gf1 = g1>0.0031308;
gf2 = g2<=0.0031308;
g1(g1>0.0031308) = 1.055*g1(g1>0.0031308).^(1/2.4)-0.055;
g1 = gf1.*g1;
g2(g2<=0.0031308) = 12.92*g2(g2<=0.0031308);
g2 = gf2.*g2;

g = g1 + g2;

end