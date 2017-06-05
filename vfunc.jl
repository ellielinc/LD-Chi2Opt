#power law
function visibility_ldhest(param, v_r)
maxk=10;
f=k->gamma(param[2]./2+2)./gamma(param[2]/2+k+2)*gamma(k+1))*(-(pi*v_r*param[1]).^2./4)^k;
V = map(f,collect(1:maxk));
#sum(V)
end

#quadratic law
function visibility_ldquad(param, v_r)
zeta = pi*v_r*param[1];
V = ((1-param[2]-param[3])*(besselj1(zeta))./zeta)
+((param[1]+2param[3])./sqrt(2/pi))*((sqrt(2/(pi*zeta))*((sin(zeta)./zeta)-cos(zeta))./zeta.^(3/2))-2param[3]*(besselj(2, zeta)./zeta.^2)./((1/2)-(param[2]./6)-(param[3]./12))
end

#uniform disk
function visibility_ud(param, v_r)
V = 2*(besselj1(pi*param*v_r))./(pi*param*v_r)
end
