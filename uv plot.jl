using FITSIO
include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")
include("vfunc.jl")

PyPlot.show()
oifitsfile= "2004-data1.oifits"
data=read_oifits(oifitsfile);
uvplot(data.uv[1,:],data.uv[2,:])
u=data.uv[1,:];
v=data.uv[2,:];
fig = figure(figsize=(10,10));
scatter(u,v)
indx=data.indx_v2
scatter(u[indx],v[indx])

r=sqrt(u[indx].^2+v[indx].^2);
fig = figure(figsize=(10,10))
xlabel("Projected Baseline")
ylabel("Visibility")
title("Graph 2")
grid("on")
scatter(r, data.v2_data)

using NLopt;
chisq=(param,g)->sum(((abs2(2*besselj1(pi*param[1]*r)./(pi*param[1]*r))-data.v2_data[indx])./data.v2_data_err).^2)/length(data.v2_data[indx])
opt = Opt(:LN_NELDERMEAD, 1);
xtol_rel!(opt,1e-4)
min_objective!(opt, chisq)
(minf,minx,ret) = optimize(opt, [3e-8])
println("got $minf at $minx (returned $ret)")

theta = minx[1]
v_model = 2(besselj1(pi*theta*r))./(pi*theta*r);
v2_model = v_model.^2;
scatter(r, v2_model)

readline()
