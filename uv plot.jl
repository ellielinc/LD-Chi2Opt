using FITSIO
include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")

PyPlot.show()
oifitsfile= "2004-data1.oifits"
data=read_oifits(oifitsfile);
data;
data.nuv;
  data.uv[Any,:]
  uvplot(datauv[1,:],datauv[2,:])
u=data.uv[1,:];
v=data.uv[2,:];
fig = figure(figsize=(10,10))
xlabel("U (M\lambda)")
ylabel("V (M\lambda)")
title("Graph 1")
scatter(u,v)
data.v2_data;
data.nuv
data.indx_v2
u;
v;
data.indx_v2
indx=data.indx_v2
u[indx];
v[indx];
scatter(u[indx],v[indx])


data.v2_data;
data.v2_data_err;
u=u[indx];
v=v[indx];
data.v2_data;
u;
v;
u.^2;
r=sqrt(u.^2+v.^2);
fig = figure(figsize=(10,10))
xlabel("Projected Baseline")
ylabel("Visibility")
title("Graph 2")
grid("on")
scatter(r, data.v2_data)
theta= 3e-8;
v_model = 2*(besselj1(pi*theta*r))./(pi*theta*r);
v2_model = v_model.^2;
scatter(r, v2_model)

x_data = data.uv[1,:];
y_model = v_model;
y_data = data.v2_data;
sigma_data = data.v2_data_err;

f= alpha->sum((y_data-(alpha[1]*x_data+alpha[2])./sigma_data).^2)
 using NLopt;

 chisq=(alpha,g)->sum((((alpha[1]*x_data+alpha[2]) - y_data)./sigma_data).^2)

 opt = Opt(:LN_NELDERMEAD, 2);
 min_objective!(opt, chisq);
 (minf,minx,ret) = optimize(opt, [1.234, 5.678])
 println("got $minf at $minx (returned $ret)")

 aplha=[,]
 y_model= (alpha[1]*x_data+alpha[2])
