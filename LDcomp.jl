using FITSIO;
using NLopt;
include("readoifits.jl")
include("setupft.jl")
include("oichi2.jl")
include("oiplot.jl")
include("vfunc.jl")
include("chi2.jl")
PyPlot.show()
oifitsfile= "2004-data1.oifits"
data=read_oifits(oifitsfile);
uvplot(data.uv[1,:],data.uv[2,:])


#visibility_ud, [1.0]
#visibility_ldpow, [1.0,0.5]
#visibility_ldquad, [1.0,0.5,0.5]


opt_chi2,opt_param,opt_diag,cvis_model = optimize_model_v2(data, visibility_ldpow, [2.0,0.5])
v2_model = cvis_to_v2(cvis_model, data.indx_v2);
v2plot_modelvsdata(data.baseline_v2,data.v2_data,data.v2_data_err, v2_model);

#Power
#chisq2(3e-8);
#v2_model2 = abs2(visibility_ldpow(minx, r));
#v2plot_modelvsdata(r,data.v2_data[indx],data.v2_data_err[indx], v2_model2)

#Quad
#chisq3(3e-8)
#v2_model3 = visibility_ldquad(minx, v_r);
#v2plot_modelvsdata(r,data.v2_data[indx],data.v2_data_err[indx], v2_model3)
