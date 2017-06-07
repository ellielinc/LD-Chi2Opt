function optimize_model_v2(data, visfunc, init_param)
  nparams=length(init_param)
  indx= data.indx_v2
  nv2 = length(data.v2_data[indx]);
  r=data.baseline_v2
  chisq=(param,g)->sum(((abs2(visfunc(param,r)-data.v2_data[indx])./data.v2_data_err[indx]).^2))/nv2;
  opt = Opt(:LN_NELDERMEAD, nparams);
  min_objective!(opt, chisq)
  (minf,minx,ret) = optimize(opt, init_param);
  println("got $minf at $minx (returned $ret)")
  #MAS = minx[1]/(4.848137*10-9);
  #println("$MAS MAS")
  #println("minx = $minx")
  cvis_model = visfunc(minx,data.baseline_v2)
  return (minf,minx,ret,cvis_model)
end
