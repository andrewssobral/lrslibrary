function data = gen_syn_data( dataname, Iobs, Inoise, Imag, Irep )

if isstruct( dataname )
    data = dataname;
else
    load( ['..\..\data\', dataname] );
end

dtemp.X = data.X;
dtemp.noise = data.mag(Imag);
dtemp.linInd = data.obs{Iobs}(:,Irep);
dtemp.noiseI = data.noise{Inoise}(Irep).ind;
dtemp.T = data.X;
dtemp.T(dtemp.noiseI) = dtemp.T(dtemp.noiseI) + data.noise{Inoise}(Irep).val(Imag);
dtemp.b = dtemp.T( dtemp.linInd );

data = dtemp;
end