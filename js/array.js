function Zero2DArray(Nx, Ny, a)
{
	for(var i=0; i<Nx; i++)
	for(var j=0; j<Ny; j++)
	{
		a[i][j] = 0.;
	}
}

function Copy2DArray(Nx, Ny, src, dest)
{
	for(var i=0; i<Nx; i++)
	for(var j=0; j<Ny; j++)
		dest[i][j] = src[i][j];
}

function Init2DArray(N, val)
{
	var a = new Array(N);
	for(var i=0; i<N; i++) a[i] = new Float64Array(N);
	
	for(var i=0; i<N; i++)
	for(var j=0; j<N; j++)
		a[i][j] = val;

	return a;
}
