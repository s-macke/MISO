"use strict";

function MISO()
{
	this.t = 0.;
	this.dt = 0.1;
	this.iter = 0;
	
	this.mat = 
	{
		Ms: 1., 
		K1: 0.00,
		A: 13e-12,
		alpha: 0.02,
	};

	this.geo = 
	{
		w: 200., // width
		h: 200., // height
		t: 20., // thickness
  
		m: 6, // power of 2 of N
		N: 1<<6 // number of elements  (1<<m)
	};

	this.Hext =
	{
		x: 0., y: 0., z: 0.
	};
	
	var N = this.geo.N;
	
	// normalized magnetization arrays
	this.mx = Init2DArray(N, 0.);
	this.my = Init2DArray(N, 1.);
	this.mz = Init2DArray(N, 0.);

	this.dmxdt = Init2DArray(N, 0.);
	this.dmydt = Init2DArray(N, 0.);
	this.dmzdt = Init2DArray(N, 0.);

	this.dmxdtold = Init2DArray(N, 0.);
	this.dmydtold = Init2DArray(N, 0.);
	this.dmzdtold = Init2DArray(N, 0.);
	
	// normalized magnetization arrays
	this.mxnew = Init2DArray(N, 0.);
	this.mynew = Init2DArray(N, 1.);
	this.mznew = Init2DArray(N, 0.);
	
	// effective field
	this.Heffx = Init2DArray(N, 0.);
	this.Heffy = Init2DArray(N, 0.);
	this.Heffz = Init2DArray(N, 0.);
	
	this.Init();
}

MISO.prototype.Init = function()
{
	this.demag = new DEMAG(this.geo);
}

MISO.prototype.Calcdmdt = function(mx, my, mz, dmxdt, dmydt, dmzdt)
{
	var i=0, j=0;

	var m0 = 4. * Math.PI * 1e-7;
	var N = this.geo.N;
	this.demag.Demag(mx, my, mz);
// ------------------------------------
	
	var Heffx = this.Heffx;
	var Heffy = this.Heffy;
	var Heffz = this.Heffz;
	var Hx = this.demag.Hx;
	var Hy = this.demag.Hy;
	var Hz = this.demag.Hz;
	
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		// External field
		Heffx[i][j] = this.Hext.x / m0 + this.mat.Ms*Hx[i][j]/m0;
		Heffy[i][j] = this.Hext.y / m0 + this.mat.Ms*Hy[i][j]/m0;
		Heffz[i][j] = this.Hext.z / m0 + this.mat.Ms*Hz[i][j]/m0;
		// Magnetocrystalline anisotropy
		var scalar = my[i][j] * 0. + my[i][j] * 1. + mz[i][j] * 0.;
		//Heffx[i][j] +=  0. * (1./mat.Ms * scalar * 2.*mat.K1);
		Heffy[i][j] +=  1. * (1./this.mat.Ms * scalar * 2.*this.mat.K1);
		//Heffz[i][j] +=  0. * (1./mat.Ms * scalar * 2.*mat.K1);
	}


	var pre = 1.e18 * 2. * this.mat.A / this.mat.Ms;  // unit A/m
	var dx = this.geo.w / this.geo.N;
	var dy = this.geo.h / this.geo.N;
	for(i=1; i<N-1; i++)
	for(j=1; j<N-1; j++)
	{
		Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
		Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
		Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );
	}

//boundary conditions, von Neuman I think

	for(j=1; j<N-1; j++)
	{
		i = 0;
		Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
		Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
		Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );
		i = N-1;
		Heffx[i][j] += pre*( (mx[i][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
		Heffy[i][j] += pre*( (my[i][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
		Heffz[i][j] += pre*( (mz[i][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );
	}

	for(i=1; i<N-1; i++)
	{
		j = 0;
		Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j])/(dy*dy) );
		Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j])/(dy*dy) );
		Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j])/(dy*dy) );
		j = N-1;
		Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
		Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
		Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );
	}

	// and finally the corners
	
	i=0; j=0;
	Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j])/(dy*dy) );
	Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j])/(dy*dy) );
	Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j])/(dy*dy) );

	i=N-1; j=0;
	Heffx[i][j] += pre*( (mx[i][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j+1] - 2.*mx[i][j] + mx[i][j])/(dy*dy) );
	Heffy[i][j] += pre*( (my[i][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j+1] - 2.*my[i][j] + my[i][j])/(dy*dy) );
	Heffz[i][j] += pre*( (mz[i][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j+1] - 2.*mz[i][j] + mz[i][j])/(dy*dy) );
	
	i=N-1; j=N-1;
	Heffx[i][j] += pre*( (mx[i][j] - 2.*mx[i][j] + mx[i-1][j])/(dx*dx) + (mx[i][j] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
	Heffy[i][j] += pre*( (my[i][j] - 2.*my[i][j] + my[i-1][j])/(dx*dx) + (my[i][j] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
	Heffz[i][j] += pre*( (mz[i][j] - 2.*mz[i][j] + mz[i-1][j])/(dx*dx) + (mz[i][j] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );

	i=0; j=N-1;
	Heffx[i][j] += pre*( (mx[i+1][j] - 2.*mx[i][j] + mx[i][j])/(dx*dx) + (mx[i][j] - 2.*mx[i][j] + mx[i][j-1])/(dy*dy) );
	Heffy[i][j] += pre*( (my[i+1][j] - 2.*my[i][j] + my[i][j])/(dx*dx) + (my[i][j] - 2.*my[i][j] + my[i][j-1])/(dy*dy) );
	Heffz[i][j] += pre*( (mz[i+1][j] - 2.*mz[i][j] + mz[i][j])/(dx*dx) + (mz[i][j] - 2.*mz[i][j] + mz[i][j-1])/(dy*dy) );

	var pre = -2.21e5/(1. + this.mat.alpha*this.mat.alpha)*1e-12;
	var pre2 = this.mat.alpha*2.21e5/(1. + this.mat.alpha*this.mat.alpha)*1e-12;	
	var jhx, jhy, jhz;
	var jax, jay, jaz;
	var l;
	
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		jhx = my[i][j]*Heffz[i][j] - mz[i][j]*Heffy[i][j];
		jhy = mz[i][j]*Heffx[i][j] - mx[i][j]*Heffz[i][j];
		jhz = mx[i][j]*Heffy[i][j] - my[i][j]*Heffx[i][j];
		
		jax = my[i][j]*jhz - mz[i][j]*jhy;
		jay = mz[i][j]*jhx - mx[i][j]*jhz;
		jaz = mx[i][j]*jhy - my[i][j]*jhx;		
		
		dmxdt[i][j] = pre * jhx - pre2*jax;
		dmydt[i][j] = pre * jhy - pre2*jay;
		dmzdt[i][j] = pre * jhz - pre2*jaz;
	}
}

MISO.prototype.Step = function()
{
	var i, j;
	var l = 0.0;
	
	var temp = null;
	var N = this.geo.N;
	/*
	temp = dmxdt;
	dmxdt = dmxdtold;
	dmxdtold = temp;
	
	temp = dmydt;
	dmydt = dmydtold;
	dmydtold = temp;
	
	temp = dmzdt;
	dmzdt = dmzdtold;
	dmzdtold = temp;
	*/
	
	var error = 1e-99;
	/*
	for(i=0; i<geo.N; i++)
	for(j=0; j<geo.N; j++)
	{
		var mxdiff = dmxdt[i][j] * dt - (1.5*dmxdt[i][j] * dt - 0.5*dmxdtold[i][j] * dt);
		var mydiff = dmydt[i][j] * dt - (1.5*dmydt[i][j] * dt - 0.5*dmydtold[i][j] * dt);
		var mzdiff = dmzdt[i][j] * dt - (1.5*dmzdt[i][j] * dt - 0.5*dmzdtold[i][j] * dt);
		
		var e = mxdiff*mxdiff + mydiff*mydiff + mzdiff*mzdiff;
		if (e > error) error = e;
	}
	dt = 1./Math.sqrt(error)*1e-8;
	console.log(dt);
	if (dt > 0.3) dt = 0.3;
*/
	var mx = this.mx;
	var my = this.my;
	var mz = this.mz;
	var mxnew = this.mxnew;
	var mynew = this.mynew;
	var mznew = this.mznew;
	var dmxdt = this.dmxdt;
	var dmydt = this.dmydt;
	var dmzdt = this.dmzdt;
	var dmxdtold = this.dmxdtold;
	var dmydtold = this.dmydtold;
	var dmzdtold = this.dmzdtold;
	var dt = this.dt;

	// Predictor	
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		mxnew[i][j] = mx[i][j] + 1.5*dmxdt[i][j] * dt - 0.5*dmxdtold[i][j] * dt;
		mynew[i][j] = my[i][j] + 1.5*dmydt[i][j] * dt - 0.5*dmydtold[i][j] * dt;
		mznew[i][j] = mz[i][j] + 1.5*dmzdt[i][j] * dt - 0.5*dmzdtold[i][j] * dt;
	}
	
	temp = dmxdt;
	dmxdt = dmxdtold;
	dmxdtold = temp;
	
	temp = dmydt;
	dmydt = dmydtold;
	dmydtold = temp;
	
	temp = dmzdt;
	dmzdt = dmzdtold;
	dmzdtold = temp;
	miso.Calcdmdt(mxnew, mynew, mznew, dmxdt, dmydt, dmzdt);
	
	// Corrector
	
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		mxnew[i][j] = mx[i][j] + 0.5*dmxdt[i][j] * dt + 0.5*dmxdtold[i][j] * dt;
		mynew[i][j] = my[i][j] + 0.5*dmydt[i][j] * dt + 0.5*dmydtold[i][j] * dt;
		mznew[i][j] = mz[i][j] + 0.5*dmzdt[i][j] * dt + 0.5*dmzdtold[i][j] * dt;
	}
	miso.Calcdmdt(mxnew, mynew, mznew, dmxdt, dmydt, dmzdt);
	
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		mx[i][j] += 0.5*dmxdt[i][j] * dt + 0.5*dmxdtold[i][j] * dt;
		my[i][j] += 0.5*dmydt[i][j] * dt + 0.5*dmydtold[i][j] * dt;
		mz[i][j] += 0.5*dmzdt[i][j] * dt + 0.5*dmzdtold[i][j] * dt;
		//mx[i][j] = mxnew[i][j];
		//my[i][j] = mynew[i][j];
		//mz[i][j] = mznew[i][j];
		
		l = Math.sqrt(mx[i][j]*mx[i][j] + my[i][j]*my[i][j] + mz[i][j]*mz[i][j]);
		mx[i][j] /= l;
		my[i][j] /= l;
		mz[i][j] /= l;
	}
	
	miso.t += miso.dt;
	miso.iter++;
}

