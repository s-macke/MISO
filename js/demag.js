
// strayfield calculation code copied from the software oommf
// Routines to calculate kernel coefficients
// See Newell et al. for details. The code below follows the
// naming conventions in that paper.

"use strict";

function DEMAG(geo)
{
	var i = 0;
	var j = 0;
	this.geo = geo;
	
	var l = 1.;
	var h = 1.;
	var e = geo.t / (geo.w/geo.N);
	//console.log(e);
	var z = 0.; // Scale size so distance between adjacent
				// cells is 1.

	var N = geo.N;
				
	this.fft = new FFT2D(geo.m+1);
	
	this.A00r = Init2DArray(geo.N*2, 0.);
	this.A00i = Init2DArray(geo.N*2, 0.);
	this.A01r = Init2DArray(geo.N*2, 0.);
	this.A01i = Init2DArray(geo.N*2, 0.);
	this.A11r = Init2DArray(geo.N*2, 0.);
	this.A11i = Init2DArray(geo.N*2, 0.);
	this.A22r = Init2DArray(geo.N*2, 0.);
	this.A22i = Init2DArray(geo.N*2, 0.);

	this.Hx = Init2DArray(geo.N, 0.);
	this.Hy = Init2DArray(geo.N, 0.);
	this.Hz = Init2DArray(geo.N, 0.);
	
	this.mxfr = Init2DArray(geo.N*2, 0.);
	this.myfr = Init2DArray(geo.N*2, 0.);
	this.mzfr = Init2DArray(geo.N*2, 0.);
	this.mxfi = Init2DArray(geo.N*2, 0.);
	this.myfi = Init2DArray(geo.N*2, 0.);
	this.mzfi = Init2DArray(geo.N*2, 0.);


	
	/*
	l = geo.w / geo.N;
	h = geo.h / geo.N;
	e = geo.t;
	*/
	
	var scale= -1. / (4. * Math.PI * l * h * e);
	
	//scaling for the implemented fft algorithm
	scale *= N * N * 2 * 2;
	
	//scale= -1. / (4. * Math.PI);
	
	// According (16) in Newell's paper, the demag field is given by
	//                        H = -N*M
	// where N is the "demagnetizing tensor," with components Nxx, Nxy,
	// etc.  With the '-1' in 'scale' we store '-N' instead of 'N',
	// so we don't have to multiply the output from the FFT + iFFT
	// by -1 in ConstMagField() below.
  
	var Nx = geo.N; // 64
	var Ny = geo.N; // 64
/*
	var NFy = 2.*geo.N;	// 128
	
	// For real inputs, the resulting complex transform is symmetric, so only the top half is computed.
	var NFx = (2.*geo.N)/2+1;	// 65  

	var Kx = 2*(NFx - 1);   // 128 Space domain kernel dimensions 
	var Ky = NFy;			// 128
*/
	var NFy = 2.*geo.N;	// 128
	var NFx = 2.*geo.N;	// 128
	var Kx = NFx;
	var Ky = NFy;

	var spacekernelr, spacekerneli;
	spacekernelr = Init2DArray(Kx, 0.);
	spacekerneli = Init2DArray(Ky, 0.);

	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		spacekernelr[i][j] = scale * this.GetSDA00(i, j, z, l, h, e);
		if (i>0) spacekernelr[Kx-i][j] = spacekernelr[i][j];
		if (j>0) spacekernelr[i][Ky-j] = spacekernelr[i][j];
		if (i>0 && j>0) spacekernelr[Kx-i][Ky-j] = spacekernelr[i][j];
	}	
	Zero2DArray(Kx, Ky, spacekerneli);
	this.fft.Forward(spacekernelr, spacekerneli, 1);
	Copy2DArray(Kx, Ky, spacekernelr, this.A00r);
	Copy2DArray(Kx, Ky, spacekerneli, this.A00i);
/*
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		//this.A00r[i][j] = Math.log(Math.abs(this.A00r[i][j]));
		this.A00r[i][j] = Math.log(Math.abs(this.A00r[i][j]));
	}
*/
	/*
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		document.write("" + i + " " + j + " " + this.A00r[i][j] + "<br>");
	}
	abort();
*/
	
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		spacekernelr[i][j] = scale * this.GetSDA11(i, j, z, l, h, e);
		if (i>0) spacekernelr[Kx-i][j] = spacekernelr[i][j];
		if (j>0) spacekernelr[i][Ky-j] = spacekernelr[i][j];
		if (i>0 && j>0) spacekernelr[Kx-i][Ky-j] = spacekernelr[i][j];
	}
	Zero2DArray(Kx, Ky, spacekerneli);
	this.fft.Forward(spacekernelr, spacekerneli, 1);
	Copy2DArray(Kx, Ky, spacekernelr, this.A11r);
	Copy2DArray(Kx, Ky, spacekerneli, this.A11i);

	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		spacekernelr[i][j] = scale*this.GetSDA22(i, j, z, l, h, e);
		if (i>0) spacekernelr[Kx-i][j] = spacekernelr[i][j];
		if (j>0) spacekernelr[i][Ky-j] = spacekernelr[i][j];
		if (i>0 && j>0) spacekernelr[Kx-i][Ky-j] = spacekernelr[i][j];
	}
	Zero2DArray(Kx, Ky, spacekerneli);
	this.fft.Forward(spacekernelr, spacekerneli, 1);	
	Copy2DArray(Kx, Ky, spacekernelr, this.A22r);
	Copy2DArray(Kx, Ky, spacekerneli, this.A22i);

	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		spacekernelr[i][j] = scale*this.GetSDA01(i,j,z,l,h,e);
		if (i>0) spacekernelr[Kx-i][j] = -1. * spacekernelr[i][j];
		if (j>0) spacekernelr[i][Ky-j] = -1. * spacekernelr[i][j];
		if (i>0 && j>0) spacekernelr[Kx-i][Ky-j] = spacekernelr[i][j];
    // A01 is odd in x and y.
	}
	Zero2DArray(Kx, Ky, spacekerneli);
	this.fft.Forward(spacekernelr, spacekerneli, 1);	
	Copy2DArray(Kx, Ky, spacekernelr, this.A01r);
	Copy2DArray(Kx, Ky, spacekerneli, this.A01i);
	/*
	for(i=0; i<Nx; i++)
	for(j=0; j<Ny; j++)
	{
		this.A22r[i][j] = Math.log(Math.abs(this.A22r[i][j]));
	}
	*/	
}


DEMAG.prototype.Demag = function(mx, my, mz)
{
	var i=0, j=0;

	var fft = this.fft;	
	var mxfr = this.mxfr;
	var mxfi = this.mxfi;
	var myfr = this.myfr;
	var myfi = this.myfi;
	var mzfr = this.mzfr;
	var mzfi = this.mzfi;
	var N = this.geo.N;

// ------------------------------------
// calculate demag field
/*
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		mx[i][j] = 0.;
		my[i][j] = 1.;
		mz[i][j] = 0.;
	}
*/
	Zero2DArray(N*2, N*2, mxfr);
	Zero2DArray(N*2, N*2, mxfi);
	Copy2DArray(N, N, mx, mxfr);
	fft.Forward(mxfr, mxfi, 1);
	
	Zero2DArray(N*2, N*2, myfr);
	Zero2DArray(N*2, N*2, myfi);
	Copy2DArray(N, N, my, myfr);
	fft.Forward(myfr, myfi, 1);	
	
	Zero2DArray(N*2, N*2, mzfr);
	Zero2DArray(N*2, N*2, mzfi);
	Copy2DArray(N, N, mz, mzfr);
	fft.Forward(mzfr, mzfi, 1);
	
	var A00r = this.A00r;
	var A00i = this.A00i;
	var A01r = this.A01r;
	var A01i = this.A01i;
	var A11r = this.A11r;
	var A11i = this.A11i;
	var A22r = this.A22r;
	var A22i = this.A22i;
	
	var tempxr = 0., tempxi = 0.;
	var tempyr = 0., tempyi = 0.;
	var tempzr = 0., tempzi = 0.;

	for(i=0; i<N*2; i++)
	for(j=0; j<N*2; j++)
	{
		tempxr = A00r[i][j] * mxfr[i][j] + A01r[i][j] * myfr[i][j] - A00i[i][j] * mxfi[i][j] - A01i[i][j] * myfi[i][j];
		tempxi = A00r[i][j] * mxfi[i][j] + A01r[i][j] * myfi[i][j] + A00i[i][j] * mxfr[i][j] + A01i[i][j] * myfr[i][j];
		tempyr = A01r[i][j] * mxfr[i][j] + A11r[i][j] * myfr[i][j] - A01i[i][j] * mxfi[i][j] - A11i[i][j] * myfi[i][j];
		tempyi = A01r[i][j] * mxfi[i][j] + A11r[i][j] * myfi[i][j] + A01i[i][j] * mxfr[i][j] + A11i[i][j] * myfr[i][j]		
		
		tempzr = A22r[i][j] * mzfr[i][j] - A22i[i][j] * mzfi[i][j];
		tempzi = A22r[i][j] * mzfi[i][j] + A22i[i][j] * mzfr[i][j];	
		mxfr[i][j] = tempxr;
		mxfi[i][j] = tempxi;
		myfr[i][j] = tempyr;
		myfi[i][j] = tempyi;
		mzfr[i][j] = tempzr;
		mzfi[i][j] = tempzi;
		
	}
	fft.Forward(mxfr, mxfi, -1);
	fft.Forward(myfr, myfi, -1);
	fft.Forward(mzfr, mzfi, -1);
	
	//Copy2DArray(N, N, A00r, this.Hz);
/*
	for(i=0; i<N; i++)
	for(j=0; j<N; j++)
	{
		document.write(i+ " " + j + " " + mzfr[i][j] + "<br>");
	}
	abort();
*/
	
	Copy2DArray(N, N, mxfr, this.Hx);
	Copy2DArray(N, N, myfr, this.Hy);
	Copy2DArray(N, N, mzfr, this.Hz);
}



DEMAG.prototype.SelfDemagDx = function(x, y, z)
{ 
// Note: Assumes x, y, and z are all >0.
  // Note 2: egcs-2.91.57 on Linux/x86 with -O1 mangles this
  //  function (produces NaN's) unless we manually group terms.
  // Here self demag Hx = -Dx Mx / (4*PI*x*y*z).

  if(x==y && y==z) return 4.*Math.PI*x*x*x/3.;  // Special case: cube

  var xsq = x*x, ysq = y*y, zsq = z*z;
  var diag = Math.sqrt(xsq + ysq + zsq);
  var arr = 0.;

  var mpxy = (x-y)*(x+y);
  var mpxz = (x-z)*(x+z);

  arr += -4.*(2.*xsq*x - ysq*y - zsq*z);
  arr +=  4.*(xsq + mpxy)*Math.sqrt(xsq + ysq);
  arr +=  4.*(xsq + mpxz)*Math.sqrt(xsq + zsq);
  arr += -4.*(ysq + zsq)*Math.sqrt(ysq + zsq);
  arr += -4.*diag*(mpxy + mpxz);

  arr += 24.*x*y*z*Math.atan(y*z/(x*diag));
  arr += 12.*(z + y)*xsq*Math.log(x);

  arr += 12.*z*ysq*Math.log((Math.sqrt(ysq+zsq)+z)/y);
  arr += -12.*z*xsq*Math.log(Math.sqrt(xsq+zsq)+z);
  arr += 12.*z*mpxy*Math.log(diag+z);
  arr += -6.*z*mpxy*Math.log(xsq+ysq);

  arr +=  12.*y*zsq*Math.log((Math.sqrt(ysq+zsq)+y)/z);
  arr += -12.*y*xsq*Math.log(Math.sqrt(xsq+ysq)+y);
  arr +=  12.*y*mpxz*Math.log(diag+y);
  arr +=  -6.*y*mpxz*Math.log(xsq+zsq);

  var Dx = arr/3.;
  return Dx;
}

DEMAG.prototype.SelfDemagNx = function(xsize, ysize, zsize)
{
  if (xsize<=0.0 || ysize<=0.0 || zsize<=0.0) return 0.;
  return SelfDemagDx(xsize, ysize, zsize) / (4.*Math.PI*xsize*ysize*zsize);
  /// 4*PI*xsize*ysize*zsize is scaling factor to convert internal
  /// Dx value to proper scaling consistent with Newell's Nxx,
  /// where Hx = -Nxx.Mx (formula (16) in Newell).
}

DEMAG.prototype.SelfDemagNy = function(xsize, ysize, zsize)
{ 
	return SelfDemagNx(ysize, zsize, xsize); 
}

DEMAG.prototype.SelfDemagNz = function(xsize, ysize, zsize)
{
	return SelfDemagNx(zsize, xsize, ysize);
}

DEMAG.prototype.f = function(x, y, z)
{ // There is mucking around here to handle case where imports
  // are near zero.  In particular, asinh(t) is written as
  // log(t+sqrt(1+t)) because the latter appears easier to
  // handle if t=y/x (for example) as x -> 0.

 // This function is even; the fabs()'s just simplify special case handling.
  x = Math.abs(x); var xsq = x*x;
  y = Math.abs(y); var ysq = y*y;
  z = Math.abs(z); var zsq = z*z;

  var R = xsq + ysq + zsq;
  if (R <= 0.0) return 0.0; else R = Math.sqrt(R);

  // f(x,y,z)
  var piece = 0;
  if (z > 0.)
  { // For 2D grids, half the calls from F1 have z==0.
    var temp1, temp2, temp3;
    piece += 2. * (2.*xsq - ysq - zsq) * R;
	temp1 = x*y*z;
    if (temp1 > 0.) piece += -12.*temp1*Math.atan2(y*z, x*R);
	temp2 = xsq + zsq;
    if (y>0. && temp2>0.)
	{
      var dummy = Math.log(((y + R)*(y + R))/temp2);
      piece +=  3.*y*zsq*dummy;
      piece += -3.*y*xsq*dummy;
    }
	temp3 = xsq+ysq;
    if (temp3 > 0.) 
	{
      var dummy = Math.log(((z + R)*(z + R))/temp3);
      piece +=  3.*z*ysq*dummy;
      piece += -3.*z*xsq*dummy;
    }
  }
  else 
  {
    // z==0
    if (x==y) 
	{
      var K = -2.45981439737106805379;
      /// K = 2*sqrt(2)-6*log(1+sqrt(2))
      piece += K*xsq*x;
    } 
	else 
	{
      piece += 2.*(2.*xsq - ysq)*R;
      if (y>0. && x>0.) piece += -6*y*xsq*Math.log((y + R)/x);
    }
  }
  return piece / 12.;
}


DEMAG.prototype.GetSDA00 = function(x, y, z, dx, dy, dz)
{ // This is Nxx*(4*PI*tau) in Newell's paper
  var result = 0.;
  if(x==0. && y==0. && z==0.) 
  {
    // Self demag term.  The base routine can handle x==y==z==0,
    // but this should be more accurate.
    result = this.SelfDemagDx(dx, dy, dz);
  } 
  else 
  {
    // Simplified (collapsed) formula based on Newell's paper.
    // This saves about half the calls to f().  There is still
    // quite a bit of redundancy from one cell site to the next,
    // but as this is an initialization-only issue speed shouldn't
    // be too critical.
    var arr = 0;
    arr += -1.*this.f(x+dx,y+dy,z+dz);
    arr += -1.*this.f(x+dx,y-dy,z+dz);
    arr += -1.*this.f(x+dx,y-dy,z-dz);
    arr += -1.*this.f(x+dx,y+dy,z-dz);
    arr += -1.*this.f(x-dx,y+dy,z-dz);
    arr += -1.*this.f(x-dx,y+dy,z+dz);
    arr += -1.*this.f(x-dx,y-dy,z+dz);
    arr += -1.*this.f(x-dx,y-dy,z-dz);

    arr +=  2.*this.f(x,y-dy,z-dz);
    arr +=  2.*this.f(x,y-dy,z+dz);
    arr +=  2.*this.f(x,y+dy,z+dz);
    arr +=  2.*this.f(x,y+dy,z-dz);
    arr +=  2.*this.f(x+dx,y+dy,z);
    arr +=  2.*this.f(x+dx,y,z+dz);
    arr +=  2.*this.f(x+dx,y,z-dz);
    arr +=  2.*this.f(x+dx,y-dy,z);
    arr +=  2.*this.f(x-dx,y-dy,z);
    arr +=  2.*this.f(x-dx,y,z+dz);
    arr +=  2.*this.f(x-dx,y,z-dz);
    arr +=  2.*this.f(x-dx,y+dy,z);

    arr += -4.*this.f(x,y-dy,z);
    arr += -4.*this.f(x,y+dy,z);
    arr += -4.*this.f(x,y,z-dz);
    arr += -4.*this.f(x,y,z+dz);
    arr += -4.*this.f(x+dx,y,z);
    arr += -4.*this.f(x-dx,y,z);

    arr +=  8.*this.f(x,y,z);
    result=arr;
  }
  return result;
  /// Multiply result by 1./(4*PI*dx*dy*dz) to get effective "demag"
  /// factor Nxx from Newell's paper.
}

DEMAG.prototype.GetSDA11 = function(x, y, z, l, h, e)
{
	return this.GetSDA00(y, z, x, h, e, l);
}

DEMAG.prototype.GetSDA22 = function(x, y, z, l, h, e)
{
	return this.GetSDA00(z, x, y, e, l, h); 
}

DEMAG.prototype.g = function(x, y, z)
{ // There is mucking around here to handle case where imports
  // are near zero.  In particular, asinh(t) is written as
  // log(t+sqrt(1+t)) because the latter appears easier to
  // handle if t=y/x (for example) as x -> 0.

  var result_sign=1.0;
  if(x<0.0) result_sign *= -1.0;  
  if(y<0.0) result_sign *= -1.0;
  x=Math.abs(x); y=Math.abs(y); z=Math.abs(z);  // This function is even in z and
  /// odd in x and y.  The fabs()'s simplify special case handling.

  var xsq = x*x,ysq = y*y,zsq = z*z;
  var R = xsq+ysq+zsq;
  if (R <= 0.) return 0.0;
  else R = Math.sqrt(R);

  // g(x,y,z)
  var piece = 0.;
  piece += -2.*x*y*R;
  if(z>0.) 
  { // For 2D grids, 1/3 of the calls from GetSDA01 have z==0.
    piece += -z*zsq*Math.atan2(x*y,z*R);
    piece += -3.*z*ysq*Math.atan2(x*z,y*R);
    piece += -3.*z*xsq*Math.atan2(y*z,x*R);

    var temp1, temp2, temp3;
    if ((temp1 = xsq+ysq)>0.)
      piece += 6.*x*y*z*Math.log((z+R)/Math.sqrt(temp1));

    if ((temp2 = ysq+zsq)>0.)
      piece += y*(3.*zsq-ysq)*Math.log((x+R)/Math.sqrt(temp2));

    if ((temp3 = xsq+zsq)>0.)
      piece += x*(3.*zsq-xsq)*Math.log((y+R)/Math.sqrt(temp3));
  } 
  else 
  {
    // z==0.
    if(y>0.) piece += -y*ysq*Math.log((x+R)/y);
    if(x>0.) piece += -x*xsq*Math.log((y+R)/x);
  }

  return (1./6.) * result_sign * piece;
}


DEMAG.prototype.GetSDA01 = function(x, y, z, l, h, e)
{ // This is Nxy*(4*PI*tau) in Newell's paper.

  // Simplified (collapsed) formula based on Newell's paper.
  // This saves about half the calls to g().  There is still
  // quite a bit of redundancy from one cell site to the next,
  // but as this is an initialization-only issure speed shouldn't
  // be too critical.
  var arr = 0.;

  arr += -1.*this.g(x-l, y-h, z-e);
  arr += -1.*this.g(x-l, y-h, z+e);
  arr += -1.*this.g(x+l, y-h, z+e);
  arr += -1.*this.g(x+l, y-h, z-e);
  arr += -1.*this.g(x+l, y+h, z-e);
  arr += -1.*this.g(x+l, y+h, z+e);
  arr += -1.*this.g(x-l, y+h, z+e);
  arr += -1.*this.g(x-l, y+h, z-e);

  arr +=  2.*this.g(x, y+h, z-e);
  arr +=  2.*this.g(x, y+h, z+e);
  arr +=  2.*this.g(x, y-h, z+e);
  arr +=  2.*this.g(x, y-h, z-e);
  arr +=  2.*this.g(x-l, y-h, z);
  arr +=  2.*this.g(x-l, y+h, z);
  arr +=  2.*this.g(x-l, y, z-e);
  arr +=  2.*this.g(x-l, y, z+e);
  arr +=  2.*this.g(x+l, y, z+e);
  arr +=  2.*this.g(x+l, y, z-e);
  arr +=  2.*this.g(x+l, y-h, z);
  arr +=  2.*this.g(x+l, y+h, z);

  arr += -4.*this.g(x-l, y, z);
  arr += -4.*this.g(x+l, y, z);
  arr += -4.*this.g(x, y, z+e);
  arr += -4.*this.g(x, y, z-e);
  arr += -4.*this.g(x, y-h, z);
  arr += -4.*this.g(x, y+h, z);

  arr +=  8.*this.g(x, y, z);
  
  return arr;
  // Multiply result by 1./(4*PI*l*h*e) to get effective "demag"
  // factor Nxy from Newell's paper.
}
