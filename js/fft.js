//   Perform a 2D FFT inplace given a complex 2D array
//   The direction dir, 1 for forward, -1 for reverse
//   The size of the array (nx,ny)
//   Return false if there are memory problems or
//      the dimensions are not powers of 2

"use strict";

function FFT2D(m)
{
	this.m = m;
	this.n = 1;
	for (var i=0;i<m;i++) this.n *= 2;
	this.nx = this.n;
	this.ny = this.n;
	this.real = new Float64Array(this.nx);
	this.imag = new Float64Array(this.nx);
}

FFT2D.prototype.Forward2 = function(cx, cy, dx, dy, dir)
{
	var i, j;
	var real = this.real;
	var imag = this.imag;

	// Transform the rows
	for (j=0;j<this.ny;j++)
	{
		for (i=0;i<this.nx;i++) 
		{
			real[i] = cx[i][j];
			imag[i] = cy[i][j];
		}
		this.FFT(dir, real, imag);
		for (i=0;i<this.nx;i++)
		{
			dx[i][j] = real[i];
			dy[i][j] = imag[i];
		}
	}

	// Transform the columns
	for (i=0;i<this.nx;i++)
	{
		this.FFT(dir, dx[i], dy[i]);
	}
}

FFT2D.prototype.Forward = function(cx, cy, dir)
{
	var i, j;
	var real = this.real;
	var imag = this.imag;

	// Transform the rows
	for (j=0;j<this.ny;j++)
	{
		for (i=0;i<this.nx;i++) 
		{
			real[i] = cx[i][j];
			imag[i] = cy[i][j];
		}
		this.FFT(dir, real, imag);
		for (i=0;i<this.nx;i++)
		{
			cx[i][j] = real[i];
			cy[i][j] = imag[i];
		}
	}

	// Transform the columns
	for (i=0;i<this.nx;i++)
	{
		for (j=0;j<this.ny;j++)
		{
			real[j] = cx[i][j];
			imag[j] = cy[i][j];
		}
		this.FFT(dir, real, imag);
		for (j=0;j<this.ny;j++)
		{
			cx[i][j] = real[j];
			cy[i][j] = imag[j];
		}
	}
}
  
/*
   This computes an in-place complex-to-complex FFT 
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform 
*/
FFT2D.prototype.FFT = function(dir, x, y)
{
	var i,i1,j,k,i2,l,l1,l2; // int
	var c1,c2,tx,ty,t1,t2,u1,u2,z; // double

	// Calculate the number of points 

	// Do the bit reversal
	i2 = this.n >> 1;
	j = 0;
	for (i=0;i<this.n-1;i++)
	{
		if (i < j)
		{
			tx = x[i];
			ty = y[i];
			x[i] = x[j];
			y[i] = y[j];
			x[j] = tx;
			y[j] = ty;
		}
		k = i2;
		while (k <= j)
		{
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	// Compute the FFT
	c1 = -1.0;
	c2 = 0.0;
	l2 = 1;
	for (l=0;l<this.m;l++) 
	{
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0; 
		u2 = 0.0;
		for (j=0;j<l1;j++)
		{
			for (i=j;i<this.n;i+=l2)
			{
				i1 = i + l1;
				t1 = u1 * x[i1] - u2 * y[i1];
				t2 = u1 * y[i1] + u2 * x[i1];
				x[i1] = x[i] - t1; 
				y[i1] = y[i] - t2;
				x[i] += t1;
				y[i] += t2;
			}
			z =  u1 * c1 - u2 * c2;
			u2 = u1 * c2 + u2 * c1;
			u1 = z;
		}
		c2 = Math.sqrt((1.0 - c1) / 2.0);
		if (dir == 1) c2 = -c2;
		c1 = Math.sqrt((1.0 + c1) / 2.0);
	}

	// Scaling for forward transform
	if (dir == 1)
		for (i=0;i<this.n;i++)
		{
			x[i] /= this.n;
			y[i] /= this.n;
		}
}
