"use strict";

function GL() // constructor
{
	this.canvas = document.getElementById("glcanvas");
	this.canvas.width = window.innerWidth;
	this.canvas.height = window.innerHeight;
	this.gl = this.Init(this.canvas);
    if (!this.gl)
	{
		alert("Could not initialise WebGL, sorry :-(");
		return;
    }
	this.gl.viewportWidth = this.canvas.width;
	this.gl.viewportHeight = this.canvas.height;
	
// vertex shader
	var fragmentshader = this.GetShader(this.gl, "shader-fs");
	var vertexshader = this.GetShader(this.gl, "shader-vs");

	this.webglprogramobject = this.gl.createProgram();
	this.gl.attachShader(this.webglprogramobject, vertexshader);
	this.gl.attachShader(this.webglprogramobject, fragmentshader);
	this.gl.linkProgram(this.webglprogramobject);
	this.gl.useProgram(this.webglprogramobject);
	
	// modelviewmatrix
	this.globalmvMatrix = new Float32Array
		([
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
		]);
	//Model.mvMatrix = mat4.identity();
	//mat4.translate(Model.mvMatrix, [0.0, 0.0, -3.0]);
	this.mvGlobalMatrixUniformID = this.gl.getUniformLocation(this.webglprogramobject, "um4GlobalModelviewMatrix");	
	this.gl.uniformMatrix4fv(this.mvGlobalMatrixUniformID, false, this.globalmvMatrix);

	this.localmvMatrix = new Float32Array
		([
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
		]);
	this.mvLocalMatrixUniformID = this.gl.getUniformLocation(this.webglprogramobject, "um4LocalModelviewMatrix");	
	this.gl.uniformMatrix4fv(this.mvLocalMatrixUniformID, false, this.localmvMatrix);
	
	this.pMatrix = this.GetPerspectiveMatrix();	
	this.pMatrixUniformID = this.gl.getUniformLocation(this.webglprogramobject, "um4PerspectiveMatrix");
	this.gl.uniformMatrix4fv(this.pMatrixUniformID, false, this.pMatrix);

	this.gl.clearColor(0.8, 0.8, 0.8, 1.0);
	this.gl.clearDepth(1.0);         // delete depth buffer value
	this.gl.enable(this.gl.DEPTH_TEST);
	this.gl.depthFunc(this.gl.LEQUAL);	//  smaller depth is front

	this.ntriangles = 0;	
	//vVertices[0] = -0.05;
	//this.gl.bufferData(this.gl.ARRAY_BUFFER, vVertices, this.gl.STATIC_DRAW);
	//this.gl.drawArrays(this.gl.TRIANGLES, 0, 3);	

}

	var palette = new Float32Array
	([
		0.0, 0.0, 0.5, 0,
		0.0, 0.0, 1.0, 0,
		0.0, 0.5, 1.0, 0,
		0.0, 1.0, 1.0, 0,
		0.5, 1.0, 0.5, 0,
		1.0, 1.0, 0.0, 0,
		1.0, 0.5, 0.0, 0,
		1.0, 0.0, 0.0, 0,
		0.5, 0.0, 0.0, 0
	]);


GL.prototype.GetShader = function(gl, id)
{  
	var shaderScript = document.getElementById(id);
	if (!shaderScript) return null;

	var str = "";
	var k = shaderScript.firstChild;
	while (k)      
	{
		if (k.nodeType == 3) str += k.textContent;
		k = k.nextSibling;
	}
	var shader;
	if (shaderScript.type == "x-shader/x-fragment")
	{
		shader = gl.createShader(gl.FRAGMENT_SHADER);
	} else
	if (shaderScript.type == "x-shader/x-vertex")
	{  
		shader = gl.createShader(gl.VERTEX_SHADER);
	}
	else
	{	  
		return null;
	}

	gl.shaderSource(shader, str);
	gl.compileShader(shader);
	if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS))
	{
		alert(gl.getShaderInfoLog(shader));
		return null;
	}
	return shader;
}
	
	
function GetColor(cf, color)
{
	if (cf < 0.) cf = 0.; else if (cf > 0.9999) cf = 0.9999;
	var idx = Math.floor(cf*8.);
	var g = cf*8. - idx;
//	for(var i = 0; i<4; i++)
	{
		color[0] = palette[idx*4 + 0]*(1.-g) + palette[idx*4 + 4]*g;
		color[1] = palette[idx*4 + 1]*(1.-g) + palette[idx*4 + 5]*g;
		color[2] = palette[idx*4 + 2]*(1.-g) + palette[idx*4 + 6]*g;
		//color[3] = palette[idx*4 + 3]*(1.-g) + palette[idx*4 + 7]*g;
	}
}

GL.prototype.Init = function(canvas)
{
	var names = [ "webgl", "experimental-webgl", "moz-webgl", "webkit-3d" ];
	var gl = null;
	for (var i=0; i<names.length; i++)
	{
		try
		{
			gl = canvas.getContext(names[i]);
			if (gl) break;
		} catch (e) { }
	}
	return gl;
}	


GL.prototype.BuildTriangle = function()
{
	this.ntriangles = 1;
	var vVertices = new Float32Array
		([ 
		0.,  1., 0.,
		-1., -1., 0.,
		1., -1., 0. 
		]);
	var vertexBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vVertices, this.gl.STATIC_DRAW);

	var vColors = new Float32Array
		([
		1.0, 0.0, 0.0, 1.0, // red
		0.0, 1.0, 0.0, 1.0, // green
		0.0, 0.0, 1.0, 1.0 // blue
		]); 
	var colorBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vColors, this.gl.STATIC_DRAW);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	var vertexAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avPosition");
	this.gl.vertexAttribPointer(vertexAttribLoc, 3, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(vertexAttribLoc);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	var colorAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avColor");
	this.gl.vertexAttribPointer(colorAttribLoc, 4, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(colorAttribLoc);
}

GL.prototype.AddQuad = function(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, c1, c2, c3, c4, vVertices, vColors, n)
{
	var color = new Array(4);
	vVertices[n*3+0] = x1;
	vVertices[n*3+1] = y1;
	vVertices[n*3+2] = z1;
	
	GetColor(c1, color);
	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;

	vVertices[n*3+0] = x2;
	vVertices[n*3+1] = y2;
	vVertices[n*3+2] = z2;

	GetColor(c2, color);
	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;
		
	vVertices[n*3+0] = x3;
	vVertices[n*3+1] = y3;
	vVertices[n*3+2] = z3;

	GetColor(c3, color);
	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;
// ----
	vVertices[n*3+0] = x3;
	vVertices[n*3+1] = y3;
	vVertices[n*3+2] = z3;

	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;

	vVertices[n*3+0] = x4;
	vVertices[n*3+1] = y4;
	vVertices[n*3+2] = z4;

	GetColor(c4, color);
	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;
		
	vVertices[n*3+0] = x1;
	vVertices[n*3+1] = y1;
	vVertices[n*3+2] = z1;

	GetColor(c1, color);
	vColors[n*4+0] = color[0];
	vColors[n*4+1] = color[1];
	vColors[n*4+2] = color[2];
	vColors[n*4+3] = 1.;
	n++;	
}

GL.prototype.AddLightQuad = function(x1, y1, z1,   x2, y2, z2,   x3, y3, z3,   x4, y4, z4,   vVertices, vColors, n)
{
	var v1x = x2 - x1;
	var v1y = y2 - y1;
	var v1z = z2 - z1;
	var v2x = x3 - x1;
	var v2y = y3 - y1;
	var v2z = z3 - z1;
	var nx = v1y*v2z - v1z*v2y;
	var ny = v1z*v2x - v1x*v2z;
	var nz = v1x*v2y - v1y*v2x;
	
	var c = (0.5 * ny/Math.sqrt(nx*nx + ny*ny + nz*nz)) + 0.5;

	vVertices[n*3+0] = x1;
	vVertices[n*3+1] = y1;
	vVertices[n*3+2] = z1;
	
	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;

	vVertices[n*3+0] = x2;
	vVertices[n*3+1] = y2;
	vVertices[n*3+2] = z2;

	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;
		
	vVertices[n*3+0] = x3;
	vVertices[n*3+1] = y3;
	vVertices[n*3+2] = z3;

	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;
// ----
	vVertices[n*3+0] = x3;
	vVertices[n*3+1] = y3;
	vVertices[n*3+2] = z3;

	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;

	vVertices[n*3+0] = x4;
	vVertices[n*3+1] = y4;
	vVertices[n*3+2] = z4;

	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;
		
	vVertices[n*3+0] = x1;
	vVertices[n*3+1] = y1;
	vVertices[n*3+2] = z1;

	vColors[n*4+0] = c;
	vColors[n*4+1] = c;
	vColors[n*4+2] = c;
	vColors[n*4+3] = 1.;
	n++;	
}



GL.prototype.Build2DMeshBox = function(width, height, thickness, Nx, Ny, map)
{
	this.ntriangles = (Nx-1) * (Ny-1) * 2;
	var vVertices = new Float32Array(this.ntriangles*3*3);
	var vColors = new Float32Array(this.ntriangles*3*4);

	var dx = width / (Nx-1);
	var dy = height/ (Ny-1);
	
	var n=0;
	for(var j=0; j<Ny-1; j++)
	for(var i=0; i<Nx-1; i++)
	{
		var x = dx * (i-(Nx-1)/2);
		var y = dy * (j-(Ny-1)/2);
/*
		this.AddQuad(
			x,    y, -thickness,
			x+dx, y, -thickness,
			x+dx, y+dy, -thickness,
			x,    y+dy, -thickness,
			map[i][j]*0.5+0.5, 
			map[i+1][j]*0.5+0.5, 
			map[i+1][j+1]*0.5+0.5, 
			map[i][j+1]*0.5+0.5, 
			vVertices, vColors, n);
			n += 6;
*/
		this.AddQuad(
			x,    y, 0.,
			x+dx, y, 0.,
			x+dx, y+dy, 0.,
			x,    y+dy, 0.,
			map[i][j]*0.5+0.5, 
			map[i+1][j]*0.5+0.5, 
			map[i+1][j+1]*0.5+0.5, 
			map[i][j+1]*0.5+0.5, 
			vVertices, vColors, n);
			n += 6;

		/*
		this.AddQuad(
			x,    y, -thickness/2,
			x+dx, y, -thickness/2, 
			x+dx, y+dy, -thickness/2, 
			x,    y+dy, -thickness/2, 
			1., 1., 1., 1, vVertices, vColors, n);
			n += 6;
		
		this.AddQuad(
			x,    y, thickness/2,
			x+dx, y, thickness/2, 
			x+dx, y+dy, thickness/2, 
			x,    y+dy, thickness/2, 
			1., 1., 1., 1, vVertices, vColors, n);
			n += 6;
		*/
	}
	
	var vertexBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vVertices, this.gl.STATIC_DRAW);
	
	var colorBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vColors, this.gl.STATIC_DRAW);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	var vertexAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avPosition");
	this.gl.vertexAttribPointer(vertexAttribLoc, 3, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(vertexAttribLoc);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	var colorAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avColor");
	this.gl.vertexAttribPointer(colorAttribLoc, 4, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(colorAttribLoc);
}

GL.prototype.BuildArrow = function()
{
	var N = 16;
	this.ntriangles = N * 4;
	var vVertices = new Float32Array(this.ntriangles*3*3);
	var vColors = new Float32Array(this.ntriangles*3*4);
	var n = 0;
	for(var i=0; i<N; i++)
	{
		var x = 0.3*Math.cos(1.*i / N * 2. * Math.PI);
		var y = 0.3*Math.sin(1.*i / N * 2. * Math.PI);
		var x2 = 0.3*Math.cos(1.*(i+1) / N * 2. * Math.PI);
		var y2 = 0.3*Math.sin(1.*(i+1) / N * 2. * Math.PI);

		this.AddLightQuad(
			-0.9, x2, y2,
			0.5, x2, y2,
			0.5, x, y,
			-0.9, x, y,
			vVertices, vColors, n);
			n += 6;

		this.AddLightQuad(
			0.5, 1.5*x, 1.5*y,
			0.5, 1.5*x2, 1.5*y2,
			1.5, 0., 0.,
			1.5, 0., 0.,
			vVertices, vColors, n);
			n += 6;
	}
	
	var vertexBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vVertices, this.gl.STATIC_DRAW);
	
	var colorBuffer = this.gl.createBuffer();
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	this.gl.bufferData(this.gl.ARRAY_BUFFER, vColors, this.gl.STATIC_DRAW);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, vertexBuffer);
	var vertexAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avPosition");
	this.gl.vertexAttribPointer(vertexAttribLoc, 3, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(vertexAttribLoc);
	
	this.gl.bindBuffer(this.gl.ARRAY_BUFFER, colorBuffer);
	var colorAttribLoc = this.gl.getAttribLocation(this.webglprogramobject, "avColor");
	this.gl.vertexAttribPointer(colorAttribLoc, 4, this.gl.FLOAT, false, 0, 0);
	this.gl.enableVertexAttribArray(colorAttribLoc);
}

var xangle = 0;
var yangle = 0;

GL.prototype.Clear = function()
{
	// clear background	and depth buffer
	this.gl.clear(this.gl.COLOR_BUFFER_BIT | this.gl.DEPTH_BUFFER_BIT);
}

GL.prototype.SetGlobalMatrix = function(xangle, yangle)
{
	this.globalmvMatrix = this.GetRotationMatrix(xangle, yangle, 0., 1.);
	this.globalmvMatrix[14] = -4.; // translation
	this.gl.uniformMatrix4fv(this.mvGlobalMatrixUniformID, false, this.globalmvMatrix);
}

GL.prototype.Draw = function()
{
	this.gl.uniformMatrix4fv(this.mvLocalMatrixUniformID, false, this.localmvMatrix);
	this.SetGlobalMatrix(xangle, yangle);
	this.gl.drawArrays(this.gl.TRIANGLES, 0, this.ntriangles * 3);
}


GL.prototype.GetPerspectiveMatrix = function()
{
	var pMatrix = new Float32Array
		([
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
		]);
		var znear = 0.1;
		var zfar = 100.;
		var fovy = 90.;
		var aspect = this.canvas.width/this.canvas.height;

		var top = znear * Math.tan(fovy * Math.PI / 360.0);
		var bottom = -top;
        
		var right = top * aspect;
		var left = -right;
		var near = znear;
		var far = zfar;	
		var rl = (right - left);
		var tb = (top - bottom);
		var fn = (far - near);
		pMatrix[0] = (near * 2) / rl;
		pMatrix[1] = 0;
		pMatrix[2] = 0;
		pMatrix[3] = 0;
		pMatrix[4] = 0;
		pMatrix[5] = (near * 2) / tb;
		pMatrix[6] = 0;
		pMatrix[7] = 0;
		pMatrix[8] = (right + left) / rl;
		pMatrix[9] = (top + bottom) / tb;
		pMatrix[10] = -(far + near) / fn;
		pMatrix[11] = -1;
		pMatrix[12] = 0;
		pMatrix[13] = 0;
		pMatrix[14] = -(far * near * 2) / fn;
		pMatrix[15] = 0;
		return pMatrix;
}

GL.prototype.GetRotationMatrix = function(ax, ay, az, s)
{
	var m = new Float32Array
		([
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
		]);
		ay = -ay;
		
	m[0] = (Math.cos(az)*Math.cos(ay) + Math.sin(ay)*Math.sin(ax)*Math.sin(az))*s;
	m[1] = (Math.cos(ay)*Math.sin(az) - Math.sin(ay)*Math.sin(ax)*Math.cos(az))*s;
	m[2] = Math.sin(ay)*Math.cos(ax)*s;
	
	m[4] = -Math.cos(ax)*Math.sin(az)*s;
	m[5] = Math.cos(ax)*Math.cos(az)*s;
	m[6] = Math.sin(ax)*s;
	
	m[8] = (-Math.sin(ay)*Math.cos(az) + Math.cos(ay)*Math.sin(ax)*Math.sin(az))*s;
	m[9] = (-Math.sin(ay)*Math.sin(az) - Math.cos(ay)*Math.sin(ax)*Math.cos(az))*s;
	m[10] = Math.cos(ay)*Math.cos(ax)*s;
	return m;
}



