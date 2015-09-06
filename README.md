# MISO
Online Micromagnetic Simulator for demonstration purposes

Just want to find out what's this project is all about? [Click here](https://s-macke.github.io/MISO/index.html)

Micromagnetism is an approximation of the magnetic behavior of ferromagnets on the nanometer scale averaging the magnetization over several nm.
It is the world of magnetic domain walls and spin waves and describe ground state properties as well as dynamics.
One of the most famous applications of this theory take place in hard drives in which the magnetization of tiny particles are switched via an external field to save bits 0 and 1.

This simulator allows you to observe the dynamics of such tiny magnetic particles. It's main purpose is for teaching and demonstrations.

You can find a detailed description of the theory on https://en.wikipedia.org/wiki/Micromagnetics
The numerical details are very similar to the OOMMF package: http://math.nist.gov/oommf/
and is based on a 2D finite difference grid. The stray field is calculated via a Fast Fourier Transformation (FFT)

### Simulator

* Main [simulator](https://s-macke.github.io/MISO/index.html) page

### LICENSE
 * The program is distributed under the terms of the MIT License. The license details can be found in the file "LICENSE.md"

### Developer
Sebastian Macke [simulationcorner.net](http://simulationcorner.net)
