RandomLinearLab
---------------------------------------------------------------------------

1. Introduction
===============

This is a Lab for testing and comparing randomized methods for solving linear systems all implemented in MATLAB, based on method described in

	Robert M. Gower and Peter Richtarik, "Randomized Iterative Methods for Linear Systems", arXiv:1506.03296

To interface with the methods try
main_pd.n               % Test methods for positive definite
main_overdet            % Test methods for overdetermined systems


The methods implemented are

Randomized Kaczmarz
Coordinate Descent for Least Squares
Coordinate Descent for Positive Definite  (with a block version)
Gaussian Kaczmarz
Gaussian for Least Squares
Gaussian for Positive Definite  (with a block version)

Additionally, you can use

main_opt_CD-pd   	% Test CD-pd with optimized probability sampling
main_opt_Kaczmarz	% Test Kaczmarz with optimiof eachzed probability sampling

2. Installation and Setup
=========================

Start Matlab and make sure that the working directory is set to the
main directory of the present package.  At the MATLAB prompt, run

  >> setup_RandomLinearLab

To test if the installation and setup have been 
completed successfully please run in the MATLAB prompt:

  >> demo_RandomLinearLab


3. License
==========

 RandomLinearLab or Random Linear Lab
 Copyright (C) 2014, Robert Gower and Peter Richtarik 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.


4. The Authors
==============

If you have any bug reports or comments, please feel free to email 

  Robert Gower <r.m.gower@sms.ed.ac.uk>
  Robert Gower <gowerrobert@gmail.com>

Robert Gower and Peter Richtarik 
18 June 2015
