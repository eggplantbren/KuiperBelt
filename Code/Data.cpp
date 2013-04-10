#include "Data.h"
#include <fstream>
#include <iostream>

using namespace std;

Data Data::instance;

Data::Data()
:loaded(false)
{

}

void Data::load(const char* filename)
{
	// Load the data
	fstream fin(filename, ios::in);
	if(!fin)
		cerr<<"Error: Can't load "<<filename<<"."<<endl;
	char hash;
	fin>>hash; fin>>numImages;
	times.resize(numImages);
	for(int i=0; i<numImages; i++)
		fin>>times[i];
	fin.ignore(1000000, '\n');
	fin>>hash>>ni>>nj>>xMin>>xMax>>yMin>>yMax; fin.ignore(1000000, '\n');

	xRange = xMax - xMin;
	yRange = yMax - yMin;
	dx = xRange/nj;
	dy = yRange/ni;
	dA = dx*dy;

	images.resize(numImages);

	for(int k=0; k<numImages; k++)
	{
		images[k].resize(ni, nj);
		for(int i=0; i<ni; i++)
			for(int j=0; j<nj; j++)
				fin>>images[k](i, j);
	}
	fin.close();
	loaded = true;

	// Assign the xc, yc arrays
	vector<double> x(nj);
	vector<double> y(ni);
	for(int j=0; j<nj; j++)
		x[j] = xMin + (j + 0.5)*dx;
	for(int i=0; i<ni; i++)
		y[i] = yMax - (i + 0.5)*dy;
	xc = images[0];
	yc = images[0];
	for(int i=0; i<ni; i++)
	{
		for(int j=0; j<nj; j++)
		{
			xc(i, j) = x[j];
			yc(i, j) = y[i];
		}
	}

}

