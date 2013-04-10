/*
    Copyright (C) 2011 Brendon J. Brewer
    This file is part of DNest, the Diffusive Nested Sampler.

    DNest is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    DNest is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with DNest.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Utils.h"
#include "RandomNumberGenerator.h"
#include "StarFieldModel.h"
#include <cmath>

using namespace DNest3;


const int StarFieldModel::maxNumStars = 200;


StarFieldModel::StarFieldModel()
:mockImage(Data::get_instance().get_ni(), Data::get_instance().get_nj())
,staleness(0)
{
	if(!Data::get_instance().isLoaded())
		std::cerr<<"WARNING: Data not loaded."<<std::endl;
}


void StarFieldModel::fromPrior()
{
	noiseSigma = exp(log(1E-3) + log(1E6)*randomU());
	noiseCoeff = exp(log(1E-3) + log(1E6)*randomU());
	background = -1000. + 2000.*randomU();

	psf.fromPrior();
	hyperparameters.fromPrior();
	numStars = randInt(maxNumStars + 1);

	stars.clear();
	for(int i=0; i<numStars; i++)
		stars.push_back(hyperparameters.generateStar());

	xc = Data::get_instance().get_xMin() + Data::get_instance().get_xRange()*randomU();
	yc = Data::get_instance().get_yMin() + Data::get_instance().get_yRange()*randomU();
	r = exp(log(0.01) + log(100.)*randomU());

	calculateMockImage();
	calculateLogLikelihood();
}


double StarFieldModel::perturb()
{
	double logH = 0.;
	int which;

	if(randomU() <= 0.05)
		which = 4;
	else
		which = randInt(6);

	if(which == 0)
		logH += perturb1();
	else if(which == 1)
		logH += perturb2();
	else if(which == 2)
		logH += perturb3();
	else if(which == 3)
		logH += perturb4();
	else if(which == 4)
		logH += perturb5();
	else if(which == 5)
		logH += perturb6();
	calculateLogLikelihood();

	return logH;
}


double StarFieldModel::perturb4()
{
	double logH = 0.;
	int which = randInt(3);
	if(which == 0)
	{
		noiseSigma = log(noiseSigma);
		noiseSigma += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		noiseSigma = mod(noiseSigma - log(1E-3), log(1E6)) + log(1E-3);
		noiseSigma = exp(noiseSigma);
	}
	else if(which == 1)
	{
		noiseCoeff = log(noiseCoeff);
		noiseCoeff += log(1E6)*pow(10., 1.5 - 6.*randomU())*randn();
		noiseCoeff = mod(noiseCoeff - log(1E-3), log(1E6)) + log(1E-3);
		noiseCoeff = exp(noiseCoeff);
	}
	else
	{
		mockImage.decrement(background);
		background += 2000.*pow(10., 1.5 - 6.*randomU())*randn();
		background = mod(background + 1000., 2000.) - 1000.;
		mockImage.increment(background);
		staleness++;
	}

	return logH;
}


double StarFieldModel::perturb5()
{
	double logH = 0.;
	logH += psf.perturb();
	calculateMockImage();
	return logH;
}



double StarFieldModel::perturb6()
{
	double logH = 0.;

	xc += Data::get_instance().get_xRange()*pow(10., 1.5 - 6.*randomU())*randn();
	xc = mod(xc - Data::get_instance().get_xMin(), Data::get_instance().get_xRange()) + Data::get_instance().get_xMin();

	yc += Data::get_instance().get_yRange()*pow(10., 1.5 - 6.*randomU())*randn();
	yc = mod(yc - Data::get_instance().get_yMin(), Data::get_instance().get_yRange()) + Data::get_instance().get_yMin();

	r = log(r);
	r += log(100.)*pow(10., 1.5 - 6.*randomU())*randn();
	r = mod(r - log(0.01), log(100.)) + log(0.01);
	r = exp(r);

	return logH;
}


double StarFieldModel::perturb1()
{
	double logH = 0.;

	// Make a proposal for the new number of stars
	int diff = static_cast<int>
			(round(maxNumStars*pow(10., 1.5 - 6.*randomU())*randn()));
	if(diff == 0)
		diff = (randomU() <= 0.5)?(-1):(1);
	int proposal = numStars + diff;
	proposal = mod(proposal, maxNumStars + 1);

	int actual_diff = proposal - numStars;

	if(actual_diff > 0)
	{
		for(int i=0; i<actual_diff; i++)
		{
			Star star = hyperparameters.generateStar();
			star.incrementImage(mockImage, psf);
			stars.push_back(star);
			numStars++;
		}
	}
	else if(actual_diff < 0)
	{
		int which;
		for(int i=0; i<-actual_diff; i++)
		{
			which = randInt(numStars);
			stars[which].decrementImage(mockImage, psf);
			stars.erase(stars.begin() + which);
			numStars--;
		}
	}

	staleness++;
	return logH;
}


double StarFieldModel::perturb2()
{
	double logH = 0.;
	int which = randInt(2);
	if(which == 0)
	{
		logH = hyperparameters.perturb1(stars);
		calculateMockImage();
	}
	else
		logH = hyperparameters.perturb2(stars);
	return logH;
}


double StarFieldModel::perturb3()
{
	double logH = 0.;

	double chance = pow(10., 0.5 - 4.*randomU());
	double scale = pow(10., 1.5 - 6.*randomU());

	int which = randInt(2);

	if(which == 0)
	{
		// Positions
		for(int i=0; i<numStars; i++)
		{
			if(randomU() <= chance)
			{
				if(chance < 1.)
					stars[i].decrementImage(mockImage, psf);

				logH += hyperparameters.perturbStar1(stars[i], scale);

				if(chance < 1.)
					stars[i].incrementImage(mockImage, psf);
			}
		}
	}
	else if(which == 1)
	{
		// Fluxes
		for(int i=0; i<numStars; i++)
		{
			if(randomU() <= chance)
			{
				if(chance < 1.)
					stars[i].decrementImage(mockImage, psf);

				logH += hyperparameters.perturbStar2(stars[i], scale);

				if(chance < 1.)
					stars[i].incrementImage(mockImage, psf);
			}
		}

	}

	if(chance < 1.)
		staleness++;
	else
		calculateMockImage();

	return logH;
}


void StarFieldModel::calculateMockImage()
{
	mockImage.set(background);
	for(size_t i=0; i<stars.size(); i++)
		stars[i].incrementImage(mockImage, psf);
	staleness = 0;
}


void StarFieldModel::calculateLogLikelihood()
{
	logL = 0.;
	double var;

	for(int k=0; k<Data::get_instance().get_numImages(); k++)
	{
		Array mock = mockImage; // This is going to have the KBO in it too

		for(int i=0; i<Data::get_instance().get_ni(); i++)
			for(int j=0; j<Data::get_instance().get_nj(); j++)
			{
				var = pow(noiseSigma, 2)
					+ noiseCoeff*(mock(i, j) - background);

				logL += -0.5*log(2*M_PI*var)
					- 0.5*pow(Data::get_instance().get_image(0)(i, j)
					- mockImage(i, j), 2)/var;
			}
	}
}


void StarFieldModel::print(std::ostream& out) const
{
	out<<numStars<<' '<<staleness<<' ';
	out<<psf<<' '<<noiseSigma<<' '<<noiseCoeff<<' '<<background<<' ';
	hyperparameters.print(out); out<<' ';
	out<<xc<<' '<<yc<<' '<<r<<' ';

	// Print x, pad with zeros
	for(int i=0; i<numStars; i++)
		out<<stars[i].x<<' ';
	for(int i=numStars; i<maxNumStars; i++)
		out<<0<<' ';

	// Print y, pad with zeros
	for(int i=0; i<numStars; i++)
		out<<stars[i].y<<' ';
	for(int i=numStars; i<maxNumStars; i++)
		out<<0<<' ';

	// Print flux, pad with zeros
	for(int i=0; i<numStars; i++)
		out<<stars[i].flux<<' ';
	for(int i=numStars; i<maxNumStars; i++)
		out<<0<<' ';

	for(int i=0; i<Data::get_instance().get_ni(); i++)
		for(int j=0; j<Data::get_instance().get_nj(); j++)
			out<<mockImage(i, j)<<' ';
}

