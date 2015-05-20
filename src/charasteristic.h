/*
	GraGLeS 2D A grain growth simulation utilizing level set approaches
    Copyright (C) 2015  Christian Miessen, Nikola Velinov

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
*/
#ifndef	__CHARACTERISTIC__
#define __CHARACTERISTIC__

/*!
 * \struct characteristics
 * \brief Structure that holds information about a neighbor of a grain.

 * This structure holds the length of the boundary between the two grains,
 * the mis-orientation between the two grains, the energy density of the boundary
 * and its mobility.
*/
struct characteristics{
    LSbox* directNeighbour;
    double length;
    double energyDensity;
    double mis_ori;
	double mobility;
    characteristics(LSbox* directNeighbour, double length, double energyDensity, double mis_ori) : directNeighbour(directNeighbour), length(length), energyDensity(energyDensity), mis_ori(mis_ori), mobility(1)
	{}
	characteristics(LSbox* directNeighbour, double length, double energyDensity, double mis_ori, double mu) : directNeighbour(directNeighbour), length(length), energyDensity(energyDensity), mis_ori(mis_ori), mobility(mu)
	{}
	characteristics( const characteristics& other ) :
		directNeighbour(other.directNeighbour), length(other.length), energyDensity(other.energyDensity), mis_ori(other.mis_ori), mobility(other.mobility)
	{}
};

#endif
