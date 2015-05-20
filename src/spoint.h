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
#ifndef __SPOINT__
#define	 __SPOINT__
#include <cmath>
/*!
 * \struct SPoint
 * \brief Structure used to represent a two dimensional point
 *
 * The point represented by this structure has coordinates of type double. No operators
 * are overloaded for this structure.
 */
struct SPoint
{
	SPoint() : x(-1), y(-1),energy(0)
	{}
	SPoint(double _x, double _y, double _energy) : x(_x), y(_y), energy(_energy)
	{}
	SPoint(const SPoint& other) : x(other.x), y(other.y), energy(other.energy)
	{}
	double squaredDistanceTo(const SPoint& other) const
	{
		return (x-other.x)*(x-other.x) + (y-other.y)*(y-other.y);
	}
	bool operator==(const SPoint &other) const
	{
		return (x == other.x) && (y == other.y);
	}
	bool operator<(const SPoint& other) const
	{
		return y>other.y || (!(other.y>y) && x<other.x);
	}
	SPoint operator+(const SPoint& other) const
	{
		SPoint result(0,0,0);
		result.x = this->x + other.x;
		result.y = this->y + other.y;
		return result;
	}
	SPoint operator-(const SPoint& other) const
	{
		SPoint result(0,0,0);
		result.x = this->x - other.x;
		result.y = this->y - other.y;
		return result;
	}
	SPoint operator*(const double other)
	{
		SPoint result;
		result.x = this->x * other;
		result.y = this->y * other;
		return result;
	}
	double dot(const SPoint& other) const
	{
		return x*other.x + y*other.y;
	}
	double cross(const SPoint& other) const
	{
		return x*other.y - y*other.x;
	}
	double len() const
	{
		return sqrt(x*x + y*y);
	}
	double x;
	double y;
	double energy;
};

#endif	//__SPOINT__
