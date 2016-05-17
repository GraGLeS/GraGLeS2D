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

struct SPoint {
	SPoint() :
		x(-1), y(-1), energy(0), mob(1) {
	}
	SPoint(double _x, double _y, double _energy, double _mob) :
		x(_x), y(_y), energy(_energy), mob(_mob) {
	}
	SPoint(const SPoint& other) :
		x(other.x), y(other.y), energy(other.energy), mob(other.mob) {
	}
	double squaredDistanceTo(const SPoint& other) const {
		return (x - other.x) * (x - other.x) + (y - other.y) * (y - other.y);
	}
	double DistanceTo(const SPoint& other) const {
			return sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y));
		}
	inline bool operator==(const SPoint &other) const {
		return (x == other.x) && (y == other.y);
	}
	inline bool operator<(const SPoint& other) const {
		return y > other.y || (!(other.y > y) && x < other.x);
	}
	inline SPoint operator+(const SPoint& other) const {
		SPoint result(0, 0, 0, 1);
		result.x = this->x + other.x;
		result.y = this->y + other.y;
		return result;
	}
	inline SPoint operator-(const SPoint& other) const {
		SPoint result(0, 0, 0, 1);
		result.x = this->x - other.x;
		result.y = this->y - other.y;
		return result;
	}
	inline SPoint operator*(const double other) {
		SPoint result;
		result.x = this->x * other;
		result.y = this->y * other;
		return result;
	}
	inline double dot(const SPoint& other) const {
		return x * other.x + y * other.y;
	}
	inline double cross(const SPoint& other) const {
		return x * other.y - y * other.x;
	}
	inline double len() const {
		return sqrt(x * x + y * y);
	}
	inline double lenSqr() const
        {
                return (x*x + y*y);
        }

	double x;
	double y;
	double energy;
	double mob;
};

#endif	//__SPOINT__
