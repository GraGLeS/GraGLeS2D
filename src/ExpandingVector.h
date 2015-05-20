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
#ifndef __EXPANDING_VECTOR__
#define __EXPANDING_VECTOR__
template<class T>
/*!
 * \class ExpandingVector
 * \brief Class that implements a vector that only expands when resized.
 */
class ExpandingVector : public std::vector<T>
{
	typedef size_t size_type;
	typedef T value_type;
public:
	/*!
	 * \brief This method is used to expand the capacity of the vector to a certain size.
	 */
	void expand(size_type n)
	{
		if(n > this->size())
		{
			this->resize(n);
		}
	}
};
#endif	//__EXPANDING_VECTOR__
