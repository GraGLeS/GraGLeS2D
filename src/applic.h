/*
	Utility math library
    Copyright (C) 2015  Markus Kuehbach

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
#ifndef _applic_h_
#define _applic_h_

#define DEBUG


#define _STR(x) _VAL(x)
#define _VAL(x) #x
#define ASSERT(cond)
#define TOLFP(x) ((1.0)-(0.99999/((double) x))) //tolerance
#define ERRTXT(text) (text" : file: "__FILE__" line:"_STR(__LINE__)"\n")

#define NOW true
#define NOTYET false

#define NEVER true


#ifndef DEBUG
        #define QUICKASSERT(cond)       ((void)0)
#else
        #define QUICKASSERT(cond)       ASSERT(cond)
#endif

void exitus (const char *s);
#endif
