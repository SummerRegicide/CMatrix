#ifndef _GENERALDEFINATION_H_
#define _GENERALDEFINATION_H_

#include <iostream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include <cmath>

namespace LAP
{

	/*Q: can not put above headers here, or std goes wrong*/


	//Define special typename inside
	namespace inSDef
	{
		//indicate which line, order
		typedef unsigned int Line;

		//indicate positive or negative
		typedef bool Symbol;
		
		//indicate number
		typedef long double RealNumValue;

		//a power is non-negative
		typedef unsigned int Pow;

		//shortcut
		typedef std::vector<RealNumValue> Vec;
		typedef std::vector<std::vector<RealNumValue> > VecSet;
		typedef std::vector<RealNumValue>::iterator VecIter;
		typedef std::vector<std::vector<RealNumValue> >::iterator VecSetIter;
		typedef std::vector<RealNumValue>::const_iterator cVecIter;
		typedef std::vector<std::vector<RealNumValue> >::const_iterator cVecSetIter;

		//save
		typedef short int PosOrNeg;

		//counter is non-negative
		typedef unsigned int Counting;


		//end of inside definition
	}
	
	//end of LAP namespace
}
#endif