/*

   Copyright (c) 2011, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#include "receptor.hpp"
#include "file.hpp"

namespace idock
{
	receptor::receptor(vector<atom>&& atoms_) : atoms(static_cast<vector<atom>&&>(atoms_))
	{
		// Assert the Rvalue constructor of vector<T> work.
		BOOST_ASSERT(atoms_.empty());
		BOOST_ASSERT(!atoms.empty());

		// Dump receptor.
		//ofile dump("receptor.csv");
		//dump << "i,x,y,z,ad,xs\n";
		//for (size_t i = 0; i < atoms.size(); ++i)
		//{
		//	const atom& a = atoms[i];
		//	dump << i << ',' << a.coordinate[0] << ',' << a.coordinate[1] << ',' << a.coordinate[2] << ',' << a.ad << ',' << a.xs << '\n';
		//}
		//dump.close();
	}
}
