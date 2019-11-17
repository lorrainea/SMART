/**
    SMART: Supermaximal Approximate Repeats Tool.
    Copyright (C) 2019 Ayad, L.A.K., Charalampopoulos, P. and Pissis, S.P. 

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
**/

#include "ReadsDB.hpp"
#include "AppConfig.hpp"
#include "compute_klcp.hpp"
#include "defs.hpp"
#include <sys/time.h>
#include <iostream>

using namespace std;


int64_t nchoosek(int64_t n, int64_t k ) 
{
	if (n < k) 
		return 0;
        if (n == k) 
		return 1;
        int s = min(k, (n - k));
        int t = n;
        for (int a = n - 1, b = 2; b <= s; a--, b++)
            t = (t * a) / b;

return t;
}

double gettime( void )
 {
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
 }

int main(int argc, char** argv)
{
	AppConfig acfg(argc, argv);
	if(!acfg.validate(std::cout))
		return 0;
	double start = gettime();

	init_substitution_score_tables ();

	// Build the reads database and get it as a string
	ReadsDB rdb(acfg.inp, acfg.ifiles, 0, acfg.a);

	if(rdb.getReadsCount() == 0) 
	{
		std::cout << "\"error\" : \"Input file not found or not in fasta/fa format\"" << std::endl;
		return 0;
	}

	compute_klcp(rdb, acfg);

	double end = gettime();
	// write time taken
	fprintf( stderr, "Elapsed time: %lf secs\n", ( end - start ) );

return 0;
}
