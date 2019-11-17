/**
    ALFRED: A practical method for alignment-free distance computation.
    Copyright (C) 2015 Thankachan, S.V., Chockalingam, S.P., Liu, Y., 
    Apostolico, A. and Aluru, S

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

#ifndef APPCONFIG_H
#define APPCONFIG_H

#include <vector>
#include <string>
#include <fstream>

void init_substitution_score_tables ( void );
double gettime( void );

// PARAMETERS   ------------------------------------------------------
struct AppConfig{
    std::string inp;
    std::string a;
    std::vector<std::string> ifiles;
    std::string outf;
    //std::string logf;
    std::string app;
    bool help;
    std::ofstream ofs;
    std::ofstream lfs;
    std::ofstream histfs;
    int64_t kv;
    int64_t l;
    int64_t r;
    int64_t t;
    int64_t method;
    bool histogram;
    std::string histf;

    void printHelp(std::ostream& ots);
    bool validate(std::ostream& ots);
    void write(std::ostream& ots);
    AppConfig(const std::vector<std::string>& files,
              const std::string& of,
              int64_t kval, int64_t ell, int64_t rmq_type, int64_t trim);
    AppConfig(int argc, char** argv);
    ~AppConfig();
};
int64_t nchoosek( int64_t n, int64_t k );

#endif /* APPCONFIG_H */
