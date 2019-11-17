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

#include "utils.hpp"
#include "AppConfig.hpp"
#include <unistd.h>

// PARAMETERS -------------------------------------------------------

void AppConfig::write(std::ostream& ots){
    ots << "\"params\"      : {" << std::endl;
    ots << " \t\"input\"   : \"" << inp << "\"," << std::endl;
    ots << " \t\"in_files\" : \"" << ifiles.size() << "\"," << std::endl;
    ots << " \t\"out_file\" : \"" << outf << "\"," << std::endl;
    ots << " \t\"histogram_file\" : \"" << histf << "\"," << std::endl;
    ots << " \t\"kvalue\"   : \"" << kv << "\"," << std::endl;
    ots << " \t\"lvalue\"   : \"" << l << "\"," << std::endl;
    ots << " \t\"rmqtype\"   : \"" << r << "\"," << std::endl;
    ots << " \t\"trim\"   : \"" << t << "\"," << std::endl;
    ots << " }," << std::endl;
}

void AppConfig::printHelp(std::ostream& ots){
    ots << "Usage " << app << " "
        << "-i <input file> "
        << "-o <output file> "
        << "[OPTIONS]" << std::endl << std::endl
        << "Available options:" << std::endl
        << "\t -i <file>  fasta/fastq input file" << std::endl
        << "\t -o <file>  output file" << std::endl
        << "\t -l <int>   minimum length of supermaximal repeats" << std::endl
        << "\t -k <int>   maximum number of mismatches" << std::endl
	<< "\t -r <int>   0 for O(n) space for rmqs. 1 for O(nlogn) space for rmqs. Default: 1" << std::endl
	<< "\t -t <int>   0 not to trim to optimize e-value and 1 to trim to optimize e-value. Default: 1" << std::endl
        << std::endl
        << "\t -h         help " << std::endl
        << std::endl;
}

AppConfig::AppConfig(int argc, char** argv){
    char c;
    const char* params = "a:i:f:o:l:k:r:t:pnhH";
    help = false;
    app = argv[0];
    inp = "";
    a = "";
    outf = "";
    kv = 1;
    l = 100;
    r = 1;
    t = 1;
    //method = 0;
    std::string fstr = "";
    histogram = false;
    histf= "";

    while ((c = getopt(argc, argv, params)) != -1) {
        switch (c) {
        case 'i':
            inp = optarg;
            break;
        case 'o':
            outf = optarg;
            break;
	case 'a':
            inp = optarg;
            break;
        case 'l':
            l = atoi(optarg);
            break;
        case 'h':
            help = true;
            break;
        case 'H':
            help = true;
            break;
        case 'k':
            kv = atoi(optarg);
            break;
        case 'n':
            method |= 1;
            break;
        case 'r':
            r = atoi(optarg);
            break;
 	case 't':
            t = atoi(optarg);
            break;
        }
    }
    //if(logf.length() == 0)
    //    logf = outf + ".log";
}

AppConfig::AppConfig(const std::vector<std::string>& files,
                     const std::string& of,
                     int64_t kval, int64_t ell, int64_t rmq_type, int64_t trim){
    app = "testApp";
    ifiles = files;
    outf = of;
    l = ell;
    kv = kval;
    r = rmq_type;
    t = trim;
    help = false;
}

bool AppConfig::validate(std::ostream& ots){
    bool validCfg = true;
    if(help){
        printHelp(ots);
        return false;
    }
    if(inp.length() == 0 && ifiles.size() == 0){
        ots << "Invalid input file " << inp << std::endl;
        validCfg = false;
    }

    if(outf.length() == 0){
        ots << "Invalid output file " << outf << std::endl;
        validCfg = false;
    }

    if(outf.length() > 0){
        ofs.open(outf, std::ofstream::out);
        if(!ofs.is_open()){
            ots << "Error opening output file : " << outf << std::endl;
            validCfg = false;
        }
    }


    if(histf.length() > 0){
        histfs.open(histf, std::ofstream::out);
        if(!histfs.is_open()){
            ots << "Error opening histogram file : " << histf << std::endl;
            validCfg = false;
        }
    }

    if(histogram){
        if(kv <= 0 ){
            ots << "Invalid k value (" << kv << std::endl;
            if(kv <= 0) kv = 1;
            ots << "Will use Default values of k = " << kv << std::endl;
        }
    }

    if(kv < 0){
        ots << "Invalid k value (" << kv
            << "). Using default value of 1" << std::endl;
        kv = 1;
    }

    if(l < 0){
        ots << "Invalid l value (" << l
            << "). Using default value of 100" << std::endl;
        l = 100;
    }


   if(r !=0 && r != 1){
        ots << "Invalid r value (" << r
            << "). Using default value of 0" << std::endl;
        r = 0;
    }

   if(t !=0 && t != 1){
        ots << "Invalid t value (" << r
            << "). Using default value of 1" << std::endl;
        t = 1;
    }
    if(!validCfg){
        printHelp(ots);
    }
    return validCfg;
}

AppConfig::~AppConfig(){
    if(ofs.is_open())
        ofs.close();
    if(histfs.is_open())
        histfs.close();
}
