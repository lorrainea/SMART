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

#ifndef READSDB_H
#define READSDB_H

#include <string>
#include <vector>

// TYPEDEFS ----------------------------------------------------------
typedef std::vector<std::string> preads_t;
typedef std::vector<preads_t>    preads_vector_t;
typedef std::vector<unsigned>    rlengths_t;
typedef std::vector<std::string> rnames_t;
typedef std::vector<rlengths_t>  rlengths_vector_t;
typedef std::vector<rnames_t>    rnames_vector_t;
//  READS DATABASE ---------------------------------------------------
//  - Encapsulates the read set as a database
//  - Provides interface to the database
//  - Constructed/Initialized by a directory containing ONLY fasta files

class ReadsDB{
    std::vector<std::string> m_readFiles; // file paths
    std::vector<unsigned>    m_readCounts; // reads cnts in each file
    std::vector<unsigned>    m_readCtsPfxSum; // prefix sum of read cnts
    rlengths_vector_t        m_readLengths; // lengths of the read
    rnames_vector_t          m_readNames;   // names of the read
    preads_vector_t          m_readsStore; // the entire read set
    unsigned                 m_readMaxLength; // length of the maximum read
    unsigned                 m_nReads; // total number of reads
    std::string              m_dbString; // whole db as a string
   
    unsigned findMax();
    unsigned findMin();
    void padStrings(const unsigned& mxLength, const char& pChar);
    void init(char padChar, std::string& alphabet);
public:

    void writeMeta(std::ostream& ots);
    void writeNames(std::ostream& ots);
    const std::string& getDBS(){ return m_dbString; };
    // read fasta file
    static void addRead(std::string& read, std::string& header,
                        preads_t& freads, rlengths_t& flengths,
                        rnames_t& fnames);
    static void readFasta(std::string& fname, preads_t& freads,
                          rlengths_t& flengths, rnames_t& fnames, std::string& alphabet);
    static void readFastq(std::string& fname, preads_t& freads,
                          rlengths_t& flengths, rnames_t& fnames);
    unsigned padReads(const char& fillChar);
    // read access functions
    unsigned getReadLength(const unsigned& fileID, const unsigned& readID);
    int getFileOfRead(const unsigned& readID) const;
    unsigned getReadsCount(){ return m_nReads; };
    unsigned getFilesCount(){ return m_readFiles.size(); };
    const std::string& getFileName(const unsigned& idx);
    const std::string& getRead(const unsigned& fileID, const unsigned& readID);

    // Functions to retrieve information, when
    //  the reads are stored as string
    void asString(std::string& total);
    void buildDBS(); // converts the database into a single string

    // getDBSXXX() functions assume that buildDBS() is completed
    unsigned getDBSFileId(const unsigned& dbStrPosition);
    unsigned getDBSReadId(const unsigned& dbStrPosition);
    unsigned getDBSReadPos(const unsigned& dbStrPosition);
    char getDBSCharBefore(const unsigned& dbStrPosition);
    char getDBSCharAt(const unsigned& dbStrPosition);
    unsigned getReadId(const unsigned& fileID, const unsigned& readID);
    const std::string& getReadById(const unsigned& grID) const;
    const std::string& getReadNameById(const unsigned& grID) const;
    unsigned getReadLengthById(const unsigned& grID);
    ReadsDB(std::vector<std::string>& fnames, char padChar = (char) 0);
    ReadsDB(std::string& inpname, char padChar = (char) 0);
    ReadsDB(std::string& inpname, std::vector<std::string>& fnames,
            char padChar = (char) 0, std::string& alphabet = (std::string&) "" );
};

#endif /* READSDB_H */
