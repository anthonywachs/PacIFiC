#ifndef _CSVWRITER_HH_
#define _CSVWRITER_HH_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// =================================================================================================
/** @brief CSV Writer for benchmark results */
class CSVWriter
{
private:
    std::string   m_filename;
    std::ofstream m_file;
    bool          m_headerWritten;

public:
    CSVWriter(const std::string& filename, bool append = false)
        : m_filename(filename)
        , m_headerWritten(false)
    {
        if(append)
            m_file.open(m_filename, std::ios::out | std::ios::app);
        else
            m_file.open(m_filename, std::ios::out);
        if(!m_file.is_open())
        {
            std::cerr << "Error: Could not open CSV file: " << m_filename << std::endl;
        }
    }

    ~CSVWriter()
    {
        if(m_file.is_open())
            m_file.close();
    }

    void writeHeader(const std::vector<std::string>& columns)
    {
        if(!m_file.is_open() || m_headerWritten)
            return;

        for(size_t i = 0; i < columns.size(); ++i)
        {
            m_file << columns[i];
            if(i < columns.size() - 1)
                m_file << ",";
        }
        m_file << "\n";
        m_file.flush();
        m_headerWritten = true;
    }

    template <typename... Args>
    void writeRow(Args... args)
    {
        if(!m_file.is_open())
            return;

        writeRowImpl(args...);
        m_file << "\n";
        m_file.flush();
    }

private:
    template <typename T>
    void writeRowImpl(T value)
    {
        m_file << value;
    }

    template <typename T, typename... Args>
    void writeRowImpl(T value, Args... args)
    {
        m_file << value << ",";
        writeRowImpl(args...);
    }
};

#endif
