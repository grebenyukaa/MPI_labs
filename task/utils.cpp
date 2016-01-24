#include "utils.h"

#include <fstream>
#include <sstream>

Logger::Logger(const Matrix::index_type nodeid)
{
    std::ostringstream oss;
    oss << "node" << nodeid << ".log";
    m_ofs = new std::ofstream(oss.str().c_str(), std::ios_base::app);
}

Logger::~Logger()
{
    m_ofs->close();
    delete m_ofs;
}
