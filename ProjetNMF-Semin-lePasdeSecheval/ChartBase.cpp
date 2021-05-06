//
//  ChartBase.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 04/05/2021.
//

#include "ChartBase.hpp"

void ChartBase::setTitle(std::string title)
{
    m_title = title;
}

void ChartBase::writeAsHTML(const std::string& file) const
{
    std::ofstream out;
    out.open(file.c_str());
    writeAsHTMLToStream(out);
    out.close();
}
