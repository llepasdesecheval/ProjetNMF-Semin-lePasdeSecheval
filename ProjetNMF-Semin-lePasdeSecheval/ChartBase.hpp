//
//  ChartBase.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 04/05/2021.
//

#if !defined CHARTBASE_HPP
#define CHARTBASE_HPP

#include "stdafx.hpp"

class ChartBase
{
public:
    /**
     Default constructor
     */
    ChartBase() : m_title("A chart") {}
    
    /**
     Constructor of ChartBase
     @param title Title of the graph
     */
    ChartBase(std::string title) : m_title(title) {}
    
    /**
     Virtual destructor
     */
    ~ChartBase() {}
    
    /**
     Title of the graph
     */
    std::string Title() const
    {
        return m_title;
    }
    
    /**
     Sets the title of the graph
     */
    void setTitle(std::string title);
    
    /**
     Writes the HTML code to generate the graph
     */
    virtual void writeAsHTMLToStream(std::ostream& out) const = 0;
    
    /**
     Writes the HTML file containing the graph
     */
    void writeAsHTML(const std::string& file) const;
    
protected:
    /**
    Title of the graph
    */
    std::string m_title;
};

#endif /* ChartBase_hpp */
