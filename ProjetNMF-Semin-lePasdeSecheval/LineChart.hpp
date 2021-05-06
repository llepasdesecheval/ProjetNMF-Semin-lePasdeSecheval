//
//  LineChart.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 04/05/2021.
//

#if !defined LINECHART_HPP
#define LINECHART_HPP

#include "ChartBase.hpp"

class LineChart : public ChartBase
{
public:
    /**
     Default constructor for LineChart
     */
    LineChart();
    
    /**
     Constructor for LineChart
     @param title Title of the chart
     */
    LineChart(std::string title);
    
    /**
     Sets the data of the chart
     @param x Input data
     @param y Output data
     */
    void setSeries(const std::vector<double>& x, const std::vector<double>& y);
    
    void writeAsHTMLToStream(std::ostream& out) const override;
    
private:
    std::vector<double> m_xSeries;
    std::vector<double> m_ySeries;
};

void testLineChart();

#endif /* LineChart_hpp */
