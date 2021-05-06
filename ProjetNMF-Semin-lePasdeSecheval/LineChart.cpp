//
//  LineChart.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 04/05/2021.
//

#include "LineChart.hpp"

LineChart::LineChart() : ChartBase() {}

LineChart::LineChart(std::string title) : ChartBase(title) {}

#pragma region Boiler plate

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Boiler plate code
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Boiler plate code for the top of a line chart
static void writeTopBoilerPlate(std::ostream& out)
{
    out << ("<html>\n");
    out << ("<head>\n");
    out << ("<script type=\"text/javascript\" src=\"https://www.google.com/jsapi\"></script>\n");
    out << ("<script type=\"text/javascript\">\n");
    out << ("google.load('visualization', '1.0', {'packages':['corechart']});\n");
    out << ("google.setOnLoadCallback(drawChart);\n");
    out << ("function drawChart() {\n");
    
}

// Boiler plate code for the bottom of a line chart
static void writeBottomBoilerPlate(std::ostream& out, const std::string& title)
{
    out<<"var options = {title: '" << title << "', 'width':800, 'height':800};\n";
    out << "var chart = new google.visualization.LineChart(document.getElementById('chart_div'));\n";
    out << "chart.draw(data, options);\n";
    out << "}\n";
    out << "</script>\n";
    out << "</head>\n";
    out << "<body>\n";
    out << "<div id='chart_div'>\n";
    out << "</body>\n";
    out << "</html>";
}

// Write the data of the line chart
static void writeLineChartData(std::ostream& out, const std::vector<double>& xSeries, const std::vector<double>& ySeries)
{
    ASSERT(xSeries.size() == ySeries.size());
    out << "var data = google.visualization.arrayToDataTable([\n";
    out << "['x values','y values'],\n";
    size_t n = xSeries.size();
    for (size_t i = 0; i<n; ++i)
    {
        double x = xSeries[i];
        double y = ySeries[i];
        out << "[" << x << ", " << y << "]";
        if (i != n - 1)
        {
            out << ",";
        }
        out << "\n";
    }
    out << "]);\n";
}

#pragma endregion

void LineChart::setSeries(const std::vector<double>& x, const std::vector<double>& y)
{
    m_xSeries = x;
    m_ySeries = y;
}

void LineChart::writeAsHTMLToStream(std::ostream& out) const
{
    writeTopBoilerPlate(out);
    writeLineChartData(out, m_xSeries, m_ySeries);
    writeBottomBoilerPlate(out, m_title);
}

static void testWriteLineChart()
{
    std::vector<double> x;
    std::vector<double> y;
    for (int i = 0; i < 3; ++i)
    {
        x.push_back(i);
        y.push_back(i * i);
    }
    std::stringstream s;
    writeLineChartData(s, x, y);
    std::string actual = s.str();
    std::string expected = std::string("")
    + "var data = google.visualization.arrayToDataTable([\n"
    + "['x values','y values'],\n"
    + "[0, 0],\n"
    + "[1, 1],\n"
    + "[2, 4]\n"
    + "]);\n";
    std::cout << actual;
    std::cout << expected;
    ASSERT(actual == expected);
}

void testLineChart()
{
    TEST(testWriteLineChart);
}
