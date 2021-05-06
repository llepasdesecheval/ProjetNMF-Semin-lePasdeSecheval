//
//  BSMModel.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined BSMMODEL_HPP
#define BSMMODEL_HPP

#include "stdafx.hpp"
#include "matlib.hpp"

class BSMModel
{
public:
    /**
     Class constructor for BSMModel
     @param stockPrice Initial stock price
     @param volatility Volatility of the stock price process
     @param rate Constant risk-free interest rate; default is 0
     @param dividend Constand dividend yield; default is 0
     @param date Date of the model; default is 0
     */
    BSMModel(double stockPrice, double volatility, double rate = 0, double dividend = 0, double date = 0);

    /**
    Initial stock price
    */
    double StockPrice() const
    {
        return m_stockPrice;
    }
    
    /**
    Stock returns' volatility
    */
    double Volatility() const
    {
        return m_volatility;
    }
    
    /**
    Constant risk-free interest rate
    */
    double InterestRate() const
    {
        return m_interestRate;
    }
    
    /**
     Stock continuous dividend yield
     */
    double DividendYield() const
    {
        return m_dividendYield;
    }
    
    /**
     Model's (start) date
     */
    double Date() const
    {
        return m_date;
    }
    
private:
    double m_stockPrice;
    double m_volatility;
    double m_interestRate;
    double m_dividendYield;
    double m_date;
};

void testBSM();

#endif /* BSMModel_hpp */
