//
//  BSMModel.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "BSMModel.hpp"

BSMModel::BSMModel(double stockPrice, double volatility, double rate, double dividend, double date) :
m_stockPrice(stockPrice), m_volatility(volatility), m_interestRate(rate), m_dividendYield(dividend), m_date(date) {}
