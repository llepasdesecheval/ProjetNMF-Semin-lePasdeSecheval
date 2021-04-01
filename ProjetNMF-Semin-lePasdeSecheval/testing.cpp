//
//  testing.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 22/03/2021.
//

#include "testing.hpp"

/**
 Whether debug messages are enabled
 */
static bool debugEnabled = false;

bool isDebugEnabled()
{
    return debugEnabled;
}

void setDebugEnabled(bool enable)
{
    debugEnabled = enable;
}
