from enum import Enum

'''
Basic logger utility module for debugging purposes.
'''

#Just basic write to file feature that nicely formats the relevant variable states and outputs them into a readable format.
#Just assume that the debugging info is already known beforehand.
#Quick and dirty logging feature. Not extensive or efficient by any means in its current state.


spaces = lambda x: ''.join('  ' for _ in range(x))

class Debug(Enum):
    MINMAX = 1 #For min/max point for parallelotopes
    POLY = 2 #For min/max polynomals calculated for bernstein expansion
    LOCAL_BOUND = 3 #Local bounds for each parallelotope for offl, offu
    GLOBAL_BOUND = 4 #Final bound after considering all parallelotopes

_Debug_Strings = {
Debug.MINMAX: 'Min/Max points for Parall {0}: {1}    {2}',
Debug.POLY:  'UPoly: {0}   LPoly: {1}',
Debug.LOCAL_BOUND: 'MaxB: {0}   MinB: {1}    for P: {2}',
Debug.GLOBAL_BOUND: 'New Offu: {0}    New Offl: {1}'
}

def write_log(*args):

    debug_cat = args[-1]
    with open('log.txt','a') as f:
        f.write('\n' + spaces(debug_cat.value) + _Debug_Strings[debug_cat].format(*args) + '\n')
