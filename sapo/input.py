import sympy as sp
from sympy import parsing.sympy_parser.parse_expr

import configparser as conf
import os.path
from functools import map

import sapo.model

class InputFileParser:

    def __init__(self, path):

        self.config = conf.ConfigParser()
        self.config.read(path)


    def _createInitBundle(self):

        var_str = self.config['VARS']['vars'].split(', ')
        dyns_str = self.config['VARS']['dyns'].split(', ')

        vars = map(var_str, lambda x: sp.Symbol(x))
        dyns = map(dyns_str, lambda dx: parse_expr(dx, vars))

        temp_str_vals = dict(self.config['DIRECTIONS']).values()
        temp_str_list = list(temp_str_vals)

        for ptope_str in temp_str_list:

            ptope_expr = ptope_str.split(' , ')

            for direct_str in ptope_expr:

                direct_A, direct_b = direct_str[1:-1].split(', ')
                direct_A = parse_expr(direct_A,  vars)

                direct = sp.Poly(direct_A, vars).coeffs()
