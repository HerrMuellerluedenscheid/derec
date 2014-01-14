import numpy as num
from collections import defaultdict

class OpticBase():
    def __init__(self, test_cases, keys={}):

        assert len(keys <= 2)
        assert map(lambda x: x.key==self.test_cases[0].key, self.test_cases)
        self.test_cases = test_cases
        self.plot_keys = keys
        self.dimensions = defaultdict(list)
        self.update_dimensions()

    def numpyrize(self):
        '''
        Transform self.test_cases into a 2-/3-dimensional numpy array.
        '''
        A = num.zeros(self.dimensions)
        y = 0
        for case in self.test_cases:
            x_val = getatr(case.event, x)
            x = self.value_to_index(self.plot_keys[x], x_val)
            
            if not y == '':
                y_val = getatr(case.event, y)
                y = self.value_to_index(self.plot_keys[y], y_val)
            
            A[x,y,:] = case.get_misfit_array()

    def update_dimensions(self, new_cases=None):
        if new_cases:
            cs = new_cases
        else: 
            cs = self.test_cases
        
        dim_mapper = defaultdict(set)

        for k in self.plot_keys: 
            for c in cs:
                dim_mapper[k].add(getatr(c,k))
            dim_mapper[k]=list(dim_mapper[k])

    def value_to_index(k, val):
        return self.dim_mapper[k].index(val)

    def add_cases(self, cases):
        self.test_cases.extend(cases)
        self.update_dimensions(cases)
