from copy import deepcopy
from itertools import pairwise
from cutgeneratingfunctionology.igp import *
import csv
import os
from cutgeneratingfunctionology.spam.basic_semialgebraic import EmptyBSA

class RepElemGenFailure(Exception):
    pass


def add_breakpoints_and_find_equiv_classes(bkpt_poly):
    """
    Takes dim k-1 breakpoint NNC polyhedron (as a :class:`BasicSemialgebraicSet_base`) and finds rep elements 
    """
    # BSAs are highly mutable, work only with copies.
    B_cap_N_b = copy(bkpt_poly)
    B_cap_N_b.add_space_dimensions_and_embed(1)
    # get new number of breakpoints
    k = B_cap_N_b.ambient_dim()
    # if k< 2:
    #     raise ValueError("bkpt_poly should have space dim at least 1.")
    model_bound_bkpts = [0]*k
    model_bound_bkpts[k-1] = 1
    # 0 < lambda_k <1
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, -1, operator.lt) # model bounds
    B_cap_N_b.add_linear_constraint(model_bound_bkpts, 0, operator.gt) # model bounds 
    bkpt_order = [0]*k
    bkpt_order[k-2] = 1
    bkpt_order[k-1] = -1
    B_cap_N_b.add_linear_constraint(bkpt_order, 0, operator.lt) # order on bkpts
    # print(B_cap_N_b)
    # rep elem list 
    rep_elems = []
    for j in range(k-1):
        for i in range(k):
            for interval_w in [0,1]:
                for line_w in [0,1]:
                    # which interval is (lambda_k,lambda_k) located in?
                    # modeled lambda_i op 2lambda_k - w < lambda_{i+1}
                    for interval_op in [operator.lt, operator.eq, operator.gt]:
                        for line_op in [operator.lt, operator.eq, operator.gt]:
                            # highly mutable objects, operater on the copy
                            B_cap_N_b_copy = copy(B_cap_N_b)
                            lhs_i = [0]*k
                            lhs_i[k-1] = -2
                            lhs_i[i] = 1
                            B_cap_N_b_copy.add_linear_constraint(lhs_i, interval_w, interval_op)
                            lhs_i_plus_1 = [0]*k
                            lhs_i_plus_1[k-1] = -2
                            if i < k-1:
                                lhs_i_plus_1[i+1] = 1            
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, interval_w, operator.gt)
                            else:
                                B_cap_N_b_copy.add_linear_constraint(lhs_i_plus_1, interval_w + 1, operator.gt)
                            if not B_cap_N_b_copy.is_empty():
                                # does the line x+y equiv lambda_k mod 1 lie on/above/below (lambda_j,lambda_j)?
                                # modeled by  2lambda_j op lambda_k + w
                                lhs_j = [0]*k
                                lhs_j[j] = 2
                                lhs_j[k-1] = -1
                                B_cap_N_b_copy.add_linear_constraint(lhs_j, -line_w, line_op)
                                try:
                                    rep_elem = B_cap_N_b_copy.find_point()
                                    rep_elems.append(tuple(rep_elem))
                                except EmptyBSA:
                                    pass
    return unique_list(rep_elems)

def nnc_poly_from_bkpt_sequence(bkpt, backend=None):
    n = len(bkpt)
    # assert(n >= 2)
    coord_names = []
    bkpt_vals = bkpt
    vals = bkpt_vals[0:n]
    bkpt_extd = list(bkpt)+[1]
    for i in range(0,n):
        coord_names.append('lambda'+str(i))
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True)
    logging.disable(logging.INFO)
    K.gens()[0] == 0
    for i in range(n-1):
        K.gens()[i] < K.gens()[i+1]
    K.gens()[n-1] < 1
    for i in range(n):
        for j in range(n):
            if bkpt[i]+bkpt[j]>= 1:
                w = 1
            else:
                w = 0
            for k in range(n):
                if bkpt_extd[k] < bkpt[i]+bkpt[j] - w and bkpt[i]+bkpt[j] - w < bkpt_extd[k+1]:
                    if k != n-1:
                        K.gens()[k] < K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < K.gens()[k+1]
                    else:
                        K.gens()[k] < K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < 1
                elif bkpt_extd[k] == bkpt[i]+bkpt[j] - w:
                    if k != n-1:
                        K.gens()[k] == K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < K.gens()[k+1]
                    else:
                        K.gens()[k] == K.gens()[i] + K.gens()[j] - w
                        K.gens()[i] + K.gens()[j] - w < 1                    
    return K._bsa


def make_rep_bkpts_with_len_n(n, k=1, bkpts=None):
    r"""
    Produce representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n.

    Note, this function does not check that the data input is correct and assume it is being used correctly.

    INPUT: n, length of breakpoint sequence, k, length of every element, bkpts, an iterable of breakpoints all of length k.

    OUTPUT: A list of representative elements of every isomorphism class of breakpoints complexes for breakpoint sequences of length n extrapolated from bkpts.
    """
    # Look into using a directed tree as an underlying data structure for generating elements.
    new_bkpts = []
    if n < 2:
        raise ValueError("n>=2")
    if k == n:
        raise ValueError("k<n")
    if k == 1 and bkpts is None:
        bkpts=[[0]]
    for bkpt in bkpts:
        new_bkpts += add_breakpoints_and_find_equiv_classes(nnc_poly_from_bkpt_sequence(bkpt).upstairs())
    new_bkpts = unique_list(new_bkpts)
    k += 1
    if k == n:
        return new_bkpts
    else:
        return make_rep_bkpts_with_len_n(n, k, new_bkpts)


def generate_assumed_symmetric_vertices_continuous(fn, f, bkpt):
    """
    Assumes the symmetry condition holds for all vertices (x,y) in bkpt's breakpoints complex
    such that x+y equiv f.
    """
    for i in range(len(bkpt)):
        x = bkpt[i]
        if x == f:
            continue
        if x < f:
            y = f - x
        else:
            y = 1 + f - x
        fn(x) + fn(y) == 1
        yield (x, y, 0, 0)


def value_nnc_polyhedron_value_cords(bkpt, f_index):
    """
    For a given ``bkpt`` seqeunce and ``f_index``, write the value polyhedron as a BSA in only the value parameters. 
    """
    # this saves a slight amount of overhead when detemrining points for the value polyhedron since the assumed
    # minimality test does not have to entierly go though parametric real field.
    n = len(bkpt)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n - 1)
    if not isinstance(bkpt, list):
        bkpt = list(bkpt)
    coord_names = []
    val = [None]*(n)
    for i in range(n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = val, mutable_values=True, big_cells=True, allow_refinement=False)
    K.gens()[0] == 0
    for i in range(1, n):
        K.gens()[i] <=1
        K.gens()[i] > 0
    h = piecewise_function_from_breakpoints_and_values(bkpt + [1], K.gens() + [0], merge=False)
    # Assumes minimality for the partially defined function.
    for vert in generate_type_1_vertices_continuous(h, operator.ge, bkpt + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, bkpt + [1]):
        vert
    for vert in generate_assumed_symmetric_vertices_continuous(h, bkpt[f_index], bkpt + [1]):
        vert
    return K._bsa
    
# def value_nnc_polyhedron_gamma_0_not_as_param(bkpt, f_index):
    # """
    # For a given ``bkpt`` seqeunce and ``f_index``, find the value polyhedron which assumes pi_(bkpt, v) is minimal.
    # """
    # n = len(bkpt)
    # assert(n >= 2)
    # assert(f_index >= 1)
    # assert(f_index <= n - 1)
    # if not isinstance(bkpt, list):
        # bkpt = list(bkpt)
    # coord_names = []
    # val = [None]*(n-1)
    # for i in range(1, n):
        # coord_names.append('gamma'+str(i))
    # logging.disable(logging.INFO)
    # K = ParametricRealField(names=coord_names, values = val, mutable_values=True, big_cells=True, allow_refinement=False)
    # for i in range(1, n):
        # K.gens()[i] <=1
        # K.gens()[i] > 0
    # h = piecewise_function_from_breakpoints_and_values(bkpt + [1], [0] + K.gens() + [0], merge=False)
    # # Assumes minimality for the partially defined function.
    # for vert in generate_type_1_vertices_continuous(h, operator.ge, bkpt + [1]):
        # vert
    # for vert in generate_type_2_vertices_continuous(h, operator.ge, bkpt + [1]):
        # vert
    # for vert in generate_assumed_symmetric_vertices_continuous(h, bkpt[f_index], bkpt + [1]):
        # vert
    # return K._bsa

def value_nnc_polyhedron(bkpt, f_index):
    """
    For a given ``bkpt`` seqeunce and ``f_index``, write a base which is the value polyhedron corrospoding in the full space of parameters.
    
    EXAMPLES::
    
    """
    n = len(bkpt)
    assert(n >= 2)
    assert(f_index >= 1)
    assert(f_index <= n)
    coord_names = []
    bkpt_vals = list(bkpt)
    vals = bkpt_vals + [None]*(n)
    for i in range(n):
        coord_names.append('lambda'+str(i))
    for i in range(n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = vals, mutable_values=True, big_cells=True, allow_refinement=False)
    # breakpoint parameters are the mesured breakpoint values. 
    for i in range(n):
        K.gens()[i] == bkpt[i]
    # necessary conditions on value parameters
    K.gens()[n] == 0 
    for i in range(1, n):
        K.gens()[i+n] <=1
        K.gens()[i+n] > 0
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n]  + [1], K.gens()[n:2*n] + [0], merge=False)
    # Assumes minimality  for the partially defined function.
    for vert in generate_type_1_vertices_continuous(h, operator.ge, K.gens()[0:n] + [1]):
        vert
    for vert in generate_type_2_vertices_continuous(h, operator.ge, K.gens()[0:n] + [1]):
        vert
    for vert in generate_assumed_symmetric_vertices_continuous(h, K.gens()[f_index-1], [0] + K.gens()[0:n] + [1]):
        vert
    return K._bsa


def bsa_of_rep_element(bkpt, vals):
    """
    Given pi_(bkpt, vals) is {minimal, not minimal}, find BSA subset of R^(2n) such that (bkpt, vals) in BSA and for all p
    in BSA, pi_p is {minimal, not minimal}.

    INPUT: (bkpt, vals) are lists or vectors of length n and bkpt is a proper breakpoints sequence and vals
    is the corresponding value parameters.

    OUTPUT: A basic semialgebraic set.
    
    EXAMPLES::
    
    
    """
    n = len(bkpt)
    assert(n>=2)
    coord_names = []
    for i in range(n):
        coord_names.append('lambda'+str(i))
    for i in range(n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = bkpt+vals, big_cells=True)
    h = piecewise_function_from_breakpoints_and_values(K.gens()[0:n] + [1], K.gens()[n:2*n] + [0], merge=False)
    minimality_test(h)
    return K.make_proof_cell().bsa


def bsa_of_rep_element_pi_of_0_not_param(bkpt, vals):
    """
    Given pi_(bkpt, vals) is {minimal, not minimal}, find BSA subset of R^(2n) such that (bkpt, vals) in BSA and for all p
    in BSA, pi_p is {minimal, not minimal}.

    INPUT: (bkpt, vals) are lists or vectors of length n and bkpt is a proper breakpoints sequence and vals
    is the corresponding value parameters.

    OUTPUT: A basic semialgebraic set.
    
    EXAMPLES::
    
    
    """
    n = len(bkpt)
    if not isinstance(bkpt, list):
        bkpt = list(bkpt)
    if not isinstance(vals, list):
        vals = list(vals)
    assert(n>=2)
    coord_names = []
    for i in range(1,n):
        coord_names.append('lambda'+str(i))
    for i in range(1,n):
        coord_names.append('gamma'+str(i))
    logging.disable(logging.INFO)
    K = ParametricRealField(names=coord_names, values = bkpt[1:]+vals[1:], big_cells=True)
    h = piecewise_function_from_breakpoints_and_values([0]+K.gens()[:n-1] + [1], [0] + K.gens()[n-1:] + [0], merge=False)
    minimality_test(h)
    return K.make_proof_cell().bsa


def find_minimal_function_reps_from_bkpts(bkpts, backend=None):
    """
    Finds representative elements of minimal functions from a given breakpoint sequence. 
    """
    rep_elems = []
    for bkpt in bkpts:
        n = len(bkpt)
        for f_index in range(1, n):
            poly_bsa = value_nnc_polyhedron_value_cords(list(bkpt), f_index)
            gammas = poly_bsa.polynomial_map()[0].parent().gens()
            try:
                test_point = poly_bsa.upstairs().find_point()
            except EmptyBSA:
                raise RepElemGenFailure("The value polyhedron {} is empty. This should not be empty. Double check inputs".format(poly_bsa))
            test_val = []
            for gamma_i in gammas:
                test_val.append(test_point[poly_bsa.v_dict()[gamma_i]])
            h = piecewise_function_from_breakpoints_and_values(list(bkpt)+[1], test_val+[0])
            if not minimality_test(h): # test bkpt doesn't make a valid minimal function, the statement still holds with original bkpt
                raise ValueError("HELP! ({}, {}) paramaterized is not a minimal function but assuming a breakpoint sequence is input, this should be minimal. GL debugging.".format(bkpt, test_val))
            rep_elems.append((bkpt, test_val))
    return rep_elems


class BreakpointComplexClassContainer:
    """
    A container for the family of breakpoint complexes for peicewise linear functions
    with at most n breakpoints.
    
    EXAMPLES::
    
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            if kwrds[backend] == "pplite":
                self._backend = "pplite"
            else:
                self._backend = None
        if "load_rep_elem_data" in kwrds.keys():
            if kwrds[load_rep_elem_data] is None:
                logging.warning("Generating representative elements. This might take a while.")
                self._data = make_rep_bkpts_with_len_n(self._n)
            else:
                file_names = kwrds["load_bkpt_data"].split(",")
                self._data = []
                for file_name in file_names:
                    file = open(file_name, "r")
                    self._data += [eval(preparse(data)) for data in list(csv.reader(file))]
                    file.close()
                if "gen_elems_from_data" in kwrds.keys():
                    if kwrds[gen_elems_from_data] == True:
                        k = len(self._data[0])
                        if k < n:
                            self._data = make_rep_bkpts_with_len_n(n, k, self._data)
        else:
            logging.warning("Generating representative elements. This might take a while.")
            self._data = make_rep_bkpts_with_len_n(self._n)

    def __repr__(self):
        return f"Container for the space breakpoint sequences of length {self._n} under equivlance of polyhedral complexes."

    def get_rep_elems(self):
        for bkpt in self._data:
            yield bkpt

    def get_nnc_poly_from_bkpt(self):
        for bkpt in self._data:
            for f_index in range(1,n):
                yield nnc_poly_from_bkpt(bkpt, f_index)

    def num_rep_elems(self):
        return len(self._data)

    def add_one_bkpt_to_all(self):
        logging.warning("Generating representative elements. This might take a while.")
        self._data = make_bkpts_with_len_n(self._n+1, self._n, self._data)

    def write_data(self, output_file_name_style=None, max_rows=None):
        """
        Writes representative element data to a `.csv` file with one column and rows of representative elements.
        Optionally, write many `.csv` files with at `most max_rows` rows per file.
        Files are named output_file_file_name_style_filenumber.csv.
        The default output_file_name_style="bkpts_of_len_n".
        """
        # TODO: Future, support writing different types of data such as polyhedra data.
        if output_file_name_style is None:
            file_name_base = "bkpts_of_len_{}".format(self._n)
        else:
            file_name_base =output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = 0
        output_file = file_name_base +".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_rows):
                try:
                    data_writer.writerow(self._data[max_rows * file_number + row])
                except IndexError:
                    break
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number+1)+".csv"


class PiMinContContainer:
    """
    A container for the space of continuous piecewise linear minimal functions with at
    most n breakpoints paramaterized by breakpoints and values using semialgebraic sets.

    TESTS::

    >>> PiMin_at_most_4_breakpoints = PiMinContContainer(4)
    >>> all([minimality_test(pi) for PiMin_at_most_4_breakpoints.get_rep_functions()])
    True
    >>>
    """
    def __init__(self, n, **kwrds):
        self._n = n
        assert(self._n >= 2)
        if "backend" in kwrds.keys():
            if kwrds[backend] == "pplite":
                self._backend = "pplite"
            else:
                self._backend = None
        if "load_bkpt_data" in kwrds.keys() and "load_rep_elem_data" not in kwrds.keys():
            file_names = kwrds["load_bkpt_data"].split(",")
            bkpts = []
            for file_name in file_names:
                with open(file_name, newline='' ) as csvfile:
                    file_reader = csv.reader(csvfile)
                    for row in file_reader:
                        bkpts.append([eval(preparse(data)) for data in row])
            self._data = find_minimal_function_reps_from_bkpts(bkpts)
        elif "load_bkpt_data" not in kwrds.keys() and "load_rep_elem_data" in kwrds.keys():
            file_names = kwrds["load_rep_elem_data"].strip(" ").split(",")
            self._data = []
            for file_name in file_names:
                with open(file_name, newline='' ) as csvfile:
                    file_reader = csv.reader(csvfile)
                    for row in file_reader:
                        self._data.append([eval(preparse(data)) for data in row])
        else:
            logging.warning("Generating representative elements. This might take a while.")
            bkpts = make_rep_bkpts_with_len_n(self._n)
            self._data = find_minimal_function_reps_from_bkpts(bkpts)

    def __repr__(self):
        return "Space of minimal functions with at most {} breakpoints parameterized by breakpoints and values using semialgebraic sets.".format(self._n)

    def get_semialgebraic_sets(self):
        for b, v in self._data:
            yield bsa_of_rep_element(list(b), list(v))

    def get_rep_elems(self):
        for b, v in self._data:
            yield (list(b), list(v))

    def get_rep_functions(self):
        for b, v in self._data:
            yield piecewise_function_from_breakpoints_and_values(list(b)+[1], list(v)+[0])

    def n(self):
        return self._n

    def covers_space(self):
        raise NotImplementedError

    def refine_space(self):
        raise NotImplementedError

    def write_data(self, output_file_name_style=None, max_rows=None):
        """
        Writes representative element data to a `.csv` file with one column and rows of representative elements.
        Optionally, write many `.csv` files with at `most max_rows` rows per file.
        Files are named output_file_file_name_style_filenumber.csv.
        The default output_file_name_style="Pi_Min_n".
        """
        # TODO: Future, support writing different types of data such as polyhedra data.
        if output_file_name_style is None:
            file_name_base = "Pi_Min_{}".format(self._n)
        else:
            file_name_base =output_file_name_style
        if max_rows is not None:
            assert(max_rows >= 1)
            num_files = len(self._data)//max_rows + 1
            file_name_base = file_name_base + "_part_0"
        if max_rows is None:
            max_rows = len(self._data)
            num_files = 1
        output_file = file_name_base +".csv"
        for file_number in range(num_files):
            out_file = open(output_file, "w")
            data_writer = csv.writer(out_file, csv.QUOTE_NONE)
            for row in range(max_rows):
                try:
                    data_writer.writerow(self._data[max_rows * file_number + row])
                except IndexError:
                    break
            out_file.close()
            output_file = file_name_base[:-1]+"{}".format(file_number+1)+".csv"


### Plotting Utilties ###