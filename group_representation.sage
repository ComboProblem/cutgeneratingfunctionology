# Make sure current directory is in path.
# That's not true while doctesting (sage -t).
if '' not in sys.path:
    sys.path = [''] + sys.path

from igp import *

def cpl3_function(r0, z1, o1, o2):
    """
    Construct a CPL3 function.

    Parameters:
        0 < r0 (real) < 1, 0 < z1 (real), 0 <= o1, o2 (real)

    EXAMPLE::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: p = cpl3_function(r0=1/7, z1=1/7, o1=1/4, o2=1/12)
        sage: q = plot(p)
    """
    if not (bool(r0 > 0) & bool(z1 > 0) & bool(o1 > 0) & bool(o2 > 0)):
        raise ValueError, "Bad parameters. Unable to construct the function."
    if not bool(((1-r0)/2-2*z1 > 0) and bool(1/2-o1-o2 > 0)):
        raise ValueError, "Conditions for a CPL-3 function are NOT satisfied."

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    values = [0, 0, o1, o1+o2, 1-(o1+o2), 1-o1, 1]

    list_of_pairs = [[(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], values[i]], [bkpt[i+1], values[i+1]])] for i in range(len(bkpt)-1)]

    return PiecewiseQuasiPeriodic(list_of_pairs)

def standard_representation_from_group_representation(fn):
    """
    Convert a standard representation 'phi' (a superadditive quasiperiodic function) from a group representation 'fn' (a subadditive periodic function).

    EXAMPLE::
        
        sage: fn = mlr_cpl3_d_3_slope(r0=1/7, z1=1/7)
        sage: phi = standard_representation_from_group_representation(fn)
        sage: q = plot(phi)
    """
    f = find_f(fn)

    bkpt = fn._end_points
    phi_at_bkpt = [bkpt[i]-fn(bkpt[i])*f for i in range(len(bkpt))]

    list_of_pairs = []
    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], phi_at_bkpt[i]], [bkpt[i+1], phi_at_bkpt[i+1]])])

    phi = PiecewiseQuasiPeriodic(list_of_pairs)

    return phi

def group_representation_from_standard_representation(phi):
    """
    Convert a group representation 'fn' (a subadditive periodic function) from a standard representation 'phi' (a superadditive quasiperiodic function).

    EXAMPLE::

        sage: phi = cpl3_function(r0=1/7, z1=1/7, o1=1/4, o2=1/12)
        sage: fn = group_representation_from_standard_representation(phi)
        sage: q = plot(fn)
    """
    bkpt = phi._end_points
    value = phi.periodic_term._values_at_end_points

    f = -min([phi(bkpt[i])-bkpt[i] for i in range(len(bkpt))])
    fn_at_bkpt = [(bkpt[i]-phi(bkpt[i]))/f for i in range(len(bkpt))]

    list_of_pairs = []
    for i in range(len(bkpt)-1):
        list_of_pairs.append([(bkpt[i], bkpt[i+1]), linear_function_through_points([bkpt[i], fn_at_bkpt[i]], [bkpt[i+1], fn_at_bkpt[i+1]])])

    fn = FastPiecewise(list_of_pairs)

    return fn

# Extreme functions (extreme_functions_in_literature.sage)
def mlr_cpl3_a_2_slope(r0=3/13, z1=3/26, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt a.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        0 < r0 < 1, 0 <= z1 < 1

    Examples:
        page 183, Fig 2, point a::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_a_2_slope(r0=3/13, z1=3/26)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 < 1) & bool(0 <= z1 < 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = r0

    bkpt = [0, r0, 1]
    slopes = [1/r0, 1/(r0-1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_b_3_slope(r0=3/26, z1=1/13, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt b.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        3*r0 + 8*z1 <= 1

    Note:
        mlr_cpl3_b_3_slope(r0,z1) is the same as drlm_backward_3_slope(f=r0,bkpt=r0+2*z1).

    Examples:
        page 183, Fig 2, point b::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h1 = mlr_cpl3_b_3_slope(r0=3/26, z1=1/13)
        sage: extremality_test(h1)
        True
        sage: h2 = drlm_backward_3_slope(f=3/26,bkpt=7/26)
        sage: extremality_test(h2)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not bool(3*r0 + 8*z1 <= 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-3*r0)/8

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (2*z1-1)/(2*z1*(1+r0)), (2*z1-1)/(2*z1*(1+r0)), 1/(1+r0), (2*z1-1)/(2*z1*(1+r0)), (2*z1-1)/(2*z1*(1+r0))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_c_3_slope(r0=5/24, z1=1/12, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt c.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        3*r0 + 4*z1 <= 1

    Note:
        mlr_cpl3_c_3_slope(r0,z1) is the same as drlm_backward_3_slope(f=r0,bkpt=r0+z1).

    Examples:
        page 183, Fig 2, point c::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h1 = mlr_cpl3_c_3_slope(r0=5/24, z1=1/12)
        sage: extremality_test(h1)
        True
        sage: h2 = drlm_backward_3_slope(f=5/24,bkpt=7/24)
        sage: extremality_test(h2)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not bool(3*r0 + 4*z1 <= 1):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-3*r0)/4

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (z1-1)/(z1*(1+r0)), 1/(1+r0), 1/(1+r0), 1/(1+r0), (z1-1)/(z1*(1+r0))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_d_3_slope(r0=1/6, z1=1/12, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt d.
        - Proven extreme p.188, thm.19.
       
    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 = 2*z1, r0 + 8*z1 <= 1

    Examples:
        p.183, Fig 2, point d1::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_d_3_slope(r0=1/6, z1=1/12)
        sage: extremality_test(h, False)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 == 2*z1) & bool(r0+8*z1 <= 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = r0/2

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (2*z1+1)/(2*z1*(r0-1)), (1-2*z1)/(2*z1*(1-r0)), 1/(r0-1), (1-2*z1)/(2*z1*(1-r0)), (2*z1+1)/(2*z1*(r0-1))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_f_2_or_3_slope(r0=1/6, z1=1/6, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 or 3; Continuous.
        - Discovered p.179, Table 3, Ext. pnt f.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 <= z1, r0 + 5*z1 = 1

    Examples:
        page 184, Fig 3, point f1 and f2::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_f_2_or_3_slope(r0=1/6, z1=1/6)
        sage: extremality_test(h)
        True
        sage: h = mlr_cpl3_f_2_or_3_slope(r0=3/23, z1=4/23)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 <= z1) & bool(r0+5*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/5

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1)), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-1)/(r0^2+4*z1+8*z1*r0-1), (r0+8*z1-2)/(r0^2+4*z1+8*z1*r0-1), (z1*r0+8*z1^2-r0-6*z1+1)/(z1*(r0^2+4*z1+8*z1*r0-1))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_g_3_slope(r0=1/12, z1=5/24, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt g.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 < z1, 2*r0 + 4*z1 = 1

    Examples:
        page 184, Fig 3, point g::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_g_3_slope(r0=1/12, z1=5/24)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 < z1) & bool(2*r0 + 4*z1 == 1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-2*r0)/4

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1), (1-3*z1)/(z1*(1-3*r0)), 3/(3*r0-1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_h_2_slope(r0=1/4, z1=1/6, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt h.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 + 4*z1 <= 1 < 2*r0 + 4*z1

    Examples:
        page 183, Fig 2, point h::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_h_2_slope(r0=1/4, z1=1/6)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 + 4*z1 <= 1 < 2*r0 + 4*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/4

    bkpt = [0, r0, r0+2*z1, 1-2*z1, 1]
    slopes = [1/r0, (4*z1-1)/(4*r0*z1), 1/r0, (4*z1-1)/(4*r0*z1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_k_2_slope(r0=7/27, z1=4/27, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt k.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 <= 2*z1, r0 + 5*z1 = 1

    Examples:
        page 185, Fig 4, point k1::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_k_2_slope(r0=7/27, z1=4/27)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 <= 2*z1) & bool(r0+5*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/5

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (3*z1-1)/(3*z1*r0), 1/r0, (2-3*r0-12*z1)/(3*r0*(1-r0-4*z1)), 1/r0, (3*z1-1)/(3*z1*r0)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_n_3_slope(r0=16/22, z1=1/22, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 3 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt n.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 > 2*z1, r0 + 8*z1 < 1

    Examples:
        page 185, Fig 4, point n2::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_n_3_slope(r0=9/25, z1=2/25)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 > 2*z1) & bool(r0+8*z1<1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/8

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (1+r0)/(r0*(r0-1)), 1/r0, 1/(r0-1), 1/r0, (1+r0)/(r0*(r0-1))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_o_2_slope(r0=3/8, z1=1/8, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt o.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 >= 2*z1, 2*r0 + 2*z1 = 1

    Examples:
        page 186, Fig 5, point o::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_o_2_slope(r0=3/8, z1=1/8)
        sage: extremality_test(h)
        True

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 >= 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-2*r0)/2

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, 2/(2*r0-1), (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 1/r0, (4*z1-4*z1*r0-1+2*r0)/(2*r0*z1*(1-2*r0)), 2/(2*r0-1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_q_2_slope(r0=5/12, z1=3/24, field=None, conditioncheck=True):
    """

    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt q.
        - Proven extreme p.188, thm.19.


    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:

        r0 > 2*z1, r0+4*z1 <= 1 < r0+5*z1

    Examples:
        page 186, Fig 5, point q::


        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_q_2_slope(r0=5/12, z1=3/24)
        sage: extremality_test(h)
        True


    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 > 2*z1) & bool(r0+4*z1 <= 1 < r0+5*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/4

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0,(r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1)), 1/r0, (r0+2*z1)/(r0*(-1+r0+2*z1))]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_r_2_slope(r0=3/7, z1=1/7, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt r.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 > 2*z1, r0+4*z1 <= 1 <= 2*r0+2*z1

    Examples:
        page 185, Fig , point r::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_r_2_slope(r0=3/7, z1=1/7)
        sage: extremality_test(h)
        True 

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 > 2*z1) & bool(r0+4*z1<=1<=2*r0+2*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-r0)/4

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, (2*z1-1)/(2*r0*z1), 1/r0, 1/r0, 1/r0, (2*z1-1)/(2*r0*z1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

# The functions 'l' and 'p' that need to be checked.
def mlr_cpl3_l_2_slope(r0=8/25, z1=4/25, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt l.

        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:

        r0 = 2*z1, 6*z1 <= 1 < 7*z1

    Examples:
        page 185, Fig 4, point l::


        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_l_2_slope(r0=8/25, z1=4/25)
        sage: extremality_test(h)
        True 


    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 == 2*z1) & bool(6*z1<=1<7*z1)):
            logging.info("Conditions for extremality are NOT satisfied.")
   	else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = r0/2

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, 2/(2*r0-1), -4*z1/(2*z1*(1-2*r0)), 2/(2*r0-1), -4*z1/(2*z1*(1-2*r0)), 2/(2*r0-1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)

def mlr_cpl3_p_2_slope(r0=5/12, z1=1/12, field=None, conditioncheck=True):
    """
    Summary:
        - The group representation of the continuous piecewise linear lifting (CPL) function.
        - Inifinite; Dim = 1; Slopes = 2 ; Continuous.
        - Discovered p.179, Table 3, Ext. pnt p.
        - Proven extreme p.188, thm.19.

    Parameters:
        0< r0 (real) < 1, 0 < z1 (real)
    Function is known to be extreme under the conditions:
        r0 > 2*z1, 2*r0 + 2*z1 = 1

    Examples:
        page 186, Fig 5, point p1 and p2::

        sage: logging.disable(logging.WARN) # Suppress warning about experimental discontinuous code.
        sage: h = mlr_cpl3_p_2_slope(r0=5/12, z1=1/12)
        sage: extremality_test(h)
        True # Fix parameters

    Check: r0=5/12, z1=1/12 (satisfiy the condition for extremality give in the paper but the function is different given in the paper) 
           h = piecewise_function_from_breakpoints_and_slopes([0,5/12,6/12,7/12,10/12,11/12,1],[12/5,-12,12/5,-68/5,12/5,-12])

    Reference:
        L. A. Miller, Y. Li, and J.-P. P. Richard, New Inequalities for Finite and Inifinite Group Problems from Approximate Lifting, Naval Research Logistics 55 (2008), no.2, 172-191, doi:10.1002/nav.20275
    """
    if conditioncheck:
        if not bool(r0 > 0):
            raise ValueError, "Bad parameters. Unable to construct the function."
        if not (bool(r0 > 2*z1) & bool(2*r0+2*z1==1)):
            logging.info("Conditions for extremality are NOT satisfied.")
        else:
            logging.info("Conditions for extremality are satisfied.")

    if z1 is None:
        z1 = (1-2*r0)/2

    bkpt = [0, r0, r0+z1, r0+2*z1, 1-2*z1, 1-z1, 1]
    slopes = [1/r0, 2/(2*r0-1), 1/r0, (2*z1-10*z1*r0+r0)/(r0*(1-2*r0)*(4*z1-1+r0)), 1/r0, 2/(2*r0-1)]

    return piecewise_function_from_breakpoints_and_slopes(bkpt, slopes, field=None)





