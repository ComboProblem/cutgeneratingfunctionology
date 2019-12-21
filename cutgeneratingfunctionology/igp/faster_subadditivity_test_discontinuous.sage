from itertools import chain

import queue as queue
import itertools
import numpy as np

class SubadditivityTestTreeNodeGeneral(object):


    def __init__(self, fn, level, intervals):

        self.intervals=tuple(tuple(I) for I in intervals)
        self.level=level
        self.function=fn
        self.vertices=verts(*intervals)
        self.projections=projections(self.vertices)
        self.left_child=None
        self.right_child=None
        self.parent=None
        self.affine_estimators=None

    def __lt__(self, other):
        return self.intervals < other.intervals
    def __eq__(self, other):
        return self.intervals == other.intervals
    def __le__(self, other):
        return self.intervals <= other.intervals
    def __ge__(self, other):
        return self.intervals >= other.intervals
    def __gt__(self, other):
        return self.intervals > other.intervals
    def __ne__(self, other):
        return self.intervals != other.intervals
    def __hash__(self):
        return hash(self.intervals)

    def I_bkpts(self):
        if hasattr(self,'_I_bkpts'):
            return self._I_bkpts
        new_I=self.projections[0]
        self._I_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_I))
        return self._I_bkpts

    def J_bkpts(self):
        if hasattr(self,'_J_bkpts'):
            return self._J_bkpts
        new_J=self.projections[1]
        self._J_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_J))
        return self._J_bkpts

    def K_bkpts(self):
        if hasattr(self,'_K_bkpts'):
            return self._K_bkpts
        new_K=self.projections[2]
        self._K_bkpts = tuple(find_all_bkpts_in_the_interval(self.function.end_points(),new_K))
        return self._K_bkpts

    def I_values(self):
        if hasattr(self,'_I_values'):
            return self._I_values
        self._I_values = {}
        for v in self.I_bkpts():
            self._I_values[v] = copy(self.function.limits(v))
        self._I_values[self.I_bkpts()[0]].pop(-1)
        self._I_values[self.I_bkpts()[-1]].pop(1)
        return self._I_values

    def J_values(self):
        if hasattr(self,'_J_values'):
            return self._J_values
        self._J_values = {}
        for v in self.J_bkpts():
            self._J_values[v] = copy(self.function.limits(v))
        self._J_values[self.J_bkpts()[0]].pop(-1)
        self._J_values[self.J_bkpts()[-1]].pop(1)
        return self._J_values

    def K_values(self):
        if hasattr(self,'_K_values'):
            return self._K_values
        self._K_values = {}
        for v in self.K_bkpts():
            self._K_values[v] = copy(self.function.limits(fractional(v)))
        self._K_values[self.K_bkpts()[0]].pop(-1)
        self._K_values[self.K_bkpts()[-1]].pop(1)
        return self._K_values

    def delta_pi_lower_bound(self):
        alpha_I=min(itertools.chain(*self.I_values().values()))
        alpha_J=min(itertools.chain(*self.J_values().values()))
        beta_K=max(itertools.chain(*self.K_values().values()))
        return alpha_I+alpha_J-beta_K

    def delta_pi_upper_bound(self):
        if hasattr(self,'_delta_pi_ub'):
            return self._delta_pi_ub
        self._delta_pi_ub=min(self.function(v[0])+self.function(v[1])-self.function(fractional(v[0]+v[1])) for v in self.vertices)
        return self._delta_pi_ub

    def delta_pi_min_of_indivisible_node(self):
        # zero dimensional face
        res=min(self.function(v[0])+self.function(v[1])-self.function(fractional(v[0]+v[1])) for v in self.vertices)

        # one dimensional face
        one_dimensional_faces=[]
        for i in range(len(self.vertices)):
            for j in range(i+1,len(self.vertices)):
                v1=self.vertices[i]
                v2=self.vertices[j]
                if v1[0] == v2[0] or v1[1] == v2[1] or v1[0]+v1[1] == v2[0]+v2[1]:
                    one_dimensional_faces.append((v1,v2))
        for edge in one_dimensional_faces:
            mid_point = find_midpoint(edge[0],edge[1])
            quarter_point_1 = find_midpoint(edge[0],mid_point)
            quarter_point_2 = find_midpoint(edge[1],mid_point)
            delta_pi_mid_point = delta_pi(self.function, mid_point[0], mid_point[1])
            delta_pi_quarter_point_1 = delta_pi(self.function, quarter_point_1[0], quarter_point_1[1])
            delta_pi_quarter_point_2 = delta_pi(self.function, quarter_point_2[0], quarter_point_2[1])
            delta_pi_left_end_point = 2*delta_pi_quarter_point_1 - delta_pi_mid_point
            delta_pi_right_end_point = 2*delta_pi_quarter_point_2 - delta_pi_mid_point
            res=min(res,delta_pi_left_end_point,delta_pi_right_end_point)

        # two dimensional face
        center = (sum(v[0] for v in self.vertices)/len(self.vertices),(sum(v[1] for v in self.vertices)/len(self.vertices)))
        delta_pi_center = delta_pi(self.function, center[0], center[1])
        for v in self.vertices:
            mid_point = find_midpoint(v,center)
            delta_pi_mid_point = delta_pi(self.function, mid_point[0], mid_point[1])
            delta_pi_v = 2*delta_pi_mid_point-delta_pi_center
            res=min(res,delta_pi_v)

        self._delta_pi_ub = res
        return self._delta_pi_ub

    def branching_direction(self):
        new_I,new_J,new_K=self.projections
        lenI=new_I[1]-new_I[0]
        lenJ=new_J[1]-new_J[0]
        lenK=new_K[1]-new_K[0]
        candidates=[]
        if len(self.I_bkpts())>2:
            candidates.append(('I',lenI))
        if len(self.J_bkpts())>2:
            candidates.append(('J',lenJ))
        if len(self.K_bkpts())>2:
            candidates.append(('K',lenK))
        if not candidates:
            return None
        else:
            ind=np.argmax([candidate[1] for candidate in candidates])
            return candidates[ind][0]

    def new_intervals(self):
        dir=self.branching_direction()
        if dir=='I':
            new_bkpt=self.I_bkpts()[floor(len(self.I_bkpts())/2)]
            return (((self.intervals[0][0],new_bkpt),self.intervals[1],self.intervals[2]),
                    ((new_bkpt,self.intervals[0][1]),self.intervals[1],self.intervals[2]))
        elif dir=='J':
            new_bkpt=self.J_bkpts()[floor(len(self.J_bkpts())/2)]
            return ((self.intervals[0],(self.intervals[1][0],new_bkpt),self.intervals[2]),
                    (self.intervals[0],(new_bkpt,self.intervals[1][1]),self.intervals[2]))
        elif dir=='K':
            new_bkpt=self.K_bkpts()[floor(len(self.K_bkpts())/2)]
            return ((self.intervals[0],self.intervals[1],(self.intervals[2][0],new_bkpt)),
                    (self.intervals[0],self.intervals[1],(new_bkpt,self.intervals[2][1])))
        else:
            raise ValueError("Indivisible Region.")

    def is_divisible(self):
        if len(self.I_bkpts())<=2 and len(self.J_bkpts())<=2 and len(self.K_bkpts())<=2:
            return False
        else:
            return True

    def generate_children(self,upper_bound=0,stop_only_if_strict=True):
        if not self.is_fathomed(upper_bound,stop_only_if_strict):
            I1,I2=self.new_intervals()
            self.left_child=SubadditivityTestTreeNodeGeneral(self.function,self.level+1,I1)
            self.right_child=SubadditivityTestTreeNodeGeneral(self.function,self.level+1,I2)
            self.left_child.parent=self
            self.right_child.parent=self

    def is_fathomed(self,value=0,stop_only_if_strict=True):
        if self.is_divisible():
            if stop_only_if_strict:
                if self.delta_pi_lower_bound()>value:
                    return True
                else:
                    return False
            else:
                # stop branching early.
                if self.delta_pi_lower_bound()>=value:
                    return True
                else:
                    return False
        else:
            return True

class SubadditivityTestTreeGeneral:

    r"""
    Class for spatial branch and bound for subadditivity testing handling discontinuous cases.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)
        sage: h = drlm_2_slope_limit()
        sage: T = SubadditivityTestTree(h)
        sage: T.minimum()
        0
    """

    def __init__(self,fn,intervals=((0,1), (0,1), (0,2)),global_upper_bound=0,objective_limit=0):
        self.function=fn
        self.intervals=intervals
        self.global_upper_bound=global_upper_bound
        self.objective_limit=objective_limit
        self.root=SubadditivityTestTreeNodeGeneral(fn,0,intervals)
        self.height=0
        self.complete_node_set=set([self.root])
        self.leaf_set=set([self.root])
        # the order of unfathomed nodes matters, like DFS or BFS
        self.unfathomed_node_list=queue.PriorityQueue()
        self.nonsubadditive_vertices=set()
        self.additive_vertices=set()
        self.maximal_additive_faces=set()

    def number_of_nodes(self):
        return len(self.complete_node_set)

    def number_of_leaves(self):
        return len(self.leaf_set)

    def node_branching(self,node,search_method='BB',find_min=True,stop_only_if_strict=True):
        if not node.left_child:
            if find_min:
                node.generate_children(self.global_upper_bound,stop_only_if_strict)
            else:
                node.generate_children(self.objective_limit,stop_only_if_strict)
            if node.left_child:
                if search_method=='BFS':
                    self.unfathomed_node_list.put((len(self.complete_node_set),node.left_child))
                    self.unfathomed_node_list.put((len(self.complete_node_set)+1,node.right_child))
                elif search_method=='DFS':
                    self.unfathomed_node_list.put((-len(self.complete_node_set),node.left_child))
                    self.unfathomed_node_list.put((-len(self.complete_node_set)-1,node.right_child))
                elif search_method=='BB':
                    if node.left_child.is_divisible():
                        self.unfathomed_node_list.put((node.left_child.delta_pi_lower_bound(),node.left_child))
                    else:
                        self.unfathomed_node_list.put((node.left_child.delta_pi_upper_bound(),node.left_child))
                    if node.right_child.is_divisible():
                        self.unfathomed_node_list.put((node.right_child.delta_pi_lower_bound(),node.right_child))
                    else:
                        self.unfathomed_node_list.put((node.right_child.delta_pi_upper_bound(),node.right_child))
                else:
                    raise ValueError("Can't recognize search_method.")
                self.height=max(self.height,node.left_child.level)
                self.complete_node_set.update({node.left_child,node.right_child})
                self.leaf_set.discard(node)
                self.leaf_set.update({node.left_child,node.right_child})

    def minimum(self,search_method='BB'):
        if search_method=='BFS' or search_method=='DFS':
            self.unfathomed_node_list.put((0,self.root))
        elif search_method=='BB':
            self.unfathomed_node_list.put((self.root.delta_pi_lower_bound(),self.root))
        else:
            raise ValueError("Can't recognize search_method.")
        while not self.unfathomed_node_list.empty():
            current_node=self.unfathomed_node_list.get()[1]
            if len(current_node.vertices)==0:
                    continue
            upper_bound=current_node.delta_pi_upper_bound()
            if upper_bound<self.global_upper_bound:
                self.global_upper_bound=upper_bound
            self.node_branching(current_node,search_method=search_method,find_min=True,stop_only_if_strict=False)
            if not current_node.is_divisible():
                upper_bound=current_node.delta_pi_min_of_indivisible_node()
                if upper_bound<self.global_upper_bound:
                    self.global_upper_bound=upper_bound
        self.min=self.global_upper_bound
        return self.min

def find_midpoint(v1,v2):
    return ((v1[0]+v2[0])/2, (v1[1]+v2[1])/2)
