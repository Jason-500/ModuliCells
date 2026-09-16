from sage.all import *
import copy as pycopy
from libs.Functions import *
from libs.CustomErrors import *

from itertools import combinations

class DFV:
    """
    class DFV: for managing Delaunay with faces and dual(Voronoi).

    D: Delaunay graph;

    F: Faces of Delaunay;
    
    V: Voronoi dual of D.
    """
    vertex_id = {i: i for i in range(1, 7)}
    edge_id = {f"x{i}": f"x{i}" for i in range(1, 13)}
    id = [vertex_id,edge_id]
    
    def __init__(self,D,F,V,is_interior=False,is_canonical=False,is_root_cell=False,root_cell=None,poly=None):
        """
        class DFV: for managing Delaunay with faces and dual(Voronoi).
        D: Delaunay graph;
        F: Faces of Delaunay;
        V: Voronoi dual of D.
        """
        self.D = D
        self.F = F
        self.V = V
        self.card_E = len(D.edges())
        self.card_F = len(F)
        self._poly = poly
        self._canonical_poly = None
        self._dim = None
        self._faces = {}
        self._canonical_faces = {}
        self._edge_labels = [label for _,_,label in D.edges()]
        self._is_root_cell = is_root_cell
        self._root_cell = root_cell

        self.is_interior = None
        if is_interior: 
            self.is_interior = True

        self.is_canonical = None
        self.canonical_image = None
        if is_canonical: 
            self.is_canonical = True
            self.canonical_image = (self,DFV.id)
        
        self._face_dict = {}
        self._all_faces_loaded = False

        self._planar_graph_automorphism_group = None

    # def remove_single_edge(self,i):
    #     """
    #     Return the subgraphs obtained by removing ith edge from given Delaunay graph D, with modification of F and V.
    #     """
    #     d = self.D.copy()
    #     ei = d.edges()[i]
    #     d.delete_edge(ei)
    #     f,v = modified_faces_with_voronoi(ei,self.F,self.V)
    #     return DFV(d,f,v)

    def remove_multiple_edges(self,indices, is_interior=True, poly = None):
        d, f, v = [pycopy.deepcopy(x) for x in (self.D, self.F, self.V)]
        e_list = self.D.edges()
        for i in indices:
            d.delete_edge(e_list[i])
            f,v = modified_faces_with_voronoi(e_list[i],f,v)
        return DFV(d,f,v,is_interior, root_cell = self, poly = poly)
    
    @property
    def is_root_cell(self):
        return self._is_root_cell
    
    @property
    def root_cell(self):
        if not self.is_root_cell:
            return self._root_cell
        else:
            return self

    @property
    def poly(self):
        if self._poly: return self._poly
        self._poly = self.canonical_poly
        return self._poly
    
    @property
    def dim(self):
        return self.poly.dim()

    @property
    def edge_labels(self):
        return self._edge_labels
    
    @property
    def canonical_poly(self):
        if self._canonical_poly: return self._canonical_poly
        self._canonical_poly = poly(self.D,self.F)
        self._dim = self._canonical_poly.dim()
        return self._canonical_poly
    
    @property
    def face_dict(self): # 暂时face_dict只做了interior face.
        if not self.is_root_cell:
            return self.root_cell.face_dict
        else:
            if self._all_faces_loaded: return self._face_dict
            self._init_faces()
            self._all_faces_loaded = True
            return self._face_dict
        
    def _init_faces(self, dim=None, solve_boundary_cells = False):
        if dim is None:
            for d in range(self.dim + 1):
                self._init_faces(d, solve_boundary_cells)
            return

        if dim in self._faces: 
            return

        raw_interior_faces = []
        if solve_boundary_cells:
            raw_boundary_faces = []

        voronoi_edge_dict = {label:(start,end) for start,end,label in self.V.edges()}
        for face in self.canonical_poly.faces(dim):
            poly = face.as_polyhedron()
            center = poly.center()
            index_of_collapsed_edges = [i for i, x in enumerate(center) if x == 1]
            labels_of_collapsed_edges = [self.D.edges()[i][2] for i in index_of_collapsed_edges]
            labels_of_remained_edges = frozenset([self.D.edges()[i][2] for i in [i for i, x in enumerate(center) if x != 1]])

            if 0 in center: # non-separating degeneration
                if solve_boundary_cells:
                    raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False,poly=poly))

            elif not Graph([voronoi_edge_dict[label] for label in labels_of_collapsed_edges],multiedges=True, loops=True).is_forest(): 
                # separating degeneration. I think it will always occur by pinching a single selfloop in V. But first we still detect the loop.
                if solve_boundary_cells:
                    raw_boundary_faces.append(self.remove_multiple_edges(index_of_collapsed_edges,is_interior=False,poly=poly))

            else:
                if self.is_root_cell: 
                    # 以下仅在生成face_dict时使用.
                    generated_face = self.remove_multiple_edges(index_of_collapsed_edges,is_interior=True,poly=poly)
                    raw_interior_faces.append(generated_face)
                    self._face_dict[labels_of_remained_edges] = generated_face
                else: 
                    raw_interior_faces.append(self.face_dict[labels_of_remained_edges])
        
        self._faces[dim]=[raw_interior_faces,[]]
        if solve_boundary_cells:
            self._faces[dim][1] = raw_boundary_faces

    def faces(self, dim=None, solve_boundary_cells = False):
        """
        Return the faces with given dimension: [interior faces, boundary faces]
        """

        if dim is None:
            self._init_faces(solve_boundary_cells = solve_boundary_cells)
            return {
                d: self.faces(d)
                for d in range(self.dim + 1)
            }
        
        if dim in self._faces: 
            return self._faces[dim]

        self._init_faces(dim, solve_boundary_cells)
        return self._faces[dim]
    
    def facets(self, solve_boundary_cells = False):
        return self.faces(self.dim - 1, solve_boundary_cells)
    
    def init_canonical_faces(self, dim = None, solve_boundary_cells = False):
        if dim is None:
            for d in range(self.dim + 1):
                self.init_canonical_faces(d, solve_boundary_cells)
            return

        if not self.is_root_cell:
            self.root_cell.init_canonical_faces(dim, solve_boundary_cells) # 顶部先确定范本, 随后为所有face指派其等价类.
            new_interior_cell_list = []

            interior_faces, boundary_faces = self.faces(dim, solve_boundary_cells)

            for interior_face in interior_faces:
                if interior_face.canonical_image[0] not in new_interior_cell_list: # 按理说这里应该已经获取到了canonical image, 因为faces的策略就是为非顶部cell返回顶部cell的对应胞腔.
                    new_interior_cell_list.append(interior_face.canonical_image[0])

            self._canonical_faces[dim] = [new_interior_cell_list,[]]

            if solve_boundary_cells:
                new_boundary_cell_list = []
                for boundary_face in boundary_faces:
                    if boundary_face.canonical_image[0] not in new_boundary_cell_list:
                        new_boundary_cell_list.append(boundary_face.canonical_image[0])

                self._canonical_faces[dim][1] = new_boundary_cell_list
            return
        
        if dim in self._canonical_faces: return

        new_interior_cell_list = []# if self.is_root_cell else self.root_cell.canonical_faces(dim)[0]

        if solve_boundary_cells:
            new_boundary_cell_list = []# if self.is_root_cell else self.root_cell.canonical_faces(dim)[1]

        interior_faces, boundary_faces = self.faces(dim, solve_boundary_cells)

        for interior_face in interior_faces:
            is_canonical = True
            for interior_cell in new_interior_cell_list:
                is_isomorphic, permutation = interior_cell.is_orientation_preserving_isomorphic_to(interior_face)
                if is_isomorphic:
                    interior_face.is_canonical = False
                    interior_face.canonical_image = (interior_cell,permutation) # may need to fix the direction of map here.
                    is_canonical = False
                    break
            if is_canonical:
                if not self.is_root_cell:
                    raise MissingCellError()
                interior_face.set_canonical()
                new_interior_cell_list.append(interior_face)

        if solve_boundary_cells:
            for boundary_face in boundary_faces:
                is_canonical = True
                for boundary_cell in new_boundary_cell_list:
                    is_isomorphic, permutation = boundary_cell.is_orientation_preserving_isomorphic_to(boundary_face)
                    if is_isomorphic:
                        boundary_face.is_canonical = False
                        boundary_face.canonical_image = (boundary_cell,permutation) # may need to fix the direction of map here.
                        is_canonical = False
                        break
                if is_canonical:
                    if not self.is_root_cell:
                        raise MissingCellError()
                    boundary_face.set_canonical()
                    new_boundary_cell_list.append(boundary_face)

        self._canonical_faces[dim] = [new_interior_cell_list,[]]
        if solve_boundary_cells:
            self._canonical_faces[dim][1] = new_boundary_cell_list

    def canonical_faces(self, dim = None, solve_boundary_cells = False):

        if dim is None:
            self._init_faces(solve_boundary_cells = solve_boundary_cells)
            return {
                d: self.canonical_faces(d)
                for d in range(self.dim + 1)
            }

        if dim in self._canonical_faces: return self._canonical_faces[dim]

        self.init_canonical_faces(dim, solve_boundary_cells)
        # if not self.is_root_cell:
        #     new_interior_cell_list = []

        #     interior_faces, boundary_faces = self.faces(dim, solve_boundary_cells)

        #     for interior_face in interior_faces:
        #         if interior_face.canonical_image[0] not in new_interior_cell_list:
        #             new_interior_cell_list.append(interior_face.canonical_image[0])

        #     self._canonical_faces[dim] = [new_interior_cell_list,[]]

        #     if solve_boundary_cells:
        #         new_boundary_cell_list = []
        #         for boundary_face in boundary_faces:
        #             if boundary_face.canonical_image[0] not in new_boundary_cell_list:
        #                 new_boundary_cell_list.append(boundary_face.canonical_image[0])

        #         self._canonical_faces[dim][1] = new_boundary_cell_list

        # new_interior_cell_list = []
        # if solve_boundary_cells:
        #     new_boundary_cell_list = []

        # interior_faces, boundary_faces = self.faces(dim, solve_boundary_cells)

        # for interior_face in interior_faces:
        #     is_canonical = True
        #     for interior_cell in new_interior_cell_list:
        #         is_isomorphic, permutation = interior_cell.is_orientation_preserving_isomorphic_to(interior_face)
        #         if is_isomorphic:
        #             interior_face.is_canonical = False
        #             interior_face.canonical_image = (interior_cell,permutation) # may need to fix the direction of map here.
        #             is_canonical = False
        #             break
        #     if is_canonical:
        #         interior_face.set_canonical()
        #         new_interior_cell_list.append(interior_face)

        # if solve_boundary_cells:
        #     for boundary_face in boundary_faces:
        #         is_canonical = True
        #         for boundary_cell in new_boundary_cell_list:
        #             is_isomorphic, permutation = boundary_cell.is_orientation_preserving_isomorphic_to(boundary_face)
        #             if is_isomorphic:
        #                 boundary_face.is_canonical = False
        #                 boundary_face.canonical_image = (boundary_cell,permutation) # may need to fix the direction of map here.
        #                 is_canonical = False
        #                 break
        #         if is_canonical:
        #             boundary_face.set_canonical()
        #             new_boundary_cell_list.append(boundary_face)

        # self._canonical_faces[dim] = [new_interior_cell_list,[]]
        # if solve_boundary_cells:
        #     self._canonical_faces[dim][1] = new_boundary_cell_list
        return self._canonical_faces[dim]

    def canonical_facets(self, solve_boundary_cells = False):
        return self.canonical_faces(self.dim - 1, solve_boundary_cells)
    
    # def solve_all_canonical_faces(self, solve_boundary_cells = False):

    def has_canonical_face(self, face_type: "DFV"): 
        return self.get_canonical_face_isomorphic_to_given_type(face_type) != None
            
    def get_canonical_face_isomorphic_to_given_type(self, face_type: "DFV"):
        # 只写了Interior部分
        for f in self.canonical_faces(face_type.dim)[0]:
            if face_type.is_orientation_preserving_isomorphic_to(f)[0]:
                return f
        return None
    
    def get_faces_isomorphic_to_given_type(self, face_type: "DFV"):
        face_list = []
        for f in self.faces(face_type.dim)[0]:
            if face_type.is_orientation_preserving_isomorphic_to(f)[0]:
                face_list.append(f)
        return face_list
    
    def get_faces_containing_given_cell(self, cell: "DFV", dim=None):
        # 没有适配boundary
        if dim is None:
            return {
                d: self.get_faces_containing_given_cell(cell,d)
                for d in range(cell.dim, self.dim + 1)
            }
        # B是C的边界当且仅当B的label全部被C包含.
        if cell.root_cell != self.root_cell:
            return None
        
        face_list = []
        for f in self.faces(dim)[0]:
            if is_subsequence(cell.edge_labels, f.edge_labels):
                face_list.append(f)

        return face_list

    def return_DFV_tuple(self):
        return (self.D,self.F,self.V)

    def is_orientation_preserving_isomorphic_to(self, dfv, return_map = True):
        return is_orientation_preserving_isomorphic(self, dfv, return_map)

    def set_canonical_image(self,dfv,v_iso,e_iso):
        self.is_canonical = False
        self.canonical_image = (dfv,[v_iso,e_iso])

    def set_canonical(self):
        self.is_canonical = True
        self.canonical_image = (self,DFV.id)

    def planar_graph_automorphism_group(self):
        if self._planar_graph_automorphism_group: return self._planar_graph_automorphism_group

        original_grp = self.D.automorphism_group()
        canonical_faces = [canonical_face(f) for f in self.F]
        group_element_list = []
        for g in original_grp:
            current_map_dict = {v: g(v) for v in self.D.vertices()}
            mapping_func = lambda a: current_map_dict[a]
            is_valid = True
            for fi in self.F:
                if canonical_image(fi,mapping_func) not in canonical_faces:
                    is_valid = False
                    break
            if is_valid:
                group_element_list.append(g)

        self._planar_graph_automorphism_group = PermutationGroup(group_element_list)
        return self._planar_graph_automorphism_group
    
    def automorphism_group(self):
        """Return the effective orientation-preserving symmetry group of this cell."""
        return self.planar_graph_automorphism_group()

    def get_equivalent_coordinates(self):
        """
        Get the equivalence between angle coordinates, labeled by edge labels.
        """
        verts = list(self.poly.vertices())
        n = self.poly.ambient_dim()

        equal_pairs = []

        for i in range(n):
            if self.root_cell.edge_labels[i] in self.edge_labels:
                for j in range(i+1, n):
                    if all(v[i] == v[j] for v in verts):
                        equal_pairs.append((self.root_cell.edge_labels[i], self.root_cell.edge_labels[j]))

        return equal_pairs

    def pointwise_stabilizer(self):
        """Return the subgroup acting trivially on the whole parameter cell."""
        element_list = []
        for g in self.planar_graph_automorphism_group():
            if self.get_fixed_point_set(g).dim() == self.dim:
                element_list.append(g)
        return PermutationGroup(element_list)

    def stabilizer(self):
        """Alias for :meth:`pointwise_stabilizer`."""
        return self.pointwise_stabilizer()

    def stablizer(self):
        """Backward-compatible alias for the former misspelled method name."""
        return self.pointwise_stabilizer()

    def get_fixed_point_set(self,g):
        edge_label_perm= get_edge_label_permutation(self,g)
        orbits = orbits_of_perm(edge_label_perm,self.edge_labels)
        eq_constraints = constraints_from_orbits(orbits, self.root_cell.edge_labels)

        fxpts = Polyhedron(eqns=eq_constraints,base_ring=QQ)
        return self.poly & fxpts
    
    def get_fixed_point_set_list(self, return_g = False):
        if return_g:
            return [(g,self.get_fixed_point_set(g)) for g in self.planar_graph_automorphism_group()]
        return [self.get_fixed_point_set(g) for g in self.planar_graph_automorphism_group()]















# methods for determining if isomorphic
def canonical_face(face):
    n = len(face)
    face_nolabel = [(e[0],e[1]) for e in face]
    rotations = [tuple(face_nolabel[i:] + face_nolabel[:i]) for i in range(n)]
    # reversed_rotations = [tuple((v2, v1) for (v1,v2,_) in reversed(r)) for r in rotations]
    # print(reversed_rotations)
    # return min(rotations + reversed_rotations)
    return min(rotations)

def flipped_canonical_face(face):
    n = len(face)
    face_nolabel_flipped = [(e[1],e[0]) for e in face]
    rotations = [tuple(face_nolabel_flipped[i:] + face_nolabel_flipped[:i]) for i in range(n)]
    return min(rotations)

def canonical_image(face,map):
    image = canonical_face([
        (map(edge[0]),map(edge[1]))
        for edge in face
    ])
    return image

def is_orientation_preserving_isomorphic(dfv1:DFV, dfv2:DFV, return_map = True):
    if not dfv1.V.is_isomorphic(dfv2.V):
        return False, None
    
    # 检查 Delaunay 图是否同构
    is_iso, iso_map = dfv1.D.is_isomorphic(dfv2.D, certificate=True)

    # 一个特判: 把0维退化胞腔丢掉.
    if dfv1.poly.is_empty():
        return True, iso_map

    if not is_iso:
        # print("not iso") 
        return False, None
    else:
        canonical_f2 = [canonical_face(f2) for f2 in dfv2.F]
        
        for phi in dfv1.D.automorphism_group():
            # 构造顶点置换dict
            current_map_dict = {v: iso_map[phi(v)] for v in dfv1.D.vertices()}
            
            mapping_func = lambda a: current_map_dict[a]
            
            is_iso = True
            for f1 in dfv1.F: # 检查所有面是否保持定向
                if all(f2 != canonical_image(f1, mapping_func) for f2 in canonical_f2):
                    is_iso = False
                    break
            
            # 如果找到了保持定向的同构，返回 True 和顶点+边映射字典
            if is_iso: 
                if return_map:
                    edge_label_map = get_edge_label_map(dfv1,dfv2,lambda x: current_map_dict[x])
                    return True, [current_map_dict, edge_label_map]
                else:
                    return True, None
        return False, None

def cyclic_match_position(source_word, target_word):
    """
    找 source_word 与 target_word 的有向循环匹配位移.
    若不存在, 返回 None.

    返回 shift, 使得:
        source_word[k] == target_word[(k + shift) % n]
    """
    n = len(source_word)

    if n != len(target_word):
        return None

    for shift in range(n):
        if all(
            source_word[k] == target_word[(k + shift) % n]
            for k in range(n)
        ):
            return shift

    return None

def get_edge_label_map(dfv1:DFV, dfv2:DFV, mapping):
    """
    给定映射 mapping: dfv1 -> dfv2 求解边映射的 dict. 
    若 mapping 不保持面集 F, 返回None.

    返回 label_perm, 使得
        label_perm[old_label] == new_label.
    """

    target_words = [
        [(u, v) for u, v, _ in face]
        for face in dfv2.F
    ]

    label_perm = {}

    for source_face in dfv1.F:
        mapped_word = [
            (mapping(u), mapping(v))
            for u, v, _ in source_face
        ]

        match = None

        for j, target_word in enumerate(target_words):
            shift = cyclic_match_position(mapped_word, target_word)

            if shift is not None:
                match = (j, shift)
                break

        if match is None:
            return None
            
        j, shift = match
        target_face = dfv2.F[j]
        n = len(source_face)

        for k, (_, _, source_label) in enumerate(source_face):
            target_label = target_face[(k + shift) % n][2]

            if source_label in label_perm:
                if label_perm[source_label] != target_label:
                    return None
            else:
                label_perm[source_label] = target_label

    return label_perm

def get_edge_label_permutation(dfv:DFV,g):
    """
    给定 DFV 和置换 g 求解边置换的 dict. 
    若 perm 不保持面集 F, 返回None.

    返回 label_perm, 使得
        label_perm[old_label] == new_label.

    这里的 label_perm 扩充过, 以适应
    """
    return get_edge_label_map(dfv,dfv,g)


def orbits_of_perm(label_perm, labels):
    """
    label_perm: dict, 表示单个 perm 的作用: label -> label
    labels: 所有 label 列表

    返回: 所有轨道
    """
    seen = []
    orbits = []

    for x in labels:
        if x in seen:
            continue

        orbit = []
        cur = x

        while cur not in orbit:
            orbit.append(cur)
            cur = label_perm[cur]

        orbits.append(orbit)
        seen += orbit

    return orbits

def delaunay_edge_orbit_of_g(dfv:DFV, g):
    """
    给定 DFV 和顶点置换 g 求解边置换的轨道.

    返回: 所有边的轨道
    """
    labels = dfv.edge_labels
    label_perm = get_edge_label_permutation(dfv,g)

    return orbits_of_perm(label_perm,labels)

def constraints_from_orbits(orbits,labels):
    """
    从 orbits 求解 Polyhedron 构造函数所需的线性条件.
    """
    n = len(labels) # amount of parameters
    eq_constraints = [[0]*(n+1)]

    for orbit in orbits:
        for a, b in combinations(orbit, 2):
            eq_constraints.append(
                [0]+[1 if a == label else -1 if b == label else 0 for label in labels]
            )
    
    return eq_constraints