import bpy
import bmesh
import numpy as np

def getBM(name):
    # Get the active mesh
    me = bpy.data.meshes[name]
    # Get a BMesh representation
    bm = bmesh.new()   # create an empty BMesh
    bm.from_mesh(me)   # fill it in from a Mesh
    return bm

def freeBM(bm, name):
    me = bpy.data.meshes[name]
    bm.free() 
    

def freeEditedBM(edited_bm,name):
    me = bpy.data.meshes[name]
    edited_bm.to_mesh(me)
    edited_bm.free() 

#find open boundary , multi boundary find the one with max number edges   BENDG
def is_include(e, exclude):
    for ex in exclude:
        if e in ex:
            return True
    return False

def findOpenEdge(edges, exclude):
    for e in edges:
        if e.is_boundary and not is_include(e, exclude):
            return e
    return None

def nextBoundaryEdge(e,preE):
    stv,nxv = e.verts
    for nxe in nxv.link_edges:
        if nxe.is_boundary and nxe != e and nxe != preE:
            return nxe
    for nxe in stv.link_edges:
        if nxe.is_boundary and nxe != e and nxe != preE:
            return nxe

def boundaryloop(e):
    print(f"start:{e}")
    boundary = [e]
    nxe = nextBoundaryEdge(e, None)
    preE = e
    while nxe != e:
        print("---")
        print(nxe,type(nxe))
        boundary.append(nxe)
        nxe = nextBoundaryEdge(boundary[-1], preE)
        preE = boundary[-1]
    return boundary

def findBoundary(bm):
    e = findOpenEdge(bm.edges, [])
    bds = []
    while e != None:
        bds.append(boundaryloop(e))
        e = findOpenEdge(bm.edges, bds)
    return bds
#find open boundary , multi boundary find the one with max number edges END

def boundaryVerts(boundaryEdges):
    stv,nxv = boundaryEdges[0].verts
    verts = [stv.index]
    for e in boundaryEdges[0:-1]:
        other = e.other_vert(stv)
        verts.append(other.index)
        stv = other
    return verts

def bdE2bdV(bdE):
    bdv = []
    for bde in bdE:
        bdv.append(boundaryVerts(bde))
    return bdv

def get_other_face_vert(v0,v1,f):
    for v in f.verts:
        if v not in [v0,v1]:
            return v

def extract(vs,edited_vs):
    vs_num = len(vs)
    connected = []
    xyz = np.zeros((vs_num,3))
    edited_xyz = np.zeros((vs_num,3))
    for v in vs:
        xyz[v.index] = v.co.x, v.co.y, v.co.z
        connected.append([[e.index,e.other_vert(v).index] for e in v.link_edges])
    for v in edited_vs:
        edited_xyz[v.index] = v.co.x, v.co.y, v.co.z
    return xyz, edited_xyz, connected

def calcWeight(eds):
    ed_num = len(eds)
    cotangent = np.ones(ed_num)
    # for e in eds:
    #     if not e.is_boundary:
    #         v0,v1 = e.verts
    #         f0,f1 = e.link_faces
    #         v2 = get_other_face_vert(v0,v1,f0)
    #         v3 = get_other_face_vert(v0,v1,f1)
    #         vec0 = v0.co-v2.co
    #         vec1 = v1.co-v2.co
    #         vec00 = v0.co-v3.co
    #         vec11 = v1.co-v3.co
    #         cotangent[e.index] = 0.5 * (vec0.dot(vec1)/vec0.cross(vec1).length + vec00.dot(vec11)/vec00.cross(vec11).length)
    return  cotangent#,local_average_area

def buildLaplacian(connected, cotangent, constrained):
    vs_num = len(connected)
    L = np.zeros((vs_num,vs_num))
    for i,con in enumerate(connected):
        if i in constrained:
            L[i,i] = 1
        else:
            w_sum = 0
            for ei,j in con:
                L[i,j] = -cotangent[ei]
                w_sum += cotangent[ei]
            L[i,i] = w_sum 
    return L

def buildB(xyz, connected, cotangent, R, constrained):
    vs_num = len(xyz)
    B = np.zeros((vs_num,3))
    for i,v in enumerate(xyz):
        if i in constrained:
            B[i] = xyz[i]
        else:
            for ei,j in connected[i]:
                eij = (v - xyz[j])
                JR = (R[i] + R[j])
                w = 0.5*cotangent[ei]
                B[i] += w*(JR@eij)
    return B
#fix  fixed vertex are fixed in position in the process
#handles  move by hand vertex  are fixed in position in the process
def buildfix(xyz, fix):
    vs_num = len(xyz)
    fix_num = len(fix)
    F = np.zeros((fix_num,vs_num))
    FB = np.zeros((fix_num,3))
    for i,f in enumerate(fix):
        F[i,f] = 1
        FB[i] = xyz[f]
    return F,FB

def localStep(xyz, fixed, R, cotangent, edited_xyz, connected):
    for i,v in enumerate(xyz):
        if not i in fixed:
            E = []
            EE = []
            for ei,j in connected[i]:
                vec = (v - xyz[j])*cotangent[ei]
                E.append(vec)
                vecc = (edited_xyz[i] - edited_xyz[j])
                EE.append(vecc)
            eij = np.array(E).transpose()
            eeij = np.array(EE)
            Si = eij@eeij
            U, _, Vh = np.linalg.svd(Si)
            rot = (U@Vh).transpose() #Vi * Ui.transpose()
            if np.linalg.det(rot) < 0:
                print('flipped...')
                minus = np.identity(3)
                minus[2,2] = -1
                R[i] = Vh.transpose()@(U@minus).transpose()
            else:
                R[i] = rot
            
def setVertexPos(constrained, edited_vs, edited_xyz):
    for v in edited_vs:
        if not v.index in constrained:
            v.co.x = edited_xyz[v.index,0]
            v.co.y = edited_xyz[v.index,1]
            v.co.z = edited_xyz[v.index,2]
    
edited_vert_index = 60
bm = getBM("Bunny_head_earR")
edited_bm = getBM("Bunny_head_earR.001")
bde = findBoundary(bm)
bdv = bdE2bdV(bde)
fixed = bdv[0]
handled = [edited_vert_index]
constrained = fixed+handled
identity = [[1.,0,0],[0,1.,0],[0,0,1.]]
xyz, edited_xyz, connected = extract(bm.verts,edited_bm.verts)
vs_num = len(bm.verts)
R = np.array(identity*vs_num).reshape(-1,3,3)
cotangent = calcWeight(bm.edges)
L = buildLaplacian(connected,cotangent,constrained)
for i in range(10):
    localStep(xyz, fixed, R, cotangent, edited_xyz, connected)
    B = buildB(xyz, connected, cotangent, R, constrained)
    for vi in handled:
        B[vi] = edited_xyz[vi]
    edited_xyz = np.linalg.solve(L,B)
setVertexPos(constrained, edited_bm.verts, edited_xyz)
freeBM(bm, "Bunny_head_earR")
freeEditedBM(edited_bm, "Bunny_head_earR.001")